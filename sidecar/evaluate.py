#!/usr/bin/env python3
"""Evaluate one new pair.

Takes a structure you supply rather than one of the benchmark entries, runs
the whole analysis on it, and optionally adds it to the dashboard.

    # a ternary complex you already have, in mmCIF or PDB
    python3 sidecar/evaluate.py my_complex.pdb

    # or pull one from the PDB by id
    python3 sidecar/evaluate.py 7XYZ --fetch

    # name the ligand explicitly if the automatic pick is wrong
    python3 sidecar/evaluate.py my_complex.cif --ligand LVY

    # override the chain assignment
    python3 sidecar/evaluate.py my_complex.cif --anchor A --partner C

    # add it to the dashboard, then rebuild
    python3 sidecar/evaluate.py my_complex.cif --add-to-dashboard
    node build.js

What this needs from you: a structure in which the two proteins are already
positioned relative to each other, with the ligand bound. That is the one
thing the pipeline cannot work out for itself.

If you have the anchor with its ligand and the partner as separate files,
you need a docking step first -- see the note printed at the end. If you have
only sequences, you need structure prediction before that. Neither is in this
repo, and both are the reason the prospective arm is unbuilt.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))

from glueforge.cavity import detect_cavity
from glueforge.cif import Ligand, Structure, read_any
from glueforge.frustration import FrustrationConfig, compute_frustration
from glueforge.geometry import (
    atom_sasa,
    contact_midpoints,
    exit_vectors,
    reach_frustrated_contacts,
)
from glueforge.interface import (
    GLUE_PROXIMAL,
    Assignment,
    assign_roles,
    label_contacts,
    ligand_chain_contacts,
)
from glueforge.metrics import compare_groups, compare_stratified, distance_gradient
from glueforge.potential import AwsemParams, load_potential

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "cache" / "raw"
OUT = ROOT / "cache" / "evaluated"
RCSB = "https://files.rcsb.org/download/{}.cif"


def fetch(entry: str) -> Path:
    RAW.mkdir(parents=True, exist_ok=True)
    dest = RAW / f"{entry.upper()}.cif"
    if dest.exists():
        return dest
    with urllib.request.urlopen(RCSB.format(entry.upper()), timeout=90) as r, \
            open(dest, "wb") as fh:
        shutil.copyfileobj(r, fh)
    return dest


def manual_assignment(
    structure: Structure, anchor: str, partner: str, ligand_comp: str | None
) -> Assignment:
    """Build an assignment from chains the user named."""
    want = ligand_comp
    if not want:
        best = structure.largest_ligand()
        if best is None:
            raise SystemExit(
                "no ligand found. Candidates in this file: "
                + (", ".join(sorted({l.comp_id for l in structure.ligands}))
                   or "none")
            )
        want = best.comp_id

    # With several copies in the asymmetric unit, take the one that actually
    # touches the chains the user named rather than the largest overall.
    copies = [l for l in structure.ligands if l.comp_id.upper() == want.upper()]
    if not copies:
        raise SystemExit(
            f"no ligand with id {want}. Candidates: "
            + (", ".join(sorted({l.comp_id for l in structure.ligands}))
               or "none")
        )
    scored = []
    for cand in copies:
        c = ligand_chain_contacts(structure, cand)
        scored.append((c.get(anchor, 0) + c.get(partner, 0), cand, c))
    scored.sort(key=lambda t: -t[0])
    total, lig, counts = scored[0]
    if total == 0:
        raise SystemExit(
            f"no copy of {want} contacts chains {anchor} and {partner}.\n"
            f"  copies of {want} are near: "
            + "; ".join(
                f"{cand.chain}:{sorted(c)}" for _, cand, c in scored[:4])
            + "\n  drop --anchor/--partner to let the roles be assigned "
              "automatically."
        )
    return Assignment(
        entry_id=structure.entry_id, ligand=lig,
        anchor_chain=anchor, partner_chain=partner,
        anchor_name=structure.name_of(anchor),
        partner_name=structure.name_of(partner),
        anchor_contacts=counts.get(anchor, 0),
        partner_contacts=counts.get(partner, 0),
        other_chains={k: v for k, v in counts.items()
                      if k not in (anchor, partner)},
    )


def evaluate(
    path: Path,
    ligand_comp: str | None = None,
    anchor: str | None = None,
    partner: str | None = None,
    n_decoys: int = 800,
    proximal: float = GLUE_PROXIMAL,
    gamma: str | None = None,
) -> dict:
    structure = read_any(path)
    if not structure.residues:
        raise SystemExit(f"{path.name}: no protein residues parsed")

    if anchor and partner:
        assign = manual_assignment(structure, anchor, partner, ligand_comp)
    else:
        try:
            assign = assign_roles(structure, ligand_comp)
        except ValueError as exc:
            raise SystemExit(
                f"{path.name}: {exc}\n"
                f"  chains present: {sorted(structure.chain_sizes())}\n"
                f"  ligands present: {sorted({l.comp_id for l in structure.ligands})}\n"
                f"  try --anchor/--partner, or --ligand to name one."
            )

    residues = [
        r for r in structure.residues
        if r.chain in (assign.anchor_chain, assign.partner_chain)
    ]
    if len(residues) < 20:
        raise SystemExit(
            f"only {len(residues)} residues in chains "
            f"{assign.anchor_chain}/{assign.partner_chain}"
        )

    pot = load_potential(gamma, AwsemParams())
    fr = compute_frustration(
        residues, pot, FrustrationConfig(n_decoys=n_decoys))
    labels = label_contacts(residues, fr.index_i, fr.index_j, assign, proximal)

    iface = labels.mask_interface
    prox = labels.mask_proximal & iface
    dist = iface & ~labels.mask_proximal
    if iface.sum() < 10:
        raise SystemExit(
            f"only {int(iface.sum())} contacts cross the interface; the two "
            f"chains may not be in contact in this structure"
        )

    rho_pair = 0.5 * (fr.density[fr.index_i] + fr.density[fr.index_j])
    contact_dlig = np.minimum(
        labels.residue_min_dist[fr.index_i], labels.residue_min_dist[fr.index_j])

    comparison = compare_groups(fr.frustration[prox], fr.frustration[dist])
    strat = compare_stratified(fr.frustration[iface], prox[iface], rho_pair[iface])
    gradient = distance_gradient(fr.frustration[iface], contact_dlig[iface])

    # Cavity and growth points.
    all_atoms = [a for r in residues for a in r.atoms]
    anchor_idx = {i for i, a in enumerate(all_atoms)
                  if a.chain == assign.anchor_chain}
    heavy = [a for a in assign.ligand.atoms if a.element.upper() != "H"]
    lig_xyz = np.array([(a.x, a.y, a.z) for a in heavy])
    cav = detect_cavity(all_atoms, lig_xyz.mean(axis=0), anchor_idx)

    recovered = (
        float((cKDTree(cav.points).query(lig_xyz, k=1)[0] < 1.2).mean())
        if len(cav.points) else 0.0
    )

    anchor_atoms = [a for r in residues if r.chain == assign.anchor_chain
                    for a in r.atoms]
    vecs = exit_vectors(assign.ligand, anchor_atoms + assign.ligand.atoms)
    mids = contact_midpoints(residues, fr.index_i, fr.index_j)
    reach = reach_frustrated_contacts(vecs, heavy, mids, fr.frustration, iface)
    sasa_binary = atom_sasa(heavy, anchor_atoms + assign.ligand.atoms)
    sasa_ternary = atom_sasa(heavy, all_atoms + assign.ligand.atoms)

    return {
        "source": str(path),
        "entry": structure.entry_id,
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "potential": pot.name,
        "potential_published": pot.is_published(),
        "assignment": assign.as_dict(),
        "chains_in_file": structure.chain_sizes(),
        "ligands_in_file": sorted({l.comp_id for l in structure.ligands}),
        "frustration": {
            "contacts_total": len(fr),
            "contacts_interface": int(iface.sum()),
            "contacts_proximal": int(prox.sum()),
            "contacts_distal": int(dist.sum()),
            "comparison": comparison.as_dict(),
            "stratified": strat.as_dict(),
            "gradient": gradient,
        },
        "cavity": {
            **cav.summary(),
            "ligand_recovered": round(recovered, 3),
            "anchor_lined_share": round(float((cav.from_anchor >= 0.6).mean()), 3)
            if len(cav.points) else None,
        },
        "presented_face": {
            "sasa_binary": round(float(sasa_binary.sum()), 1),
            "sasa_ternary": round(float(sasa_ternary.sum()), 1),
            "buried_by_partner": round(
                float(1 - sasa_ternary.sum() / max(sasa_binary.sum(), 1e-6)), 3),
            "growth_points": len(reach),
            "top_vectors": [v.as_dict() for v in reach[:5]],
        },
    }


def report(r: dict) -> None:
    a, f, c, p = (r["assignment"], r["frustration"], r["cavity"],
                  r["presented_face"])
    print(f"\n{'=' * 64}")
    print(f"{r['entry']}   ligand {a['ligand']}   potential {r['potential']}"
          + ("" if r["potential_published"] else "  (STAND-IN)"))
    print(f"{'=' * 64}")
    print(f"  anchor   {a['anchor_chain']}  {a['anchor_name'][:44]}")
    print(f"  partner  {a['partner_chain']}  {a['partner_name'][:44]}")
    print(f"  the ligand puts {100 * a['anchor_share']:.0f}% of its protein "
          f"contacts into the anchor")

    print(f"\n  FRUSTRATION")
    cmp_ = f["comparison"]
    print(f"    interface contacts      {f['contacts_interface']} "
          f"({f['contacts_proximal']} glue-proximal, {f['contacts_distal']} distal)")
    print(f"    proximal vs distal      {cmp_['difference']:+.3f}  "
          f"(p = {cmp_['p_value']})")
    print(f"    burial-controlled       "
          f"{f['stratified']['pooled_difference']:+.3f}"
          if f["stratified"]["pooled_difference"] is not None else
          "    burial-controlled       n/a")
    g = f["gradient"]
    if g.get("spearman_rho") is not None:
        print(f"    gradient toward ligand  {g['spearman_rho']:+.3f} "
              f"(p = {g['p_value']})")
        verdict = ("consistent with the benchmark"
                   if g["spearman_rho"] > 0 and cmp_["difference"] < 0
                   else "does NOT follow the benchmark pattern")
        print(f"    -> {verdict}")

    print(f"\n  CAVITY")
    if c.get("points"):
        print(f"    volume                  {c['volume']:.0f} A^3 "
              f"({c['points']} grid points)")
        print(f"    enclosure quality       {c['enclosure_quality']:.2f}  "
              + ("an enclosed pocket" if c["enclosure_quality"] >= 0.4
                 else "closer to a surface groove"))
        print(f"    anchor-lined            {100 * c['anchor_lined_share']:.0f}%")
        print(f"    mean electrostatic      {c['mean_electrostatic']:+.2f}  "
              + ("wants a cation or donor" if c["mean_electrostatic"] < 0
                 else "wants an anion"))
        print(f"    ligand recovered        {100 * c['ligand_recovered']:.0f}%"
              + ("" if c["ligand_recovered"] > 0.8
                 else "   <- low; check the site is really a pocket"))
    else:
        print("    none detected at this site")

    print(f"\n  PRESENTED FACE AND GROWTH POINTS")
    print(f"    face on the anchor      {p['sasa_binary']:.0f} A^2, "
          f"{100 * p['buried_by_partner']:.0f}% buried by the partner")
    print(f"    growth points           {p['growth_points']}")
    for v in p["top_vectors"][:3]:
        print(f"      {v['atom_name']:>5s} ({v['element']})  reaches "
              f"{v['n_reachable']:2d} frustrated contacts, "
              f"addressable {v['addressable']:.2f}")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Evaluate one protein pair with a bound ligand.")
    ap.add_argument("structure", help="a .cif/.pdb file, or a PDB id with --fetch")
    ap.add_argument("--fetch", action="store_true",
                    help="treat the argument as a PDB id and download it")
    ap.add_argument("--ligand", default=None,
                    help="component id of the glue, if the automatic pick is wrong")
    ap.add_argument("--anchor", default=None, help="anchor chain id")
    ap.add_argument("--partner", default=None, help="partner chain id")
    ap.add_argument("--decoys", type=int, default=800)
    ap.add_argument("--proximal", type=float, default=GLUE_PROXIMAL)
    ap.add_argument("--gamma", default=None,
                    help="path to gamma tables, or 'builtin' for the stand-in")
    ap.add_argument("--add-to-dashboard", action="store_true",
                    help="copy into cache/raw so export_web.py picks it up")
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()

    if bool(args.anchor) != bool(args.partner):
        raise SystemExit("--anchor and --partner must be given together")

    path = fetch(args.structure) if args.fetch else Path(args.structure)
    if not path.exists():
        raise SystemExit(f"{path} not found")

    result = evaluate(path, args.ligand, args.anchor, args.partner,
                      args.decoys, args.proximal, args.gamma)
    report(result)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / f"{result['entry']}.json"
    dest.write_text(json.dumps(result, indent=2))
    print(f"\n  wrote {dest}")

    if args.add_to_dashboard:
        RAW.mkdir(parents=True, exist_ok=True)
        target = RAW / f"{result['entry']}{path.suffix}"
        if path.resolve() != target.resolve():
            shutil.copy(path, target)
        print(f"  copied to {target}")
        print(f"  now run:  python3 sidecar/export_web.py && node build.js")

    print(
        "\n  Note: this assumes the two chains are already positioned "
        "correctly.\n"
        "  If you have the anchor-plus-ligand and the partner as separate\n"
        "  structures, dock them first - the cavity holds up to about 2 A of\n"
        "  partner error and falls apart by 3 A (see the cavity page).\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
