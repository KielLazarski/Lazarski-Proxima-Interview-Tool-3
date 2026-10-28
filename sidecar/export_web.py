#!/usr/bin/env python3
"""Write the frozen cache the dashboard reads.

Two kinds of output land in cache/web:

* one file per system, holding everything a page needs to draw it -- residue
  positions, the contact list with its frustration values, and the interface
  and proximity labels;
* `crossval.json`, a small fixture of inputs plus the Python engine's answers,
  which the JavaScript engine is checked against. That file is the reason the
  architecture's "validates" arrow is a real check rather than a claim.

Parallel arrays rather than arrays of objects, because the contact lists run
to a few thousand entries per system and the difference in file size is
roughly fourfold.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from glueforge.cavity import detect_cavity
from glueforge.cif import read_structure
from glueforge.frustration import (
    FrustrationConfig,
    compute_frustration,
    representative_coords,
)
from glueforge.geometry import (
    atom_sasa,
    contact_midpoints,
    exit_vectors,
    reach_frustrated_contacts,
    secondary_structure,
)
from glueforge.interface import GLUE_PROXIMAL, assign_roles, label_contacts
from glueforge.metrics import (
    compare_groups,
    compare_stratified,
    distance_gradient,
)
from glueforge.potential import AwsemParams, load_potential

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "cache" / "raw"
WEB = ROOT / "cache" / "web"

#: Which biological system each entry belongs to, so the dashboard can group
#: replicates and avoid presenting them as independent evidence.
SYSTEM_GROUPS = {
    "6TD3": ("CDK12-cyclin K x DDB1", "Faust/Thoma form"),
    "8BU5": ("CDK12-cyclin K x DDB1", "Faust/Thoma form"),
    "8BU7": ("CDK12-cyclin K x DDB1", "Faust/Thoma form"),
    "8BU9": ("CDK12-cyclin K x DDB1", "Faust/Thoma form"),
    "8BUA": ("CDK12-cyclin K x DDB1", "Faust/Thoma form"),
    "8BUJ": ("CDK12-cyclin K x DDB1", "Faust/Thoma form"),
    "8BUK": ("CDK12-cyclin K x DDB1", "Faust/Thoma form"),
    "8BUL": ("CDK12-cyclin K x DDB1", "Faust/Thoma form"),
    "6Q0W": ("DCAF15 x RBM39", "Faust form"),
    "6Q0R": ("DCAF15 x RBM39", "Faust form"),
    "6Q0V": ("DCAF15 x RBM39", "Faust form"),
    "6UD7": ("DCAF15 x RBM39", "Bussiere form"),
    "6H0F": ("CRBN x C2H2 zinc finger", "IKZF1 ZF2"),
    "5FQD": ("CRBN x C2H2 zinc finger", "CK1a"),
}


# _cavity_blob moved to glueforge/viewblob.py, so that the benchmark
# export and a live docking run build the same cavity shape from the
# same code. The page has one cavity view and it reads one set of keys.
from glueforge.viewblob import cavity_blob as _cavity_blob  # noqa: E402


def export_system(path: Path, pot, n_decoys: int, proximal: float) -> dict:
    structure = read_structure(path)
    assign = assign_roles(structure)
    residues = [
        r for r in structure.residues
        if r.chain in (assign.anchor_chain, assign.partner_chain)
    ]

    cfg = FrustrationConfig(n_decoys=n_decoys, flavour="configurational")
    fr = compute_frustration(residues, pot, cfg)
    labels = label_contacts(residues, fr.index_i, fr.index_j, assign, proximal)

    xyz = representative_coords(residues)
    centre = xyz.mean(axis=0)
    xyz = xyz - centre

    iface = labels.mask_interface
    prox = labels.mask_proximal & iface
    dist = iface & ~labels.mask_proximal
    rho_pair = 0.5 * (fr.density[fr.index_i] + fr.density[fr.index_j])

    comparison = compare_groups(fr.frustration[prox], fr.frustration[dist])
    strat = compare_stratified(
        fr.frustration[iface], prox[iface], rho_pair[iface]
    )
    contact_dlig = np.minimum(
        labels.residue_min_dist[fr.index_i], labels.residue_min_dist[fr.index_j]
    )
    gradient = distance_gradient(fr.frustration[iface], contact_dlig[iface])

    lig = assign.ligand
    heavy = [a for a in lig.atoms if a.element.upper() != "H"]
    lig_xyz = np.asarray(
        [(a.x, a.y, a.z) for a in heavy], dtype=np.float64
    ) - centre

    group, form = SYSTEM_GROUPS.get(structure.entry_id, ("unassigned", ""))

    def r3(a) -> list:
        arr = np.asarray(a, dtype=np.float64)
        # NaN and inf are not valid JSON, and this file is parsed by the page.
        arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
        return np.round(arr, 3).tolist()

    # --- drawing extras: CA trace, secondary structure, interface sidechains
    ca = []
    for r in residues:
        a = r.atom("CA") or r.representative()
        ca.append((a.x, a.y, a.z))
    ca_xyz = np.asarray(ca, dtype=np.float64) - centre
    ss = secondary_structure(residues)

    iface_res = np.zeros(len(residues), dtype=bool)
    for k in np.nonzero(iface)[0]:
        iface_res[fr.index_i[k]] = True
        iface_res[fr.index_j[k]] = True

    side_idx, side_x, side_y, side_z, side_name = [], [], [], [], []
    backbone = {"N", "CA", "C", "O", "OXT"}
    for ri in np.nonzero(iface_res)[0]:
        for atom in residues[ri].atoms:
            if atom.element.upper() == "H" or atom.name in backbone:
                continue
            side_idx.append(int(ri))
            side_x.append(atom.x - centre[0])
            side_y.append(atom.y - centre[1])
            side_z.append(atom.z - centre[2])
            side_name.append(atom.name)

    # --- step 4 scoring half: the composite cavity and its fields ---------
    cav = detect_cavity(
        [a for r in residues for a in r.atoms],
        np.asarray([(x.x, x.y, x.z) for x in heavy], dtype=np.float64).mean(axis=0),
        {i for i, at in enumerate(a for r in residues for a in r.atoms)
         if at.chain == assign.anchor_chain},
    )

    # --- step 2 and step 6: exposed face, exit vectors, what they reach ----
    # Exposure is computed against the ANCHOR ALONE, not the ternary complex.
    # In the ternary the glue is sandwiched between both proteins and is
    # buried by definition, which would report no growth points at all. The
    # face that matters is the one presented to solvent once the ligand sits
    # in the anchor's pocket, because that is the only part a second protein
    # can ever touch and the only place chemistry can grow from.
    anchor_atoms = [
        a for r in residues if r.chain == assign.anchor_chain for a in r.atoms
    ]
    binary_context = anchor_atoms + lig.atoms
    vecs = exit_vectors(lig, binary_context)
    mids = contact_midpoints(residues, fr.index_i, fr.index_j)
    reach = reach_frustrated_contacts(
        vecs, heavy, mids, fr.frustration, iface
    )
    lig_sasa = atom_sasa(heavy, binary_context)
    # How much of that presented face the partner then buries -- the
    # structural statement that a glue contributes only a sliver of interface.
    ternary_context = [a for r in residues for a in r.atoms] + lig.atoms
    lig_sasa_ternary = atom_sasa(heavy, ternary_context)

    return {
        "entry": structure.entry_id,
        "group": group,
        "form": form,
        "assignment": assign.as_dict(),
        "potential": pot.name,
        "potential_published": pot.is_published(),
        "n_decoys": n_decoys,
        "proximal_cutoff": proximal,
        "residues": {
            "chain": [r.chain for r in residues],
            "seq": [r.seq_id for r in residues],
            "one": "".join(r.one for r in residues),
            "x": r3(xyz[:, 0]),
            "y": r3(xyz[:, 1]),
            "z": r3(xyz[:, 2]),
            "rho": r3(fr.density),
            "dlig": r3(np.minimum(labels.residue_min_dist, 99.0)),
            "cax": r3(ca_xyz[:, 0]),
            "cay": r3(ca_xyz[:, 1]),
            "caz": r3(ca_xyz[:, 2]),
            "ss": ss,
            "iface": iface_res.astype(np.uint8).tolist(),
        },
        "sidechains": {
            "note": "heavy sidechain atoms, interface residues only",
            "res": side_idx,
            "x": r3(side_x),
            "y": r3(side_y),
            "z": r3(side_z),
            "name": side_name,
        },
        "ligand": {
            "comp_id": lig.comp_id,
            "x": r3(lig_xyz[:, 0]),
            "y": r3(lig_xyz[:, 1]),
            "z": r3(lig_xyz[:, 2]),
            "element": [a.element for a in heavy],
            "name": [a.name for a in heavy],
            "sasa": r3(lig_sasa),
            "sasa_ternary": r3(lig_sasa_ternary),
            "sasa_total": round(float(lig_sasa.sum()), 1),
            "sasa_total_ternary": round(float(lig_sasa_ternary.sum()), 1),
            "buried_by_partner": round(
                float(1.0 - lig_sasa_ternary.sum() / max(lig_sasa.sum(), 1e-6)), 3
            ),
        },
        "cavity": _cavity_blob(cav, centre, r3),
        "exit_vectors": {
            "note": (
                "Exposed ligand heavy atoms and their outward normals, with "
                "the frustrated interface contacts each could reach. Reach is "
                "a distance and angle feasibility test, not a docking result."
            ),
            "exposure_context": "anchor chain only (binary complex)",
            "max_reach": 12.0,
            "max_angle_deg": 75.0,
            "frustrated_below": -0.3,
            "vectors": [v.as_dict() for v in reach],
        },
        "contacts": {
            "i": fr.index_i.tolist(),
            "j": fr.index_j.tolist(),
            "d": r3(fr.distance),
            "e": r3(fr.energy),
            "f": r3(fr.frustration),
            "interface": iface.astype(np.uint8).tolist(),
            "proximal": labels.mask_proximal.astype(np.uint8).tolist(),
        },
        "result": {
            "contacts_total": len(fr),
            "contacts_interface": int(iface.sum()),
            "contacts_proximal": int(prox.sum()),
            "contacts_distal": int(dist.sum()),
            "comparison": comparison.as_dict(),
            "stratified": strat.as_dict(),
            "gradient": gradient,
            "mean_density_proximal": round(float(rho_pair[prox].mean()), 3)
            if prox.any() else None,
            "mean_density_distal": round(float(rho_pair[dist].mean()), 3)
            if dist.any() else None,
            "whole_complex": fr.summary(),
        },
    }


def build_crossval(path: Path, pot, n_residues: int = 120) -> dict:
    """A small fixture the JavaScript engine must reproduce exactly."""
    structure = read_structure(path)
    assign = assign_roles(structure)
    residues = [
        r for r in structure.residues
        if r.chain in (assign.anchor_chain, assign.partner_chain)
    ]
    # Take residues closest to the ligand so the subset is a real interface
    # rather than an arbitrary slice of two proteins.
    from glueforge.interface import min_distance_to_ligand

    d = min_distance_to_ligand(residues, assign.ligand)
    keep = np.argsort(d)[:n_residues]
    keep.sort()
    subset = [residues[k] for k in keep]

    cfg = FrustrationConfig(n_decoys=200, flavour="configurational")
    fr = compute_frustration(subset, pot, cfg)
    cfg_m = FrustrationConfig(n_decoys=200, flavour="mutational")
    fm = compute_frustration(subset, pot, cfg_m)

    xyz = representative_coords(subset)
    return {
        "source": structure.entry_id,
        "note": (
            "Inputs plus the Python engine's outputs. The JavaScript engine "
            "must reproduce these to within 1e-9; both draw decoys from the "
            "same mulberry32 stream in the same order."
        ),
        "potential": pot.as_dict(),
        "config": {
            "n_decoys": 200,
            "seed": cfg.seed,
            "min_sep_contact": cfg.min_sep_contact,
            "contact_cutoff": cfg.contact_cutoff,
        },
        "input": {
            "chain": [r.chain for r in subset],
            "one": "".join(r.one for r in subset),
            "x": np.round(xyz[:, 0], 6).tolist(),
            "y": np.round(xyz[:, 1], 6).tolist(),
            "z": np.round(xyz[:, 2], 6).tolist(),
        },
        "expected": {
            "density": np.round(fr.density, 9).tolist(),
            "i": fr.index_i.tolist(),
            "j": fr.index_j.tolist(),
            "energy": np.round(fr.energy, 9).tolist(),
            "configurational": np.round(fr.frustration, 9).tolist(),
            "mutational": np.round(fm.frustration, 9).tolist(),
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("entries", nargs="*")
    ap.add_argument("--gamma", default=None)
    ap.add_argument("--decoys", type=int, default=800)
    ap.add_argument("--proximal", type=float, default=GLUE_PROXIMAL)
    args = ap.parse_args()

    pot = load_potential(args.gamma, AwsemParams())
    WEB.mkdir(parents=True, exist_ok=True)

    paths = (
        [RAW / f"{e.upper()}.cif" for e in args.entries]
        if args.entries else sorted(RAW.glob("*.cif"))
    )

    index = []
    for p in paths:
        if not p.exists():
            continue
        try:
            rec = export_system(p, pot, args.decoys, args.proximal)
        except Exception as exc:  # noqa: BLE001
            print(f"  {p.stem}: FAILED {type(exc).__name__}: {exc}")
            continue
        dest = WEB / f"{rec['entry']}.json"
        blob = json.dumps(rec, separators=(",", ":"), allow_nan=False)
        dest.write_text(blob)
        index.append({
            "entry": rec["entry"],
            "group": rec["group"],
            "form": rec["form"],
            "ligand": rec["assignment"]["ligand"],
            "anchor": rec["assignment"]["anchor_name"],
            "partner": rec["assignment"]["partner_name"],
            "anchor_share": rec["assignment"]["anchor_share"],
            "result": rec["result"]["comparison"],
            "stratified": rec["result"]["stratified"]["pooled_difference"],
            "gradient_rho": rec["result"]["gradient"]["spearman_rho"],
            "gradient_p": rec["result"]["gradient"]["p_value"],
            "sweep_stable": rec["result"]["gradient"].get("sweep_all_same_sign"),
            "n_interface": rec["result"]["contacts_interface"],
            "n_proximal": rec["result"]["contacts_proximal"],
            "cavity_volume": rec["cavity"].get("summary", {}).get("volume"),
            "enclosure_quality": rec["cavity"].get("summary", {}).get("enclosure_quality"),
            "anchor_lined_share": rec["cavity"].get("anchor_lined_share"),
            "n_distal": rec["result"]["contacts_distal"],
            "bytes": dest.stat().st_size,
        })
        print(f"  {rec['entry']}: {dest.stat().st_size/1024:.0f} kB")

    sens_path = ROOT / "cache" / "cavity" / "cavity.json"
    sensitivity = {}
    if sens_path.exists():
        blob = json.loads(sens_path.read_text())
        for row in blob.get("sensitivity", []):
            sensitivity[row["entry"]] = row
        for entry, row in sensitivity.items():
            f = WEB / f"{entry}.json"
            if f.exists():
                rec = json.loads(f.read_text())
                rec.setdefault("cavity", {})["sensitivity"] = row
                f.write_text(json.dumps(rec, separators=(",", ":"), allow_nan=False))
        print(f"  merged cavity sensitivity for {len(sensitivity)} systems")

    (WEB / "index.json").write_text(json.dumps({
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "arm": "oracle",
        "potential": pot.name,
        "potential_published": pot.is_published(),
        "potential_provenance": pot.provenance,
        "n_decoys": args.decoys,
        "proximal_cutoff": args.proximal,
        "systems": index,
    }, indent=2))

    if paths:
        cv = build_crossval(paths[0], pot)
        (WEB / "crossval.json").write_text(json.dumps(cv, separators=(",", ":")))
        print(f"  crossval.json: {(WEB/'crossval.json').stat().st_size/1024:.0f} kB "
              f"({len(cv['expected']['i'])} contacts from {cv['source']})")

    print(f"\nwrote {len(index)} systems to {WEB}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
