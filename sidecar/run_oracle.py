#!/usr/bin/env python3
"""Phase 1: the oracle run.

For each system, takes the crystal ternary complex, deletes the glue, computes
contact frustration across the anchor-partner interface, and asks whether the
contacts where the glue used to sit are more frustrated than the rest of that
interface.

This is the upper bound, not a prediction: the interface geometry comes from
the structure that already contains the answer. It tests whether the
hypothesis is true at all, which is a different question from whether the
pipeline can find it blind. The prospective arm is the other half.

Usage:
    python run_oracle.py                     # every system in the cache
    python run_oracle.py 6TD3 6Q0W           # named entries only
    python run_oracle.py --gamma path.json   # use the published AWSEM tables
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from glueforge.cif import read_structure
from glueforge.frustration import FrustrationConfig, compute_frustration
from glueforge.interface import GLUE_PROXIMAL, assign_roles, label_contacts
from glueforge.metrics import compare_groups, compare_stratified
from glueforge.potential import AwsemParams, load_potential

RAW = Path(__file__).resolve().parent.parent / "cache" / "raw"
OUT = Path(__file__).resolve().parent.parent / "cache" / "oracle"


def run_one(
    path: Path,
    pot,
    n_decoys: int,
    flavours: tuple[str, ...] = ("configurational", "mutational"),
    proximal: float = GLUE_PROXIMAL,
) -> dict:
    t0 = time.time()
    structure = read_structure(path)
    assign = assign_roles(structure)

    residues = [
        r for r in structure.residues
        if r.chain in (assign.anchor_chain, assign.partner_chain)
    ]

    record: dict = {
        "entry": structure.entry_id,
        "source_file": path.name,
        "assignment": assign.as_dict(),
        "n_residues": len(residues),
        "potential": pot.name,
        "potential_published": pot.is_published(),
        "proximal_cutoff": proximal,
        "flavours": {},
    }

    for flavour in flavours:
        cfg = FrustrationConfig(n_decoys=n_decoys, flavour=flavour)
        fr = compute_frustration(residues, pot, cfg)
        labels = label_contacts(
            residues, fr.index_i, fr.index_j, assign, proximal
        )

        iface = labels.mask_interface
        prox = labels.mask_proximal & iface
        dist = iface & ~labels.mask_proximal

        f_prox = fr.frustration[prox]
        f_dist = fr.frustration[dist]

        comparison = compare_groups(f_prox, f_dist)

        # Burial control: the glue sits in the interface core, so compare
        # within bands of comparable contact density.
        rho_pair = 0.5 * (fr.density[fr.index_i] + fr.density[fr.index_j])
        strat = compare_stratified(
            fr.frustration[iface], prox[iface], rho_pair[iface]
        )

        record["flavours"][flavour] = {
            "contacts_total": len(fr),
            "contacts_interface": int(iface.sum()),
            "contacts_proximal": int(prox.sum()),
            "contacts_distal": int(dist.sum()),
            "interface_summary": {
                "mean_frustration": round(float(fr.frustration[iface].mean()), 4),
                "fraction_highly_frustrated": round(
                    float((fr.frustration[iface] < -1.0).mean()), 4
                ),
            },
            "whole_complex_summary": fr.summary(),
            "comparison": comparison.as_dict(),
            "stratified": strat.as_dict(),
            "mean_density_proximal": round(float(rho_pair[prox].mean()), 3)
            if prox.any() else None,
            "mean_density_distal": round(float(rho_pair[dist].mean()), 3)
            if dist.any() else None,
        }

    record["seconds"] = round(time.time() - t0, 1)
    return record


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("entries", nargs="*", help="PDB ids; default is all cached")
    ap.add_argument("--gamma", default=None,
                    help="path to AWSEM_2015.json; omit to use the stand-in")
    ap.add_argument("--decoys", type=int, default=800)
    ap.add_argument("--proximal", type=float, default=GLUE_PROXIMAL)
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()

    pot = load_potential(args.gamma, AwsemParams())
    print(f"potential: {pot.name}  ({'published' if pot.is_published() else 'STAND-IN'})")
    print(f"           {pot.provenance}\n")

    if args.entries:
        paths = [RAW / f"{e.upper()}.cif" for e in args.entries]
    else:
        paths = sorted(RAW.glob("*.cif"))

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for p in paths:
        if not p.exists():
            print(f"  {p.name}: missing")
            continue
        try:
            rec = run_one(p, pot, args.decoys, proximal=args.proximal)
        except Exception as exc:  # noqa: BLE001 - report and continue
            print(f"  {p.stem}: FAILED {type(exc).__name__}: {exc}")
            continue
        results.append(rec)
        a = rec["assignment"]
        c = rec["flavours"]["configurational"]
        cmp_ = c["comparison"]
        print(
            f"  {rec['entry']} [{a['ligand']}] "
            f"{a['anchor_name'][:22]:22s} x {a['partner_name'][:22]:22s}"
        )
        print(
            f"      interface {c['contacts_interface']:4d} contacts "
            f"({c['contacts_proximal']} proximal) | "
            f"F prox {cmp_['mean_proximal']:+.3f} vs distal "
            f"{cmp_['mean_distal']:+.3f} | diff {cmp_['difference']:+.3f} "
            f"| rb {cmp_['rank_biserial']:+.3f} | p {cmp_['p_value']} "
            f"| {rec['seconds']}s"
        )

    payload = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "potential": pot.as_dict(),
        "proximal_cutoff": args.proximal,
        "n_decoys": args.decoys,
        "arm": "oracle",
        "note": (
            "Interface geometry taken from the crystal ternary complex with "
            "the glue deleted. Upper bound on what the hypothesis can show; "
            "not a prospective prediction."
        ),
        "systems": results,
    }
    dest = out_dir / "oracle.json"
    dest.write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {dest}  ({len(results)} systems)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
