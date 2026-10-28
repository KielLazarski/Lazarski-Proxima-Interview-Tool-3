#!/usr/bin/env python3
"""Composite cavity analysis, and how much docking error it survives.

Answers a question the frustration work could not: if the design objective is
the cavity at the glued interface rather than frustration, how accurately does
the partner have to be placed before that cavity stops being the same cavity?

The anchor-first setup makes this better conditioned than it sounds. The
ligand's pose in the anchor comes from a co-crystal, so most of the pocket is
known exactly; only the partner-lined part moves when docking is wrong. This
script measures that split.

    python3 sidecar/run_cavity.py                 # the whole cached set
    python3 sidecar/run_cavity.py 6TD3 --sensitivity
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))

from glueforge.cavity import (
    cavity_overlap,
    detect_cavity,
    field_correlation,
    rigid_perturb,
)
from glueforge.cif import read_structure
from glueforge.interface import assign_roles

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "cache" / "raw"
OUT = ROOT / "cache" / "cavity"

SHIFTS = (0.5, 1.0, 2.0, 3.0, 5.0)


def load(entry: str):
    s = read_structure(RAW / f"{entry}.cif")
    a = assign_roles(s)
    res = [r for r in s.residues if r.chain in (a.anchor_chain, a.partner_chain)]
    anchor = [at for r in res if r.chain == a.anchor_chain for at in r.atoms]
    partner = [at for r in res if r.chain == a.partner_chain for at in r.atoms]
    lig = np.array(
        [(x.x, x.y, x.z) for x in a.ligand.atoms if x.element.upper() != "H"]
    )
    return s, a, anchor, partner, lig


def describe(entry: str) -> dict:
    s, a, anchor, partner, lig = load(entry)
    centre = lig.mean(axis=0)
    cav = detect_cavity(anchor + partner, centre, set(range(len(anchor))))
    rec = (
        float((cKDTree(cav.points).query(lig, k=1)[0] < 1.2).mean())
        if len(cav.points) else 0.0
    )
    out = cav.summary()
    out.update({
        "entry": s.entry_id,
        "ligand": a.ligand.comp_id,
        "ligand_recovered": round(rec, 3),
        "anchor_lined_share": round(float((cav.from_anchor >= 0.6).mean()), 3),
        "partner_lined_share": round(float((cav.from_anchor <= 0.4).mean()), 3),
    })
    return out


def sensitivity(entry: str, n_trials: int = 3, seed: int = 7) -> dict:
    """Perturb the partner rigidly and watch the cavity degrade."""
    rng = np.random.default_rng(seed)
    _, _, anchor, partner, lig = load(entry)
    centre = lig.mean(axis=0)
    anc_idx = set(range(len(anchor)))
    ref = detect_cavity(anchor + partner, centre, anc_idx)
    if not len(ref.points):
        return {"entry": entry, "error": "no reference cavity"}

    anchor_side = ref.from_anchor >= 0.6
    partner_side = ref.from_anchor <= 0.4

    rows = []
    for mag in SHIFTS:
        acc = {k: [] for k in
               ("jaccard", "coverage", "volume_ratio", "elec", "donor",
                "anchor_keep", "partner_keep")}
        for _ in range(n_trials):
            t = rng.normal(size=3)
            t = t / np.linalg.norm(t) * mag
            moved = rigid_perturb(partner, t, mag * 2.0, int(rng.integers(1e6)))
            cav = detect_cavity(anchor + moved, centre, anc_idx)
            if not len(cav.points):
                continue
            ov = cavity_overlap(ref, cav)
            fc = field_correlation(ref, cav)
            acc["jaccard"].append(ov["jaccard"])
            acc["coverage"].append(ov["coverage"])
            acc["volume_ratio"].append(ov["volume_ratio"])
            if fc.get("electrostatic") is not None:
                acc["elec"].append(fc["electrostatic"])
            if fc.get("donor") is not None:
                acc["donor"].append(fc["donor"])
            tb = cKDTree(cav.points)
            for name, mask in (("anchor_keep", anchor_side),
                               ("partner_keep", partner_side)):
                if mask.sum() >= 10:
                    d, _ = tb.query(ref.points[mask], k=1)
                    acc[name].append(float((d <= 1.0).mean()))

        rows.append({
            "shift_angstrom": mag,
            **{k: (round(float(np.mean(v)), 3) if v else None)
               for k, v in acc.items()},
        })

    return {
        "entry": entry,
        "reference": ref.summary(),
        "anchor_lined_share": round(float(anchor_side.mean()), 3),
        "partner_lined_share": round(float(partner_side.mean()), 3),
        "rows": rows,
        "note": (
            "The partner is displaced rigidly by the given magnitude with a "
            "matched rotation; the anchor and the ligand stay put. Jaccard is "
            "shape agreement with the unperturbed cavity, elec and donor are "
            "field correlations over the overlap, and anchor_keep / "
            "partner_keep are the share of each region that survives."
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("entries", nargs="*")
    ap.add_argument("--sensitivity", action="store_true")
    ap.add_argument("--trials", type=int, default=3)
    args = ap.parse_args()

    entries = args.entries or sorted(p.stem for p in RAW.glob("*.cif"))
    OUT.mkdir(parents=True, exist_ok=True)

    described, sens = [], []
    print(f"{'entry':6s} {'lig':>5s} {'vol A3':>7s} {'encl':>5s} {'comp':>5s} "
          f"{'anchor%':>8s} {'elec':>6s} {'recov':>6s}")
    print("-" * 56)
    for e in entries:
        if not (RAW / f"{e}.cif").exists():
            continue
        try:
            d = describe(e)
        except Exception as exc:  # noqa: BLE001
            print(f"{e}: FAILED {type(exc).__name__}: {exc}")
            continue
        described.append(d)
        print(f"{d['entry']:6s} {d['ligand']:>5s} {d.get('volume', 0):7.0f} "
              f"{d.get('enclosure_quality', 0):5.2f} "
              f"{d.get('composite_fraction', 0):5.2f} "
              f"{d['anchor_lined_share']:8.2f} "
              f"{d.get('mean_electrostatic', 0):+6.2f} "
              f"{100 * d['ligand_recovered']:5.0f}%")

        if args.sensitivity:
            t0 = time.time()
            s = sensitivity(e, args.trials)
            sens.append(s)
            print(f"       sensitivity: ", end="")
            for r in s["rows"]:
                print(f"{r['shift_angstrom']}A J={r['jaccard']} ", end="")
            print(f" ({time.time() - t0:.0f}s)")

    payload = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "systems": described,
        "sensitivity": sens,
    }
    (OUT / "cavity.json").write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {OUT / 'cavity.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
