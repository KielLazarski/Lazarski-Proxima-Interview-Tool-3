#!/usr/bin/env python3
"""Calibrate the gate against real negatives, not just against known glues.

The assessment screen has to answer one question: does this interface have a
pocket worth designing chemistry into? Calibrating that on the 14 solved glue
systems alone gives thresholds that say "this resembles a known glue" and
cannot say "this is not one", because every example is a positive. A gate
built that way passes whatever it is shown and feels informative while being
vacuous.

Negatives are available without needing any new structures. A decoy pose of
the partner -- a plausible complex, in contact, not interpenetrating, in the
wrong place -- produces a real interface with a real cavity at the ligand
site. It is exactly the thing the gate has to reject. The hardest negatives
are the decoys that bury the presented face best, so those are the ones used.

What comes out: for each metric, how well it separates the crystal pose from
its own hardest decoys, as the area under the ROC curve. 0.5 is no signal at
all; 1.0 is perfect separation. A metric below about 0.7 has no business
being in a verdict, however sensible it sounds.

    python3 sidecar/run_gate.py                   # all 14, 12 decoys each
    python3 sidecar/run_gate.py 6TD3 --decoys 30
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from glueforge.cif import read_structure  # noqa: E402
from glueforge.face import presented_face  # noqa: E402
from glueforge.gate import cavity_metrics  # noqa: E402
from glueforge.interface import assign_roles  # noqa: E402
from glueforge.pose import (  # noqa: E402
    PoseScorer,
    coords,
    heavy,
    random_surface_poses,
)

RAW = HERE.parent / "cache" / "raw"
OUT = HERE.parent / "cache" / "gate"

#: Metrics the gate could use, and which direction is good.
METRICS = (
    ("volume", +1, "cavity volume", "Å³"),
    ("enclosure_quality", +1, "enclosure quality", ""),
    ("composite_fraction", +1, "lined by both proteins", ""),
    ("mean_enclosure", +1, "mean enclosure", ""),
    ("ligand_recovery", +1, "cavity covers the ligand", ""),
    ("face_burial", +1, "presented face buried", ""),
)


def youden(pos: np.ndarray, neg: np.ndarray, sign: int) -> dict:
    """The threshold that separates best, and how well it does.

    Using an endpoint of the positive distribution instead -- the weakest
    validated system, which is the obvious choice -- lets one unusual system
    drag the bar down for everything, and the hardest decoys then clear it.
    The Youden point maximises (sensitivity + specificity - 1), so it is set
    by where the two distributions actually part company rather than by how
    far one of them happens to reach.
    """
    cand = np.unique(np.concatenate([pos, neg]))
    best = None
    for t in cand:
        tpr = float((pos * sign >= t * sign).mean())
        fpr = float((neg * sign >= t * sign).mean())
        j = tpr - fpr
        if best is None or j > best[0]:
            best = (j, float(t), tpr, fpr)
    j, t, tpr, fpr = best
    return {"threshold": round(t, 4), "youden_j": round(j, 4),
            "sensitivity": round(tpr, 4), "specificity": round(1 - fpr, 4)}


def auc(pos: np.ndarray, neg: np.ndarray) -> float:
    """Area under the ROC, by rank. Ties count a half."""
    if not len(pos) or not len(neg):
        return float("nan")
    wins = 0.0
    for p in pos:
        wins += float((neg < p).sum()) + 0.5 * float((neg == p).sum())
    return wins / (len(pos) * len(neg))


def analyse(entry: str, n_decoys: int, n_keep: int, seed: int) -> dict | None:
    path = RAW / f"{entry}.cif"
    if not path.exists():
        print(f"  {entry}: not in cache/raw, skipping")
        return None

    t0 = time.time()
    st = read_structure(path)
    asg = assign_roles(st)
    anchor = [a for r in st.residues if r.chain == asg.anchor_chain
              for a in r.atoms]
    partner = [a for r in st.residues if r.chain == asg.partner_chain
               for a in r.atoms]
    face = presented_face(st, asg.ligand, anchor_chains=[asg.anchor_chain])
    if not face.exposed:
        print(f"  {entry}: no exposed face")
        return None

    lig_xyz = coords(heavy(asg.ligand.atoms))
    scorer = PoseScorer(anchor, asg.ligand.atoms, face)

    native = cavity_metrics(heavy(anchor), heavy(partner), lig_xyz)
    native["face_burial"] = scorer.score(partner, "crystal").face_burial

    # Generate plenty, keep the ones that bury the face best. Those are the
    # decoys the gate would actually be asked about, because the docking run
    # would have ranked them highest.
    poses = random_surface_poses(anchor, partner, n_decoys, seed=seed, gap=2.6)
    scored = []
    p_heavy = heavy(partner)
    for p in poses:
        moved = p.apply_atoms(p_heavy)
        s = scorer.score(moved, p.label)
        if s.viable:
            scored.append((s.face_burial, moved, s))
    scored.sort(key=lambda t: -t[0])
    hardest = scored[:n_keep]

    decoys = []
    for burial, moved, s in hardest:
        m = cavity_metrics(heavy(anchor), moved, lig_xyz)
        m["face_burial"] = burial
        decoys.append(m)

    row = {
        "entry": entry,
        "ligand": asg.ligand.comp_id,
        "buried_fraction": round(face.buried_fraction, 4),
        "directionality": round(face.directionality, 3),
        "native": {k: round(v, 4) for k, v in native.items()},
        "decoys": [{k: round(v, 4) for k, v in d.items()} for d in decoys],
        "n_decoys_scored": len(scored),
        "n_decoys_kept": len(decoys),
        "seconds": round(time.time() - t0, 1),
    }

    beat = {}
    for key, sign, _, _ in METRICS:
        vals = np.array([d[key] for d in decoys])
        n = native[key]
        beat[key] = int((vals * sign >= n * sign).sum())
    row["decoys_matching_or_beating_native"] = beat

    print("  %-6s %-5s  native vol %6.0f encl %.2f comp %.2f rec %.2f  "
          "| %2d decoys, beaten by %s  [%.0fs]" % (
              entry, asg.ligand.comp_id, native["volume"],
              native["enclosure_quality"], native["composite_fraction"],
              native["ligand_recovery"], len(decoys),
              "/".join(str(beat[k]) for k, _, _, _ in METRICS),
              row["seconds"]))
    return row


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("entries", nargs="*")
    ap.add_argument("--decoys", type=int, default=120,
                    help="random poses generated per system")
    ap.add_argument("--keep", type=int, default=12,
                    help="hardest decoys kept and given a cavity")
    ap.add_argument("--seed", type=int, default=20261027)
    args = ap.parse_args()

    entries = [e.upper() for e in args.entries] or sorted(
        p.stem for p in RAW.glob("*.cif"))
    if not entries:
        print("cache/raw is empty. Run sidecar/fetch.py first.")
        return 1

    print(f"\nGate calibration: {len(entries)} systems, {args.decoys} poses "
          f"each, cavity for the {args.keep} that bury the face best")
    print("beaten-by counts are in metric order: " +
          ", ".join(k for k, _, _, _ in METRICS) + "\n")

    rows = []
    for e in entries:
        try:
            r = analyse(e, args.decoys, args.keep, args.seed)
        except Exception as exc:  # noqa: BLE001
            print(f"  {e}: FAILED {type(exc).__name__}: {exc}")
            continue
        if r:
            rows.append(r)
    if not rows:
        return 1

    print("\n" + "=" * 72)
    print("HOW WELL DOES EACH METRIC SEPARATE THE RIGHT POSE FROM THE WRONG ONES?")
    print("=" * 72)
    print("  AUC over all systems pooled. 0.5 is no signal; below about 0.7 a")
    print("  metric should not be in a verdict however sensible it sounds.\n")

    pooled = {}
    print("  %-26s %7s %9s %11s %10s %7s" % (
        "metric", "AUC", "native", "decoy p90", "threshold", "spec"))
    for key, sign, label, unit in METRICS:
        pos = np.array([r["native"][key] for r in rows]) * sign
        neg = np.array([d[key] for r in rows for d in r["decoys"]]) * sign
        a = auc(pos, neg)
        pv, nv = pos * sign, neg * sign
        pooled[key] = {
            "auc": round(float(a), 4),
            "native_median": round(float(np.median(pv)), 4),
            "native_min": round(float(pv.min()), 4),
            "native_p10": round(float(np.percentile(pv, 10 if sign > 0 else 90)), 4),
            "decoy_p90": round(float(np.percentile(nv, 90 if sign > 0 else 10)), 4),
            "decoy_median": round(float(np.median(nv)), 4),
            "decoy_max": round(float(nv.max() if sign > 0 else nv.min()), 4),
            "direction": sign,
            "label": label,
            "unit": unit,
            **youden(pos, neg, sign),
        }
        print("  %-26s %7.3f %9.3f %11.3f %10.3f %6.0f%% %s" % (
            label, a, np.median(pos * sign), pooled[key]["decoy_p90"],
            pooled[key]["threshold"], 100 * pooled[key]["specificity"],
            "" if a >= 0.7 else "  <- too weak to use"))

    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / "gate.json"
    dest.write_text(json.dumps({
        "note": (
            "Gate calibration against decoy poses. Positives are the crystal "
            "placement of the partner; negatives are the plausible random "
            "placements that bury the presented face best, which are the ones "
            "a docking run would rank highest and the gate would actually be "
            "asked about. AUC 0.5 means the metric carries no signal."
        ),
        "n_systems": len(rows),
        "decoys_per_system": args.keep,
        "seed": args.seed,
        "metrics": pooled,
        "systems": rows,
    }, indent=2), encoding="utf-8")
    print(f"\nwrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
