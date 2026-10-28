#!/usr/bin/env python3
"""Does the presented-face restraint carry any signal?

This is the question the whole tool rests on, and it can be answered without
a docking program. The claim is that requiring the partner to bury the
ligand's exposed face narrows the search to the right part of the surface. If
that claim is false, no amount of docking quality rescues it: a restraint that
does not separate the native pose from plausible wrong ones is a restraint
that does nothing except make you think you have constrained something.

The test. For each system, take the crystal complex, keep the anchor with its
ligand, and throw the partner away. Then put the partner back in several
hundred random orientations, each pushed in along a random direction until the
two just touch. Every decoy is therefore a *plausible* complex -- in contact,
not interpenetrating, not floating in space -- which is the only kind worth
comparing against. Score each on how much of the presented face it buries, and
ask where the crystal pose lands.

Two numbers come out, and they say different things:

  **rank**        where the native pose sits among the decoys. Rank 1 means
                  no random placement buried the face better than the real
                  one.
  **enrichment**  how much the restraint narrows the search. If 2% of random
                  plausible poses clear the burial threshold, then insisting
                  on it is a 50-fold reduction in what a docking run has to
                  search.

Enrichment is the one that matters operationally. A native pose can rank first
while half the decoys are nearly as good, and that restraint has told you
almost nothing.

    python3 sidecar/run_decoys.py                 # all 14, 400 decoys each
    python3 sidecar/run_decoys.py 6TD3 8BU7       # named entries
    python3 sidecar/run_decoys.py --decoys 2000   # tighter tails
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
from glueforge.interface import assign_roles  # noqa: E402
from glueforge.pose import (  # noqa: E402
    Pose,
    PoseScorer,
    coords,
    heavy,
    random_surface_poses,
)

RAW = HERE.parent / "cache" / "raw"
OUT = HERE.parent / "cache" / "decoys"

#: Burial thresholds to report enrichment at. 0.5 is the one the gate would
#: plausibly use; the others are there so the threshold's arbitrariness is
#: visible rather than hidden in a single number.
THRESHOLDS = (0.25, 0.50, 0.75)


def analyse(entry: str, n_decoys: int, seed: int,
            gap: float = 2.6) -> dict | None:
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
        print(f"  {entry}: no exposed face, nothing to restrain on")
        return None

    scorer = PoseScorer(anchor, asg.ligand.atoms, face)

    # How close the two proteins actually come in this crystal. It bounds
    # how close a decoy may fairly be pushed: if the native may press in to
    # a hydrogen bond while decoys are held at van der Waals contact, the
    # native gets face burial it has not earned and the whole comparison is
    # tilted. Reported so the choice of `gap` can be judged rather than
    # trusted.
    from scipy.spatial import cKDTree
    native_xyz = coords(heavy(partner))
    native_nearest = float(
        cKDTree(coords(heavy(anchor))).query(native_xyz, k=1)[0].min()
    )
    native = scorer.score(partner, label="crystal",
                          reference_xyz=native_xyz)

    decoys = random_surface_poses(anchor, partner, n_decoys, seed=seed,
                                  gap=gap)
    scores = []
    for p in decoys:
        moved = p.apply_atoms(heavy(partner))
        scores.append(scorer.score(moved, label=p.label,
                                   reference_xyz=native_xyz))

    # Only plausible decoys count. One that drives the partner through the
    # anchor can bury the whole face and is not a competing hypothesis.
    viable = [s for s in scores if s.viable]
    burials = np.array([s.face_burial for s in viable]) if viable else np.zeros(0)

    better = int((burials > native.face_burial).sum())
    rank = better + 1
    pct = float((burials <= native.face_burial).mean() * 100) if len(burials) else float("nan")

    enrich = {}
    for t in THRESHOLDS:
        frac = float((burials >= t).mean()) if len(burials) else float("nan")
        enrich[str(t)] = {
            "decoy_fraction": round(frac, 5),
            "fold": (round(1.0 / frac, 1) if frac > 0 else None),
            "native_clears": bool(native.face_burial >= t),
        }

    # Does burial track being in the right place at all? If the correlation
    # with RMSD is flat, high burial is not evidence of a near-native pose.
    rms = np.array([s.rmsd_to_reference for s in viable]) if viable else np.zeros(0)
    if len(rms) > 8:
        from scipy.stats import spearmanr
        rho, pval = spearmanr(burials, rms)
    else:
        rho, pval = float("nan"), float("nan")

    out = {
        "entry": entry,
        "ligand": asg.ligand.comp_id,
        "anchor_chain": asg.anchor_chain,
        "partner_chain": asg.partner_chain,
        "presented_area": round(face.presented_area, 1),
        "buried_fraction": round(face.buried_fraction, 4),
        "n_exposed": len(face.exposed),
        "n_rim": len(face.rim),
        "n_decoys": n_decoys,
        "decoy_gap": gap,
        "n_viable": len(viable),
        "native_nearest_approach": round(native_nearest, 2),
        "native": native.as_dict(),
        "decoy_burial": {
            "mean": round(float(burials.mean()), 4) if len(burials) else None,
            "median": round(float(np.median(burials)), 4) if len(burials) else None,
            "p90": round(float(np.percentile(burials, 90)), 4) if len(burials) else None,
            "max": round(float(burials.max()), 4) if len(burials) else None,
        },
        "native_rank": rank,
        "native_percentile": round(pct, 1) if pct == pct else None,
        "decoys_better": better,
        "enrichment": enrich,
        "burial_vs_rmsd": {
            "spearman_rho": round(float(rho), 4) if rho == rho else None,
            "p": float(pval) if pval == pval else None,
        },
        "seconds": round(time.time() - t0, 1),
    }

    print(
        f"  {entry} {asg.ligand.comp_id:>4}  native {native.face_burial:5.1%}  "
        f"decoys {len(viable):4d} viable, median {np.median(burials):5.1%} "
        f"max {burials.max():5.1%}" if len(burials) else
        f"  {entry}: no viable decoys"
    )
    print(
        f"           rank {rank}/{len(viable) + 1}  "
        f"({better} better)   "
        f"enrichment at 50%: "
        + (f"{enrich['0.5']['fold']}x" if enrich["0.5"]["fold"] else "no decoy clears it")
        + f"   rho(burial, rmsd) = "
        + (f"{rho:+.2f}" if rho == rho else "n/a")
        + f"   [{out['seconds']}s]"
    )
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("entries", nargs="*", help="entries, default all cached")
    ap.add_argument("--decoys", type=int, default=400)
    ap.add_argument("--seed", type=int, default=20261027)
    ap.add_argument(
        "--gap", type=float, default=2.6,
        help=(
            "closest a decoy's heavy atoms may come to the anchor's, in "
            "angstroms. The default matches how close the real partners come "
            "in these crystals (2.41-2.75 A, median 2.66): holding decoys "
            "further out than the native sits would hand the native burial "
            "it has not earned."
        ),
    )
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    entries = [e.upper() for e in args.entries] or sorted(
        p.stem for p in RAW.glob("*.cif")
    )
    if not entries:
        print("cache/raw is empty. Run sidecar/fetch.py first.")
        return 1

    print(f"\nRestraint informativeness: {len(entries)} systems, "
          f"{args.decoys} random surface poses each")
    print(f"decoys held at {args.gap} A nearest heavy-atom approach "
          f"(the real partners sit at 2.41-2.75)\n")

    rows = []
    for e in entries:
        try:
            r = analyse(e, args.decoys, args.seed, args.gap)
        except Exception as exc:  # noqa: BLE001 - one bad entry must not stop the run
            print(f"  {e}: FAILED {type(exc).__name__}: {exc}")
            continue
        if r:
            rows.append(r)

    if not rows:
        print("\nnothing analysed")
        return 1

    print("\n" + "=" * 74)
    print("SUMMARY")
    print("=" * 74)

    ranks = [r["native_rank"] for r in rows]
    first = sum(1 for r in ranks if r == 1)
    top5 = sum(1 for r in ranks if r <= 5)
    folds = [r["enrichment"]["0.5"]["fold"] for r in rows
             if r["enrichment"]["0.5"]["fold"]]
    nat50 = sum(1 for r in rows if r["enrichment"]["0.5"]["native_clears"])
    rhos = [r["burial_vs_rmsd"]["spearman_rho"] for r in rows
            if r["burial_vs_rmsd"]["spearman_rho"] is not None]

    print(f"  native pose ranks first          {first}/{len(rows)}")
    print(f"  native pose in the top five      {top5}/{len(rows)}")
    print(f"  native buries >=50% of the face  {nat50}/{len(rows)}")
    # Most systems have no decoy clearing the threshold at all, so a median
    # over the few that do is a median over a biased handful. Report the
    # count first and the fold only for the systems where it is finite.
    clean = len(rows) - len(folds)
    print(f"  no decoy cleared 50% burial      {clean}/{len(rows)} systems"
          f"  (enrichment there is > the viable decoy count, ~200x)")
    if folds:
        print(f"  where some decoy did clear it    {len(folds)} system(s), "
              f"enrichment " + ", ".join(f"{f:.0f}x" for f in sorted(folds)))
    if rhos:
        print(f"  rho(burial, RMSD)                median {np.median(rhos):+.2f}"
              f"  (negative = burial rises as the pose gets closer)")

    # The operationally useful question: can you tell, from the upload screen
    # alone, whether the restraint is going to be worth anything? Buried
    # fraction costs nothing to compute and is known before any docking.
    bf = np.array([r["buried_fraction"] for r in rows])
    rk = np.array([r["native_rank"] for r in rows], dtype=float)
    if len(rows) >= 6:
        from scipy.stats import spearmanr
        rho_rank, p_rank = spearmanr(bf, rk)
        print()
        print("  Can the restraint's usefulness be predicted before docking?")
        print(f"    rho(buried fraction, native rank) = {rho_rank:+.2f}"
              f"   p = {p_rank:.4f}   n = {len(rows)}")
        print("    negative means a more deeply buried ligand gives a better rank")

        deep = [r for r in rows if r["buried_fraction"] >= 0.55]
        shallow = [r for r in rows if r["buried_fraction"] < 0.55]
        for name, grp in (("pocket  (>=55% buried)", deep),
                          ("groove  ( <55% buried)", shallow)):
            if not grp:
                continue
            first = sum(1 for r in grp if r["native_rank"] == 1)
            med = float(np.median([r["native_rank"] for r in grp]))
            folds = [r["enrichment"]["0.5"]["fold"] for r in grp
                     if r["enrichment"]["0.5"]["fold"]]
            print(f"    {name}  n={len(grp):2d}   rank 1 in {first}/{len(grp)}"
                  f"   median rank {med:.0f}"
                  f"   no decoy cleared 50% in {len(grp)-len(folds)}/{len(grp)}")

    print()
    print("  %-6s %-5s %7s %8s %8s %8s %8s %9s" % (
        "entry", "lig", "buried", "native", "median", "best", "rank",
        "enrich50"))
    for r in sorted(rows, key=lambda x: -x["buried_fraction"]):
        d = r["decoy_burial"]
        f = r["enrichment"]["0.5"]["fold"]
        print("  %-6s %-5s %6.0f%% %7.1f%% %7.1f%% %7.1f%% %4d/%-3d %9s" % (
            r["entry"], r["ligand"],
            100 * r["buried_fraction"],
            100 * r["native"]["face_burial"],
            100 * (d["median"] or 0), 100 * (d["max"] or 0),
            r["native_rank"], r["n_viable"] + 1,
            (f"{f:.0f}x" if f else ">%dx" % r["n_viable"])))

    dest = Path(args.out or (OUT / f"decoys-gap{args.gap:.1f}.json"))
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps({
        "note": (
            "Whether the presented-face restraint separates the crystal pose "
            "from plausible random placements. No docking program involved: "
            "decoys are uniform random orientations of the partner, each "
            "pushed in along a random direction until contact."
        ),
        "n_decoys": args.decoys,
        "seed": args.seed,
        "gap": args.gap,
        "thresholds": list(THRESHOLDS),
        "systems": rows,
    }, indent=2), encoding="utf-8")
    print(f"\nwrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
