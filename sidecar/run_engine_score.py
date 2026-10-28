#!/usr/bin/env python3
"""Can the engine's scoring function recognise the right answer when handed it?

Recognising a pose and finding one are different abilities, they fail for
different reasons, and a recovery rate confuses them. This asks only the
first question, across the whole benchmark: score the crystal placement on
the engine's own grid, score a set of plausible wrong placements the same
way, and report where the crystal one lands among them.

If the crystal pose is at the top, the scoring function is sound and any
recovery failure is the search not reaching it -- which is fixed by sampling
orientations more finely, and is expensive but straightforward. If the
crystal pose is *not* at the top, no amount of searching will help and the
grid parameters are wrong. These need telling apart before spending days on
either, and on 6TD3 the answer was already clear: the crystal pose scored
+190 where forty wrong placements scored between -22 and +23.

The decoys are the same `random_surface_poses` used to calibrate the gate:
each is a random orientation brought in along a random direction until the
two proteins just touch. They are plausible complexes, not obviously-wrong
placements -- a decoy set of absurdities would make any score look
discriminating.

    python3 sidecar/run_engine_score.py
    python3 sidecar/run_engine_score.py 6TD3 --decoys 100

Needs nothing but numpy and scipy, and runs in a couple of minutes.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from glueforge import fftdock as F  # noqa: E402
from glueforge.cif import read_structure  # noqa: E402
from glueforge.face import presented_face  # noqa: E402
from glueforge.interface import assign_roles  # noqa: E402
from glueforge.pose import coords, heavy, random_surface_poses  # noqa: E402

RAW = HERE.parent / "cache" / "raw"
OUT = HERE.parent / "cache" / "recovery"


def say(msg: str = "") -> None:
    print(msg, flush=True)


def one(entry: str, n_decoys: int, spacing: float) -> dict | None:
    from scipy import fft as sfft

    path = RAW / f"{entry}.cif"
    if not path.exists():
        return None
    st = read_structure(path)
    asg = assign_roles(st)
    face = presented_face(st, asg.ligand, anchor_chains=[asg.anchor_chain])
    if not face.exposed:
        return None

    anchor = [a for r in st.residues if r.chain == asg.anchor_chain
              for a in r.atoms]
    partner = [a for r in st.residues if r.chain == asg.partner_chain
               for a in r.atoms]

    rec = coords(heavy(anchor) + heavy(asg.ligand.atoms))
    par = coords(heavy(partner))
    lig = coords(heavy(asg.ligand.atoms))
    face_xyz = lig[[a.index for a in face.exposed]]
    reach = float(np.linalg.norm(par - par.mean(axis=0), axis=1).max())

    # The same box the engine would build, so the score is the engine's.
    half = 2 * reach + F.FACE_CONTACT + F.CORE_RADIUS + F.SKIN + spacing
    lo = face_xyz.mean(axis=0) - half
    hi = face_xyz.mean(axis=0) + half
    dims = tuple(int(sfft.next_fast_len(
        int(np.ceil((h - l) / spacing)) + 1)) for l, h in zip(lo, hi))
    lat = F._Lattice(lo, dims, spacing)
    inside = rec[np.all((rec >= lo) & (rec <= hi), axis=1)]
    core = lat.occupancy(inside, F.CORE_RADIUS)
    shell = lat.occupancy(inside, F.CORE_RADIUS + F.SKIN) & ~core
    grid = np.where(core, -F.PENALTY, 0.0) + shell.astype(np.float64)

    def score(xyz: np.ndarray) -> float:
        return float(grid[lat.occupancy(xyz, F.CORE_RADIUS)].sum())

    native = score(par)
    decoys = [score(p.apply(par))
              for p in random_surface_poses(anchor, partner, n_decoys,
                                            seed=7, gap=3.0)]
    better = sum(1 for d in decoys if d >= native)
    out = {
        "entry": entry,
        "ligand": asg.ligand.comp_id,
        "partner_atoms": int(len(par)),
        "buried_fraction": round(face.buried_fraction, 4),
        "native": round(native, 1),
        "decoy_best": round(max(decoys), 1),
        "decoy_median": round(float(np.median(decoys)), 1),
        "decoys_at_or_above_native": int(better),
        "n_decoys": len(decoys),
        "recognised": bool(better == 0),
    }
    say("  %-6s %-5s %7d %6.0f%% %9.0f %9.0f %9.0f   %s" % (
        entry, asg.ligand.comp_id, out["partner_atoms"],
        100 * out["buried_fraction"], native, max(decoys),
        float(np.median(decoys)),
        "yes" if out["recognised"] else f"NO ({better}/{len(decoys)} beat it)"))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("entries", nargs="*")
    ap.add_argument("--decoys", type=int, default=60)
    ap.add_argument("--spacing", type=float, default=1.4)
    args = ap.parse_args()

    entries = [e.upper() for e in args.entries] or sorted(
        p.stem for p in RAW.glob("*.cif"))
    if not entries:
        say("cache/raw is empty. Run sidecar/fetch.py first.")
        return 1

    say(f"\n  Crystal pose against {args.decoys} plausible wrong placements, "
        f"on the engine's own grid at {args.spacing} A.\n")
    say("  %-6s %-5s %7s %7s %9s %9s %9s   %s" % (
        "entry", "lig", "atoms", "buried", "crystal", "best decoy",
        "median", "recognised"))

    rows = []
    for e in entries:
        try:
            r = one(e, args.decoys, args.spacing)
        except Exception as exc:  # noqa: BLE001
            say(f"  {e}: FAILED {type(exc).__name__}: {exc}")
            continue
        if r:
            rows.append(r)

    if not rows:
        say("\n  nothing measured")
        return 1

    ok = sum(1 for r in rows if r["recognised"])
    say("")
    say(f"    {ok}/{len(rows)} systems: the crystal pose outscores every "
        "wrong placement tried.")
    say("")
    if ok == len(rows):
        say("    So the scoring function is not what is losing poses. Every")
        say("    recovery failure is the search failing to reach a pose this")
        say("    grid would have ranked first.")
    else:
        bad = [r["entry"] for r in rows if not r["recognised"]]
        say("    The grid does not recognise the right answer for: "
            + ", ".join(bad) + ".")
        say("    For those, searching harder cannot help -- the parameters")
        say("    (CORE_RADIUS, SKIN, PENALTY) are what need changing.")

    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / "engine_score.json"
    dest.write_text(json.dumps(
        {"spacing": args.spacing, "n_decoys": args.decoys, "systems": rows},
        indent=2), encoding="utf-8")
    say(f"\n  wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
