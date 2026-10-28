#!/usr/bin/env python3
"""How often does the built-in docking engine find the right answer?

The same question `run_recovery.py` asks of LightDock, asked of the engine
that ships with GlueForge, and asked the same way so the two numbers mean the
same thing: take a solved ternary complex, throw the partner away, measure
the ligand's presented face against the anchor alone, aim the search at it,
and ask how far the best returned pose is from where the partner really was.

Two things about what this measures, both of which make the number here
**better than real use**, and both of which have to be said next to it.

**It is bound-bound docking.** Both halves come from the same crystal, so
both are already in the conformation they adopt on binding. Real use starts
from two separately solved structures whose side chains, and sometimes
backbones, are in the wrong place. Published docking benchmarks separate
these cases for good reason: bound-bound is substantially the easier one,
and a rigid-body method has nothing to say about the difference.

**Fourteen systems is a small benchmark**, and they were chosen as systems
with a ternary complex to measure a presented face on, not as a
representative sample of docking difficulty.

    python3 sidecar/run_builtin_recovery.py                  # quick, all
    python3 sidecar/run_builtin_recovery.py --preset standard
    python3 sidecar/run_builtin_recovery.py 6TD3 8BU7        # named entries

Needs nothing but numpy and scipy, which is the entire point of it.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from glueforge.cif import read_structure  # noqa: E402
from glueforge.face import presented_face  # noqa: E402
from glueforge.fftdock import PRESETS, dock  # noqa: E402
from glueforge.interface import assign_roles  # noqa: E402
from glueforge.pose import PoseScorer, coords, heavy  # noqa: E402

RAW = HERE.parent / "cache" / "raw"
OUT = HERE.parent / "cache" / "recovery"

#: A hit for docking purposes.
HIT = 5.0
#: Close enough for the cavity at the interface to still be the same cavity.
DESIGNABLE = 2.0


def say(msg: str = "") -> None:
    print(msg, flush=True)


def one(entry: str, settings, top: int) -> dict | None:
    path = RAW / f"{entry}.cif"
    if not path.exists():
        say(f"  {entry}: not in cache/raw, skipping")
        return None

    t0 = time.time()
    st = read_structure(path)
    asg = assign_roles(st)
    face = presented_face(st, asg.ligand, anchor_chains=[asg.anchor_chain])
    if not face.exposed:
        say(f"  {entry}: no exposed face, nothing to aim at")
        return None

    anchor_atoms = [a for r in st.residues if r.chain == asg.anchor_chain
                    for a in r.atoms]
    partner_atoms = [a for r in st.residues if r.chain == asg.partner_chain
                     for a in r.atoms]
    native_xyz = coords(heavy(partner_atoms))

    say(f"  {entry} {asg.ligand.comp_id:>4}  "
        f"face {face.presented_area:6.1f} A^2, "
        f"{100*face.buried_fraction:3.0f}% buried, "
        f"partner {len(native_xyz)} atoms ... ")

    poses = dock(anchor_atoms, asg.ligand.atoms, partner_atoms, face, settings)
    if not poses:
        say("    the search returned nothing")
        return {"entry": entry, "preset": settings.key,
                "error": "no poses returned"}

    scorer = PoseScorer(anchor_atoms, asg.ligand.atoms, face)
    rows = []
    for p in (poses[:top] if top else poses):
        s = scorer.score(p.apply_atoms(partner_atoms), label=p.label,
                         reference_xyz=native_xyz)
        rows.append(s.as_dict())

    viable = [r for r in rows if r["viable"]]
    pool = viable or rows
    best_rmsd = min(r["rmsd_to_reference"] for r in pool)
    by_face = max(pool, key=lambda r: r["face_burial"])
    # `rows` is in the engine's own order, which after refinement is by the
    # continuous score it optimised. Recorded separately from the face-burial
    # pick because the two disagree, and which of them to trust is the open
    # question once the search is finding the right answer at all.
    by_engine = pool[0]

    out = {
        "entry": entry,
        "ligand": asg.ligand.comp_id,
        "preset": settings.key,
        "engine": "builtin",
        "n_poses": len(rows),
        "n_viable": len(viable),
        "partner_atoms": int(len(native_xyz)),
        "presented_area": round(face.presented_area, 1),
        "buried_fraction": round(face.buried_fraction, 4),
        "native_face_burial": round(
            scorer.score(partner_atoms, "crystal").face_burial, 4),
        "best_rmsd": round(best_rmsd, 2),
        "best_by_face_burial": by_face,
        # Every scored pose, in the engine's own order. Kept because
        # choosing one pose from the output is a separate question from
        # finding it, and compare_rankings.py can only replay that choice
        # over poses it can see -- without this it had the answer but not
        # the working, and reported that there was nothing scored.
        "poses": rows,
        "engine_top_rmsd": round(by_engine["rmsd_to_reference"], 2),
        "recovered": bool(best_rmsd <= HIT),
        "designable": bool(best_rmsd <= DESIGNABLE),
        "top_pose_rmsd": round(by_face["rmsd_to_reference"], 2),
        "seconds": round(time.time() - t0, 1),
    }
    say(f"    best {best_rmsd:6.2f} A   "
        f"by face {by_face['rmsd_to_reference']:6.2f} A   "
        f"by score {by_engine['rmsd_to_reference']:6.2f} A   "
        f"{'HIT' if out['recovered'] else 'miss'}"
        f"{'  (designable)' if out['designable'] else ''}   "
        f"[{out['seconds']:.0f}s, {len(rows)} poses]")
    return out


def report(rows: list[dict], settings) -> None:
    good = [r for r in rows if "error" not in r]
    if not good:
        say("\n  every system failed")
        return
    rec = sum(1 for r in good if r["recovered"])
    des = sum(1 for r in good if r["designable"])
    rms = [r["best_rmsd"] for r in good]

    say("")
    say("  " + "=" * 68)
    say(f"  built-in engine, {settings.label}: {settings.rotations} "
        f"orientations at {settings.spacing} A")
    say("  " + "=" * 68)
    say(f"    recovered (<={HIT:.0f} A)     {rec}/{len(good)}"
        f"   = {100*rec/len(good):.0f}%")
    say(f"    designable (<={DESIGNABLE:.0f} A)    {des}/{len(good)}"
        f"   = {100*des/len(good):.0f}%")
    say(f"    median best RMSD       {np.median(rms):.2f} A")

    say("")
    eng = [r for r in good if r.get("engine_top_rmsd") is not None]
    if eng:
        hit_e = sum(1 for r in eng if r["engine_top_rmsd"] <= HIT)
        hit_f = sum(1 for r in eng if r["top_pose_rmsd"] <= HIT)
        say(f"    top pose right (score)  {hit_e}/{len(eng)}")
        say(f"    top pose right (face)   {hit_f}/{len(eng)}")

    say("")
    say("    %-6s %-5s %7s %8s %9s %9s %9s %7s" % (
        "entry", "lig", "buried", "atoms", "best", "by face", "by score",
        "secs"))
    for r in rows:
        if "error" in r:
            say("    %-6s %-5s   %s" % (r["entry"], "", r["error"]))
            continue
        say("    %-6s %-5s %6.0f%% %8d %9.2f %9.2f %9.2f %7.0f%s" % (
            r["entry"], r["ligand"], 100 * r["buried_fraction"],
            r["partner_atoms"], r["best_rmsd"], r["top_pose_rmsd"],
            r.get("engine_top_rmsd", float("nan")),
            r["seconds"], "  HIT" if r["recovered"] else ""))

    say("")
    say("    'best' is the closest pose anywhere in the output, which is what")
    say("    the search is worth. 'by face' and 'by score' are the two ways of")
    say("    choosing one from it -- face burial, and the engine's own")
    say("    continuous score. The gap between 'best' and either of those is")
    say("    the ranking's failure, not the search's, and they are separate")
    say("    problems with separate fixes.")
    say("")
    say("    Both halves of every pair come from the same crystal, so both")
    say("    are already in their bound conformation. Real use is harder.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("entries", nargs="*")
    ap.add_argument("--preset", default="quick", choices=sorted(PRESETS))
    # Everything the engine returned, not a prefix of it. The engine now
    # returns each pose and its refinement, ordered by a score -- so a
    # prefix drops the lower-scoring half, which is exactly the half the
    # refinement stage was just fixed for discarding. Cutting here would
    # have reinstated the bug one level up and measured the fix as a
    # failure.
    ap.add_argument("--top", type=int, default=0,
                    help="score only the first N poses; 0 means all")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    entries = [e.upper() for e in args.entries] or sorted(
        p.stem for p in RAW.glob("*.cif"))
    if not entries:
        say("cache/raw is empty. Run sidecar/fetch.py first.")
        return 1

    settings = PRESETS[args.preset]
    lo, hi = settings.minutes_estimate
    say(f"\nbuilt-in engine, {settings.label}: {len(entries)} systems, "
        f"up to {lo}-{hi} minutes each\n")
    say("  The estimate is for the largest partner in the benchmark. The box "
        "is\n  sized by the partner's reach and the cost goes as its cube, so "
        "a small\n  partner is far quicker than that and a large one is not.\n")

    rows = []
    for e in entries:
        try:
            r = one(e, settings, args.top)
        except Exception as exc:  # noqa: BLE001
            say(f"  {e}: FAILED {type(exc).__name__}: {exc}")
            rows.append({"entry": e, "preset": settings.key,
                         "error": f"{type(exc).__name__}: {exc}"})
            continue
        if r:
            rows.append(r)

    report(rows, settings)

    OUT.mkdir(parents=True, exist_ok=True)
    dest = Path(args.out) if args.out else OUT / "builtin.json"
    dest.write_text(json.dumps({
        "machine": {"platform": platform.platform(),
                    "python": platform.python_version()},
        "engine": "builtin",
        "bound_bound": True,
        "cutoff": HIT,
        "designable_cutoff": DESIGNABLE,
        "settings": settings.as_dict(),
        "systems": rows,
    }, indent=2), encoding="utf-8")
    say(f"\nwrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
