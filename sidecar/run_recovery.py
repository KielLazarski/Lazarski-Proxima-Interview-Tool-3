#!/usr/bin/env python3
"""How often does a restrained docking run find the right answer?

This is the number the three presets are supposed to carry, and it is the one
thing in GlueForge that cannot be measured without LightDock. Run it on a
machine that has LightDock and send the output back; until then the presets
show their run times and admit they do not know their recovery rates, which
is better than showing a made-up one.

What it does, for each benchmark system: take the solved ternary complex,
throw the partner away, measure the ligand's presented face against the
anchor alone, restrain a docking run to the rim of that face, and then ask
how far the best returned pose is from where the partner really was.

"Recovered" means the best pose is within `--cutoff` angstroms of the crystal
placement, measured without superposition -- superposing would remove exactly
the error being measured. The default of 5 A is not arbitrary: the cavity at
the interface tolerates about 2 A of partner error and falls apart by 3, so
a 5 A recovery is a hit for docking purposes and still too coarse to design
chemistry into. Both numbers are reported.

    python3 sidecar/run_recovery.py                     # standard preset, all 14
    python3 sidecar/run_recovery.py --preset quick      # the fast one
    python3 sidecar/run_recovery.py 6TD3 8BU7           # named entries
    python3 sidecar/run_recovery.py --all-presets       # all three, one after another

Expect this to take hours. It is a benchmark, not a step in the workflow.
"""

from __future__ import annotations

import argparse
import json
import platform
import shutil
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from glueforge.cif import read_any, read_structure  # noqa: E402
from glueforge.face import presented_face  # noqa: E402
from glueforge.interface import assign_roles  # noqa: E402
from glueforge.pose import PoseScorer, coords, heavy, rmsd  # noqa: E402
from glueforge.restraints import (  # noqa: E402
    PRESETS,
    cluster_command,
    conformations_command,
    count_swarms,
    missing_tools,
    rank_command,
    restraints_from_face,
    run,
    run_command,
    setup_command,
    split_for_docking,
)

RAW = HERE.parent / "cache" / "raw"
OUT = HERE.parent / "cache" / "recovery"

#: A hit for docking purposes.
HIT = 5.0
#: Close enough for the cavity at the interface to still be the same cavity.
DESIGNABLE = 2.0


def say(msg: str = "") -> None:
    print(msg, flush=True)


def one(entry: str, preset, work_root: Path, cutoff: float,
        cores: int, top: int) -> dict | None:
    path = RAW / f"{entry}.cif"
    if not path.exists():
        say(f"  {entry}: not in cache/raw, skipping")
        return None

    t0 = time.time()
    work = work_root / f"{entry}-{preset.key}"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    st = read_structure(path)
    asg = assign_roles(st)
    face = presented_face(st, asg.ligand, anchor_chains=[asg.anchor_chain])
    if not face.exposed:
        say(f"  {entry}: no exposed face, nothing to restrain on")
        return None

    rst = restraints_from_face(face)
    files = split_for_docking(
        st, [asg.anchor_chain], [asg.partner_chain], work, ligand=asg.ligand)
    (work / "restraints.list").write_text(rst.text(), encoding="ascii")

    rec, par = "receptor.pdb", "partner.pdb"
    log_path = work / "lightdock.log"
    log_file = log_path.open("w", encoding="utf-8")

    def log(line: str) -> None:
        log_file.write(line + "\n")
        log_file.flush()

    say(f"  {entry} {asg.ligand.comp_id:>4}  "
        f"face {face.presented_area:6.1f} A^2, {100*face.buried_fraction:3.0f}% buried, "
        f"{len(rst.active)}+{len(rst.passive)} restraints ... ", )

    code = run(setup_command(Path(rec), Path(par),
                             Path("restraints.list"), preset), work, log)
    if code != 0:
        log_file.close()
        say(f"    setup failed (exit {code}); see {log_path}")
        return {"entry": entry, "preset": preset.key, "error": f"setup exit {code}"}

    n_swarms = count_swarms(work)
    if n_swarms == 0:
        log_file.close()
        say("    no swarm satisfied the restraints")
        return {"entry": entry, "preset": preset.key, "error": "no swarms",
                "restraints": rst.as_dict()}

    code = run(run_command(preset, cores), work, log)
    if code != 0:
        log_file.close()
        say(f"    docking failed (exit {code}); see {log_path}")
        return {"entry": entry, "preset": preset.key, "error": f"dock exit {code}"}

    for i in range(n_swarms):
        d = work / f"swarm_{i}"
        if d.is_dir():
            run(conformations_command(Path("../" + rec), Path("../" + par),
                                      preset), d, None)
            run(cluster_command(preset), d, None)
    run(rank_command(n_swarms, preset), work, log)
    log_file.close()

    # --- measure ----------------------------------------------------------
    anchor_atoms = [a for r in st.residues if r.chain == asg.anchor_chain
                    for a in r.atoms]
    partner_atoms = [a for r in st.residues if r.chain == asg.partner_chain
                     for a in r.atoms]
    native_xyz = coords(heavy(partner_atoms))
    scorer = PoseScorer(anchor_atoms, asg.ligand.atoms, face)

    poses = sorted(work.glob("swarm_*/lightdock_*.pdb"))[:top]
    rows = []
    for p in poses:
        try:
            got = read_any(p)
        except Exception as exc:  # noqa: BLE001
            log(f"could not read {p}: {exc}")
            continue
        atoms = heavy([a for r in got.residues for a in r.atoms])
        if len(atoms) != len(native_xyz):
            # LightDock's output should hold the same atoms it was given. If
            # it does not, RMSD against the crystal is meaningless and a
            # silently wrong number here would be the worst outcome of the
            # whole harness.
            log(f"{p.name}: {len(atoms)} atoms, expected {len(native_xyz)}")
            continue
        s = scorer.score(atoms, label=p.stem, reference_xyz=native_xyz)
        rows.append({"file": str(p.relative_to(work)), **s.as_dict()})

    if not rows:
        say("    no comparable poses came back")
        return {"entry": entry, "preset": preset.key,
                "error": "no comparable poses", "n_swarms": n_swarms}

    viable = [r for r in rows if r["viable"]]
    pool = viable or rows
    best_rmsd = min(r["rmsd_to_reference"] for r in pool)
    by_face = max(pool, key=lambda r: r["face_burial"])

    out = {
        "entry": entry,
        "ligand": asg.ligand.comp_id,
        "preset": preset.key,
        "n_swarms": n_swarms,
        "n_poses": len(rows),
        "n_viable": len(viable),
        "presented_area": round(face.presented_area, 1),
        "buried_fraction": round(face.buried_fraction, 4),
        "restraints": rst.as_dict(),
        "native_face_burial": round(
            scorer.score(partner_atoms, "crystal").face_burial, 4),
        "best_rmsd": round(best_rmsd, 2),
        "best_by_face_burial": by_face,
        "recovered": bool(best_rmsd <= cutoff),
        "designable": bool(best_rmsd <= DESIGNABLE),
        "top_pose_rmsd": round(by_face["rmsd_to_reference"], 2),
        "seconds": round(time.time() - t0, 1),
    }
    say(f"    best {best_rmsd:5.2f} A   "
        f"top-by-face {by_face['rmsd_to_reference']:5.2f} A   "
        f"{'HIT' if out['recovered'] else 'miss'}"
        f"{'  (designable)' if out['designable'] else ''}   "
        f"[{out['seconds']:.0f}s, {n_swarms} swarms, {len(rows)} poses]")
    return out


def report(rows: list[dict], preset, cutoff: float) -> None:
    good = [r for r in rows if "error" not in r]
    if not good:
        say("\n  every system failed; see the logs")
        return
    rec = sum(1 for r in good if r["recovered"])
    des = sum(1 for r in good if r["designable"])
    rms = [r["best_rmsd"] for r in good]

    say("")
    say("  " + "=" * 66)
    say(f"  {preset.label}: {preset.swarms} swarms, {preset.glowworms} "
        f"glowworms, {preset.steps} steps")
    say("  " + "=" * 66)
    say(f"    recovered (<={cutoff:.0f} A)    {rec}/{len(good)}"
        f"   = {100*rec/len(good):.0f}%")
    say(f"    designable (<={DESIGNABLE:.0f} A)   {des}/{len(good)}"
        f"   = {100*des/len(good):.0f}%")
    say(f"    median best RMSD      {np.median(rms):.2f} A")
    if len(rows) != len(good):
        say(f"    failed outright       {len(rows) - len(good)}")

    say("")
    say("    %-6s %-5s %7s %9s %11s %8s" % (
        "entry", "lig", "buried", "best RMSD", "top-by-face", "swarms"))
    for r in rows:
        if "error" in r:
            say("    %-6s %-5s   %s" % (r["entry"], "", r["error"]))
            continue
        say("    %-6s %-5s %6.0f%% %8.2f %10.2f %8d%s" % (
            r["entry"], r["ligand"], 100 * r["buried_fraction"],
            r["best_rmsd"], r["top_pose_rmsd"], r["n_swarms"],
            "  HIT" if r["recovered"] else ""))

    say("")
    say("    'best RMSD' is the closest pose anywhere in the output; "
        "'top-by-face'")
    say("    is the one this tool would actually have picked. The gap "
        "between them")
    say("    is how much the face-burial score costs you.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("entries", nargs="*")
    ap.add_argument("--preset", default="standard", choices=sorted(PRESETS))
    ap.add_argument("--all-presets", action="store_true")
    ap.add_argument("--cutoff", type=float, default=HIT)
    ap.add_argument("--cores", type=int, default=0)
    ap.add_argument("--top", type=int, default=100)
    ap.add_argument("--work", default=None)
    args = ap.parse_args()

    missing = missing_tools()
    if missing:
        say("\n  LightDock is not usable: " + ", ".join(missing) +
            " not on PATH.")
        say("  This harness is the one part of GlueForge that cannot run "
            "without it.")
        say("  Install with:  pip install lightdock")
        say("  Then:          python3 sidecar/dock.py --check")
        return 2

    entries = [e.upper() for e in args.entries] or sorted(
        p.stem for p in RAW.glob("*.cif"))
    if not entries:
        say("cache/raw is empty. Run sidecar/fetch.py first.")
        return 1

    work_root = Path(args.work or (OUT / "work"))
    work_root.mkdir(parents=True, exist_ok=True)
    presets = ([PRESETS[k] for k in ("quick", "standard", "thorough")]
               if args.all_presets else [PRESETS[args.preset]])

    all_out = {
        "machine": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "cores": args.cores or None,
        },
        "cutoff": args.cutoff,
        "designable_cutoff": DESIGNABLE,
        "presets": {},
    }

    for preset in presets:
        lo, hi = preset.minutes_estimate
        say(f"\n{preset.label}: {len(entries)} systems, roughly "
            f"{lo}-{hi} minutes each\n")
        rows = []
        for e in entries:
            try:
                r = one(e, preset, work_root, args.cutoff, args.cores, args.top)
            except Exception as exc:  # noqa: BLE001
                say(f"  {e}: FAILED {type(exc).__name__}: {exc}")
                rows.append({"entry": e, "preset": preset.key,
                             "error": f"{type(exc).__name__}: {exc}"})
                continue
            if r:
                rows.append(r)
        report(rows, preset, args.cutoff)
        all_out["presets"][preset.key] = {
            "preset": preset.as_dict(), "systems": rows,
        }

    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / "recovery.json"
    dest.write_text(json.dumps(all_out, indent=2), encoding="utf-8")
    say(f"\nwrote {dest}")
    say("\nSend that file back and the measured recovery rates go into the "
        "presets,\nwhere the docking screen shows them beside the run times.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
