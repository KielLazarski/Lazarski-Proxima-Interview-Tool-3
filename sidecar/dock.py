#!/usr/bin/env python3
"""Run a restrained docking, from two structures to scored poses.

    python3 sidecar/dock.py A.pdb B.pdb --ligand RC8
    python3 sidecar/dock.py A.cif B.cif --preset quick
    python3 sidecar/dock.py --check                 # is LightDock usable?
    python3 sidecar/dock.py A.pdb B.pdb --dry-run   # print the plan, run nothing

What happens, in order: protein A is read and its presented face measured
against A alone; the rim of that face becomes LightDock residue restraints;
LightDock places B; every returned pose is re-scored here on how much of the
presented face it buries and whether it overlaps the ligand.

That last step is not optional. LightDock is given protein only -- its scoring
function has no parameters for an arbitrary small molecule -- so the ligand's
pocket is empty as far as the docking run is concerned, and a pose that
settles into it buries the presented face completely by occupying it. Those
poses are the best-scoring ones in an unfiltered run and they are nonsense.

This file is the half of the pipeline that needs LightDock installed. The
measurement of whether the restraint is worth anything in the first place is
`run_decoys.py`, which needs nothing.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from glueforge.cif import read_any  # noqa: E402
from glueforge.face import default_anchor_chains, presented_face  # noqa: E402
from glueforge.pose import PoseScorer, coords, heavy  # noqa: E402
from glueforge.restraints import (  # noqa: E402
    DEFAULT_PRESET,
    PRESETS,
    TOOLS,
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


def say(msg: str = "") -> None:
    print(msg, flush=True)


def check() -> int:
    missing = missing_tools()
    say("\nLightDock")
    for t in TOOLS:
        where = shutil.which(t)
        say("  %-32s %s" % (t, where or "NOT FOUND"))
    if missing:
        say(
            "\n  %d of %d not on PATH.\n"
            "  Install with:  pip install lightdock\n"
            "  If pip says it is already installed, its scripts are not on\n"
            "  your PATH -- on Windows that is usually the Scripts folder of\n"
            "  the Python you installed it into." % (len(missing), len(TOOLS))
        )
        return 1
    say("\n  all present")
    return 0


def pick_ligand(structure, wanted: str | None):
    if wanted:
        for l in structure.ligands:
            if l.comp_id.upper() == wanted.upper():
                return l
        raise SystemExit(
            f"no ligand {wanted} in protein A. Found: "
            + (", ".join(sorted({l.comp_id for l in structure.ligands}))
               or "none")
        )
    best = structure.largest_ligand()
    if best is None:
        raise SystemExit(
            "protein A has no ligand. The method starts from a molecule "
            "already bound to A; without one there is nothing to aim at."
        )
    return best


def run_builtin(sa, sb, lig, anchor_chains, partner_chains, face, rst,
                work, args) -> int:
    """The same pipeline without LightDock, using the engine in `fftdock`.

    The re-scoring step that follows a LightDock run is not needed in the
    same way here -- the ligand is inside the engine's own receptor grid, so
    a pose that settles into its pocket is penalised while it is generated --
    but it is still run, because face burial is what the caller wants
    reported and `PoseScorer` is the measured way to get it.
    """
    from glueforge.fftdock import PRESETS as FFT_PRESETS
    from glueforge.fftdock import dock as fft_dock

    settings = FFT_PRESETS[args.preset]
    anchor_atoms = [a for r in sa.residues if r.chain in set(anchor_chains)
                    for a in r.atoms]
    partner_atoms = [a for r in sb.residues if r.chain in set(partner_chains)
                     for a in r.atoms]
    if not partner_atoms:
        raise SystemExit(
            "protein B has none of the chains asked for: "
            + ", ".join(partner_chains)
        )

    n_par = len([a for a in partner_atoms if a.element.upper() != "H"])
    lo, hi = settings.minutes_estimate
    say(f"\n  engine            built-in (numpy and scipy only)")
    say(f"  preset            {settings.label} — "
        f"{settings.rotations} orientations at {settings.spacing} A, "
        f"then {settings.focus} revisited")
    say(f"  rough time        up to {lo}-{hi} minutes")
    # Two warnings, both from the same cause and both measured rather than
    # guessed: this engine finds the right answer when the search space is
    # genuinely narrow, and loses it when it is not.
    if n_par > 2500:
        say(f"\n  ! This partner has {n_par} heavy atoms, which is large for "
            "this engine.\n    The search can only sample orientations so "
            "finely, and a large body\n    swings further for the same "
            "angular error -- on the benchmark's\n    6500-atom partners the "
            "best pose found was still 11 A out. LightDock,\n    via WSL on "
            "Windows, is the better tool for a pair this size.")
    if face.buried_fraction < 0.55:
        say(f"\n  ! The ligand is only {100*face.buried_fraction:.0f}% buried, "
            "so the face it presents is\n    broad and flat rather than a "
            "pocket mouth. That weakens the restraint\n    for any engine, "
            "and this one leans on the restraint harder than\n    LightDock "
            "does. Expect the search to spread.")

    last = [-1]

    def progress(done: int, total: int) -> None:
        pct = int(100 * done / max(1, total))
        if pct >= last[0] + 10:
            last[0] = pct - pct % 10
            say(f"    {last[0]:3d}%")

    say("\n  [searching]")
    poses = fft_dock(anchor_atoms, lig.atoms, partner_atoms, face,
                     settings, progress=progress)
    if not poses:
        say("\n  Nothing satisfied the restraint. That usually means the "
            "presented\n  face is too small or too enclosed for this partner "
            "to reach.")
        return 4

    say(f"\n  [re-scoring against the ligand]  {len(poses)} poses")
    scorer = PoseScorer(anchor_atoms, lig.atoms, face)
    rows = [(p, scorer.score(p.apply_atoms(partner_atoms), label=p.label))
            for p in (poses[: args.top] if args.top else poses)]
    # Left in the engine's own order rather than re-sorted by face burial.
    #
    # Face burial is the right way to rank a set of poses that have not
    # already been filtered on it -- it is how a glue interface is told from
    # a merely favourable one, and that is what it was validated for. But
    # every pose here has already cleared a hard face-contact filter, and
    # ranking a pre-filtered set by the filter's own criterion sorts it by
    # how many atoms each pose piles onto the face rather than by whether it
    # is right. On 6H0F, the one benchmark system where this engine's search
    # did find the answer, face burial chose a 16.7 A pose and the engine's
    # order chose the 2.6 A one.
    #
    # Across all six the two rules are nearly level (compare_rankings.py puts
    # them at 14.4 and 14.9 A of median regret), so this is not a large
    # claim -- but where it mattered it mattered a lot, and it costs nothing.
    kept = [(p, s) for p, s in rows if s.viable]

    say(f"  {len(rows)} scored, {len(kept)} viable "
        f"({len(rows) - len(kept)} rejected: clashing or not touching)")
    say("\n  %-16s %8s %8s %7s %7s" % ("pose", "face", "contacts", "clash", "lig"))
    for p, s in (kept or rows)[:12]:
        say("  %-16s %7.1f%% %8d %7d %7d" % (
            p.label, 100 * s.face_burial, s.contacts, s.clashes,
            s.ligand_clashes))

    # The poses themselves, so they can be looked at rather than taken on
    # trust. Written as PDB because that is what every viewer opens.
    from glueforge.restraints import write_pdb
    out_dir = work / "poses"
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for p, s in (kept or rows)[:20]:
        dest = out_dir / f"{p.label}.pdb"
        write_pdb(p.apply_atoms(partner_atoms), dest)
        written.append(str(dest.relative_to(work)))

    # And the thing to look at them against. A pose file on its own is
    # protein B floating in space: the whole content of a docking result is
    # where it sits relative to protein A, and without the receptor beside
    # it there is nothing in the folder that shows that. The LightDock path
    # gets this from `split_for_docking`, which the built-in path returns
    # before reaching, so it was shipping twenty files nobody could open.
    #
    # The ligand goes in with the anchor, because it is what the partner was
    # aimed at and leaving it out would hide the only part of protein A the
    # result is really about.
    receptor = work / "receptor.pdb"
    write_pdb(anchor_atoms + lig.atoms, receptor)
    (work / "face.json").write_text(
        json.dumps(face.as_dict(), indent=2), encoding="utf-8")
    say(f"\n  wrote {receptor.name} — open it with any pose in "
        f"{out_dir.name}/ to see the placement")

    (work / "scored.json").write_text(json.dumps({
        "engine": "builtin",
        "settings": settings.as_dict(),
        "face": face.as_dict(),
        "restraints": rst.as_dict(),
        "n_poses": len(rows),
        "poses": [
            {"label": p.label, **s.as_dict()} for p, s in rows
        ],
        "written": written,
        "receptor": "receptor.pdb",
    }, indent=2), encoding="utf-8")
    say(f"\n  wrote {work / 'scored.json'} and {len(written)} poses in "
        f"{out_dir}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("a", nargs="?", help="protein A, with the ligand bound")
    ap.add_argument("b", nargs="?", help="protein B")
    ap.add_argument("--ligand", help="ligand code, if the automatic pick is wrong")
    ap.add_argument("--anchor", nargs="*", help="chains of A to treat as the anchor")
    ap.add_argument("--partner", nargs="*", help="chains of B to dock")
    ap.add_argument("--engine", default="auto",
                    choices=["auto", "lightdock", "builtin"],
                    help="auto uses LightDock when it is installed")
    ap.add_argument("--preset", default=DEFAULT_PRESET, choices=sorted(PRESETS))
    ap.add_argument("--cores", type=int, default=0)
    ap.add_argument("--out", default=None, help="working directory")
    ap.add_argument("--check", action="store_true", help="is LightDock usable?")
    ap.add_argument("--dry-run", action="store_true",
                    help="write the inputs and print the commands, run nothing")
    # All of them by default. The engine returns each pose and its
    # refinement ordered by its own score, so a prefix drops the lower half
    # of that pairing -- and on the benchmark the closest pose to the right
    # answer was frequently not in the top fifty. Re-scoring is cheap next
    # to the search that produced them.
    ap.add_argument("--top", type=int, default=0,
                    help="re-score only the first N poses; 0 means all")
    args = ap.parse_args()

    if args.check:
        return check()
    if not args.a or not args.b:
        ap.error("give protein A and protein B, or --check")

    preset = PRESETS[args.preset]
    work = Path(args.out or f"runs/dock-{time.strftime('%Y%m%d-%H%M%S')}")
    work.mkdir(parents=True, exist_ok=True)

    # --- read and measure ------------------------------------------------
    say(f"\n  reading {args.a}")
    sa = read_any(args.a)
    say(f"  reading {args.b}")
    sb = read_any(args.b)

    lig = pick_ligand(sa, args.ligand)
    anchor_chains = args.anchor or default_anchor_chains(sa, lig)
    if not anchor_chains:
        raise SystemExit("the ligand in A touches no protein chain")
    partner_chains = args.partner or sorted(sb.chain_sizes())

    face = presented_face(sa, lig, anchor_chains=anchor_chains)
    for w in face.warnings:
        say(f"  ! {w}")
    if not face.exposed:
        raise SystemExit(
            "No part of the ligand is exposed once A is in place. A fully "
            "enclosed ligand cannot be reached by a second protein, so this "
            "pair has no glue to make and docking it would not tell you "
            "anything."
        )

    say(f"\n  ligand            {lig.comp_id} ({len(heavy(lig.atoms))} heavy atoms)")
    say(f"  anchor chains     {', '.join(anchor_chains)}")
    say(f"  presented face    {face.presented_area:.1f} A^2 over "
        f"{len(face.exposed)} atoms, {100*face.buried_fraction:.0f}% buried")

    # The one thing worth warning about before spending hours: a ligand that
    # is barely buried presents a large flat face, which restrains a docking
    # run much more weakly than a ligand in a pocket.
    if face.buried_fraction < 0.55:
        say(
            "\n  ! This ligand sits in a shallow groove rather than a pocket\n"
            f"    ({100*face.buried_fraction:.0f}% buried). On the benchmark, "
            "the restraint is much\n"
            "    weaker in that case -- many wrong placements also cover a "
            "flat face.\n"
            "    See `run_decoys.py` for what that costs."
        )

    rst = restraints_from_face(face)
    say(f"  restraints        {len(rst.active)} active, {len(rst.passive)} passive")

    # --- which engine ------------------------------------------------------
    engine = args.engine
    if engine == "auto":
        engine = "lightdock" if not missing_tools() else "builtin"
        if engine == "builtin":
            say("\n  LightDock is not installed, so the built-in engine "
                "will be used.\n  Its measured recovery is in "
                "cache/recovery/builtin.json -- read that before\n  trusting "
                "a pose. The two engines have not been compared here, so "
                "this\n  is not a fallback to something known to be worse, "
                "only to the one that\n  has been measured. "
                "`dock.py --check` says how to get LightDock.")

    if engine == "builtin":
        return run_builtin(sa, sb, lig, anchor_chains, partner_chains,
                           face, rst, work, args)

    # --- write the inputs -------------------------------------------------
    files = split_for_docking(
        sa, anchor_chains, partner_chains, work,
        ligand=lig, partner_structure=sb,
    )

    (work / "restraints.list").write_text(rst.text(), encoding="ascii")
    (work / "face.json").write_text(
        json.dumps(face.as_dict(), indent=2), encoding="utf-8")
    (work / "plan.json").write_text(json.dumps({
        "preset": preset.as_dict(),
        "restraints": rst.as_dict(),
        "files": files,
        "ligand": lig.comp_id,
        "anchor_chains": anchor_chains,
        "partner_chains": partner_chains,
    }, indent=2), encoding="utf-8")

    rec = Path(files["receptor"]).name
    par = Path(files["partner"]).name
    steps = [
        ("setup", setup_command(Path(rec), Path(par),
                                Path("restraints.list"), preset)),
        ("dock", run_command(preset, args.cores)),
    ]

    lo, hi = preset.minutes_estimate
    say(f"\n  preset            {preset.label} "
        f"({preset.swarms} swarms, {preset.glowworms} glowworms, "
        f"{preset.steps} steps)")
    say(f"  rough time        {lo}-{hi} minutes on 8 cores")
    say(f"  working directory {work}")

    if args.dry_run:
        say("\n  commands (not run):")
        for _, cmd in steps:
            say("    " + " ".join(cmd))
        say("    cd swarm_N && " + " ".join(
            conformations_command(Path("../" + rec), Path("../" + par), preset)))
        say("    cd swarm_N && " + " ".join(cluster_command(preset)))
        say("    " + " ".join(rank_command(preset.swarms, preset)))
        return 0

    missing = missing_tools()
    if missing:
        say("\n  LightDock is not usable: " + ", ".join(missing) +
            " not on PATH.\n  The inputs above are written and valid; run "
            "`dock.py --check` for the fix.")
        return 2

    # --- run --------------------------------------------------------------
    log_path = work / "lightdock.log"
    log_file = log_path.open("w", encoding="utf-8")

    def log(line: str) -> None:
        log_file.write(line + "\n")
        log_file.flush()
        say("    " + line)

    for name, cmd in steps:
        say(f"\n  [{name}]")
        code = run(cmd, work, log)
        if code != 0:
            log_file.close()
            say(f"\n  {name} failed (exit {code}). Full output: {log_path}")
            return 3

    n_swarms = count_swarms(work)
    say(f"\n  {n_swarms} swarms survived the restraints")
    if n_swarms == 0:
        say("  No swarm satisfied the restraints. That usually means the "
            "active set is too strict for this pair; try --preset thorough "
            "or widen ACTIVE_WITHIN.")
        log_file.close()
        return 4

    say("\n  [conformations and clustering]")
    for i in range(n_swarms):
        d = work / f"swarm_{i}"
        if not d.is_dir():
            continue
        run(conformations_command(Path("../" + rec), Path("../" + par), preset),
            d, None)
        run(cluster_command(preset), d, None)
    say(f"  done, {n_swarms} swarms")

    say("\n  [ranking]")
    run(rank_command(n_swarms, preset), work, log)
    log_file.close()

    # --- re-score ---------------------------------------------------------
    say("\n  [re-scoring against the ligand]")
    anchor_atoms = [a for r in sa.residues if r.chain in set(anchor_chains)
                    for a in r.atoms]
    scorer = PoseScorer(anchor_atoms, lig.atoms, face)

    found = sorted(work.glob("swarm_*/lightdock_*.pdb"))
    poses = found[: args.top] if args.top else found
    if not poses:
        say("  no pose files found; look in " + str(work))
        return 5

    rows = []
    for p in poses:
        try:
            st = read_any(p)
        except Exception as exc:  # noqa: BLE001
            say(f"  ! could not read {p.name}: {exc}")
            continue
        atoms = [a for r in st.residues for a in r.atoms]
        rows.append((p, scorer.score(atoms, label=p.stem)))

    rows.sort(key=lambda t: -t[1].face_burial)
    kept = [(p, s) for p, s in rows if s.viable]

    say(f"\n  {len(rows)} poses scored, {len(kept)} viable "
        f"({len(rows) - len(kept)} rejected: clashing or not touching)")
    in_pocket = sum(1 for _, s in rows if s.ligand_clashes > 2)
    if in_pocket:
        say(f"  {in_pocket} of them had settled into the ligand's own pocket, "
            "which is what the re-scoring is for")

    say("\n  %-28s %8s %8s %7s %7s" % ("pose", "face", "contacts", "clash", "lig"))
    for p, s in (kept or rows)[:12]:
        say("  %-28s %7.1f%% %8d %7d %7d" % (
            p.stem[:28], 100 * s.face_burial, s.contacts, s.clashes,
            s.ligand_clashes))

    (work / "scored.json").write_text(json.dumps({
        "preset": preset.as_dict(),
        "face": face.as_dict(),
        "restraints": rst.as_dict(),
        "n_swarms": n_swarms,
        "poses": [
            {"file": str(p.relative_to(work)), **s.as_dict()} for p, s in rows
        ],
    }, indent=2), encoding="utf-8")
    say(f"\n  wrote {work / 'scored.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
