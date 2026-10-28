#!/usr/bin/env python3
"""Which way of choosing one pose from the output is least wrong?

The search and the ranking fail separately and have to be fixed separately.
On the first benchmark run the search put a 1.67 A pose in its output for
6H0F and the tool then picked an 18 A one, so the search was working and the
ranking was not.

This reads `cache/recovery/builtin.json`, which holds every scored pose of
every system, and replays the choice under several rules. No docking is
re-run, so comparing ten rules costs the same as comparing one, and the rules
are being compared on the same poses rather than on separate runs that
differ by chance as well as by rule.

What it cannot tell you: whether a rule that wins here wins anywhere else.
Six systems is a small sample and these are its own training data. A rule
that only just wins is not a result.

    python3 sidecar/compare_rankings.py
    python3 sidecar/compare_rankings.py --json other/run.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DEFAULT = HERE.parent / "cache" / "recovery" / "builtin.json"

HIT = 5.0


def say(msg: str = "") -> None:
    print(msg, flush=True)


#: Each rule takes the list of scored poses and returns the one it would
#: choose. `engine order` is the order the poses arrive in, which after
#: refinement is the continuous score the engine optimised.
RULES: dict[str, object] = {
    "engine order": lambda rows: rows[0],
    "face burial": lambda rows: max(rows, key=lambda r: r["face_burial"]),
    "contacts": lambda rows: max(rows, key=lambda r: r["contacts"]),
    "contacts less clashes": lambda rows: max(
        rows, key=lambda r: r["contacts"] - 4 * r["clashes"]),
    "restraints satisfied": lambda rows: max(
        rows, key=lambda r: (r.get("restraint_satisfaction", 0.0),
                             r["contacts"])),
    "face burial, then contacts": lambda rows: max(
        rows, key=lambda r: (round(r["face_burial"], 2), r["contacts"])),
    "buried area": lambda rows: max(rows, key=lambda r: r["buried_area"]),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", default=str(DEFAULT))
    ap.add_argument("--cutoff", type=float, default=HIT)
    args = ap.parse_args()

    path = Path(args.json)
    if not path.exists():
        say(f"no {path}. Run sidecar/run_builtin_recovery.py first.")
        return 1
    data = json.loads(path.read_text(encoding="utf-8"))
    systems = [s for s in data.get("systems", []) if "error" not in s
               and s.get("poses")]
    if not systems:
        say("nothing scored in that file.")
        return 1

    say(f"\n  {len(systems)} systems, choosing one pose from each\n")

    # The ceiling: how often the right answer was in the output at all. No
    # ranking rule can beat this, and the gap between a rule and this number
    # is what better ranking could still buy.
    reachable = sum(1 for s in systems
                    if min(p["rmsd_to_reference"] for p in s["poses"])
                    <= args.cutoff)
    say(f"    {'in the output at all':<28} {reachable}/{len(systems)}"
        "   <- nothing can beat this")
    say("")

    # The hit rate alone is nearly useless here: when the search tops out
    # above the cutoff, every rule scores zero and they look identical. What
    # separates them is how much worse than the best available pose each one
    # chooses -- a rule that picks 20 A when 7 A was in the list is failing
    # even though neither would have counted as a hit.
    bests = [min(p["rmsd_to_reference"] for p in s["poses"])
             for s in systems]

    results = []
    for name, rule in RULES.items():
        picks = []
        for s in systems:
            viable = [p for p in s["poses"] if p["viable"]] or s["poses"]
            picks.append(rule(viable)["rmsd_to_reference"])
        hits = sum(1 for r in picks if r <= args.cutoff)
        gaps = [p - b for p, b in zip(picks, bests)]
        results.append((hits, float(np.median(gaps)), name, picks, gaps))

    # Ordered by what the choice costs, not by the hit rate, for the reason
    # above.
    results.sort(key=lambda t: (t[1], -t[0]))
    say(f"    {'rule':<28} {'hits':>5}  {'median miss':>11}   picks")
    for hits, gap, name, picks, gaps in results:
        say(f"    {name:<28} {hits}/{len(systems)}   {gap:8.1f} A   [" +
            " ".join(f"{r:5.1f}" for r in picks) + "]")
    say("")
    say("    'median miss' is how much further from the right answer the")
    say("    rule's choice is than the best pose that was available to it.")
    say("    Zero would mean the ranking costs nothing and every remaining")
    say("    failure belongs to the search.")
    say("")
    say("    best available:              " + " " * 14 +
        "[" + " ".join(f"{b:5.1f}" for b in bests) + "]")

    say("")
    say("    Columns are the systems in file order: "
        + " ".join(s["entry"] for s in systems))
    say("")
    say("    These six are the data the rules are being compared on, so a")
    say("    rule that only just wins has not been shown to be better. The")
    say("    number worth acting on is the gap between the best rule and")
    say("    'in the output at all'.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
