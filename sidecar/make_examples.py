#!/usr/bin/env python3
"""Prepare the worked examples that ship with GlueForge.

Shipping raw PDB entries would mean shipping tens of megabytes and then
asking the user to do the trimming that the tool itself has to teach them
anyway. These are prepared instead: each example is the chains that matter,
already split into protein A with its ligand and protein B on its own, so
one click gives a pair that is ready to run.

The measurements in each example's description are produced here, by running
the same code the tool runs. Nothing is written by hand, so an example whose
numbers would change cannot keep describing itself the old way -- which is
the failure mode of curated demos, and the reason this is a script rather
than a JSON file somebody edits.

    python3 sidecar/make_examples.py           # rebuild cache/examples
    python3 sidecar/make_examples.py --check   # verify without writing
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from glueforge.cavity import detect_cavity  # noqa: E402
from glueforge.cif import read_any  # noqa: E402
from glueforge.face import presented_face  # noqa: E402
from glueforge.gate import assess as judge  # noqa: E402
from glueforge.gate import cavity_metrics  # noqa: E402
from glueforge.interface import assign_roles  # noqa: E402
from glueforge.pose import PoseScorer, coords, heavy  # noqa: E402
from glueforge.restraints import write_pdb  # noqa: E402

RAW = HERE.parent / "cache" / "raw"
OUT = HERE.parent / "cache" / "examples"

#: The examples, and what each is here to show. `kind` decides how the tool
#: offers it: "ternary" splits into A and B so it can be assessed as given or
#: re-docked; "binary" has no partner and is a starting point for docking
#: something new.
#:
#: Chosen for what they demonstrate and for being different from each other.
#: The benchmark holds eight CDK12 entries and four DCAF15 ones; shipping
#: them all would be a longer list showing less.
EXAMPLES = [
    {
        "id": "6H0F",
        "kind": "ternary",
        "title": "CRBN with lenalidomide, recruiting IKZF1",
        "shows": (
            "The cleanest case in the set. The molecule sits deep in "
            "cereblon's pocket and presents a small, strongly directional "
            "face, and the zinc finger covers almost all of it -- 94% of "
            "the presented face buried, the highest of the fourteen. It is "
            "also the only system where the built-in docking engine "
            "recovered the crystal placement, so it is the one example "
            "where every stage of the tool can be run end to end."
        ),
    },
    {
        "id": "5FQD",
        "kind": "ternary",
        "title": "CRBN with lenalidomide, recruiting CK1\u03b1",
        "shows": (
            "The same ligase and the same chemistry recruiting an entirely "
            "different protein. That is the premise of molecular glue in "
            "one pair of files, and it is worth opening beside 6H0F: the "
            "presented face is nearly identical, 74.8 against 72.6 square "
            "angstroms, and what arrives to sit on it is not."
        ),
    },
    {
        "id": "8BU5",
        "kind": "ternary",
        "title": "CDK12\u2013cyclin K glued to DDB1",
        "shows": (
            "A glue with no ligase pocket involved at all: the molecule "
            "binds CDK12 and DDB1 is recruited straight onto it. Useful as "
            "a check that the gate is measuring interfaces rather than "
            "having learned to recognise cereblon, since nothing about "
            "this system resembles the first two."
        ),
    },
    {
        "id": "6UD7",
        "kind": "ternary",
        "title": "DCAF15 with E7820, recruiting RBM39",
        "shows": (
            "The largest cavity in the benchmark at 1274 cubic angstroms, "
            "which makes it the clearest illustration of what the cavity "
            "view is for -- that is the space chemistry would have to "
            "fill. One caveat is visible in its own numbers: the face has "
            "a directionality of 0.34, barely above the threshold where "
            "the tool stops trusting the approach direction at all."
        ),
    },
    {
        "id": "6Q0W",
        "kind": "ternary",
        "title": "DCAF15 with indisulam, recruiting RBM39 \u2014 the hard case",
        "shows": (
            "The same two proteins as 6UD7 with a different ligand, and a "
            "different verdict: borderline rather than pass. The ligand "
            "lies in a shallow groove rather than a pocket, so it presents "
            "343 square angstroms of broad flat face at 36% buried, and a "
            "face like that restrains a docking search weakly -- many "
            "wrong placements cover it. Included on purpose. A tool worth "
            "trusting has to be seen saying \"not sure\", and the pair "
            "with 6UD7 shows the verdict following the chemistry rather "
            "than the proteins."
        ),
    },
]


def measure(path: Path) -> dict:
    """Everything the description quotes, computed the way the tool does."""
    st = read_any(path)
    a = assign_roles(st)
    res = [r for r in st.residues
           if r.chain in (a.anchor_chain, a.partner_chain)]
    anchor = [x for r in res if r.chain == a.anchor_chain for x in r.atoms]
    partner = [x for r in res if r.chain == a.partner_chain for x in r.atoms]
    face = presented_face(st, a.ligand, anchor_chains=[a.anchor_chain])

    lig_xyz = coords(heavy(a.ligand.atoms))
    metrics = cavity_metrics(heavy(anchor), heavy(partner), lig_xyz)
    scorer = PoseScorer(anchor, a.ligand.atoms, face)
    metrics["face_burial"] = scorer.score(partner, "crystal").face_burial
    verdict = judge(metrics, face=face.as_dict())

    allat = [x for r in res for x in r.atoms]
    cav = detect_cavity(
        allat, lig_xyz.mean(axis=0),
        {i for i, x in enumerate(allat) if x.chain == a.anchor_chain},
    )
    return {
        "entry": st.entry_id,
        "ligand": a.ligand.comp_id,
        "anchor_chain": a.anchor_chain,
        "partner_chain": a.partner_chain,
        "anchor_name": a.anchor_name,
        "partner_name": a.partner_name,
        "anchor_residues": len([r for r in res if r.chain == a.anchor_chain]),
        "partner_residues": len([r for r in res if r.chain == a.partner_chain]),
        "ligand_heavy_atoms": len(heavy(a.ligand.atoms)),
        "presented_area": round(face.presented_area, 1),
        "buried_fraction": round(face.buried_fraction, 3),
        "directionality": round(face.directionality, 3),
        "face_burial_by_partner": round(float(metrics["face_burial"]), 3),
        "cavity_volume": round(float(cav.summary()["volume"]), 1),
        "cavity_points": int(len(cav.points)),
        "verdict": (verdict.get("verdict") if isinstance(verdict, dict)
                    else getattr(verdict, "verdict", None)),
    }


def build(spec: dict, check: bool) -> dict | None:
    src = RAW / f"{spec['id']}.cif"
    if not src.exists():
        print(f"  {spec['id']}: not in cache/raw, skipped")
        return None

    m = measure(src)
    st = read_any(src)

    # One file holding the whole ternary complex, not two pre-split ones.
    #
    # Splitting them first looked tidier and was wrong. Assessing a complex
    # is the validated part of this tool -- the gate, the cavity, the face,
    # the growth points -- and it needs both proteins in one structure. An
    # example that arrives pre-split cannot be assessed at all until it has
    # been through the docking engine, which is the part that recovers one
    # pose in six. That would route every first-time user through the
    # weakest stage to reach the strongest one.
    #
    # Shipped assembled, the chain picker then makes the other route
    # available anyway: drop the partner chain and dock something else at
    # the same face.
    both = [x for r in st.residues
            if r.chain in (m["anchor_chain"], m["partner_chain"])
            for x in r.atoms]
    for lig in st.ligands:
        if lig.chain in (m["anchor_chain"], m["partner_chain"]):
            both += list(lig.atoms)
    a_atoms, b_atoms = both, []

    out = OUT / spec["id"]
    record = {
        **{k: v for k, v in spec.items()},
        **m,
        "files": {"A": "complex.pdb", "B": None},
    }
    if check:
        return record

    out.mkdir(parents=True, exist_ok=True)
    write_pdb(both, out / "complex.pdb")
    (out / "example.json").write_text(
        json.dumps(record, indent=2), encoding="utf-8")
    size = sum(f.stat().st_size for f in out.glob("*"))
    print(f"  {spec['id']}: {len(a_atoms)} atoms "
          f"(chains {m['anchor_chain']}+{m['partner_chain']}), "
          f"{size // 1024} kB  [{m['verdict']}]")
    return record


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true",
                    help="measure and print, write nothing")
    args = ap.parse_args()

    if not args.check and OUT.exists():
        shutil.rmtree(OUT)
    print(f"\n  {len(EXAMPLES)} examples\n")
    built = [r for r in (build(s, args.check) for s in EXAMPLES) if r]
    if not built:
        print("\n  nothing built; run sidecar/fetch.py first")
        return 1

    if not args.check:
        (OUT / "index.json").write_text(
            json.dumps({"examples": built}, indent=2), encoding="utf-8")
        total = sum(f.stat().st_size for f in OUT.rglob("*") if f.is_file())
        print(f"\n  wrote {OUT} ({total // 1024} kB)")

    print("\n  %-6s %-5s %8s %7s %8s %9s  %s" % (
        "entry", "lig", "face", "buried", "cavity", "verdict", "shows"))
    for r in built:
        print("  %-6s %-5s %7.1fÅ² %6.0f%% %7.0fÅ³ %9s  %s" % (
            r["entry"], r["ligand"], r["presented_area"],
            100 * r["buried_fraction"], r["cavity_volume"],
            r["verdict"], r["title"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
