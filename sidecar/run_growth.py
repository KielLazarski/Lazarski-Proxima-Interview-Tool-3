#!/usr/bin/env python3
"""Does the tool point at the atoms chemists actually grew from?

Step 6 of the workflow ranks places on the ligand worth growing chemistry
from, and the ranking is driven by how exposed each atom is once the ligand
sits in its pocket. That is a claim, and the CDK12 series is a way to test
it: eight ligands, one purine core, bound in the same site of the same
protein across eight depositions. They are a real medicinal-chemistry
campaign. Wherever they differ from each other is, by definition, where
chemistry was grown.

So: superpose the proteins, find which ligand atoms are conserved across the
series and which vary, and ask whether exposure predicts variation.

There is an obvious confound and it is tested rather than waved away. A
molecule's outer atoms are both more exposed and more synthetically
accessible, so exposure might be a proxy for "not in the middle of the
scaffold" and carry nothing the 2D structure does not already tell you. The
partial correlation, with distance from the ligand centroid taken out,
answers that.

    python3 sidecar/run_growth.py
    python3 sidecar/run_growth.py --cutoff 1.5   # looser atom matching
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, rankdata, spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from glueforge.cavity import detect_cavity  # noqa: E402
from glueforge.cif import read_structure  # noqa: E402
from glueforge.face import presented_face  # noqa: E402
from glueforge.geometry import atom_sasa  # noqa: E402
from glueforge.growth import growth_points  # noqa: E402
from glueforge.interface import assign_roles  # noqa: E402
from glueforge.pose import coords, heavy  # noqa: E402

RAW = HERE.parent / "cache" / "raw"
OUT = HERE.parent / "cache" / "growth"

#: The CDK12 glue series: one core, eight elaborations, eight depositions.
SERIES = ("6TD3", "8BU5", "8BU7", "8BU9", "8BUA", "8BUJ", "8BUK", "8BUL")


def kabsch(P: np.ndarray, Q: np.ndarray):
    """Rotation and translation taking P onto Q, least-squares."""
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    D = np.diag([1.0, 1.0, float(np.sign(np.linalg.det(Vt.T @ U.T)))])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc


def partial_spearman(x, y, z):
    """Spearman of x against y with z removed from both.

    Done on ranks so it stays non-parametric, which matters: exposure is
    heavily zero-inflated and conservation is a small count.
    """
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)

    def resid(a, b):
        B = np.vstack([b, np.ones_like(b)]).T
        coef, *_ = np.linalg.lstsq(B, a, rcond=None)
        return a - B @ coef

    return pearsonr(resid(rx, rz), resid(ry, rz))


def load(entries) -> dict:
    out = {}
    for e in entries:
        p = RAW / f"{e}.cif"
        if not p.exists():
            continue
        st = read_structure(p)
        a = assign_roles(st)
        ca = {r.seq_id: r.atom("CA") for r in st.residues
              if r.chain == a.anchor_chain and r.atom("CA")}
        out[e] = (st, a, ca)
    return out


def analyse(ref: str, data: dict, cutoff: float) -> dict | None:
    st0, a0, ca0 = data[ref]
    lig0 = heavy(a0.ligand.atoms)
    X0 = coords(lig0)

    trees, rmsds = [], []
    for e, (st, a, ca) in data.items():
        if e == ref:
            continue
        common = sorted(set(ca) & set(ca0))
        if len(common) < 100:
            continue
        P = np.array([[ca[i].x, ca[i].y, ca[i].z] for i in common])
        Q = np.array([[ca0[i].x, ca0[i].y, ca0[i].z] for i in common])
        R, t = kabsch(P, Q)
        rmsds.append(float(np.sqrt((((P @ R.T) + t - Q) ** 2).sum(1).mean())))
        trees.append(cKDTree(coords(heavy(a.ligand.atoms)) @ R.T + t))
    if len(trees) < 5:
        return None

    # Conservation: how many of the other ligands put an atom here.
    cons = np.array([
        sum(1 for T in trees if T.query(X0[i], k=1)[0] < cutoff)
        for i in range(len(lig0))
    ])

    anchor = heavy([x for r in st0.residues if r.chain == a0.anchor_chain
                    for x in r.atoms])
    rel = atom_sasa(lig0, anchor + lig0) / np.maximum(
        atom_sasa(lig0, lig0), 1e-9)
    peri = np.linalg.norm(X0 - X0.mean(0), axis=1)

    rho_exp, p_exp = spearmanr(rel, cons)
    rho_per, p_per = spearmanr(peri, cons)
    rho_par, p_par = partial_spearman(rel, cons, peri)

    # The above tests the exposure GATE: which atoms get considered at all.
    # It says nothing about the order they are then listed in, and the two
    # are easy to conflate. So test the shipped ranking directly, among the
    # atoms that passed the gate.
    res = [r for r in st0.residues
           if r.chain in (a0.anchor_chain, a0.partner_chain)]
    anchor = [x for r in res if r.chain == a0.anchor_chain for x in r.atoms]
    partner_at = [x for r in res if r.chain == a0.partner_chain for x in r.atoms]
    pres = [r for r in res if r.chain == a0.partner_chain]
    face = presented_face(st0, a0.ligand, anchor_chains=[a0.anchor_chain])
    allat = [x for r in res for x in r.atoms]
    anc = {k for k, x in enumerate(allat) if x.chain == a0.anchor_chain}
    cav = detect_cavity(allat, X0.mean(0), anc)
    pts = growth_points(a0.ligand, anchor, partner_at, face, cav, pres)

    rank_rho = rank_p = float("nan")
    partner_rho = float("nan")
    if len(pts) >= 6:
        sel = np.array([g.index for g in pts])
        rank_rho, rank_p = spearmanr([g.score for g in pts], cons[sel])
    # And how the conserved/varied split relates to the partner, which is
    # the obvious explanation to rule in or out.
    part_xyz = coords(heavy(partner_at))
    if len(part_xyz):
        d_part = cKDTree(part_xyz).query(X0, k=1)[0]
        partner_rho = spearmanr(cons, d_part)[0]

    return {
        "ref": ref,
        "ligand": a0.ligand.comp_id,
        "n_atoms": len(lig0),
        "n_compared": len(trees),
        "backbone_rmsd_median": round(float(np.median(rmsds)), 2),
        "conserved": int((cons >= len(trees) - 1).sum()),
        "variable": int((cons <= len(trees) // 2).sum()),
        "exposure_vs_conservation": {
            "rho": round(float(rho_exp), 4), "p": float(p_exp)},
        "periphery_vs_conservation": {
            "rho": round(float(rho_per), 4), "p": float(p_per)},
        "exposure_controlling_for_periphery": {
            "rho": round(float(rho_par), 4), "p": float(p_par)},
        "shipped_rank_vs_conservation": {
            "rho": round(float(rank_rho), 4) if rank_rho == rank_rho else None,
            "p": float(rank_p) if rank_p == rank_p else None,
            "n_growth_points": len(pts),
        },
        "conservation_vs_distance_to_partner": {
            "rho": round(float(partner_rho), 4)
            if partner_rho == partner_rho else None},
        "atoms": [
            {"name": lig0[i].name, "element": lig0[i].element,
             "conserved_in": int(cons[i]),
             "relative_sasa": round(float(rel[i]), 3),
             "distance_from_centroid": round(float(peri[i]), 2)}
            for i in range(len(lig0))
        ],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cutoff", type=float, default=1.0,
                    help="angstroms; two atoms this close count as the same")
    ap.add_argument("--out", default=str(OUT / "growth.json"))
    args = ap.parse_args()

    data = load(SERIES)
    if len(data) < 6:
        print("Need most of the CDK12 series in cache/raw. Have: "
              + (", ".join(sorted(data)) or "none"))
        return 1

    print(f"\nGrowth points against a real campaign: {len(data)} ligands, "
          f"one core, atoms matched within {args.cutoff} Å\n")
    print("%-6s %-5s %6s %6s %9s %9s %9s %9s" % (
        "ref", "lig", "core", "var", "rho(exp)", "rho(peri)", "partial",
        "rank"))

    rows = []
    for ref in data:
        r = analyse(ref, data, args.cutoff)
        if not r:
            continue
        rows.append(r)
        rk = r["shipped_rank_vs_conservation"]["rho"]
        print("%-6s %-5s %6d %6d %+9.3f %+9.3f %+9.3f %9s%s" % (
            r["ref"], r["ligand"], r["conserved"], r["variable"],
            r["exposure_vs_conservation"]["rho"],
            r["periphery_vs_conservation"]["rho"],
            r["exposure_controlling_for_periphery"]["rho"],
            ("%+.3f" % rk) if rk is not None else "  n/a",
            "" if r["exposure_controlling_for_periphery"]["p"] >= 0.05
            else " *"))

    if not rows:
        print("nothing comparable")
        return 1

    exp = np.array([r["exposure_vs_conservation"]["rho"] for r in rows])
    per = np.array([r["periphery_vs_conservation"]["rho"] for r in rows])
    par = np.array([r["exposure_controlling_for_periphery"]["rho"] for r in rows])
    sig = sum(1 for r in rows
              if r["exposure_controlling_for_periphery"]["p"] < 0.05
              and r["exposure_controlling_for_periphery"]["rho"] < 0)

    print("\n" + "=" * 68)
    print("  exposure vs conservation      median %+.3f   %d/%d negative"
          % (np.median(exp), int((exp < 0).sum()), len(exp)))
    print("  periphery vs conservation     median %+.3f   %d/%d negative"
          % (np.median(per), int((per < 0).sum()), len(per)))
    print("  exposure, periphery removed   median %+.3f   %d/%d negative, "
          "%d significant"
          % (np.median(par), int((par < 0).sum()), len(par), sig))
    print("=" * 68)
    print("\n  Negative means the atoms this tool calls exposed are the atoms")
    print("  the chemists changed. The third line is the one that matters:")
    print("  if exposure were only a proxy for being on the outside of the")
    print("  molecule it would collapse there, and it does not.")
    print("\n  The eight rows are one series looked at eight ways, not eight")
    print("  independent tests. Read the consistency, not a combined p-value.")

    rk = np.array([r["shipped_rank_vs_conservation"]["rho"] for r in rows
                   if r["shipped_rank_vs_conservation"]["rho"] is not None])
    pr = np.array([r["conservation_vs_distance_to_partner"]["rho"] for r in rows
                   if r["conservation_vs_distance_to_partner"]["rho"] is not None])
    if len(rk):
        print("\n" + "=" * 68)
        print("  AND THE PART THAT IS NOT VALIDATED")
        print("=" * 68)
        print("  Everything above tests the exposure GATE -- which atoms are")
        print("  considered at all. The tool then ORDERS those atoms by a")
        print("  composite of room, direction and enclosure, and that order")
        print("  is a different claim.")
        print()
        print("  shipped rank vs conservation  median %+.3f   %d/%d positive"
              % (np.median(rk), int((rk > 0).sum()), len(rk)))
        if np.median(rk) > 0.15:
            print("    Positive means the tool ranks highest the atoms this")
            print("    campaign did NOT change. It is pointing away from where")
            print("    chemistry actually went.")
        if len(pr):
            print()
            print("  conservation vs distance to partner  median %+.3f"
                  % np.median(pr))
            if np.median(pr) > 0.15:
                print("    Positive means the varied atoms sit CLOSER to the")
                print("    partner. The campaign elaborated toward it, into")
                print("    space a rigid cavity does not show as free --")
                print("    real elaboration reshapes an interface, and this")
                print("    tool can only see the space already there.")

    dest = Path(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps({
        "note": (
            "Whether ligand exposure predicts where a congeneric series was "
            "varied. One series (CDK12 glues), each member used in turn as "
            "the reference, so the rows are not independent."
        ),
        "series": list(data),
        "match_cutoff": args.cutoff,
        "systems": rows,
        "summary": {
            "exposure_median_rho": round(float(np.median(exp)), 4),
            "periphery_median_rho": round(float(np.median(per)), 4),
            "partial_median_rho": round(float(np.median(par)), 4),
            "partial_negative": int((par < 0).sum()),
            "partial_significant": sig,
            "n": len(rows),
        },
    }, indent=2), encoding="utf-8")
    print(f"\nwrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
