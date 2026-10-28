"""Turn a complex into the shape the page's viewer already draws.

`export_web.py` builds this blob for the fourteen benchmark systems, but it
builds it as a by-product of the frustration analysis -- eight hundred decoys
and a published potential, minutes of work per system. A docked pose needs
none of that to be *drawn*. The drawing wants coordinates, a backbone trace,
secondary structure, which residues are at the interface, the ligand, and the
cavity; all of which are geometry.

So this is the drawing half on its own, computed from coordinates alone and
fast enough to run while someone waits. What it deliberately does not
produce is the frustration analysis: no per-contact energies, no exit
vectors, no comparison of proximal against distal. Those are measurements
with a method behind them and a published potential underneath, and
inventing plausible-looking numbers for them so that every style in the
viewer has something to draw would be the worst thing this file could do.
The viewer is told what is missing and offers only what it can honestly
show.

The contract is `cache/web/<entry>.json`, which the validation browser reads,
and `test_viewblob.py` checks this against a real one of those rather than
against a description of one.
"""

from __future__ import annotations

import numpy as np

from .cavity import detect_cavity
from .cif import Atom, Residue
from .geometry import atom_sasa, secondary_structure

#: Heavy-atom distance within which two residues are in contact. The same 5 A
#: used everywhere else here, so "interface" means one thing across the tool.
CONTACT = 5.0

#: Backbone atoms, which the sidechain style does not draw.
BACKBONE = frozenset({"N", "CA", "C", "O", "OXT"})


def _r3(a) -> list:
    arr = np.asarray(a, dtype=np.float64)
    # NaN and inf are not valid JSON and this is parsed by the page.
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    return np.round(arr, 3).tolist()



def cavity_blob(cav, centre, r3=_r3, max_points: int = 1400) -> dict:
    """The cavity, thinned for the browser.

    A cavity runs to a few thousand grid points and the page only needs
    enough to draw the shape and its fields, so this subsamples on a stride
    rather than shipping the lot.

    Lives here rather than in `export_web.py`, where it was written, because
    the benchmark export and a live run have to produce byte-identical shapes
    for the same cavity -- the page has one cavity view and it reads one set
    of keys. Two copies of this would drift, and the drift would show up as a
    view that works for the fourteen cached systems and not for anything the
    user docked themselves.
    """
    n = len(cav.points)
    if n == 0:
        return {"points": 0, "note": "no cavity detected at this site"}
    stride = max(1, n // max_points)
    sel = slice(None, None, stride)
    pts = cav.points[sel] - centre
    return {
        "summary": cav.summary(),
        "points": int(len(pts)),
        "stride": stride,
        "total_points": n,
        "spacing": cav.spacing,
        "x": r3(pts[:, 0]), "y": r3(pts[:, 1]), "z": r3(pts[:, 2]),
        "enclosure": [int(v) for v in cav.enclosure[sel]],
        "electrostatic": r3(cav.electrostatic[sel]),
        "donor": r3(cav.donor[sel]),
        "acceptor": r3(cav.acceptor[sel]),
        "hydrophobic": r3(cav.hydrophobic[sel]),
        "from_anchor": r3(cav.from_anchor[sel]),
        "anchor_lined_share": round(float((cav.from_anchor >= 0.6).mean()), 3),
        "partner_lined_share": round(float((cav.from_anchor <= 0.4).mean()), 3),
        "note": (
            "Electrostatics are formal charges on ionisable groups with a "
            "distance-dependent dielectric. That gets the acid-or-base "
            "character of a cavity right and nothing finer."
        ),
    }


def _heavy(atoms: list[Atom]) -> list[Atom]:
    return [a for a in atoms if a.element.upper() != "H"]


def _residue_min_distances(
    residues: list[Residue], target: np.ndarray
) -> np.ndarray:
    """Closest heavy-atom approach from each residue to a set of points."""
    from scipy.spatial import cKDTree

    if not len(target):
        return np.full(len(residues), 99.0)
    tree = cKDTree(target)
    out = np.empty(len(residues))
    for i, r in enumerate(residues):
        xyz = np.array([(a.x, a.y, a.z) for a in _heavy(r.atoms)])
        out[i] = tree.query(xyz, k=1)[0].min() if len(xyz) else 99.0
    return np.minimum(out, 99.0)


def build(
    residues: list[Residue],
    ligand,
    anchor_chain: str,
    partner_chain: str,
    anchor_name: str = "",
    partner_name: str = "",
    entry: str = "run",
    label: str = "",
) -> dict:
    """The viewer blob for one complex.

    `residues` holds both proteins, already positioned -- for a docking
    result that means the partner's residues carry the pose's coordinates,
    not the ones it was uploaded with.
    """
    from scipy.spatial import cKDTree

    anchor_res = [r for r in residues if r.chain == anchor_chain]
    partner_res = [r for r in residues if r.chain == partner_chain]
    if not anchor_res:
        raise ValueError("no residues on the anchor chain %r" % anchor_chain)

    heavy_lig = _heavy(ligand.atoms)
    lig_xyz_raw = np.array([(a.x, a.y, a.z) for a in heavy_lig])

    xyz_raw = np.array([
        (a.x, a.y, a.z) for r in residues
        for a in [r.representative()]
    ])
    centre = xyz_raw.mean(axis=0)
    xyz = xyz_raw - centre

    ca_raw = np.array([
        (a.x, a.y, a.z) for r in residues
        for a in [r.atom("CA") or r.representative()]
    ])
    ca = ca_raw - centre

    # --- who is at the interface, by geometry ----------------------------
    # The viewer's `iface` is 0, 1 or 2: not at the interface, at it, or at
    # it and close to the ligand. The benchmark blob gets the third state
    # from the proximal/distal split that the frustration work defines; here
    # it is simply distance to the ligand, which is what that split measures
    # anyway and is honest to compute without a potential.
    partner_xyz = np.array([
        (a.x, a.y, a.z) for r in partner_res for a in _heavy(r.atoms)
    ]) if partner_res else np.zeros((0, 3))
    anchor_xyz = np.array([
        (a.x, a.y, a.z) for r in anchor_res for a in _heavy(r.atoms)
    ])

    iface = np.zeros(len(residues), dtype=np.uint8)
    if len(partner_xyz):
        a_tree, p_tree = cKDTree(anchor_xyz), cKDTree(partner_xyz)
        for i, r in enumerate(residues):
            rx = np.array([(a.x, a.y, a.z) for a in _heavy(r.atoms)])
            if not len(rx):
                continue
            other = p_tree if r.chain == anchor_chain else a_tree
            if other.query(rx, k=1)[0].min() <= CONTACT:
                iface[i] = 1

    dlig = _residue_min_distances(residues, lig_xyz_raw)
    iface[(iface == 1) & (dlig <= 6.0)] = 2

    # --- sidechains, interface residues only -----------------------------
    side_res, side_x, side_y, side_z, side_name = [], [], [], [], []
    for ri, r in enumerate(residues):
        if not iface[ri]:
            continue
        for atom in r.atoms:
            if atom.element.upper() == "H" or atom.name in BACKBONE:
                continue
            side_res.append(int(ri))
            side_x.append(atom.x - centre[0])
            side_y.append(atom.y - centre[1])
            side_z.append(atom.z - centre[2])
            side_name.append(atom.name)

    # --- the ligand, and how much of it each protein covers ---------------
    anchor_atoms = [a for r in anchor_res for a in r.atoms]
    binary = anchor_atoms + list(ligand.atoms)
    ternary = [a for r in residues for a in r.atoms] + list(ligand.atoms)
    sasa_binary = atom_sasa(heavy_lig, binary)
    sasa_ternary = atom_sasa(heavy_lig, ternary)

    # --- the cavity at the interface --------------------------------------
    all_atoms = [a for r in residues for a in r.atoms]
    anchor_idx = {i for i, at in enumerate(all_atoms) if at.chain == anchor_chain}
    cav = detect_cavity(all_atoms, lig_xyz_raw.mean(axis=0), anchor_idx)

    # How much of the ligand each protein touches. Recomputed rather than
    # carried over, because for a docked pose the partner has moved and the
    # count it came with describes where it used to be.
    lig_tree = cKDTree(lig_xyz_raw)
    n_anchor_lig = int(sum(
        len(c) for c in lig_tree.query_ball_point(anchor_xyz, r=CONTACT)))
    n_partner_lig = int(sum(
        len(c) for c in lig_tree.query_ball_point(partner_xyz, r=CONTACT))
    ) if len(partner_xyz) else 0

    cavity = cavity_blob(cav, centre)

    return {
        "entry": entry,
        "label": label,
        # Said in the blob rather than inferred from a missing key, so the
        # viewer never has to guess why a style is unavailable.
        "drawing_only": True,
        # `residues.rho` is the local density the AWSEM potential defines,
        # not a geometric neighbour count, and the page's live frustration
        # engine would read it as the former. A plausible number in its
        # place would be wrong in a way nothing on screen could show, so it
        # is absent and declared absent.
        "missing": ["contacts", "exit_vectors", "result", "residues.rho"],
        "missing_note": (
            "Built from coordinates alone. The frustration analysis is not "
            "run for a docking result, so contact energies, local densities "
            "and exit vectors are not available, and the styles that draw "
            "them are not offered."
        ),
        # The same keys `Assignment.as_dict` produces, so a run and a
        # benchmark system are interchangeable to every page that reads one.
        # The contact counts are recomputed here from the coordinates given,
        # which for a docked pose is the only correct source: the partner is
        # somewhere new, so the counts it arrived with are stale.
        "assignment": {
            "entry": entry,
            "ligand": ligand.comp_id,
            "ligand_chain": getattr(ligand, "chain", ""),
            "ligand_seq": getattr(ligand, "seq_id", 0),
            "ligand_heavy_atoms": len(heavy_lig),
            "anchor_chain": anchor_chain,
            "anchor_name": anchor_name or ("chain " + anchor_chain),
            "partner_chain": partner_chain,
            "partner_name": partner_name or ("chain " + partner_chain),
            "anchor_ligand_contacts": n_anchor_lig,
            "partner_ligand_contacts": n_partner_lig,
            "anchor_share": round(
                n_anchor_lig / max(n_anchor_lig + n_partner_lig, 1), 4),
            "other_contacting_chains": {},
        },
        "residues": {
            "chain": [r.chain for r in residues],
            "seq": [r.seq_id for r in residues],
            "one": "".join(r.one for r in residues),
            "x": _r3(xyz[:, 0]), "y": _r3(xyz[:, 1]), "z": _r3(xyz[:, 2]),
            "cax": _r3(ca[:, 0]), "cay": _r3(ca[:, 1]), "caz": _r3(ca[:, 2]),
            "ss": secondary_structure(residues),
            "iface": iface.tolist(),
            "dlig": _r3(dlig),
        },
        "sidechains": {
            "note": "heavy sidechain atoms, interface residues only",
            "res": side_res,
            "x": _r3(side_x), "y": _r3(side_y), "z": _r3(side_z),
            "name": side_name,
        },
        "ligand": {
            "comp_id": ligand.comp_id,
            "x": _r3(lig_xyz_raw[:, 0] - centre[0]),
            "y": _r3(lig_xyz_raw[:, 1] - centre[1]),
            "z": _r3(lig_xyz_raw[:, 2] - centre[2]),
            "element": [a.element for a in heavy_lig],
            "name": [a.name for a in heavy_lig],
            "sasa": _r3(sasa_binary),
            "sasa_ternary": _r3(sasa_ternary),
            "sasa_total": round(float(sasa_binary.sum()), 1),
            "sasa_total_ternary": round(float(sasa_ternary.sum()), 1),
            "buried_by_partner": round(float(
                1.0 - sasa_ternary.sum() / max(sasa_binary.sum(), 1e-6)), 3),
        },
        "cavity": cavity,
    }
