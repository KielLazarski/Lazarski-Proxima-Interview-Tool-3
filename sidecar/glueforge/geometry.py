"""Solvent accessibility, ligand exit vectors, and what they can reach.

This is step 2 and step 6 of the workflow. Step 2 asks which part of the
ligand is actually presented to the solvent once it sits in the anchor's
pocket -- that exposed face is the only part that can ever touch a second
protein, and the only place chemistry can grow from. Step 6 asks, for each of
those growth points, which frustrated interface contacts are within reach.

Nothing here needs a chemistry toolkit. Exposure is Shrake-Rupley, exit
vectors are local outward normals, and reach is a distance-and-angle
feasibility test. RDKit would be needed to turn a chosen vector into an
actual molecule; it is not needed to decide which vector is worth growing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .cif import Atom, Ligand, Residue
from .residues import element_radius

PROBE_RADIUS = 1.4

#: Golden-spiral sphere points. Deterministic, unlike a random sampling, so
#: two runs give the same accessibility to the last decimal.
_SPHERE_CACHE: dict[int, np.ndarray] = {}


def sphere_points(n: int = 192) -> np.ndarray:
    if n in _SPHERE_CACHE:
        return _SPHERE_CACHE[n]
    k = np.arange(n, dtype=np.float64) + 0.5
    phi = np.arccos(1.0 - 2.0 * k / n)
    golden = np.pi * (1.0 + 5.0 ** 0.5)
    theta = golden * k
    pts = np.stack(
        [np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)],
        axis=1,
    )
    _SPHERE_CACHE[n] = pts
    return pts


def atom_sasa(
    target: list[Atom],
    context: list[Atom],
    n_points: int = 192,
    probe: float = PROBE_RADIUS,
    include_hydrogens: bool = False,
) -> np.ndarray:
    """Per-atom solvent-accessible surface of `target`, in the presence of
    `context` (which should include the target itself).

    Shrake-Rupley: roll a probe over each atom's expanded sphere and count the
    fraction of sample points no other atom occludes.

    **Hydrogens are discarded from the context by default, and this is not a
    detail.** Whether a deposition models hydrogens is a choice of the
    depositor, not a property of the molecule: 9 of our 14 benchmark
    structures have them and 5 do not. A ligand's own hydrogens shield its
    own heavy atoms heavily, so keeping them changes its apparent exposed
    area by 1.7x to 4.3x -- which means absolute areas computed with
    hydrogens present cannot be compared against areas computed without
    them. Discarding them always costs a little physical accuracy on the
    structures that have them and buys comparability across every structure,
    which is the trade this tool needs.

    The protein's hydrogens turn out not to matter either way (dropping them
    alone moves the answer by about 2%); it is the ligand's own that do the
    work. Relative measures -- buried fraction, per-atom relative accessibility
    -- are stable across the two conventions to within about 0.05, because the
    shielding largely cancels between the bound and free states. Prefer those
    when comparing one system against another.

    The returned array stays aligned with `target` as given: a hydrogen in
    `target` gets 0.0 rather than being dropped, so callers can index straight
    back into their own atom list.
    """
    if not target:
        return np.zeros(0)

    if include_hydrogens:
        keep_t = list(range(len(target)))
        ctx = context
    else:
        keep_t = [i for i, a in enumerate(target) if a.element.upper() != "H"]
        ctx = [a for a in context if a.element.upper() != "H"]
        if not keep_t or not ctx:
            return np.zeros(len(target))

    sub = [target[i] for i in keep_t]

    tx = np.array([(a.x, a.y, a.z) for a in sub], dtype=np.float64)
    tr = np.array([element_radius(a.element) for a in sub], dtype=np.float64) + probe

    cx = np.array([(a.x, a.y, a.z) for a in ctx], dtype=np.float64)
    cr = np.array([element_radius(a.element) for a in ctx], dtype=np.float64) + probe

    pts = sphere_points(n_points)
    out = np.zeros(len(target), dtype=np.float64)
    max_reach = float(cr.max() + tr.max())

    for k, i in enumerate(keep_t):
        centre, radius = tx[k], tr[k]
        d = np.linalg.norm(cx - centre, axis=1)
        near = np.nonzero((d < max_reach) & (d > 1e-6))[0]
        if near.size == 0:
            out[i] = 4.0 * np.pi * radius * radius
            continue
        probes = centre + pts * radius
        gap = np.linalg.norm(
            probes[:, None, :] - cx[None, near, :], axis=2
        ) - cr[None, near]
        accessible = (gap > 0).all(axis=1)
        out[i] = (
            4.0 * np.pi * radius * radius * accessible.sum() / float(n_points)
        )
    return out


@dataclass
class ExitVector:
    """One place the ligand could be grown from, and which way it points."""

    atom_index: int
    atom_name: str
    element: str
    position: np.ndarray
    direction: np.ndarray
    sasa: float
    relative_sasa: float
    #: Frustrated interface contacts this vector could plausibly address.
    reachable: list[int]
    addressable: float
    mean_reach_distance: float

    def as_dict(self) -> dict:
        return {
            "atom_index": int(self.atom_index),
            "atom_name": self.atom_name,
            "element": self.element,
            "position": [round(float(v), 3) for v in self.position],
            "direction": [round(float(v), 4) for v in self.direction],
            "sasa": round(float(self.sasa), 2),
            "relative_sasa": round(float(self.relative_sasa), 3),
            "n_reachable": len(self.reachable),
            "reachable": [int(k) for k in self.reachable],
            "addressable": round(float(self.addressable), 3),
            # null, not NaN: bare NaN is not valid JSON and the dashboard
            # parses this file.
            "mean_reach_distance": (
                round(float(self.mean_reach_distance), 2)
                if math.isfinite(self.mean_reach_distance) else None
            ),
        }


def exit_vectors(
    ligand: Ligand,
    context_atoms: list[Atom],
    min_relative_sasa: float = 0.08,
) -> list[tuple[int, np.ndarray, np.ndarray, float, float]]:
    """Growth points on the ligand: exposed heavy atoms and their outward normals.

    The direction is the local outward normal, taken as the vector from the
    mean of an atom's bonded neighbours to the atom itself. That is where a
    substituent would project, which is what matters for growing a molecule.
    Atoms with no neighbour within bonding distance fall back to the vector
    from the ligand centroid.
    """
    heavy = [a for a in ligand.atoms if a.element.upper() != "H"]
    if not heavy:
        return []

    sasa = atom_sasa(heavy, context_atoms)
    free = atom_sasa(heavy, heavy)
    rel = np.divide(sasa, np.maximum(free, 1e-6))

    pos = np.array([(a.x, a.y, a.z) for a in heavy], dtype=np.float64)
    centroid = pos.mean(axis=0)

    out = []
    for i, atom in enumerate(heavy):
        if rel[i] < min_relative_sasa:
            continue
        d = np.linalg.norm(pos - pos[i], axis=1)
        bonded = np.nonzero((d > 0.4) & (d < 1.95))[0]
        if bonded.size:
            direction = pos[i] - pos[bonded].mean(axis=0)
        else:
            direction = pos[i] - centroid
        norm = np.linalg.norm(direction)
        if norm < 1e-6:
            continue
        out.append((i, pos[i], direction / norm, float(sasa[i]), float(rel[i])))
    return out


def reach_frustrated_contacts(
    vectors,
    heavy_atoms: list[Atom],
    contact_midpoints: np.ndarray,
    contact_frustration: np.ndarray,
    contact_mask: np.ndarray,
    max_reach: float = 12.0,
    max_angle_deg: float = 75.0,
    frustrated_below: float = -0.3,
) -> list[ExitVector]:
    """For each growth point, which frustrated interface contacts it could address.

    A contact is reachable when it sits within `max_reach` angstroms of the
    growth point and within `max_angle_deg` of the exit direction. The
    distance bound is roughly what three to five rotatable bonds span; the
    angle bound keeps the chemistry from having to fold back over the ligand.

    `addressable` sums the frustration it could reach, discounted by how far
    and how far off-axis each contact is, so a vector pointing straight at a
    cluster of frustrated contacts scores above one that merely has them
    somewhere in front of it.
    """
    out: list[ExitVector] = []
    cos_limit = np.cos(np.deg2rad(max_angle_deg))

    candidates = np.nonzero(
        contact_mask & (contact_frustration < frustrated_below)
    )[0]
    if candidates.size == 0:
        candidates = np.nonzero(contact_mask)[0]

    mids = contact_midpoints[candidates]
    fr = contact_frustration[candidates]

    for idx, position, direction, sasa, rel in vectors:
        delta = mids - position
        dist = np.linalg.norm(delta, axis=1)
        ok = (dist > 1e-6) & (dist <= max_reach)
        if not ok.any():
            out.append(ExitVector(
                atom_index=idx, atom_name=heavy_atoms[idx].name,
                element=heavy_atoms[idx].element, position=position,
                direction=direction, sasa=sasa, relative_sasa=rel,
                reachable=[], addressable=0.0, mean_reach_distance=float("nan"),
            ))
            continue

        unit = delta[ok] / dist[ok][:, None]
        cos = unit @ direction
        within = cos >= cos_limit
        sel = np.nonzero(ok)[0][within]

        if sel.size == 0:
            out.append(ExitVector(
                atom_index=idx, atom_name=heavy_atoms[idx].name,
                element=heavy_atoms[idx].element, position=position,
                direction=direction, sasa=sasa, relative_sasa=rel,
                reachable=[], addressable=0.0, mean_reach_distance=float("nan"),
            ))
            continue

        d_sel = dist[sel]
        c_sel = (unit[within] @ direction)
        # Linear falloff with distance, and with how far off-axis it sits.
        weight = (1.0 - d_sel / max_reach) * np.clip(
            (c_sel - cos_limit) / max(1.0 - cos_limit, 1e-6), 0.0, 1.0
        )
        score = float(np.sum(np.abs(np.minimum(fr[sel], 0.0)) * weight))

        out.append(ExitVector(
            atom_index=idx,
            atom_name=heavy_atoms[idx].name,
            element=heavy_atoms[idx].element,
            position=position,
            direction=direction,
            sasa=sasa,
            relative_sasa=rel,
            reachable=[int(candidates[k]) for k in sel],
            addressable=score,
            mean_reach_distance=float(d_sel.mean()),
        ))

    out.sort(key=lambda v: -v.addressable)
    return out


def contact_midpoints(
    residues: list[Residue], index_i: np.ndarray, index_j: np.ndarray
) -> np.ndarray:
    """Midpoint between the representative atoms of each contact pair."""
    xyz = np.array(
        [[r.representative().x, r.representative().y, r.representative().z]
         for r in residues],
        dtype=np.float64,
    )
    return 0.5 * (xyz[index_i] + xyz[index_j])


def secondary_structure(residues: list[Residue]) -> str:
    """Rough helix / strand / coil assignment from CA geometry alone.

    A full DSSP needs backbone hydrogen bonding and therefore N, C and O.
    This uses the i to i+3 and i to i+4 CA distances, which separate helix
    from extended strand cleanly enough to draw a ribbon that reads correctly.
    It is for display only and is labelled as such wherever it surfaces.
    """
    ca = []
    for r in residues:
        a = r.atom("CA") or r.representative()
        ca.append((a.x, a.y, a.z))
    pos = np.array(ca, dtype=np.float64)
    n = len(residues)
    out = ["C"] * n
    if n < 5:
        return "".join(out)

    chains = [r.chain for r in residues]
    for i in range(n - 4):
        if chains[i] != chains[i + 4]:
            continue
        d3 = np.linalg.norm(pos[i + 3] - pos[i])
        d4 = np.linalg.norm(pos[i + 4] - pos[i])
        if 4.6 <= d3 <= 6.6 and 5.2 <= d4 <= 7.4:
            for k in range(i, i + 4):
                out[k] = "H"
        elif d3 >= 8.8 and d4 >= 11.5:
            for k in range(i, i + 3):
                if out[k] == "C":
                    out[k] = "E"
    return "".join(out)
