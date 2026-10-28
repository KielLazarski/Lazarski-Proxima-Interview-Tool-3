"""The presented face: the part of the ligand a second protein could touch.

This is the single idea the whole tool rests on. A ligand sitting in its
anchor's pocket is mostly buried. Whatever sliver of it still sees solvent is
the only surface a partner protein can ever reach, so it is both the thing a
glue has to be grown from and the thing a docking run has to be aimed at.

The measurement that matters is deliberately narrow: ligand accessibility
**against the anchor chain alone**. Computing it in the assembled ternary
complex instead is the obvious mistake and it is self-defeating -- in a real
complex the partner has already covered the face, so every atom reads as
buried and there is nothing to grow from and nothing to aim at. Against the
anchor alone, CR8 in DDB1-CDK12 presents about 59 square angstroms from a
ligand that is 97 percent buried, and that sliver is exactly where cyclin K
binds. The number is small by construction; it being small is the finding,
not a failure of the measurement.

Everything downstream reads from here. The rim residues this returns become
the LightDock restraint in phase 2, and the exposed atoms become the growth
points in phase 4.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .cif import Atom, Ligand, Residue, Structure
from .geometry import atom_sasa
from .residues import NON_LIGAND

#: An atom keeping this share of its unhindered surface counts as presented.
#: Low on purpose: the whole face is often only a handful of atoms, and at
#: 0.15 several published glues present nothing at all.
MIN_RELATIVE_SASA = 0.08

#: How far an anchor residue may sit from an exposed ligand atom and still be
#: said to line the face. 8 angstroms is about one residue's reach past van
#: der Waals contact, which is the distance across which a partner residue
#: packing onto the ligand would also touch this one.
RIM_CUTOFF = 8.0

#: Below this, the exposed face wraps the ligand rather than pointing one
#: way, its normal is arbitrary, and anything ranked by that normal is noise.
#:
#: Placed in an observed gap: the three groove-bound ligands in the benchmark
#: score 0.215-0.274 and the eleven pocket-bound ones 0.343-0.877. That is a
#: clean separation but it is n=14, so treat this as a soft warning rather
#: than a law -- what it changes is whether active restraints are set at all,
#: not whether the pair is worth looking at.
MIN_DIRECTIONALITY = 0.30

#: An anchor residue has to be solvent-facing itself to be part of the rim.
#: Without this the set fills up with residues buried *behind* the ligand,
#: which a partner protein can never reach and which would pull a docking
#: run into the anchor's core.
MIN_RIM_SASA = 5.0


def _heavy(atoms) -> list[Atom]:
    return [a for a in atoms if a.element.upper() != "H"]


def _coords(atoms: list[Atom]) -> np.ndarray:
    if not atoms:
        return np.zeros((0, 3))
    return np.array([(a.x, a.y, a.z) for a in atoms], dtype=np.float64)


@dataclass
class RimResidue:
    chain: str
    comp_id: str
    seq_id: int
    insertion: str
    one: str
    min_distance: float
    sasa: float
    #: How far this residue sits along the face's outward normal, from the
    #: face centroid. Positive is on the side a partner approaches from.
    #:
    #: This, not `min_distance`, is what predicts whether the real partner
    #: touches a rim residue. Measured across the 14 benchmark systems,
    #: distance to the nearest exposed ligand atom carries no signal at all
    #: (r = -0.07, p = 0.22 over 324 rim residues) while normal projection
    #: carries a lot (r = +0.57, p < 1e-4). The reason is geometric: a
    #: partner comes in along the normal, so two residues equally close to
    #: the ligand are not equally reachable -- one is on the approach side
    #: and the other is behind it.
    normal_projection: float = 0.0

    def as_dict(self) -> dict:
        return {
            "chain": self.chain,
            "comp_id": self.comp_id,
            "seq_id": int(self.seq_id),
            "insertion": self.insertion,
            "one": self.one,
            "min_distance": round(float(self.min_distance), 2),
            "sasa": round(float(self.sasa), 1),
            "normal_projection": round(float(self.normal_projection), 2),
            "label": self.label,
        }

    @property
    def label(self) -> str:
        return "%s.%s.%d%s" % (
            self.chain, self.comp_id, self.seq_id, self.insertion.strip()
        )


@dataclass
class ExposedAtom:
    index: int
    name: str
    element: str
    position: np.ndarray
    sasa: float
    relative_sasa: float

    def as_dict(self) -> dict:
        return {
            "index": int(self.index),
            "name": self.name,
            "element": self.element,
            "position": [round(float(v), 3) for v in self.position],
            "sasa": round(float(self.sasa), 2),
            "relative_sasa": round(float(self.relative_sasa), 3),
        }


@dataclass
class PresentedFace:
    """What protein A offers the solvent at its ligand site."""

    ligand_comp: str
    ligand_chain: str
    ligand_seq: int
    anchor_chains: list[str]
    ligand_heavy_atoms: int

    #: Ligand surface with nothing else present.
    free_sasa: float
    #: Ligand surface once the anchor is there. The difference is the pocket.
    bound_sasa: float
    #: Sum over atoms that clear MIN_RELATIVE_SASA -- the usable face.
    presented_area: float

    exposed: list[ExposedAtom]
    rim: list[RimResidue]

    centroid: np.ndarray
    normal: np.ndarray

    #: How much the exposed face points one way, 0 to 1. The area-weighted
    #: mean of the unit vectors from the ligand centroid to each exposed
    #: atom: near 1 and the face is a single opening, near 0 and it wraps the
    #: ligand. Everything that uses `normal` is only as meaningful as this.
    directionality: float = 0.0

    warnings: list[str] = field(default_factory=list)

    @property
    def directional(self) -> bool:
        """Is the normal worth ranking anything by?"""
        return self.directionality >= MIN_DIRECTIONALITY

    @property
    def buried_fraction(self) -> float:
        if self.free_sasa <= 0:
            return 0.0
        return 1.0 - (self.bound_sasa / self.free_sasa)

    def as_dict(self) -> dict:
        return {
            "ligand": self.ligand_comp,
            "ligand_chain": self.ligand_chain,
            "ligand_seq": int(self.ligand_seq),
            "anchor_chains": list(self.anchor_chains),
            "ligand_heavy_atoms": int(self.ligand_heavy_atoms),
            "free_sasa": round(float(self.free_sasa), 1),
            "bound_sasa": round(float(self.bound_sasa), 1),
            "buried_fraction": round(float(self.buried_fraction), 4),
            "presented_area": round(float(self.presented_area), 1),
            "n_exposed": len(self.exposed),
            "exposed": [a.as_dict() for a in self.exposed],
            "rim": [r.as_dict() for r in self.rim],
            "n_rim": len(self.rim),
            "centroid": [round(float(v), 3) for v in self.centroid],
            "normal": [round(float(v), 4) for v in self.normal],
            "directionality": round(float(self.directionality), 3),
            "directional": bool(self.directional),
            "warnings": list(self.warnings),
        }

    def restraint_labels(self, limit: int = 40) -> list[str]:
        """The rim as LightDock wants it: `chain.RES.seq`, nearest first."""
        return [r.label for r in self.rim[:limit]]


def default_anchor_chains(
    structure: Structure, ligand: Ligand, cutoff: float = 4.5
) -> list[str]:
    """Chains in actual contact with the ligand, most contacts first.

    When a file holds the whole ternary complex this returns both partners,
    and the caller has to say which one is the anchor. When it holds protein
    A alone -- the normal case here -- there is only one answer.
    """
    lig = _coords(_heavy(ligand.atoms))
    if not len(lig):
        return []
    counts: dict[str, int] = {}
    for res in structure.residues:
        if res.comp_id in NON_LIGAND:
            continue
        xyz = _coords(_heavy(res.atoms))
        if not len(xyz):
            continue
        d = np.linalg.norm(xyz[:, None, :] - lig[None, :, :], axis=2)
        n = int((d < cutoff).sum())
        if n:
            counts[res.chain] = counts.get(res.chain, 0) + n
    return [c for c, _ in sorted(counts.items(), key=lambda kv: -kv[1])]


def presented_face(
    structure: Structure,
    ligand: Ligand,
    anchor_chains: list[str] | None = None,
    min_relative_sasa: float = MIN_RELATIVE_SASA,
    rim_cutoff: float = RIM_CUTOFF,
    min_rim_sasa: float = MIN_RIM_SASA,
) -> PresentedFace:
    """Measure the ligand's exposed face against the anchor, and who lines it.

    `anchor_chains` defaults to every chain touching the ligand. Pass it
    explicitly when the file holds a complete ternary complex and you want the
    face as protein A alone would present it -- which is the only way the
    answer means anything, since the partner would otherwise be sitting on
    the very surface being measured.
    """
    warnings: list[str] = []

    lig_atoms = _heavy(ligand.atoms)
    if not lig_atoms:
        raise ValueError("ligand %s has no heavy atoms" % ligand.comp_id)

    if anchor_chains is None:
        anchor_chains = default_anchor_chains(structure, ligand)
        if len(anchor_chains) > 1:
            warnings.append(
                "The ligand touches %d chains (%s). The face was measured "
                "against all of them, so it reads as the face in the "
                "assembled complex, not the face protein A presents on its "
                "own. Name the anchor chain to measure it properly."
                % (len(anchor_chains), ", ".join(anchor_chains))
            )
    if not anchor_chains:
        raise ValueError(
            "ligand %s is not in contact with any protein chain"
            % ligand.comp_id
        )

    wanted = set(anchor_chains)
    anchor_res = [
        r for r in structure.residues
        if r.chain in wanted and r.comp_id not in NON_LIGAND
    ]
    if not anchor_res:
        raise ValueError(
            "no protein residues in chain(s) %s" % ", ".join(sorted(wanted))
        )
    anchor_atoms = _heavy([a for r in anchor_res for a in r.atoms])

    # -- the face itself ---------------------------------------------------
    # Context is the anchor plus the ligand: the ligand occludes itself, and
    # leaving it out of its own context inflates every atom's accessibility.
    free = atom_sasa(lig_atoms, lig_atoms)
    bound = atom_sasa(lig_atoms, anchor_atoms + lig_atoms)
    rel = np.divide(bound, np.maximum(free, 1e-6))

    pos = _coords(lig_atoms)
    keep = np.nonzero(rel >= min_relative_sasa)[0]

    exposed = [
        ExposedAtom(
            index=int(i),
            name=lig_atoms[i].name,
            element=lig_atoms[i].element,
            position=pos[i],
            sasa=float(bound[i]),
            relative_sasa=float(rel[i]),
        )
        for i in keep
    ]
    exposed.sort(key=lambda a: -a.sasa)
    presented_area = float(bound[keep].sum()) if keep.size else 0.0

    if not exposed:
        warnings.append(
            "No part of the ligand is exposed once the anchor is in place. "
            "A fully enclosed ligand cannot be reached by a second protein "
            "and cannot be grown from, so this pair has no glue to make."
        )

    # -- where it faces ----------------------------------------------------
    # Area-weighted so one large opening is not outvoted by several pinholes.
    directionality = 0.0
    if exposed:
        lig_centre = pos.mean(axis=0)
        w = np.array([a.sasa for a in exposed], dtype=np.float64)
        vec = np.array([a.position for a in exposed]) - lig_centre
        mag = np.linalg.norm(vec, axis=1, keepdims=True)
        unit = vec / np.maximum(mag, 1e-9)
        directionality = float(np.linalg.norm(
            (unit * w[:, None]).sum(axis=0) / w.sum()
        ))

    if exposed:
        w = np.array([a.sasa for a in exposed], dtype=np.float64)
        pts = np.array([a.position for a in exposed], dtype=np.float64)
        centroid = (pts * w[:, None]).sum(axis=0) / w.sum()
        direction = centroid - pos.mean(axis=0)
        norm = float(np.linalg.norm(direction))
        if norm < 1e-3:
            # The face wraps the ligand rather than pointing anywhere -- the
            # ligand is sitting in a groove, not a pocket. Fall back to the
            # direction away from the anchor's bulk.
            away = centroid - _coords(anchor_atoms).mean(axis=0)
            n2 = float(np.linalg.norm(away))
            normal = away / n2 if n2 > 1e-6 else np.array([0.0, 0.0, 1.0])
            warnings.append(
                "The exposed atoms wrap around the ligand rather than facing "
                "one way. The ligand is probably in a surface groove, which "
                "makes the docking restraint much less directional."
            )
        else:
            normal = direction / norm
    else:
        centroid = pos.mean(axis=0)
        normal = np.array([0.0, 0.0, 1.0])

    # -- who lines it ------------------------------------------------------
    rim = _rim_residues(
        anchor_res, anchor_atoms, exposed, rim_cutoff, min_rim_sasa,
        centroid, normal,
    )
    if exposed and directionality < MIN_DIRECTIONALITY:
        warnings.append(
            "The exposed face wraps the ligand rather than pointing one way "
            "(directionality %.2f, below %.2f). This ligand is in a surface "
            "groove, not a pocket. Docking will still be restrained to the "
            "region, but none of the residues can be required, because "
            "there is no approach direction to pick them by."
            % (directionality, MIN_DIRECTIONALITY)
        )

    if exposed and not rim:
        warnings.append(
            "The exposed face has no solvent-facing anchor residues around "
            "it, so there is nothing to restrain a docking run with."
        )

    return PresentedFace(
        ligand_comp=ligand.comp_id,
        ligand_chain=ligand.chain,
        ligand_seq=ligand.seq_id,
        anchor_chains=list(anchor_chains),
        ligand_heavy_atoms=len(lig_atoms),
        free_sasa=float(free.sum()),
        bound_sasa=float(bound.sum()),
        presented_area=presented_area,
        exposed=exposed,
        rim=rim,
        centroid=centroid,
        normal=normal,
        directionality=directionality,
        warnings=warnings,
    )


def _rim_residues(
    anchor_res: list[Residue],
    anchor_atoms: list[Atom],
    exposed: list[ExposedAtom],
    rim_cutoff: float,
    min_rim_sasa: float,
    face_centroid: np.ndarray,
    face_normal: np.ndarray,
) -> list[RimResidue]:
    """Anchor residues that both surround the exposed face and face out.

    Both halves of that are load-bearing. Proximity alone picks up the pocket
    floor underneath the ligand, which a partner protein can never touch;
    exposure alone picks up the whole surface of the protein.
    """
    if not exposed:
        return []

    face = np.array([a.position for a in exposed], dtype=np.float64)

    # One SASA pass over the whole anchor, then summed per residue.
    sasa = atom_sasa(anchor_atoms, anchor_atoms)
    by_atom = {id(a): float(s) for a, s in zip(anchor_atoms, sasa)}

    out: list[RimResidue] = []
    for res in anchor_res:
        heavy = _heavy(res.atoms)
        if not heavy:
            continue
        xyz = _coords(heavy)
        d = np.linalg.norm(xyz[:, None, :] - face[None, :, :], axis=2)
        dmin = float(d.min())
        if dmin > rim_cutoff:
            continue
        res_sasa = sum(by_atom.get(id(a), 0.0) for a in heavy)
        if res_sasa < min_rim_sasa:
            continue
        out.append(RimResidue(
            chain=res.chain,
            comp_id=res.comp_id,
            seq_id=res.seq_id,
            insertion=res.insertion,
            one=res.one,
            min_distance=dmin,
            sasa=res_sasa,
            normal_projection=float(
                np.dot(xyz.mean(axis=0) - face_centroid, face_normal)
            ),
        ))

    out.sort(key=lambda r: r.min_distance)
    return out
