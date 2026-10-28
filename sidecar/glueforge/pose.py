"""Rigid-body poses, and how to tell a good one from a bad one.

A docking result is a placement of protein B relative to protein A. This
module is everything needed to represent one, measure how far it is from
another, and score it -- in particular on the thing this tool cares about and
general docking does not: **how much of the ligand's presented face the
partner buries.**

That is the scoring decision worth stating plainly. An ordinary docking score
asks whether the interface is favourable. A glue interface is often not
especially favourable; what makes it a glue interface is that the ligand is in
it. A pose that buries the whole presented face is a candidate however
ordinary its interface energy, and a pose with a beautiful interface on the
far side of the protein is not a candidate at all. So face burial is scored
separately and reported next to the energy rather than folded into it, because
collapsing them into one number hides exactly the disagreement you want to
see.

Nothing here needs LightDock. It scores poses from any source: a docking run,
a crystal structure, or the random placements used to ask whether the
restraint carries any signal in the first place.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .cif import Atom
from .face import PresentedFace
from .geometry import atom_sasa

#: Heavy-atom contact distance. The conventional 5 A for interface definition.
CONTACT = 5.0

#: Below this the two proteins are interpenetrating, not docked. Measured as
#: the count of heavy-atom pairs closer than this.
CLASH = 2.5


def heavy(atoms: list[Atom]) -> list[Atom]:
    return [a for a in atoms if a.element.upper() != "H"]


def coords(atoms: list[Atom]) -> np.ndarray:
    if not atoms:
        return np.zeros((0, 3))
    return np.array([(a.x, a.y, a.z) for a in atoms], dtype=np.float64)


def rotation_from_axis_angle(axis: np.ndarray, radians: float) -> np.ndarray:
    """Rodrigues. Returns a 3x3 rotation matrix."""
    axis = np.asarray(axis, dtype=np.float64)
    n = float(np.linalg.norm(axis))
    if n < 1e-12:
        return np.eye(3)
    axis = axis / n
    K = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    return (np.eye(3) + np.sin(radians) * K
            + (1.0 - np.cos(radians)) * (K @ K))


def random_rotation(rng: np.random.Generator) -> np.ndarray:
    """Uniform over SO(3).

    Via a normalised random quaternion rather than three random Euler angles,
    which concentrate near the poles and would make the decoy set quietly
    non-uniform -- and a non-uniform decoy set is how a restraint gets
    credited with signal it does not have.
    """
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
    ])


@dataclass
class Pose:
    """A rigid placement of the partner, as a rotation about its own centroid
    followed by a translation. Identity is the pose as given."""

    rotation: np.ndarray = field(default_factory=lambda: np.eye(3))
    translation: np.ndarray = field(default_factory=lambda: np.zeros(3))
    #: Where the rotation is applied about; the partner's own centroid.
    centre: np.ndarray | None = None
    label: str = ""

    def apply(self, xyz: np.ndarray) -> np.ndarray:
        c = self.centre if self.centre is not None else xyz.mean(axis=0)
        return (xyz - c) @ self.rotation.T + c + self.translation

    def apply_atoms(self, atoms: list[Atom]) -> list[Atom]:
        moved = self.apply(coords(atoms))
        return [
            Atom(a.serial, a.name, a.element, a.comp_id, a.chain, a.seq_id,
                 float(p[0]), float(p[1]), float(p[2]), a.occupancy,
                 a.hetatm, a.insertion)
            for a, p in zip(atoms, moved)
        ]

    @staticmethod
    def identity(atoms: list[Atom]) -> "Pose":
        return Pose(centre=coords(atoms).mean(axis=0), label="crystal")


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """Positional RMSD with no superposition.

    Deliberately not superposed: the question is how far this placement is
    from the right one, and superposing would remove exactly the error being
    measured. This is "ligand RMSD" in docking terms, with protein B playing
    the ligand.
    """
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch: {a.shape} vs {b.shape}")
    if not len(a):
        return float("nan")
    d = a - b
    return float(np.sqrt((d * d).sum() / len(a)))


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


@dataclass
class PoseScore:
    """Everything measured about one placement.

    `face_burial` is the number this tool exists to compute; the rest are
    there so that a high burial achieved by driving the two proteins through
    each other is visibly not a result.
    """

    label: str
    face_burial: float          # 0-1, share of the presented face covered
    buried_area: float          # A^2 of the face covered
    contacts: int               # heavy-atom pairs within CONTACT
    clashes: int                # heavy-atom pairs within CLASH
    interface_residues: int
    rmsd_to_reference: float | None = None
    #: Partner atoms driven into the ligand itself.
    ligand_clashes: int = 0
    #: Share of the ACTIVE restraint residues this pose actually contacts.
    #: This is what `lgd_filter_restraints.py --fnat` tests, so it measures
    #: the restraint directly rather than through face burial as a proxy.
    restraint_satisfaction: float = 0.0

    @property
    def viable(self) -> bool:
        """Touching, not interpenetrating, and not sitting on the ligand.

        The protein clash budget is not zero because rigid-body docking of two
        separately-solved structures always leaves a few: side chains that
        would move on binding. More than a handful means the backbones
        overlap.

        The ligand budget is tighter and it matters more than it looks.
        LightDock docks protein to protein; the ligand is not in the receptor
        file it is given, so nothing in the docking run stops a pose from
        settling into the ligand's own pocket. Such a pose buries the
        presented face perfectly -- by occupying it -- and would otherwise
        score as the best result in the run.
        """
        return (self.contacts >= 20 and self.clashes <= 15
                and self.ligand_clashes <= 2)

    def as_dict(self) -> dict:
        return {
            "label": self.label,
            "face_burial": round(float(self.face_burial), 4),
            "buried_area": round(float(self.buried_area), 1),
            "contacts": int(self.contacts),
            "clashes": int(self.clashes),
            "interface_residues": int(self.interface_residues),
            "ligand_clashes": int(self.ligand_clashes),
            "restraint_satisfaction": round(float(self.restraint_satisfaction), 4),
            "viable": bool(self.viable),
            "rmsd_to_reference": (
                round(float(self.rmsd_to_reference), 2)
                if self.rmsd_to_reference is not None else None
            ),
        }


class PoseScorer:
    """Scores many poses of one partner against one fixed anchor.

    Built once per system because the expensive parts -- the anchor's atoms,
    the ligand's free surface, the exposed-atom subset -- do not change
    between poses, and a decoy run scores thousands.
    """

    def __init__(
        self,
        anchor_atoms: list[Atom],
        ligand_atoms: list[Atom],
        face: PresentedFace,
        contact: float = CONTACT,
        clash: float = CLASH,
        restraint_labels: list[str] | None = None,
    ) -> None:
        self.anchor = heavy(anchor_atoms)
        self.ligand = heavy(ligand_atoms)
        self.face = face
        self.contact = contact
        self.clash = clash

        self.anchor_xyz = coords(self.anchor)
        self.ligand_xyz = coords(self.ligand)

        from scipy.spatial import cKDTree
        self._anchor_tree = cKDTree(self.anchor_xyz)
        self._ligand_tree = cKDTree(self.ligand_xyz)

        # Only the exposed atoms are scored, not the whole ligand. A decoy
        # scan runs thousands of poses and accessibility is the expensive
        # part, so scoring 11 atoms instead of 32 is most of the difference
        # between a minute and twenty. `face.exposed[i].index` indexes the
        # heavy-atom list, which is what `self.ligand` is.
        self.face_atoms = [self.ligand[a.index] for a in face.exposed]
        self.base_face_area = (
            float(atom_sasa(self.face_atoms, self.anchor + self.ligand).sum())
            if self.face_atoms else 0.0
        )
        self.face_centre = (
            coords(self.face_atoms).mean(axis=0) if self.face_atoms
            else self.ligand_xyz.mean(axis=0)
        )

        # The restraint residues, as atom coordinates grouped per residue, so
        # a pose can be asked the same question LightDock's own filter asks:
        # how many of these did you actually land on.
        if restraint_labels is None:
            from .restraints import restraints_from_face
            restraint_labels = restraints_from_face(face).active
        wanted = set(restraint_labels)
        groups: dict[str, list[int]] = {}
        for i, a in enumerate(self.anchor):
            key = "%s.%s.%d" % (a.chain, a.comp_id, a.seq_id)
            if key in wanted:
                groups.setdefault(key, []).append(i)
        self._restraint_groups = list(groups.values())
        self.n_restraints = len(self._restraint_groups)

    def score(
        self,
        partner_atoms: list[Atom],
        label: str = "",
        reference_xyz: np.ndarray | None = None,
    ) -> PoseScore:
        part = heavy(partner_atoms)
        pxyz = coords(part)
        ref_rmsd = (rmsd(pxyz, reference_xyz)
                    if reference_xyz is not None else None)

        if not self.face_atoms or not len(pxyz):
            return PoseScore(label, 0.0, 0.0, 0, 0, 0, ref_rmsd)

        lig_clash = sum(len(c) for c in self._ligand_tree.query_ball_point(
            pxyz, r=self.clash))

        satisfied = 0.0
        if self._restraint_groups:
            from scipy.spatial import cKDTree
            ptree = cKDTree(pxyz)
            hit = 0
            for idx in self._restraint_groups:
                if any(ptree.query_ball_point(self.anchor_xyz[i],
                                              r=self.contact)
                       for i in idx):
                    hit += 1
            satisfied = hit / len(self._restraint_groups)

        # Contacts and clashes, from the anchor tree built once in __init__.
        pairs = self._anchor_tree.query_ball_point(pxyz, r=self.contact)
        contacts = sum(len(p) for p in pairs)
        clashes = sum(
            len(c) for c in self._anchor_tree.query_ball_point(
                pxyz, r=self.clash)
        )
        iface = len({
            part[i].chain + str(part[i].seq_id)
            for i, lst in enumerate(pairs) if lst
        })

        # Face burial: re-run accessibility with the partner present and see
        # how much of the exposed face it took away. Only partner atoms that
        # could possibly reach the face are in the context -- a decoy on the
        # far side of the anchor then costs nothing at all.
        reach = 20.0
        near = pxyz[np.linalg.norm(pxyz - self.face_centre, axis=1) < reach]
        if not len(near):
            return PoseScore(label, 0.0, 0.0, int(contacts), int(clashes),
                             int(iface), ref_rmsd, int(lig_clash), satisfied)
        near_atoms = [
            a for a, p in zip(part, pxyz)
            if np.linalg.norm(p - self.face_centre) < reach
        ]

        remaining = float(atom_sasa(
            self.face_atoms, self.anchor + self.ligand + near_atoms
        ).sum())
        buried = max(0.0, self.base_face_area - remaining)
        burial = buried / self.base_face_area if self.base_face_area else 0.0

        return PoseScore(
            label=label,
            face_burial=float(burial),
            buried_area=buried,
            contacts=int(contacts),
            clashes=int(clashes),
            interface_residues=int(iface),
            rmsd_to_reference=ref_rmsd,
            ligand_clashes=int(lig_clash),
            restraint_satisfaction=float(satisfied),
        )


# ---------------------------------------------------------------------------
# Decoys
# ---------------------------------------------------------------------------


def random_surface_poses(
    anchor_atoms: list[Atom],
    partner_atoms: list[Atom],
    n: int,
    seed: int = 0,
    gap: float = 3.0,
) -> list[Pose]:
    """Random placements of the partner on the anchor's surface.

    Each is a uniform random orientation plus an approach along a uniform
    random direction, pushed in until the two just touch -- `gap` being the
    closest heavy-atom centres are allowed to come. 3 A is van der Waals
    contact; pressing in to 1.5 would interpenetrate, and would hand each
    decoy more contact and more face burial than a real complex could have. That makes every
    decoy a *plausible* complex -- touching, not interpenetrating, not
    floating in space -- which is the only kind worth comparing against. A
    decoy set of obviously-wrong placements would make any score look
    discriminating.
    """
    rng = np.random.default_rng(seed)
    a_xyz = coords(heavy(anchor_atoms))
    p_atoms = heavy(partner_atoms)
    p_xyz = coords(p_atoms)
    a_centre = a_xyz.mean(axis=0)
    p_centre = p_xyz.mean(axis=0)

    from scipy.spatial import cKDTree
    a_tree = cKDTree(a_xyz)

    # Radii bound the search: start far enough out that nothing touches.
    a_reach = float(np.linalg.norm(a_xyz - a_centre, axis=1).max())
    p_reach = float(np.linalg.norm(p_xyz - p_centre, axis=1).max())
    start = a_reach + p_reach + 5.0

    out: list[Pose] = []
    for i in range(n):
        R = random_rotation(rng)
        rotated = (p_xyz - p_centre) @ R.T + p_centre

        d = rng.normal(size=3)
        d /= np.linalg.norm(d)

        # Binary search along the approach for the closest non-clashing offset.
        lo, hi = 0.0, start
        for _ in range(28):
            mid = 0.5 * (lo + hi)
            shift = a_centre + d * mid - p_centre
            nearest = a_tree.query(rotated + shift, k=1)[0].min()
            if nearest < gap:
                lo = mid           # too close: move out
            else:
                hi = mid           # clear: move in
        shift = a_centre + d * hi - p_centre
        out.append(Pose(rotation=R, translation=shift, centre=p_centre,
                        label=f"decoy-{i:04d}"))
    return out


def perturbed_poses(
    partner_atoms: list[Atom],
    offsets: list[float],
    seed: int = 0,
    rotation_per_angstrom: float = 2.0,
) -> list[Pose]:
    """Poses at a known displacement from where the partner already is.

    For checking that a score degrades the way it should as a pose gets
    worse, without needing a docking program to produce the bad poses.
    """
    rng = np.random.default_rng(seed)
    centre = coords(heavy(partner_atoms)).mean(axis=0)
    out = []
    for off in offsets:
        if off == 0.0:
            out.append(Pose(centre=centre, label="0.0A"))
            continue
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        axis = rng.normal(size=3)
        R = rotation_from_axis_angle(
            axis, np.deg2rad(off * rotation_per_angstrom)
        )
        out.append(Pose(rotation=R, translation=d * off, centre=centre,
                        label=f"{off:.1f}A"))
    return out
