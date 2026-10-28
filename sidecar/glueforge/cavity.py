"""Composite cavity detection and pharmacophore field mapping.

The cavity at a glued interface is the one a molecule has to fill: part of it
is the anchor's pocket, part of it is completed by the partner. This module
finds it on a grid and characterises it with the fields a generative model
would be conditioned on -- shape, enclosure, electrostatics, hydrogen-bond
capability and hydrophobicity.

Method is LIGSITE-style: rasterise the protein, mark free points, keep the
ones enclosed along several axes, take the connected component at the site.
Deliberately not fpocket -- not because fpocket is wrong, but because the
interesting quantity here is how the cavity *moves* when the partner moves,
and that needs the detector to be inside the loop rather than a subprocess.

The electrostatic map uses formal charges on ionisable groups with a
distance-dependent dielectric. That is a crude field: it captures the
dominant charge pattern that decides whether a cavity wants an acid or a
base, and nothing finer. Anything quantitative needs real partial charges
and a Poisson-Boltzmann solver.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from .cif import Atom
from .residues import element_radius

PROBE = 1.4

#: How far a charge can be and still count toward a cavity point's
#: electrostatic character, in angstroms. See the note where it is used:
#: without it the sum reports the whole protein's net charge rather than the
#: site's local environment.
ELECTROSTATIC_CUTOFF = 12.0

#: Formal charges placed on the atoms that carry them at pH 7.
CHARGED_ATOMS = {
    ("ASP", "OD1"): -0.5, ("ASP", "OD2"): -0.5,
    ("GLU", "OE1"): -0.5, ("GLU", "OE2"): -0.5,
    ("LYS", "NZ"): 1.0,
    ("ARG", "NH1"): 0.5, ("ARG", "NH2"): 0.5,
    ("HIS", "ND1"): 0.05, ("HIS", "NE2"): 0.05,
}
#: Termini carry charge too, but chain ends in a crystal are often disordered
#: and their contribution is local, so they are left out rather than guessed.

DONORS = {
    ("SER", "OG"), ("THR", "OG1"), ("TYR", "OH"), ("CYS", "SG"),
    ("ASN", "ND2"), ("GLN", "NE2"), ("TRP", "NE1"), ("HIS", "ND1"),
    ("HIS", "NE2"), ("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"),
    ("ARG", "NH2"),
}
ACCEPTORS = {
    ("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2"),
    ("ASN", "OD1"), ("GLN", "OE1"), ("SER", "OG"), ("THR", "OG1"),
    ("TYR", "OH"), ("HIS", "ND1"), ("HIS", "NE2"),
}
APOLAR = {"C"}


@dataclass
class Cavity:
    """A connected pocket on the grid, with its fields."""

    points: np.ndarray            # (n, 3) cartesian coordinates
    spacing: float
    enclosure: np.ndarray         # 0-3 axes enclosed, per point
    electrostatic: np.ndarray     # kcal/mol/e, sign is the field a + probe feels
    donor: np.ndarray             # proximity score to protein donors
    acceptor: np.ndarray
    hydrophobic: np.ndarray
    from_anchor: np.ndarray       # fraction of nearby protein atoms in the anchor
    centroid: np.ndarray = field(init=False)

    def __post_init__(self):
        self.centroid = (
            self.points.mean(axis=0) if len(self.points) else np.zeros(3)
        )

    @property
    def volume(self) -> float:
        return float(len(self.points) * self.spacing ** 3)

    def enclosure_quality(self) -> float:
        """Share of the cavity enclosed on all three axes.

        Near 1 is a buried pocket of the sort generative models are mostly
        trained on. Near 0 is a surface groove, which is what several glue
        sites actually are -- and a much harder target to fill.
        """
        if not len(self.points):
            return 0.0
        return float((self.enclosure >= 3).mean())

    def composite_fraction(self) -> float:
        """How much of the cavity is lined by both proteins rather than one.

        A value near 0.5 means the cavity genuinely straddles the interface,
        which is what makes it a glue site rather than an ordinary pocket on
        one protein.
        """
        if not len(self.points):
            return 0.0
        mixed = (self.from_anchor > 0.15) & (self.from_anchor < 0.85)
        return float(mixed.mean())

    def summary(self) -> dict:
        if not len(self.points):
            return {"volume": 0.0, "points": 0}
        return {
            "points": int(len(self.points)),
            "volume": round(self.volume, 1),
            "spacing": self.spacing,
            "mean_enclosure": round(float(self.enclosure.mean()), 2),
            "enclosure_quality": round(self.enclosure_quality(), 3),
            "composite_fraction": round(self.composite_fraction(), 3),
            "mean_electrostatic": round(float(self.electrostatic.mean()), 3),
            "electrostatic_range": [
                round(float(self.electrostatic.min()), 2),
                round(float(self.electrostatic.max()), 2),
            ],
            "mean_hydrophobic": round(float(self.hydrophobic.mean()), 2),
            "donor_rich_fraction": round(float((self.donor > 0.5).mean()), 3),
            "acceptor_rich_fraction": round(float((self.acceptor > 0.5).mean()), 3),
        }


def _grid(centre: np.ndarray, half_width: float, spacing: float):
    n = int(2 * half_width / spacing) + 1
    ax = [np.linspace(c - half_width, c + half_width, n) for c in centre]
    return ax, n


def detect_cavity(
    protein: list[Atom],
    site_centre: np.ndarray,
    anchor_atoms: set[int] | None = None,
    half_width: float = 14.0,
    spacing: float = 0.8,
    min_enclosure: int = 2,
    probe: float = PROBE,
) -> Cavity:
    """Find the enclosed free volume around `site_centre`.

    `anchor_atoms` holds indices into `protein` that belong to the anchor
    chain, which lets the cavity report how much of its lining each protein
    contributes. Indices refer to `protein` as passed in; the hydrogen
    filtering below is invisible to the caller.

    Hydrogens are discarded for the same reason `atom_sasa` discards them:
    whether a structure has them is a choice of the depositor, and a modelled
    hydrogen projects about 1.1 angstroms past its parent heavy atom, so
    keeping them shrinks every cavity in a structure that happens to have
    them and leaves cavity volume incomparable between depositions. Since
    cavity volume is one of the numbers the gate is calibrated on, that
    incomparability would propagate straight into the verdict.
    """
    if any(a.element.upper() == "H" for a in protein):
        keep = [i for i, a in enumerate(protein) if a.element.upper() != "H"]
        remap = {old: new for new, old in enumerate(keep)}
        if anchor_atoms is not None:
            anchor_atoms = {remap[i] for i in anchor_atoms if i in remap}
        protein = [protein[i] for i in keep]

    if not protein:
        raise ValueError("no heavy atoms to build a cavity from")

    coords = np.array([(a.x, a.y, a.z) for a in protein], dtype=np.float64)
    radii = np.array([element_radius(a.element) for a in protein])

    ax, n = _grid(site_centre, half_width, spacing)
    gx, gy, gz = np.meshgrid(*ax, indexing="ij")
    pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)

    tree = cKDTree(coords)
    # A grid point is occupied if it falls inside any atom's vdW sphere.
    nearest_d, nearest_i = tree.query(pts, k=1)
    occupied = nearest_d < radii[nearest_i]
    occ = occupied.reshape(n, n, n)

    # LIGSITE-style enclosure: along each axis, is there protein on both
    # sides? Points enclosed on several axes sit in a pocket rather than on
    # an open surface.
    enclosure = np.zeros((n, n, n), dtype=np.int8)
    for axis in (0, 1, 2):
        fwd = np.maximum.accumulate(occ, axis=axis)
        rev = np.flip(np.maximum.accumulate(np.flip(occ, axis=axis), axis=axis), axis=axis)
        enclosure += (fwd & rev).astype(np.int8)

    free = (~occ) & (enclosure >= min_enclosure)
    # Points must also be far enough from protein to hold a probe.
    clearance = (nearest_d - radii[nearest_i]).reshape(n, n, n)
    free &= clearance > probe * 0.65
    # Bridge the thin necks that fragment a shallow groove into dozens of
    # pieces. Without this, DCAF15's site comes back as 13 cubic angstroms.
    free = ndimage.binary_closing(free, np.ones((3, 3, 3))) & (~occ)

    if not free.any():
        # Six arrays, not five: enclosure, electrostatic, donor, acceptor,
        # hydrophobic, from_anchor.
        empty = np.zeros((0, 3))
        return Cavity(empty, spacing, *[np.zeros(0) for _ in range(6)])

    labels, n_lab = ndimage.label(free)
    # Keep the component nearest the site centre.
    centre_idx = tuple(
        int(np.clip(round((c - (s - half_width)) / spacing), 0, n - 1))
        for c, s in zip(site_centre, site_centre)
    )
    want = labels[centre_idx]
    if want == 0:
        sizes = ndimage.sum(free, labels, range(1, n_lab + 1))
        centres = ndimage.center_of_mass(free, labels, range(1, n_lab + 1))
        best, best_d = 1, 1e9
        for li, com in enumerate(centres, start=1):
            if sizes[li - 1] < 20:
                continue
            world = np.array([ax[k][int(round(com[k]))] for k in range(3)])
            d = np.linalg.norm(world - site_centre)
            if d < best_d:
                best_d, best = d, li
        want = best

    sel = (labels == want).reshape(-1)
    cav_pts = pts[sel]
    cav_enc = enclosure.reshape(-1)[sel]

    fields = _fields(cav_pts, protein, coords, tree, anchor_atoms)
    return Cavity(cav_pts, spacing, cav_enc, *fields)


def _fields(
    pts: np.ndarray,
    protein: list[Atom],
    coords: np.ndarray,
    tree: cKDTree,
    anchor_atoms: set[int] | None,
):
    """Electrostatic, H-bond and hydrophobic maps on the cavity points."""
    n = len(pts)
    elec = np.zeros(n)
    donor = np.zeros(n)
    acceptor = np.zeros(n)
    phobic = np.zeros(n)
    frac_anchor = np.zeros(n)

    charges, charge_xyz = [], []
    don_xyz, acc_xyz, apolar_xyz = [], [], []
    for i, a in enumerate(protein):
        key = (a.comp_id.upper(), a.name)
        q = CHARGED_ATOMS.get(key)
        if q:
            charges.append(q)
            charge_xyz.append((a.x, a.y, a.z))
        if key in DONORS:
            don_xyz.append((a.x, a.y, a.z))
        if key in ACCEPTORS:
            acc_xyz.append((a.x, a.y, a.z))
        if a.element.upper() in APOLAR:
            apolar_xyz.append((a.x, a.y, a.z))

    # Coulomb with a distance-dependent dielectric, eps = 4r, cut off at
    # ELECTROSTATIC_CUTOFF. Crude, but it is the term that decides whether a
    # cavity wants an acid or a base.
    #
    # The cutoff is not a performance shortcut, it is the physics. Summed
    # over every charge in both proteins this reports the proteins' overall
    # acidity as though it were a property of the site: at 6TD3, 86% of the
    # field came from charges more than 12 angstroms away, and at 6H0F the
    # unbounded sum had the opposite sign to the local environment. Real
    # solvent and counterions screen long-range fields almost completely,
    # and a distance-dependent dielectric only partly stands in for that.
    if charges:
        q = np.array(charges)
        cxyz = np.array(charge_xyz)
        for start in range(0, n, 4000):
            chunk = pts[start:start + 4000]
            d = np.linalg.norm(chunk[:, None, :] - cxyz[None, :, :], axis=2)
            d = np.maximum(d, 1.5)
            contrib = 332.0 * q[None, :] / (4.0 * d * d)
            contrib[d > ELECTROSTATIC_CUTOFF] = 0.0
            elec[start:start + 4000] = contrib.sum(axis=1)

    def proximity(xyz, cutoff=4.5):
        if not xyz:
            return np.zeros(n)
        t = cKDTree(np.array(xyz))
        out = np.zeros(n)
        for i, idx in enumerate(t.query_ball_point(pts, cutoff)):
            if idx:
                d = np.linalg.norm(np.array(xyz)[idx] - pts[i], axis=1)
                out[i] = float(np.sum(1.0 - d / cutoff))
        return out

    donor = proximity(don_xyz)
    acceptor = proximity(acc_xyz)
    phobic = proximity(apolar_xyz, 5.0)

    if anchor_atoms:
        for i, idx in enumerate(tree.query_ball_point(pts, 6.0)):
            if idx:
                frac_anchor[i] = np.mean([j in anchor_atoms for j in idx])

    return elec, donor, acceptor, phobic, frac_anchor


def cavity_overlap(a: Cavity, b: Cavity, tol: float = 1.0) -> dict:
    """How much two cavities describe the same volume.

    Jaccard on occupied grid cells, plus how much of the smaller one the
    larger covers, which is the more forgiving measure and the relevant one
    when a perturbation shrinks the pocket rather than moving it.
    """
    if not len(a.points) or not len(b.points):
        return {"jaccard": 0.0, "coverage": 0.0,
                "centroid_shift": float("nan"),
                "volume_ratio": 0.0}

    ta, tb = cKDTree(a.points), cKDTree(b.points)
    in_b = np.array([len(x) > 0 for x in ta.query_ball_tree(tb, tol)])
    in_a = np.array([len(x) > 0 for x in tb.query_ball_tree(ta, tol)])
    inter = 0.5 * (in_b.sum() + in_a.sum())
    union = len(a.points) + len(b.points) - inter

    return {
        "jaccard": round(float(inter / union), 4) if union else 0.0,
        "coverage": round(float(in_b.mean()), 4),
        "centroid_shift": round(float(np.linalg.norm(a.centroid - b.centroid)), 2),
        "volume_ratio": round(float(b.volume / max(a.volume, 1e-6)), 3),
    }


def field_correlation(a: Cavity, b: Cavity, tol: float = 1.0) -> dict:
    """Do two cavities agree about their chemistry where they overlap?

    Shape stability is necessary but not sufficient: a generative model is
    conditioned on the fields, so the fields have to be stable too.
    """
    if not len(a.points) or not len(b.points):
        return {}
    tb = cKDTree(b.points)
    d, idx = tb.query(a.points, k=1)
    ok = d <= tol
    if ok.sum() < 20:
        return {"n_overlap": int(ok.sum())}

    out = {"n_overlap": int(ok.sum())}
    for name in ("electrostatic", "donor", "acceptor", "hydrophobic"):
        va = getattr(a, name)[ok]
        vb = getattr(b, name)[idx[ok]]
        if va.std() < 1e-9 or vb.std() < 1e-9:
            out[name] = None
        else:
            out[name] = round(float(np.corrcoef(va, vb)[0, 1]), 4)
    return out


def rigid_perturb(
    atoms: list[Atom], translation: np.ndarray, rotation_deg: float, seed: int
) -> list[Atom]:
    """Copy a set of atoms with a rigid-body displacement applied.

    Used to ask how much docking error a cavity can absorb before it stops
    being the same cavity.
    """
    rng = np.random.default_rng(seed)
    xyz = np.array([(a.x, a.y, a.z) for a in atoms], dtype=np.float64)
    centre = xyz.mean(axis=0)

    axis = rng.normal(size=3)
    axis /= np.linalg.norm(axis)
    th = np.deg2rad(rotation_deg)
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    R = np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)

    moved = (xyz - centre) @ R.T + centre + translation
    out = []
    for a, p in zip(atoms, moved):
        out.append(Atom(a.serial, a.name, a.element, a.comp_id, a.chain,
                        a.seq_id, float(p[0]), float(p[1]), float(p[2]),
                        a.occupancy, a.hetatm, a.insertion))
    return out
