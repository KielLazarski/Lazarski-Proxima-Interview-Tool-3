"""Rigid-body docking by FFT correlation, with the presented face as a filter.

This exists because LightDock cannot be installed on Windows -- no wheel for
any version, and upstream does not claim to support it -- and asking people to
set up a Linux environment to run one stage of a tool is a poor trade. What
is here needs numpy and scipy and nothing else.

The method is Katchalski-Katzir correlation docking, which is the 1992 idea
underneath FTDock and ZDOCK. Both proteins are drawn onto the same lattice.
The receptor is drawn twice over: a thin shell at its surface worth +1 a
cell, and its interior worth -PENALTY a cell. The partner is drawn once, as
solid body. The sum of the product over every relative position is then a
score that rewards surface touching surface and punishes one protein inside
the other -- and because that sum over all translations is a cross-correlation,
one FFT pair computes every translation at once. The search is therefore only
over orientations, which is three dimensions instead of six.

Two things here are not in the textbook version, and both come from what this
tool is for.

**The ligand is part of the receptor.** LightDock is given protein only, so
the ligand's pocket is empty as far as the docking run is concerned and poses
that settle into it score best. `dock.py` has to detect and discard those
afterwards. Here the ligand is simply drawn into the receptor grid along with
the anchor, so a pose that tries to occupy the pocket is penalised while it is
being generated rather than filtered out after.

**The presented face is a hard filter, not a scoring bias.** A placement is
kept only if some partner atom comes within contact of an exposed ligand
atom. Most of the translation map cannot satisfy that however well it scores,
so it is excluded before anything is ranked -- which is both the restraint
`dock.py` hands LightDock as a residue list and, incidentally, most of what
makes the search affordable. The lattice narrows it; the surviving peaks are
then checked against real coordinates rather than cells, because a filter
that is the whole point of the method should not be quantised.

The search runs in two sweeps. The first is a uniform set of orientations
and only has to get the right neighbourhood into its upper ranks; the second
spends the rest of the budget inside the neighbourhoods the first liked. That
is worth the complication because of how sharply the fit dies with
orientation error -- see the note above PRESETS, which is measured rather
than borrowed from the docking literature.

What this does not do: side chains do not move, there is no desolvation or
electrostatic term anywhere in it, and nothing models conformational change
on binding. It is a coarse search whose job is to put a near-native pose
somewhere in its output, after which `PoseScorer` -- the part validated
against 210 adversarial decoys -- does the ranking. Whether it succeeds at
that job is measured by `run_builtin_recovery.py`, and the answer is a
number, not a claim.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

import numpy as np

from .cif import Atom
from .face import PresentedFace
from .pose import Pose, coords, heavy, rotation_from_axis_angle

# --- the grid ---------------------------------------------------------------

#: Cells within this of a receptor heavy atom are its interior. Roughly a
#: carbon's van der Waals radius less a little, so that two surfaces in van
#: der Waals contact do not read as interpenetrating.
CORE_RADIUS = 1.9

#: Interior cells out to here are the surface shell that scores. 1.4 A is the
#: radius of a water, so the shell is about where a partner atom would sit.
SKIN = 1.4

#: What a cell of partner inside the receptor costs, against +1 for a cell in
#: the surface shell. The classic value. Large enough that no amount of
#: surface contact pays for burying one protein in the other.
PENALTY = 15.0

#: A partner cell must come this close to an exposed ligand atom to count as
#: burying the face. Heavy-atom contact distance, as everywhere else here.
FACE_CONTACT = 5.0


@dataclass(frozen=True)
class Settings:
    """How hard to search. Time scales with `rotations`, and everything else
    is fixed by the geometry."""

    key: str
    label: str
    rotations: int
    spacing: float = 1.2
    #: Translations kept per orientation before the global cut. Measured,
    #: not guessed: at the best orientation a 2000-point set offers on 6TD3,
    #: the peak closest to the right answer was the nineteenth. Keeping ten
    #: threw it away, and the transform that found it had already been paid
    #: for -- pulling more peaks out of a map already in memory is nearly
    #: free, so this is generous on purpose.
    per_rotation: int = 30
    #: Poses handed to the scorer. This cut ranks by the lattice score,
    #: which is not the same thing as being close to the right answer, so it
    #: can discard the best pose in the run -- and because searching harder
    #: adds candidates that compete for the same slots, a larger search can
    #: return a worse best-of than a smaller one. Generous for that reason.
    keep: int = 400
    #: Two kept translations for one orientation must differ by this much,
    #: or they are the same pose reported twice.
    separation: float = 6.0
    #: Walk the shortlist downhill afterwards. The lattice cannot resolve
    #: better than the gap between its orientations, so without this stage
    #: the output is neighbourhoods rather than placements.
    refine: bool = True
    #: Orientations to revisit closely after the first sweep, and how many
    #: neighbours to try around each. Zero turns the second sweep off.
    focus: int = 80
    focus_each: int = 30

    @property
    def minutes_estimate(self) -> tuple[int, int]:
        # Timed on two cores at 1.4 A against a 6500-atom partner, which is
        # the largest in the benchmark and therefore the slowest: 0.155 s an
        # orientation, scaling with the cells in the box. Most machines have
        # more cores than that and most partners are smaller, so the low end
        # is what to expect and the high end is the bad case.
        n = self.rotations + (self.focus * self.focus_each if self.focus else 0)
        per = n * (1.4 / self.spacing) ** 3 / 1000.0
        return (max(1, int(round(1.0 * per))), max(2, int(round(3.0 * per))))

    def as_dict(self) -> dict:
        lo, hi = self.minutes_estimate
        return {
            "key": self.key,
            "label": self.label,
            "engine": "builtin",
            "rotations": self.rotations,
            "focus": [self.focus, self.focus_each],
            "spacing": self.spacing,
            "keep": self.keep,
            "minutes_estimate": [lo, hi],
        }


# How many orientations is not a free choice, and the first guess at it was
# far too low. Measured on 6TD3 by rotating the crystal partner by a known
# angle and taking the best score any translation could then reach:
#
#     off native    0d     3d     5d     8d    12d    20d
#     best score  +190   +126   +136   +104   -291   nothing
#
# The fit does not degrade gently, it dies: by twelve degrees no placement
# anywhere near the right one scores above zero, and the search goes off and
# returns a confident miss somewhere else. Eight degrees is the working
# limit, and a 400-orientation set -- which is what this shipped with first
# -- leaves gaps of fifteen. That is the whole of why it was returning 28 A.
#
# Sets of 2000, 5000 and 12000 points leave gaps of roughly 8, 5 and 4
# degrees. There is no cheaper way around it: the orientation is the one
# dimension the transform cannot search for free.
PRESETS: dict[str, Settings] = {
    "quick": Settings("quick", "Quick", rotations=3000, spacing=1.6,
                      per_rotation=24, keep=200),
    "standard": Settings("standard", "Standard", rotations=8000, spacing=1.4,
                         keep=400),
    "thorough": Settings("thorough", "Thorough", rotations=20000, spacing=1.2,
                         per_rotation=36, keep=600),
}


# --- orientations -----------------------------------------------------------


def super_fibonacci(n: int) -> np.ndarray:
    """`n` unit quaternions spread evenly over SO(3).

    Alexa's Super-Fibonacci spiral. Deterministic, and low-discrepancy rather
    than random: a random set of the same size leaves holes and doubles up
    elsewhere, and a hole is a pose the search cannot reach however long it
    runs. There is no seed here on purpose -- a docking run that found the
    answer should find it again.
    """
    phi = math.sqrt(2.0)
    psi = 1.533751168755204288118041  # the real root of x^4 = x + 4
    i = np.arange(n, dtype=np.float64) + 0.5
    s = i / n
    r, R = np.sqrt(s), np.sqrt(1.0 - s)
    alpha = 2.0 * math.pi * i / phi
    beta = 2.0 * math.pi * i / psi
    return np.stack([r * np.sin(alpha), r * np.cos(alpha),
                     R * np.sin(beta), R * np.cos(beta)], axis=1)


def covering_degrees(n: int) -> float:
    """Roughly how far any orientation can be from the nearest of `n`.

    SO(3) has finite volume, so the gap shrinks as the cube root of the
    count. Calibrated against the measured distance from the identity to the
    nearest member of sets of 400, 1200 and 4000 points, which came out at
    14.9, 9.1 and 5.2 degrees.
    """
    return 14.9 * (400.0 / max(1, n)) ** (1.0 / 3.0)


def neighbours(q: np.ndarray, n: int, spread_degrees: float) -> np.ndarray:
    """`n` orientations near `q`, within about `spread_degrees` of it.

    Used for the second sweep: the first one says which neighbourhoods are
    worth looking at, and this fills them in. Spending the whole budget on a
    uniform set instead means most of the transforms are spent on
    orientations the first sweep already ruled out, and the one dimension
    the transform cannot search for free is the one that needed the budget.
    """
    # A low-discrepancy set on the sphere of small rotations, composed with
    # the base orientation. Deterministic, like the first sweep.
    i = np.arange(n, dtype=np.float64) + 0.5
    z = 1.0 - 2.0 * i / n
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, 1.0))
    phi = i * math.pi * (3.0 - math.sqrt(5.0))
    axes = np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)
    # Angles spread through the ball rather than all at its surface, so the
    # set covers the neighbourhood instead of ringing it.
    ang = np.radians(spread_degrees) * np.cbrt((i / n))

    half = ang / 2.0
    w = np.cos(half)[:, None]
    v = axes * np.sin(half)[:, None]
    small = np.concatenate([w, v], axis=1)
    return np.array([_quat_mul(s, q) for s in small])


def _quat_mul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return np.array([
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    ])


def quaternion_to_matrix(q: np.ndarray) -> np.ndarray:
    w, x, y, z = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
    ])


# --- lattice ----------------------------------------------------------------


def _stamp(radius: float, spacing: float) -> np.ndarray:
    """Lattice offsets within `radius` of a point, as an (n, 3) int array."""
    k = int(math.ceil(radius / spacing))
    rng = np.arange(-k, k + 1)
    gx, gy, gz = np.meshgrid(rng, rng, rng, indexing="ij")
    d2 = (gx * spacing) ** 2 + (gy * spacing) ** 2 + (gz * spacing) ** 2
    m = d2 <= radius * radius
    return np.stack([gx[m], gy[m], gz[m]], axis=1).astype(np.int64)


def _mark(shape: tuple[int, int, int], idx: np.ndarray,
          stamp: np.ndarray) -> np.ndarray:
    """Boolean grid with every cell within the stamp of every point set.

    One vectorised pass per stamp offset rather than per atom: a stamp of a
    couple of hundred offsets over a few thousand atoms is a few hundred numpy
    operations, where the obvious loop would be a few million.
    """
    out = np.zeros(shape, dtype=bool)
    nx, ny, nz = shape
    for dx, dy, dz in stamp:
        i = idx[:, 0] + dx
        j = idx[:, 1] + dy
        k = idx[:, 2] + dz
        ok = ((i >= 0) & (i < nx) & (j >= 0) & (j < ny) & (k >= 0) & (k < nz))
        out[i[ok], j[ok], k[ok]] = True
    return out


class _Lattice:
    """A shared coordinate frame for both proteins and the face."""

    def __init__(self, origin: np.ndarray, shape: tuple[int, int, int],
                 spacing: float) -> None:
        self.origin = origin
        self.shape = shape
        self.spacing = spacing

    def index(self, xyz: np.ndarray) -> np.ndarray:
        return np.rint((xyz - self.origin) / self.spacing).astype(np.int64)

    def occupancy(self, xyz: np.ndarray, radius: float) -> np.ndarray:
        if not len(xyz):
            return np.zeros(self.shape, dtype=bool)
        return _mark(self.shape, self.index(xyz),
                     _stamp(radius, self.spacing))


# --- the search -------------------------------------------------------------


@dataclass
class Candidate:
    """One placement, before anything physical has been said about it."""

    rotation: np.ndarray
    translation: np.ndarray
    shape_score: float
    #: Closest approach to an exposed ligand atom, in angstroms. Kept because
    #: a pose that only grazes the face is a different thing from one that
    #: sits on it, and the number says which without re-deriving it.
    face_distance: float


def _top_separated(values: np.ndarray, where: np.ndarray, n: int,
                   min_sep_cells: int) -> list[tuple[tuple[int, int, int], float]]:
    """The n highest cells among `where`, no two within `min_sep_cells`.

    `where` is a flat index array, not a boolean grid, and that is the point:
    the allowed region does not change between orientations, so it is found
    once and only those cells are ever looked at. Masking the whole grid
    instead built a multi-million-cell array on every orientation, which on
    the benchmark cost about as much as the transform it was feeding on.

    The separation matters as much as the ranking. Without it every one of
    the n results is the same pose moved by a lattice step, and the output
    looks like many answers when it is one.
    """
    if not len(where):
        return []
    vals = values.ravel()[where]
    take = int(min(vals.size, max(n * 60, 1500)))
    cand = np.argpartition(vals, -take)[-take:]
    cand = cand[np.argsort(vals[cand])[::-1]]
    coords_ = np.stack(np.unravel_index(where[cand], values.shape), axis=1)

    picked: list[tuple[tuple[int, int, int], float]] = []
    pts = np.empty((n, 3), dtype=np.int64)
    k = 0
    for c in range(len(cand)):
        p = coords_[c]
        if k and np.abs(pts[:k] - p).max(axis=1).min() < min_sep_cells:
            continue
        picked.append((tuple(int(x) for x in p), float(vals[cand[c]])))
        pts[k] = p
        k += 1
        if k >= n:
            break
    return picked


# --- refinement -------------------------------------------------------------
#
# The lattice scan cannot do better than the gap between its orientations,
# and that gap is eight degrees at best. Eight degrees is several angstroms
# of error on a protein-sized body, so the scan's own output is a shortlist
# of neighbourhoods rather than a set of answers. Every FFT docking program
# has a stage like this one after it, for the same reason.
#
# The scan's score cannot be used here: it is quantised to the lattice, so it
# is flat over a cell and then steps, and a local optimiser on a staircase
# stops at the first tread. This is the same measurement made continuously,
# on coordinates.

#: Closer than this and two heavy atoms are overlapping.
SOFT_CLASH = 3.2
#: Out to here they are in contact.
SOFT_CONTACT = 5.0
#: What one angstrom of overlap costs, against +1 for a contacting atom.
SOFT_PENALTY = 8.0
#: What burying one exposed ligand atom is worth, over and above the contact
#: it already counts as. Zero, and that is a correction rather than a
#: default.
#:
#: It was 6.0, on the reasoning that the face is the point of the exercise
#: and deserves weighting. The reasoning was wrong twice over. The face is
#: already a hard filter -- every pose here has cleared it, so rewarding it
#: again only ranks poses by how much they pile onto it. And the score this
#: one stands in for, the lattice score that was measured and separates the
#: crystal pose from sixty wrong ones, has no face term at all. Adding one
#: here meant refinement descended a different objective from the one that
#: found the pose, and on 6UD7 it dragged the best placement from 7.85 A to
#: 16.79 A while reporting a better number for it.
FACE_WEIGHT = 0.0


def _soft_score(partner: np.ndarray, rec_tree, face_tree) -> float:
    """Shape complementarity and face burial, continuous in the coordinates.

    Deliberately not an energy. It is the same thing the lattice counts --
    atoms in contact, atoms overlapping, face covered -- written so that it
    has a slope everywhere, which is all a local optimiser needs.
    """
    d, _ = rec_tree.query(partner, k=1, workers=-1)
    overlap = np.clip(SOFT_CLASH - d, 0.0, None)
    contact = np.clip(
        np.minimum(d - SOFT_CLASH, SOFT_CONTACT + 1.5 - d) / 1.5, 0.0, 1.0)
    shape = float(contact.sum() - SOFT_PENALTY * overlap.sum())

    fd, _ = face_tree.query(partner, k=1, workers=-1)
    buried = float(np.clip((SOFT_CONTACT + 1.5 - fd) / 1.5, 0.0, 1.0).sum())
    return shape + FACE_WEIGHT * buried


def refine(
    pose: Pose,
    partner_xyz: np.ndarray,
    rec_tree,
    face_tree,
    max_degrees: float = 12.0,
    max_shift: float = 6.0,
) -> tuple[Pose, float]:
    """Walk one pose downhill in six dimensions. Returns it and its score.

    Bounded on purpose. An unbounded local search from a shortlisted pose
    will happily slide off to a better interface somewhere else on the
    protein, which turns a list of distinct neighbourhoods into many copies
    of whichever one has the deepest well.
    """
    from scipy.optimize import minimize

    centre = pose.centre if pose.centre is not None else partner_xyz.mean(0)
    base = pose.apply(partner_xyz)
    base_centre = base.mean(axis=0)
    local = base - base_centre

    def place(p: np.ndarray) -> np.ndarray:
        rot = rotation_from_axis_angle(
            p[:3], float(np.linalg.norm(p[:3])) or 0.0)
        return local @ rot.T + base_centre + p[3:]

    def cost(p: np.ndarray) -> float:
        return -_soft_score(place(p), rec_tree, face_tree)

    lim = math.radians(max_degrees)
    bounds = [(-lim, lim)] * 3 + [(-max_shift, max_shift)] * 3
    res = minimize(cost, np.zeros(6), method="Powell", bounds=bounds,
                   options={"maxiter": 20, "xtol": 0.05, "ftol": 0.05})

    p = np.asarray(res.x, dtype=np.float64)
    extra = rotation_from_axis_angle(p[:3], float(np.linalg.norm(p[:3])) or 0.0)

    # Fold the adjustment into the pose it started from, so the result is one
    # rigid transform of the original coordinates rather than two stacked.
    # The translation is read back off the coordinates the optimiser actually
    # settled on rather than derived, which keeps this correct whatever
    # `pose.centre` happens to be -- an earlier version quietly assumed it
    # was the partner's centroid and would have been wrong for any pose where
    # it was not.
    final = place(p)
    combined = extra @ pose.rotation
    rotated = (partner_xyz - centre) @ combined.T
    translation = (final - rotated - centre).mean(axis=0)
    return (Pose(rotation=combined, translation=translation, centre=centre,
                 label=pose.label),
            float(-res.fun))


def _distinct(candidates: list["Candidate"], local: np.ndarray, keep: int,
              separation: float) -> list["Candidate"]:
    """The best `keep` candidates, no two of them the same placement.

    The cut that used to be here was `sorted_by_score[:keep]`, and it had a
    hole big enough to lose the answer through. `_top_separated` keeps
    distinct translations *within* one orientation, but nothing kept distinct
    placements *across* orientations -- and the focused second sweep
    deliberately produces thirty near-identical orientations around each
    promising one, every one of which reports its own version of the same
    peak. So the global list filled with one well reported many times, and
    genuinely different placements further down were cut.

    Sameness is judged on where the placement actually puts the protein, by
    moving a handful of spread-out probe atoms and comparing them. Centroid
    alone would call two opposite orientations about the same point one pose.
    """
    if not candidates:
        return []

    # Probes: the extremes along each axis plus the centre, which between
    # them move under any rotation or translation that matters.
    n = len(local)
    idx = [0, n - 1]
    for ax in range(3):
        idx += [int(np.argmin(local[:, ax])), int(np.argmax(local[:, ax]))]
    probes = local[sorted(set(idx))]

    # Only the strongest candidates are worth the comparison; past twenty
    # times the quota nothing survives the greedy pass anyway.
    pool = candidates[: max(keep * 20, 2000)]
    moved = np.einsum("pj,cij->cpi",
                      probes, np.array([c.rotation for c in pool])) \
        + np.array([c.translation for c in pool])[:, None, :]

    taken: list[int] = []
    acc = np.empty((keep, moved.shape[1], 3))
    for i in range(len(pool)):
        if taken:
            d = np.abs(acc[: len(taken)] - moved[i]).max(axis=(1, 2))
            if d.min() < separation:
                continue
        acc[len(taken)] = moved[i]
        taken.append(i)
        if len(taken) >= keep:
            break
    return [pool[i] for i in taken]


def _closest(points: np.ndarray, targets: np.ndarray) -> float:
    """Smallest distance between any point and any target.

    Direct rather than through a tree: the target set is a dozen exposed
    ligand atoms, and building a tree for a dozen points costs more than the
    distances do.
    """
    if not len(points) or not len(targets):
        return float("inf")
    d2 = ((points[:, None, :] - targets[None, :, :]) ** 2).sum(axis=2)
    return float(np.sqrt(d2.min()))


def dock(
    anchor_atoms: list[Atom],
    ligand_atoms: list[Atom],
    partner_atoms: list[Atom],
    face: PresentedFace,
    settings: Settings | str = "standard",
    progress: Callable[[int, int], None] | None = None,
    workers: int = -1,
) -> list[Pose]:
    """Place the partner against the anchor, aimed at the presented face.

    Returns poses in descending order of shape complementarity. That order is
    *not* the answer -- it is a shortlist for `PoseScorer`, which is the part
    of this pipeline that was checked against decoys.
    """
    from scipy import fft as sfft

    if isinstance(settings, str):
        settings = PRESETS[settings]
    sp = settings.spacing

    # An anchor-less run would dock the partner against the bare ligand, which
    # is not a smaller version of this method but a different and meaningless
    # one. Refused rather than returned empty, because an empty result reads
    # as "nothing scored" and this is "nothing was asked".
    if not heavy(anchor_atoms):
        raise ValueError(
            "docking needs the anchor protein; the partner cannot be placed "
            "against the ligand alone"
        )
    rec = coords(heavy(anchor_atoms) + heavy(ligand_atoms))
    par = coords(heavy(partner_atoms))
    if not len(rec) or not len(par):
        return []

    lig = coords(heavy(ligand_atoms))
    face_xyz = lig[[a.index for a in face.exposed]] if face.exposed else lig

    rec_centre, par_centre = rec.mean(axis=0), par.mean(axis=0)
    rec_local, par_local = rec - rec_centre, par - par_centre
    rec_reach = float(np.linalg.norm(rec_local, axis=1).max())
    par_reach = float(np.linalg.norm(par_local, axis=1).max())

    # Which body turns, and how big the box has to be. These are one
    # decision, because the box is sized by whichever body moves.
    #
    # The correlation is cyclic, so anything that wraps round an edge is
    # scored against the opposite side of the other protein, and the box has
    # to be big enough that no pose worth keeping ever does. The restraint
    # bounds that tightly, but differently depending on which body moves,
    # because the presented face belongs to the receptor:
    #
    #   partner moves   the face stands still, so the partner's centre is
    #                   within (its reach + contact) of the face and the
    #                   partner reaches one more `reach` past that.
    #   receptor moves  the face travels with it and may meet the partner
    #                   anywhere on the partner's surface, so the bound
    #                   carries the partner's reach as well.
    #
    # The second is usually the larger of the two and sometimes much smaller,
    # when the anchor is a small protein and the partner a large one. Both
    # are computed and the cheaper is taken. Nothing about accuracy turns on
    # this: a relative orientation error of eight degrees is eight degrees
    # whichever body was turned to produce it.
    half_fixed_face = 2.0 * par_reach + FACE_CONTACT
    half_moving_face = par_reach + FACE_CONTACT + 2.0 * rec_reach
    swapped = half_moving_face < half_fixed_face

    if swapped:
        mob_local, mob_centre, mob_reach = rec_local, rec_centre, rec_reach
        stat, stat_centre, stat_reach = par, par_centre, par_reach
        face_local = face_xyz - rec_centre      # the face travels with it
        centre_on, half = par_centre, half_moving_face
    else:
        mob_local, mob_centre, mob_reach = par_local, par_centre, par_reach
        stat, stat_centre, stat_reach = rec, rec_centre, rec_reach
        face_local = None
        centre_on, half = face_xyz.mean(axis=0), half_fixed_face

    half += CORE_RADIUS + SKIN + sp
    lo, hi = centre_on - half, centre_on + half
    dims = tuple(int(sfft.next_fast_len(int(math.ceil((h - l) / sp)) + 1))
                 for l, h in zip(lo, hi))
    lat = _Lattice(lo, dims, sp)

    # Atoms of the fixed body outside the box cannot be reached by any pose
    # the restraint allows, so they are dropped. Not an approximation: a
    # placement that cannot bring the partner within FACE_CONTACT of the face
    # is absent from the output either way.
    stat = stat[np.all((stat >= lo) & (stat <= hi), axis=1)]
    if not len(stat):
        return []

    core = lat.occupancy(stat, CORE_RADIUS)
    shell = lat.occupancy(stat, CORE_RADIUS + SKIN) & ~core
    S = np.where(core, -PENALTY, 0.0) + shell.astype(np.float64)
    FS = sfft.rfftn(S, workers=workers)
    sep_cells = max(1, int(round(settings.separation / sp)))

    axes = []
    for d in dims:
        a = np.arange(d)
        axes.append(np.where(a > d // 2, a - d, a).astype(np.float64) * sp)
    shifts = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1)

    def allowed_for(target: np.ndarray, radius: float) -> np.ndarray:
        """Flat indices of translations that put `target` within `radius`.

        A necessary condition only -- it asks whether the two could touch at
        the face, not whether they do. Every survivor is checked exactly
        afterwards, against coordinates rather than cells.
        """
        d2 = ((axes[0] - target[0]) ** 2)[:, None, None] \
            + ((axes[1] - target[1]) ** 2)[None, :, None] \
            + ((axes[2] - target[2]) ** 2)[None, None, :]
        return np.flatnonzero((d2 <= radius * radius).ravel())

    # When the partner moves, the face stands still and the allowed region is
    # the same for every orientation, so it is found once. When the receptor
    # moves the face goes with it and the region has to be rebuilt each time.
    fixed_allowed = None
    if not swapped:
        fixed_allowed = allowed_for(face_xyz.mean(axis=0) - par_centre,
                                    par_reach + FACE_CONTACT)
        if not len(fixed_allowed):
            return []

    out: list[Candidate] = []

    def sweep(quats: np.ndarray) -> list[float]:
        """One pass over a set of orientations. Returns each one's best
        score, which is what the second pass chooses from."""
        best: list[float] = []
        for n, q in enumerate(quats):
            if progress is not None and n % 25 == 0:
                progress(done + n, total)

            M = quaternion_to_matrix(q)
            moved = mob_local @ M.T + mob_centre
            L = lat.occupancy(moved, CORE_RADIUS).astype(np.float64)
            corr = sfft.irfftn(FS * np.conj(sfft.rfftn(L, workers=workers)),
                               s=dims, workers=workers)

            if swapped:
                face_moved = face_local @ M.T + mob_centre
                allowed = allowed_for(stat_centre - face_moved.mean(axis=0),
                                      stat_reach + FACE_CONTACT)
            else:
                face_moved, allowed = None, fixed_allowed
            if not len(allowed):
                best.append(-np.inf)
                continue

            top = -np.inf
            for (i, j, k), v in _top_separated(corr, allowed,
                                               settings.per_rotation,
                                               sep_cells):
                if v <= 0.0:
                    continue
                s = shifts[i, j, k]
                # The exact test the lattice only approximated: with this
                # placement, does the partner reach the presented face.
                near = (_closest(face_moved + s, stat) if swapped
                        else _closest(moved + s, face_xyz))
                if near > FACE_CONTACT:
                    continue
                top = max(top, v)
                out.append(Candidate(M, s.copy(), float(v), near))
            best.append(top)
        return best

    # Two sweeps. The first is uniform and only has to get the right
    # neighbourhood into its upper ranks; the second spends the rest of the
    # budget inside the neighbourhoods the first liked.
    #
    # This is worth the complication because of how sharply the fit dies
    # with orientation error. Reaching four degrees uniformly takes twenty
    # thousand transforms; reaching it in the places that matter takes three
    # thousand and then eighty neighbourhoods of thirty. The gamble is that
    # the first sweep ranks the right neighbourhood inside its top eighty --
    # on 6TD3 the near-native orientation's best peak was +76 against a
    # global best of about +95, which clears that bar comfortably.
    coarse = super_fibonacci(settings.rotations)
    n_focus = min(settings.focus, len(coarse)) if settings.focus else 0
    total = settings.rotations + n_focus * settings.focus_each
    done = 0

    scores = sweep(coarse)
    done = settings.rotations

    if n_focus and settings.focus_each:
        # Half again the typical gap, because `covering_degrees` is
        # calibrated on the distance to the *nearest* member of a set, which
        # is the typical case rather than the worst one. A neighbourhood
        # narrower than the gap it is meant to fill leaves the hole it was
        # added to close.
        gap = 1.5 * covering_degrees(settings.rotations)
        order = np.argsort(scores)[::-1][:n_focus]
        for idx in order:
            if not np.isfinite(scores[idx]):
                continue
            sweep(neighbours(coarse[idx], settings.focus_each, gap))
            done += settings.focus_each

    if progress is not None:
        progress(total, total)

    out.sort(key=lambda c: -c.shape_score)
    kept = _distinct(out, par_local if not swapped else mob_local,
                     settings.keep, settings.separation)

    if not swapped:
        poses = [
            Pose(rotation=c.rotation, translation=c.translation,
                 centre=par_centre, label="fft-%04d" % i)
            for i, c in enumerate(kept)
        ]
    else:
        # The search moved the receptor; the caller wants the partner moved.
        # Undoing the receptor's transform and applying its inverse to the
        # partner leaves the relative geometry identical. Pinned by
        # `test_a_swapped_result_is_the_same_relative_placement`, because
        # getting this wrong produces poses that are plausible and wrong.
        poses = []
        for i, c in enumerate(kept):
            M, t = c.rotation, c.translation
            poses.append(Pose(
                rotation=M.T,
                translation=(par_centre - rec_centre - t) @ M
                + rec_centre - par_centre,
                centre=par_centre,
                label="fft-%04d" % i,
            ))

    if not settings.refine or not poses:
        return poses

    from scipy.spatial import cKDTree
    rec_tree = cKDTree(rec)
    face_tree = cKDTree(face_xyz)

    # Both the pose and its refinement, never the refinement alone.
    #
    # Measured on 6UD7, where the grid ranks the crystal pose at four times
    # the best wrong placement:
    #
    #     keep 200, no refinement     best pose 7.85 A
    #     keep 200, refined           best pose 16.79 A
    #
    # Refinement was not polishing the answer, it was walking away from it:
    # the continuous score it descends is not the lattice score that found
    # the pose, and where the two disagree it drags a good placement to a
    # worse one and the original is gone. Raising the cut to 20000 moved
    # 7.85 to 7.32, so the cut was not the loss -- this was.
    #
    # Keeping both makes the stage incapable of losing ground. It is the
    # same rule the focused sweep already had to learn, and the third time
    # tonight that a stage added to improve the result made it worse by
    # replacing rather than adding.
    scored = []
    for i, p in enumerate(poses):
        if progress is not None and i % 10 == 0:
            progress(settings.rotations, settings.rotations)
        scored.append((p, _soft_score(p.apply(par), rec_tree, face_tree)))
        scored.append(refine(p, par, rec_tree, face_tree))
    scored.sort(key=lambda t: -t[1])

    # No second deduplication here, and the reason is the whole lesson of
    # this stage.
    #
    # There used to be one, on the sensible-sounding grounds that refinement
    # can pull two distinct starts into the same well and reporting one
    # answer nine times makes a shortlist look like a consensus. But it
    # merged each refined pose with the parent it came from, and then kept
    # whichever of the two scored higher -- which is always the refined one,
    # because that is what the optimiser just maximised. So the parent was
    # discarded by construction, and whenever refinement had walked away
    # from the better placement, the better placement went with it.
    #
    # That is not fixable by weighting the score differently; it was tried.
    # No score is a perfect proxy for being near the right answer, so
    # "replace the pose with a higher-scoring neighbour" can always lose,
    # and the only stage that cannot lose is one that does not replace.
    # `_distinct` already ran before refinement, so what arrives here is
    # distinct; at worst this returns each answer twice, which costs a line
    # in a table and no correctness at all.
    return [
        Pose(rotation=p.rotation, translation=p.translation,
             centre=p.centre, label="fft-%04d" % i)
        for i, (p, _) in enumerate(scored)
    ]
