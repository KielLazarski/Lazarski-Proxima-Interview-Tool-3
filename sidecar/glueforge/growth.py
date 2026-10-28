"""Where to grow from, which way, and into what.

Step 6. The presented face says which ligand atoms a second protein can
reach; this says what is in front of each of them.

The ranking is exposure-led, and that is not an assumption -- `run_growth.py`
tests it against the CDK12 campaign, eight analogues of one purine core
solved in the same site. Exposure predicts which atoms that campaign varied
at a median rank correlation of -0.50, and -0.43 with distance from the
molecule's centroid partialled out, so it is not merely restating that outer
atoms are easier to reach with a synthesis.

What this module adds is the second half of the question. Knowing an atom is
exposed tells you chemistry *can* go there. What decides whether it *should*
is what the cone in front of it contains: how much room, how enclosed, what
charge, whether it wants a donor or an acceptor, whether it is greasy. A
growth point facing 20 cubic angstroms of nothing is not the same
opportunity as one facing an enclosed, negatively charged pocket, and a
ranking by exposure alone cannot tell them apart.

No chemistry toolkit is needed here. RDKit is needed to turn a chosen vector
into an actual molecule; it is not needed to decide which vector is worth
the attempt.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .cavity import CHARGED_ATOMS, ELECTROSTATIC_CUTOFF, Cavity
from .cif import Atom, Ligand, Residue
from .face import PresentedFace
from .geometry import atom_sasa
from .residues import NON_LIGAND

#: Half-angle of the cone in front of a growth point, in degrees.
#:
#: A substituent projects along the exit vector and swings, but not nearly as
#: far as the first version of this allowed. At 50 degrees and 7.5 angstroms
#: the cone spans about nine angstroms across its mouth, which is a chunk of
#: interface rather than a place a substituent goes. The field averaged over
#: that is not the field a group attached here would feel, however correctly
#: each individual value is computed.
CONE_DEGREES = 38.0

#: How far along the vector to look. Beyond this a substituent would be
#: several bonds out and the geometry here stops predicting anything.
CONE_REACH = 6.5

#: Off-axis falloff. Points in the cone are weighted by a Gaussian in their
#: perpendicular distance from the exit vector, so the reported field is
#: dominated by where a substituent actually sits rather than by the widest
#: part of the cone, which holds the most grid points purely by geometry.
CONE_SIGMA = 2.2

# All three of these numbers were chosen by reasoning about what a
# substituent can reach, not measured against anything, so how much the
# answer depends on them is worth knowing. Across six systems, varying the
# angle between 30 and 45 degrees and the reach between 5 and 8 angstroms:
#
#   Kendall tau against the shipped order   0.81 to 0.96
#   the top-ranked atom is unchanged        in 4 to 5 of 6
#   top-three overlap                       0.67 to 1.00
#   sigma has no effect on the order at all (tau 1.000)
#
# So the best one or two points are robust to these choices and the exact
# order below them is not. Taken with the fact that the order is unvalidated
# against the one campaign available, read the top of the list and ignore the
# precision of the rest.

#: An atom keeping less than this share of its unhindered surface has no
#: room to grow from, whatever else is true of it.
MIN_RELATIVE_SASA = 0.08


def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    return v / n if n > 1e-9 else np.array([0.0, 0.0, 1.0])


@dataclass
class GrowthPoint:
    """One place on the ligand, and what lies in front of it."""

    index: int
    name: str
    element: str
    position: np.ndarray
    direction: np.ndarray

    sasa: float
    relative_sasa: float

    #: Cavity points inside the cone, turned into a volume.
    reachable_volume: float
    #: Mean enclosure of those points. Open solvent scores near zero.
    enclosure: float
    #: Field means over the cone. Sign conventions follow cavity.py.
    electrostatic: float
    donor: float
    acceptor: float
    hydrophobic: float
    #: Share of the cone's lining contributed by the partner rather than the
    #: anchor. A growth point facing only the anchor is growing into the
    #: protein the ligand is already bound to, which buys no new interface.
    toward_partner: float

    #: Partner residues within reach of the cone, nearest first.
    partner_contacts: list[str] = field(default_factory=list)
    #: The charged residues actually driving the field here, with the sign of
    #: each one's contribution.
    #:
    #: Listing only the nearest residues alongside a field value invites a
    #: reading the numbers do not support: 6TD3's best growth point sits
    #: beside an arginine and wants a basic group, which looks like a
    #: contradiction until you see that two aspartates and a glutamate
    #: further out dominate the sum. Showing what drives the number next to
    #: the number is the difference between a result and a puzzle.
    field_drivers: list[dict] = field(default_factory=list)
    n_cavity_points: int = 0

    @property
    def score(self) -> float:
        """How much room there is here, and whose it is.

        **This is a description of available space, not a prediction of where
        to grow, and the distinction is load-bearing.** The exposure gate
        above is validated against the CDK12 campaign: how exposed an atom is
        predicts whether that campaign changed it, at -0.50 and -0.43 with
        peripherality partialled out. This ordering is not. Checked the same
        way it comes out at +0.37 -- it puts highest the atoms the chemists
        left alone.

        The reason looks structural. The atoms that campaign did change sit
        closer to the partner, in space a rigid cavity does not show as free:
        real elaboration remodels an interface and this measures the room
        already in it. Nothing here will fix that, so the honest thing is to
        label what this is rather than tune the weights until the number
        comes out the other way on a sample of one campaign.
        """
        room = min(1.0, self.reachable_volume / 60.0)
        return float(
            0.45 * room
            + 0.35 * self.toward_partner
            + 0.20 * min(1.0, self.enclosure / 3.0)
        )

    def suggests(self) -> str:
        """What kind of group the cone in front of this atom would suit.

        A description of the environment, not a prescription, and a coarse
        one: the electrostatics are a Coulomb sum over formal charges with a
        distance-dependent dielectric and no explicit solvent, which is
        enough to say whether a region leans positive or negative and not
        nearly enough to rank two real compounds against each other.

        Whether such a group can be attached here, and whether it would be
        stable or synthesisable, is RDKit's problem and then a chemist's.
        """
        if self.reachable_volume < 8.0:
            return "no room — anything here would clash"
        bits = []
        if abs(self.electrostatic) > 1.5:
            bits.append("a basic group" if self.electrostatic < 0
                        else "an acidic group")
        if self.donor > 0.25:
            bits.append("an acceptor")
        if self.acceptor > 0.25:
            bits.append("a donor")
        if self.hydrophobic > 1.0 and not bits:
            bits.append("something greasy")
        if not bits:
            bits.append("no strong preference — size and shape decide")
        room = ("tight" if self.reachable_volume < 25
                else "roomy" if self.reachable_volume > 60 else "moderate")
        return "%s; %s" % (", ".join(bits), room)

    def as_dict(self) -> dict:
        return {
            "index": int(self.index),
            "name": self.name,
            "element": self.element,
            "position": [round(float(v), 3) for v in self.position],
            "direction": [round(float(v), 4) for v in self.direction],
            "sasa": round(float(self.sasa), 2),
            "relative_sasa": round(float(self.relative_sasa), 3),
            "reachable_volume": round(float(self.reachable_volume), 1),
            "enclosure": round(float(self.enclosure), 3),
            "electrostatic": round(float(self.electrostatic), 3),
            "donor": round(float(self.donor), 3),
            "acceptor": round(float(self.acceptor), 3),
            "hydrophobic": round(float(self.hydrophobic), 3),
            "toward_partner": round(float(self.toward_partner), 3),
            "partner_contacts": list(self.partner_contacts),
            "field_drivers": list(self.field_drivers),
            "n_cavity_points": int(self.n_cavity_points),
            "score": round(float(self.score), 4),
            "suggests": self.suggests(),
        }


def exit_direction(pos: np.ndarray, i: int, all_pos: np.ndarray) -> np.ndarray:
    """Which way a substituent on atom i would project.

    The vector from the mean of its bonded neighbours to the atom itself.
    Atoms with no neighbour in bonding range fall back to pointing away from
    the molecule's centroid.
    """
    d = np.linalg.norm(all_pos - all_pos[i], axis=1)
    bonded = np.nonzero((d > 0.4) & (d < 1.95))[0]
    if bonded.size:
        return _unit(all_pos[i] - all_pos[bonded].mean(axis=0))
    return _unit(all_pos[i] - all_pos.mean(axis=0))


def growth_points(
    ligand: Ligand,
    anchor_atoms: list[Atom],
    partner_atoms: list[Atom],
    face: PresentedFace,
    cavity: Cavity | None = None,
    partner_residues: list[Residue] | None = None,
    cone_degrees: float = CONE_DEGREES,
    cone_reach: float = CONE_REACH,
    min_relative_sasa: float = MIN_RELATIVE_SASA,
) -> list[GrowthPoint]:
    """Rank the places on the ligand worth growing from, best first."""
    heavy = [a for a in ligand.atoms if a.element.upper() != "H"]
    if not heavy:
        return []

    pos = np.array([(a.x, a.y, a.z) for a in heavy], dtype=np.float64)
    anchor_h = [a for a in anchor_atoms if a.element.upper() != "H"]

    # Exposure against the anchor alone, same convention as the face. An
    # atom the partner happens to be covering is still a growth point --
    # that is the whole mechanism.
    sasa = atom_sasa(heavy, anchor_h + heavy)
    free = atom_sasa(heavy, heavy)
    rel = sasa / np.maximum(free, 1e-9)

    cav_pts = (np.asarray(cavity.points, dtype=np.float64)
               if cavity is not None and len(cavity.points) else
               np.zeros((0, 3)))
    spacing = 0.8
    cos_cut = math.cos(math.radians(cone_degrees))

    # Charged atoms of both proteins, for explaining the field.
    charges = []
    for a in list(anchor_h) + [x for x in partner_atoms
                               if x.element.upper() != "H"]:
        q = CHARGED_ATOMS.get((a.comp_id, a.name))
        if q:
            charges.append((np.array([a.x, a.y, a.z]), q, a))

    part_h = [a for a in partner_atoms if a.element.upper() != "H"]
    part_pos = (np.array([(a.x, a.y, a.z) for a in part_h], dtype=np.float64)
                if part_h else np.zeros((0, 3)))

    out: list[GrowthPoint] = []
    for i, atom in enumerate(heavy):
        if rel[i] < min_relative_sasa:
            continue
        direction = exit_direction(pos, i, pos)

        vol = encl = elec = don = acc = phob = 0.0
        toward = 0.0
        n_pts = 0
        sel = np.zeros(0, dtype=int)
        w = np.zeros(0)
        if len(cav_pts):
            v = cav_pts - pos[i]
            dist = np.linalg.norm(v, axis=1)
            near = dist < cone_reach
            if near.any():
                cosang = (v[near] @ direction) / np.maximum(dist[near], 1e-9)
                sel = np.nonzero(near)[0][cosang >= cos_cut]
                n_pts = int(sel.size)
                if n_pts:
                    vol = n_pts * spacing ** 3
                    # Weight by perpendicular distance from the exit vector.
                    # An unweighted mean is dominated by the mouth of the
                    # cone, which holds the most grid points for purely
                    # geometric reasons and is the part a substituent is
                    # least likely to occupy.
                    along = v[sel] @ direction
                    perp = np.sqrt(np.maximum(
                        dist[sel] ** 2 - along ** 2, 0.0))
                    w = np.exp(-0.5 * (perp / CONE_SIGMA) ** 2)
                    wsum = float(w.sum()) or 1.0

                    def mean(arr):
                        return float((np.asarray(arr)[sel] * w).sum() / wsum)

                    encl = mean(cavity.enclosure)
                    elec = mean(cavity.electrostatic)
                    don = mean(cavity.donor)
                    acc = mean(cavity.acceptor)
                    phob = mean(cavity.hydrophobic)
                    # from_anchor is the share of nearby lining atoms that
                    # belong to the anchor, so 1 - it is the partner's share.
                    toward = 1.0 - mean(cavity.from_anchor)

        # Which charged residues make the field what it is.
        #
        # Computed over the same cavity points with the same weights as
        # `elec` itself, so the parts sum to the whole. An explanation
        # evaluated at a single representative point does not: the potential
        # along one of these vectors can change sign within three angstroms,
        # so a one-point decomposition of a many-point average can come out
        # with the opposite sign to the number it claims to explain, which is
        # worse than offering no explanation at all.
        drivers: dict[str, float] = {}
        if charges and n_pts:
            pts = cav_pts[sel]
            for xyz, q, at in charges:
                d = np.maximum(np.linalg.norm(pts - xyz, axis=1), 1.5)
                per = 332.0 * q / (4.0 * d * d)
                per[d > ELECTROSTATIC_CUTOFF] = 0.0
                contrib = float((per * w).sum() / (w.sum() or 1.0))
                if abs(contrib) < 1e-9:
                    continue
                key = "%s%d" % (at.comp_id, at.seq_id)
                drivers[key] = drivers.get(key, 0.0) + contrib
        # The tail is not negligible and pretending otherwise inverts the
        # story. At 6TD3's best growth point the five largest contributors
        # sum to +1.5 while the field is -3.0: two arginines sit close and
        # push hard, and thirty-odd acidic residues further out quietly
        # outvote them. A top-five list alone would show a positive
        # explanation under a negative number. So whatever is left over is
        # listed as its own row, and the rows always sum to the field.
        ordered = sorted(drivers.items(), key=lambda kv: -abs(kv[1]))
        top = ordered[:5]
        rest = sum(v for _, v in ordered[5:])
        field_drivers = [
            {"residue": k, "contribution": round(v, 2),
             "sign": "+" if v > 0 else "-"}
            for k, v in top if abs(v) > 0.15
        ]
        if abs(rest) > 0.15:
            field_drivers.append({
                "residue": "%d others further out" % max(0, len(ordered) - 5),
                "contribution": round(rest, 2),
                "sign": "+" if rest > 0 else "-",
                "remainder": True,
            })

        contacts: list[str] = []
        if len(part_pos) and partner_residues:
            tip = pos[i] + direction * (cone_reach * 0.6)
            for res in partner_residues:
                if res.comp_id in NON_LIGAND:
                    continue
                xyz = np.array([(a.x, a.y, a.z) for a in res.atoms
                                if a.element.upper() != "H"],
                               dtype=np.float64)
                if not len(xyz):
                    continue
                d = float(np.linalg.norm(xyz - tip, axis=1).min())
                if d < cone_reach:
                    contacts.append((d, "%s%d" % (res.comp_id, res.seq_id)))
            contacts.sort()
            contacts = [c[1] for c in contacts[:8]]

        out.append(GrowthPoint(
            index=i, name=atom.name, element=atom.element,
            position=pos[i], direction=direction,
            sasa=float(sasa[i]), relative_sasa=float(rel[i]),
            reachable_volume=vol, enclosure=encl,
            electrostatic=elec, donor=don, acceptor=acc, hydrophobic=phob,
            toward_partner=toward, partner_contacts=contacts,
            field_drivers=field_drivers, n_cavity_points=n_pts,
        ))

    out.sort(key=lambda g: -g.score)
    return out


def summarise(points: list[GrowthPoint]) -> dict:
    """What the set of growth points says about the site as a whole."""
    if not points:
        return {
            "n": 0,
            "usable": 0,
            "note": (
                "No part of the ligand is exposed enough to grow from. A "
                "fully enclosed ligand cannot be elaborated toward a second "
                "protein, whatever else the interface looks like."
            ),
        }
    usable = [g for g in points if g.reachable_volume >= 8.0]
    toward = [g for g in usable if g.toward_partner >= 0.25]
    return {
        "n": len(points),
        "usable": len(usable),
        "facing_partner": len(toward),
        "best_score": round(points[0].score, 3),
        "total_reachable_volume": round(
            sum(g.reachable_volume for g in usable), 1),
        "note": (
            "%d exposed atoms, %d with room to grow into, %d of those facing "
            "the partner rather than back into the anchor."
            % (len(points), len(usable), len(toward))
        ),
    }
