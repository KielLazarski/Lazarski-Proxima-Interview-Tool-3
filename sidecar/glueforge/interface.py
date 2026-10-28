"""Working out which protein is the anchor, which is the partner, and where
the glue sits -- without being told.

The roles fall out of the geometry. A glue, by definition, is buried in one
protein's pocket and presents a small face to the other. So the chain that
buries more of the ligand is the anchor (protein A) and the chain that
contacts the presented face is the partner (protein B). That rule also hands
back, for free, the buried-surface asymmetry the proposal rests on.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cif import Ligand, Residue, Structure
from .residues import element_radius

#: A protein atom is in contact with a ligand atom within this distance.
LIGAND_CONTACT = 4.5
#: An interface contact counts as glue-proximal within this distance of any
#: ligand heavy atom. Chosen as one solvation shell beyond direct contact.
GLUE_PROXIMAL = 6.0


@dataclass
class Assignment:
    """One biological copy of a ternary complex, with roles assigned."""

    entry_id: str
    ligand: Ligand
    anchor_chain: str
    partner_chain: str
    anchor_name: str
    partner_name: str
    anchor_contacts: int
    partner_contacts: int
    other_chains: dict[str, int]

    @property
    def asymmetry(self) -> float:
        """Anchor contacts as a share of all ligand-protein contacts."""
        total = self.anchor_contacts + self.partner_contacts
        return self.anchor_contacts / total if total else 0.0

    def as_dict(self) -> dict:
        return {
            "entry": self.entry_id,
            "ligand": self.ligand.comp_id,
            "ligand_chain": self.ligand.chain,
            "ligand_seq": self.ligand.seq_id,
            "ligand_heavy_atoms": self.ligand.heavy_atom_count,
            "anchor_chain": self.anchor_chain,
            "anchor_name": self.anchor_name,
            "partner_chain": self.partner_chain,
            "partner_name": self.partner_name,
            "anchor_ligand_contacts": self.anchor_contacts,
            "partner_ligand_contacts": self.partner_contacts,
            "anchor_share": round(self.asymmetry, 4),
            "other_contacting_chains": self.other_chains,
        }


def _heavy_coords(atoms) -> np.ndarray:
    pts = [(a.x, a.y, a.z) for a in atoms if a.element.upper() != "H"]
    return np.asarray(pts, dtype=np.float64)


def _residue_coords(residues: list[Residue]) -> list[np.ndarray]:
    return [
        np.asarray(
            [(a.x, a.y, a.z) for a in r.atoms if a.element.upper() != "H"],
            dtype=np.float64,
        )
        for r in residues
    ]


def ligand_chain_contacts(
    structure: Structure, ligand: Ligand, cutoff: float = LIGAND_CONTACT
) -> dict[str, int]:
    """Count protein heavy atoms within `cutoff` of the ligand, by chain."""
    lig = _heavy_coords(ligand.atoms)
    if lig.size == 0:
        return {}
    counts: dict[str, int] = {}
    centre = lig.mean(axis=0)
    reach = np.linalg.norm(lig - centre, axis=1).max() + cutoff + 1.0

    for res in structure.residues:
        rep = res.representative()
        if abs(rep.x - centre[0]) > reach + 12:
            continue
        pts = np.asarray(
            [(a.x, a.y, a.z) for a in res.atoms if a.element.upper() != "H"],
            dtype=np.float64,
        )
        if pts.size == 0:
            continue
        if np.linalg.norm(pts.mean(axis=0) - centre) > reach + 12:
            continue
        d = np.sqrt(((pts[:, None, :] - lig[None, :, :]) ** 2).sum(-1))
        hits = int((d < cutoff).sum())
        if hits:
            counts[res.chain] = counts.get(res.chain, 0) + hits
    return counts


def assign_roles(
    structure: Structure,
    ligand_comp: str | None = None,
    cutoff: float = LIGAND_CONTACT,
) -> Assignment:
    """Pick the best copy of the complex and assign anchor / partner roles.

    With several copies in the asymmetric unit, the copy whose ligand has the
    most protein contact is taken; that is the one with the least disorder
    around the site.
    """
    if ligand_comp:
        candidates = [
            l for l in structure.ligands
            if l.comp_id.upper() == ligand_comp.upper()
        ]
    else:
        best = structure.largest_ligand()
        if best is None:
            raise ValueError(f"{structure.entry_id}: no ligand found")
        candidates = [
            l for l in structure.ligands if l.comp_id == best.comp_id
        ]
    if not candidates:
        raise ValueError(f"{structure.entry_id}: ligand {ligand_comp} not found")

    scored = []
    for lig in candidates:
        counts = ligand_chain_contacts(structure, lig, cutoff)
        scored.append((sum(counts.values()), lig, counts))
    scored.sort(key=lambda t: -t[0])
    _, ligand, counts = scored[0]

    if len(counts) < 2:
        raise ValueError(
            f"{structure.entry_id}: ligand {ligand.comp_id} touches "
            f"{len(counts)} chain(s); not a glued interface"
        )

    ranked = sorted(counts.items(), key=lambda kv: -kv[1])
    (anchor, n_anchor), (partner, n_partner) = ranked[0], ranked[1]

    return Assignment(
        entry_id=structure.entry_id,
        ligand=ligand,
        anchor_chain=anchor,
        partner_chain=partner,
        anchor_name=structure.name_of(anchor),
        partner_name=structure.name_of(partner),
        anchor_contacts=n_anchor,
        partner_contacts=n_partner,
        other_chains=dict(ranked[2:]),
    )


def min_distance_to_ligand(
    residues: list[Residue], ligand: Ligand
) -> np.ndarray:
    """Closest approach of each residue's heavy atoms to the ligand."""
    lig = _heavy_coords(ligand.atoms)
    out = np.full(len(residues), np.inf)
    for i, pts in enumerate(_residue_coords(residues)):
        if pts.size == 0:
            continue
        d = np.sqrt(((pts[:, None, :] - lig[None, :, :]) ** 2).sum(-1))
        out[i] = d.min()
    return out


@dataclass
class InterfaceLabels:
    """Per-contact labels for the A-B interface, with the glue deleted."""

    mask_interface: np.ndarray
    mask_proximal: np.ndarray
    contact_min_dist: np.ndarray
    residue_min_dist: np.ndarray

    @property
    def n_interface(self) -> int:
        return int(self.mask_interface.sum())

    @property
    def n_proximal(self) -> int:
        return int((self.mask_interface & self.mask_proximal).sum())


def label_contacts(
    residues: list[Residue],
    index_i: np.ndarray,
    index_j: np.ndarray,
    assignment: Assignment,
    proximal: float = GLUE_PROXIMAL,
) -> InterfaceLabels:
    """Mark which contacts cross the A-B interface and which sit by the glue.

    A contact's distance to the glue is the smaller of its two residues'
    distances, so a contact counts as proximal if either partner reaches the
    ligand.
    """
    chains = np.array([r.chain for r in residues])
    a, b = assignment.anchor_chain, assignment.partner_chain
    ci, cj = chains[index_i], chains[index_j]
    mask_interface = ((ci == a) & (cj == b)) | ((ci == b) & (cj == a))

    res_d = min_distance_to_ligand(residues, assignment.ligand)
    contact_d = np.minimum(res_d[index_i], res_d[index_j])

    return InterfaceLabels(
        mask_interface=mask_interface,
        mask_proximal=contact_d <= proximal,
        contact_min_dist=contact_d,
        residue_min_dist=res_d,
    )
