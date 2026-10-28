"""Growing the molecule: attaching fragments and judging what comes back.

This is the only part of GlueForge that needs a chemistry toolkit, and it is
split so that the part needing one is as small as possible. `score_product`
is numpy and takes coordinates; everything RDKit touches is confined to
`build`. That is not tidiness for its own sake -- RDKit would not install
where this was written, so the scoring half is tested and the building half
is not, and the boundary is where the tested work stops.

**What is scored, and what deliberately is not.**

Products are judged on whether they fit, whether they fill the cavity, and
whether they reach the partner. Those are geometry and they are reliable.

They are *not* ranked on electrostatic complementarity, and the temptation to
do so is the reason for this note. The field here is a Coulomb sum over
formal charges with a distance-dependent dielectric and no explicit solvent.
It is good enough to say a region leans positive or negative -- and it only
became that good after a cutoff was added, because without one it was
reporting each protein's net charge as a property of the site. Turning that
into a number that orders two real compounds would be inventing precision.
So the field character is reported next to each product and left out of the
ranking, the same way cavity volume is measured and kept out of the gate's
verdict.

A chemist reading "fills 40 cubic angstroms of the cavity, reaches four
partner residues, and sits in a region that leans acidic" has something they
can act on. A chemist reading "score 0.73" has a number whose provenance they
would have to take on trust, and it would not deserve it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .cavity import Cavity
from .cif import Atom
from .growth import GrowthPoint

#: Heavy-atom separation below which a product atom is clashing with protein.
#:
#: Set from the benchmark rather than guessed. Across the 14 crystal
#: structures the tightest contact a real bound ligand makes with its protein
#: is 2.47 angstroms, and the median of those minima is 2.88 -- these are
#: hydrogen bonds, not clashes. A threshold of 2.9 would reject four of the
#: fourteen real ligands as impossible, which is a strange standard to hold a
#: substituent to when the molecule it is attached to cannot meet it.
CLASH = 2.4

#: A product atom this close to a cavity point counts as filling it.
FILL = 1.3

#: Reach for counting new partner residues touched by the added atoms.
CONTACT = 4.5

#: A tiny starter library, as SMILES with a single attachment point marked
#: by [*]. Enough to exercise the machinery before anything is downloaded.
#: Chosen to span the obvious directions -- small and greasy, polar donor,
#: polar acceptor, basic, acidic, aromatic -- rather than to be good.
#:
#: Names are the substituent, and specific where the attachment position
#: changes the chemistry: `pyridin-4-yl` rather than `pyridyl`, because a
#: chemist reading the latter cannot tell which isomer this is and the 2-,
#: 3- and 4- isomers behave differently. `carboxylic acid` rather than
#: `carboxylate`, because the SMILES is the neutral acid -- at pH 7 it would
#: be deprotonated, and nothing here models the ligand's own protonation.
#: These names reach the interface and any exported SDF, so a wrong one
#: travels.
STARTER_FRAGMENTS = (
    ("methyl", "[*]C"),
    ("ethyl", "[*]CC"),
    ("cyclopropyl", "[*]C1CC1"),
    ("hydroxyl", "[*]O"),
    ("hydroxymethyl", "[*]CO"),
    ("methoxy", "[*]OC"),
    ("amino", "[*]N"),
    ("methylamino", "[*]NC"),
    ("dimethylamino", "[*]N(C)C"),
    ("acetamido", "[*]NC(C)=O"),
    ("carboxamide", "[*]C(N)=O"),
    ("carboxylic acid", "[*]C(O)=O"),
    ("sulfonamide", "[*]S(N)(=O)=O"),
    ("morpholin-4-yl", "[*]N1CCOCC1"),
    ("piperazin-1-yl", "[*]N1CCNCC1"),
    ("pyrrolidin-1-yl", "[*]N1CCCC1"),
    ("imidazol-5-yl", "[*]c1cnc[nH]1"),
    ("phenyl", "[*]c1ccccc1"),
    ("pyridin-4-yl", "[*]c1ccncc1"),
    ("fluoro", "[*]F"),
    ("chloro", "[*]Cl"),
    ("trifluoromethyl", "[*]C(F)(F)F"),
    ("nitrile", "[*]C#N"),
    ("methylsulfonyl", "[*]S(C)(=O)=O"),
)


class NeedsRDKit(RuntimeError):
    """Raised where a chemistry toolkit is genuinely required."""


def have_rdkit() -> bool:
    try:
        import rdkit  # noqa: F401
        return True
    except Exception:  # noqa: BLE001
        return False


# ---------------------------------------------------------------------------
# Scoring: no toolkit, and tested
# ---------------------------------------------------------------------------


@dataclass
class ProductScore:
    """What was measured about one grown molecule.

    `rank_score` combines only the geometric terms. The field numbers travel
    alongside so a chemist can see what the site wants; they are not folded
    in, because the electrostatics are not good enough to order compounds and
    a combined number would hide that.
    """

    label: str
    n_added: int
    clashes: int
    filled_volume: float
    fill_fraction: float
    partner_contacts: list[str] = field(default_factory=list)
    mean_electrostatic: float = 0.0
    mean_hydrophobic: float = 0.0
    buried_fraction: float = 0.0

    @property
    def viable(self) -> bool:
        return self.clashes == 0 and self.n_added > 0

    @property
    def rank_score(self) -> float:
        if not self.viable:
            return 0.0
        fill = min(1.0, self.filled_volume / 45.0)
        reach = min(1.0, len(self.partner_contacts) / 5.0)
        return float(0.5 * fill + 0.3 * reach + 0.2 * self.buried_fraction)

    def character(self) -> str:
        """What the region these atoms landed in is like."""
        bits = []
        if self.mean_electrostatic < -1.5:
            bits.append("leans negative, so a basic group would suit it")
        elif self.mean_electrostatic > 1.5:
            bits.append("leans positive, so an acidic group would suit it")
        if self.mean_hydrophobic > 1.0:
            bits.append("greasy")
        return "; ".join(bits) or "no strong character"

    def as_dict(self) -> dict:
        return {
            "label": self.label,
            "n_added": int(self.n_added),
            "clashes": int(self.clashes),
            "filled_volume": round(float(self.filled_volume), 1),
            "fill_fraction": round(float(self.fill_fraction), 3),
            "partner_contacts": list(self.partner_contacts),
            "n_partner_contacts": len(self.partner_contacts),
            "mean_electrostatic": round(float(self.mean_electrostatic), 3),
            "mean_hydrophobic": round(float(self.mean_hydrophobic), 3),
            "buried_fraction": round(float(self.buried_fraction), 3),
            "viable": bool(self.viable),
            "rank_score": round(float(self.rank_score), 4),
            "character": self.character(),
        }


def score_product(
    added_xyz: np.ndarray,
    protein: list[Atom],
    partner_residues: list,
    cavity: Cavity | None,
    label: str = "",
    clash: float = CLASH,
) -> ProductScore:
    """Judge the atoms a fragment added. Pure geometry.

    `added_xyz` holds only the new atoms -- the parent ligand is already
    where it was and scoring it again would reward every product equally for
    work none of them did.
    """
    added_xyz = np.asarray(added_xyz, dtype=np.float64).reshape(-1, 3)
    if not len(added_xyz):
        return ProductScore(label=label, n_added=0, clashes=0,
                            filled_volume=0.0, fill_fraction=0.0)

    from scipy.spatial import cKDTree

    heavy = [a for a in protein if a.element.upper() != "H"]
    pxyz = np.array([(a.x, a.y, a.z) for a in heavy], dtype=np.float64)
    clashes = 0
    if len(pxyz):
        d = cKDTree(pxyz).query(added_xyz, k=1)[0]
        clashes = int((d < clash).sum())

    filled = 0
    fill_fraction = 0.0
    elec = phob = 0.0
    buried = 0.0
    if cavity is not None and len(cavity.points):
        cpts = np.asarray(cavity.points, dtype=np.float64)
        tree = cKDTree(cpts)
        # Which cavity points this product occupies.
        hit = set()
        for i, idx in enumerate(tree.query_ball_point(added_xyz, r=FILL)):
            hit.update(idx)
        filled = len(hit)
        fill_fraction = filled / max(1, len(cpts))
        if hit:
            sel = sorted(hit)
            elec = float(np.mean(np.asarray(cavity.electrostatic)[sel]))
            phob = float(np.mean(np.asarray(cavity.hydrophobic)[sel]))
        # How many added atoms are inside the cavity at all, rather than
        # sticking out into solvent where they buy nothing.
        near = tree.query(added_xyz, k=1)[0]
        buried = float((near < FILL).mean())

    contacts: list[str] = []
    for res in partner_residues or []:
        xyz = np.array([(a.x, a.y, a.z) for a in res.atoms
                        if a.element.upper() != "H"], dtype=np.float64)
        if not len(xyz):
            continue
        if float(np.linalg.norm(
                xyz[:, None, :] - added_xyz[None, :, :], axis=2).min()) < CONTACT:
            contacts.append("%s%d" % (res.comp_id, res.seq_id))

    return ProductScore(
        label=label,
        n_added=len(added_xyz),
        clashes=clashes,
        filled_volume=filled * 0.8 ** 3,
        fill_fraction=fill_fraction,
        partner_contacts=contacts[:10],
        mean_electrostatic=elec,
        mean_hydrophobic=phob,
        buried_fraction=buried,
    )


# ---------------------------------------------------------------------------
# Building: needs RDKit, and is not tested
# ---------------------------------------------------------------------------
#
# Everything below is written against RDKit's documented interface on a
# machine where RDKit could not be installed. It is deliberately short and it
# fails loudly: a wrong call raises rather than returning something
# plausible, because a silently mangled molecule that still scores is the one
# failure mode that would waste a chemist's time rather than a minute of it.


def load_fragments(path=None) -> list[tuple[str, str]]:
    """Fragments as (name, SMILES) pairs.

    With no path, the built-in starter set. With one, an SDF read through
    RDKit. Either way each fragment must carry exactly one attachment point,
    marked as a dummy atom.
    """
    if path is None:
        return list(STARTER_FRAGMENTS)
    if not have_rdkit():
        raise NeedsRDKit(
            "Reading an SDF needs RDKit. Without it the built-in starter set "
            "is available by passing no path."
        )
    from rdkit import Chem

    out = []
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    for i, mol in enumerate(supplier):
        if mol is None:
            continue
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        if len(dummies) != 1:
            # Skipped rather than guessed. Choosing an attachment point for
            # someone else's fragment is a chemistry decision and not one
            # this code is entitled to make.
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else "frag%d" % i
        out.append((name, Chem.MolToSmiles(mol)))
    if not out:
        raise ValueError(
            "No usable fragments in %s. Each needs exactly one attachment "
            "point, written as a dummy atom ([*] in SMILES)." % path
        )
    return out


def build(
    ligand_block: str,
    growth: GrowthPoint,
    fragment_smiles: str,
    seed: int = 0xC0FFEE,
):
    """Attach one fragment at one growth point and embed the product.

    `ligand_block` is the parent ligand as a molblock, so its coordinates
    come in and the core stays where the crystal put it. The product is
    embedded with every parent atom constrained, which is the whole point:
    a product whose core has drifted is not a product of this ligand.

    Returns an RDKit Mol with 3D coordinates, and the indices of the atoms
    that were added.
    """
    if not have_rdkit():
        raise NeedsRDKit("Building molecules needs RDKit.")

    from rdkit import Chem
    from rdkit.Chem import AllChem

    parent = Chem.MolFromMolBlock(ligand_block, removeHs=False)
    if parent is None:
        raise ValueError("could not read the ligand as a molblock")

    frag = Chem.MolFromSmiles(fragment_smiles)
    if frag is None:
        raise ValueError("could not read fragment %r" % fragment_smiles)

    dummies = [a.GetIdx() for a in frag.GetAtoms() if a.GetAtomicNum() == 0]
    if len(dummies) != 1:
        raise ValueError(
            "fragment %r has %d attachment points; it needs exactly one"
            % (fragment_smiles, len(dummies))
        )
    dummy = dummies[0]
    neighbours = frag.GetAtomWithIdx(dummy).GetNeighbors()
    if len(neighbours) != 1:
        raise ValueError("the attachment point is not singly bonded")
    frag_attach = neighbours[0].GetIdx()

    n_parent = parent.GetNumAtoms()
    combo = Chem.RWMol(Chem.CombineMols(parent, frag))

    # The growth point indexes heavy atoms of the ligand as the rest of the
    # pipeline sees them; map that onto this molecule's atom order.
    heavy_idx = [a.GetIdx() for a in parent.GetAtoms() if a.GetAtomicNum() > 1]
    if growth.index >= len(heavy_idx):
        raise ValueError(
            "growth point %d is past the end of this ligand (%d heavy atoms)"
            % (growth.index, len(heavy_idx))
        )
    anchor_idx = heavy_idx[growth.index]

    # A substituent replaces a hydrogen. Without removing one the valence is
    # wrong and sanitisation throws, which is the right behaviour but a
    # confusing error to hand a user.
    hydrogens = [n.GetIdx() for n in combo.GetAtomWithIdx(anchor_idx).GetNeighbors()
                 if n.GetAtomicNum() == 1]
    if not hydrogens:
        raise ValueError(
            "%s has no hydrogen to replace, so nothing can be attached there"
            % growth.name
        )

    combo.AddBond(anchor_idx, n_parent + frag_attach, Chem.BondType.SINGLE)
    for idx in sorted([hydrogens[0], n_parent + dummy], reverse=True):
        combo.RemoveAtom(idx)

    product = combo.GetMol()
    Chem.SanitizeMol(product)

    # Constrain every surviving parent atom to where it already is.
    kept = [i for i in range(n_parent) if i != hydrogens[0]]
    conf = parent.GetConformer()
    coord_map = {}
    for new_i, old_i in enumerate(kept):
        coord_map[new_i] = conf.GetAtomPosition(old_i)

    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = False
    if AllChem.EmbedMolecule(product, coordMap=coord_map,
                             randomSeed=seed) != 0:
        raise ValueError(
            "could not place %s here without moving the core" % fragment_smiles
        )
    AllChem.MMFFOptimizeMolecule(product, maxIters=400)

    added = list(range(len(kept), product.GetNumAtoms()))
    return product, added
