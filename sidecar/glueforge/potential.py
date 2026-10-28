"""Contact potentials, in the AWSEM functional form.

The functional form here follows the AWSEM contact term as implemented in the
`frustratometer` package:

    theta(r)    = 0.25 * (1 + tanh(eta * (r - r_min))) * (1 + tanh(eta * (r_max - r)))
    thetaII(r)  = same, with r_minII / r_maxII
    rho_i       = sum over j, |i - j| >= 2, of theta(r_ij)
    sigma_w     = 0.25 * (1 - tanh(eta_sigma * (rho_i - rho_0)))
                       * (1 - tanh(eta_sigma * (rho_j - rho_0)))
    sigma_p     = 1 - sigma_w
    E_ij        = -k_contact * ( theta * g_direct
                                 + sigma_w * thetaII * g_water
                                 + sigma_p * thetaII * g_protein )

Distances are between CB atoms, with CA standing in for glycine.

Two sets of gamma tables can be used:

* `awsem` -- the published AWSEM_2015 tables, loaded from the JSON shipped
  with the `frustratometer` package. This is the real method.
* `builtin` -- a transparent stand-in assembled from published residue scales
  (Kyte-Doolittle hydropathy, formal charge, hydrogen-bond capability). It
  exists so the pipeline runs and is testable without the AWSEM file. It is
  NOT the published potential and anything computed with it is a statement
  about this stand-in, not about AWSEM. Results carry the tag so the
  distinction survives into the cache and onto the page.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .residues import (
    AA_INDEX,
    AA_ORDER,
    HB_ACCEPTOR,
    HB_DONOR,
    KD_NORMALISED,
    SIDECHAIN_CHARGE,
)

N_AA = len(AA_ORDER)


@dataclass(frozen=True, slots=True)
class AwsemParams:
    """Switching-function constants. Defaults are the published values."""

    r_min: float = 4.5
    r_max: float = 6.5
    r_min_ii: float = 6.5
    r_max_ii: float = 9.5
    eta: float = 5.0
    eta_sigma: float = 7.0
    rho_0: float = 2.6
    k_contact: float = 4.184
    min_sequence_separation_rho: int = 2
    min_sequence_separation_contact: int = 0
    #: Pairs beyond this are skipped outright; thetaII is ~0 well before it.
    distance_cutoff: float = 12.0
    #: Burial term: three density shells with a tanh switch, published values.
    burial_kappa: float = 4.0
    burial_ro_min: tuple[float, float, float] = (0.0, 3.0, 6.0)
    burial_ro_max: tuple[float, float, float] = (3.0, 6.0, 9.0)


@dataclass
class Potential:
    """The gamma tables plus the switching constants.

    `burial` is 20 x 3: one value per residue type per density shell. It is
    optional because the stand-in has no burial term; when it is present the
    engine adds the burial energy of both residues to each contact, which is
    what AWSEM does.
    """

    name: str
    direct: list[list[float]]
    water: list[list[float]]
    protein: list[list[float]]
    params: AwsemParams
    provenance: str
    burial: list[list[float]] | None = None

    def is_published(self) -> bool:
        return self.name.startswith("awsem")

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "provenance": self.provenance,
            "published": self.is_published(),
            "alphabet": AA_ORDER,
            "params": {
                "r_min": self.params.r_min,
                "r_max": self.params.r_max,
                "r_minII": self.params.r_min_ii,
                "r_maxII": self.params.r_max_ii,
                "eta": self.params.eta,
                "eta_sigma": self.params.eta_sigma,
                "rho_0": self.params.rho_0,
                "k_contact": self.params.k_contact,
                "min_sequence_separation_rho":
                    self.params.min_sequence_separation_rho,
                "min_sequence_separation_contact":
                    self.params.min_sequence_separation_contact,
                "distance_cutoff": self.params.distance_cutoff,
            },
            "direct": self.direct,
            "water": self.water,
            "protein": self.protein,
            "burial": self.burial,
            "burial_params": {
                "kappa": self.params.burial_kappa,
                "ro_min": list(self.params.burial_ro_min),
                "ro_max": list(self.params.burial_ro_max),
            } if self.burial else None,
        }


def _symmetric(mat: list[list[float]]) -> list[list[float]]:
    out = [[0.0] * N_AA for _ in range(N_AA)]
    for i in range(N_AA):
        for j in range(N_AA):
            out[i][j] = 0.5 * (mat[i][j] + mat[j][i])
    return out


# ---------------------------------------------------------------------------
# The transparent stand-in
# ---------------------------------------------------------------------------

#: Weights for the stand-in. Chosen so the three terms land on a comparable
#: scale and the resulting gammas sit in roughly the same numeric range as
#: published contact potentials. They are not fitted to any data.
_W_HYDROPHOBIC = 1.60
_W_ELECTROSTATIC = 0.90
_W_HBOND = 0.35
_W_POLAR_SOLVATION = 1.10


def builtin_potential(params: AwsemParams | None = None) -> Potential:
    """Assemble the stand-in gammas from published residue scales."""
    p = params or AwsemParams()
    direct = [[0.0] * N_AA for _ in range(N_AA)]
    water = [[0.0] * N_AA for _ in range(N_AA)]
    protein = [[0.0] * N_AA for _ in range(N_AA)]

    for a in AA_ORDER:
        ia = AA_INDEX[a]
        na, qa = KD_NORMALISED[a], SIDECHAIN_CHARGE[a]
        da, aa_ = HB_DONOR[a], HB_ACCEPTOR[a]
        for b in AA_ORDER:
            ib = AA_INDEX[b]
            nb, qb = KD_NORMALISED[b], SIDECHAIN_CHARGE[b]
            db, ab = HB_DONOR[b], HB_ACCEPTOR[b]

            # Gamma is the favourability of the pairing: larger is better,
            # because the energy expression negates the whole sum.
            hydrophobic = _W_HYDROPHOBIC * na * nb
            electrostatic = -_W_ELECTROSTATIC * qa * qb
            hbond = _W_HBOND * (da * ab + aa_ * db)
            direct[ia][ib] = hydrophobic + electrostatic + hbond

            # A water-mediated contact rewards polar pairs and screens
            # charge, so the electrostatic term is damped rather than absent.
            polar = _W_POLAR_SOLVATION * (1.0 - na) * (1.0 - nb)
            water[ia][ib] = polar + 0.3 * electrostatic + 0.5 * hbond

            # A buried mediated contact behaves like a weakened direct one.
            protein[ia][ib] = 0.55 * direct[ia][ib]

    return Potential(
        name="builtin-v1",
        direct=_symmetric(direct),
        water=_symmetric(water),
        protein=_symmetric(protein),
        params=p,
        provenance=(
            "Stand-in assembled from Kyte-Doolittle hydropathy, formal "
            "sidechain charge and hydrogen-bond capability. Not a published "
            "potential; use AWSEM_2015 for any result meant to be compared "
            "with the literature."
        ),
    )


# ---------------------------------------------------------------------------
# The published tables
# ---------------------------------------------------------------------------


def load_awsem_json(path: str | Path, params: AwsemParams | None = None) -> Potential:
    """Load the published AWSEM gamma tables.

    The file stores one flat `gamma_array` plus a `segment_definition` that
    says how to slice it, and uses AWSEM's own residue ordering rather than
    an alphabetical one. Both are handled here; see data/PROVENANCE.md.
    """
    path = Path(path)
    with open(path, "r", encoding="utf-8") as fh:
        blob = json.load(fh)

    for key in ("gamma_array", "segment_definition", "alphabet"):
        if key not in blob:
            raise KeyError(
                f"{path.name} is missing {key!r}; expected the AWSEM gamma "
                f"format with keys gamma_array, segment_definition, alphabet. "
                f"Found: {sorted(blob)}"
            )

    flat = [float(v) for v in blob["gamma_array"]]
    segs = {k: tuple(v) for k, v in blob["segment_definition"].items()}
    src_alphabet = [str(a) for a in blob["alphabet"]]
    if len(src_alphabet) != N_AA:
        raise ValueError(f"expected a {N_AA}-letter alphabet, got {len(src_alphabet)}")

    expected = sum(int(np.prod(shape)) for shape in segs.values())
    if expected != len(flat):
        raise ValueError(
            f"{path.name}: segment definition implies {expected} values but "
            f"gamma_array holds {len(flat)}"
        )

    # Segments are laid out in the order the definition lists them.
    cursor = 0
    blocks: dict[str, np.ndarray] = {}
    for name, shape in segs.items():
        size = int(np.prod(shape))
        blocks[name] = np.array(flat[cursor:cursor + size]).reshape(shape)
        cursor += size

    # Reorder from the file's alphabet onto this package's.
    pos = {a: i for i, a in enumerate(src_alphabet)}
    missing = [a for a in AA_ORDER if a not in pos]
    if missing:
        raise ValueError(f"{path.name}: alphabet is missing {missing}")
    order = np.array([pos[a] for a in AA_ORDER])

    def pair_block(name: str) -> list[list[float]]:
        m = blocks[name]
        m = m.reshape(m.shape[-2], m.shape[-1])
        return m[np.ix_(order, order)].tolist()

    direct = pair_block("Direct")
    water = pair_block("Water")
    protein = pair_block("Protein")

    burial = None
    if "Burial" in blocks:
        # Stored as (3 shells, 20 residues); the engine wants (20, 3).
        burial = blocks["Burial"].reshape(3, N_AA)[:, order].T.tolist()

    return Potential(
        name="awsem-2015",
        direct=_symmetric(direct),
        water=_symmetric(water),
        protein=_symmetric(protein),
        params=params or AwsemParams(),
        provenance=(
            f"Published AWSEM-MD (2015) gamma tables from {path.name} "
            f"({blob.get('description', 'no description')}), reordered from "
            f"the AWSEM alphabet onto this package's. Contact and burial "
            f"terms only; the electrostatic term is not implemented."
        ),
        burial=burial,
    )


#: The published tables ship with this package, so the real potential is the
#: default and the stand-in has to be asked for.
BUNDLED_AWSEM = Path(__file__).resolve().parent.parent / "data" / "AWSEM_2015.json"


def load_potential(
    awsem_json: str | Path | None = None,
    params: AwsemParams | None = None,
) -> Potential:
    """Load the published tables, falling back to the stand-in only if absent.

    Pass the string "builtin" to force the stand-in, which is useful for
    checking how much a result depends on the choice of potential.
    """
    if awsem_json == "builtin":
        return builtin_potential(params)
    candidate = Path(awsem_json) if awsem_json else BUNDLED_AWSEM
    if candidate.exists():
        return load_awsem_json(candidate, params)
    return builtin_potential(params)
