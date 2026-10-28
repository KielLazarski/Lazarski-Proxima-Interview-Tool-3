"""The frustration engine -- Python reference implementation.

Computes contact-level frustration indices in the AWSEM form. The matching
JavaScript engine in `src/js/03-frustration.js` implements the same thing;
`sidecar/tests/test_cross_engine.py` checks the two agree to within floating
point, which is only meaningful because both draw decoys from the same
deterministic PRNG in the same order.

Sign convention, which is the thing most easily got backwards:

    F = (mean(E_decoy) - E_native) / std(E_decoy)

A native pairing much *lower* in energy than its decoys gives a large positive
F and is minimally frustrated. A native pairing that scores worse than the
decoy average gives negative F and is frustrated. The conventional cut-offs
are F > 0.78 for minimally frustrated and F < -1.0 for highly frustrated.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .cif import Residue
from .potential import N_AA, Potential
from .residues import AA_INDEX

MASK32 = np.uint64(0xFFFFFFFF)

MINIMALLY_FRUSTRATED = 0.78
HIGHLY_FRUSTRATED = -1.0


# ---------------------------------------------------------------------------
# Vectorised PRNG
# ---------------------------------------------------------------------------


def mulberry32_stream(seed: int, n: int) -> np.ndarray:
    """n draws of mulberry32, identical to calling `Mulberry32.next_float()`.

    mulberry32 advances its state by a fixed increment, so the state after k
    calls is seed + k*0x6D2B79F5 and the whole stream can be produced at once
    instead of in a Python loop. `test_prng.py` pins this against the scalar
    implementation.
    """
    if n <= 0:
        return np.zeros(0, dtype=np.float64)
    k = np.arange(1, n + 1, dtype=np.uint64)
    a = (np.uint64(seed) + k * np.uint64(0x6D2B79F5)) & MASK32
    t = a
    t = ((t ^ (t >> np.uint64(15))) * (np.uint64(1) | t)) & MASK32
    inner = ((t ^ (t >> np.uint64(7))) * (np.uint64(61) | t)) & MASK32
    t = (((t + inner) & MASK32) ^ t) & MASK32
    return (((t ^ (t >> np.uint64(14))) & MASK32).astype(np.float64)
            / 4294967296.0)


# ---------------------------------------------------------------------------
# Configuration and results
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class FrustrationConfig:
    n_decoys: int = 1000
    seed: int = 20261026
    #: Pairs closer in sequence than this are not reported as contacts. The
    #: AWSEM default is 1; 2 is used here so backbone neighbours do not fill
    #: the contact list. Inter-chain pairs are always kept.
    min_sep_contact: int = 2
    #: Pairs beyond this distance (CB-CB, angstroms) are not contacts. The
    #: mediated well is essentially zero past 9.5.
    contact_cutoff: float = 9.5
    flavour: str = "configurational"  # or "mutational"


@dataclass
class FrustrationResult:
    index_i: np.ndarray
    index_j: np.ndarray
    distance: np.ndarray
    energy: np.ndarray
    decoy_mean: np.ndarray
    decoy_std: np.ndarray
    frustration: np.ndarray
    density: np.ndarray
    flavour: str
    potential_name: str
    config: FrustrationConfig
    residue_keys: list[tuple[str, int, str]] = field(default_factory=list)
    residue_seq: str = ""

    def __len__(self) -> int:
        return int(self.index_i.size)

    def classify(self) -> np.ndarray:
        """-1 frustrated, 0 neutral, +1 minimally frustrated, per contact."""
        out = np.zeros(self.frustration.shape, dtype=np.int8)
        out[self.frustration > MINIMALLY_FRUSTRATED] = 1
        out[self.frustration < HIGHLY_FRUSTRATED] = -1
        return out

    def summary(self) -> dict:
        cls = self.classify()
        n = max(len(self), 1)
        return {
            "contacts": len(self),
            "flavour": self.flavour,
            "potential": self.potential_name,
            "minimally_frustrated": int((cls == 1).sum()),
            "neutral": int((cls == 0).sum()),
            "highly_frustrated": int((cls == -1).sum()),
            "fraction_frustrated": float((cls == -1).sum() / n),
            "mean_frustration": float(self.frustration.mean()) if len(self) else 0.0,
        }


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


def representative_coords(residues: list[Residue]) -> np.ndarray:
    """CB coordinates, with CA standing in for glycine and truncated sidechains."""
    xyz = np.empty((len(residues), 3), dtype=np.float64)
    for i, r in enumerate(residues):
        a = r.representative()
        xyz[i] = (a.x, a.y, a.z)
    return xyz


def sequence_indices(residues: list[Residue]) -> np.ndarray:
    """Residue identities as indices into the potential's 20-letter alphabet."""
    return np.array(
        [AA_INDEX.get(r.one, AA_INDEX["A"]) for r in residues], dtype=np.int64
    )


def _theta(d: np.ndarray, r_min: float, r_max: float, eta: float) -> np.ndarray:
    return 0.25 * (1.0 + np.tanh(eta * (d - r_min))) * (
        1.0 + np.tanh(eta * (r_max - d))
    )


def local_density(
    dist: np.ndarray, same_chain: np.ndarray, seq_gap: np.ndarray, pot: Potential
) -> np.ndarray:
    """rho_i: the AWSEM local density, summed over the direct-well window."""
    p = pot.params
    th = _theta(dist, p.r_min, p.r_max, p.eta)
    mask = (~same_chain) | (seq_gap >= p.min_sequence_separation_rho)
    np.fill_diagonal(mask, False)
    return (th * mask).sum(axis=1)


# ---------------------------------------------------------------------------
# The engine
# ---------------------------------------------------------------------------


def burial_energy(
    pot: Potential, aa: np.ndarray, rho: np.ndarray, gb: np.ndarray | None
) -> np.ndarray:
    """AWSEM burial term: how well a residue type suits its local density.

    Three density shells, each a tanh window. Returns zero when the potential
    carries no burial table, which is the case for the stand-in.
    """
    if gb is None:
        return np.zeros_like(rho, dtype=np.float64)
    p = pot.params
    out = np.zeros(rho.shape, dtype=np.float64)
    for k in range(3):
        indicator = (
            np.tanh(p.burial_kappa * (rho - p.burial_ro_min[k]))
            + np.tanh(p.burial_kappa * (p.burial_ro_max[k] - rho))
        )
        out += gb[aa, k] * indicator
    return -0.5 * p.k_contact * out


def _pair_energy(
    pot: Potential,
    ai: np.ndarray,
    aj: np.ndarray,
    theta: np.ndarray,
    theta_ii: np.ndarray,
    sigma_w: np.ndarray,
    gd: np.ndarray,
    gw: np.ndarray,
    gp: np.ndarray,
    gb: np.ndarray | None = None,
    rho_i: np.ndarray | None = None,
    rho_j: np.ndarray | None = None,
) -> np.ndarray:
    """Contact energy, plus the burial energy of both residues.

    E = -k * (theta*g_direct + sw*thetaII*g_water + sp*thetaII*g_protein)
        + E_burial(i) + E_burial(j)
    """
    sigma_p = 1.0 - sigma_w
    e = -pot.params.k_contact * (
        theta * gd[ai, aj]
        + sigma_w * theta_ii * gw[ai, aj]
        + sigma_p * theta_ii * gp[ai, aj]
    )
    if gb is not None and rho_i is not None and rho_j is not None:
        e = e + burial_energy(pot, ai, rho_i, gb) + burial_energy(pot, aj, rho_j, gb)
    return e


def compute_frustration(
    residues: list[Residue],
    pot: Potential,
    config: FrustrationConfig | None = None,
) -> FrustrationResult:
    """Contact-level frustration for one set of residues.

    Pass the residues of both chains together to get interface contacts as
    well as intra-chain ones; pass one chain to get that chain alone.
    """
    cfg = config or FrustrationConfig()
    p = pot.params
    n = len(residues)
    if n < 2:
        raise ValueError("need at least two residues")

    xyz = representative_coords(residues)
    aa = sequence_indices(residues)

    diff = xyz[:, None, :] - xyz[None, :, :]
    dist = np.sqrt((diff * diff).sum(-1))

    chain_ids = np.array(
        [hash(r.chain) & 0x7FFFFFFF for r in residues], dtype=np.int64
    )
    same_chain = chain_ids[:, None] == chain_ids[None, :]

    order = np.arange(n)
    seq_gap = np.abs(order[:, None] - order[None, :])

    rho = local_density(dist, same_chain, seq_gap, pot)

    # Contact list: upper triangle, within cut-off, far enough in sequence.
    sep_ok = (~same_chain) | (seq_gap >= cfg.min_sep_contact)
    close = dist < cfg.contact_cutoff
    keep = np.triu(sep_ok & close, k=1)
    ii, jj = np.nonzero(keep)
    m = ii.size
    if m == 0:
        raise ValueError("no contacts found; check the chain selection")

    d = dist[ii, jj]
    theta = _theta(d, p.r_min, p.r_max, p.eta)
    theta_ii = _theta(d, p.r_min_ii, p.r_max_ii, p.eta)
    sw = 0.25 * (1.0 - np.tanh(p.eta_sigma * (rho[ii] - p.rho_0))) * (
        1.0 - np.tanh(p.eta_sigma * (rho[jj] - p.rho_0))
    )

    gd = np.asarray(pot.direct, dtype=np.float64)
    gw = np.asarray(pot.water, dtype=np.float64)
    gp = np.asarray(pot.protein, dtype=np.float64)
    gb = np.asarray(pot.burial, dtype=np.float64) if pot.burial else None

    e_native = _pair_energy(
        pot, aa[ii], aa[jj], theta, theta_ii, sw, gd, gw, gp,
        gb, rho[ii], rho[jj],
    )

    # --- decoys -----------------------------------------------------------
    nd = cfg.n_decoys
    configurational = cfg.flavour == "configurational"
    draws_per = 5 if configurational else 2
    stream = mulberry32_stream(cfg.seed, m * nd * draws_per)
    stream = stream.reshape(m, nd, draws_per)

    da = (stream[:, :, 0] * n).astype(np.int64)
    db = (stream[:, :, 1] * n).astype(np.int64)
    np.clip(da, 0, n - 1, out=da)
    np.clip(db, 0, n - 1, out=db)
    ia, ib = aa[da], aa[db]

    if configurational:
        # Also resample the environment: a distance drawn from the structure's
        # own contact distances, and two densities drawn from its own rho.
        ri = (stream[:, :, 2] * m).astype(np.int64)
        np.clip(ri, 0, m - 1, out=ri)
        d_dec = d[ri]
        qi = (stream[:, :, 3] * n).astype(np.int64)
        qj = (stream[:, :, 4] * n).astype(np.int64)
        np.clip(qi, 0, n - 1, out=qi)
        np.clip(qj, 0, n - 1, out=qj)
        th_dec = _theta(d_dec, p.r_min, p.r_max, p.eta)
        th2_dec = _theta(d_dec, p.r_min_ii, p.r_max_ii, p.eta)
        sw_dec = 0.25 * (1.0 - np.tanh(p.eta_sigma * (rho[qi] - p.rho_0))) * (
            1.0 - np.tanh(p.eta_sigma * (rho[qj] - p.rho_0))
        )
        # Burial follows the resampled environment too.
        rho_i_dec, rho_j_dec = rho[qi], rho[qj]
    else:
        # Mutational: identity changes, the contact's geometry does not.
        th_dec = np.repeat(theta[:, None], nd, axis=1)
        th2_dec = np.repeat(theta_ii[:, None], nd, axis=1)
        sw_dec = np.repeat(sw[:, None], nd, axis=1)
        rho_i_dec = np.repeat(rho[ii][:, None], nd, axis=1)
        rho_j_dec = np.repeat(rho[jj][:, None], nd, axis=1)

    e_dec = _pair_energy(pot, ia, ib, th_dec, th2_dec, sw_dec, gd, gw, gp,
                         gb, rho_i_dec, rho_j_dec)

    mu = e_dec.mean(axis=1)
    sd = e_dec.std(axis=1)
    sd_safe = np.where(sd < 1e-12, 1e-12, sd)
    frustration = (mu - e_native) / sd_safe

    return FrustrationResult(
        index_i=ii,
        index_j=jj,
        distance=d,
        energy=e_native,
        decoy_mean=mu,
        decoy_std=sd,
        frustration=frustration,
        density=rho,
        flavour=cfg.flavour,
        potential_name=pot.name,
        config=cfg,
        residue_keys=[r.key for r in residues],
        residue_seq="".join(r.one for r in residues),
    )
