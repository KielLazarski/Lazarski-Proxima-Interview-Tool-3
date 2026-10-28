"""Statistics for the benchmark.

The headline test is metric 4 from the proposal: at a glued interface, with
the glue deleted, are the contacts near where the glue sits more frustrated
than the rest of the interface?

Two things make that question easy to answer badly, so both are handled here.

First, frustration is not normally distributed, so the comparison is
rank-based (Mann-Whitney) with a rank-biserial effect size, and the p-value
comes from a permutation test rather than a parametric assumption.

Second, and more important: a glue sits in the *core* of an interface, and
core contacts differ from rim contacts for reasons that have nothing to do
with frustration. A raw comparison would partly measure burial. So the same
comparison is also run stratified by local density, which asks the question
within bands of comparable burial.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import stats


@dataclass
class GroupComparison:
    """One frustrated-vs-rest comparison."""

    n_a: int
    n_b: int
    mean_a: float
    mean_b: float
    median_a: float
    median_b: float
    difference: float
    rank_biserial: float
    cohens_d: float
    p_value: float
    test: str

    def as_dict(self) -> dict:
        return {
            "n_proximal": self.n_a,
            "n_distal": self.n_b,
            "mean_proximal": round(self.mean_a, 4),
            "mean_distal": round(self.mean_b, 4),
            "median_proximal": round(self.median_a, 4),
            "median_distal": round(self.median_b, 4),
            "difference": round(self.difference, 4),
            "rank_biserial": round(self.rank_biserial, 4),
            "cohens_d": round(self.cohens_d, 4),
            "p_value": float(f"{self.p_value:.4g}"),
            "test": self.test,
        }


def compare_groups(
    a: np.ndarray,
    b: np.ndarray,
    n_permutations: int = 20000,
    seed: int = 20261026,
) -> GroupComparison:
    """Compare two samples of frustration indices.

    `a` is the glue-proximal group, `b` the rest of the interface. A negative
    difference means the proximal contacts are more frustrated, which is the
    direction the hypothesis predicts.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.size < 3 or b.size < 3:
        return GroupComparison(
            n_a=int(a.size), n_b=int(b.size),
            mean_a=float(a.mean()) if a.size else float("nan"),
            mean_b=float(b.mean()) if b.size else float("nan"),
            median_a=float(np.median(a)) if a.size else float("nan"),
            median_b=float(np.median(b)) if b.size else float("nan"),
            difference=float("nan"), rank_biserial=float("nan"),
            cohens_d=float("nan"), p_value=float("nan"),
            test="insufficient data",
        )

    u, _ = stats.mannwhitneyu(a, b, alternative="two-sided")
    rank_biserial = 1.0 - (2.0 * u) / (a.size * b.size)

    pooled = np.sqrt(
        ((a.size - 1) * a.var(ddof=1) + (b.size - 1) * b.var(ddof=1))
        / max(a.size + b.size - 2, 1)
    )
    d = (a.mean() - b.mean()) / pooled if pooled > 1e-12 else 0.0

    rng = np.random.default_rng(seed)
    perm = stats.permutation_test(
        (a, b),
        lambda x, y, axis=-1: np.mean(x, axis=axis) - np.mean(y, axis=axis),
        permutation_type="independent",
        alternative="two-sided",
        n_resamples=n_permutations,
        random_state=rng,
        vectorized=True,
    )

    return GroupComparison(
        n_a=int(a.size),
        n_b=int(b.size),
        mean_a=float(a.mean()),
        mean_b=float(b.mean()),
        median_a=float(np.median(a)),
        median_b=float(np.median(b)),
        difference=float(a.mean() - b.mean()),
        rank_biserial=float(rank_biserial),
        cohens_d=float(d),
        p_value=float(perm.pvalue),
        test=f"permutation ({n_permutations}) + Mann-Whitney U",
    )


@dataclass
class StratifiedComparison:
    """The same comparison run inside bands of comparable burial."""

    bands: list[dict] = field(default_factory=list)
    pooled_difference: float = float("nan")
    pooled_weight: float = 0.0

    def as_dict(self) -> dict:
        return {
            "bands": self.bands,
            "pooled_difference": (
                round(self.pooled_difference, 4)
                if np.isfinite(self.pooled_difference) else None
            ),
            "note": (
                "Difference in mean frustration, proximal minus distal, "
                "computed inside density bands and pooled by band size. "
                "Controls for a glue sitting in the interface core."
            ),
        }


def compare_stratified(
    values: np.ndarray,
    proximal: np.ndarray,
    covariate: np.ndarray,
    n_bands: int = 4,
) -> StratifiedComparison:
    """Compare proximal against distal within quantile bands of a covariate."""
    values = np.asarray(values, dtype=np.float64)
    proximal = np.asarray(proximal, dtype=bool)
    covariate = np.asarray(covariate, dtype=np.float64)

    if values.size < 20:
        return StratifiedComparison()

    edges = np.quantile(covariate, np.linspace(0, 1, n_bands + 1))
    edges[-1] += 1e-9
    out = StratifiedComparison()
    weighted, total = 0.0, 0.0

    for k in range(n_bands):
        sel = (covariate >= edges[k]) & (covariate < edges[k + 1])
        a = values[sel & proximal]
        b = values[sel & ~proximal]
        entry = {
            "band": k + 1,
            "covariate_range": [round(float(edges[k]), 3),
                                round(float(edges[k + 1]), 3)],
            "n_proximal": int(a.size),
            "n_distal": int(b.size),
        }
        if a.size >= 3 and b.size >= 3:
            diff = float(a.mean() - b.mean())
            entry["mean_proximal"] = round(float(a.mean()), 4)
            entry["mean_distal"] = round(float(b.mean()), 4)
            entry["difference"] = round(diff, 4)
            w = float(a.size + b.size)
            weighted += diff * w
            total += w
        out.bands.append(entry)

    if total > 0:
        out.pooled_difference = weighted / total
        out.pooled_weight = total
    return out


def enrichment_rank(scores: np.ndarray, true_index: int) -> dict:
    """Where the known answer lands in a ranked list. Lower rank is better."""
    scores = np.asarray(scores, dtype=np.float64)
    order = np.argsort(-scores)
    rank = int(np.nonzero(order == true_index)[0][0]) + 1
    return {
        "rank": rank,
        "of": int(scores.size),
        "percentile": round(100.0 * (1.0 - (rank - 1) / max(scores.size - 1, 1)), 2),
    }


def distance_gradient(
    frustration: np.ndarray,
    ligand_distance: np.ndarray,
    cutoffs: tuple[float, ...] = (4.0, 5.0, 6.0, 8.0, 10.0),
) -> dict:
    """Is the effect a gradient, or an artefact of where the cut-off was drawn?

    Splitting interface contacts at one distance and comparing the two halves
    cannot distinguish "frustration rises toward the ligand" from "the
    threshold happened to land somewhere flattering". Two things settle it:
    the rank correlation between frustration and distance, which uses no
    threshold at all, and the same difference recomputed at several cut-offs.

    A positive rho means frustration rises as the ligand is approached, which
    is the direction the hypothesis predicts.
    """
    f = np.asarray(frustration, dtype=np.float64)
    d = np.asarray(ligand_distance, dtype=np.float64)
    if f.size < 10:
        return {"n": int(f.size), "spearman_rho": None, "p_value": None, "sweep": {}}

    rho, p = stats.spearmanr(d, f)
    sweep = {}
    for c in cutoffs:
        a, b = f[d <= c], f[d > c]
        sweep[str(c)] = (
            round(float(a.mean() - b.mean()), 4)
            if a.size > 2 and b.size > 2 else None
        )
    stable = [v for v in sweep.values() if v is not None]
    return {
        "n": int(f.size),
        "spearman_rho": round(float(rho), 4),
        "p_value": float(f"{p:.4g}"),
        "sweep": sweep,
        "sweep_all_same_sign": bool(stable) and (
            all(v < 0 for v in stable) or all(v > 0 for v in stable)
        ),
        "note": (
            "rho is the rank correlation between a contact's frustration and "
            "its distance to the ligand, over interface contacts only. "
            "Positive means frustration rises toward the ligand, with no "
            "threshold involved. The sweep repeats the group difference at "
            "several cut-offs."
        ),
    }
