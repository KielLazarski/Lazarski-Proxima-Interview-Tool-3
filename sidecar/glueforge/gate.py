"""The gate: is this interface worth doing chemistry into?

This is a deliverable in its own right. It answers the question without any
chemistry being done, and the answer is often no -- which is the point. Glue
discovery is mostly pair selection, and a tool that says "yes, proceed" to
everything has moved the work rather than done it.

Three principles decide what goes in a verdict.

**A metric has to have earned its place.** Every threshold here comes from
`run_gate.py`, which scores the crystal pose of each benchmark system against
its own hardest decoys -- plausible wrong placements of the same partner, the
ones a docking run would have ranked highest. A metric that cannot separate
those two has no business in a verdict however sensible it sounds, and
cavity volume, which is the first thing anyone reaches for, is one of them.

**Calibrating on positives alone is the easy mistake.** Fourteen solved glue
systems tell you what a glue site looks like and nothing at all about what
one does not look like. Thresholds drawn from them would pass whatever they
were shown while feeling informative. The decoys are what make this a test.

**The verdict says what it is made of.** Every number is shown next to the
range the validated systems occupy and the range the decoys occupy, so a
borderline call is visibly borderline rather than arriving as a word.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .cavity import detect_cavity

#: A metric must separate crystal poses from their hardest decoys at least
#: this well to be allowed into a verdict. 0.5 is a coin flip.
MIN_AUC = 0.70

CALIBRATION = Path(__file__).resolve().parent.parent.parent / "cache" / "gate" / "gate.json"

PASS, BORDERLINE, FAIL = "pass", "borderline", "fail"

#: The grid spacing detect_cavity uses, for turning a point count into a
#: volume. Kept here so the two cannot drift apart.
SPACING = 0.8


def cavity_metrics(anchor, partner, ligand_heavy_xyz) -> dict:
    """Cavity at the ligand site for one placement of the partner."""
    from scipy.spatial import cKDTree

    all_atoms = anchor + partner
    anchor_idx = set(range(len(anchor)))
    centre = ligand_heavy_xyz.mean(axis=0)
    cav = detect_cavity(all_atoms, centre, anchor_idx)

    if not len(cav.points):
        return {
            "volume": 0.0, "enclosure_quality": 0.0, "composite_fraction": 0.0,
            "mean_enclosure": 0.0, "ligand_recovery": 0.0,
        }

    enc = np.asarray(cav.enclosure, dtype=float)
    frac = np.asarray(cav.from_anchor, dtype=float)
    rec = float(
        (cKDTree(cav.points).query(ligand_heavy_xyz, k=1)[0] < 1.2).mean()
    )
    return {
        "volume": float(len(cav.points) * SPACING ** 3),
        # How enclosed the average cavity point is, normalised the way the
        # dashboard reports it.
        "enclosure_quality": float((enc >= 3).mean()),
        "mean_enclosure": float(enc.mean()),
        # A glue pocket is lined by BOTH proteins. One lined almost entirely
        # by the anchor is a pocket in the anchor, not an interface pocket.
        "composite_fraction": float(np.minimum(frac, 1.0 - frac).mean() * 2.0),
        "ligand_recovery": rec,
    }


@dataclass
class Check:
    """One metric, measured and placed against what is known."""

    key: str
    label: str
    unit: str
    value: float
    verdict: str
    auc: float
    native_median: float
    native_min: float
    decoy_median: float
    decoy_p90: float
    direction: int
    note: str = ""

    def as_dict(self) -> dict:
        return {
            "key": self.key,
            "label": self.label,
            "unit": self.unit,
            "value": round(float(self.value), 4),
            "verdict": self.verdict,
            "auc": round(float(self.auc), 3),
            "native_median": round(float(self.native_median), 4),
            "native_min": round(float(self.native_min), 4),
            "decoy_median": round(float(self.decoy_median), 4),
            "decoy_p90": round(float(self.decoy_p90), 4),
            "direction": int(self.direction),
            "note": self.note,
        }


@dataclass
class Verdict:
    verdict: str
    headline: str
    checks: list[Check] = field(default_factory=list)
    excluded: list[dict] = field(default_factory=list)
    reasons: list[str] = field(default_factory=list)
    caveats: list[str] = field(default_factory=list)
    calibrated: bool = True
    #: Weighted agreement across the checks, -1 to +1.
    score: float = 0.0
    #: The two metrics that discriminate best, and what they said.
    strongest: list[dict] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "verdict": self.verdict,
            "headline": self.headline,
            "checks": [c.as_dict() for c in self.checks],
            "excluded": list(self.excluded),
            "reasons": list(self.reasons),
            "caveats": list(self.caveats),
            "calibrated": self.calibrated,
            "n_pass": sum(1 for c in self.checks if c.verdict == PASS),
            "n_fail": sum(1 for c in self.checks if c.verdict == FAIL),
            "n_checks": len(self.checks),
            "score": round(float(self.score), 3),
            "strongest": list(self.strongest),
        }


def load_calibration(path: Path | None = None) -> dict | None:
    p = path or CALIBRATION
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def judge(value: float, spec: dict) -> tuple[str, str]:
    """Place one measurement against where the two distributions part.

    The cut is the Youden point from the calibration, not an endpoint of the
    validated range: taking the weakest validated system as the bar lets one
    unusual system drag it down until the hardest decoys clear it.

    The middle band is the honest one. Where the validated and decoy
    distributions overlap, no threshold can separate them, and saying so is
    more useful than rounding to a verdict on the strength of which side of
    an arbitrary line a value fell.
    """
    sign = spec["direction"]
    v = value * sign
    cut = spec.get("threshold", spec["native_min"]) * sign

    # The overlap zone: from where the weakest validated systems sit to where
    # the decoys reach. Inside it, the evidence does not separate.
    native_p10 = spec.get("native_p10", spec["native_min"]) * sign
    decoy_p90 = spec["decoy_p90"] * sign
    lo, hi = min(native_p10, decoy_p90), max(native_p10, decoy_p90)
    if lo <= v <= hi:
        return BORDERLINE, (
            "inside the range where validated systems and wrong poses "
            "overlap, so this measurement cannot tell them apart"
        )
    if v >= cut:
        return PASS, "clear of where wrong poses reach"
    return FAIL, "on the wrong-pose side of the separating threshold"


def assess(metrics: dict, calibration: dict | None = None,
           face: dict | None = None) -> Verdict:
    """Turn measured metrics into a verdict with its reasons.

    `metrics` holds the same keys `run_gate.cavity_metrics` produces, plus
    `face_burial`. Anything missing is skipped rather than guessed.
    """
    cal = calibration if calibration is not None else load_calibration()

    if not cal or "metrics" not in cal:
        return Verdict(
            verdict=BORDERLINE,
            headline="Not calibrated",
            calibrated=False,
            caveats=[
                "No gate calibration on disk, so there is nothing to compare "
                "these numbers against. Run sidecar/run_gate.py to build it. "
                "The measurements below are real; only the verdict is missing."
            ],
            checks=[],
        )

    checks: list[Check] = []
    excluded: list[dict] = []

    for key, spec in cal["metrics"].items():
        if key not in metrics:
            continue
        if spec.get("auc", 0.0) < MIN_AUC:
            # Measured, shown, and kept out of the verdict. Silently dropping
            # it would hide that the obvious metric is the useless one.
            excluded.append({
                "key": key,
                "label": spec.get("label", key),
                "value": round(float(metrics[key]), 4),
                "auc": round(float(spec.get("auc", 0.0)), 3),
                "why": (
                    "separates right poses from wrong ones at AUC %.2f, which "
                    "is close enough to a coin flip that including it would "
                    "add confidence without adding evidence"
                    % spec.get("auc", 0.0)
                ),
            })
            continue

        v, note = judge(float(metrics[key]), spec)
        checks.append(Check(
            key=key,
            label=spec.get("label", key),
            unit=spec.get("unit", ""),
            value=float(metrics[key]),
            verdict=v,
            auc=float(spec["auc"]),
            native_median=float(spec["native_median"]),
            native_min=float(spec["native_min"]),
            decoy_median=float(spec["decoy_median"]),
            decoy_p90=float(spec["decoy_p90"]),
            direction=int(spec["direction"]),
            note=note,
        ))

    if not checks:
        return Verdict(
            verdict=BORDERLINE,
            headline="Nothing measurable",
            calibrated=True,
            caveats=["None of the calibrated metrics could be computed here."],
        )

    n_fail = sum(1 for c in checks if c.verdict == FAIL)
    n_pass = sum(1 for c in checks if c.verdict == PASS)

    # Counting checks equally would let four mediocre metrics outvote two
    # excellent ones. On the calibration, face burial separates right poses
    # from wrong ones at AUC 0.995 and composite fraction at 0.972, while the
    # enclosure measures sit at 0.75-0.83. Weighting by how much better than
    # chance each one is (AUC - 0.5) makes the strong ones carry the verdict,
    # which is what they earned.
    weight = {c.key: max(0.0, c.auc - 0.5) for c in checks}
    total = sum(weight.values()) or 1.0
    sign = {PASS: 1.0, BORDERLINE: 0.0, FAIL: -1.0}
    score = sum(weight[c.key] * sign[c.verdict] for c in checks) / total

    if score >= 0.60:
        verdict, headline = PASS, "This looks like a glue site"
    elif score <= -0.20:
        verdict, headline = FAIL, "This does not look like a glue site"
    else:
        verdict, headline = BORDERLINE, "Mixed"

    reasons = []
    for c in sorted(checks, key=lambda c: (c.verdict != FAIL, -c.auc)):
        arrow = "above" if c.direction > 0 else "below"
        reasons.append(
            "%s is %.3g, %s the weakest validated system (%.3g) — %s"
            % (c.label, c.value,
               arrow if c.verdict == PASS else ("below" if c.direction > 0 else "above"),
               c.native_min, c.note)
        )

    caveats = [
        "Thresholds come from 14 solved glue systems and the plausible wrong "
        "poses of those same partners. Fourteen is a small set and they are "
        "not a random sample of anything, so treat a borderline call as "
        "'not separated by this evidence' rather than as a weak yes."
    ]
    if excluded:
        caveats.append(
            "Measured but deliberately kept out of the verdict: "
            + ", ".join(e["label"] for e in excluded)
            + ". Each failed to separate right poses from wrong ones."
        )
    if face and not face.get("directional", True):
        caveats.append(
            "The ligand sits in a surface groove rather than a pocket, so the "
            "docking restraint that produced this pose was weaker than usual "
            "and the pose is correspondingly less trustworthy."
        )

    out = Verdict(verdict=verdict, headline=headline, checks=checks,
                  excluded=excluded, reasons=reasons, caveats=caveats)
    out.score = round(float(score), 3)
    out.strongest = sorted(
        ({"label": c.label, "auc": round(c.auc, 3), "verdict": c.verdict}
         for c in checks), key=lambda d: -d["auc"]
    )[:2]
    return out
