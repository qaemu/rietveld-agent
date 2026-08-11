"""Verdicts: governed decisions from calibration resolution + hypothesis.

M1 fingerprint-only evidence. Thresholds are explicit constants (they will
become policy entries once refinement-backed evidence lands); the
calibration requirement from the analysis policy
(``governance.calibration_requirement: released-only``) is a hard stop:
an unresolved or ambiguous instrument means ``abstain``, never a guess.

Identification is POSITION-BASED: the fingerprint compares d-space peak
positions, so counting statistics do not gate the identification claim.
The counting-statistics gate (spike 05, envelope L3: 1e5 peak counts) now
guards the VERIFICATION stage instead (core.verification.verifier,
VERIFY_MIN_PEAK_COUNTS): refinements against weak data are documented as
statistics-below-gate instead of silently refusing to identify.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

from core.calibration import Resolution, ResolutionStatus
from core.hypothesis import CandidateMatch, HypothesisRanking

#: M1 fingerprint-only thresholds (pre-registered in
#: governance/policies/m1_fingerprint.policy.json -- eval-backed by spike 05;
#: to be promoted into full policy records when refinement-backed evidence lands)
MIN_TOP_SIMILARITY = 0.35
MIN_MARGIN = 0.10                # vs the best OTHER phase family
POLICY_VERSION = "policy/0.1.0"


@dataclass
class Verdict:
    status: str                     # "supported" | "abstain"
    primary: Optional[CandidateMatch] = None
    reasons: List[str] = field(default_factory=list)
    thresholds: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {"status": self.status,
                "primary": self.primary.to_dict() if self.primary else None,
                "reasons": list(self.reasons),
                "thresholds": dict(self.thresholds)}


def decide(ranking: HypothesisRanking, resolution: Resolution,
           min_top_similarity: float = MIN_TOP_SIMILARITY,
           min_margin: float = MIN_MARGIN,
           policy_version: str = POLICY_VERSION) -> Verdict:
    thresholds = {"min_top_similarity": min_top_similarity,
                  "min_margin": min_margin,
                  "policy_version": policy_version}
    reasons: List[str] = []

    # --- hard stop: instrument must be uniquely resolved to a released
    # calibration before any further analysis may proceed
    if resolution.status != ResolutionStatus.RESOLVED:
        reasons.append(f"calibration not resolved ({resolution.status.value}): "
                       f"{resolution.reason}")
        return Verdict(status="abstain", reasons=reasons, thresholds=thresholds)

    if not ranking.ranked:
        reasons.append("no candidates ranked")
        return Verdict(status="abstain", reasons=reasons, thresholds=thresholds)

    # NOTE: no counting-statistics gate here -- the fingerprint is a
    # d-space position signature, so identification is counts-independent
    # (the gate lives in the verification stage, where refinements against
    # weak data actually depend on counting statistics).

    top = ranking.ranked[0]
    if top.similarity < min_top_similarity:
        reasons.append(
            f"top candidate similarity {top.similarity:.3f} < {min_top_similarity}")
        return Verdict(status="abstain", reasons=reasons, thresholds=thresholds)

    margin = ranking.family_margin
    if margin < min_margin:
        reasons.append(
            f"family margin {margin:.3f} < {min_margin} "
            f"(top vs best other-phase-family candidate "
            f"{_best_other(top, ranking)})")
        return Verdict(status="abstain", reasons=reasons, thresholds=thresholds)

    reasons.append(f"calibration resolved ({resolution.record['id']}); "
                   f"top candidate separated from the best other phase family "
                   f"by margin {margin:.3f}")
    return Verdict(status="supported", primary=top,
                   reasons=reasons, thresholds=thresholds)


def _best_other(top: CandidateMatch, ranking: HypothesisRanking) -> str:
    others = [c for c in ranking.ranked[1:] if c.phase_family != top.phase_family]
    if not others:
        return "none"
    return max(others, key=lambda c: c.similarity).material_id