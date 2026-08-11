"""Verdict stage: governed supported/abstain decisions.

Identification is position-based (d-space fingerprint), so no counting-
statistics gate lives here; that gate guards the verification stage
(core.verification.VERIFY_MIN_PEAK_COUNTS).
"""
from .verdict import (
    MIN_MARGIN,
    MIN_TOP_SIMILARITY,
    POLICY_VERSION,
    Verdict,
    decide,
)

__all__ = ["MIN_MARGIN", "MIN_TOP_SIMILARITY",
           "POLICY_VERSION", "Verdict", "decide"]