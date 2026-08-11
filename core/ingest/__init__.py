"""Public API for the ingest layer: XRDML reading + fingerprinting."""
from .fingerprint import (
    InstrumentFingerprint,
    SampleFingerprint,
    find_peaks,
    instrument_equivalent,
    instrument_match,
    match_peaks,
    profile_similarity,
    sample_fingerprint,
    tth_to_d,
)
from .models import InstrumentParams, PowderPattern
from .xrdml import XRDMLError, parse_xrdml

__all__ = [
    "InstrumentFingerprint",
    "InstrumentParams",
    "PowderPattern",
    "SampleFingerprint",
    "XRDMLError",
    "find_peaks",
    "instrument_equivalent",
    "instrument_match",
    "match_peaks",
    "parse_xrdml",
    "profile_similarity",
    "sample_fingerprint",
    "tth_to_d",
]