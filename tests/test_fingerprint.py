"""Fingerprint tests: instrument identity, d-space cross-instrument matching,
and material discrimination (see data/unit02/results for measured values).
"""
import numpy as np
import pytest

from core.ingest import (
    InstrumentFingerprint,
    instrument_equivalent,
    instrument_match,
    match_peaks,
    parse_xrdml,
    profile_similarity,
    sample_fingerprint,
    tth_to_d,
)

FIX = "tests/fixtures/xrdml"


def load(name: str):
    return parse_xrdml(f"{FIX}/{name}")


PBS_CU = lambda: sample_fingerprint(load("cu_PbSO4.xrdml"))
PBS_FE = lambda: sample_fingerprint(load("fe_PbSO4.xrdml"))
SIO2_CU = lambda: sample_fingerprint(load("cu_quartz.xrdml"))


# --------------------------------------------------------------- d-spacing math

def test_tth_to_d_sanity():
    # PbSO4 strong reflection: 2th ~29.69 deg for Cu Kalpha1
    d = float(tth_to_d(29.6851, 1.5405))
    assert d == pytest.approx(3.0069, abs=1e-3)


def test_cross_instrument_d_consistency():
    """The same reflection must give the same d on Cu and Fe anodes."""
    dc = PBS_CU().peaks[0].d
    df = PBS_FE().peaks[0].d
    assert abs(dc - df) < 0.01


# ------------------------------------------------------------ same sample match

def test_profile_similarity_self_is_one():
    f = PBS_CU()
    assert profile_similarity(f, f) > 0.999


def test_same_material_across_anodes_matches():
    """PbSO4 measured on Cu and Fe must be recognized as the same material."""
    s = profile_similarity(PBS_CU(), PBS_FE())
    assert s > 0.75, f"cross-anode same-sample similarity too low: {s:.3f}"


def test_peak_counts_are_not_selective_alone():
    """Honest design check: greedy d-matching over-counts across materials;
    the profile metric is the discriminator (finding from unit 02)."""
    n_ss = match_peaks(PBS_CU(), PBS_FE())[0]      # same material
    n_dd = match_peaks(PBS_CU(), SIO2_CU())[0]     # different material
    assert n_ss >= 5
    assert n_dd >= n_ss * 0.4  # just documents the overlap behaviour


# ------------------------------------------------------- material discrimination

def test_different_material_discriminated():
    s_ss = profile_similarity(PBS_CU(), SIO2_CU())
    s_ff = profile_similarity(PBS_FE(), SIO2_CU())
    assert s_ss < 0.30
    assert s_ff < 0.30
    # and better than the cross-anode same-material score
    assert profile_similarity(PBS_CU(), PBS_FE()) > s_ss


# ----------------------------------------------------------------- instruments

def test_instrument_fingerprint_same_setup():
    a = InstrumentFingerprint.from_pattern(load("cu_PbSO4.xrdml"))
    b = InstrumentFingerprint.from_pattern(load("cu_quartz.xrdml"))
    assert instrument_match(a, b)
    assert instrument_equivalent(a.params, b.params)
    assert a.id == b.id


def test_instrument_fingerprint_anode_discrimination():
    a = InstrumentFingerprint.from_pattern(load("cu_PbSO4.xrdml"))
    b = InstrumentFingerprint.from_pattern(load("fe_PbSO4.xrdml"))
    assert not instrument_match(a, b)
    assert not instrument_equivalent(a.params, b.params)
    assert a.id != b.id


def test_fingerprint_id_stable_across_encodings():
    a = InstrumentFingerprint.from_pattern(load("cu_PbSO4.xrdml"))
    b = InstrumentFingerprint.from_pattern(load("base64_PbSO4.xrdml"))
    assert a.id == b.id


def test_peak_line_list_is_json_clean():
    peaks = PBS_CU().peak_line_list()
    assert isinstance(peaks, list) and peaks
    for p in peaks:
        assert isinstance(p["tth"], float)
        assert isinstance(p["d"], float)
        assert isinstance(p["height"], float)
        assert isinstance(p["fwhm"], float)
    assert peaks[0]["d"] == pytest.approx(3.01, abs=0.02)