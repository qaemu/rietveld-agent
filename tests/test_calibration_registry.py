"""Calibration registry tests: PRM parsing, schema compliance, released-only
resolution (resolved/ambiguous/unknown), persistence, XYE default.
"""
import json
import re
from pathlib import Path

import jsonschema
import pytest

from core.calibration import (
    CalibrationRegistry,
    PrmError,
    ResolutionStatus,
    build_record,
    parse_prm,
    load_schema,
)
from core.ingest import InstrumentParams, parse_xrdml

ROOT_REPO = __file__.rsplit("/", 2)[0]  # repo root (tests/ is one level deep)
CU_PRM = f"{ROOT_REPO}/data/spike/input/INST_XRY.PRM"
FE_PRM = f"{ROOT_REPO}/data/spike2/input/INST_FE.PRM"
FIX = f"{ROOT_REPO}/tests/fixtures/xrdml"
EVIDENCE = f"{ROOT_REPO}/data/spike/results/spike_report.md"

APPROVED = {"reviewed_by": "tester", "reviewed_at": "2026-08-09T00:00:00+00:00",
            "status": "approved", "evidence_ref": EVIDENCE}


def pattern(name: str):
    return parse_xrdml(f"{FIX}/{name}")


# -------------------------------------------------------------------- PRM parse

def test_parse_prm_cu():
    p = parse_prm(CU_PRM)
    assert p.anode == "CuKa"
    assert p.wavelengths == pytest.approx((1.5405, 1.5443), abs=1e-4)
    assert p.instrument_type == "PXC"
    assert p.profile_function == "3"
    assert len(p.profile_coefficients) >= 8
    assert re.fullmatch(r"[0-9a-f]{64}", p.sha256)


def test_parse_prm_fe():
    p = parse_prm(FE_PRM)
    assert p.anode == "FeKa"
    assert p.wavelengths == pytest.approx((1.9360, 1.9399), abs=1e-4)


def test_parse_prm_missing_file():
    with pytest.raises(PrmError):
        parse_prm("/nonexistent.prm")


def test_parse_prm_missing_icons(tmp_path):
    bad = tmp_path / "nokey.prm"
    bad.write_text("INS   HTYPE   PXC\nINS  1 IRAD     3\n")
    with pytest.raises(PrmError):
        parse_prm(str(bad))


# -------------------------------------------------------------- schema contract

def test_schema_loaded():
    s = load_schema()
    assert s["$id"].endswith("calibration.schema.json")
    assert s["title"] == "CalibrationRecord"


def test_build_record_validates_against_schema():
    prm = parse_prm(CU_PRM)
    rec = build_record(prm, name="CNEA-Cu-test", kind="XML-FINGERPRINTED",
                       state="released", approval=APPROVED,
                       fingerprint_rules=[{"field": "tube_anode_material",
                                           "expected": "CuKa"}])
    jsonschema.validate(rec, load_schema())  # must not raise
    assert rec["id"].startswith("cal-") and len(rec["id"]) >= 8
    assert rec["wavelength"]["ratio_kalpha2_kalpha1"] == pytest.approx(
        1.5443 / 1.5405, rel=1e-6)


def test_add_record_rejects_schema_violation(tmp_path):
    reg = CalibrationRegistry()
    prm = parse_prm(CU_PRM)
    rec = build_record(prm, name="bad", kind="XML-FINGERPRINTED", state="draft")
    del rec["approval"]  # approval is required by the schema
    with pytest.raises(jsonschema.ValidationError):
        reg.add_record(rec)


def test_duplicate_id_rejected(tmp_path):
    reg = CalibrationRegistry()
    rid = reg.add_prm(CU_PRM, name="a", kind="XYE-DEFAULT", state="draft")
    with pytest.raises(ValueError):
        reg.add_prm(CU_PRM, name="b", kind="XYE-DEFAULT", state="draft")
    assert rid in reg.ids


# ----------------------------------------------------------------- resolution

def _released_registry() -> CalibrationRegistry:
    reg = CalibrationRegistry()
    reg.add_prm(CU_PRM, name="CNEA-Cu-XYE", kind="XYE-DEFAULT", state="released",
                approval=APPROVED)
    reg.add_prm(CU_PRM, name="CNEA-Cu-XML", kind="XML-FINGERPRINTED", state="released",
                approval=APPROVED,
                fingerprint_rules=[{"field": "tube_anode_material", "expected": "CuKa"}])
    reg.add_prm(FE_PRM, name="CNEA-Fe-XML", kind="XML-FINGERPRINTED", state="released",
                approval=APPROVED,
                fingerprint_rules=[{"field": "tube_anode_material", "expected": "FeKa"}])
    return reg


def test_lookup_resolves_cu_and_fe():
    reg = _released_registry()
    r_cu = reg.lookup(pattern("cu_PbSO4.xrdml"))
    assert r_cu.status == ResolutionStatus.RESOLVED
    assert "Cu" in r_cu.record["name"]
    r_q = reg.lookup(pattern("cu_quartz.xrdml"))
    assert r_q.status == ResolutionStatus.RESOLVED and "Cu" in r_q.record["name"]
    r_fe = reg.lookup(pattern("fe_PbSO4.xrdml"))
    assert r_fe.status == ResolutionStatus.RESOLVED and "Fe" in r_fe.record["name"]
    # the XYE-DEFAULT record must never win a fingerprint lookup
    assert r_cu.record["kind"] == "XML-FINGERPRINTED"


def test_lookup_unknown_instrument():
    reg = _released_registry()
    mo = InstrumentParams(anode="MoKa", wavelengths=(0.7093, 0.7136),
                          scan_axis="2Theta/Theta")
    assert reg.lookup(mo).status == ResolutionStatus.UNKNOWN


def test_lookup_unknown_on_plain_params():
    """An InstrumentParams without scan axis still resolves by physics."""
    reg = _released_registry()
    cu = InstrumentParams(anode="Cu", wavelengths=(1.5405, 1.5443))
    assert reg.lookup(cu).status == ResolutionStatus.RESOLVED


def test_released_only_gating():
    reg = _released_registry()
    # a *draft* duplicate (explicit id; content-addressed ids reject dupes)
    # must not influence resolution...
    prm = parse_prm(CU_PRM)
    draft = build_record(prm, name="CNEA-Cu-draft", kind="XML-FINGERPRINTED",
                         state="draft", approval=APPROVED,
                         fingerprint_rules=[{"field": "tube_anode_material",
                                             "expected": "CuKa"}])
    draft["id"] = "cal-draft-copy"
    reg.add_record(draft)
    assert reg.lookup(pattern("cu_PbSO4.xrdml")).status == ResolutionStatus.RESOLVED
    # ...and must not resolve when the released one is removed
    cu_xml = [r for r in reg.records(state="released")
              if r["name"] == "CNEA-Cu-XML"][0]
    reg.remove(cu_xml["id"])
    assert reg.lookup(pattern("cu_PbSO4.xrdml")).status == ResolutionStatus.UNKNOWN


def test_ambiguous_when_two_released_match():
    reg = _released_registry()
    base = [r for r in reg.records(state="released") if r["name"] == "CNEA-Cu-XML"][0]
    alt = dict(base)
    alt["id"], alt["name"] = "cal-ambig-test", "CNEA-Cu-alt"
    reg.add_record(alt)
    res = reg.lookup(pattern("cu_PbSO4.xrdml"))
    assert res.status == ResolutionStatus.AMBIGUOUS
    assert set(res.matches) == {base["id"], "cal-ambig-test"}


# ------------------------------------------------------------- xye default

def test_default_xye_released_single():
    reg = _released_registry()
    d = reg.default_xye()
    assert d is not None and d["kind"] == "XYE-DEFAULT"


def test_default_xye_none_without_released():
    reg = _released_registry()
    d = reg.default_xye()
    reg.remove(d["id"])
    assert reg.default_xye() is None


# --------------------------------------------------------------- persistence

def test_save_load_roundtrip(tmp_path):
    reg = _released_registry()
    path = reg.save(str(tmp_path / "registry.json"))
    reg2 = CalibrationRegistry.load(path)
    assert sorted(reg2.ids) == sorted(reg.ids)
    for rid in reg.ids:
        assert reg2.get(rid) == reg.get(rid)
    # resolutions unchanged after reload
    assert reg2.lookup(pattern("fe_PbSO4.xrdml")).status == ResolutionStatus.RESOLVED


def test_saved_records_validate_against_schema(tmp_path):
    reg = _released_registry()
    payload = json.loads(Path(reg.save(str(tmp_path / "r.json"))).read_text())
    schema = load_schema()
    for rec in payload["records"]:
        jsonschema.validate(rec, schema)