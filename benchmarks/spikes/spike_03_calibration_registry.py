"""Spike 03: calibration registry end-to-end on the golden corpus.

Registers trusted calibrations (Cu + Fe, from the validated PRM files) as
schema-conformant CalibrationRecords, resolves every golden XRDML fixture
through core.ingest fingerprints, and demonstrates the three governed
outcomes: resolved / ambiguous / unknown.

Run:  python spike_03_calibration_registry.py
"""
from __future__ import annotations

import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES_DIR = os.path.join(ROOT, "data", "spike3", "results")
FIX_DIR = os.path.join(ROOT, "tests", "fixtures", "xrdml")

CU_PRM = os.path.join(ROOT, "data", "spike", "input", "INST_XRY.PRM")
FE_PRM = os.path.join(ROOT, "data", "spike2", "input", "INST_FE.PRM")
SPIKE1_EVIDENCE = os.path.join(ROOT, "data", "spike", "results", "spike_report.md")

sys.path.insert(0, ROOT)

from core.calibration import CalibrationRegistry, ResolutionStatus
from core.calibration.prm import parse_prm
from core.ingest import InstrumentParams, parse_xrdml


def approved(reviewed_by: str, evidence_ref: str) -> dict:
    now = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    return {"reviewed_by": reviewed_by, "reviewed_at": now,
            "status": "approved", "evidence_ref": evidence_ref}


def main() -> None:
    os.makedirs(RES_DIR, exist_ok=True)
    report: dict = {"mode": "simulate", "started_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")}

    reg = CalibrationRegistry()

    # --- register trusted calibrations ------------------------------------
    cu_rules = [
        {"field": "tube_anode_material", "expected": "CuKa"},
        {"field": "scan_axis", "expected": "2Theta/Theta"},
    ]
    fe_rules = [
        {"field": "tube_anode_material", "expected": "FeKa"},
        {"field": "scan_axis", "expected": "2Theta/Theta"},
    ]
    reg.add_prm(CU_PRM, name="CNEA-Cu-XRD-default", kind="XYE-DEFAULT",
                state="released", approval=approved("cnca-lab-admin", SPIKE1_EVIDENCE),
                notes="default for plain XYE input")
    reg.add_prm(CU_PRM, name="CNEA-Cu-PANalytical-XPert3", kind="XML-FINGERPRINTED",
                state="released", approval=approved("cnca-lab-admin", SPIKE1_EVIDENCE),
                notes="Cu Kalpha Bragg-Brentano", fingerprint_rules=cu_rules)
    reg.add_prm(FE_PRM, name="CNEA-Fe-PANalytical-XPert3", kind="XML-FINGERPRINTED",
                state="released", approval=approved("cnca-lab-admin", SPIKE1_EVIDENCE),
                notes="Fe Kalpha Bragg-Brentano", fingerprint_rules=fe_rules)

    # print parsed PRM content (audit)
    for label, path in (("cu", CU_PRM), ("fe", FE_PRM)):
        p = parse_prm(path)
        print(f"[prm:{label}] anode={p.anode} wl={tuple(round(w, 6) for w in p.wavelengths)} "
              f"type={p.instrument_type} sha={p.sha256[:12]}")
        report.setdefault("prms", {})[label] = {
            "anode": p.anode, "wavelengths": list(p.wavelengths),
            "instrument_type": p.instrument_type, "sha256": p.sha256,
            "profile_function": p.profile_function,
            "profile_coefficients": list(p.profile_coefficients),
        }

    # --- resolve golden fixtures -------------------------------------------
    cases = [
        ("cu_PbSO4.xrdml", "Cu"),
        ("cu_quartz.xrdml", "Cu"),
        ("fe_PbSO4.xrdml", "Fe"),
    ]
    resolutions: dict = {}
    for fname, expect in cases:
        pattern = parse_xrdml(os.path.join(FIX_DIR, fname))
        res = reg.lookup(pattern)
        ok = res.status == ResolutionStatus.RESOLVED and expect in res.record["name"]
        print(f"[resolve:{fname}] {res.status.value:9s} -> "
              f"{res.record['name'] if res.record else '-'}  ok={ok}")
        resolutions[fname] = {"expect": expect, "status": res.status.value,
                              "record_id": res.record["id"] if res.record else None,
                              "reason": res.reason}
        if not ok:
            raise SystemExit(f"resolution mismatch for {fname}")

    # --- governed outcomes: unknown & ambiguous -----------------------------
    unknown = InstrumentParams(anode="MoKa", wavelengths=(0.7093, 0.7136),
                               scan_axis="2Theta/Theta")
    res_u = reg.lookup(unknown)
    print(f"[resolve:unknown] {res_u.status.value}  ({res_u.reason})")
    resolutions["unknown_MoKa"] = {"expect": "unknown", "status": res_u.status.value,
                                   "reason": res_u.reason}
    assert res_u.status == ResolutionStatus.UNKNOWN

    reg2 = CalibrationRegistry()
    base = reg.get(sorted(reg.ids)[0])  # any released record
    for r in reg.records():
        reg2.add_record(r)
    alt = dict(base)
    alt["id"] = "cal-ambig-demo"
    alt["name"] = "CNEA-Cu-alt-duplicate"
    reg2.add_record(alt)  # second RELEASED Cu record with same fingerprint
    res_a = reg2.lookup(parse_xrdml(os.path.join(FIX_DIR, "cu_PbSO4.xrdml")))
    print(f"[resolve:ambiguous] {res_a.status.value} matches={res_a.matches}")
    resolutions["ambiguity_demo"] = {"expect": "ambiguous", "status": res_a.status.value,
                                     "matches": res_a.matches}
    assert res_a.status == ResolutionStatus.AMBIGUOUS

    # --- XYE default ---------------------------------------------------------
    dflt = reg.default_xye()
    print(f"[xye-default] {dflt['name'] if dflt else None} id={dflt['id'] if dflt else None}")
    assert dflt is not None and dflt["kind"] == "XYE-DEFAULT"

    # --- persist + roundtrip --------------------------------------------------
    reg_path = reg.save(os.path.join(RES_DIR, "registry.json"))
    reg_rt = CalibrationRegistry.load(reg_path)
    assert sorted(reg_rt.ids) == sorted(reg.ids)
    print(f"[persist] saved {len(reg.ids)} records -> {reg_path}; load roundtrip ok")

    report["registry_path"] = reg_path
    report["records"] = reg.records()
    report["resolutions"] = resolutions
    report["summary"] = {
        "n_records": len(reg.ids),
        "resolved": sum(1 for v in resolutions.values() if v["status"] == "resolved"),
        "unknown": 1, "ambiguous": 1,
        "governance": "released-only; schema calibration/v0",
    }

    with open(os.path.join(RES_DIR, "spike03_report.json"), "w") as fh:
        json.dump(report, fh, indent=2, default=str)

    md = [
        "# Spike 03: calibration registry",
        "",
        f"- records: {len(reg.ids)} (Cu XYE-DEFAULT, Cu XML-FINGERPRINTED, Fe XML-FINGERPRINTED)",
        f"- schema: governance/schemas/calibration.schema.json (calibration/v0), validated on insert",
        f"- policy: calibration_requirement = released-only (analysis_policy.example.json)",
        "",
        "## PRM audit",
        "",
        "| prm | anode | kalpha1 | kalpha2 | profile fn | sha256 |",
        "|---|---|---|---|---|---|",
    ]
    for label, p in report["prms"].items():
        md.append(f"| {label} | {p['anode']} | {p['wavelengths'][0]:.4f} | "
                  f"{p['wavelengths'][1]:.4f} | {p['profile_function']} | {p['sha256'][:12]} |")
    md += ["", "## Resolution of golden fixtures", "",
           "| fixture | expected | outcome | record |", "|---|---|---|---|"]
    for fname, r in resolutions.items():
        md.append(f"| {fname} | {r['expect']} | {r['status']} | {r.get('record_id') or r.get('matches') or '-'} |")
    md += [
        "",
        "## Findings",
        "- Resolution is physical (anode + kalpha wavelengths + scan axis), not "
        "filename/extension based; scan-grid is soft context (real scans vary).",
        "- ambiguous and unknown are hard stops: the agent must abstain or ask the "
        "administrator to register the instrument; never guess a PRM.",
        "- evidence_ref ties every released calibration to its verification run "
        "(spike01 report) -> audit trail is part of the record.",
        "- XYE-DEFAULT is the single admin-approved fallback for plain XYE input; "
        "metadata-rich XRDML input always resolves via fingerprints.",
        "",
        "## Verdict",
        "- [ ] all 3 golden fixtures resolve to the correct released calibration",
        "- [ ] unknown instrument detected as unknown (Mo variant)",
        "- [ ] duplicate released calibrations detected as ambiguous",
        "- [ ] registry persists and round-trips identically",
        "- [ ] every record validates against calibration.schema.json",
    ]
    with open(os.path.join(RES_DIR, "spike03_report.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    print("wrote", os.path.join(RES_DIR, "spike03_report.json"))


if __name__ == "__main__":
    main()