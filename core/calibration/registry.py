"""Calibration registry: trusted instrument calibrations (PRM references) with
governed resolution semantics.

Implements the contract in governance/schemas/calibration.schema.json
(CalibrationRecord v0):

* every stored record is JSON-schema validated on insertion;
* resolution is *released-only* (policy: "calibration_requirement": "released-only");
* a lookup returns exactly one of three outcomes:
    RESOLVED  -- exactly one released calibration matches the observation;
    AMBIGUOUS -- more than one released calibration matches (must stop);
    UNKNOWN   -- no released calibration matches (must stop / register);
  so the agent can never silently proceed on a mis-assigned instrument.

Matching is physical: anode material + kalpha wavelengths + scan axis
(metadata-rich XRDML). The scan grid is treated as soft context since real
measurements vary in range/step.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Iterable, Optional

import jsonschema

from core.calibration.prm import PrmContent, parse_prm
from core.ingest import InstrumentParams, PowderPattern

_SCHEMA_REL = Path(__file__).resolve().parent.parent.parent / "governance" / "schemas" / "calibration.schema.json"


class ResolutionStatus(str, Enum):
    RESOLVED = "resolved"
    AMBIGUOUS = "ambiguous"
    UNKNOWN = "unknown"


@dataclass
class Resolution:
    status: ResolutionStatus
    record: Optional[dict] = None
    matches: list = field(default_factory=list)
    reason: str = ""
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return {"status": self.status.value, "record": self.record,
                "matches": self.matches, "reason": self.reason, "notes": self.notes}


def load_schema(path: Optional[str] = None) -> dict:
    p = Path(path or os.environ.get("RIETVELD_AGENT_CALIBRATION_SCHEMA") or _SCHEMA_REL)
    if not p.exists():
        raise FileNotFoundError(f"calibration schema not found: {p}")
    return json.loads(p.read_text())


def _record_id(prm_sha256: str, kind: str) -> str:
    h = hashlib.sha1(f"{prm_sha256}:{kind}".encode()).hexdigest()[:12]
    return f"cal-{h}"


def build_record(prm: PrmContent, *, name: str, kind: str, state: str = "draft",
                 approval: Optional[dict] = None, notes: str = "",
                 fingerprint_rules: Optional[list[dict]] = None,
                 evidence_ref: str = "", created_at: Optional[str] = None) -> dict:
    """Build a CalibrationRecord dict from parsed PRM content."""
    r_k2_k1 = 1.0
    wl = prm.wavelengths
    if len(wl) >= 2:
        r_k2_k1 = wl[1] / wl[0]
    created_at = created_at or time.strftime("%Y-%m-%dT%H:%M:%S%z")
    record = {
        "schema": "calibration/v0",
        "id": _record_id(prm.sha256, kind),
        "name": name,
        "kind": kind,
        "geometry": "Bragg-Brentano",
        "wavelength": {"alpha1": wl[0], "alpha2": wl[1] if len(wl) > 1 else wl[0],
                       "ratio_kalpha2_kalpha1": r_k2_k1},
        "state": state,
        "approval": approval or {
            "reviewed_by": "unset", "reviewed_at": created_at,
            "status": "needs_review", "evidence_ref": evidence_ref,
        },
        "created_at": created_at,
        "content": {"gsasii_instprm": {"path": prm.source_path, "sha256": prm.sha256,
                                       "format": "GSAS-II instrument parameter file"},
                    "notes": notes},
    }
    if fingerprint_rules:
        record["fingerprint"] = {"rules": fingerprint_rules, "accept_incomplete": False}
    return record


def _normalize_anode(value: str) -> str:
    """Normalize anode spellings: 'CuKa' == 'Cu', 'FeKa' == 'Fe', ..."""
    v = value.strip().upper()
    for suffix in ("KALPHA", "KA", "K"):
        if v.endswith(suffix):
            v = v[: -len(suffix)]
            break
    return v


class CalibrationRegistry:
    """Schema-validated, released-only instrument calibration registry."""

    def __init__(self, schema: Optional[dict] = None):
        self._schema = schema or load_schema()
        self._records: dict[str, dict] = {}

    # ------------------------------------------------------------------ records

    def add_record(self, record: dict) -> str:
        jsonschema.validate(record, self._schema)   # raises ValidationError
        rid = record["id"]
        if rid in self._records:
            raise ValueError(f"duplicate calibration id: {rid}")
        self._records[rid] = record
        return rid

    def add_prm(self, path: str, *, name: str, kind: str, state: str = "draft",
                approval: Optional[dict] = None, notes: str = "",
                fingerprint_rules: Optional[list[dict]] = None,
                evidence_ref: str = "") -> str:
        prm = parse_prm(path)
        record = build_record(prm, name=name, kind=kind, state=state,
                              approval=approval, notes=notes,
                              fingerprint_rules=fingerprint_rules,
                              evidence_ref=evidence_ref)
        return self.add_record(record)

    @property
    def ids(self) -> list[str]:
        return sorted(self._records)

    def records(self, state: Optional[str] = None) -> list[dict]:
        recs = list(self._records.values())
        if state:
            recs = [r for r in recs if r.get("state") == state]
        return recs

    def get(self, rid: str) -> Optional[dict]:
        return self._records.get(rid)

    def remove(self, rid: str) -> None:
        self._records.pop(rid, None)

    # ------------------------------------------------------------------ storage

    def save(self, path: str) -> str:
        payload = {"schema": "calibration-registry/v0",
                   "updated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                   "records": [self._records[k] for k in sorted(self._records)]}
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text(json.dumps(payload, indent=2))
        return path

    @classmethod
    def load(cls, path: str, schema: Optional[dict] = None) -> "CalibrationRegistry":
        payload = json.loads(Path(path).read_text())
        reg = cls(schema=schema)
        for rec in payload.get("records", []):
            reg.add_record(rec)
        return reg

    # ---------------------------------------------------------------- resolution

    def _primary_match(self, params: InstrumentParams, record: dict,
                       wl_tol: float = 1e-4) -> bool:
        if record.get("kind") != "XML-FINGERPRINTED":
            return False
        wl = record["wavelength"]
        rec_wl = (wl["alpha1"], wl["alpha2"])
        if len(params.wavelengths) != len(rec_wl):
            return False
        if any(abs(a - b) > wl_tol for a, b in zip(params.wavelengths, rec_wl)):
            return False
        expected_anode = (record.get("fingerprint") or {}).get("rules")
        anode_ok = True
        if expected_anode:
            for rule in expected_anode:
                if rule.get("field") == "tube_anode_material":
                    anode_ok = _normalize_anode(params.anode) == _normalize_anode(
                        str(rule["expected"]))
        return anode_ok

    def lookup(self, params_or_pattern: InstrumentParams | PowderPattern,
               state: str = "released", wl_tol: float = 1e-4) -> Resolution:
        params = params_or_pattern.instrument if isinstance(params_or_pattern, PowderPattern) \
            else params_or_pattern
        matches = [
            r for r in self.records(state=state)
            if self._primary_match(params, r, wl_tol)
        ]
        if len(matches) == 1:
            return Resolution(ResolutionStatus.RESOLVED, record=matches[0],
                              matches=[m["id"] for m in matches],
                              reason=f"unique released calibration for "
                                     f"{params.anode} ({params.wavelengths})")
        if len(matches) > 1:
            return Resolution(ResolutionStatus.AMBIGUOUS, matches=[m["id"] for m in matches],
                              reason=f"{len(matches)} released calibrations match "
                                     f"{params.anode} ({params.wavelengths}); analysis must stop")
        return Resolution(ResolutionStatus.UNKNOWN, matches=[],
                          reason=f"no released calibration matches anode={params.anode} "
                                 f"wavelengths={tuple(round(w, 6) for w in params.wavelengths)}")

    def default_xye(self) -> Optional[dict]:
        """The single administratively approved XYE-DEFAULT calibration."""
        recs = self.records(state="released")
        defaults = [r for r in recs if r.get("kind") == "XYE-DEFAULT"]
        if len(defaults) == 1:
            return defaults[0]
        return None