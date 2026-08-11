"""RunBundle report assembly (governance/schemas/run_bundle.schema.json).

Every analysis produces one immutable, schema-validated bundle: inputs and
their hashes, controlled versions, environment, plan approval, stage
checkpoints, verdicts, evidence separation, and artifact pointers.
"""
from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from pathlib import Path
from typing import Optional

import jsonschema

from core.calibration import Resolution, ResolutionStatus
from core.hypothesis import HypothesisRanking
from core.verdict import Verdict

_SCHEMA_REL = (Path(__file__).resolve().parent.parent.parent
               / "governance" / "schemas" / "run_bundle.schema.json")


def load_bundle_schema() -> dict:
    return json.loads(_SCHEMA_REL.read_text())


def _sha256s(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def normalize_env() -> dict:
    import numpy
    return {"python": sys.version.split()[0],
            "platform": platform.platform(),
            "package_versions": {"numpy": numpy.__version__}}


def build_run_bundle(
    *,
    input_path: str,
    raw_sha256: str,
    sample_name: Optional[str],
    resolution: Resolution,
    ranking: HypothesisRanking,
    verdict: Verdict,
    policy_version: str = "policy/0.1.0",
    gsasii_version: str = "5.6.3",
    library_manifest: str = "",
    recipe_key: str = "fingerprint-only-m1",
    artifact_path: str = "",
    verification: Optional[dict] = None,
) -> dict:
    """Assemble (and schema-validate) the immutable run bundle."""
    run_id = "run-" + _sha256s(f"{input_path}:{raw_sha256}:{time.time_ns()}")[:14]

    if verdict.status == "supported":
        status = "completed"
    else:
        status = "held"            # governed stop: admin/expert action required

    verdict_items = []
    if verdict.primary is not None and verdict.status == "supported":
        verdict_items.append({
            "phase_family": verdict.primary.phase_family or verdict.primary.name,
            "cod_ids": [],
            "status": "supported",
            "context_only": False,
             "evidence": {"evidence_level": "fingerprint-only",
                          "top_similarity": round(ranking.top_similarity, 6),
                          "margin": round(ranking.family_margin, 6),
                          "raw_margin": round(ranking.margin, 6),
                          "matched_peaks": verdict.primary.matched_peaks,
                          "top_peaks": verdict.primary.top_peaks[:5]},
        })

    plan_text = json.dumps({"recipe_key": recipe_key,
                            "policy_version": policy_version}, sort_keys=True)
    checkpoints = [
        {"stage": "ingest", "gpx_path": "", "metrics": {"parsed": True}},
        {"stage": "calibration", "gpx_path": "",
         "metrics": {"resolution": resolution.status.value,
                     "calibration_id": resolution.record["id"] if resolution.record else None}},
        {"stage": "hypothesis", "gpx_path": "",
         "metrics": {"n_candidates": ranking.n_candidates,
                     "top_similarity": round(ranking.top_similarity, 6),
                     "margin": round(ranking.family_margin, 6),
                     "raw_margin": round(ranking.margin, 6)}},
        {"stage": "verdict", "gpx_path": "",
         "metrics": {"verdict": verdict.status,
                     "reasons": list(verdict.reasons)}},
    ]
    bundle = {
        "schema": "run_bundle/v0",
        "run_id": run_id,
        "status": status,
        "inputs": {
            "raw_data": {"path": input_path, "original_sha256": raw_sha256,
                         "format": "XRDML"},
            "sample_name": sample_name,
            "sample_name_sha256": _sha256s(sample_name) if sample_name else None,
        },
        "controlled_versions": {
            "calibration_id": resolution.record["id"] if resolution.record else None,
            "calibration_sha256": None,
            "catalog_version": library_manifest,
            "policy_version": policy_version,
            "gsasii_version": gsasii_version,
            "vocabulary_version": None,
        },
        "environment": normalize_env(),
        "plan": {
            "plan_sha256": _sha256s(plan_text),
            "recipe_key": recipe_key,
            "approval": {"mode": "policy-delegated", "policy_version": policy_version,
                         "expert": None, "approved_at": None},
        },
        "checkpoints": checkpoints,
        "verdicts": verdict_items,
        "evidence_separation": {
            "measured": "", "deterministic_checks": "", "ai_interpretation": "",
        },
        "artifacts": [{"kind": "report-json", "path": artifact_path}],
    }
    if verification is not None:
        bundle["verification"] = verification
    jsonschema.validate(bundle, load_bundle_schema())
    return bundle