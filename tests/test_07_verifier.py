"""Unit 07 tests: bounded Rietveld verification.

Exercises core/verification/verifier.py on the two e2e cases (anglesite vs
calcite; GaAsO4 homeotype vs SiO2 quartz), the refinement policy record, and
the additive run-bundle schema upgrade (verification evidence optional; old
bundles keep validating).
"""
import json
import os
import sys
import unittest
from pathlib import Path

import jsonschema

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from benchmarks.eval.sim import ensure_gsasii                        # noqa: E402
from core.report import build_run_bundle, load_bundle_schema         # noqa: E402
from core.verification import (confirmed_by_policy,                  # noqa: E402
                               load_refinement_policy, refine_candidate,
                               verify_case)

POLICY = ROOT / "governance" / "policies" / "refinement-budget.v1.json"
POLICY_SCHEMA = ROOT / "governance" / "schemas" / "refinement_policy.schema.json"
CID_CIF = ROOT / "data" / "unit06" / "input" / "cod"
RELEASE = ROOT / "data" / "catalog" / "releases" / "catalog_0.1.1.json"
WORK = ROOT / "data" / "unit07" / "work"
RES = ROOT / "data" / "unit07" / "results"
FIX = ROOT / "tests" / "fixtures" / "xrdml"
RELEASE_STATE = ROOT / "data" / "unit04" / "results" / "bundles"


class TestRefinementPolicy(unittest.TestCase):
    def test_policy_loads_and_validates(self):
        jsonschema.Draft7Validator.check_schema(json.loads(POLICY_SCHEMA.read_text()))
        policy = load_refinement_policy(str(POLICY))
        self.assertEqual(policy["version"], "1.0")
        self.assertEqual(policy["recipe"], "refinement-verify-v1")
        self.assertIn("Atoms", policy["prohibited"])
        self.assertIn("Mustrain", policy["prohibited"])
        self.assertEqual(policy["confirm"]["max_rwp"], 0.35)


class TestBoundedVerification(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.policy = load_refinement_policy(str(POLICY))
        cls.prm = ensure_gsasii(str(ROOT), str(ROOT / ".vendor" / "GSAS-II"), "")
        cls.release = json.loads(RELEASE.read_text())
        cls.by_id = {e["cod_id"]: e for e in cls.release["entries"]}

    def _candidates(self, true_id, comp_ids):
        out = [(true_id, self.by_id[true_id]["family"], str(CID_CIF / f"{true_id}.cif"))]
        for cid in comp_ids:
            out.append((cid, self.by_id[cid]["family"], str(CID_CIF / f"{cid}.cif")))
        return out

    def test_anglesite_lowest_rwp_and_confirmed(self):
        # case A: cu_PbSO4 -> anglesite 1010950, competitor calcite 1010928
        out = verify_case("cu_PbSO4.xrdml", str(CID_CIF / "1010950.cif"),
                          self.by_id[1010950]["family"],
                          self._candidates(1010950, [1010928]),
                          str(WORK), self.prm, self.policy)
        ranked = out.sorted_results()
        self.assertEqual(ranked[0].family, "PbSO4 (anglesite)")
        self.assertLess(ranked[0].rwp, ranked[1].rwp)
        self.assertLessEqual(ranked[0].rwp, self.policy["confirm"]["max_rwp"])
        self.assertGreaterEqual(out.separation,
                                self.policy["confirm"]["separation_min"])
        self.assertTrue(confirmed_by_policy(out, self.policy))

    def test_gaaso4_lowest_rwp_and_confirmed(self):
        # case B: cu_quartz -> GaAsO4 homeotype 1009000, competitor SiO2
        # alpha-quartz 9009666 (release 0.1.1; the compressed 9012601 from
        # release 0.1.0 is no longer a release entry)
        out = verify_case("cu_quartz.xrdml", str(CID_CIF / "1009000.cif"),
                          self.by_id[1009000]["family"],
                          self._candidates(1009000, [9009666]),
                          str(WORK), self.prm, self.policy)
        ranked = out.sorted_results()
        self.assertEqual(ranked[0].family, "GaAsO4 (quartz homeotype)")
        self.assertLess(ranked[0].rwp, ranked[1].rwp)
        self.assertTrue(confirmed_by_policy(out, self.policy))

    def test_wrong_phase_inverts_confirmation(self):
        # truth is anglesite; if we ask the verifier to start from the WRONG
        # phase as truth table, the outcome must NOT confirm calcite
        out = verify_case("wrong", str(CID_CIF / "1010950.cif"),
                          "CaCO3 (calcite)",
                          self._candidates(1010950, [1010928]),
                          str(WORK), self.prm, self.policy)
        self.assertFalse(confirmed_by_policy(out, self.policy),
                         "calcite must not be confirmed against an anglesite "
                         "observation")

    def test_determinism_across_two_runs(self):
        # same inputs -> bit-identical Rwp/GoF (cached observed pattern)
        tth, yobs, wy = None, None, None
        from core.verification.verifier import make_observed
        tth, yobs, wy = make_observed(str(CID_CIF / "1010950.cif"), str(WORK),
                                      self.prm, self.policy)
        a = refine_candidate(1010950, "PbSO4 (anglesite)",
                             str(CID_CIF / "1010950.cif"), tth, yobs, wy,
                             str(WORK), self.prm, self.policy)
        b = refine_candidate(1010950, "PbSO4 (anglesite)",
                             str(CID_CIF / "1010950.cif"), tth, yobs, wy,
                             str(WORK), self.prm, self.policy)
        self.assertEqual(a.rwp, b.rwp)
        self.assertEqual(a.gof, b.gof)


class TestBundleSchemaAdditive(unittest.TestCase):
    def test_existing_bundle_still_valid(self):
        # pre-upgrade bundles carry no "verification" key; they must still
        # validate against the additive schema
        schema = load_bundle_schema()
        bundles_dir = ROOT / "data" / "unit04" / "results" / "bundles"
        paths = list(bundles_dir.glob("*.bundle.json"))
        self.assertGreaterEqual(len(paths), 3)
        for bpath in paths:
            jsonschema.validate(json.loads(bpath.read_text()), schema)

    def test_verification_evidence_embedded_in_bundle(self):
        # verification outcome from unit 07 embeds into a real run bundle
        report = json.loads((RES / "unit07_report.json").read_text())
        case = report["cases"][0]
        ranked = case["candidates"]
        evidence = {
            "recipe": report["policy"]["recipe"],
            "policy_version": report["policy"]["version"],
            "evidence_level": "fingerprint + refinement",
            "observed_from": case["observed"],
            "confirmed_family": ranked[0]["family"],
            "separation": round(ranked[1]["rwp"] - ranked[0]["rwp"], 6),
            "candidates": ranked,
        }
        from core.calibration import CalibrationRegistry
        from core.hypothesis import load_library, rank_candidates
        from core.ingest import parse_xrdml, sample_fingerprint
        from core.verdict import decide
        lib = json.loads((ROOT / "data" / "candidates" / "library.json").read_text())
        reg = CalibrationRegistry.load(str(ROOT / "data" / "unit03" / "results"
                                           / "registry.json"))
        q = sample_fingerprint(parse_xrdml(str(FIX / "cu_PbSO4.xrdml")))
        res = reg.lookup(parse_xrdml(str(FIX / "cu_PbSO4.xrdml")).instrument)
        ranking = rank_candidates(q, load_library(lib["materials"]),
                                  names={m["id"]: m["name"]
                                         for m in lib["materials"]})
        verdict = decide(ranking, res)
        bundle = build_run_bundle(
            input_path=str(FIX / "cu_PbSO4.xrdml"), raw_sha256="0" * 64,
            sample_name="cu_PbSO4", resolution=res, ranking=ranking,
            verdict=verdict, library_manifest=lib["manifest_sha256"],
            verification=evidence)
        jsonschema.validate(bundle, load_bundle_schema())
        self.assertEqual(bundle["verification"]["evidence_level"],
                         "fingerprint + refinement")
        self.assertEqual(bundle["verification"]["confirmed_family"],
                         "PbSO4 (anglesite)")


if __name__ == "__main__":
    unittest.main(verbosity=2)