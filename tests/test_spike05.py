"""Spike 05 tests: noise-model determinism, family-aware margin semantics,
statistics gate, policy-record consistency, and eval report integrity.
"""
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

from benchmarks.eval.noise import amorphous_pattern, perturb  # noqa: E402
from core.calibration import CalibrationRegistry  # noqa: E402
from core.hypothesis import load_library, rank_candidates  # noqa: E402
from core.ingest import parse_xrdml, sample_fingerprint  # noqa: E402
from core.verdict import (MIN_MARGIN, MIN_TOP_SIMILARITY,  # noqa: E402
                          decide)
from core.verification.verifier import (VERIFY_MIN_PEAK_COUNTS,  # noqa: E402
                                        VERIFY_MIN_PEAK_COUNTS as MIN_PEAK_COUNTS)

FIX = ROOT / "tests" / "fixtures" / "xrdml"
REGISTRY = ROOT / "data" / "spike3" / "results" / "registry.json"
LIBRARY = ROOT / "data" / "candidates" / "library.json"
REPORT = ROOT / "data" / "spike5" / "results" / "spike05_report.json"
POLICY = ROOT / "governance" / "policies" / "m1_fingerprint.policy.json"
POLICY_SCHEMA = ROOT / "governance" / "schemas" / "policy.schema.json"


class TestNoiseModel(unittest.TestCase):
    def test_perturb_is_deterministic_given_seed(self):
        clean = parse_xrdml(str(FIX / "cu_PbSO4.xrdml"))
        a = perturb(clean, counts=1e5, s_mm=0.05, bg_fraction=0.02, seed=7)
        b = perturb(clean, counts=1e5, s_mm=0.05, bg_fraction=0.02, seed=7)
        np.testing.assert_array_equal(a.intensity, b.intensity)

    def test_perturb_is_seed_sensitive(self):
        clean = parse_xrdml(str(FIX / "cu_PbSO4.xrdml"))
        a = perturb(clean, counts=1e5, s_mm=0.05, bg_fraction=0.02, seed=7)
        b = perturb(clean, counts=1e5, s_mm=0.05, bg_fraction=0.02, seed=8)
        self.assertFalse(np.allclose(a.intensity, b.intensity))

    def test_counts_never_negative(self):
        clean = parse_xrdml(str(FIX / "fe_PbSO4.xrdml"))
        for counts in (1e3, 1e6, None):
            p = perturb(clean, counts=counts, s_mm=0.5, bg_fraction=0.2, seed=3)
            self.assertGreaterEqual(p.intensity.min(), 0.0)

    def test_amorphous_pattern_has_broad_humps(self):
        cu_params = parse_xrdml(str(FIX / "cu_PbSO4.xrdml")).instrument
        tth = parse_xrdml(str(FIX / "cu_PbSO4.xrdml")).tth
        pat = amorphous_pattern(tth, counts=1e5, s_mm=0.0, bg_fraction=0.0,
                                seed=1, instrument=cu_params)
        self.assertEqual(pat.sample_name, "amorphous-control")
        self.assertGreater(pat.peak_max, 0.0)


class TestFamilyMargin(unittest.TestCase):
    """Spike 05 core finding: differently-realized entries of the SAME phase
    family must not compete. margin = top vs second overall;
    family_margin = top vs best different-family candidate."""

    @classmethod
    def setUpClass(cls):
        cls.lib = json.loads(LIBRARY.read_text())
        cls.library = load_library(cls.lib["materials"])
        cls.names = {m["id"]: m["name"] for m in cls.lib["materials"]}
        cls.families = {m["id"]: m["phase_family"] for m in cls.lib["materials"]}

    def _rank(self, fixture):
        q = sample_fingerprint(parse_xrdml(str(FIX / fixture)))
        return rank_candidates(q, self.library, names=self.names,
                               families=self.families)

    def test_fe_margin_is_family_aware(self):
        """fe query: raw margin (~0.11, vs the Cu twin) is knife-edge;
        family margin (vs sio2/nacl) must be >> 0.5."""
        r = self._rank("fe_PbSO4.xrdml")
        self.assertEqual(r.ranked[0].phase_family, "PbSO4 (anglesite)")
        self.assertLess(r.margin, 0.15)              # twin-entry competition
        self.assertGreater(r.family_margin, 0.5)     # vs genuinely other phases

    def test_same_material_entries_share_family(self):
        r = self._rank("cu_PbSO4.xrdml")
        fams = {c.material_id: c.phase_family for c in r.ranked}
        self.assertEqual(fams["mat-pbso4"], "PbSO4 (anglesite)")
        self.assertEqual(fams["mat-pbso4-fe"], "PbSO4 (anglesite)")

    def test_family_margin_falls_back_to_plain_margin_without_families(self):
        q = sample_fingerprint(parse_xrdml(str(FIX / "cu_PbSO4.xrdml")))
        r = rank_candidates(q, self.library, names=self.names)
        self.assertEqual(r.family_margin, r.margin)

    def test_verdict_uses_family_margin_in_reasons(self):
        reg = CalibrationRegistry.load(str(REGISTRY))
        r = self._rank("fe_PbSO4.xrdml")
        res = reg.lookup(parse_xrdml(str(FIX / "fe_PbSO4.xrdml")).instrument)
        v = decide(r, res)
        self.assertEqual(v.status, "supported")
        self.assertIn("best other phase family", v.reasons[-1])


class TestVerdictCountingIndependence(unittest.TestCase):
    """The counting-statistics gate lives in the verification stage
    (core.verification.VERIFY_MIN_PEAK_COUNTS), not in the verdict stage:
    identification is position-based (d-space fingerprint).  These tests
    pin that contract: the verdict does not depend on peak counts and never
    reasons about counting statistics (see the core/verdict docstring and
    cli/analyze.py module contract).
    """

    def _bundle_parts(self):
        reg = CalibrationRegistry.load(str(REGISTRY))
        lib = json.loads(LIBRARY.read_text())
        library = load_library(lib["materials"])
        names = {m["id"]: m["name"] for m in lib["materials"]}
        families = {m["id"]: m["phase_family"] for m in lib["materials"]}
        q = sample_fingerprint(parse_xrdml(str(FIX / "cu_PbSO4.xrdml")))
        r = rank_candidates(q, library, names=names, families=families)
        res = reg.lookup(parse_xrdml(str(FIX / "cu_PbSO4.xrdml")).instrument)
        return r, res

    def test_verdict_supports_strong_evidence_regardless_of_peak_counts(self):
        r, res = self._bundle_parts()
        v = decide(r, res)
        self.assertEqual(v.status, "supported")
        self.assertFalse(any("counting" in reason for reason in v.reasons))

    def test_gate_floor_lives_in_verification(self):
        # the calibrated floor (spike-05 envelope L3) guards the refinement
        # stage; it must not gate the identification verdict
        self.assertEqual(VERIFY_MIN_PEAK_COUNTS, 100_000.0)
        self.assertNotIn("peak_max", decide.__code__.co_varnames)


class TestPolicyRecord(unittest.TestCase):
    def test_record_validates_against_policy_schema(self):
        import jsonschema
        rec = json.loads(POLICY.read_text())
        sch = json.loads(POLICY_SCHEMA.read_text())
        jsonschema.validate(rec, sch)

    def test_record_consistent_with_verdict_defaults(self):
        rec = json.loads(POLICY.read_text())["thresholds"]["m1_fingerprint"]
        self.assertEqual(rec["min_top_similarity"], MIN_TOP_SIMILARITY)
        self.assertEqual(rec["min_margin"], MIN_MARGIN)
        self.assertEqual(rec["min_peak_counts"], MIN_PEAK_COUNTS)
        self.assertEqual(rec["margin_kind"], "family_margin")
        self.assertIn("spike05_report.json", rec["evidence"])


class TestEvalReport(unittest.TestCase):
    def test_report_exists_and_passed(self):
        rep = json.loads(REPORT.read_text())
        self.assertTrue(rep["passed"])
        self.assertIn("evidence_pass", rep)
        self.assertIn("fused_pass", rep)
        self.assertIn("statistics_gate", rep)
        self.assertEqual(rep["statistics_gate"]["min_peak_counts"],
                         MIN_PEAK_COUNTS)

    def test_report_fp_rate_zero(self):
        rep = json.loads(REPORT.read_text())
        for a in rep["evidence_pass"]["amorphous_control"]:
            self.assertEqual(a["fp_rate"], 0.0)

    def test_report_markdown_exists(self):
        self.assertTrue((REPORT.parent / "spike05_report.md").exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)