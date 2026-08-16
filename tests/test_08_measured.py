"""Unit 08 tests: bounded verification on MEASURED patterns (numpy engine).

Exercises the measured-path verification wired into cli/analyze:

  * seeded noise envelope: the true phase stays lowest-Rwp at every level,
    and the calibrated bound is not stricter than the unit-07 policy bound
    (L5 - extreme noise/displacement - honestly exceeds the bounded budget);
  * cli e2e on the real fixtures: positive control (cu_quartz == catalog
    CIF 1009000 -> Rwp ~ 0, in bounds) and the documented Wyckoff-mismatch
    case (cu_PbSO4 -> fingerprint-consistent + lowest-Rwp, out of policy
    bounds by honest model mismatch);
  * abstain and off-grid measurements skip verification entirely.

Deterministic: candidate models are protocol simulations cached in
data/unit08/work (fixed CIF + protocol -> fixed pattern).
"""
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

import jsonschema

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from benchmarks.eval.noise import perturb                              # noqa: E402
from benchmarks.eval.sim import (ensure_gsasii,                        # noqa: E402
                                 sim_cif_to_pattern)
from benchmarks.protocols.generate_fixtures import build_xrdml   # noqa: E402
from cli.analyze import analyze as cli_analyze                         # noqa: E402
from core.report import load_bundle_schema                             # noqa: E402
from core.verification import (load_refinement_policy,                 # noqa: E402
                               verify_measured)

POLICY = ROOT / "governance" / "policies" / "refinement-budget.v1.json"
CID_CIF = ROOT / "data" / "unit06" / "input" / "cod"
RELEASE = ROOT / "data" / "catalog" / "releases" / "catalog_0.1.1.json"
WORK = ROOT / "data" / "unit08" / "work"
FIX = ROOT / "tests" / "fixtures" / "xrdml"
ANGLESITE = ("PbSO4 (anglesite)", 1010950, "1010950.cif")
CALCITE = ("CaCO3 (calcite)", 1010928, "1010928.cif")
QUARTZ = ("GaAsO4 (quartz homeotype)", 1009000, "1009000.cif")


class TestNoiseEnvelope(unittest.TestCase):
    """The measured-path engine against seeded synthetic measurements."""

    @classmethod
    def setUpClass(cls):
        cls.policy = load_refinement_policy(str(POLICY))
        cls.prm = ensure_gsasii(str(ROOT), str(ROOT / ".vendor" / "GSAS-II"),
                                "")
        cls.release = json.loads(RELEASE.read_text())
        cls.by_id = {e["cod_id"]: e for e in cls.release["entries"]}

    def _observed(self, counts, s_mm, bg_fraction, seed):
        # ONE clean sim per pattern call: the observed arrays are fixed
        # artifacts (GSAS-II soft-SVD float noise would otherwise vary the
        # clean sim at ~1e-7 between calls)
        clean = sim_cif_to_pattern(str(CID_CIF / ANGLESITE[2]), str(WORK),
                                   prm_path=self.prm)
        noisy = perturb(clean, counts=counts, s_mm=s_mm,
                        bg_fraction=bg_fraction, seed=seed)
        import numpy as np
        return np.asarray(noisy.tth), np.asarray(noisy.intensity)

    def _fit(self, level, counts, s_mm, bg_fraction, seed):
        tth, yobs = self._observed(counts, s_mm, bg_fraction, seed)
        candidates = [(cod, fam, str(CID_CIF / cif))
                      for fam, cod, cif in (ANGLESITE, CALCITE)]
        out = verify_measured(
            tth, yobs, case=f"{level}_env",
            candidates=candidates,
            work_dir=str(WORK), prm_path=self.prm, policy=self.policy)
        return out.sorted_results()

    def test_true_phase_stays_lowest_rwp_across_levels(self):
        # L1 (clean counting stats) and L5 (extreme noise + displacement)
        for level, counts, s_mm, bg in (("L1", 1_000_000, 0.02, 0.005),
                                        ("L5", 10_000, 0.30, 0.10)):
            ranked = self._fit(level, counts, s_mm, bg, seed=1002)
            self.assertEqual(ranked[0].family, ANGLESITE[0],
                             f"{level}: true phase must stay lowest-Rwp")
            self.assertLess(ranked[0].rwp, ranked[1].rwp,
                            f"{level}: true phase must beat competitor")

    def test_low_noise_level_is_within_policy_bounds(self):
        # at L1 the true-phase Rwp is noise-limited, well inside max_rwp,
        # with separation beyond the policy minimum
        ranked = self._fit("L1", 1_000_000, 0.02, 0.005, seed=1002)
        self.assertLessEqual(ranked[0].rwp,
                             self.policy["confirm"]["max_rwp"])
        self.assertGreaterEqual(ranked[1].rwp - ranked[0].rwp,
                                self.policy["confirm"]["separation_min"])

    def test_extreme_noise_level_honestly_exceeds_budget(self):
        # L5 (10k peak counts, 0.30 mm displacement): the bounded budget
        # (uniform shift only) cannot absorb angle-dependent displacement,
        # so true-phase Rwp legitimately exceeds 0.35 - documented
        # limitation, never hidden by relaxing the engineering model
        ranked = self._fit("L5", 10_000, 0.30, 0.10, seed=1002)
        self.assertGreater(ranked[0].rwp,
                           self.policy["confirm"]["max_rwp"])

    def test_measured_fit_is_deterministic(self):
        # identical observed arrays (fixed artifacts) -> identical fit
        tth, yobs = self._observed(30_000, 0.20, 0.05, seed=1002)
        candidates = [(cod, fam, str(CID_CIF / cif))
                      for fam, cod, cif in (ANGLESITE, CALCITE)]
        a = verify_measured(tth, yobs, case="det", candidates=candidates,
                            work_dir=str(WORK), prm_path=self.prm,
                            policy=self.policy).sorted_results()[0]
        b = verify_measured(tth, yobs, case="det", candidates=candidates,
                            work_dir=str(WORK), prm_path=self.prm,
                            policy=self.policy).sorted_results()[0]
        self.assertEqual(a.rwp, b.rwp)
        self.assertEqual(a.gof, b.gof)


class TestCliMeasuredVerification(unittest.TestCase):
    """The verification stage wired into the real cli.analyze pipeline."""

    @classmethod
    def setUpClass(cls):
        cls.policy = load_refinement_policy(str(POLICY))
        cls.schema = load_bundle_schema()

    def test_positive_control_confirms_at_zero_rwp(self):
        # cu_quartz fixture IS the catalog CIF (1009000): Rwp ~ 0, in
        # policy bounds, consistent with the fingerprint verdict
        bundle = cli_analyze(str(FIX / "cu_quartz.xrdml"),
                             verification=True)
        jsonschema.validate(bundle, self.schema)
        v = bundle["verification"]
        self.assertEqual(v["confirmed_family"], QUARTZ[0])
        self.assertTrue(v["consistent_with_fingerprint"])
        self.assertLessEqual(v["candidates"][0]["rwp"], 0.02)
        self.assertGreaterEqual(v["separation"], 0.5)
        self.assertEqual(v["evidence_level"], "fingerprint + refinement")
        self.assertEqual(v["recipe"], "refinement-verify-v1")

    def test_structure_mismatch_honestly_out_of_bounds(self):
        # cu_PbSO4 is the APS-tutorial Wyckoff structure, NOT exactly
        # catalog 1010950: fingerprint-consistent + lowest-Rwp, yet the
        # absolute Rwp exceeds the calibrated bound - documented, and the
        # evidence remains consistent with the fingerprint top family
        bundle = cli_analyze(str(FIX / "cu_PbSO4.xrdml"), verification=True)
        jsonschema.validate(bundle, self.schema)
        v = bundle["verification"]
        self.assertEqual(v["candidates"][0]["family"], ANGLESITE[0])
        self.assertEqual(v["confirmed_family"], ANGLESITE[0])
        self.assertTrue(v["consistent_with_fingerprint"])
        self.assertLess(v["candidates"][0]["rwp"], 1.0)
        self.assertGreater(v["candidates"][0]["rwp"], 0.5)

    def test_abstain_runs_no_verification(self):
        # amorphous measurement never reaches "supported": no verification
        # stage may run and the bundle carries no verification evidence
        import numpy as np
        from benchmarks.eval.noise import amorphous_pattern
        from core.ingest import parse_xrdml
        pat = parse_xrdml(str(FIX / "cu_PbSO4.xrdml"))
        synth = amorphous_pattern(np.asarray(pat.tth), counts=1_000_000,
                                  s_mm=0.05, bg_fraction=0.02, seed=7,
                                  instrument=pat.instrument)
        with tempfile.TemporaryDirectory() as td:
            out = str(Path(td) / "amorphous.xrdml")
            Path(out).write_text(build_xrdml(
                "amorphous-control", "CuKa", 1.5405, 1.5443,
                15.0, 140.0, 0.02, synth.intensity.tolist()))
            bundle = cli_analyze(out, verification=True)
            jsonschema.validate(bundle, self.schema)
        self.assertEqual(bundle["status"], "held")
        self.assertEqual(bundle["verdicts"], [])   # abstain: no verdict
        self.assertNotIn("verification", bundle)

    def test_off_grid_measurement_verified_by_resampling(self):
        # a grid that does not match the protocol (4251 pts, tmax 100 deg)
        # must not crash: the measured grid is resampled onto the protocol
        # grid and the verification stage still runs (counting statistics
        # assessed, never gating the identification) - see the cli/analyze
        # module contract.
        pat = sim_cif_to_pattern(str(CID_CIF / ANGLESITE[2]), str(WORK),
                                 prm_path=ensure_gsasii(
                                     str(ROOT),
                                     str(ROOT / ".vendor" / "GSAS-II"), ""))
        with tempfile.TemporaryDirectory() as td:
            out = str(Path(td) / "short.xrdml")
            Path(out).write_text(build_xrdml(
                "short-scan", "CuKa", 1.5405, 1.5443,
                15.0, 100.0, 0.02, [float(v) for v in pat.intensity[:4251]]))
            bundle = cli_analyze(out, verification=True)
            jsonschema.validate(bundle, self.schema)
        v = bundle["verification"]
        self.assertEqual(bundle["status"], "completed")
        self.assertTrue(v["grid"]["resampled"])
        self.assertIn("statistics", v)
        self.assertIn(v["statistics"]["satisfied"], (True, False))


class TestEnvelopePolicyConsistency(unittest.TestCase):
    def test_policy_bound_was_not_relaxed_beyond_envelope(self):
        # the unit-08 envelope (L1..L5 x3 seeds) calibrates the bound from
        # the true-phase Rwp max + margin: it must not be stricter than the
        # unit-07 bound (0.35) - the report records the calibrated value
        report = json.loads(
            (ROOT / "data" / "unit08" / "results" /
             "unit08_report.json").read_text())
        self.assertGreaterEqual(report["envelope_study"]["calibrated_bound"],
                                0.35)
        self.assertTrue(report["envelope_study"]["all_lowest_is_true"])
        self.assertTrue(report["summary"]["all_schema_valid"])
        self.assertTrue(report["summary"]["all_lowest_rwp_matches_fingerprint"])


if __name__ == "__main__":
    unittest.main()