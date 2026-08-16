"""Unit 04 tests: report/bundle assembly + full analyze pipeline e2e.

Also extends hypothesis tests with library-driven ranking behavior
(same-anode vs cross-anode similarity expectations).
"""
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from core.calibration import CalibrationRegistry   # noqa: E402
from core.hypothesis import load_library, rank_candidates  # noqa: E402
from core.ingest import parse_xrdml, sample_fingerprint   # noqa: E402
from core.report import build_run_bundle, load_bundle_schema  # noqa: E402

FIX = ROOT / "tests" / "fixtures" / "xrdml"
REGISTRY = ROOT / "data" / "unit03" / "results" / "registry.json"
LIBRARY = ROOT / "data" / "candidates" / "library.json"


class TestHypothesisLibrary(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lib_payload = json.loads(LIBRARY.read_text())
        cls.library = load_library(cls.lib_payload["materials"])

    def _rank(self, fixture):
        q = sample_fingerprint(parse_xrdml(str(FIX / fixture)))
        mats = self.lib_payload["materials"]
        names = {m["id"]: m["name"] for m in mats}
        families = {m["id"]: m["phase_family"] for m in mats}
        return rank_candidates(q, self.library, names=names, families=families)

    def test_same_anode_pbso4_ranks_pbso4_first(self):
        r = self._rank("cu_PbSO4.xrdml")
        self.assertEqual(r.ranked[0].material_id, "mat-pbso4")
        self.assertGreater(r.top_similarity, 0.99)
        self.assertGreater(r.family_margin, 0.5)   # vs other phase families
        self.assertGreater(r.margin, 0.09)

    def test_fe_pbso4_ranks_fe_entry_first(self):
        r = self._rank("fe_PbSO4.xrdml")
        self.assertEqual(r.ranked[0].material_id, "mat-pbso4-fe")
        # cross-anode match to the Cu entry must still be high and separated
        self.assertGreater(r.ranked[1].similarity, 0.85)
        self.assertLess(r.ranked[1].similarity, r.ranked[0].similarity)
        self.assertGreater(r.family_margin, 0.5)
        self.assertGreater(r.margin, 0.09)

    def test_cu_quartz_ranks_gaaso4_homeotype_first(self):
        # unit-06 chemistry validation: the M1 'quartz' fixture is GaAsO4
        # (quartz homeotype, COD 1009000); it must top its own family and
        # stay well separated from real SiO2 alpha-quartz (cod-9009666,
        # release 0.1.1) and every other family
        r = self._rank("cu_quartz.xrdml")
        self.assertEqual(r.ranked[0].material_id, "cod-1009000")
        self.assertEqual(r.ranked[0].phase_family, "GaAsO4 (quartz homeotype)")
        self.assertGreater(r.top_similarity, 0.99)
        self.assertGreater(r.family_margin, 0.5)

    def test_library_is_catalog_backed(self):
        ids = {m["id"] for m in self.lib_payload["materials"]}
        self.assertGreaterEqual(len(self.lib_payload["materials"]), 8)
        for want in ("cod-1000041", "cod-9009666", "cod-1010950", "mat-pbso4"):
            self.assertIn(want, ids)
        # superseded mislabelled entries are gone
        self.assertNotIn("mat-nacl", ids)
        self.assertNotIn("mat-sio2", ids)
        self.assertRegex(self.lib_payload["manifest_sha256"], r"^[0-9a-f]{64}$")


class TestRunBundle(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.schema = load_bundle_schema()
        cls.reg = CalibrationRegistry.load(str(REGISTRY))
        cls.lib_payload = json.loads(LIBRARY.read_text())

    def _resolve(self, anode):
        from core.ingest import InstrumentParams
        return self.reg.lookup(InstrumentParams(anode=anode,
                                                wavelengths=(1.5405, 1.5443)))

    def _bundle(self, fixture="cu_PbSO4.xrdml"):
        q = sample_fingerprint(parse_xrdml(str(FIX / fixture)))
        resolution = self._resolve("CuKa")
        names = {m["id"]: m["name"] for m in self.lib_payload["materials"]}
        ranking = rank_candidates(q, load_library(self.lib_payload["materials"]),
                                  names=names)
        from core.verdict import decide
        verdict = decide(ranking, resolution)
        return build_run_bundle(
            input_path=str(FIX / fixture), raw_sha256="0" * 64,
            sample_name=parse_xrdml(str(FIX / fixture)).sample_name,
            resolution=resolution, ranking=ranking, verdict=verdict,
            library_manifest=self.lib_payload["manifest_sha256"],
            artifact_path=str(FIX / f"{fixture}.bundle.json"))

    def test_bundle_validates_against_schema(self):
        bundle = self._bundle()
        self.assertIsNone(bundle.get("boom"))  # sanity: dict is what we built
        from jsonschema import validate
        validate(bundle, self.schema)          # raises on any violation

    def test_bundle_required_shape(self):
        b = self._bundle()
        self.assertEqual(b["schema"], "run_bundle/v0")
        self.assertEqual(b["status"], "completed")
        self.assertRegex(b["plan"]["plan_sha256"], r"^[0-9a-f]{64}$")
        self.assertRegex(b["inputs"]["raw_data"]["original_sha256"],
                         r"^[0-9a-f]{64}$")
        self.assertEqual([c["stage"] for c in b["checkpoints"]],
                         ["ingest", "calibration", "hypothesis", "verdict"])
        self.assertEqual(b["artifacts"][0]["kind"], "report-json")

    def test_unknown_calibration_abstains(self):
        from core.ingest import InstrumentParams
        from core.verdict import decide
        res = self.reg.lookup(InstrumentParams(anode="MoKa",
                                               wavelengths=(0.7093, 0.7136)))
        self.assertEqual(res.status.value, "unknown")
        q = sample_fingerprint(parse_xrdml(str(FIX / "cu_PbSO4.xrdml")))
        empty = rank_candidates(q, {})
        verdict = decide(empty, res)
        self.assertEqual(verdict.status, "abstain")
        bundle = build_run_bundle(
            input_path="unknown.xrdml", raw_sha256="0" * 64,
            sample_name=None, resolution=res, ranking=empty, verdict=verdict,
            library_manifest=self.lib_payload["manifest_sha256"])
        self.assertEqual(bundle["status"], "held")
        self.assertEqual(bundle["verdicts"], [])

    def test_e2e_bundles_exist_for_all_fixtures(self):
        from cli.analyze import analyze
        for fname in ("cu_PbSO4.xrdml", "cu_quartz.xrdml", "fe_PbSO4.xrdml"):
            out = str(FIX / f"{fname}.e2e.json")
            b = analyze(str(FIX / fname), str(REGISTRY), str(LIBRARY), out)
            from jsonschema import validate
            validate(b, self.schema)
            self.assertEqual(b["status"], "completed")
            self.assertEqual(b["verdicts"][0]["status"], "supported")


if __name__ == "__main__":
    unittest.main(verbosity=2)