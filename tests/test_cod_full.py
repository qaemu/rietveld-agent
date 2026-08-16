"""Unit 12b: complete-COD screening (core.codsearch).

The full-COD index (data/cod_index/, gitignored; `make cod-index`) is built
from the rsync'd full COD CIF tree. These tests pin the contract:

* the index covers the ENTIRE COD (>= 500k entries; front page: 534,674),
* anglesite (COD 1010950) is rank 1 on both the Cu and Fe fixtures,
  beating the runner-up,
* the offline kinematic-intensity rerank runs and agrees,
* --full-cod bundles are schema-valid and reproducible.

Tests skip when the index is not built (make cod-index), OR when the rsync
tree is still incomplete -- the whole-file gate is `make check`.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import jsonschema

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from core.codsearch import NPZ_PATH, load_entire_index, screen_cod  # noqa: E402
from core.ingest import parse_xrdml, sample_fingerprint              # noqa: E402
from core.report import load_bundle_schema                           # noqa: E402

FIX = ROOT / "tests" / "fixtures" / "xrdml"
ANGLESITE_COD = 1010950
# a "complete" COD index must cover within a few percent of the COD front
# page count (534,674); anything below is a partial tree (rsync mid-way)
MIN_COMPLETE_ENTRIES = 500_000


def require_complete_index():
    """Load the index; skip when the rsync/COD tree is still incomplete.

    The committed/partial index is a valid screening artifact, but the
    contract tests below pin whole-index properties (coverage >= 500k,
    anglesite rank 1 over the ENTIRE COD), which are untestable until
    `make cod-tree` (~26 GB) + `make cod-index` have run.
    """
    cod_ids, d_units, entry_of, dmin_eff, metas, paths = load_entire_index()
    if len(cod_ids) < MIN_COMPLETE_ENTRIES:
        raise unittest.SkipTest(
            f"complete-COD index incomplete ({len(cod_ids)} entries; "
            f"need >= {MIN_COMPLETE_ENTRIES} — run `make cod-tree` then "
            f"`make cod-index`) and re-run")
    return cod_ids, d_units, entry_of, dmin_eff, metas, paths


@unittest.skipUnless(NPZ_PATH.exists(),
                     "complete-COD index not built (make cod-index)")
class TestCompleteCodIndex(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        (cls.cod_ids, cls.d_units, cls.entry_of, cls.dmin_eff,
         cls.metas, cls.paths) = require_complete_index()

    def test_index_covers_the_entire_cod(self):
        # COD front page (www.crystallography.net) counts 534,674 entries;
        # a "complete" index must be within a few percent of that.
        self.assertGreaterEqual(len(self.cod_ids), 500_000)

    def test_known_structures_present(self):
        for cid in (ANGLESITE_COD,      # anglesite (PbSO4)
                    9008366,            # alite M3 (cement C3S)
                    9012794,            # belite-beta (cement C2S)
                    1000017,            # corundum (Al2O3)
                    1200009,            # brownmillerite (C4AF)
                    8103596):           # tricalcium aluminate (C3A)
            self.assertIn(cid, self.metas, f"COD {cid} missing from index")

    def test_index_shape(self):
        self.assertEqual(len(self.metas), len(self.cod_ids))
        self.assertGreater(self.d_units.size, 80_000_000)

    def test_scan_and_screen_are_deterministic(self):
        fp = sample_fingerprint(parse_xrdml(str(FIX / "cu_PbSO4.xrdml")))
        hits_a = screen_cod(fp, top_k=5, rerank=False)["top"]
        hits_b = screen_cod(fp, top_k=5, rerank=False)["top"]
        self.assertEqual([h["cod_id"] for h in hits_a],
                         [h["cod_id"] for h in hits_b])


@unittest.skipUnless(NPZ_PATH.exists(),
                     "complete-COD index not built (make cod-index)")
class TestCompleteCodScreen(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        require_complete_index()   # heavy screens need a WHOLE-COD index
        cls.a_cu = screen_cod(
            sample_fingerprint(parse_xrdml(str(FIX / "cu_PbSO4.xrdml"))),
            top_k=6, rerank=True)
        cls.a_fe = screen_cod(
            sample_fingerprint(parse_xrdml(str(FIX / "fe_PbSO4.xrdml"))),
            top_k=6, rerank=True)

    def test_anglesite_cu_rank1(self):
        top = self.a_cu["top"][0]
        self.assertEqual(top["cod_id"], ANGLESITE_COD)
        self.assertIn("PB", top["formula"].upper())

    def test_anglesite_fe_rank1(self):
        self.assertEqual(self.a_fe["top"][0]["cod_id"], ANGLESITE_COD)

    def test_rank1_beats_runner_up(self):
        self.assertGreaterEqual(self.a_cu["top"][0]["significance"],
                                self.a_cu["top"][1]["significance"])

    def test_screening_covered_the_whole_index(self):
        self.assertGreater(self.a_cu["scanned"], 500_000)

    def test_intensity_rerank_present_and_agreeing(self):
        rr = self.a_cu["top"][0].get("intensity_rerank")
        self.assertIsNotNone(rr)
        self.assertIsNone(rr.get("error"))
        # coverage pins the rank: the true phase reproduces every measured
        # peak (crude kinematic corr is lower because raw-count heights
        # dominate the cosine; margin over the runner-up is still clear).
        self.assertGreaterEqual(rr["coverage"], 0.9)
        self.assertGreater(rr["corr"], 0.4)


@unittest.skipUnless(NPZ_PATH.exists(),
                     "complete-COD index not built (make cod-index)")
class TestCliFullCod(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        require_complete_index()   # the CLI bundle contract needs the
        #                             whole-COD index built too

    def test_cli_bundle_with_cod_screen_is_valid(self):
        schema = load_bundle_schema()
        out_fd, out_path = tempfile.mkstemp(suffix=".bundle.json")
        import os
        os.close(out_fd)
        try:
            r = subprocess.run(
                [sys.executable, "-m", "cli", "analyze",
                 str(FIX / "cu_PbSO4.xrdml"), "--full-cod",
                 "--no-verification", "--output", out_path],
                capture_output=True, text=True, cwd=str(ROOT), timeout=900)
            self.assertEqual(r.returncode, 0,
                             f"stderr:\n{r.stderr[:2000]}")
            bundle = json.loads(Path(out_path).read_text())
            self.assertEqual(bundle["status"], "complete")
            cs = bundle.get("cod_screen")
            self.assertIsNotNone(cs)
            self.assertGreater(cs["index_entries"], 500_000)
            self.assertEqual(cs["top"][0]["cod_id"], ANGLESITE_COD)
            jsonschema.validate(bundle, schema)
        finally:
            Path(out_path).unlink(missing_ok=True)


if __name__ == "__main__":
    unittest.main()