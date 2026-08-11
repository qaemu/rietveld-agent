"""Spike 06: COD-backed catalogue release + rebuilt runtime library.

Checks the artefacts written by benchmarks/spikes/spike_06_catalog_build.py:
  data/catalog/schemas/catalog_release.v0.json
  data/catalog/releases/catalog_0.1.1.json
  data/candidates/library.json
  data/spike6/results/spike06_report.json

Pinned criteria (from the spike 06 gate):
  * release: >= 8 entries, all schema-valid, CIF sha256 verified on disk,
    cif_validation_rate >= 0.98 (12/12 after the pyrite close-out), >= 5
    newly fetched CODs, pyrite present as a validated FeS2/Pa-3 entry (COD
    text search "pyrite" is contaminated; formula search "FeS2" resolves),
    manifest deterministic.
  * library: catalog-backed cod-* entries + legacy PbSO4 fixtures kept,
    M1 mislabels (mat-nacl = Ag0.5Bi0.5S, mat-sio2 = GaAsO4) superseded.
"""

import hashlib
import json
import os
import unittest

import jsonschema

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(TESTS_DIR)
RELEASE = os.path.join(ROOT, "data", "catalog", "releases", "catalog_0.1.1.json")
SCHEMA = os.path.join(ROOT, "governance", "schemas",
                      "catalog_release.schema.json")
LIBRARY = os.path.join(ROOT, "data", "candidates", "library.json")
REPORT = os.path.join(ROOT, "data", "spike6", "results",
                      "spike06_report.json")
IN_DIR = os.path.join(ROOT, "data", "spike6", "input", "cod")


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _load(path):
    with open(path) as fh:
        return json.load(fh)


class TestCatalogRelease(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.schema = _load(SCHEMA)
        cls.release = _load(RELEASE)
        cls.library = _load(LIBRARY)
        cls.report = _load(REPORT)

    # --- schema + gate ----------------------------------------------------

    def test_schema_is_valid(self):
        jsonschema.Draft7Validator.check_schema(self.schema)

    def test_release_against_schema(self):
        jsonschema.validate(self.release, self.schema)

    def test_release_gate_minimums(self):
        entries = self.release["entries"]
        self.assertGreaterEqual(len(entries), 8)
        self.assertGreaterEqual(self.release["validation"]
                                ["cif_validation_rate"], 0.98)
        self.assertEqual(len(self.release["entries"])
                         + len(self.release["validation"]["rejected"]), 12)

    def test_newly_fetched_minimum(self):
        self.assertGreaterEqual(self.report["n_newly_fetched"], 5)

    def test_cif_sha256_matches_files_on_disk(self):
        for e in self.release["entries"]:
            path = os.path.join(IN_DIR, f"{e['cod_id']}.cif")
            with open(path, "rb") as fh:
                digest = sha256_text(fh.read().decode("utf-8", "replace"))
            self.assertEqual(digest, e["cif_sha256"],
                             f"CIF content changed for COD {e['cod_id']}")

    def test_pyrite_entry_validated(self):
        # G0 close-out: COD *text* search "pyrite" yields only arsenopyrite /
        # chalcopyrite / marcasite; the formula search "FeS2" gives true
        # pyrites (Pa-3), so the target must now be a validated entry
        pyrite = [e for e in self.release["entries"]
                  if e["family"] == "FeS2 (pyrite)"]
        self.assertEqual(len(pyrite), 1)
        e = pyrite[0]
        self.assertIn("Fe", e["formula"])
        self.assertIn("S", e["formula"])
        sg = e["space_group"].lower().replace(":", "")
        self.assertTrue(sg.startswith("p a -3"),
                        f"pyrite must be Pa-3, got {e['space_group']}")
        self.assertEqual(e["validation"]["status"], "pass")
        # no pyrite rejection should remain
        reasons = " | ".join(r["reason"]
                             for r in self.release["validation"]["rejected"])
        self.assertNotIn("pyrite", reasons)
        self.assertTrue(all(r.get("cod_id") and isinstance(r["cod_id"], int)
                            and r.get("reason") for r in
                            self.release["validation"]["rejected"]),
                        "every rejected entry needs an integer cod_id + reason")

    def test_manifest_is_deterministic(self):
        # scientific content only: provenance timestamps + manifest itself
        # are excluded from the manifest hash
        clone = {k: v for k, v in self.release.items()
                 if k not in ("manifest_sha256", "built_at")}
        clone["source"] = {k: v for k, v in self.release["source"].items()
                           if k != "accessed_at"}
        recomputed = sha256_text(json.dumps(clone, sort_keys=True))
        self.assertEqual(recomputed, self.release["manifest_sha256"])

    # --- chemistry sanity of curated entries ------------------------------

    def test_entry_chemistry_matches_family(self):
        want = {
            "NaCl (halite)":         (["Cl", "Na"], 2),
            "AgBiS2 (schapbachite)": (["Ag", "S"], 1),
            "GaAsO4 (quartz homeotype)": (["As", "Ga"], 2),
            "TiO2 (rutile)":         (["Ti"], 1),
            "PbSO4 (anglesite)":     (["S", "O", "Pb"], 1),
            "SiO2 (quartz)":         (["Si", "O"], 2),
            "CaCO3 (calcite)":       (["Ca", "C"], 1),
            "Al2O3 (corundum)":      (["Al"], 1),
            "CaF2 (fluorite)":       (["Ca", "F"], 2),
            "MgO (periclase)":       (["Mg", "O"], 2),
            "FeS2 (pyrite)":         (["Fe", "S"], 1),
        }
        for e in self.release["entries"]:
            fam = e["family"]
            if fam not in want:
                continue
            elems, _ = want[fam]
            for el in elems:
                self.assertIn(el, e["formula"],
                              f"{fam}: missing {el} in {e['formula']}")
            self.assertTrue(e["space_group"] and e["space_group"] != "?")

    # --- clustering + validation block ------------------------------------

    def test_clustering_grouping(self):
        fam = self.release["clustering"]["families"]
        self.assertIn("PbSO4 (anglesite)", fam)
        self.assertIn("SiO2 (quartz)", fam)
        self.assertTrue(all(len(v) >= 1 for v in fam.values()))

    def test_cif_validation_rate_matches_count(self):
        n = len(self.release["entries"])
        total = n + len(self.release["validation"]["rejected"])
        self.assertAlmostEqual(self.release["validation"]["cif_validation_rate"],
                               round(n / total, 3))


class TestRuntimeLibrary(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.schema = _load(SCHEMA)
        cls.release = _load(RELEASE)
        cls.library = _load(LIBRARY)

    def test_references_release(self):
        self.assertEqual(self.library["release_ref"]["version"], "0.1.1")
        self.assertEqual(
            self.library["release_ref"]["manifest_sha256"],
            self.release["manifest_sha256"])

    def test_catalog_backed_and_legacy_kept(self):
        mats = self.library["materials"]
        ids = {m["id"] for m in mats}
        release_ids = {f"cod-{e['cod_id']}" for e in self.release["entries"]}
        self.assertTrue(release_ids <= ids, "every released entry must be in "
                                            "the library")
        self.assertIn("mat-pbso4", ids)
        self.assertIn("mat-pbso4-fe", ids)
        self.assertNotIn("mat-nacl", ids)
        self.assertNotIn("mat-sio2", ids)
        self.assertGreaterEqual(len(mats), len(release_ids) + 2)

    def test_every_entry_has_fingerprint_and_provenance(self):
        for m in self.library["materials"]:
            self.assertTrue(m["fingerprint"], f"{m['id']} lost fingerprint")
            self.assertIn("catalog_ref", m)

    def test_library_manifest_is_deterministic(self):
        mats = sorted(self.library["materials"], key=lambda e: e["id"])
        recomputed = sha256_text(json.dumps(mats, sort_keys=True))
        self.assertEqual(recomputed, self.library["manifest_sha256"])

    def test_legacy_pbso4_matches_catalog_anglesite(self):
        release_by_id = {e["cod_id"]: e for e in self.release["entries"]}
        mats = {m["id"]: m for m in self.library["materials"]}
        for leg in ("mat-pbso4", "mat-pbso4-fe"):
            ref = mats[leg]["catalog_ref"]
            self.assertEqual(ref["release"], "0.1.1")
            entry = release_by_id.get(ref["cod_id"])
            self.assertIsNotNone(entry)
            self.assertEqual(entry["family"], "PbSO4 (anglesite)")


if __name__ == "__main__":
    unittest.main()