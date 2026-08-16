# Gate G0 — Phase-0 completion checklist

Status: **all items met** (2026-08-09). Re-verify anytime with the commands
at the bottom; every claim below has a live report artifact behind it.

Phase-0 is the engineering foundation: GSAS-II pinning, ingestion,
calibration registry, instrument-aware simulation, threshold evidence,
and a chemistry-validated COD catalog. Gate G0 does **not** claim any
scientific conclusions — it certifies that the deterministic core,
controlled inputs, and evidence trail are trustworthy enough to start
verification work (unit 07+, Rietveld-backed).

## Criteria

### C1 — Catalog release (unit 06)

- [x] `data/catalog/releases/catalog_0.1.0.json` validates against
      `governance/schemas/catalog_release.schema.json`
      (schema versioned, `additionalProperties: false`).
- [x] >= 8 entries: **12** (halite, schapbachite, GaAsO4 quartz homeotype,
      rutile, periclase, magnetite-family, anglesite, calcite, corundum,
      pyrite, fluorite, quartz).
- [x] Every entry carries `cod_id`, `cif_sha256` (verified against the CIF
      on disk in tests), reduced `formula`, `mineral_name`, `space_group`,
      `cell`, `validation.status`, and `source` with COD URL + citation.
- [x] >= 5 entries newly fetched from COD for this project: **9**
      (unit-06-introduced structures; 3 pre-existing fixture CODs excluded
      by definition).
- [x] `cif_validation_rate >= 0.98`: **1.0** (12/12 after the pyrite
      close-out; no rejections outstanding).
- [x] Rejected list is honest: recording keeps per-candidate reasons; the
      pyrite rejection observed in the first runs was a COD *text-search*
      artifact (arsenopyrite/chalcopyrite/marcasite) — resolved by the
      formula search "FeS2" (true pyrite Pa-3, COD 1544891).
- [x] Deterministic manifests: `manifest_sha256` covers scientific content
      only (timestamps excluded) — two consecutive builds give identical
      hashes; the test suite recomputes the manifest from the payload and
      from the library materials and matches both.
- [x] Clustering/family grouping recorded; chemistry sanity per family
      asserted in tests.
- [x] Licensing notice: COD is CC0/public domain, attribution via per-entry
      citations.

Evidence: `data/unit06/results/unit06_report.{md,json}`.

### C2 — Runtime library (unit 06, rebuild)

- [x] `data/candidates/library.json` rebuilt **from the release**: every
      `cod-<id>` material points at `catalog_ref {cod_id, release}` and
      embeds the release manifest (controlled-input chain: release -> library).
- [x] Legacy M1 fixture entries (`mat-pbso4`, `mat-pbso4-fe`) kept with
      provenance-corrected pointers (COD anglesite reference); superseded
      mislabels `mat-nacl` (actually Ag0.5Bi0.5S) and `mat-sio2` (actually
      GaAsO4) dropped.
- [x] `library.manifest_sha256` deterministic (recomputed in tests);
      `release_ref` matches the release manifest.
- [x] 14 materials (12 cod-backed + 2 legacy fixtures).

Evidence: `data/candidates/library.json`,
`data/unit06/results/unit06_report.md` (M1 provenance audit table).

### C3 — End-to-end pipeline + bundles (unit 04)

- [x] `cu_PbSO4` -> `supported`, primary family **PbSO4 (anglesite)**,
      sim 1.000, family margin ~0.91.
- [x] `cu_quartz` -> `supported`, primary family **GaAsO4 (quartz
      homeotype)** (chemistry-validated identity; real SiO2 quartz sits in
      the library as its own entry, homeotype vs quartz separated).
- [x] `fe_PbSO4` -> `supported` (cross-anode Fe vs Cu-only fingerprints in
      d-space; family margin ~0.92 while the raw twinning margin is ~0.11 —
      the family-aware margin is what carries the verdict, per unit 05).
- [x] Unknown instrument (MoKa, no calibration record) -> **hard-stop
      abstain** before any hypothesis work (policy: released-only
      conclusions from unknown calibrations).
- [x] Every bundle validates against
      `governance/schemas/run_bundle.schema.json` at write time and in the
      test suite; bundles carry `phase_family` (real family, not entry
      name), evidence separation, controlled versions, plan hash.

Evidence: `data/unit04/results/unit04_report.{md,json}` + bundles in
`data/unit04/results/bundles/`.

### C4 — Threshold evidence (unit 05) + statistics gate

- [x] Constants `min_top_similarity = 0.35`, `min_margin = 0.10` evaluated
      against 50 deterministic seeded replicates x 10 perturbation levels x
      2 anodes x 4 materials: verdict **KEEP** (family-aware margin).
- [x] First majority-support flip at perturbation level L4 for
      sio2-cu / nacl-cu; thresholds are therefore **not** tightened at G0
      (documented in `data/unit05/results/unit05_report.md`).
- [x] Counting-statistics gate `min_peak_counts` derived (100000) and
      enforced: below-gate abstains before verdict; # gate boundary never
      overrides evidence.
- [x] Evidence level remains honestly `fingerprint-only`; refinement-backed
      evidence is the declared next milestone (unit 07), not a G0 claim.

Evidence: `data/unit05/results/unit05_report.{md,json}`.

### C5 — Calibration registry (unit 03)

- [x] Registry schema-versioned, manifests hashed, review/approval
      machinery exercised in docs; Cu and Fe calibrations resolve by
      fingerprint.
- [x] `XYE-DEFAULT` is the sole admin-approved fallback for plain XYE
      input; XRDML always resolves via fingerprints.
- [x] Mo (and any unregistered) instrument -> governed hard-stop abstain
      (deliberate at G0; an Mo calibration entry is an optional post-G0
      item, lab-dependent).

Evidence: `data/unit03/results/unit03_report.{md,json}`.

### C6 — GSAS-II pinning + instrument-aware simulation (units 01/02)

- [x] GSAS-II installed via official channels, never bundled; environment
      pinned (version + python captured per bundle), scriptable API unit
      recorded.
- [x] Instrument-aware simulator (Cu/Fe anodes, Cu-Ka/Fe-Ka doublets,
      Bragg-Brentano, peak shapes, counting noise) reproduces the fixture
      patterns; seeded and deterministic (unit 05 replicates).
- [x] Fixtures: PbSO4 (Cu + Fe), quartz homeotype (Cu), NaCl (Cu) live in
      `tests/fixtures/xrdml/` and underpin C3/C4.

Evidence: `data/unit{,2}/results/unit*_report.md`.

### C7 — Test suite

- [x] `python -m pytest tests/ -q` -> **all pass** (79 tests), including
      `tests/test_unit06.py` (release schema, gate minimums, CIF sha256 on
      disk, pyrite entry Pa-3, manifest determinism, library provenance).

## Re-verify commands

```bash
python benchmarks/protocols/06_catalog_build.py   # exit 0, rate 1.0
python benchmarks/protocols/04_analyze.py         # all cases ok=True
python benchmarks/protocols/05_eval_thresholds.py # exit 0, verdict KEEP
python benchmarks/protocols/08_measured_verify.py # exit 0, report written
python benchmarks/protocols/09_multiphase.py      # exit 0, report written
python benchmarks/protocols/10_real_data.py       # exit 0, report written
python -m pytest tests/ -q                           # all pass
```

Run twice is by design: build scripts are idempotent and must reproduce
identical manifests on cache.

## Known open items (post-G0, not gate failures)

1. ~~Refinement-backed verification stage~~ **DELIVERED (unit 07/08)**:
   bounded verification is wired into `cli/analyze` (runs after a supported
   fingerprint verdict; skips on abstain/held/unsupported; `--no-verification`
   to disable). Two engines, both policy-bounded (bg + shift + scale; cell
   fixed in the measured engine — the fingerprint stage already constrained
   the family):
   - sim-observed path (`refine_candidate`, unit 07): GSAS-II refinement of
     each candidate against the catalog-CIF observation; e2e re-confirms the
     fingerprint top family (anglesite Rwp≈0.00 vs calcite 0.93; GaAsO4 ≈0.00
     vs SiO2 0.98); report in `data/unit07/results/unit07_report.{md,json}`.
   - measured path (`refine_measured_candidate` via `verify_measured`,
     unit 08): deterministic numpy bounded fit on the ACTUAL measured grid
     with counting weights; candidate models are protocol GSAS-II simulations
     (cached, deterministic). Noise envelope (L1..L5 × 3 seeds, anglesite vs
     calcite): true phase lowest-Rwp at every level; worst-case true Rwp
     0.50 (L5, 10k counts + 0.30 mm displacement) exceeds `confirm.max_rwp`
     0.35 — honest bounded-budget limitation at extreme noise/displacement,
     bound left at 0.35 (envelope not stricter). Measured e2e: cu_quartz
     (fixture == catalog CIF 1009000) confirms in-bounds (Rwp ≈ 0); cu_PbSO4
     (APS Wyckoff structure ≠ catalog 1010950) stays fingerprint-consistent +
     lowest-Rwp but out of policy bounds — strength-of-model evidence, not
     identity proof; report in
     `data/unit08/results/unit08_report.{md,json}`, tests in
     `tests/test_unit08.py`.
2. Mo-anode calibration entry (only if the lab measures Mo).
3. Multi-phase policy (2-phase cap is the Phase-0 policy field; enabled by
   the catalog but unexercised — unit 09 measures the gap empirically:
   a weak 2nd phase (~2% peak intensity) is absorbed by the single-phase
   bounded fit and passes in-bounds; an equipeak 50/50 mixture honestly
   ABSTAINS at the fingerprint stage (family margin 0.054 < 0.10, bundle
   held). Sources + reproductions in `data/unit09/` and
   `benchmarks/protocols/09_multiphase.py`).
4. `docs/` formalization of release process (review/approval of catalog
    releases, schema-validated) — draft state remains `draft` until then.
5. **Out-of-sample validation (unit 10)**: real RRUFF published patterns
   (quartz/calcite/corundum/rutile, identities fixed by REFINE v3.0 cell
   refinement) abstain on the 1e5 counting gate (archive fast scans,
   1e3-6e3 peak counts) and 2/4 fingerprint rankings match (calcite 0.92,
   corundum 0.64; quartz fails on the 3.34 A vs anglesite collision,
   rutile on weak counts). Fingerprint/verification are calibrated on
   noiseless protocol sims only — real-data discrimination remains
   unvalidated; sources + numbers in `data/unit10/` and
   `benchmarks/protocols/10_real_data.py`.