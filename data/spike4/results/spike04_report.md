# Spike 04: M1 vertical slice -- analyze

- library: data/candidates/library.json (manifest 340127a8343ea8f0...)
- registry: data/spike3/results/registry.json

## Cases

| fixture | expected | got | top sim | bundle |
|---|---|---|---|---|
| cu_PbSO4.xrdml | supported | supported | 1.0 | data/spike4/results/bundles/cu_PbSO4.xrdml.bundle.json |
| cu_quartz.xrdml | supported | supported | 1.0 | data/spike4/results/bundles/cu_quartz.xrdml.bundle.json |
| fe_PbSO4.xrdml | supported | supported | 1.0 | data/spike4/results/bundles/fe_PbSO4.xrdml.bundle.json |
| unknown_MoKa(no file, params only) | unknown_calibration | unknown | - | - |

## Findings
- Fe measurement resolves against a Cu-only fingerprint library in d-space (mat-pbso4 cross-anode), so calibration (not library) carries the instrument coupling.
- Abstention is total-precedence: calibration unknown/ambiguous stops before any hypothesis work (policy: released-only).
- M1 evidence is fingerprint-only; thresholds (top sim >= 0.35, margin >= 0.10) are explicit constants in core/verdict pending policy promotion when refinement-backed evidence lands.
- SUPERSEDED by spike 05: the 0.108 cross-anode margin shown in M1 was an exact-match artifact (query == its own library entry). The margin is now family-aware (against the best DIFFERENT phase family); see data/spike5/results/spike05_report.md. M1 verdicts unchanged, same-anode and cross-anode both hold.
- Catalogue-backing (spike 06): the M1 sampling of quartz-family patterns was chemistry-validated as GaAsO4 (quartz homeotype) and the 'NaCl' spot as Ag0.5Bi0.5S; the library and these expectations now reflect the validated identities (data/catalog/releases/catalog_0.1.0.json).
- Every bundle validates against run_bundle.schema.json at write time and in the test suite.

## Verdict
- [ ] cu_PbSO4 -> supported PbSO4
- [ ] cu_quartz -> supported SiO2
- [ ] fe_PbSO4 -> supported PbSO4 (cross-anode)
- [ ] unknown MoKa -> hard-stop abstain (calibration unknown)
- [ ] all bundles schema-valid
