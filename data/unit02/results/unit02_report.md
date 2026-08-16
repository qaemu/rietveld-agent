# Unit 02: XRDML ingest + fingerprinting

## Fixtures (GSAS-II v5.6.3 simulated, noiseless ycalc, rounded to 0.01)

| label | anode | material | top peak 2th (deg) | sha256 |
|---|---|---|---|---|
| cu_PbSO4 | Cu | PbSO4 | 29.68 | 9cb9084c7b76 |
| fe_PbSO4 | Fe | PbSO4 | 37.56 | ba81b8e1dcd9 |
| cu_quartz | Cu | SiO2 | 25.84 | d3906c9ede93 |
| base64_PbSO4 | - | - | - | f6da72cde811 |
| nonamespace_quartz | - | - | - | bd5b312168d0 |

## Verification
- `--verify` rerun: 5 fixtures unchanged, 0 mismatches

## Physics checks (fingerprinting)
- PbSO4 Cu: top peak at 29.68 deg (unit01 measured 29.68)
- PbSO4 Fe: top peak at 37.56 deg (kinematic expectation ~37.6 deg for the same d = 3.005 A: d-space invariant)
- SiO2 (COD 1009000, quartz-family): top peak at 25.84 deg, d = 3.44 A (distinct from PbSO4 in d-space)

## Findings (agent design)
- Greedy peak-by-peak d matching over-counts across different materials (22/27 coincidental overlaps at tol=0.02 A); the d-space profile cosine is the reliable discriminator (0.908 same material cross-anode vs 0.018 different material). Peak line lists remain for reporting only.
- Instrument matching is exact (fingerprint id): pairs observations with the PRM in the calibration registry keyed on anode+wavelengths+grid.
- Fixtures embed the noiseless GSAS-II ycalc (rounded to 0.01 counts): deterministic goldens; lab data will carry Poisson noise which the profile metric tolerates.
- XRDML parser handles namespaced/plain XML and ASCII/base64 payloads; validation must reject empty and non-numeric payloads explicitly (numpy parses 'nan' as a float).

## Verdict
- [ ] goldens regenerate deterministically (hash-stable)
- [ ] xrdml parser round-trips all variants (ns / no-ns / base64)
- [ ] instrument fingerprint discriminates anodes; sample fingerprint matches across anodes in d-space and discriminates materials
