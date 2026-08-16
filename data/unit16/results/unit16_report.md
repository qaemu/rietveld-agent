# Unit 16: validation of the unit-15 RQPA protocol

harness: benchmarks/protocols/validate.py | protocol: docs/rqpa_protocol.md

## A. Reproducibility

- recorded content md5: f43d8be2c932420676d612242dd049a5
- self-consistent: True
- rerun: done
- prior recorded md5 matches legacy (timing-included) hash: False (pre-canonicalization reports only)
- prior canonical md5 (committed run, timing stripped): f43d8be2c932420676d612242dd049a5
- fresh md5 (recorded): f43d8be2c932420676d612242dd049a5
- cross-run canonical identical: True
- rerun rc: 0 (900.5 s)
- verdict: REPRODUCIBLE

## B. Synthetic ground-truth recovery (known-answer test)

injected = published wt% (normalized to 100); Poisson peak-max 200000 counts (seed 16); bands: +-1.5 wt% majors / +-1.0 wt% minors

### Clinker_Nist_CuKalpha1_R1.xrdml  (pass)

wR = 0.45%, rwp_norm = 0.0045, converged = True, materialize_bkg = True
| phase | injected | recovered | |diff| | band | pass |
|---|---|---|---|---|---|
| alite-M3 | 66.0 | 66.03 | 0.030000000000001137 | 1.5 | PASS |
| aluminate-cub | 0.7 | 0.63 | 0.06999999999999995 | 1.0 | PASS |
| aluminate-ort | 1.2 | 1.19 | 0.010000000000000009 | 1.0 | PASS |
| aphthitalite | 0.8 | 0.79 | 0.010000000000000009 | 1.0 | PASS |
| belite-alphaH | 2.7 | 2.7 | 0.0 | 1.0 | PASS |
| belite-beta | 13.5 | 13.56 | 0.0600000000000005 | 1.5 | PASS |
| ferrite-C4AF | 11.1 | 11.1 | 0.0 | 1.5 | PASS |
| periclase | 4.0 | 4.0 | 0.0 | 1.0 | PASS |
- KAT score: 8/8 (87.4 s)

### Silicate_enriched_residue_Nist_CuKalpha1_R1.xrdml  (pass)

wR = 0.56%, rwp_norm = 0.0056, converged = True, materialize_bkg = True
| phase | injected | recovered | |diff| | band | pass |
|---|---|---|---|---|---|
| alite-M3 | 78.7 | 78.8 | 0.09999999999999432 | 1.5 | PASS |
| belite-alphaH | 2.9 | 2.88 | 0.020000000000000018 | 1.0 | PASS |
| belite-beta | 13.4 | 13.34 | 0.0600000000000005 | 1.5 | PASS |
| periclase | 5.0 | 4.97 | 0.03000000000000025 | 1.5 | PASS |
- KAT score: 4/4 (55.5 s)

### aluminate_enriched_residue_clinkerNIST_180718_R1.xrdml  (pass)

wR = 0.38%, rwp_norm = 0.00377, converged = True, materialize_bkg = True
| phase | injected | recovered | |diff| | band | pass |
|---|---|---|---|---|---|
| aluminate-cub | 5.1657 | 5.22 | 0.05429999999999957 | 1.5 | PASS |
| aluminate-ort | 7.6023 | 7.87 | 0.2677000000000005 | 1.5 | PASS |
| aphthitalite | 2.4366 | 2.33 | 0.1065999999999998 | 1.0 | PASS |
| ferrite-C4AF | 68.0312 | 67.8 | 0.23120000000000118 | 1.5 | PASS |
| periclase | 16.7641 | 16.78 | 0.015900000000002024 | 1.5 | PASS |
- KAT score: 5/5 (40.8 s)

### Clinker_Synchrotron.dat  (pass)

wR = 0.46%, rwp_norm = 0.00461, converged = True, materialize_bkg = True
| phase | injected | recovered | |diff| | band | pass |
|---|---|---|---|---|---|
| alite-M3 | 65.3935 | 65.65 | 0.2565000000000026 | 1.5 | PASS |
| aluminate | 1.9898 | 1.8 | 0.18979999999999997 | 1.0 | PASS |
| aphthitalite | 0.5699 | 0.56 | 0.009899999999999909 | 1.0 | PASS |
| belite-alphaH | 2.9997 | 3.01 | 0.010299999999999976 | 1.0 | PASS |
| belite-beta | 13.7986 | 13.76 | 0.038600000000000634 | 1.5 | PASS |
| ferrite-C4AF | 11.5988 | 11.58 | 0.018800000000000594 | 1.5 | PASS |
| periclase | 3.6496 | 3.64 | 0.009599999999999831 | 1.0 | PASS |
- KAT score: 7/7 (322.0 s)

## C. Gate scoring (real data, fresh report)

| sample | model | wR | wR gate | rwp | rwp gate | wt gate | rank tau | operational | publication |
|---|---|---|---|---|---|---|---|---|---|
| Clinker_Nist_CuKalpha1_R1.xrdml | M3 | 15.07 | FAIL (6.5) | 0.15071 | PASS (1.0) | FAIL | 0.3571 | yes | FAIL |
| Clinker_Nist_CuKalpha1_R1.xrdml | T1 | 20.91 | FAIL (6.5) | 0.20909 | PASS (1.0) | FAIL | 0.2857 | yes | FAIL |
| Silicate_enriched_residue_Nist_CuKalpha1_R1.xrdml | M3 | 20.18 | FAIL (6.5) | 0.20175 | PASS (1.0) | FAIL | 0.6667 | yes | FAIL |
| Silicate_enriched_residue_Nist_CuKalpha1_R1.xrdml | T1 | 27.11 | FAIL (6.5) | 0.27113 | PASS (1.0) | FAIL | 0.6667 | no | FAIL |
| aluminate_enriched_residue_clinkerNIST_180718_R1.xrdml | M3 | 13.56 | FAIL (6.5) | 0.13564 | PASS (1.0) | FAIL | -0.2 | yes | FAIL |
| Clinker_Synchrotron.dat | M3 | 9.76 | FAIL (5.0) | 0.26437 | PASS (0.5) | FAIL | 0.619 | yes | FAIL |

- operational: 5/6 | publication: 0/6

Constrained trace phases (fixed-composition reinsertion at the normalized published share):
- aluminate_enriched_residue_clinkerNIST_180718_R1.xrdml: aphthitalite

## Summary

- reproducibility: identical content hash across independent runs
- KAT recovery: 4/4 samples within band
- publication gate: 0/6 (wR 9.76-27.11 vs <=6.5/<=5) -- protocol operational, fits not yet publication-grade
