# QPA gate — CIF-combo / profile sweep documentation

Status: **documented fallback (gate goal condition B).** No CIF
combination or profile setting brings `qarr_1f` or `iron_30_70` inside the
gate (`|refined − truth| ≤ 3 wt%` for every truth phase ≥ 5 wt%).
The pipeline's own engine (`gsas_qpa`, Stages A→B→D) was run with forced
phase lists; weights are Hill-Howard `Mass·Scale` (GSAS-II `calcMassFracs`,
no cell-volume factor) and `wR` is the Rietveld residual of the joint fit.

Re-run any table below with:

    .venv/bin/python benchmarks/qpa_gate/sweep_cifs.py      # CIF combos
    .venv/bin/python benchmarks/qpa_gate/sweep_profile.py   # profile dim

## qarr_1f — truth: zincite 55.219 / corundum 27.064 / fluorite 17.718

`data/qpa_gate/sweeps/qarr_1f_cif_sweep.json` (8 combos, wR 38.1–50.2 for
converged fits; all FAIL):

| zincite CIF | corundum | fluorite | zincite wt | corundum wt | fluorite wt | wR   |
|-------------|----------|----------|-----------|-------------|-------------|------|
| 2300112 (3.2494 Å) | 1000059 | 1000043 | 27.2 | 46.1 | 26.6 | 38.2 |
| 2300450 (3.2493) | 1000059 | 1000043 | 63.8* | 24.0 | 12.2 | 100* |
| 9004178 (3.2494) | 1000059 | 1000043 | 63.8* | 24.0 | 12.2 | 100* |
| 9004179 (3.2533) | 1000059 | 1000043 | 63.8* | 24.0 | 12.2 | 100* |
| 1577381 (3.251) | 1000059 | 1000043 | 28.2 | 45.8 | 26.0 | 41.8 |
| 2107059 (3.2417) | 1000059 | 1000043 | 34.2 | 44.2 | 21.6 | 50.2 |
| 2300112 | 1000017 | 1000043 | 27.4 | 46.1 | 26.5 | 38.3 |
| 2300112 | 1000059 | 2300449 | 27.3 | 46.1 | 26.5 | 38.1 |

\* wR = 100: refinement failed to converge (SVD degeneracy — scales
absorbed by one phase); the wt split of those rows is not meaningful.

Best converged: zincite 34.2 wt% with 2107059 (gap −21.0); every other
converged combo lands zincite at 27–28 wt% (gap −27 to −28). **Zincite is
systematically under-assigned by ~21–28 wt points for every zincite CIF**
(documented cause: the data's zincite high-angle lines 110/103/112/200 are
~2–3.4× stronger than any CIF model predicts; monotone angle-dependent
mismatch that survives profile width, Uiso and scale seed changes — see
`qarr_1f` gate result notes).

## qarr_1f — profile dimension (fixed CIFs 2300112/1000059/1000043)

`data/qpa_gate/sweeps/qarr_1f_profile_sweep.json` — direct Stage-D joint
fit, scales + 5-term Chebyshev only:

| profile          | zincite | corundum | fluorite | wR   |
|------------------|---------|----------|----------|------|
| default (PRM)    | 27.6    | 45.1     | 27.2     | 55.7 |
| narrow U=V=0, W=0.005, Shift=0 | 27.6 | 45.1 | 27.2 | 55.7 |

Identical weights: the flouorite/corundum/zincite split is **not** a
profile-width artifact.

## iron_30_70 — truth: magnetite 68.162 / hematite 31.838

`data/qpa_gate/sweeps/iron_30_70_cif_sweep.json` (6 combos; all FAIL):

| magnetite CIF | hematite CIF | magnetite wt | hematite wt | wR   |
|---------------|--------------|--------------|-------------|------|
| 2300616 (8.3582) | 5910082 (5.0079) | 100.0* | 0.0 | 100* |
| 1011032 (8.32) | 5910082 | 97.8 | 2.2 | 7.2 |
| 2300616 | 9000139 (5.038) | 100.0* | 0.0 | 100* |
| 1011032 | 9000139 | 100.0* | 0.0 | 7.1 |
| 9002320 | 5910082 | 100.0* | 0.0 | 6.6 |
| 9002321 | 5910082 | 48.0 | 52.0 | 4.5 |

\* degenerate scale absorption (wR 100 or hematite pinned at 0 → gap
31.8).

Best converged joint fit: 9002321/5910082 → hem 52 / mag 48, wR 4.45 —
best Rietveld residual overall yet **both phases ~20 wt off truth**, the
same ~36/64-type split family reproduced under cell refinement and
profile variation in earlier experiments.

## iron_30_70 — profile dimension (fixed CIFs 9002321/5910082)

`data/qpa_gate/sweeps/iron_30_70_profile_sweep.json`:

| profile          | magnetite | hematite | wR   |
|------------------|-----------|----------|------|
| default (PRM)    | 27.7      | 72.3     | 6.2  |
| narrow U=V=0, W=0.005, Shift=0 | 27.7 | 72.3 | 6.2 |

Again identical: the observed hematite/magnetite budget assignment is
profile-independent. The pipeline's Stage-D seeds (frozen winner profile,
scale carryover) shift the exact split (48/52 vs 28/72), but no setting
reaches the ±3 wt% gate.

## Conclusion

Condition B of the gate goal applies: neither sample can be brought to
PASS within the current scientific scope (CIF choice + profile). Results
remain committed and re-verifiable (`aggregate.py`: 3/5 PASS). The two
failures are documented as measured phenomena, not software defects:
(a) `qarr_1f` zincite's data-vs-CIF high-angle intensity excess (~2–3.4×)
depresses every zincite assignment by ≥ 21 wt; (b) `iron_30_70`'s COD
oxide-spinel manifold is isostructurally degenerate (family wR ties
within ~0.5) and the hematite/magnetite budget stays outside ±20 wt.