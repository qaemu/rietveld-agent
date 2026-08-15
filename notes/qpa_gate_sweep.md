# QPA gate: CIF-combination sweep for the two failing samples

Date: 2026-08-15 — continuation of the gate work (iron_30_70, qarr_1f).

Gate contract (`benchmarks/qpa_gate/qpa_gate.py::gate`): every truth phase
>= 5 wt% must be identified and |refined wt% - truth wt%| <= 3.

Method: GSAS-II joint refinements (Hill-Howard) on the raw patterns with a
narrow data-scale profile (Caglioti U,V,W = 0, 0, 0.005; X = Y = 0; sample
shift 0), chebyshev background (5 coeffs), scales + bg refined, 60 cycles.
Molar masses for wt%: ZnO 81.38, CaF2 78.07, Al2O3 101.96,
Fe2O3 159.69, Fe3O4 231.54.

Reproduce: `python benchmarks/qpa_gate/sweep_cifs.py` (~40 min, writes
project gpx/lst under `data/qpa_gate/sweep/`).

## qarr_1f (truth: zincite 55.2 / fluorite 17.7 / corundum 27.1)

Fluorite fixed to 1000043, corundum to 1000059; zincite swept. Profile as
above.

| zincite CIF | wR | zincite | fluorite | corundum | Δmax |
|---|---|---|---|---|---|
| 2300112 (pipeline pick) | 58.25 | 23.6 | 35.9 | 40.5 | 31.6 |
| 2300450 | 58.23 | 73.5 | 12.5 | 14.0 | 18.3 |
| 9004178 | 58.40 | 74.6 | 12.0 | 13.5 | 19.4 |
| 9004179 | 65.62 | 68.0 | 15.0 | 17.0 | 12.8 |
| 1577381 | 62.02 | 20.6 | 37.3 | 42.1 | 34.6 |
| 2107059 | 68.06 | 14.2 | 40.5 | 45.2 | 41.0 |
| 2300450 + March-Dollase PO (axis 001) | 58.23 | 73.5 | 12.5 | 14.0 | 18.3 |

Finding: the zincite CIF choice swings the refined zincite wt% from 14 to
75 (truth 55.2), i.e. no library zincite reproduces the data's relative
intensities across the full angular range; preferred orientation (MD, axis
001) does not change the fit. The data's zincite high-angle lines
(110/103/112/200) run 2-3.4x stronger than the model's CIF-based
intensities (observed in the joint-fit residual: ycalc ~0.5x yobs at
56.6/62.9/68.0°), a monotonic angle-dependent mismatch that resists
profile/Uiso/seeds/PO. No combination passes (±3).

## iron_30_70 (truth: hematite 31.8 / magnetite 68.2)

Joint fits over the local library magnetite (a = 8.25-8.36) x hematite
CIFs:

| magnetite | hematite | wR | hematite | magnetite | Δmax |
|---|---|---|---|---|---|
| 2300616 (a=8.3582) | 5910082 | 5.68 | 38.0 | 62.0 | 6.2 |
| 2300616 | 9000139 | 5.98 | 0.0* | 100.0* | 31.8 |
| 2300616 | 9015065 | 5.36 | 45.3 | 54.7 | 13.5 |
| 1011032 (a=8.32) | 5910082 | 6.56 | 45.4 | 54.6 | 13.6 |
| 1011032 | 9000139 | 6.59 | 39.0 | 61.0 | 7.2 |
| 1011032 | 9015065 | 6.33 | 55.7 | 44.3 | 23.9 |
| 9002320 (a=8.3122) | 5910082 | 6.90 | 85.5 | 14.5 | 53.7 |
| 9002320 | 9000139 | 6.86 | 82.4 | 17.6 | 50.6 |
| 9002320 | 9015065 | 6.50 | 90.5 | 9.5 | 58.7 |

(* = hematite scale pinned to zero by the LSQ; parameter dropped warning.)

Finding: the data's magnetite is ambient (a = 8.38-8.40 from 220/400/440/
511 positions) but every library Fe3O4 CIF is a high-pressure variant
(a = 8.25-8.36); the best achievable joint fit (2300616 + 5910082,
wR 5.68) stops at 38.0/62.0 — 6.2 wt points off, outside ±3. Cell
refinement of magnetite in the joint fit moves 38/62 -> 36/64 (worse), so
the cell mismatch is not the binding term — the hematite/magnetite scale
split is set by the CIFs' structure factors vs the data's intensities.

## Conclusion

With the current COD library (data/spike12/work/cifs, COD-pinned) and
fixed profile stages, neither sample can pass the ±3 gate: qarr_1f fails
because no zincite CIF reproduces the data's angle-dependent intensities;
iron_30_70 fails because no ambient-magnetite CIF exists in the library
and the F2-consistent split lands 6+ points off. Fix options for future
work: ingest an ambient magnetite CIF (a ~ 8.39) from the COD web tree,
and model zincite with a specimen-consistent intensity model (PO + Uiso
refinement in Stage D, or a measured zincite CIF).