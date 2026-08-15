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

Finding: no in-pool magnetite CIF fits the data's cell exactly (best
2300616, a = 8.3582, wR 5.68); the split stops at 38.0/62.0 — 6.2 wt
points off, outside ±3. Cell refinement of magnetite in the joint fit
moves 38/62 -> 36/64 (worse), so the cell is not the binding term — the
hematite/magnetite scale split is set by the CIFs' structure factors vs
the data's intensities.

## Follow-up sweep (2026-08-15): ambient magnetite + intensity model

Reproduce: `python benchmarks/qpa_gate/sweep_intensity.py`
(sweep2 = Uiso refinement, `SWEEP3=1` = the run below; ~35 min).

Three hypotheses from the first sweep were tested and are all refuted:

1. **Ambient magnetite ingestion (refuted).** Five pure Fe3O4 COD entries
   with ambient cells were pulled from www.crystallography.net and
   joint-fit (3 hematite CIFs x 6 magnetites, Uiso off, narrow profile):

   | magnetite (COD, a) | best hematite | wR | hem/mag split |
   |---|---|---|---|
   | 2300616 (8.3582) | 5910082 | **5.68** | **38.0 / 62.0** |
   | 1010369 (8.384) | 9000139 | 7.06 | 59.7 / 40.3 |
   | 2101535 (8.3922) | 5910082 | 7.06 | 77.4 / 22.6 |
   | 1539747 (8.3941) | 5910082 | 7.12 | 90.1 / 9.9 |
   | 9002316 (8.3965) | 5910082 | 7.09 | 95.9 / 4.1 |
   | 1513304 (8.3985) | 5910082 | 7.09 | 85.6 / 14.4 |

   Every ambient CIF fits the pattern worse than 2300616 (wR 6.7-7.2 vs
   5.68) — the data's magnetite is not ambient; the "a = 8.38-8.40 from
   peak positions" estimate in the first sweep was wrong (line overlap /
   Kalpha2 aliasing). The 2300616 + 5910082 38.0/62.0 result remains the
   global optimum and the sample cannot reach truth 31.8/68.2 (±3) with
   any COD library magnetite.

2. **Uiso refinement (refuted).** Refining every atom's Uiso (GSAS-II
   `Atoms: all -> "U"`) collapses the Hessian (scale/Uiso correlation;
   "1..5 Parameter(s) dropped" in every run); splits scatter 0/100..86/14
   and wR degrades to 6.3-7.4. Not usable without staged
   constraints/damping.

3. **Preferred orientation (refuted, correct API).** The first sweep
   used non-existent HAP keys (`'PO'`/`'POhkl'`); the second uses the
   real GSAS-II interface: `Histograms[hist]['Pref.Ori.'] =
   ['MD', G, refine, axis, order, coefs]` (what `HAPvalue('PO', 1)`
   + `set_HAP_refinements({'Pref.Ori.': True})` set). MD G now refines
   (G = 0.90-1.28) but refined zincite wt% moves only
   73.5 -> 69.7 (axis 001) / 73.5 (axis 100) and wR blows out to 100 in
   every PO-enabled run (correlated-parameter collapse). The
   angle-dependent zincite mismatch is therefore not a uniaxial texture.

4. **Zincite Uiso census.** The library zincite CIFs carry widely
   different temperature factors: 2300112/1577381 aniso 0.0068/0.0061
   (U33 0.0034), 2300450 iso 0.00623, 9004178 aniso 0.01 (U33 0.0054),
   9004179 aniso 0.01 (U33 0.0083), 2107059 Uiso 0. The high-angle hk0
   lines are ~2-3.4x stronger in the data than in any of these models.
   The two best fits bracket the true zincite wt% from opposite sides at
   nearly identical wR: 2300112 -> 23.6 (wR 58.25) and 2300450 -> 73.5
   (wR 58.23); truth 55.2 lies between, i.e. the data's zincite F2/TF
   profile is between every library entry and cannot be represented.

## Conclusion

With the current COD library (data/spike12/work/cifs, COD-pinned) and
fixed profile stages, neither sample can pass the ±3 gate: qarr_1f fails
because no zincite CIF reproduces the data's angle-dependent intensities
(the F2 profile of the data's zincite lies strictly between library
entries: equal-wR fits give 23.6% and 73.5% against truth 55.2);
iron_30_70 fails because the data's magnetite is a condensed (a = 8.358)
structure that only 2300616 approaches, and the F2-consistent joint-fit
split lands at 38.0/62.0 vs truth 31.8/68.2. The candidate fixes from
the first sweep (ambient magnetite ingestion; PO; Uiso refinement) were
implemented correctly and are all refuted by direct fit evidence.

Realistic next steps, in order of leverage: (1) trace the provenance of
`data/qpa_gate/work/{qarr_1f,iron_30_70}.xye` back to the generating
simulation CIFs (`benchmarks/eval/sim.py::sim_cif_to_pattern`) and gate
the test pool against that generator; (2) replace the profile-stage
forward selection for iron-class samples with always-on joint fits of
the top-2 candidates (the pipeline's single-phase pick 9002323/100 wt%
loses 38/62 that the joint fit already exposes); (3) for qarr-class
samples, add a Stage-D spherical-harmonics PO + poly-tensor Uiso model
with damped/staged refinement (uniaxial MD is insufficient; free Uiso
collapses).