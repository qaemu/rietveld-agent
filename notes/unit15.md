# Unit 15: publication-grade RQPA runner -- run log and locked decisions

Runner: `benchmarks/protocols/rqpa_protocol.py`
Result artifacts: `data/unit15/results/unit15_report.{json,md}` (md5
recorded in the JSON header).
Protocol doc: `docs/rqpa_protocol.md` (section 7 is the normative summary
of this page).

## What the unit delivers

Deterministic, staged-ladder RQPA of the four SRM 2686a patterns with the
unit-14 published polymorph set, dual alite model (M3 / T1), Hill-Howard
normalization, per-sample tiering, and published-inventory comparison.

## Empirical probes that drove the final configuration

All probes ran against the real patterns with GSAS-II (vendored,
scriptable), work files in `data/unit15/work/` (`.gpx` ignored by git,
reproducible from the runner).

### 1. Aphthitalite on the aluminate residue (the long one)

- `alu_n1/alu_n2` (first probe): with aphthitalite (9007639) in the phase
  set the refinement aborts at cycle 0 (`wR=100 cycles=0`); without it,
  wR=14.87%, converged, sane scales. → aphthitalite is the killer on the
  alu window.
- The "Symmetry element ... not matched in GSAS-II setting" console
  warning is a red herring from other phases (also appears on runs that
  converge); the imported SGData for 9007639 is verified correct
  (P-3m1, trigonal, 3m1 Laue, 12 ops with inversion, 6 atoms).
- Scale probe matrix (alu, all with published-prior scale init unless
  noted): `alu_noinit` (GSAS defaults, aphth scale refined) → wR=100,
  cycles=0; `alu_aphfixed` (aphth scale pinned at the prior) → wR=15.49
  but only in the 4-phase set without aluminate-ort; `alu_tiny`
  (aphth scale from 0.01, refined) → aborts, GSAS-II drives it to 0;
  `alu_rel` (release the pinned scale after a converged pass) → cycle-0
  abort; `alu_ort4` (ferrite, periclase, cub, ort -- no aphth) → wR=14.87,
  scales 6.0/81/0.56/3.4 (healthy); `alu_all5` (same + pinned aphth) →
  wR=16.71 but the other scales explode to ~1e13 and wt% garbage.
- Conclusion: on the alu data window the aphthitalite scale column is
  numerically degenerate (SVD-zeroed from any start) AND its mere
  presence corrupts the Hessian even pinned. Deterministic fix: the
  runner's `DROP_FALLBACK` retries the ladder without aphthitalite and
  records `phases_dropped` in the report row; the residue report carries
  the renormalization caveat. Clinker/silicate/sync (where published
  aphthitalite is 0.05-0.43 wt%) refine it normally.

### 2. Published-prior scale init (`_init_scales` / `INIT_SCALES`)

Helps only where GSAS-II's default 1.0/1e-12 trace starts make the first
Marquardt step diverge (alu residue). It *broke* the clinker and
silicate Cu ladders at the "+belite cells" stage (regression reproduced
twice), so it now applies to the alu sample only.

### 3. Staged ladder vs single pass

Single-pass budget: belite cell metric blew up on the Cu patterns.
Ladder (scales → +alite cell → +belite cells → +minor/periclase cells)
restarts the LM from the previous converged solution and is stable.

### 4. Alu cell refinement removed

Refining the ferrite/orthorhombic-aluminate cells from the
trace-dominated residue window finds a spurious supercell (volumes x10,
scales x1e12). Cells stay at the unit-14 published values.

### 5. Sync pattern: T1 variant excluded

alite-T1 + belite-β scale columns are jointly degenerate on the 0.825 Å
window (pin alite → GSAS-II drops belite-β to zero; unpinned → cycle-0
abort). The published sync fit is M3-only; sync is reported on the M3
model (wR = 9.76%).

### 6. Variant-substitution bug (fixed mid-session)

Attempt-loop reset `PHASESETS[fname] = saved` used the pre-substitution
list, so T1 runs silently used the M3 CIF and died at stage 2 with
KeyError `'alite-T1'`. Now `base`/`variant` are separated and
`PHASESETS[fname]` is only reset to `base` after each attempt.

## Final table

| sample [model] | wR % | rwp_norm | tier | worst |phase|diff| |
|---|---|---|---|---|---|
| clinker Cu [M3] | 15.07 | 0.1507 | ok | 29.3 (alite) |
| clinker Cu [T1] | 20.91 | 0.2091 | ok | 10.5 (C3A) |
| silicate [M3] | 20.18 | 0.2018 | ok | 19.8 (alite) |
| silicate [T1] | 27.11 | 0.2711 | wR-over | 5.6 (alite) |
| aluminate residue | 14.87 | 0.1487 | ok | 61.4 (ferrite) |
| clinker sync [M3] | 9.76 | 0.2644 | ok | 32.1 (alite) |

The large diffs concentrate on strongly absorbing phases (alite M3,
ferrite) — attributed to microabsorption, addressed by Brindley
corrections in unit 17; unit 16 validates reproducibility/rank
structure with a synthetic-pattern harness.

Open items for units 16-17: gate wR <= 6.5 (Cu) / <= 5 (sync) vs
the published full-budget fits; ferrite/alite bias after Brindley;
bedding/porosity on the residue patterns.