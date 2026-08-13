# Known-answer validation of an automated Rietveld QPA engine

**Reproducibility, synthetic recovery, and gate scoring** — Paper P-2 of
the docs/papers series.

**Thesis.** The scientific validity of an automated Rietveld quantitative
phase analysis (RQPA) pipeline can be established with three independent
checks: (A) **bit-level reproducibility** across independent runs, (B)
**known-answer recovery** of injected ground-truth fractions on
noise-corrupted synthetic patterns, and (C) **honest, scored acceptance
gates** on real data.

**Status.** Working paper, v0.2.0, August 2026. All numbers are from
`data/spike16/results/spike16_report.md` (reproduce with
`python3 benchmarks/spikes/spike_16_validate.py`). Constructed following
the fifteen-step writing framework of Drake & Han (2025)
[doi:10.1371/journal.pcbi.1013505](https://doi.org/10.1371/journal.pcbi.1013505).

---

## 1. Introduction

A refinement pipeline that runs unattended must earn trust differently
from an interactive analyst session. There is no human to notice a
plausible-but-wrong fit, and no session log to audit.

**Gap.** Published RQPA validation usually stops at "converges and looks
reasonable". What is missing is a three-pronged, executable proof: the
results reproduce, the pipeline recovers known truth, and the remaining
distance to publication quality is measured and attributed.

**Response.** Spike 16 of the `rietveld-agent` repository [1] implements
exactly that proof for the protocol of Paper P-1 [2]. Its harness,
`benchmarks/spikes/spike_16_validate.py`, runs three independent checks
and always completes; verdicts are data, not exit codes.

## 2. Methods

### 2.1 A — Reproducibility check

The full spike-15 suite (six refinements, ≈15 min on the development
hardware) is re-run via subprocess. The report's canonical content hash —
computed over all fields except wall-clock timing, which legitimately
varies — is compared with the fresh run. An identical hash across
independent runs is the reproducibility verdict.

### 2.2 B — Synthetic known-answer test (KAT)

For each sample (M3 model), the published wt% are **injected** as ground
truth, normalized to sum to 100 (the raw Table-3 aluminate-residue row
sums to 102.6 owing to replicate-mean rounding; the normalized reference
`published_norm()` is the shared comparison baseline of the repository).
The injection procedure:

1. per-phase scales started from the published-fraction prior via
   Hill–Howard inversion [3], S ∼ w/(M·V);
2. only the Chebyschev background refined (three cycles, scales pinned),
   so the materialized calculated pattern equals the published-fraction
   phase signal plus a smooth window-adapted background — no fit-residual
   leakage;
3. Poisson counting noise added at a 200 000-count peak maximum, seed 16
   (`benchmarks/eval/noise.py`);
4. the **unmodified** protocol ladder of Paper P-1 runs on the synthetic
   pattern from its own starts, with its own constraint/fallback logic;
5. recovery scored vs. injected truth: |Δ| ≤ 1.5 wt% for published
   majors (≥ 5 wt%), ≤ 1.0 wt% for minors.

### 2.3 C — Gate scoring on real data

Each of the six real refinements is scored against publication-style
gates: wR ≤ 6.5% (Cu) / ≤ 5% (synchrotron); rwp_norm < 1.0 / 0.5;
per-phase wt% within the same bands; rank agreement (Kendall τ).
*Operational* = converged, no error markers, wR < 25. *Publication* =
all gates pass. Every instrument, band and gate definition is recorded in
`docs/rqpa_protocol.md` and in the spike-16 report, together with the
probe classification used during development: **wheel** = procedural
choice; **counter** = counterexample probe; **pin** = pinned constant;
**gate** = guard on a comparison.

## 3. Results

### 3.1 A — Reproducibility

Two independent full-suite runs produced bit-identical payloads:

| run | canonical md5 |
|---|---|
| committed run | `f43d8be2c932420676d612242dd049a5` |
| fresh rerun | `f43d8be2c932420676d612242dd049a5` |

`cross_run_identical = True`, rerun rc = 0 (900.5 s), fresh report
self-consistent. **Verdict: REPRODUCIBLE.**

One harness bug is worth recording: an earlier version crashed in
`write_report` (KeyError `prior_canonical_md5`) when a killed subprocess
was still writing the shared work directory. The writer now renders a
failed/skipped rerun explicitly — the verdict is data and the harness
always completes.

### 3.2 B — Known-answer recovery

All four samples pass; **24/24 phase rows within band** (Table 1);
synthetic fits land at the noise floor (wR 0.38–0.56%).

**Table 1.** Synthetic known-answer recovery (peak 200 000 counts, seed
16; bands ±1.5 majors / ±1.0 minors wt%).

| sample | KAT score | wR (%) | max \|Δ\| (wt%) |
|---|---|---|---|
| clinker Cu | 8/8 | 0.45 | 0.07 (aphthitalite recovered freely) |
| silicate residue | 4/4 | 0.56 | 0.10 |
| aluminate residue | 5/5 | 0.38 | 0.27 (ferrite 67.8 vs. 68.0) |
| clinker synchrotron | 7/7 | 0.46 | 0.27 (aluminate 1.80 vs. 1.99) |

The recovery ceiling is ≈0.3 wt% — the size of the Poisson noise floor at
this counting statistics. Three fallouts matter:

- **aphthitalite** recovers freely on the synthetic aluminate window
  (2.33 vs. 2.44 injected, within band), proving the degeneracy observed
  on the real pattern is *real-data specific* (peak overlap /
  microabsorption), not intrinsic to the window;
- the **aluminate ferrite** deficit is a normalization artifact of the
  *reference*, not of the pipeline: 67.8 vs. 68.03 normalized injected,
  |Δ| = 0.23, PASS (the earlier 4/5 score of 67.82 vs. 69.80 raw was
  scored against a reference that sums to 102.6);
- the **real-data alu** fit gained from the cell + breadth ladder: wR
  14.9 → 13.56, ferrite 8.4 → 17.4 wt% (remainder of the gap is
  microabsorption — spike 17).

### 3.3 C — Gate scoring

**Table 2.** Real-data gate scores (rwp_norm gate passes everywhere; wt
gate fails; the aluminate residue reports 5/5 phases with aphthitalite as
fixed-composition constraint).

| sample | model | wR (%) | wR gate | rank τ | operational |
|---|---|---|---|---|---|
| clinker Cu | M3 | 15.07 | FAIL (6.5) | 0.36 | yes |
| clinker Cu | T1 | 20.91 | FAIL (6.5) | 0.29 | yes |
| silicate residue | M3 | 20.18 | FAIL (6.5) | 0.67 | yes |
| silicate residue | T1 | 27.11 | FAIL (6.5) | 0.67 | no |
| aluminate residue | M3 | 13.56 | FAIL (6.5) | −0.20 | yes |
| clinker synchrotron | M3 | 9.76 | FAIL (5.0) | 0.62 | yes |

**Operational: 5/6; publication: 0/6**, as expected. Every analysis
completes, every phase carries a wt%, and the protocol runs
deterministically — but the fits are not publication-grade, and the
report says so. The dominant residual is microabsorption, which is the
target of the Brindley-correction spike 17.

## 4. Discussion and closing

Three results carry the validity argument:

1. **Reproducibility is a machine-checked fact**, not a claim — the
   canonical content hash is identical across independent full-suite
   runs.
2. **The engine recovers known fractions to the noise floor** of the
   data it was not trained on — the KAT design reuses the unmodified
   production ladder, so the recovery is a property of the pipeline.
3. **The gates are honest**: the publication gap is attributed to physics
   (microabsorption) with a concrete remedy (spike 17), and the
   normalization trap — scoring against a reference that sums to 102.6% —
   was found by the harness and removed by adopting `published_norm()` as
   the shared comparison baseline.

Nothing in this paper claims accuracy the gates do not certify. Validity
here means: *reproduces, recovers known truth, and reports its remaining
distance truthfully*. As in Papers P-1 [2] and P-3 [4], the composition
follows the fifteen-step framework [5].

## References

1. qaemu (2026). rietveld-agent: a deterministic Rietveld QPA engine for scientific agent runtimes. https://github.com/qaemu/rietveld-agent-spikes (Apache-2.0)
2. qaemu (2026). *Reproducible Rietveld quantitative phase analysis under a bounded budget* — Paper P-1 of this series. [paper1_protocol.md](paper1_protocol.md)
3. Hill, R. J., & Howard, C. J. (1987). Quantitative phase analysis from neutron powder diffraction data using the Rietveld method. *J. Appl. Cryst.*, 20, 467–474. doi:10.1107/S0021889887086199
4. qaemu (2026). *Deploying rietveld-agent on OpenCode, Claude Code, and Codex* — Paper P-3 of this series. [paper3_deployment.md](paper3_deployment.md)
5. Drake, J. M., & Han, B. A. (2025). How to write a scientific paper in fifteen steps. *PLoS Comput. Biol.*, 21(9), e1013505. doi:10.1371/journal.pcbi.1013505 (PMC12459795)

---

*License: Apache-2.0. Full results: `data/spike16/results/spike16_report.{json,md}`.*