# Decision criteria (contracts/0.1.0, draft)

## Status semantics (per phase family)

| Status | Meaning |
|---|---|
| `supported` | Diffraction-based and refinement-based evidence jointly satisfy the strict criteria below. |
| `inconclusive` | Plausible but not proven; alternatives remain. Always the ceiling for sample-name-context-only candidates. |
| `not_selected` | Not chosen during bounded search. **This is not equivalent to absent.** |
| `out_of_domain` | Material/data characteristics outside declared capability. |
| `held` | Numerical or scientific validation requires expert review. |
| `failed` | Invalid input, missing calibration, unsupported format, numerical failure, unavailable dependency. |

## `supported` requires ALL of (conjunction, no exceptions)

1. **Independent retrieval**: the family was retrieved from the pattern alone
   (blind retrieval), without use of the sample name.
2. **Improvement under equal budget**: adding it improves the comparable model
   under identical refinement budgets, exceeding the policy threshold
   `min_delta_gof_improvement` (benchmark-derived, not ad hoc).
3. **Positive, stable, plausible scale/fraction**: scale/fraction positive and
   stable across starts; within `fraction_bounds`.
4. **Diagnostic reflections present** or clearly documented as overlapped
   (overlap computed, not asserted).
5. **No decisive expected reflections absent** (computed against the family's
   expected reflection set).
6. **Survives alternative-model testing**: alternatives within
   `tolerance_chi2_ratio` must also be reported (`inconclusive`), not hidden.
7. **Remaining unexplained peaks do not invalidate the phase set** (within
   `max_unexplained_peak_budget`).
8. **Numerical and covariance checks pass**: converged shifts, sane
   correlations, no physical-bound violations.

## Abstention and out-of-domain triggers (minimum set)

- Organic materials, major amorphous fractions (signal fraction above
  `max_amorphous_signal_fraction`), nanocrystalline-dominated patterns
  (estimated crystallite size below `max_nanocrystallite_size_nm`).
- Strongly disordered/mixed-layer clays; catalog-absent phases evinced by
  strong unmatched decisive peaks.
- Synchrotron, neutron, TOF, capillary, non-Bragg-Brentano data.
- Highly multiphase patterns (phase count beyond `phase_cap`).
- Geographic sample-name tokens without a validated geological reference
  source (non-compositional context).

## `held` triggers (requires expert review)

- Preferred-orientation indicators, plausible microabsorption.
- Unstable scales/fractions across refinement starts.
- Correlation-matrix alarms or covariance warnings.
- Borderline alternative models within statistical tolerance.
- Unresolved mid-strength peaks; catalog entries flagged low-confidence.

## `failed` triggers

- Unreadable/unsupported input format; missing or ambiguous calibration
  (exactly-one default violated for XYE; non-unique/unknown fingerprint match
  for XRDML); missing catalog/policy release; numeric divergence;
  GSAS-II version mismatch or missing dependency.

## Hard rules

- **Never call a phase "absent" or "unsupported"** without a phase-specific
  detection-limit study (injection/recovery test).
- **No phase purity, structure solution, publication readiness, or robust QPA
  claims from Rwp alone.**
- Context-only candidates **never** become `supported`, even if they
  superficially improve a fit.