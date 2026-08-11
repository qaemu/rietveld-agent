# Reporting rules (contracts/0.1.0, draft)

## Mandatory separation

Every result must visibly separate three kinds of content:

| Kind | Source | Label |
|---|---|---|
| Measured observations | Raw data, GSAS-II outputs (observed/calc/background/residual, R factors, esds, covariances) | `measured` |
| Deterministic checks | Calibration match, retrieval, overlap computation, criteria application | `deterministic` |
| AI interpretation | Any language-model-generated explanation | `ai` |

AI-generated interpretation must never be interleaved into the deterministic
scientific records; it is appended in a clearly labeled section.

## Evidence-bundle minimum contents (expert review surface)

1. Supported, inconclusive, and competing phase-family hypotheses, with
   context-only candidates flagged.
2. Blind-retrieval vs name-context candidate provenance.
3. Peak attribution: matched diagnostic reflections, expected-but-missing,
   documented overlaps, important unexplained regions.
4. Observed, calculated, background, and residual profiles at readable scale.
5. Refinement sequence: released/fixed parameters, checkpoints, R factors,
   shifts, uncertainties, correlations, covariance warnings, physical-bound
   warnings.
6. Provenance block: input SHA-256, calibration identity+hash, catalog
   release, policy release, GSAS-II version + environment, model version if
   any, contracts revision.
7. Plain-language abstention / hold / out-of-domain reason.
8. A way for the expert to approve, reject, or document an alternative
   interpretation **without rewriting the historical evidence**.

## Warnings

- Results produced with unreviewed, obsolete, incompatible, or experimental
  controlled inputs must carry a prominent banner.
- "Automated preliminary result" vs "expert-reviewed conclusion" distinction
  must be unambiguous in the artifact itself.

## Privacy

- Raw diffraction arrays remain local unless the user explicitly requests
  their inclusion in an export.
- No telemetry. No cloud refinement. Exports are opt-in and list exactly
  which arrays are included.