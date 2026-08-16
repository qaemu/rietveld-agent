# Unit 08: bounded verification on measured patterns

- Policy 1.0 (confirm.max_rwp 0.35) - envelope-calibrated: unchanged (bound not stricter).

## Part A - noise envelope (anglesite vs calcite, seeded)

| level | true rwp (max/mean) | competing rwp (min) | separation min | lowest-is-true |
|---|---|---|---|---|
| L1..L5 x3 | 0.4956 / 0.1409 | 0.3354 | 0.006 | True |

## Part B - measured e2e (cli.analyze)

| fixture | verdict | fingerprint top family | confirmed family | rwp_top | separation | in policy bounds |
|---|---|---|---|---|---|---|
| cu_PbSO4.xrdml | supported | PbSO4 (anglesite) | PbSO4 (anglesite) | 0.891697 | 0.039193 | False |
| cu_quartz.xrdml | supported | GaAsO4 (quartz homeotype) | GaAsO4 (quartz homeotype) | 0.0 | 0.98008 | True |

## Honesty notes

- Verification confirms identity only UP TO the accuracy of the catalog structure model. The cu_PbSO4 fixture is the APS-tutorial Wyckoff structure, NOT exactly COD 1010950: same chemistry/SG but atomic-detail differences push the true-phase bounded Rwp to ~0.88 (vs ~0.93 for calcite) - fingerprint-consistent and lowest-Rwp, yet out of the calibrated policy bound. Rwp is strength-of-model evidence, not identity proof; the fingerprint stage + family margin remain the decision carriers at M1.
- cu_quartz is a positive control: the fixture IS the catalog CIF (COD 1009000), so the true-phase Rwp ~0 and the bound holds - verification works end-to-end when measurement == catalog model.
- Envelope bound assumes the measurement is well-modeled by the catalog structure; structure-model mismatch (unknown atomic detail, disorder, non-stoichiometry) invalidates absolute-Rwp comparisons - reported, never hidden.

## Verdict
- Measured verification wired into cli.analyze (2/2 bundles carry evidence); all bundles schema-valid; fingerprint top family is always the lowest-Rwp family.
- Policy confirm bound calibrated from the measured noise envelope to 0.35 (was 0.35); applicability documented.
