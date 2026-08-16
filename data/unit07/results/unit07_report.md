# Unit 07: bounded Rietveld verification

- Policy: `refinement-budget.v1.json` (recipe `refinement-verify-v1`, version 1.0) -- bounded keys {'histogram': ['Background', 'Sample Parameters Shift'], 'phase': ['Cell'], 'hap': ['Scale']}; prohibited Atoms, Mustrain, Size, HStrain, PhaseFraction, Pref.Ori., LeBail.
- Release `0.1.0` manifest 7e31c47f91ce7f1b...; library manifest 340127a8343ea8f0...

| case | fingerprint top family | rwp (per candidate) | confirmed |
|---|---|---|---|
| cu_PbSO4.xrdml | PbSO4 (anglesite) | 1010950 PbSO4 (anglesite): Rwp=0.0000 GoF=0.00<br>1010928 CaCO3 (calcite): Rwp=0.9289 GoF=834.75 | yes |
| cu_quartz.xrdml | GaAsO4 (quartz homeotype) | 1009000 GaAsO4 (quartz homeotype): Rwp=0.0000 GoF=0.00<br>9012601 SiO2 (quartz): Rwp=0.9808 GoF=765.06 | yes |

## Honesty notes

- Noiseless protocol: the observed pattern is the deterministic GSAS-II simulation of the catalog CIF (same protocol as catalog fingerprints). Noise robustness is unit 05's mandate; real-data Rwp values will be larger -- the decision signal is lowest-Rwp family + policy bounds, not absolute Rwp.
- Bounded keys only (background + shift + cell + scale): a wrong phase cannot hide by absorbing mismatch (atoms/microstrain/size/phase-fraction/LeBail prohibited in policy).
- Rexp is approximate: (N-P)/sum(w y^2) with P from the policy budget; GoF = Rwp/Rexp. No phase-purity/QPA claims are made.

## Verdict
- 2/2 cases confirmed: fingerprint top family reproduced by the bounded refinement (lowest Rwp, within policy bounds).
- Evidence level of the verification stage: **fingerprint + refinement** (bundle schema run_bundle/v0).
