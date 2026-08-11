# Spike 03: calibration registry

- records: 3 (Cu XYE-DEFAULT, Cu XML-FINGERPRINTED, Fe XML-FINGERPRINTED)
- schema: governance/schemas/calibration.schema.json (calibration/v0), validated on insert
- policy: calibration_requirement = released-only (analysis_policy.example.json)

## PRM audit

| prm | anode | kalpha1 | kalpha2 | profile fn | sha256 |
|---|---|---|---|---|---|
| cu | CuKa | 1.5405 | 1.5443 | 3 | e59413059bc3 |
| fe | FeKa | 1.9360 | 1.9399 | 3 | c64443a939bf |

## Resolution of golden fixtures

| fixture | expected | outcome | record |
|---|---|---|---|
| cu_PbSO4.xrdml | Cu | resolved | cal-967c5da10ad2 |
| cu_quartz.xrdml | Cu | resolved | cal-967c5da10ad2 |
| fe_PbSO4.xrdml | Fe | resolved | cal-9bba7cdaeb39 |
| unknown_MoKa | unknown | unknown | - |
| ambiguity_demo | ambiguous | ambiguous | ['cal-967c5da10ad2', 'cal-ambig-demo'] |

## Findings
- Resolution is physical (anode + kalpha wavelengths + scan axis), not filename/extension based; scan-grid is soft context (real scans vary).
- ambiguous and unknown are hard stops: the agent must abstain or ask the administrator to register the instrument; never guess a PRM.
- evidence_ref ties every released calibration to its verification run (spike01 report) -> audit trail is part of the record.
- XYE-DEFAULT is the single admin-approved fallback for plain XYE input; metadata-rich XRDML input always resolves via fingerprints.

## Verdict
- [ ] all 3 golden fixtures resolve to the correct released calibration
- [ ] unknown instrument detected as unknown (Mo variant)
- [ ] duplicate released calibrations detected as ambiguous
- [ ] registry persists and round-trips identically
- [ ] every record validates against calibration.schema.json
