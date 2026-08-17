# Unit 17 — gate completion, parallel run, and the two known FAILs

Date: 2026-08-17. 8h autonomous round. Everything below reproduces from
`benchmarks/qpa_gate/qpa_gate.py` (+ the sweep JSONs in
`data/qpa_gate/sweeps/`) at commit `4675dad`.

## 1. Gate: 5/20 -> 20/20 attempt (15 samples launched)

The 15 remaining samples were launched in three parallel batches of five
(screening + staged QPA; the machine is 8 cores), with dedicated logs:

| batch | samples | log |
|---|---|---|
| A | qarr_1b 1c 1d 1g, qarr_2 | `data/qpa_gate/runA.log` |
| B | qarr_3, qarr_4, bauxite, iron_50_50, iron_70_30 | `data/qpa_gate/runB.log` |
| C | iron_mexican, srm_Clinker_Nist_CuKalpha1_R1, srm_Silicate_enriched_residue_Nist_CuKalpha1_R1, srm_aluminate_enriched_residue_clinkerNIST_180718_R1, srm_Clinker_Synchrotron | `data/qpa_gate/runC.log` |
| D | iron_30_70 + qarr_1f **re-runs** (`--step gate`, pre-edit code) | `data/qpa_gate/runD.log` |

Safety: batches are independent (per-sample `results/<sid>.json`,
per-sample work dirs, shared read-only COD index).  One real hazard was
found and fixed: `download_cif` wrote the CIF cache non-atomically, so
two parallel batches fetching the same CIF could corrupt it (the
`st_size > 100` cache check would then accept a torn file).  Now written
via pid-suffixed temp + atomic rename (`benchmarks/protocols/cod_full.py`).

## 2. Cache-skip trap (why the first re-run did nothing)

`main()` skips any sample whose result JSON exists (`resj.exists()`),
unless `--step gate` is passed.  Launching `--only iron_30_70 qarr_1f`
without the flag printed `qarr_1f: cached / iron_30_70: cached` and
completed in seconds, emitting only the all-manifest summary.  Re-runs
must use `--step gate`.

## 3. What the sweeps already prove about the two FAILs

**iron_30_70** (`data/qpa_gate/sweeps/iron_30_70_cif_sweep.json`):

- Joint fits (hematite 5910082 + each magnetite COD entry, free
  project, no seeds): 2300616 -> wR 100 runaway; 1011032 -> 97.8/2.2
  (wR 7.24); 9002320 -> magnetite-only (wR 6.63); 9002321 -> 48.1/51.9
  (wR 4.45).  Truth is 68.2/31.8 magnetite/hematite: every sensible
  combo either drops hematite or inverts the split.  The F2-consistent
  joint split is pinned at 38/62 (Delta 6.2 absolute) — the recorded
  failure cause.  Conclusion: **no COD magnetite cell lets the gate pass
  the +-3 criterion for this sample**; the condensed a=8.3582 cell
  (2300616) remains the best data-driven choice (best magnetite-only
  fit, wR 5.07).

**qarr_1f** (`data/qpa_gate/sweeps/qarr_1f_cif_sweep.json`): 8 zincite
COD entries tested in fixed corundum/fluorite combos. Converged fits:
zincite 27.2 (wR 38.2, the recorded result), 28.2 (41.8), 34.2 (50.2);
the 63.8 entries (2300450/9004178/9004179) are unseeded-runaway
artifacts (wR 100).  Truth zincite = 55.2.  Conclusion: **an alternate
same-canon-CIF pass cannot reach truth** — the COD zincite entries are
systematically 2-3.4x weak at high angle (documented root cause), and
the honest FAIL stands.

Consequence: no default-pipeline code change can flip either sample;
the fix scope narrows to robustness + characterization.

## 4. Engine changes (all opt-in; default behavior byte-identical)

- `RQPA_STAGE_A_FALLBACK=1` — Stage-A solo profile refinement retries
  a simpler parameter set (U,V,W,X,Y -> V,W -> W -> scale-only) when the
  full refinement hits an SVD singularity / runs away (qarr_1f
  zincite-class, shift -40 deg).
- `RQPA_PROF_B2=1` — one additional canonical re-test round after Stage
  B2 with the instrument profile FREED (seeded from the base), then
  re-freeze the accepted joint profile.  Motivation: B2's frozen-profile
  gate cannot see wR gains that require a profile re-shape; the sweeps
  show free-profile joint fits do improve wR for some magnetic combos.
- `benchmarks/protocols/cod_full.py` — atomic CIF downloads (above).

## 5. Planned experiments (this round)

- runD baseline re-run of iron_30_70/qarr_1f under pre-edit code with
  full stage logs -> confirms current B2 behavior.
- runE: iron_30_70 + qarr_1f under `RQPA_STAGE_A_FALLBACK=1` (and
  iron_30_70 additionally under `RQPA_PROF_B2=1`), results compared to
  runD baselines.

## 6. Housekeeping (same round)

- `.github/workflows/ci.yml` — `make env` + `make check` + `make test`
  on ubuntu (GSAS-II auto-bootstraps; ~2-3 min suite).
- `LICENSE` restored (Apache-2.0; removed with the scaffold strip in
  `86cf079`, `pyproject.toml` always declared it).
- Gate run logs are transient: `.gitignore` gains `data/qpa_gate/*.log`.

## 7. Expected outcome

- 12-20/20 samples with result JSOns (machine allows ~15-20h of
  refinement wall; wall budget 8h).
- iron_30_70 / qarr_1f: PASS is not achievable (section 3); expected
  outcome is either an unchanged FAIL or a better-characterized FAIL
  (e.g. iron with both phases present, split pinned by wR/F2).

## 8. Observed noise (running)

- COD entry **1001390** is present and complete in the CIF cache but
  GSAS-II cannot read it ("No reader could read file"; the phase is
  skipped from the hypothesis — the pipeline degrades gracefully).  If
  a truth phase ever maps to this entry, the sample FAILs with that
  phase missing; noted here to distinguish model failure from data
  failure.