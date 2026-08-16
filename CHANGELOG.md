# Changelog

All notable changes to this repository are documented here, in
[Keep a Changelog](https://keepachangelog.com/) format; the repository
does not use semantic version tags yet (see `pyproject.toml` for the
package version).

## [Unreleased]

### Changed
- `README.md` / `README.en.md` rewritten project-first: hero + results
  preview figures from `assets/`, template-scaffolding vocabulary
  removed, entry-point and gate tables kept.
- `governance/schemas/policy.schema.json` — description strings that
  pointed into the removed `skills/contracts/` now point at
  `docs/rqpa_protocol.md`.

### Removed
- Template scaffolding from the Light-skills restyle, superseded by the
  project's own docs: `skills/contracts/` (agent-facing criteria, now
  covered by `AGENTS.md` + `docs/rqpa_protocol.md`), `docs/G0-CHECKLIST.md`
  (phase-0 milestone doc, now covered by `docs/index.md`),
  `notes/qpa_gate_sweep.md` (analysis preserved in
  `benchmarks/qpa_gate/sweeps/README.md` + `data/qpa_gate/sweeps/`).
- Boilerplate files that shipped with the restyle or pointed at removed
  paths: `CONTRIBUTING.md`, `SECURITY.md`, `LICENSE`, `CITATION.cff`
  (bibliographic metadata remains in `pyproject.toml` and the README
  citation block).

### Changed
- Unit-18/19/20 unified into one maintained module: the full-COD
  screening gate now lives in `benchmarks/qpa_gate/qpa_gate.py` with
  `benchmarks/qpa_gate/aggregate.py` (results under `data/qpa_gate/`);
  `notes/unit18.md` removed (its codsearch design notes live with
  `core/codsearch.py`); the repository landing pages (`README.md`,
  `README.en.md`) no longer expose the unit numbering.

### Added
- `benchmarks/qpa_gate/qpa_gate.py` — full-COD screening and QPA gate
  harness: 20-sample manifest (qarr 1a-1h + 2/3/4 + bauxite,
  iron oxide 30/70-50/50-70/30 + Mexican magnetite, SRM 2686a clinker
  suite), strip-based screening, staged GSAS-II refinement, verifiable
  gate verdicts; `benchmarks/qpa_gate/aggregate.py` re-verifies the
  result jsons against the current gate policy.
- `core/codsearch.py` — line-index search over the complete COD CIF tree
  (`make cod-tree` / `make cod-index`).
- `tests/test_cod_full.py` — regression tests for the full-COD path.
- `data/benchmark/` — benchmark pattern corpus (qarr synthetic QPA set,
  iron oxide mixtures, SRM 2686a inputs).
- Repository restyled to the Light-skills layout: `_shared/`,
  `assets/`, `examples/`, `projects/`, `scripts/` plus `CHANGELOG.md`,
  `CONTRIBUTING.md`, `SECURITY.md`, `README.en.md`.

### Fixed
- `benchmarks/protocols/unit_20_fullcod_qpa.py`: `canon_of` gained a
  formula-element fallback (COD records whose `mineral`/`chemname` are
  empty, e.g. 2300112 `- O Zn -` → zincite, 2300616 `- Fe3 O4 -` →
  magnetite); `canon_recall` and the screening call sites pass the
  formula too, so empty-metadata entries are no longer lost to the
  gate (`zincite` was previously reported `MISSING` for qarr_1f even
  though it was the fitted phase at 27.2 wt%).
- `benchmarks/protocols/unit_20_fullcod_qpa.py` manifest: qarr sample ids
  aligned with `QARR_TRUTH` keys (`qarr_1a`..`qarr_1h`); SRM inputs
  resolve to their real extensions (`.xrdml` / `.dat`).
- Complete-COD runs now execute with the repository virtualenv
  (`.venv/bin/python`, Makefile `PY`) so `gemmi`-dependent CIF reading
  (`cif_calc_lines`) works and the strong-line exclusion set `sel_d` is
  populated; previously a missing `gemmi` in the ambient interpreter
  silently disabled the overlap/“distinctive”/mask defenses during
  screening.

### Removed
- Dead or non-working scratch code: `periclase_diag.py`; un-numbered
  diagnostic scripts that no longer run against the current engine.

### Fixed
- README/README.en unit-20 gate badge claimed "20/20 PASS" while only 3
  of 20 samples were verified; badge now shows the honest verified state
  (3/20 PASS, 2026-08-14) and `README.en.md` documents the per-sample
  verdicts, including the known single-phase-100% convergence failure
  mode (e.g. `qarr_1f`, `iron_30_70`).

## [0.2.0] - 2026

### Added
- `core/codsearch.py` groundwork (line index schema, sharded build).
- Unit-20 harness files (see Unreleased — historical entry kept for
  provenance).

## [0.1.0] - 2025

### Added
- Deterministic RQPA engine: ingest, calibration registry, structure
  catalog, hypothesis ranking, verdict and reporting
  (`core/`, `cli/`).
- Unit-by-unit research log (units 01-16) with reproducible
  experiments, regression tests (`make test`) and the md5-locked
  protocol + validation reports (unit 15, unit 16).
- Paper series P-1..P-3 and manuscript (`paper/`, `make paper`).
- Governance schemas and policies for controlled scientific inputs
  (`governance/`), agent-facing contracts (`skills/contracts/`).