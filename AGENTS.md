# AGENTS.md — orientation for AI agent runtimes

You are operating in a deterministic scientific repository. Read this
before running anything.

## What this is

`rietveld-agent` is an open, deterministic **Rietveld quantitative phase
analysis (RQPA)** engine demonstrated end-to-end on the four NIST SRM
2686a reference patterns of García-Maté et al. (2024). GSAS-II is the
numerical authority; this repository provides the COD-pinned structure
set, the staged budget-bounded refinement protocol, the evaluation
harness and honest reporting. **Every reported number is reproducible
with `make`** and hash-locked (canonical content md5
`f43d8be2c932420676d612242dd049a5`).

## Entry points

| task | command |
|---|---|
| environment | `make env` |
| syntax/structure checks | `make check` |
| regression suite | `make test` |
| full protocol rerun (6 refinements, GSAS-II ~45 min) | `make report` |
| complete-COD mirror (~26 GB, resumable) + line index | `make cod-tree`, `make cod-index` |
| validation: reproducibility + known-answer + gates | `python3 benchmarks/protocols/validate.py [--skip-rerun\|--skip-synth]` |
| figures, papers | `make figures`, `make paper` |

## Results live in

- `data/unit15/results/unit15_report.json` — refinements
- `data/unit16/results/unit16_report.{json,md}` — validation evidence
- `data/structures/catalog.json` — md5-recorded structure catalogue
- `paper/` — manuscript and paper series as PDFs (P-1 protocol, P-2
  validation, P-3 deployment; `make paper`)

## Invariants — do not violate

1. **Scientific inputs are controlled.** Catalog, calibrations, policies
   and models are hashed, versioned and reviewable. Never silently modify
   `data/structures/catalog.json`, pinned calibrations, or scientific
   policy files; route changes through the documented governance
   (`governance/`).
2. **Local first.** No telemetry, no cloud refinement, no API calls with
   data. Apache-2.0.
3. **Weak context only.** The sample name may add candidate phases; it
   may never remove candidates, improve scores, or confirm a phase.
4. **Defensible uncertainty.** Per-phase statuses are
   `supported` / `inconclusive` / `not_selected` / `out_of_domain` /
   `held` / `failed`. Never produce universal accuracy claims; never
   claim publication readiness from Rwp alone.
5. **Interchangeable hosts.** You are a front-end around the same
   deterministic engine. Your runtime must not change results — verify
   with the content hash rather than by re-running interactively.

## Documentation map

- `README.md` — overview, results, roadmap
- `docs/rqpa_protocol.md` — the normative protocol
- `docs/installation.md` — install on macOS/Linux/Windows via OpenCode /
  Claude Code / Codex
- `paper/` — manuscript and paper series (PDF; `make paper`)
- `notes/` — unit-by-unit research log (unit15: protocol decisions,
  unit16: validation)