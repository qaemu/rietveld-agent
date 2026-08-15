# rietveld-agent (English mirror)

[![License: Apache 2.0](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)
[![Reference data (Zenodo)](https://img.shields.io/badge/data-Zenodo%2010.5281%2Fzenodo.1318501-lightgrey.svg)](https://doi.org/10.5281/zenodo.1318501)
[![OpenCode ready](https://img.shields.io/badge/OpenCode-ready-22C55E.svg)](docs/installation.md)
[![Claude Code ready](https://img.shields.io/badge/Claude%20Code-ready-8AA0FF.svg)](docs/installation.md)
[![Codex ready](https://img.shields.io/badge/Codex-ready-FFA63D.svg)](docs/installation.md)
[![Gate 3/20](https://img.shields.io/badge/full--COD%20gate-3%2F20%20PASS-yellow.svg)](benchmarks/spikes/spike_20_fullcod_qpa.py)

> English mirror of [README.md](README.md). The canonical
> documentation is the English `README.md`; this page keeps the
> reference-layout entry points of the Light-skills repository style
> (badges, quickstart, entry-point table, skill map) for readers who
> arrive from either README.

## What this is

A deterministic Rietveld quantitative phase analysis (RQPA) engine for
agent runtimes (OpenCode / Claude Code / Codex), demonstrated
end-to-end on the NIST SRM 2686a cement-clinker reference patterns:

- **GSAS-II is the numerical authority** — pinned in `.vendor/`,
  bootstrapped on first use, never bundled as a dependency.
- **COD-pinned structures** — `data/structures/catalog.json`,
  md5-recorded, governed (see `governance/`).
- **Machine-checked reproducibility** — canonical content hash
  `f43d8be2c932420676d612242dd049a5`; every table in the README
  rebuilds with `make report`.

## Quick start

```bash
git clone https://github.com/qaemu/rietveld-agent-spikes
cd rietveld-agent-spikes
make env       # virtualenv (numpy/scipy/matplotlib/pytest)
make check     # syntax + structure checks
make test      # regression suite
```

Per-platform / per-runtime instructions: [`docs/installation.md`](docs/installation.md).

## Entry points

| current need | entry point |
|---|---|
| I want to run an analysis | `python -m cli analyze <pattern>` |
| I want the full-COD screening gate (20 samples) | `python benchmarks/spikes/spike_20_fullcod_qpa.py` |
| I want the validation evidence | `python benchmarks/spikes/spike_16_validate.py` |
| I want the reproducible protocol run | `make report` (~45 min) |
| I want the papers | `make paper` (needs tectonic) |
| I want orientation for an agent runtime | [`AGENTS.md`](AGENTS.md) |

## Skill map (repository layout)

| module | contents |
|---|---|
| `core/` | deterministic engine: ingest, calibration, catalog, codsearch, hypothesis, verdict, reporting |
| `benchmarks/` | evaluation harness (`eval/`) and numbered experiments (`spikes/`) |
| `cli/` | operator CLI |
| `governance/` + `skills/contracts/` | schemas/policies for controlled scientific inputs and agent decision criteria |
| `data/` | benchmark corpus, COD index, structure set, protocol/validation results |
| `docs/` `paper/` `notes/` | protocol, installation, manuscript + paper series, research log |
| `tests/` | regression suite (`make test`) |

## Status of the full-COD gate (spike 20)

The spike-20 harness ranks candidate phases against a line-index search
of the complete COD (524,948 entries), strips the winners, and refines a
staged GSAS-II QPA for each accepted hypothesis. See
[`benchmarks/spikes/spike_20_fullcod_qpa.py`](benchmarks/spikes/spike_20_fullcod_qpa.py)
and `notes/spike18.md` for the screening design; the gate result table
lands in `data/spike20/results/`.

Verified state (re-verified 2026-08-14 with
`benchmarks/spikes/spike_20_aggregate.py`; 5 of 20 samples run, 15 not
run):

| sample | verdict | wR | inferred phases (wt%) |
|---|---|---|---|
| `qarr_1a` | PASS | 42.85 | fluorite 94.7, zincite 3.7, corundum 1.6 |
| `qarr_1e` | PASS | 38.50 | corundum 57.6, fluorite 28.2, zincite 14.2 |
| `qarr_1h` | PASS | 38.57 | corundum 37.3, fluorite 33.8, zincite 28.9 |
| `qarr_1f` | FAIL | 50.75 | single phase 2300112 @ 100.0 — corundum/zincite/fluorite all MISSING |
| `iron_30_70` | FAIL | 5.07 | magnetite 100.0 — hematite MISSING (truth 31.8) |

Known pipeline failure mode: samples can converge on a single candidate
phase at 100 wt% (screen ranks non-mineral CIFs above the true phases;
the multi-phase loop then never adds them).

## License

Apache-2.0 for original code; GSAS-II external under its own Apache-2.0
license; structure files CC0 from the Crystallography Open Database.
See [LICENSE](LICENSE) and [CITATION.cff](CITATION.cff).