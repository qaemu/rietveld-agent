# rietveld-agent (English mirror)

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)
[![Reference data (Zenodo)](https://img.shields.io/badge/data-Zenodo%2010.5281%2Fzenodo.1318501-lightgrey.svg)](https://doi.org/10.5281/zenodo.1318501)
[![OpenCode ready](https://img.shields.io/badge/OpenCode-ready-22C55E.svg)](docs/installation.md)
[![Claude Code ready](https://img.shields.io/badge/Claude%20Code-ready-8AA0FF.svg)](docs/installation.md)
[![Codex ready](https://img.shields.io/badge/Codex-ready-FFA63D.svg)](docs/installation.md)
[![Gate 3/5](https://img.shields.io/badge/full--COD%20gate-3%2F5%20PASS-yellow.svg)](benchmarks/qpa_gate/qpa_gate.py)

> English mirror of [README.md](README.md). The canonical
> documentation is the English `README.md`; this page keeps the
> entry-point tables and the per-sample gate verdicts for readers who
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
git clone https://github.com/qaemu/rietveld-agent
cd rietveld-agent
make env       # virtualenv (numpy/scipy/matplotlib/pytest)
make check     # syntax + structure checks
make test      # regression suite
```

Per-platform / per-runtime instructions: [`docs/installation.md`](docs/installation.md).

## Entry points

| current need | entry point |
|---|---|
| I want to run an analysis | `python -m cli analyze <pattern>` |
| I want the full-COD screening gate (20 samples) | `python benchmarks/qpa_gate/qpa_gate.py` |
| I want the validation evidence | `make test` — full known-answer harness documented in `AGENTS.md` |
| I want the reproducible protocol run | `make report` (~45 min) |
| I want the papers | `make paper` (needs tectonic) |
| I want orientation for an agent runtime | [`AGENTS.md`](AGENTS.md) |

## Repository layout

| module | contents |
|---|---|
| `core/` | deterministic engine: ingest, calibration, catalog, codsearch, hypothesis, verdict, reporting |
| `benchmarks/` | evaluation harness (`eval/`), QPA gate (`qpa_gate/`) |
| `cli/` | operator CLI |
| `governance/` | schemas/policies for controlled scientific inputs |
| `data/` | benchmark corpus, COD index, structure set, protocol/validation results |
| `docs/` `paper/` `notes/` | protocol, installation, manuscript + paper series, research log |
| `tests/` | regression suite (`make test`) |

## Status of the full-COD gate

The gate harness ranks candidate phases against a line-index search of
the complete COD (524,948 entries), strips the winners, and refines a
staged GSAS-II QPA for each accepted hypothesis. See
[`benchmarks/qpa_gate/qpa_gate.py`](benchmarks/qpa_gate/qpa_gate.py)
(entries indexed by `core/codsearch`); the gate result table lands in
`data/qpa_gate/results/`.

Verified state (re-verified 2026-08-15 with
`benchmarks/qpa_gate/aggregate.py`; 5 of 20 samples run, 15 not run):

| sample | verdict | wR | inferred phases (wt%) |
|---|---|---|---|
| `qarr_1a` | PASS | 42.85 | fluorite 94.7, zincite 3.7, corundum 1.6 |
| `qarr_1e` | PASS | 38.50 | corundum 57.6, fluorite 28.2, zincite 14.2 |
| `qarr_1h` | PASS | 38.57 | corundum 37.3, fluorite 33.8, zincite 28.9 |
| `qarr_1f` | FAIL | 38.23 | corundum 46.1, zincite 27.2, fluorite 26.6 (truth 27.1/55.2/17.7) |
| `iron_30_70` | FAIL | 5.07 | magnetite 100.0 — hematite MISSING (truth 31.8) |

Documented failure causes (root-cause analysis with quantified
CIF-combination and profile sweeps in
[`benchmarks/qpa_gate/sweeps/README.md`](benchmarks/qpa_gate/sweeps/README.md),
reproducible via `benchmarks/qpa_gate/sweep_cifs.py` and
`benchmarks/qpa_gate/sweep_intensity.py`):
`iron_30_70` — the data's magnetite matches only the condensed
a=8.3582 CIF (2300616; five ambient Fe3O4 COD entry cells all fit
worse), and the F2-consistent joint-fit split is pinned at 38/62 vs
truth 32/68 (Δ6.2); further, the pipeline's profile-stage forward
selection drops hematite entirely (single phase, magnetite 100%);
`qarr_1f` — the zincite CIF's high-angle lines are ~2–3.4× weaker than
the data's (monotonic angle-dependent intensity mismatch; the zincite
CIF choice alone swings the refined zincite wt% from 14 to 75 around
truth 55.2, with equal-wR fits bracketing the truth: 23.6% and 73.5%),
and free Uiso / March-Dollase PO (applied correctly via GSAS-II
`Pref.Ori.`) do not fix it — the Stage-A profile refinement also runs
away (U,V,W,X,Y SVD singularity, shift −40°).

## License

Apache-2.0 for original code; GSAS-II external under its own Apache-2.0
license; structure files CC0 from the Crystallography Open Database.
