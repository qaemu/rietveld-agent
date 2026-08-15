# rietveld-agent

[![License: Apache 2.0](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)
[![Reference data (Zenodo)](https://img.shields.io/badge/data-Zenodo%2010.5281%2Fzenodo.1318501-lightgrey.svg)](https://doi.org/10.5281/zenodo.1318501)
[![OpenCode ready](https://img.shields.io/badge/OpenCode-ready-22C55E.svg)](docs/installation.md)
[![Claude Code ready](https://img.shields.io/badge/Claude%20Code-ready-8AA0FF.svg)](docs/installation.md)
[![Codex ready](https://img.shields.io/badge/Codex-ready-FFA63D.svg)](docs/installation.md)
[![Gate 3/5](https://img.shields.io/badge/full--COD%20gate-3%2F5%20PASS-yellow.svg)](benchmarks/qpa_gate/qpa_gate.py)

Deterministic Rietveld quantitative phase analysis (RQPA) of laboratory
powder X-ray diffraction (PXRD) data, demonstrated end-to-end on the four
NIST SRM 2686a cement clinker reference patterns of García-Maté et al.
(2024). GSAS-II performs the refinement; this repository supplies the
COD-pinned structure set, a staged and budget-bounded refinement
protocol, the evaluation harness, and hash-locked reporting.

Every reported number rebuilds from the raw patterns with `make report`,
and the validation harness re-runs the whole pipeline and
compares canonical content hashes — reproducibility is a machine-checked
property, not a claim. The deterministic engine is runtime-agnostic:
OpenCode, Claude Code, and Codex are interchangeable front-ends around it
([`AGENTS.md`](AGENTS.md)).

## Quick start

```bash
git clone https://github.com/qaemu/rietveld-agent-qpa
cd rietveld-agent
make env          # virtualenv with numpy/scipy/matplotlib/pytest
make check        # syntax checks
make test         # regression suite
```

GSAS-II is pinned in `.vendor/GSAS-II` and bootstrapped on first use
(official installer, Apache-2.0, never bundled as a dependency) — no
manual GSAS-II installation. Full per-platform, per-runtime instructions
(macOS / Linux / Windows, OpenCode / Claude Code / Codex):
[`docs/installation.md`](docs/installation.md).

## Results

The protocol runner implements the staged refinement; all six refinements converge
deterministically (same inputs → same outputs, result-JSON md5 recorded):

| sample | model | wR (%) | note |
|---|---|---|---|
| clinker Cu | M3 | 15.07 | alite over-absorbed (95.3 wt%) |
| clinker Cu | **T1** | 20.91 | **alite 63.5 vs 66.0 published** |
| silicate residue | M3 | 20.18 | alite 98.5 wt% |
| silicate residue | T1 | 27.11 | alite 84.3 vs 78.7 |
| aluminate residue | M3 | 13.56 | all 5 phases reported (aphthitalite constraint) |
| clinker synchrotron | M3 | **9.76** | best wR; alite absorbs |

The validation harness checks the suite: two independent full runs produce
bit-identical payloads (canonical md5 `f43d8be2c932420676d612242dd049a5`),
and four synthetic known-answer patterns recover every injected phase
within band — **4/4 samples, 24/24 phase rows** — no refinement fails,
no phase is left indeterminate. The residual distance to the
publication-grade targets (wR ≤ 6.5% Cu, ≤ 5% sync) is attributed to
microabsorption (planned Brindley-correction work).

The full-COD screening gate ranks candidate phases by strip-iteration
against the complete COD line index (524,948 entries), refines the
accepted hypotheses with staged GSAS-II QPA, and reports a verified
phase set per benchmark sample (`make gate`; state: **3/5 PASS — 2 FAIL
honestly documented** — `iron_30_70` and `qarr_1f`, with root causes in
`data/qpa_gate/results/`; see [README.en.md](README.en.md) for the
reference-style entry table).

## Papers

Four working papers in [`paper/`](paper/) (PDF, rebuilt with `make
paper`), written following the fifteen-step framework of Drake & Han
(2025), *How to write a scientific paper in fifteen steps*
(doi:10.1371/journal.pcbi.1013505):

| paper | contents |
|---|---|
| [`paper/main.pdf`](paper/main.pdf) | manuscript: protocol, results, reproducibility |
| [`paper/paper1_protocol.pdf`](paper/paper1_protocol.pdf) | P-1 — the bounded-budget protocol: how the software works |
| [`paper/paper2_validation.pdf`](paper/paper2_validation.pdf) | P-2 — known-answer validation: why the results are valid |
| [`paper/paper3_deployment.pdf`](paper/paper3_deployment.pdf) | P-3 — deployment on OpenCode, Claude Code, Codex |

## Repository layout

```
core/                 deterministic scientific engine (ingest, calibration,
                      catalog, retrieval, hypothesis, verdict, reporting)
benchmarks/eval/      instrument-aware noise model, evaluation harness
benchmarks/qpa_gate/  full-COD screening + QPA gate (qpa_gate.py, aggregate.py)
cli/                  operator CLI (python -m cli analyze ...)
governance/           schemas and policies for controlled scientific inputs
skills/contracts/     agent-facing decision criteria and parameter allowlists
_shared/              shared conventions for reuse across skills/experiments
scripts/              repository entry points (e.g. gate_runner.sh)
examples/ e2e-qpa/    runnable single-sample demo
projects/ qarr-qpa/   project landing page (composition demo, milestones)
assets/               repository preview figures
data/structures/      COD-pinned RQPA structure set (md5-recorded)
data/qpa_gate/        full-COD gate: results/ (tracked)
docs/                 protocol specification, installation reference
paper/                manuscript + paper series (PDF, make paper)
notes/                research log
tests/                regression suite (make test)
```

## Documentation

- [`docs/rqpa_protocol.md`](docs/rqpa_protocol.md) — the normative protocol
  (data, structure models, refinement budget, staged ladder,
  quantification, acceptance criteria)
- [`docs/installation.md`](docs/installation.md) — installation per
  platform and agent runtime
- [`AGENTS.md`](AGENTS.md) — orientation and invariants for agent runtimes
- [`notes/`](notes/) — research log (protocol decisions, validation)
- [`governance/`](governance/) — schemas and policies for controlled
  scientific inputs

## Reproducing the results

```bash
make report     # rerun every RQPA refinement (GSAS-II, ~45 min)
make figures    # regenerate the manuscript figures
make paper      # rebuild all paper PDFs (tectonic)
make test       # regression suite (pytest)
```

Reference patterns are the four SRM 2686a files published with
García-Maté et al. (Zenodo 10.5281/zenodo.1318501); the ingest parsers
re-emit them as two-column data with intensities untouched (sha256
recorded).

## Citation

If you use this repository, cite the reference study and the repository
itself ([`CITATION.cff`](CITATION.cff)):

```bibtex
@article{garcia-mate2024srms,
  author  = {Garc{\'i}a-Mat{\'e}, M. and De la Torre, {\'A}. G. and
             Le{\'o}n-Reina, L. and Aranda, M. A. G.},
  title   = {Reproducibility and accuracy of Rietveld quantitative phase
             analysis of {NIST} {SRM} 2686a cement clinker},
  journal = {Cement and Concrete Research},
  volume  = {180},
  pages   = {107506},
  year    = {2024},
  doi     = {10.1016/j.cemconres.2021.106376}
}

@misc{rietveld-agent,
  author = {qaemu},
  title  = {{rietveld-agent}: a deterministic Rietveld {QPA} engine for
            scientific agent runtimes},
  year   = {2026},
  url    = {https://github.com/qaemu/rietveld-agent-qpa}
}
```

## Roadmap

| deliverable | state |
|---|---|
| PXRD ingest + COD catalog + synthetic simulator | done |
| multiphase Rietveld QPA (5-phase baseline) | done |
| published structure set, md5-locked | done |
| publication-grade protocol runner (bounded budget) | done |
| validation harness (reproducibility, synthetic recovery, gates) | done |
| complete-COD screening (`core.codsearch` line index over the full COD
  CIF tree, `make cod-tree` ~26 GB + `make cod-index`; `--full-cod` in
  `cli analyze`; candidate filter only, not identification) | in progress |
| full-COD screening + 20-sample QPA gate (`make gate`; 5 samples run,
  3/5 PASS, 2 FAIL documented: iron_30_70, qarr_1f) | in progress |
| Brindley microabsorption corrections | planned |

## License

Apache-2.0 for original code ([LICENSE](LICENSE)). GSAS-II is an external
dependency under its own Apache-2.0 license. Structure files under
`data/structures/` are CC0 from the Crystallography Open Database
(attribution preserved in catalog releases).
