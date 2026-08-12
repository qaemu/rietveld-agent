# rietveld-agent

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.10%2B-informational.svg)](pyproject.toml)
[![Manuscript](https://img.shields.io/badge/manuscript-PDF-brightgreen.svg)](paper/main.pdf)
[![Data](https://img.shields.io/badge/data-Zenodo%2010.5281%2Fzenodo.1318501-lightgrey.svg)](https://doi.org/10.5281/zenodo.1318501)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.0000000.svg)](https://doi.org/10.5281/zenodo.0000000)

An open, deterministic toolkit for validating, refining and reporting
**Rietveld quantitative phase analysis (RQPA)** of laboratory powder
X-ray diffraction (PXRD) data — demonstrated end-to-end on the four
**NIST SRM 2686a** reference patterns of García-Maté et al. (2024).
GSAS-II is the numerical authority; this repository provides the
COD-pinned structure set, the staged, budget-bounded refinement
protocol, the evaluation harness, and an honest-uncertainty reporting
layer. Every number in the paper is reproducible with `make`.

> The reference study and the associated reproducibility benchmark:
> García-Maté, M., De la Torre, Á. G., León-Reina, L., & Aranda, M. A. G.
> (2024). *Reproducibility and accuracy of Rietveld quantitative phase
> analysis of NIST SRM 2686a cement clinker.* Cement and Concrete
> Research, 180, 107506.
> doi: [10.1016/j.cemconres.2021.106376](https://doi.org/10.1016/j.cemconres.2021.106376)

## Table of contents

- [Research artifacts](#research-artifacts)
- [Repository layout](#repository-layout)
- [Results](#results)
- [Installation](#installation)
- [Reproducing the results](#reproducing-the-results)
- [Data and structure models](#data-and-structure-models)
- [Methodology](#methodology)
- [Design posture](#design-posture)
- [Roadmap](#roadmap)
- [License](#license)

## Research artifacts

| artifact | where |
|---|---|
| manuscript (this work) | [`paper/main.pdf`](paper/main.pdf) — TeX source in [`paper/`](paper/) |
| protocol specification | [`docs/rqpa_protocol.md`](docs/rqpa_protocol.md) |
| structure catalogue | [`data/structures/catalog.json`](data/structures/catalog.json) |
| full numeric results (md5-locked) | [`data/spike15/results/spike15_report.json`](data/spike15/results/spike15_report.json) |
| spike log & locked decisions | [`notes/spike15.md`](notes/spike15.md) |

Please cite the reference study above (and this repository, via
[`CITATION.cff`](CITATION.cff)) if you use these materials:

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
```

## Repository layout

```
benchmarks/eval/      instrument-aware noise model, evaluation harness
benchmarks/spikes/    numbered, self-contained experiments (spikes 01-15)
core/                 deterministic scientific engine (ingest, calibration,
                      catalog, retrieval, hypothesis, verdict, reporting)
cli/                  operator + expert CLI
admin/                administrator stubs (calibrations, catalog releases)
governance/           schemas and policies for controlled scientific inputs
data/structures/      COD-pinned RQPA structure set (md5-recorded, spike 14)
data/spike15/         RQPA protocol run: results/ (tracked) + work/ (ignored)
data/catalog/         pinned COD catalog release (catalog_0.1.0)
docs/                 protocol specification, phase-0 checklist
paper/                TeX manuscript + figure generation (make paper)
references.bib        bibliography shared by the paper and the docs
notes/                spike-by-spike research log
tests/                regression suite
```

## Results

Spike 15 implements the protocol; all six refinements converge
deterministically (same inputs → same outputs, result-JSON md5 recorded):

| sample | model | wR (%) | note |
|---|---|---|---|
| clinker Cu | M3 | 15.1 | alite over-absorbed (95.3 wt%) |
| clinker Cu | **T1** | 20.9 | **alite 63.5 vs 66.0 published** |
| silicate residue | M3 | 20.2 | alite 98.5 wt% |
| silicate residue | T1 | 27.1 | alite 84.3 vs 78.7 |
| aluminate residue | — | 14.9 | aphthitalite below detection (see protocol) |
| clinker synchrotron | M3 | **9.8** | wR best; alite absorbs |

The remaining distance to the publication-grade targets (wR ≤ 6.5 % Cu,
≤ 5 % sync) is attributed to microabsorption and is the target of the
roadmap's Brindley-correction spike. Figure reproduction:
[`paper/figures/`](paper/figures/).

## Installation

```bash
git clone https://github.com/qaemu/rietveld-agent
cd rietveld-agent
make env          # creates .venv with numpy/scipy/matplotlib
```

GSAS-II is pinned in `.vendor/GSAS-II` and bootstrapped by
`benchmarks/eval/sim.py:ensure_gsasii` (official installer; Apache-2.0,
never bundled as a dependency). `tectonic` is required for the paper
(`brew install tectonic` or your distro's package).

## Reproducing the results

```bash
make report     # rerun every RQPA refinement (GSAS-II, ~45 min)
make figures    # regenerate the manuscript figures
make paper      # compile paper/main.pdf (tectonic)
make check      # syntax + structure checks
make cite       # show the recommended citation
```

The reference patterns are the four SRM 2686a files published with
García-Maté et al. (Zenodo 10.5281/zenodo.1318501); spike-11/12
parsers re-emit them as .xye with intensities untouched (sha256
recorded).

## Data and structure models

- **Patterns** — NIST SRM 2686a: Cu Kα₁ clinker, silicate and
  aluminate residues, and 0.82543 Å synchrotron clinker (ALBA).
- **Structure set** — the published polymorph inventory, re-sourced
  from COD (CC0) and validated by an independent parser: alite M3,
  alite T1, belite β, belite α′H, cubic and orthorhombic aluminate,
  brownmillerite ferrite, periclase, aphthitalite. Every file is
  md5-recorded in `data/structures/catalog.json` (see the README in
  `data/structures/` for provenance and substitutions).

## Methodology

The full specification is [`docs/rqpa_protocol.md`](docs/rqpa_protocol.md):
staged refinement ladder (scales → alite cell → belite cells → minor
cells), bounded budget (Chebyschev background, sample shift, per-phase
scale, major-phase cells), Hill–Howard normalisation
(W_i = S_i·M_i·V_i / Σ S·M·V), acceptance criteria, and the locked
empirical decisions (scale priors on the aluminate residue, aphthitalite
handling, sync T1 exclusion). Spike notes in [`notes/`](notes/) document
every probe behind those decisions.

## Design posture

- **Hybrid scientific workflow, not an LLM over GSAS-II.** GSAS-II is the
  numerical authority for Rietveld refinement; the deterministic core plans and
  executes bounded, policy-approved refinement recipes.
- **Raw data + sample name in; bounded, evidence-rich result out.** The sample
  name is weak context only (it may add `context_only` candidates, never remove
  candidates, never improve scores, never confirm a phase).
- **Defensible uncertainty over false certainty.** Statuses per phase family:
  `supported`, `inconclusive`, `not_selected`, `out_of_domain`, `held`, `failed`.
  No "absent" verdicts without detection-limit studies. No universal accuracy
  numbers. No claims of phase purity, structure solution, publication
  readiness, or robust QPA from Rwp alone.
- **Controlled scientific inputs** (catalog, calibrations, sample-name
  vocabulary, scientific policies, optional models) are hashed, versioned,
  reviewed, replayable, and rollback-able. No update may silently alter the
  meaning of an existing analysis.
- **Local-first, no telemetry, no cloud refinement, Apache-2.0.**
- **AI hosts are interchangeable assistants** (Codex, Claude Code, OpenCode)
  around one deterministic engine, governed by shared versioned contracts.

## Roadmap

| spike | deliverable | state |
|---|---|---|
| 11-12 | PXRD ingest + COD catalog + synthetic simulator | done |
| 13 | multiphase Rietveld QPA (5-phase baseline) | done |
| 14 | published structure set (8 phases + T1), md5-locked | done |
| 15 | publication-grade protocol runner (bounded budget) | done |
| 16 | validation harness (reproducibility, synthetic recovery, gates) | in progress |
| 17 | Brindley microabsorption corrections | planned |

## License

Apache-2.0 for original code ([LICENSE](LICENSE)). GSAS-II is an external
dependency under its own license (Apache-2.0). The reference catalog is a
curated, revisioned subset of the Crystallography Open Database (CC0 /
public domain, attribution preserved in catalog releases). Structure files
under `data/structures/` are CC0 from COD.