# qarr-qpa — quantitative phase analysis project

Landing page for the quantitative-phase-analysis project archived in
this repository (reference-style `projects/` entry).

## Goal

Deterministic, verifiable Rietveld quantitative phase analysis (RQPA)
of laboratory PXRD data driven by agent runtimes, with every reported
number machine-checkable.

## Milestones

| milestone | artifact | status |
|---|---|---|
| Protocol (bounded-budget refinements) | `data/unit15/results/` — [protocol paper](paper/paper1_protocol.pdf) | done |
| Validation (reproducibility + known-answer + gates) | `data/unit16/results/` — [validation paper](paper/paper2_validation.pdf); hash `f43d8be2c932420676d612242dd049a5` | done |
| Full-COD screening + 20-sample QPA gate | `benchmarks/qpa_gate/qpa_gate.py`, results in `data/qpa_gate/results/` | in progress — 5/20 run, 3 PASS (verified 2026-08-15) |
| Brindley microabsorption corrections | planned | planned |

## Composition demo

| pattern | truth (wt%, normalized) | gate |
|---|---|---|
| qarr_1a..1h (corundum/zincite/fluorite ternary) | 1.15–94.81 per mix | see `data/qpa_gate/results/` |
| qarr_2 (ternary + brucite), qarr_3 (+ glass), qarr_4 (corundum/magnetite/zircon) | | |
| bauxite (7-phase oxide/hydroxide suite) | | |
| iron oxides 30/70, 50/50, 70/30, Mexican magnetite | | |
| SRM 2686a clinker suite (Cu + synchrotron) | García-Maté et al. 2024 tables | |

## How to reproduce

```bash
make env && make check && make test
./scripts/gate_runner.sh          # full 20-sample gate
```

## Files of record

- `README.md` / `README.en.md` — project overview
- `docs/rqpa_protocol.md` — normative protocol
- `governance/`, `skills/contracts/` — controlled-input policies
- `notes/` — research log