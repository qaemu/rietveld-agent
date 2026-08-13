# Documentation

Scientific documentation of the repository: the normative protocol, the
paper series, and the installation reference.

| file | contents |
|---|---|
| `rqpa_protocol.md` | The normative Rietveld QPA protocol (data, structure models, refinement budget, staged ladder, quantification, acceptance criteria, and the spike-15 empirical decisions in §7). |
| `G0-CHECKLIST.md` | Phase-0 gate checklist for the deterministic engine (catalog, calibrations, ingest, verdict, governance). |
| `installation.md` | Installation on macOS / Linux / Windows from OpenCode, Claude Code, or Codex; verification; troubleshooting. |
| `papers/` | The scientific paper series — [P-1 protocol](papers/paper1_protocol.md) (how the software works), [P-2 validation](papers/paper2_validation.md) (why it is scientifically valid), [P-3 deployment](papers/paper3_deployment.md) (how to run it). Index: [`papers/README.md`](papers/README.md). |

The manuscript (TeX source in `../paper/`, compiled with `make paper`)
restates the protocol and results in publication form;
`../references.bib` is the shared bibliography. The spike-by-spike
research log lives in `../notes/`.

## Reading order

1. `../README.md` — what this is and how to reproduce everything.
2. `rqpa_protocol.md` — the protocol itself.
3. `papers/paper1_protocol.md` — how the software works.
4. `papers/paper2_validation.md` — why the results are valid.
5. `papers/paper3_deployment.md` + `installation.md` — how to run it.
6. `../notes/spike15.md` — the empirical probes behind the locked
   decisions; `../notes/spike16.md` — the validation.