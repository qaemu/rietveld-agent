# Documentation

Scientific documentation of the repository: the normative protocol, the
paper series, and the installation reference.

| file | contents |
|---|---|
| `rqpa_protocol.md` | The normative Rietveld QPA protocol (data, structure models, refinement budget, staged ladder, quantification, acceptance criteria, and the unit-15 empirical decisions in §7). |
| `G0-CHECKLIST.md` | Phase-0 gate checklist for the deterministic engine (catalog, calibrations, ingest, verdict, governance). |
| `installation.md` | Installation on macOS / Linux / Windows from OpenCode, Claude Code, or Codex; verification; troubleshooting. |

The paper series lives in [`../paper/`](../paper/) as PDFs (rebuilt with
`make paper`): the manuscript `main.pdf`, P-1 protocol, P-2 validation,
P-3 deployment. The unit-by-unit research log lives in `../notes/`;
`../references.bib` is the shared bibliography.

## Reading order

1. `../README.md` — what this is and how to reproduce everything.
2. `rqpa_protocol.md` — the protocol itself.
3. `../paper/main.pdf` — the manuscript.
4. `../paper/paper1_protocol.pdf` — how the software works.
5. `../paper/paper2_validation.pdf` — why the results are valid.
6. `../paper/paper3_deployment.pdf` + `installation.md` — how to run it.
7. `../notes/unit15.md` — the empirical probes behind the locked
   decisions; `../notes/unit16.md` — the validation.