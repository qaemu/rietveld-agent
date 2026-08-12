# Documentation

This directory holds the scientific documentation of the repository.

| file | contents |
|---|---|
| `rqpa_protocol.md` | The normative Rietveld QPA protocol (data, structure models, refinement budget, staged ladder, quantification, acceptance criteria, and the spike-15 empirical decisions in §7). |
| `G0-CHECKLIST.md` | Phase-0 gate checklist for the deterministic engine (catalog, calibrations, ingest, verdict, governance). |

The manuscript (TeX source in `../paper/`, compiled with
`make paper`) restates the protocol and results in publication form;
`../references.bib` is the shared bibliography. The spike-by-spike
research log lives in `../notes/`.

## Reading order

1. `../README.md` — what this is and how to reproduce everything.
2. `rqpa_protocol.md` — the protocol itself.
3. `../paper/main.pdf` — the manuscript.
4. `../notes/spike15.md` — the empirical probes behind the locked decisions.
