# e2e QPA demo

Minimal end-to-end run of a single quantitative phase analysis from a
raw pattern to the verified gate verdict, using the full-COD gate
harness (screening over the complete COD line index + staged GSAS-II
QPA).

## Run (one sample)

```bash
./scripts/gate_runner.sh --only qarr_1a
```

`qarr_1a` is a synthetic laboratory pattern of corundum + zincite +
fluorite (truth: 1.15 / 4.04 / 94.81 wt%, normalized). Expected flow:

1. `screen_and_rank` — strip-iteration screening against the COD index
   (524,948 entries), zincite + fluorite family recovered first;
2. staged GSAS-II QPA of the accepted hypothesis;
3. gate verdict printed for the sample and written to
   `data/qpa_gate/results/qarr_1a.json`.

## Expected output shape

```
--- qarr_1a ---
qarr_1a: [('corundum', 1.2), ('zincite', 4.0), ('fluorite', 94.6)] rwp=… → PASS
```

Phase names are canonical mineral keys (see `PHASE_CANON` in
`benchmarks/qpa_gate/qpa_gate.py` and the contracts under
`skills/contracts/`); normalized wt% are reported.

## Prerequisites

- `make env` (repository virtualenv; gemmi available there)
- COD metadata + line index cached under `data/cod_index/`
  (`make cod-index` if absent)
- GSAS-II bootstraps itself into `.vendor/` on first use (network, one
  time; see `docs/installation.md`)