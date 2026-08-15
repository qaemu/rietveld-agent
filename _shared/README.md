# _shared

Reference-style home for resources shared across the repository's
skills and experiments (mirrors the `_shared/` layout of
Light0305/Light-skills).

## Contents

| path | shared resource |
|---|---|
| `skills/contracts/` | agent-facing decision criteria, reporting rules, parameter allowlists (canonical home; linked from here) |
| `core/codindex.py` | COD line-index units (`DQ`/`D_UNIT` conventions) shared by `codsearch`, the spikes and the CLI |
| `core/codsearch.py` | complete-COD screening primitives used by the QPA gate, `cli analyze --full-cod` and the run bundle |

## Conventions

- **Determinism**: anything "shared" must be hash-locked or
  regenerable; do not edit `data/structures/catalog.json` or
  governance payloads in place (see `AGENTS.md` invariants).
- **Reuse, don't fork**: new experiments should import from `core/`
  (and from `benchmarks/eval/` for the synthetic pattern model) rather
  than copy behavior.