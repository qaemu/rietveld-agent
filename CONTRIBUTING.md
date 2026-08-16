# Contributing

Thanks for considering a contribution to `rietveld-agent`. This is a
deterministic scientific repository: every number in the README is
reproduced by `make`, and the invariants in `AGENTS.md` are normative.
Please read them before opening an issue or a pull request.

## Ground rules

1. **Never silently modify scientific inputs.** `data/structures/catalog.json`,
   pinned calibrations, `governance/` policies and schemas are
   hash-locked and reviewed. Route changes through the documented
   governance workflow; do not edit them as a side effect of a feature.
2. **Keep results reproducible.** New features must keep
   `make check` and `make test` green. If a feature changes scientific
   payloads, the canonical content hash in `docs/` and the validation
   evidence must be regenerated and committed together with the change.
3. **Local first.** No telemetry, no cloud refinement, no calls that ship
   user data anywhere. Offline operation after `make env` is a design
   requirement.
4. **Weak context only.** Sample-name context may add candidate phases;
   it may never remove candidates, improve scores, or confirm a phase.

## Workflow

1. Open an issue describing the problem or the unit you intend to run.
2. Fork or branch; commit in small, self-contained steps (see
   [CHANGELOG.md](CHANGELOG.md) for the format).
3. Add or update a regression test under `tests/` and — for engine
   changes — a numbered unit entry with a log under `notes/`.
4. Run locally:
   ```bash
   make env
   make check
   make test
   ```
   For protocol-level changes also run `make report` (GSAS-II, ~45 min)
   and the validation harness
   (`python3 benchmarks/protocols/validate.py`).
5. Open a pull request; the diff must not touch hash-locked payloads
   unless the change is intentional and documented.

## Code style

- Python 3.10+, `from __future__ import annotations`, no external
  dependencies beyond those in `requirements.txt`.
- Keep the engine runtime-agnostic: `core/` must not import agent
  runtime SDKs; `cli/` and `benchmarks/` are the only places that may
  touch execution specifics.
- Type hints on public functions; one logical change per commit.

## Reporting issues

Include: the command you ran, the exact error, the repository commit
(`git rev-parse HEAD`), and whether you reproduced it after
`make env`. Security-sensitive findings: see [SECURITY.md](SECURITY.md).