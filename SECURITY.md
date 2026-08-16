# Security Policy

## Scope

This repository is research software (Apache-2.0) with no network
service of its own. The trust boundary worth caring about:

- **Scientific integrity** — the engine is deterministic and
  hash-locked; a vulnerability that lets inputs alter reported phases or
  quantities without changing the hashes would break the research
  contract. The invariants in `AGENTS.md` and the governance schemas in
  `governance/` are part of the security model.
- **Supply chain** — dependencies are pinned in `requirements.txt` and
  GSAS-II is pinned under `.vendor/` (bootstrapped from the official
  installer); do not relax pins without review.
- **Prompt/context injection** — the engine is designed to be driven by
  agent runtimes that read project files and logs. Content from
  untrusted documents must never be able to relax policies or alter
  verdict gates (weak-context rule, `AGENTS.md`).

## Reporting a vulnerability

Do not open a public issue for security problems. Contact the
maintainer privately through the GitHub repository owner profile
(`qaemu`, see `CITATION.cff`), or open a draft advisory / private
vulnerability report via GitHub's security tab.

Please include: repository commit, affected inputs, expected vs actual
behavior, and a minimal reproduction. Acknowledgment within 7 days;
mitigation and disclosure coordinated with you before public release.

## Supported status

| component | supported |
|---|---|
| `core/` engine, `cli/` operator CLI | supported |
| `benchmarks/` units and harness | maintained; `make check` must pass |
| `paper/`, `docs/`, `notes/` | documentation (non-executable content has no security surface) |
| `.vendor/GSAS-II` | pinned external tool; update only with review |

## Data handling

Reference patterns and structure files are either shipped with the
repository or fetched from the zenodo-pinned dataset; COD entries are
queried read-only. The repository performs no telemetry and makes no
outbound calls with user data. If a contribution introduces an outbound
request, it must be optional, documented, and reviewed.