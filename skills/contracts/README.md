# Skills contracts (shared, versioned)

These documents are the **shared deterministic contracts** referenced by the
agent skills (`setup-rietveld`, `analyze-raw-pxrd`, `review-rietveld-run`).
The same versioned contracts apply identically to Codex, Claude Code, and
OpenCode. **Skills must not duplicate these documents** — they reference them.

| Document | What it defines |
|---|---|
| `parameter-allowlist.json` | Exact set of GSAS-II refinement parameters a plan may touch, with rationale and bounds. |
| `decision-criteria.md` | Status semantics and the strict `supported`-call criteria. |
| `reporting-rules.md` | Evidence-bundle structure and the measured/checks/AI separation. |
| `safety-constraints.md` | Non-negotiable constraints for skills and AI hosts. |

Versioning: all changes to these files must bump the `contracts_version`
below and be reviewed by a scientist before release. Nothing in a published
analysis may reference a contracts revision other than the pinned one.

Current revision: `contracts/0.1.0` (draft — Phase 0 spike).