# Deploying rietveld-agent on OpenCode, Claude Code, and Codex

**Installation, verification, and governance** — Paper P-3 of the
docs/papers series.

**Thesis.** The scientific engine of `rietveld-agent` is
runtime-agnostic: the same deterministic pipeline (Paper P-1 [1]) and the
same validation evidence (Paper P-2 [2]) are available from any
agent-aware terminal — OpenCode, Claude Code, or Codex — on macOS, Linux,
or Windows (native or WSL), with a three-command installation, scripted
dependency bootstrapping (GSAS-II pinned, never bundled), and a checkable
verification path.

**Status.** Working paper, v0.2.0, August 2026. The complete installation
reference is [`docs/installation.md`](../installation.md). Constructed
following the fifteen-step writing framework of Drake & Han (2025)
[doi:10.1371/journal.pcbi.1013505](https://doi.org/10.1371/journal.pcbi.1013505).

---

## 1. Introduction

A repository is only as portable as its installation story. Scientific
repositories have an additional requirement: whoever runs them — human or
agent — must be able to verify that what they got is what was measured.

**Gap.** Agent runtimes (OpenCode, Claude Code, Codex) are
terminal-native, but most scientific repositories assume interactive GUI
sessions, ad-hoc environments, or cloud execution.

**Response.** The `rietveld-agent` repository [3] is built for
non-interactive execution from the start: every result is rebuilt by
`make` targets, GSAS-II [4] is pinned in `.vendor/` and bootstrapped
deterministically, and an `AGENTS.md` orientation file points any agent
at the contracts it must respect. This paper documents installation,
verification, and governance for the three runtimes.

## 2. Installation matrix

Requirements in all cases: **Python 3.10+** and **git**. GSAS-II
(Apache-2.0) is pinned in `.vendor/GSAS-II` and installed by
`benchmarks/eval/sim.py:ensure_gsasii` on first use; it is never bundled
as a dependency. `tectonic` is required only for building the manuscript
PDF (`make paper`), not for running analyses.

### OpenCode

| platform | install |
|---|---|
| macOS | `brew install opencode` or `npm i -g opencode-ai` |
| Linux | `curl -fsSL https://opencode.ai/install | bash` or `npm i -g opencode-ai` |
| Windows | `npm i -g opencode-ai` (native), or WSL + the Linux route |

Verify with `opencode --version`, then start in the repository root:
`cd rietveld-agent && opencode`.

### Claude Code

| platform | install |
|---|---|
| macOS / Linux | `npm install -g @anthropic-ai/claude-code` (or the native installer script) |
| Windows | native installer, or WSL + the Linux route |

Authenticate with `claude` (browser login) or an `ANTHROPIC_API_KEY`;
verify with `claude --version`; start with `cd rietveld-agent && claude`.

### Codex (OpenAI)

| platform | install |
|---|---|
| macOS | `brew install codex` or `npm i -g @openai/codex` |
| Linux | `curl -fsSL https://codex.openai.com/install.sh | bash` or `npm i -g @openai/codex` |
| Windows | WSL + the Linux route (codex is a terminal-native CLI) |

Authenticate with `codex login` (or `OPENAI_API_KEY`); verify with
`codex --version`; start with `cd rietveld-agent && codex`.

## 3. First run and verification

```bash
make env       # virtualenv with numpy/scipy/matplotlib
make check     # syntax + structure checks (fast, no network)
make report    # rerun every RQPA refinement (GSAS-II, ~45 min)
make paper     # rebuild the manuscript PDF (tectonic; optional)
```

The full validation evidence of Paper P-2 is reproduced with

```bash
python3 benchmarks/spikes/spike_16_validate.py               # everything (~32 min)
python3 benchmarks/spikes/spike_16_validate.py --skip-rerun  # KAT + gates only
python3 benchmarks/spikes/spike_16_validate.py --skip-synth  # rerun + gates only
```

Any agent can therefore verify the engine from scratch without human
intervention, and the canonical content hash of the report
(`f43d8be2…`) makes the verification objective.

## 4. Governance for agent-driven scientific work

The repository enforces five invariants that make agent operation safe:

1. **Controlled scientific inputs** — catalog, calibrations, policies,
   models — are hashed, versioned and reviewable; no update may silently
   change the meaning of an existing analysis.
2. **Local first** — no telemetry, no cloud refinement, Apache-2.0.
3. **Weak context only** — the sample name may add candidate phases,
   never remove them or improve scores.
4. **Defensible uncertainty** — per-phase statuses (`supported`,
   `inconclusive`, `not_selected`, `out_of_domain`, `held`, `failed`)
   instead of unearned accuracy claims.
5. **Interchangeable hosts** — the three runtimes are thin front-ends
   around the same engine, so results do not depend on which agent drove
   them.

`AGENTS.md` at the repository root is the single orientation point for
any runtime.

## 5. Closing

Installation, verification, and governance together make
`rietveld-agent` usable by humans and agents alike on all major
platforms, with objective reproducibility as the anchor. The only
external requirement is Python 3.10+; everything else is pinned or
scripted. See Paper P-1 [1] for the protocol and Paper P-2 [2] for the
validation evidence.

## References

1. qaemu (2026). *Reproducible Rietveld quantitative phase analysis under a bounded budget* — Paper P-1 of this series. [paper1_protocol.md](paper1_protocol.md)
2. qaemu (2026). *Known-answer validation of an automated Rietveld QPA engine* — Paper P-2 of this series. [paper2_validation.md](paper2_validation.md)
3. qaemu (2026). rietveld-agent: a deterministic Rietveld QPA engine for scientific agent runtimes. https://github.com/qaemu/rietveld-agent-spikes (Apache-2.0)
4. Toby, B. H., & Von Dreele, R. B. (2013). GSAS-II: the genesis of a modern open-source all purpose crystallography software package. *J. Appl. Cryst.*, 46, 544–549. doi:10.1107/S0021889813003531
5. Drake, J. M., & Han, B. A. (2025). How to write a scientific paper in fifteen steps. *PLoS Comput. Biol.*, 21(9), e1013505. doi:10.1371/journal.pcbi.1013505 (PMC12459795)

---

*License: Apache-2.0 for original code; GSAS-II under its own Apache-2.0
license; catalog and structure files CC0 from COD (attribution
preserved).*