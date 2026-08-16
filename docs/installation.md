# Installation

How to install and run the `rietveld-agent` engine — on any platform,
from any of the three supported agent runtimes.

## Requirements

| requirement | version | notes |
|---|---|---|
| Python | 3.10+ | only hard dependency for analyses |
| git | any | to clone the repository |
| tectonic | optional | only to rebuild the manuscript PDF (`make paper`) |

GSAS-II (Apache-2.0) is **never installed by hand**: it is pinned in
`.vendor/GSAS-II` and bootstrapped by `benchmarks/eval/sim.py:ensure_gsasii`
on first use — the first refinement clones the official GSAS-II
repository into `.vendor/` (one-time network access, then offline).

COD screening: the repository ships a **candidate-restricted** COD
metadata export (`data/candidates/cod_entries.csv`, ~3.6k entries, CC0)
so the whole gate pipeline runs out of the box. The full-COD metadata
export (~528k entries, 368 MB) is optional: drop
`data/cod_index/cod_metadata_full.csv` in place (see `make cod-tree`,
`make cod-index`) to upgrade screening from candidate-restricted to
whole-database.

## 1. Clone

```bash
git clone https://github.com/qaemu/rietveld-agent
cd rietveld-agent
```

## 2. Choose your runtime

One of the three agent runtimes below. All are terminal-native,
local-first CLIs; the engine itself is runtime-agnostic (see
[Paper P-3](papers/paper3_deployment.md)).

### OpenCode

| platform | command |
|---|---|
| macOS | `brew install opencode` — or — `npm i -g opencode-ai` |
| Linux | `curl -fsSL https://opencode.ai/install | bash` — or — `npm i -g opencode-ai` |
| Windows | `npm i -g opencode-ai` (native); or use WSL and the Linux route |

Verify: `opencode --version` · Start: `opencode` (in the repo root)

### Claude Code (Anthropic)

| platform | command |
|---|---|
| macOS / Linux | `npm install -g @anthropic-ai/claude-code` (or the native installer script) |
| Windows | native installer; or use WSL and the Linux route |

Authenticate: run `claude` once and log in, or export
`ANTHROPIC_API_KEY`. · Verify: `claude --version` · Start: `claude`

### Codex (OpenAI)

| platform | command |
|---|---|
| macOS | `brew install codex` — or — `npm i -g @openai/codex` |
| Linux | `curl -fsSL https://codex.openai.com/install.sh | bash` — or — `npm i -g @openai/codex` |
| Windows | use WSL and the Linux route (codex is a terminal-native CLI) |

Authenticate: `codex login` (or export `OPENAI_API_KEY`). · Verify:
`codex --version` · Start: `codex`

All three runtimes read `AGENTS.md` at the repository root for
orientation on their first run.

## 3. Bootstrap the environment

```bash
make env      # creates .venv with numpy/scipy/matplotlib
make check    # syntax + structure checks (no network)
```

## 4. Run and verify

```bash
make report    # rerun every RQPA refinement (GSAS-II auto-bootstraps, ~45 min)
make figures   # regenerate the manuscript figures
make paper     # rebuild all paper PDFs (needs tectonic; optional)
make test      # regression suite (pytest)
make cite      # show the recommended citation
```

Validation evidence (reproducibility + known-answer recovery + gate
scoring, see [Paper P-2](papers/paper2_validation.md)):

```bash
python3 benchmarks/protocols/validate.py               # everything (~32 min)
python3 benchmarks/protocols/validate.py --skip-rerun  # KAT + gates only
python3 benchmarks/protocols/validate.py --skip-synth  # rerun + gates only
python3 benchmarks/protocols/validate.py --skip-rerun --sample <name>.xrdml  # one sample
```

## 5. Where the results land

| artifact | path |
|---|---|
| protocol refinements (md5-locked) | `data/unit15/results/unit15_report.json` |
| validation evidence (md5-locked) | `data/unit16/results/unit16_report.{json,md}` |
| structure catalogue (md5-recorded) | `data/structures/catalog.json` |
| manuscript + paper series (PDF) | `paper/main.pdf`, `paper/paper{1,2,3}_*.pdf` |
| research log | [`notes/`](../notes/) |

## Troubleshooting

- **`ensure_gsasii` needs network on first `make report`** — expected:
  it downloads the official GSAS-II installer once into `.vendor/`.
- **`tectonic: command not found`** — only `make paper` needs it
  (`brew install tectonic`); analyses and validation do not.
- **Windows without WSL** — OpenCode and Claude Code install natively
  via npm; use a POSIX shell (git-bash / PowerShell + npm) for `make`
  targets, or use WSL for the full experience.
- **Reproducibility mismatch after a code change** — the canonical
  content hash is computed over the scientific payload only; wall-clock
  timing is excluded by design (`elapsed_s`). A changed hash means a
  changed payload: check what you changed before claiming
  reproducibility.