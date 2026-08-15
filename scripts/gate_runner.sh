#!/usr/bin/env bash
# Run the full-COD screening + QPA gate (spike 20) with the repository
# virtualenv and stream a summary. Usage:
#   scripts/gate_runner.sh                 # whole 20-sample manifest
#   scripts/gate_runner.sh --only qarr_1a  # single sample
set -euo pipefail
cd "$(dirname "$0")/.."

PY="${PY:-.venv/bin/python}"
if [[ ! -x "$PY" ]]; then
  echo "virtualenv missing — run: make env" >&2
  exit 1
fi

LOG="data/spike20/gate.log"
mkdir -p data/spike20
echo "== spike-20 gate: $(date +%FT%T) ==" | tee -a "$LOG"
"$PY" benchmarks/spikes/spike_20_fullcod_qpa.py "$@" 2>&1 | tee -a "$LOG"
status=${PIPESTATUS[0]}
echo "== gate exit: $status ==" | tee -a "$LOG"
exit "$status"