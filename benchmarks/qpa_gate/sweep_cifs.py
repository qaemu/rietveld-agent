# -*- coding: utf-8 -*-
"""CIF-combo sweep driver for the QPA gate real-data samples.

For every combination, the *unmodified* pipeline engine (``gsas_qpa``:
Stage A model selection -> B forward selection -> D Hill-Howard wt%) runs
with a FORCED phase list (the combos below).  Output is written as JSON
lines including per-phase wt and Rietveld wR, so the sweep is a quantified,
reproducible record of "which CIF choice drives the gate outcome" --
the configured fallback documentation artifact for the gate.

Usage:
    .venv/bin/python benchmarks/qpa_gate/sweep_cifs.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import benchmarks.qpa_gate.qpa_gate as g  # noqa: E402

CIFS = Path("data/unit12/work/cifs")
OUT = Path("data/qpa_gate/sweeps")

# sample -> ordered list of forced CIF combos (cod, canon)
COMBOS = {
    "qarr_1f": [
        [("2300112", "zincite"), ("1000059", "corundum"), ("1000043", "fluorite")],
        [("2300450", "zincite"), ("1000059", "corundum"), ("1000043", "fluorite")],
        [("9004178", "zincite"), ("1000059", "corundum"), ("1000043", "fluorite")],
        [("9004179", "zincite"), ("1000059", "corundum"), ("1000043", "fluorite")],
        [("1577381", "zincite"), ("1000059", "corundum"), ("1000043", "fluorite")],
        [("2107059", "zincite"), ("1000059", "corundum"), ("1000043", "fluorite")],
        [("2300112", "zincite"), ("1000017", "corundum"), ("1000043", "fluorite")],
        [("2300112", "zincite"), ("1000059", "corundum"), ("2300449", "fluorite")],
    ],
    "iron_30_70": [
        [("2300616", "magnetite"), ("5910082", "hematite")],
        [("1011032", "magnetite"), ("5910082", "hematite")],
        [("2300616", "magnetite"), ("9000139", "hematite")],
        [("1011032", "magnetite"), ("9000139", "hematite")],
        [("9002320", "magnetite"), ("5910082", "hematite")],
        [("9002321", "magnetite"), ("5910082", "hematite")],
    ],
}


def run_sweep(sid: str) -> None:
    g.build_manifest()
    entry = [m for m in g.MANIFEST if m[0] == sid][0]
    _, fpath, instr, truth, (lo, hi) = entry
    pat = g.load_any(fpath, instr)
    if pat is None:
        raise SystemExit(f"{sid}: pattern load failed")
    a1, a2, ratio = g._wl(instr)
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for combo in COMBOS[sid]:
        phases = [
            {"name": f"p{cod}", "cif": CIFS / f"{cod}.cif",
             "canon": canon}
            for cod, canon in combo
        ]
        tag = "-".join(cod for cod, _ in combo)
        t0 = time.time()
        try:
            q = g.gsas_qpa(pat, phases, OUT / sid, f"{sid}_{tag}", lo, hi,
                           a1, a2, ratio, sync=(a2 == 0.0))
            row = {
                "combo": [cod for cod, _ in combo],
                "phases": [{"cod": w["cod"], "canon": w["canon"],
                            "wt": w["wt"]} for w in q["wt"]],
                "wR": round(q["wR"], 3),
                "converged": q.get("converged"),
                "time_s": round(time.time() - t0, 1),
            }
        except Exception as exc:  # noqa: BLE001 -- sweep is a record
            row = {"combo": [cod for cod, _ in combo],
                   "error": f"{type(exc).__name__}: {exc}"}
        rows.append(row)
        print(sid, row["combo"], "->", row.get("phases"),
              "wR", row.get("wR"), flush=True)
    fp = OUT / f"{sid}_cif_sweep.json"
    fp.write_text(json.dumps(rows, indent=1) + "\n")
    print(f"[sweep] {sid}: {len(rows)} combos -> {fp}", flush=True)


if __name__ == "__main__":
    only = sys.argv[1:] or list(COMBOS)
    for s in only:
        run_sweep(s)