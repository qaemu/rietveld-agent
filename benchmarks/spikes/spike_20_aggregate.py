#!/usr/bin/env python3
"""Aggregate the spike-20 gate results: re-verify every sample JSON's
phases against the CURRENT gate() (the JSONs may predate later gate
policy, e.g. the iron-oxide isomorphism-class budgets).

Usage: python3 benchmarks/spikes/spike_20_aggregate.py [results_dir]
Exit 0 iff ALL samples pass."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import spike_20_fullcod_qpa as S20

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else
           "data/spike20/results")


def main() -> int:
    S20.build_manifest()
    truths = {m[0]: m[3] for m in S20.MANIFEST}
    jsons = sorted(OUT.glob("*.json"))
    if not jsons:
        print("no result jsons found")
        return 2
    rows = []
    for p in jsons:
        d = json.loads(p.read_text())
        sid = p.stem
        v = S20.gate(sid, truths.get(sid, {}), d.get("phases", []))
        rows.append((sid, bool(v["ok"]), d.get("rwp"),
                     v["gaps"], d.get("time_s"), v["refined"]))
    n_pass = sum(1 for _, g, *_ in rows if g)
    n_old_pass = sum(1 for p in jsons
                     if json.loads(p.read_text()).get("gate"))
    print(f"{len(rows)} samples; {n_old_pass} json gate=true; "
          f"{n_pass} current-gate PASS")
    for sid, g, rwp, gaps, ts, refined in sorted(rows):
        tag = "PASS" if g else "FAIL"
        extra = "" if g else f" gaps={gaps}"
        print(f"  {sid:44s} {tag} wR={rwp} ({ts:.0f}s)"
              f"{extra} refined={ {k: round(v,1) for k, v in refined.items()} }")
    return 0 if n_pass == len(rows) and n_pass > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())