#!/usr/bin/env python3
"""Audit the shared CIF cache after a (possibly parallel) gate run.

Two failure classes are distinguished:
  - torn:   file too small or lacking a `data_` header -> a truncated
            download (pre-atomic-write races could produce these)
  - broken: complete COD file that gemmi/RQPA cannot read (e.g. 1001390)

Usage:
  python3 benchmarks/qpa_gate/check_cifs.py [--recover]

With --recover, torn files are deleted so the next gate run re-downloads
them atomically.  Broken entries are never touched (they are COD-side
data; re-downloading yields the same file).

Exit 0 if the cache is clean, 1 otherwise (report only)."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "benchmarks" / "protocols"))

from core.codindex import cif_calc_lines          # noqa: E402

CIF_DIR = ROOT / "data" / "unit12" / "work" / "cifs"


def main() -> int:
    recover = "--recover" in sys.argv[1:]
    if not CIF_DIR.exists():
        print(f"no CIF cache at {CIF_DIR}")
        return 0
    cifs = sorted(CIF_DIR.glob("*.cif"))
    torn, broken = [], []
    for p in cifs:
        try:
            txt = p.read_text(errors="ignore")
        except Exception:
            torn.append((p, "unreadable file"))
            continue
        if len(txt) < 100 or "data_" not in txt[:400]:
            torn.append((p, f"short/headerless ({len(txt)} bytes)"))
            continue
        try:
            cif_calc_lines(str(p), dmin=1.0, dmax=22.0)
        except Exception as e:
            broken.append((p, str(e)[:80]))
    print(f"{len(cifs)} cached CIFs: {len(torn)} torn, {len(broken)} broken")
    for p, why in torn:
        print(f"  TORN    {p.name}: {why}")
        if recover:
            p.unlink(missing_ok=True)
            print(f"          -> deleted (will re-download atomically)")
    for p, why in broken:
        print(f"  BROKEN  {p.name}: {why}")
    if recover and torn:
        print("torn files deleted; re-run the gate to re-fetch them")
    return 1 if (torn or broken) and not recover else 0


if __name__ == "__main__":
    raise SystemExit(main())