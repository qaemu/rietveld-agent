"""Spike 04 (a): build the M1 candidate fingerprint library.

Library = deterministic d-space fingerprints of known materials:
  - PbSO4 (Cu + Fe)   from the spike02 golden fixtures (GSAS-II simulated)
  - SiO2 quartz-family from the spike02 golden fixture (Cu)
  - NaCl halite       simulated now from COD 9011025 (Cu) via GSAS-II

Runtime analysis never needs GSAS-II: the library is pure JSON.

Run:  python spike_04_build_library.py
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIX_DIR = os.path.join(ROOT, "tests", "fixtures", "xrdml")
SPIKE4_DIR = os.path.join(ROOT, "data", "spike4")
IN_DIR = os.path.join(SPIKE4_DIR, "input")
WORK_DIR = os.path.join(SPIKE4_DIR, "work")
LIB_PATH = os.path.join(ROOT, "data", "candidates", "library.json")
VENDOR = os.path.join(ROOT, ".vendor", "GSAS-II")
CU_PRM = os.path.join(ROOT, "data", "spike", "input", "INST_XRY.PRM")

sys.path.insert(0, ROOT)

from core.ingest import InstrumentParams, PowderPattern, parse_xrdml, sample_fingerprint  # noqa: E402


def _sim_halite(work_dir: str) -> dict:
    """GSAS-II: simulate the halite CIF (Cu Kalpha, same protocol as spike02)."""
    sys.path.insert(0, VENDOR)
    sys.path.insert(0, os.path.join(VENDOR, "GSASII"))
    sys.path.insert(0, os.path.join(VENDOR, "bin"))
    os.environ.setdefault("GSASIIDATA", os.path.join(VENDOR, "GSASII"))
    from GSASII.GSASIIscriptable import G2Project

    cif = os.path.join(IN_DIR, "halite_9011025.cif")
    proj = G2Project(newgpx=os.path.join(work_dir, "halite_init.gpx"))
    proj.add_phase(cif, phasename="NaCl", fmthint="CIF")
    hist = proj.add_simulated_powder_histogram(
        "sim_nacl", CU_PRM, 15.0, 140.0, Tstep=0.02, scale=50000.0,
        phases=proj.phases())
    hist.set_refinements(
        {"Background": {"type": "chebyschev-1", "no. coeffs": 3, "refine": True}})
    proj.phase(0).set_HAP_refinements({"Scale": True})
    proj.data["Controls"]["data"]["max cyc"] = 3
    proj.refine()
    ycalc = hist.getdata("ycalc")
    x = hist.getdata("x")
    proj.save(os.path.join(work_dir, "halite.gpx"))
    imax = int(ycalc.argmax())
    return {"tth": x, "ycalc": ycalc, "top_peak_tth": float(x[imax]),
            "sha256": hashlib.sha256(open(cif, "rb").read()).hexdigest()}


def _entry(mid: str, name: str, phase_family: str, source_cif: str,
           provenance: str, simulated_with: dict, fp) -> dict:
    return {"id": mid, "name": name, "phase_family": phase_family,
            "source_cif": source_cif, "provenance": provenance,
            "simulated_with": simulated_with,
            "fingerprint": fp.to_dict()}


def main() -> None:
    os.makedirs(IN_DIR, exist_ok=True)
    os.makedirs(WORK_DIR, exist_ok=True)
    os.makedirs(os.path.dirname(LIB_PATH), exist_ok=True)

    cu_sim = {"anode": "CuKa", "wavelengths": [1.5405, 1.5443]}
    fe_sim = {"anode": "FeKa", "wavelengths": [1.9360, 1.9399]}

    entries = [
        _entry("mat-pbso4", "PbSO4", "PbSO4", "PbSO4-Wyckoff.cif",
               "APS GSAS-II tutorial", cu_sim,
               sample_fingerprint(parse_xrdml(os.path.join(FIX_DIR, "cu_PbSO4.xrdml")))),
        _entry("mat-pbso4-fe", "PbSO4", "PbSO4", "PbSO4-Wyckoff.cif",
               "APS GSAS-II tutorial", fe_sim,
               sample_fingerprint(parse_xrdml(os.path.join(FIX_DIR, "fe_PbSO4.xrdml")))),
        _entry("mat-sio2", "SiO2 (quartz-family)", "SiO2 (quartz-family)", "quartz_1009000.cif",
               "COD 1009000 (public domain)", cu_sim,
               sample_fingerprint(parse_xrdml(os.path.join(FIX_DIR, "cu_quartz.xrdml")))),
    ]
    print("[halite] simulating COD 9011025 (Cu Ka)...")
    sim = _sim_halite(WORK_DIR)
    hal = PowderPattern(sample_name="NaCl", source="sim_nacl",
                        tth=sim["tth"], intensity=sim["ycalc"],
                        instrument=InstrumentParams(
                            anode="CuKa", wavelengths=(1.5405, 1.5443),
                            tmin=15.0, tmax=140.0, step=0.02, npts=len(sim["tth"])))
    entries.append(_entry("mat-nacl", "NaCl (halite)", "NaCl (halite)", "halite_9011025.cif",
                          "COD 9011025 (public domain)", cu_sim,
                          sample_fingerprint(hal)))
    print(f"[halite] top peak @{sim['top_peak_tth']:.2f} deg (expect ~31.7 NaCl 200)")

    entries.sort(key=lambda e: e["id"])
    manifest = hashlib.sha256(
        json.dumps(entries, sort_keys=True).encode()).hexdigest()
    library = {"schema": "candidate-library/v0",
               "built_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
               "manifest_sha256": manifest,
               "materials": entries}
    with open(LIB_PATH, "w") as fh:
        json.dump(library, fh, indent=2)
    print(f"[library] {len(entries)} entries -> {LIB_PATH}")
    print(f"[manifest] {manifest[:16]}...")
    for e in entries:
        print(f"  - {e['id']}: {e['name']} ({e['simulated_with']['anode']}) "
              f"peaks={len(e['fingerprint']['peaks'])}")


if __name__ == "__main__":
    main()