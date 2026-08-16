# -*- coding: utf-8 -*-
"""Profile sweep for the QPA gate real-data samples (documented fallback).

Direct GSAS-II joint fits (pipeline Stage-D pattern: scales + 5-term
Chebyshev background only, frozen instrument profile) for the best CIF
combos from ``sweep_cifs.py``, at three profile settings:

  * default : instrument profile exactly as the calibrated PRM provides
  * narrow  : U=V=0, W~0.005 (FWHM ~0.13-0.14 deg, matches the data),
              Shift = 0

Output: ``data/qpa_gate/sweeps/<sample>_profile_sweep.json`` (phase wt
via Hill-Howard Mass*Scale, Rietveld wR).

Usage:
    .venv/bin/python benchmarks/qpa_gate/sweep_profile.py
"""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402

import benchmarks.qpa_gate.qpa_gate as g  # noqa: E402

CIFS = Path("data/unit12/work/cifs")
OUT = Path("data/qpa_gate/sweeps")

COMBOS = {
    "qarr_1f": [("2300112", "zincite"), ("1000059", "corundum"),
                ("1000043", "fluorite")],
    "iron_30_70": [("9002321", "magnetite"), ("5910082", "hematite")],
}
W_NARROW = 0.005

BAD_LST = ("SVD    singular", "on the verge of being singular",
           "convergence to a d  minimum")


def build_proj(sid, combo, w_val, work: Path, tag: str, lo: float, hi: float,
               a1: float, a2: float, ratio: float, pat):
    from benchmarks.eval.sim import ensure_gsasii
    ensure_gsasii(str(Path('.').resolve()),
                  str(Path('.').resolve() / '.vendor' / 'GSAS-II'), "")
    from GSASII.GSASIIscriptable import G2Project

    tth, y = pat.tth, pat.intensity
    m = (tth >= lo) & (tth <= hi)
    sig = np.sqrt(np.maximum(y[m], 1.0))
    xye = work / f"{tag}.xye"
    np.savetxt(xye, np.column_stack([tth[m], y[m], sig]),
               fmt="%.5f %.3f %.3f")
    prm = g._clone_prm(tag, work, tag, a1, a2, ratio)
    proj = G2Project(newgpx=str(work / f"{tag}.gpx"))
    names = []
    canon_of = {}
    for cod, canon in combo:
        nm = f"p{cod}"
        proj.add_phase(str(CIFS / f"{cod}.cif"), phasename=nm, fmthint="CIF")
        names.append(nm)
        canon_of[nm] = canon
    h = proj.add_powder_histogram(str(xye), iparams=str(prm), phases=names)
    h.set_refinements({"Background": {"type": "chebyschev-1",
                                      "no. coeffs": 5, "refine": True}})
    h.set_refinements({"Sample Parameters": ["Shift"]})
    if w_val is not None:
        ip = h.data["Instrument Parameters"][0]
        for k, v in (("U", 0.0), ("V", 0.0), ("W", w_val), ("X", 0.0),
                     ("Y", 0.0)):
            if k in ip and isinstance(ip[k], list) and len(ip[k]) >= 2:
                ip[k][0] = v
        h.data["Sample Parameters"]["Shift"] = [0.0, 0.0005, False]
    proj.data["Controls"]["data"]["max cyc"] = 40
    for i in range(len(names)):
        proj.phase(i).set_HAP_refinements({"Scale": True})
    proj.refine(makeBack=False)
    lstp = work / f"{tag}.lst"
    txt = lstp.read_text(errors="ignore") if lstp.exists() else ""
    wr = (re.findall(r"Final refinement wR =\s*([\d.]+)", txt)
          or [None])[-1]
    wr = float(wr) if wr is not None else None
    conv = bool(("Refinement successful" in txt)
                or ("Final refinement" in txt))
    bad = any(mm in txt for mm in BAD_LST)
    per_phase = []
    for nm in names:
        pd = proj.data["Phases"][nm]
        mass = float(pd["General"]["Mass"])
        try:
            scale = float(pd["Histograms"][h.name]["Scale"][0])
        except (KeyError, TypeError):
            scale = 1.0
        per_phase.append({"cod": nm[1:], "canon": canon_of[nm],
                          "scale": round(scale, 6), "mass": mass})
    smv = {q["cod"]: q["scale"] * q["mass"] for q in per_phase}
    tot = sum(smv.values())
    for q in per_phase:
        q["wt"] = round(100.0 * smv[q["cod"]] / tot, 2) if tot else 0.0
    return {"wR": wr, "converged": conv, "bad": bad,
            "phases": [{"cod": q["cod"], "canon": q["canon"],
                        "wt": q["wt"]} for q in per_phase]}


def main() -> None:
    g.build_manifest()
    OUT.mkdir(parents=True, exist_ok=True)
    for sid, combo in COMBOS.items():
        entry = [m for m in g.MANIFEST if m[0] == sid][0]
        _, fpath, instr, truth, (lo, hi) = entry
        a1, a2, ratio = g._wl(instr)
        pat = g.load_any(fpath, instr)
        rows = []
        for wname, w_val in (("default", None), ("narrow_w005", W_NARROW)):
            tag = f"{sid}_prof_{wname}"
            work = OUT / sid
            work.mkdir(parents=True, exist_ok=True)
            t0 = time.time()
            r = build_proj(sid, combo, w_val, work, tag, lo, hi, a1, a2,
                           ratio, pat)
            r["profile"] = wname
            rows.append(r)
            print(sid, wname, "->", r.get("phases"), "wR", r.get("wR"),
                  "conv", r.get("converged"), "in",
                  round(time.time() - t0), "s", flush=True)
        fp = OUT / f"{sid}_profile_sweep.json"
        fp.write_text(json.dumps(rows, indent=1) + "\n")
        print(f"[profile] {sid} -> {fp}", flush=True)


if __name__ == "__main__":
    main()