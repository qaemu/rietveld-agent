# -*- coding: utf-8 -*-
"""qarr_1f preferred-orientation (March-Dollase) experiment.

Data-zincite high-angle lines (110/103/112/200) are ~2-3.4x stronger than
any COD CIF models; every CIF/profile combination under-assigns zincite by
21-28 wt points.  Zincite is uniaxial (wurtzite, P63mc): a c-axis March-
Dollase texture with ratio r < 1 boosts reflections perpendicular to c
(110 is equatorial) -- a physically plausible explanation for the excess.

Protocol: direct GSAS-II joint fit (pipeline Stage-D pattern) of
qarr_1f with CIFs 2300112 (zincite)/1000059 (corundum)/1000043 (fluorite),
scales + 5-term Chebyshev background refined, with March-Dollase on the
zincite phase: axis [0,0,1], ratio grid (fixed), then one run refining the
ratio.  Output: data/qpa_gate/sweeps/qarr_1f_po_sweep.json.

Usage:
    .venv/bin/python benchmarks/qpa_gate/exp_po.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402

import benchmarks.qpa_gate.qpa_gate as g  # noqa: E402
from benchmarks.eval.sim import ensure_gsasii  # noqa: E402

CIFS = Path("data/unit12/work/cifs")
OUT = Path("data/qpa_gate/sweeps")
COMBO = [("2300112", "zincite"), ("1000059", "corundum"),
         ("1000043", "fluorite")]
AXIS = [0, 0, 1]
RATIOS = [1.0, 0.8, 0.65, 0.5, 0.35, "refine"]

BAD_LST = ("SVD    singular", "on the verge of being singular",
           "convergence to a d  minimum")


def joint_fit(combo, ratio, tag: str, work: Path, lo: float, hi: float,
              a1: float, a2: float, ratio2: float, pat, zincite_cod: str):
    tth, y = pat.tth, pat.intensity
    m = (tth >= lo) & (tth <= hi)
    sig = np.sqrt(np.maximum(y[m], 1.0))
    xye = work / f"{tag}.xye"
    np.savetxt(xye, np.column_stack([tth[m], y[m], sig]),
               fmt="%.5f %.3f %.3f")
    prm = g._clone_prm(tag, work, tag, a1, a2, ratio2)
    from GSASII.GSASIIscriptable import G2Project

    proj = G2Project(newgpx=str(work / f"{tag}.gpx"))
    names, canon_of = [], {}
    for cod, canon in combo:
        nm = f"p{cod}"
        proj.add_phase(str(CIFS / f"{cod}.cif"), phasename=nm, fmthint="CIF")
        names.append(nm)
        canon_of[nm] = canon
    h = proj.add_powder_histogram(str(xye), iparams=str(prm), phases=names)
    h.set_refinements({"Background": {"type": "chebyschev-1",
                                      "no. coeffs": 5, "refine": True}})
    proj.data["Controls"]["data"]["max cyc"] = 60
    for i in range(len(names)):
        proj.phase(i).set_HAP_refinements({"Scale": True})

    # March-Dollase on the zincite phase only.  GSASIIstrIO consumes it from
    # the PHASE's per-histogram entry 'Phases[phase]['Histograms'][h]
    # ['Pref.Ori.'] = ['MD', ratio, refineFlag, axis, SH order, ...].
    zn = f"p{zincite_cod}"
    for i in range(len(names)):
        ph = proj.phase(i)
        if ph.name == zn:
            ref = (ratio == "refine")
            r = 1.0 if ref else float(ratio)
            # phase-level (GUI/bookkeeping) + phase-per-histogram (engine)
            ph.data["Pref.Ori."] = ["MD", r, ref, AXIS, 0, {}, [], 0.1]
            ph.data["Histograms"][h.name]["Pref.Ori."] = [
                "MD", r, ref, AXIS, 0, {}, [], 0.1]

    proj.refine(makeBack=False)
    lstp = work / f"{tag}.lst"
    txt = lstp.read_text(errors="ignore") if lstp.exists() else ""
    wr = (re.findall(r"Final refinement wR =\s*([\d.]+)", txt)
          or [None])[-1]
    wr = float(wr) if wr is not None else None
    conv = bool(("Refinement successful" in txt)
                or ("Final refinement" in txt))
    bad = any(mm in txt for mm in BAD_LST)
    # report the refined MD ratio back (float or 'n/a')
    md = proj.data["Phases"][zn]["Pref.Ori."]
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
            "md_ratio": None if ratio != "refine" else round(float(md[1]), 4),
            "axis": list(md[3]),
            "phases": [{"cod": q["cod"], "canon": q["canon"],
                        "wt": q["wt"]} for q in per_phase]}


def main() -> None:
    ensure_gsasii(str(Path('.').resolve()),
                  str(Path('.').resolve() / '.vendor' / 'GSAS-II'), "")
    g.build_manifest()
    entry = [m for m in g.MANIFEST if m[0] == "qarr_1f"][0]
    _, fpath, instr, truth, (lo, hi) = entry
    a1, a2, ratio2 = g._wl(instr)
    pat = g.load_any(fpath, instr)
    OUT.mkdir(parents=True, exist_ok=True)
    work = OUT / "qarr_1f"
    rows = []
    for r in RATIOS:
        tag = f"qarr_1f_po_{'refine' if r == 'refine' else str(r).replace('.', 'p')}"
        try:
            row = joint_fit(COMBO, r, tag, work, lo, hi, a1, a2, ratio2,
                            pat, "2300112")
            row["md_setting"] = str(r)
            rows.append(row)
            print(row["md_setting"], "->",
                  {p["canon"]: p["wt"] for p in row["phases"]},
                  "wR", row["wR"], "md_refined", row["md_ratio"],
                  "conv", row["converged"], flush=True)
        except Exception as exc:  # noqa: BLE001 -- experiment record
            rows.append({"md_setting": str(r),
                         "error": f"{type(exc).__name__}: {exc}"})
            print("ERROR", r, exc, flush=True)
    fp = OUT / "qarr_1f_po_sweep.json"
    fp.write_text(json.dumps(rows, indent=1) + "\n")
    print(f"[po] {len(rows)} rows -> {fp}", flush=True)


if __name__ == "__main__":
    main()