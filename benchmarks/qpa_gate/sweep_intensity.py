import sys, re, shutil, os
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))
CIF = ROOT / "data/spike12/work/cifs"

from benchmarks.qpa_gate.qpa_gate import _clone_prm
from benchmarks.eval.sim import ensure_gsasii
ensure_gsasii(str(ROOT), str(ROOT / '.vendor' / 'GSAS-II'), "")
from GSASII.GSASIIscriptable import G2Project

M_Q = {"z": 81.38, "f": 78.07, "c": 101.96}
M_I = {"h": 159.69, "m": 231.54}


def run(tag, sample, pid_map, mass_map, uiso=True, po_axis=None, maxcyc=80):
    WORK = ROOT / "data/qpa_gate" / "sweep2" / tag
    shutil.rmtree(WORK, ignore_errors=True)
    WORK.mkdir(parents=True)
    xye_src = ROOT / f"data/qpa_gate/work/{sample}/{sample}.xye"
    d = np.genfromtxt(xye_src, usecols=(0, 1, 2))
    xye = WORK / f"{tag}.xye"
    np.savetxt(xye, d, fmt="%.5f %.3f %.3f")
    prm = _clone_prm(tag, WORK, tag, 1.54056, 1.54439, 0.5)
    proj = G2Project(newgpx=str(WORK / f"{tag}.gpx"))
    names = []
    for pid in pid_map:
        proj.add_phase(str(CIF / f"{pid}.cif"), phasename=f"p{pid}", fmthint="CIF")
        names.append(f"p{pid}")
    h = proj.add_powder_histogram(str(xye), iparams=str(prm), phases=names)
    ip = h.data["Instrument Parameters"][0]
    ip["U"][0], ip["V"][0], ip["W"][0], ip["X"][0], ip["Y"][0] = (0.0, 0.0, 0.005, 0.0, 0.0)
    h.set_refinements({"Limits": {"low": 5.0, "high": 90.0},
                       "Background": {"type": "chebyschev-1", "no. coeffs": 5, "refine": True}})
    for i, n in enumerate(names):
        proj.phase(i).set_HAP_refinements({"Scale": True})
        if uiso:
            proj.phase(i).set_refinements({"Atoms": {"all": "U"}})
        if po_axis:
            po = proj.data["Phases"][n]["Histograms"][h.name]["Pref.Ori."]
            po[0], po[1], po[2], po[3] = "MD", 1.0, True, list(po_axis)
    proj.data["Controls"]["data"]["max cyc"] = maxcyc
    proj.refine(makeBack=False)
    txt = (WORK / f"{tag}.lst").read_text(errors="ignore")
    wr = (re.findall(r"Final refinement wR =\s*([\d.]+)", txt) or ["nan"])[-1]
    scales = {n: proj.data["Phases"][n]["Histograms"][h.name]["Scale"][0] for n in names}
    tot = sum(scales[n] * mass_map[pid_map[n.split("p")[1]]] for n in names)
    wt = {pid_map[n.split("p")[1]]: round(100 * scales[n] * mass_map[pid_map[n.split("p")[1]]] / tot, 1)
          for n in names}
    extras = ""
    if po_axis:
        zn = names[0]
        g = proj.data["Phases"][zn]["Histograms"][h.name].get("Pref.Ori.", [None, None])[1]
        extras = f" G={g}"
    print(f"[{tag}] wR={float(wr):.2f} wt={wt}{extras}", flush=True)


if __name__ == "__main__" and os.environ.get("SWEEP3") != "1":
    print("== sweep2 iron_30_70: ambient-magnetite CIFs x hematite, Uiso refined ==", flush=True)
    for mid in ("2300616", "2101535", "1539747", "9002316", "1513304", "1010369"):
        for hid in ("5910082", "9000139"):
            run(f"ir_{mid}_{hid}", "iron_30_70", {hid: "h", mid: "m"}, M_I)
    print("== sweep2 qarr_1f: zincite CIFs, Uiso refined ==", flush=True)
    for zid in ("2300112", "2300450", "9004178"):
        run(f"q_{zid}", "qarr_1f", {zid: "z", "1000043": "f", "1000059": "c"}, M_Q)
    print("== sweep2 qarr_1f: best zincite + March-Dollase PO axis 001 ==", flush=True)
    run("q_po_001", "qarr_1f", {"2300450": "z", "1000043": "f", "1000059": "c"}, M_Q,
        po_axis=(0, 0, 1))
    print("done", flush=True)


if os.environ.get("SWEEP3") == "1":
    print("== sweep3 iron: ambient magic mags x hems, NO Uiso ==", flush=True)
    for mid in ("2300616", "2101535", "1539747", "9002316", "1513304", "1010369"):
        for hid in ("5910082", "9000139", "9015065"):
            run(f"s3_{mid}_{hid}", "iron_30_70", {hid: "h", mid: "m"}, M_I, uiso=False)
    print("== sweep3 qarr: zincite x PO (MD) via proper Pref.Ori. key, no Uiso ==", flush=True)
    for zid in ("2300112", "2300450", "9004178"):
        run(f"s3q_{zid}_po001", "qarr_1f", {zid: "z", "1000043": "f", "1000059": "c"}, M_Q,
            uiso=False, po_axis=(0, 0, 1))
        run(f"s3q_{zid}_po100", "qarr_1f", {zid: "z", "1000043": "f", "1000059": "c"}, M_Q,
            uiso=False, po_axis=(1, 0, 0))
    print("done", flush=True)
