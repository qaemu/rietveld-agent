import sys, re, shutil
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
ZINC = {"2300112": "z", "2300450": "z", "9004178": "z", "9004179": "z",
        "1577381": "z", "2107059": "z"}
FLU = {"1000043": "f"}
COR = {"1000059": "c"}


def run(tag, sample, pid_map, mass_map, prof, po_axis=None, maxcyc=60):
    WORK = ROOT / "data/qpa_gate" / "sweep" / tag
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
    ip["U"][0], ip["V"][0], ip["W"][0], ip["X"][0], ip["Y"][0] = prof
    h.set_refinements({"Limits": {"low": 5.0, "high": 90.0},
                       "Background": {"type": "chebyschev-1", "no. coeffs": 5, "refine": True}})
    for i, n in enumerate(names):
        proj.phase(i).set_HAP_refinements({"Scale": True})
        if po_axis and pid_map.get(n.split('p')[1]):
            # raw HAP PO entry: [h, k, l, G(March-Dollase), refine]
            proj.data["Phases"][n]["Histograms"][h.name]["PO"] = [*po_axis, 1.0, 1]
            proj.data["Phases"][n]["Histograms"][h.name]["POhkl"] = [*po_axis, 1.0, 1]
    proj.data["Controls"]["data"]["max cyc"] = maxcyc
    proj.refine(makeBack=False)
    txt = (WORK / f"{tag}.lst").read_text(errors="ignore")
    wr = (re.findall(r"Final refinement wR =\s*([\d.]+)", txt) or ["nan"])[-1]
    scales = {n: proj.data["Phases"][n]["Histograms"][h.name]["Scale"][0] for n in names}
    tot = sum(scales[n] * mass_map[pid_map[n.split('p')[1]]] for n in names)
    wt = {n.split('p')[1]: round(100 * scales[n] * mass_map[pid_map[n.split('p')[1]]] / tot, 1)
          for n in names}
    print(f"[{tag}] wR={float(wr):.2f} wt={wt}", flush=True)


print("== sweep qarr_1f: zincite-CIF (fluorite 1000043, corundum 1000059), W=0.005 ==")
for zid in ZINC:
    run(f"q1f_z{zid}", "qarr_1f", {zid: "z", **FLU, **COR}, M_Q, (0.0, 0.0, 0.005, 0.0, 0.0))
print("== sweep iron_30_70: magnetite x hematite, W=0.005 ==")
for mid in ("2300616", "1011032", "9002320"):
    for hid in ("5910082", "9000139", "9015065"):
        run(f"ir_{mid}_{hid}", "iron_30_70",
            {hid: "h", mid: "m"}, {"h": 159.69, "m": 231.54},
            (0.0, 0.0, 0.005, 0.0, 0.0))
print("== sweep qarr_1f PO (March-Dollase, axis 001) ==")
run("q1f_po_md", "qarr_1f", {"2300450": "z", **FLU, **COR}, M_Q,
    (0.0, 0.0, 0.005, 0.0, 0.0), po_axis=(0, 0, 1))
print("done", flush=True)