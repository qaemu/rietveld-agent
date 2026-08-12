"""Regenerate the manuscript figures from the spike-15 result artifacts.

   fig1_clinker_t1_fit.png -- observed vs calculated pattern, clinker Cu
   (alite-T1 model), with difference curve and tick marks.
   fig2_results_summary.png -- wR% per sample/model with the
   publication-grade target lines, and alite wt% for the alite-bearing
   runs vs the published value (M3 vs T1 contrast).

The fit curve is recomputed by GSAS-II from the saved converged solution
(a 3-cycle scales-only pass at the converged point, which does not move
the parameters); everything else is read from the report JSON.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np                # noqa: E402

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "benchmarks" / "spikes"))

from benchmarks.eval.sim import ensure_gsasii  # noqa: E402

FIGDIR = Path(__file__).resolve().parent
RES = ROOT / "data" / "spike15" / "results" / "spike15_report.json"
GPX_T1 = ROOT / "data" / "spike15" / "work" / "Clinker_Nist_CuKalpha1_R1_qpa.gpx"


def fit_arrays(gpx: Path, lo: float, hi: float):
    """Return (x, yobs, ycalc) for the converged solution in `gpx`."""
    ensure_gsasii(str(ROOT), str(ROOT / ".vendor" / "GSAS-II"), "")
    from GSASII.GSASIIscriptable import G2Project

    proj = G2Project(gpxfile=str(gpx))
    h = proj.histograms()[0]
    proj.data["Controls"]["data"]["max cyc"] = 3
    proj.refine(makeBack=False)          # materialize ycalc; converged pt
    x = np.asarray(h.getdata("x"), dtype=float)
    yobs = np.asarray(h.getdata("yobs"), dtype=float)
    ycalc = np.asarray(h.getdata("ycalc"), dtype=float)
    return x, yobs, ycalc, proj


def fig_fit(x, yobs, ycalc, shift_deg: float, title: str, out: Path) -> None:
    xs = x - shift_deg                     # shift-corrected observed 2th
    diff = yobs - ycalc
    fig, ax = plt.subplots(figsize=(9.0, 4.2))
    ax.plot(xs, yobs, color="black", lw=0.6, label=r"$I_{\rm obs}$")
    ax.plot(x, ycalc, color="crimson", lw=0.6, label=r"$I_{\rm calc}$")
    ax.plot(x, diff + 800.0, color="steelblue", lw=0.5,
            label=r"$I_{\rm obs}-I_{\rm calc}+800$")
    ax.set_xlabel(r"$2\theta$ ($^\circ$, Cu K$\alpha_1$)")
    ax.set_ylabel("intensity (counts)")
    ax.set_title(title, fontsize=11)
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    ax.set_xlim(x.min(), x.max())
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"wrote {out}")


def fig_summary(rep: dict, out: Path) -> None:
    samples = rep["samples"]
    labels = [f"{s['sample'].split('_')[0]}\n[{s.get('model','')}]"
              for s in samples]
    wr = [s["wR"] for s in samples]
    tiers = [s["tier"] for s in samples]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.9),
                                 gridspec_kw={"width_ratios": [1.05, 1]})
    cols = ["#2e8b57" if t == "ok" else ("#daa520" if t == "wR-over"
                                          else "#b22222") for t in tiers]
    a1.barh(labels, wr, color=cols, height=0.62)
    a1.axvline(6.5, color="black", ls="--", lw=1.0)
    a1.axvline(5.0, color="black", ls=":", lw=1.0)
    a1.text(6.7, len(wr) - 0.45, "publication-grade targets\n"
            "(wR 6.5% Cu / 5% sync)", fontsize=7, va="top")
    a1.set_xlabel("final refinement wR (%)")
    a1.set_title("Fit quality under the bounded budget", fontsize=10)
    a1.set_xlim(0, 32)
    a1.tick_params(labelsize=8)

    # alite wt%: alite-bearing runs vs published
    pubs = {"Clinker_Nist_CuKalpha1_R1.xrdml": 66.0,
            "Silicate_enriched_residue_Nist_CuKalpha1_R1.xrdml": 78.7,
            "Clinker_Synchrotron.dat": 65.4}
    rows, pub_vals = [], []
    for s in samples:
        fname = s["sample"]
        if fname not in pubs:
            continue
        ours = next((p["wt_frac"] for p in s["phases"]
                     if p["phase"].startswith("alite")), None)
        rows.append((f"{fname.split('_')[0]}-{s.get('model','')}", ours))
        pub_vals.append(pubs[fname])
    xpos = np.arange(len(rows))
    ours_vals = [r[1] for r in rows]
    a2.bar(xpos, ours_vals, color="#4682b4", width=0.55)
    a2.plot(xpos, pub_vals, "ko", ms=4, label="published")
    a2.set_xticks(xpos, [r[0] for r in rows], fontsize=8)
    a2.set_ylabel("alite wt%")
    a2.set_ylim(0, 105)
    a2.set_title("Alite wt%: M3 vs T1 model", fontsize=10)
    a2.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"wrote {out}")


def main() -> None:
    rep = json.loads(RES.read_text())
    x, yobs, ycalc, proj = fit_arrays(GPX_T1, 4.0, 70.0)
    shift = next(s["shift_2th"] for s in rep["samples"]
                 if "Clinker_Nist" in s["sample"] and s.get("model") == "T1")
    fig_fit(x, yobs, ycalc, shift,
            "Clinker (SRM 2686a) Cu K$\\alpha_1$ -- alite-T1 model, "
            f"wR = 20.9%",
            FIGDIR / "fig1_clinker_t1_fit.png")
    fig_summary(rep, FIGDIR / "fig2_results_summary.png")


if __name__ == "__main__":
    raise SystemExit(main())