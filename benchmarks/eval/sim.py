"""Shared GSAS-II simulated-pattern helper (Unit 06 and friends).

Deterministic protocol: one phase from CIF, Cu Ka, 15-140 deg 2theta,
0.02 deg step, scale 50000, 3-cycle background+scale refine, return ycalc.
Same recipe as the unit02/unit04 halite sims.
"""
from __future__ import annotations

import os
import sys

from core.ingest import InstrumentParams, PowderPattern


def sim_cif_to_pattern(cif_path: str, work_dir: str, *,
                       anode: str = "CuKa",
                       wavelengths=(1.5405, 1.5443),
                       tmin: float = 15.0, tmax: float = 140.0,
                       step: float = 0.02, scale: float = 50000.0,
                       prm_path: str = "") -> PowderPattern:
    """Simulate a powder pattern from a CIF via GSAS-II scriptable.

    Returns the PowderPattern (tth, ycalc, instrument) with provenance in
    ``metadata``. Needs the vendored GSAS-II on sys.path (callers handle it).
    """
    import numpy as np
    from GSASII.GSASIIscriptable import G2Project

    os.makedirs(work_dir, exist_ok=True)
    tag = os.path.splitext(os.path.basename(cif_path))[0]
    proj = G2Project(newgpx=os.path.join(work_dir, f"{tag}_sim.gpx"))
    proj.add_phase(cif_path, phasename=tag, fmthint="CIF")
    hist = proj.add_simulated_powder_histogram(
        f"sim_{tag}", prm_path or "", tmin, tmax, Tstep=step, scale=scale,
        phases=proj.phases())
    hist.set_refinements(
        {"Background": {"type": "chebyschev-1", "no. coeffs": 3,
                        "refine": True}})
    proj.phase(0).set_HAP_refinements({"Scale": True})
    proj.data["Controls"]["data"]["max cyc"] = 3
    proj.refine()
    x = np.asarray(hist.getdata("x"))
    y = np.asarray(hist.getdata("ycalc"))
    proj.save(os.path.join(work_dir, f"{tag}_sim_final.gpx"))
    return PowderPattern(
        sample_name=tag, source=f"sim:{tag}",
        tth=x, intensity=y,
        instrument=InstrumentParams(anode=anode, wavelengths=wavelengths,
                                    tmin=tmin, tmax=tmax, step=step,
                                    npts=int(x.size)),
        metadata={"sim_protocol": "unit06-v1",
                  "simulated_from": os.path.basename(cif_path),
                  "scale": scale, "cycles": 3})


def _vendor_gsasii(vendor: str) -> None:
    """Bootstrap the vendored GSAS-II copy (network, one-time) if absent.

    Fresh clones never ship .vendor/ (gitignored); first run downloads the
    official GSAS-II repository into ``.vendor/GSAS-II`` so every refinement
    path works out of the box (see docs/installation.md).
    """
    marker = os.path.join(vendor, "GSASII", "GSASIIscriptable.py")
    if os.path.exists(marker):
        return
    os.makedirs(vendor, exist_ok=True)
    os.environ.setdefault("GSAS_VENDOR_URL",
                          "https://github.com/GSAS-II/GSAS-II.git")
    url = os.environ["GSAS_VENDOR_URL"]
    print(f"[vendor] GSAS-II missing at {vendor}; cloning {url} ...",
          flush=True)
    import subprocess
    subprocess.run(["git", "clone", "--depth", "1", url, vendor],
                   check=True)
    if not os.path.exists(marker):
        raise RuntimeError(f"GSAS-II clone did not produce {marker}")


def ensure_gsasii(root: str, vendor: str, prm: str) -> str:
    """Put vendored GSAS-II on sys.path; returns the PRM path to use."""
    _vendor_gsasii(vendor)
    sys.path.insert(0, root)
    sys.path.insert(0, vendor)
    sys.path.insert(0, os.path.join(vendor, "GSASII"))
    sys.path.insert(0, os.path.join(vendor, "bin"))
    os.environ.setdefault("GSASIIDATA", os.path.join(vendor, "GSASII"))
    if not prm or not os.path.exists(prm):
        return os.path.join(root, "data", "unit00", "input", "INST_XRY.PRM")
    return prm