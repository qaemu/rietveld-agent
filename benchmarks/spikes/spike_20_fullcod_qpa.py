"""Spike 20: FULL-COD exhaustive QPA — 20 published multi-phase datasets.

End-to-end: raw powder data -> fingerprint -> exhaustive COD geometry
screen (the complete ~530k-entry COD index, no curated guess list, no
chemistry priors) -> per-candidate fast full-pattern intensity rank (CIFs
fetched online, generated via gemmi structure factors) -> top-N phase
hypothesis -> staged GSAS-II Rietveld QPA (Hill-Howard) -> vs weighed /
published truth.

Datasets (all externally published, raw diffractometer data):
  1-8    IUCr CPD QARR 1a-1h (corundum/zincite/fluorite weighed)
  9      IUCr CPD QARR sample 2 (corundum/zincite/fluorite/brucite)
  10     IUCr CPD QARR sample 3 (corundum/zincite/fluorite + 29.5% glass)
  11     IUCr CPD QARR sample 4 (corundum/magnetite/zircon)
  12     IUCr CPD QARR synthetic bauxite (7-phase)
  13-16  NIST SRM 2686a clinker (4 patterns, published wt%)
  17-19  Fe2O3/Fe3O4 30-70, 50-50, 70-30 (Stoops/Praetz; weigh-in basis)
  20     Natural magnetite ore "MexicanMagnetite"

Gate: every hypothesis phase present in truth at >=5 wt% must be
identified, and |refined wt% - truth wt%| <= 3 (absolute wt% points),
QARR-consistent, truth-normalized for the glass sample.
"""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "benchmarks" / "spikes"))

from core.codindex import (              # noqa: E402
    D_MIN, D_MAX, D_UNIT, HKL_CAP, N_WORKERS,
    build_index, cif_calc_lines, intensity_match, load_index, parse_cod_csv,
    save_index, screen_fingerprint,
)
from core.ingest import (                # noqa: E402
    InstrumentParams, PowderPattern, parse_xrdml, sample_fingerprint,
)
from spike_12_cod_full import (          # noqa: E402
    CIF_DIR, download_cif, ensure_cod_csv, ensure_gsasii, ensure_index,
    load_pattern,
)

WORK = ROOT / "data" / "spike20" / "work"
OUT = ROOT / "data" / "spike20" / "results"
CSV = ROOT / "data" / "cod_index" / "cod_metadata_full.csv"
IDX_NPZ = ROOT / "data" / "cod_index" / "cod_full_v1.npz"
IDX_JSON = ROOT / "data" / "cod_index" / "cod_full_v1.meta.json"

STRUCT = ROOT / "data" / "structures"

# ----------------------------------------------------------------------------
# manifest: (sample id, file, instrument meta, truth label dict)
# ----------------------------------------------------------------------------
QARR = ROOT / "data" / "benchmark" / "qarr" / "lhpm"
IRON = ROOT / "data" / "benchmark" / "ironox"
SRM_IN = ROOT / "data" / "spike11" / "input"

CU = InstrumentParams(anode="CU", wavelengths=(1.54056, 1.54439),
                       scan_axis="2Theta/Theta")
FE = InstrumentParams(anode="CU", wavelengths=(1.540598, 1.544426),
                      scan_axis="2Theta/Theta")
SYNC = InstrumentParams(anode="SYNC", wavelengths=(0.82543,),
                        scan_axis="2Theta/Theta")


def _wl(instr):
    wl = instr.wavelengths
    return (float(wl[0]), float(wl[1]) if len(wl) > 1 else 0.0,
            0.5 if len(wl) > 1 else 0.0)

QARR_TRUTH = {
    "qarr_1a": {"corundum": 1.15, "zincite": 4.04, "fluorite": 94.81},
    "qarr_1b": {"corundum": 94.31, "zincite": 1.36, "fluorite": 4.33},
    "qarr_1c": {"corundum": 5.04, "zincite": 93.59, "fluorite": 1.36},
    "qarr_1d": {"corundum": 13.53, "zincite": 32.89, "fluorite": 53.58},
    "qarr_1e": {"corundum": 55.12, "zincite": 15.25, "fluorite": 29.62},
    "qarr_1f": {"corundum": 27.06, "zincite": 55.22, "fluorite": 17.72},
    "qarr_1g": {"corundum": 31.37, "zincite": 34.21, "fluorite": 34.42},
    "qarr_1h": {"corundum": 35.12, "zincite": 30.19, "fluorite": 34.69},
    "qarr_2": {"corundum": 21.27, "zincite": 19.94, "fluorite": 22.53,
               "brucite": 36.26},
    "qarr_3": {"corundum": 30.79, "zincite": 19.68, "fluorite": 20.06,
               "glass": 29.47},   # amorphous; crystalline-normalize for QPA
    "qarr_4": {"corundum": 50.46, "magnetite": 19.64, "zircon": 29.90},
    "bauxite": {"quartz": 5.16, "boehmite": 14.93, "anatase": 2.00,
                "goethite": 9.98, "kaolinite": 3.02, "gibbsite": 54.90,
                "hematite": 10.00},
}

IRON_TRUTH = {
    "iron_30_70": {"hematite": 31.8, "magnetite": 68.2},
    "iron_50_50": {"hematite": 50.6, "magnetite": 49.4},
    "iron_70_30": {"hematite": 70.5, "magnetite": 29.5},
    "iron_mexican": {"hematite": 73.0, "magnetite": 27.0},
}

SRM_TRUTH = {  # normalized published wt% (García-Maté et al. 2024, Tables)
    "Clinker_Nist_CuKalpha1_R1": {
        "alite": 66.0, "belite": 13.5 + 2.7, "ferrite": 11.1,
        "periclase": 4.0, "aluminate": 0.7 + 1.2, "aphthitalite": 0.8},
    "Silicate_enriched_residue_Nist_CuKalpha1_R1": {
        "alite": 78.7, "belite": 13.4 + 2.9, "periclase": 5.0},
    "aluminate_enriched_residue_clinkerNIST_180718_R1": {
        "ferrite": 68.0, "periclase": 16.8, "aluminate": 7.6 + 5.2,
        "aphthitalite": 2.4},
    "Clinker_Synchrotron": {
        "alite": 65.4, "belite": 13.8 + 3.0, "ferrite": 11.6,
        "periclase": 3.65, "aluminate": 1.99 + 0.57},
}

MANIFEST = []   # filled below: (id, path, instr, truth, fit (lo, hi))

# ----------------------------------------------------------------------------
# phase tag (canonical mineral name) -> (truth key, [cod mineral ^ regexes])
# ----------------------------------------------------------------------------
PHASE_CANON = {
    "corundum":     (r"corundum|aluminum oxide|aluminium oxide|al2o3"),
    "zincite":      (r"zincite|zinc oxide|zno\b"),
    "fluorite":     (r"fluorite|calcium fluoride|caf2"),
    "brucite":      (r"brucite|magnesium hydroxide|mg\(?oh\)?2"),
    "magnetite":    (r"magnetite|fe3o4"),
    "hematite":     (r"hematite|alpha.?fe2o3|fe2o3"),
    "zircon":       (r"zircon|zirconium silicate|zrsio4|zrsi"),
    "quartz":       (r"quartz|alpha.?sio2|low sio2|\bsio2\b"),
    "boehmite":     (r"boehmite|aluminum oxyhydroxide|aluminum oxide hydroxide"),
    "goethite":     (r"goethite|alpha.?feooh|feo\(oh\)"),
    "gibbsite":     (r"gibbsite|aluminum hydroxide|al\(oh\)3|aluminium hydroxide"),
    "kaolinite":    (r"kaolinite|al2si2o5\(oh\)4"),
    "anatase":      (r"anatase|titanium dioxide|tio2"),
    "alite":        (r"alite|tricalcium silicate|ca3sio5|hatrurite"),
    "belite":       (r"belite|dicalcium silicate|ca2sio4|lamite|bredigite|calcium lanthanum"),
    "aluminate":    (r"aluminate|tricalcium aluminate|ca3al2o6|c3a\b"),
    "ferrite":      (r"ferrite|brownmillerite|c4af|dicalcium aluminoferrite|ferrite \("),
    "periclase":    (r"periclase|magnesium oxide|\bmgo\b"),
    "aphthitalite": (r"aphthitalite|glaserite|k3na\(so4\)2|aphthitalite"),
}


def canon_of(cod_mineral: str, cod_chem: str) -> str | None:
    """Map a COD mineral/chemistry string to a canonical phase key."""
    hay = f"{cod_mineral} {cod_chem}".lower()
    for canon, pat in PHASE_CANON.items():
        if re.search(pat, hay):
            return canon
    return None


def sample_id(path: Path) -> str:
    return path.stem


# ----------------------------------------------------------------------------
# loading
# ----------------------------------------------------------------------------
def parse_lhpm(path: Path):
    """LHPM/RIET7: header, range line `start step end MeasureDateTime...`,
    then raw counts (10 per line).  Returns (meta, tth, counts)."""
    txt = path.read_text(errors="replace").splitlines()
    meta = {}
    i = 0
    while i < len(txt):
        line = txt[i].rstrip("\r")
        m = re.match(r"^\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+MeasureDateTime",
                     line)
        if m:
            start, step, end = (float(m.group(g)) for g in (1, 2, 3))
            n = int(round((end - start) / step)) + 1
            i += 1
            vals = []
            while len(vals) < n and i < len(txt):
                for tok in re.findall(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?",
                                      txt[i]):
                    vals.append(float(tok))
                i += 1
            tth = start + step * np.arange(n)
            return meta, tth, np.asarray(vals[:n], dtype=float)
        for key in ("SampleIdent", "Anode", "Alpha1", "Alpha2", "Ratio"):
            if key in line:
                meta[key] = line
        i += 1
    raise ValueError(f"{path.name}: no LHPM range line")


def load_any(path: Path, instr=InstrumentParams(anode="CU",
                                                wavelengths=(1.54056,))):
    """Return (pattern) from .dat (LHPM or 2-col), .xrdml."""
    head = path.read_text(errors="replace")[:200]
    if "MeasureDateTime" in head or "SampleIdent" in head:
        meta, tth, y = parse_lhpm(path)
        return PowderPattern(sample_name=path.stem, source="qarr",
                             tth=tth, intensity=y, instrument=instr)
    pat = load_pattern(path)
    if pat is None:
        # 2-col / 3-col ascii fallback
        a = np.loadtxt(path)
        pat = PowderPattern(sample_name=path.stem, source="ascii",
                            tth=a[:, 0], intensity=a[:, 1], instrument=instr)
    return pat


# ----------------------------------------------------------------------------
# step 1: exhaustive screen + per-candidate structure-factor rank
# ----------------------------------------------------------------------------
def _dgrid(tth, wl):
    st = np.sin(np.deg2rad(tth) / 2.0)
    with np.errstate(divide="ignore"):
        d = wl / (2.0 * st)
    return d


def lp_factor(d, wl, mono2th=26.6):
    """Lorentz-polarization factor for a powder diffractometer; graphite
    monochromator angle default (26.6 deg 2theta). mono2th<=0 -> no
    monochromator term."""
    st = wl / (2.0 * d)
    st = np.clip(st, 1e-9, 1.0)
    cth = np.sqrt(np.clip(1.0 - st * st, 0.0, 1.0))
    t2 = np.arcsin(st) * 2.0
    c2 = np.cos(t2)
    if mono2th > 0.0:
        cm = np.cos(np.radians(mono2th))
        lp = (1.0 + c2 * c2 * cm * cm) / (st * st * cth)
    else:
        lp = (1.0 + c2 * c2) / (st * st * cth)
    return lp


def fit_candidate(y, tth, wl, du, I, dshifts, la2=0.0,
                  ratio=0.0, lb=0.0, kb=0.12, c_w=0.00105, mono2th=26.6):
    """1-parameter full-pattern profile fit of calc lines (du, I) to
    (tth, y): each line is broadened to a Gaussian of FWHM ~ fwhm_d in
    d-space (peak-width emulation), a small relative d-shift is scanned
    (cell tolerance), and an LSQ scale is fitted per shift. Returns
    (best scale, relative chi2 = chi2/sum(y^2), best dshift)."""
    la1 = wl
    d = _dgrid(tth, wl)
    m = (d > 1.05) & (d < 24.0)
    d, y0 = d[m], y[m]
    dcalc = du / 100.0
    I = I * lp_factor(dcalc, la1, mono2th)
    I = I / (I.max() or 1.0)
    # physical widths: FWHM_d(d) ~ c_w * d^2 (constant angular width)
    w2 = (c_w * dcalc) ** 2
    best = (None, np.inf, 0.0)
    for sh in dshifts:
        dc = dcalc * (1.0 + sh)
        if dc.size > 4000:          # huge-line organics: coarse sample
            keep = np.linspace(0, dc.size - 1, 4000).astype(int)
            dc, I = dc[keep], I[keep]
        D = (d[:, None] - dc[None, :]) ** 2 / (w2 * (1.0 + sh) ** 2)
        Ic = np.exp(-D).dot(I)
        # alpha2 doublet (d_obs = d * la2/la1) + K-beta ghost
        if la2 > 0.0:
            dc2 = dc * la2 / la1
            D2 = (d[:, None] - dc2[None, :]) ** 2 / (w2 * (la2 / la1) ** 2)
            Ic = Ic + ratio * np.exp(-D2).dot(I)
        if lb > 0.0:
            dck = dc * lb / la1
            Dk = (d[:, None] - dck[None, :]) ** 2 / (w2 * (lb / la1) ** 2)
            Ic = Ic + kb * np.exp(-Dk).dot(I)
        ones = np.ones_like(Ic)
        A = np.vstack([Ic, ones]).T
        # least squares [s, b]: A @ [s,b] ~ y0
        M = A.T @ A
        try:
            coef = np.linalg.solve(M, A.T @ y0)
        except np.linalg.LinAlgError:
            continue
        s0, b0 = float(coef[0]), float(coef[1])
        if s0 < 0.0:
            s0 = 0.0
        res = y0 - s0 * Ic - b0
        r2 = float(np.dot(res, res)) / (float(np.dot(y0, y0)) or 1.0)
        if r2 < best[1]:
            best = (s0, r2, sh, b0)
    return best


def strip_round(pat_y, tth, wl, cod_ids, d_units, entry_of, metas_by_id,
                pool_k, dshifts, wl2=0.0, wlr=0.0, wlb=0.0, mono2th=26.6,
                n_strong=3, prom=0.002, exclude=()):
    """Screen current pattern against the FULL index, rank the pool by the
    1-parameter profile fit (d-shift scanned). Returns (pick, hits_sorted).

    ``exclude``: d-spacings of already-picked phases' strong lines; pool
    entries whose line set is essentially a subset of the excluded lines
    (>=85% overlap) are near-isostructural decoys and are skipped without
    a CIF download (silicon / mngse / naybf class).
    ``prom``: fingerprint prominence; 0.002 keeps weak true minors (zincite
    at ~4 wt% is rank ~500-600 in the flu-subtracted residual) in the pool.
    """
    fp = sample_fingerprint(PowderPattern(
        sample_name="residual", source="qarr", tth=tth, intensity=pat_y,
        instrument=InstrumentParams(anode="CU", wavelengths=(wl,))),
        prominence=prom)
    hits = screen_fingerprint(fp, cod_ids, d_units, entry_of, metas_by_id,
                              top_k=pool_k)
    # per-entry index line lookup (entry_of is grouped by d; build the
    # entry -> line-range map once per call)
    order = np.argsort(entry_of, kind="stable")
    seo = entry_of[order]
    sdu = d_units[order]
    n_ent = len(cod_ids)
    bounds = np.searchsorted(seo, np.arange(n_ent + 1), side="left")
    if exclude:
        sel = sorted(exclude)
        scored_pool = []
        for h in hits:
            cid = int(h["cod_id"])
            i = int(np.searchsorted(cod_ids, cid))
            if i >= n_ent or int(cod_ids[i]) != cid:
                scored_pool.append(h)
                continue
            lo, hi = bounds[i], bounds[i + 1]
            if hi - lo < 4:
                scored_pool.append(h)
                continue
            ds = sdu[lo:hi] / 100.0
            n_ov = 0
            for d in ds:
                j = np.searchsorted(sel, d)
                if (j < len(sel) and abs(sel[j] - d) <= 0.02) or \
                        (j > 0 and abs(sel[j - 1] - d) <= 0.02):
                    n_ov += 1
            if n_ov / (hi - lo) >= 0.85:
                continue
            scored_pool.append(h)
        hits = scored_pool
    scored = []
    n_dl = 0
    for h in hits:
        cod_id = int(h["cod_id"])
        cif = CIF_DIR / f"{cod_id}.cif"
        if not cif.exists() or cif.stat().st_size < 200:
            if n_dl >= 400:          # bound downloads on cache-miss
                continue
            try:
                cif = download_cif(cod_id, delay=0.0)
                n_dl += 1
            except Exception:
                continue
            if cif is None or not cif.exists():
                continue
        try:
            du, I = cif_calc_lines(str(cif))
        except Exception:
            continue
        if (I >= 0.15 * I.max()).sum() < n_strong:
            continue
        got = fit_candidate(pat_y, tth, wl, du, I, dshifts,
                           la2=wl2, ratio=wlr, lb=wlb, mono2th=mono2th)
        if got[0] is None:
            continue
        s0, r2, sh, b0 = got
        meta = metas_by_id.get(cod_id, {}) or {}
        scored.append({"cod_id": cod_id, "r2": r2, "scale": s0,
                       "bkg": b0, "dshift": sh,
                       "mineral": meta.get("mineral", ""),
                       "chem": meta.get("chemname", ""),
                       "canon": canon_of(meta.get("mineral", ""),
                                         meta.get("chemname", ""))})
    scored.sort(key=lambda r: (r["r2"], -r["scale"]))
    return scored


def screen_and_rank(pat, entries, cod_ids, d_units, entry_of, metas_by_id,
                    pool_k=4000, top_n=20, rounds=3, sample=""):
    """Sequential phase stripping over the FULL COD index. Each round:
    fingerprint (prominence 0.002 -> weak true minors stay in the pool) ->
    exhaustive screen -> full-pattern fit rank -> subtract the winning calc
    from the pattern.  Pool entries whose lines are >=85% a subset of the
    already-picked phases' strong lines (near-isostructural decoys: silicon,
    mngse, naybf, sphalerite class) are skipped.  Union of winners +
    runners-up of ALL rounds is the hypothesis set (the true phases,
    isostructural decoys and all)."""
    tth, y = pat.tth, pat.intensity
    wl = float(pat.instrument.wavelengths[0])
    wl2 = float(pat.instrument.wavelengths[1]) if len(
        pat.instrument.wavelengths) > 1 else 0.0
    wlr = 0.5 if wl2 > 0.0 else 0.0
    wlb = 1.3922 if wl2 > 0.0 else 0.0     # Cu K-beta ghost
    mono2th = 0.0 if wl2 <= 0.0 else 26.6  # sync: no monochromator
    dshifts = tuple(d / 1000.0 for d in (-2.5, -2, -1.5, -1, -0.5, 0,
                                         0.5, 1, 1.5, 2, 2.5))
    resid = y.astype(float).copy()
    picks, runners = [], []
    sel_d = set()          # strong lines (I>=0.3) of the picked phases
    t0 = time.time()
    for rnd in range(rounds):
        scored = strip_round(resid, tth, wl, cod_ids, d_units, entry_of,
                             metas_by_id, pool_k, dshifts, wl2=wl2,
                             wlr=wlr, wlb=wlb, mono2th=mono2th,
                             exclude=sel_d)
        if not scored:
            break
        best = scored[0]
        picks.append(best)
        runners.extend(scored[1:])
        # subtract best phase lines (scaled by the fitted scale)
        cif = CIF_DIR / f"{best['cod_id']}.cif"
        try:
            du, I = cif_calc_lines(str(cif))
            for d, i in zip(du, I):        # extend the exclusion line set
                if i >= 0.3:
                    sel_d.add(round(float(d) / 100.0, 3))
            dc0 = du / 100.0
            I = I * lp_factor(dc0, wl, mono2th)
            I = I / (I.max() or 1.0)
            dc = dc0 * (1.0 + best["dshift"])
            d = _dgrid(tth, wl)
            w2 = (0.004 / 2.3548) ** 2
            Ic = np.exp(-((d[:, None] - dc[None, :]) ** 2) / w2).dot(I)
            if wl2 > 0.0:
                Ic += wlr * np.exp(-((d[:, None] - dc[None, :] * wl2 / wl)
                                     ** 2) / w2).dot(I)
                Ic += 0.12 * np.exp(-((d[:, None] - dc[None, :] * wlb / wl)
                                      ** 2) / w2).dot(I)
            resid -= best["scale"] * Ic
            resid -= best.get("bkg", 0.0)
            resid = np.clip(resid, 0.0, None)
        except Exception:
            pass
        print(f"    round {rnd + 1}: pick {best['cod_id']} "
              f"r2={best['r2']:.4f} scale={best['scale']:.0f} "
              f"({best['mineral'][:16]} {best['chem'][:16]})",
              flush=True)
    print(f"    strip done in {time.time() - t0:.0f}s", flush=True)
    # hypothesis: winners (all rounds) + next round-1 runners, dedupe by id
    hyp = []
    seen = set()
    for r in picks + runners:
        c = r["cod_id"]
        if c in seen or len(hyp) >= top_n:
            continue
        seen.add(c)
        hyp.append(r)
    return hyp


# ----------------------------------------------------------------------------
# step 2: staged GSAS-II RQPA (spike-15 protocol: scales -> cells -> shapes,
# lst-parsed convergence, Hill-Howard extraction)
# ----------------------------------------------------------------------------
BKG_COEFFS = 8
BAD_LST = ("Refinement failed", "Invalid metric tensor",
           "Refinement appears to be stuck", "singular")


def cell_volume(cell6) -> float:
    a, b, c, al, be, ga = (float(x) for x in cell6[:6])
    ca, cb, cg = np.cos(np.radians((al, be, ga)))
    return a * b * c * np.sqrt(1.0 - ca ** 2 - cb ** 2 - cg ** 2
                               + 2.0 * ca * cb * cg)


def _clone_prm(sample: str, work: Path, tag: str, a1: float, a2: float,
               ratio: float, sync: bool = False) -> Path:
    """Clone the spike-16 protocol PRM for the sample geometry and patch the
    wavelengths/ratio (GSAS-II keeps the rest of the calibration)."""
    src = (ROOT / "data" / "spike16" / "work" / "INST_SYNC_PROTOCOL.PRM"
           if sync else ROOT / "data" / "spike16" / "work" /
           "INST_CU_PROTOCOL.PRM")
    txt = src.read_text()
    txt = re.sub(
        r"INS\s+1 ICONS.*",
        f"INS  1 ICONS  {a1:.6f} {a2:.6f} {ratio:.4f}  0     0.7    0     0.5",
        txt)
    prm = work / f"{tag}.prm"
    prm.write_text(txt)
    return prm


def gsas_qpa(pat, phases, work: Path, tag: str, lo: float, hi: float,
             a1: float, a2: float, ratio: float, sync: bool = False,
             maxcyc: int = 40, dwr_gate: float = 0.5, n_resid: int = 3,
             resid_sig: float = 2.5, resid_dwr_gate: float = 0.05,
             n_resid_rounds: int = 3, resid_d_tol: float = 0.005,
             cod_ids=None, d_units=None, entry_of=None, metas_by_id=None):
    """Four-stage GSAS-II RQPA.

    Stage A (model selection): refine EVERY hypothesis phase STAND-ALONE
    (scale + background + shift + U,V,W,X,Y profile). The Rietveld wR is the
    gold-standard discriminator between isostructural phases; the winner's
    refined profile becomes the FROZEN profile for all later stages.

    Stage B (forward selection, frozen profile): add the best remaining
    candidate one at a time; accept only if wR improves by > ``dwr_gate``
    (0.5).  A frozen profile keeps wrong isostructural phases from absorbing
    scale by distorting the peak shape (empirically silicon collapses from
    ~17% to ~4%).

    Stage C (residual screening): yobs - ycalc of the clean base fit;
    fingerprint the residual and screen the full COD index for weak TRUE
    minors (e.g. corundum at ~1%) that Stage B cannot distinguish from
    decoys by wR.  Pre-filter: skip hypothesis members, near-isostructural
    decoys (>=90% line overlap) and candidates without a distinctive
    low-angle (d >= 2.0 A) line matching a residual peak; then verify each
    remaining candidate by a frozen-profile delta-wR fit (gate 0.05).

    Stage D (final): joint fit with the frozen profile, Hill-Howard wt%
    via Mass*Scale (GSAS-II calcMassFracs; no cell-volume factor).
    """
    from copy import deepcopy
    work.mkdir(parents=True, exist_ok=True)
    prm = _clone_prm(tag, work, tag, a1, a2, ratio, sync=sync)
    from benchmarks.eval.sim import ensure_gsasii
    ensure_gsasii(str(ROOT), str(ROOT / '.vendor' / 'GSAS-II'),
                  str(prm))
    from GSASII.GSASIIscriptable import G2Project

    xye = work / f"{tag}.xye"
    tth, y = pat.tth, pat.intensity
    m = (tth >= lo) & (tth <= hi)
    sig = np.sqrt(np.maximum(y[m], 1.0))
    np.savetxt(xye, np.column_stack([tth[m], y[m], sig]),
               fmt="%.5f %.3f %.3f")

    def _mkproj(phases_sub, tag2, inst=None):
        gpx = work / f"{tag2}.gpx"
        for sfx in (".gpx", "_final.gpx", ".lst", ".bak0.gpx"):
            st = str(gpx).replace(".gpx", sfx)
            if Path(st).exists():
                Path(st).unlink()
        proj = G2Project(newgpx=str(gpx))
        for p in phases_sub:
            proj.add_phase(str(p["cif"]), phasename=p["name"],
                           fmthint="CIF")
        h = proj.add_powder_histogram(str(xye), iparams=str(prm),
                                      phases=[p["name"] for p in phases_sub])
        if inst is not None:
            h.data['Instrument Parameters'][0] = deepcopy(inst)
        h.set_refinements({"Limits": {"low": lo, "high": hi}})
        h.set_refinements({"Background": {"type": "chebyschev-1",
                                          "no. coeffs": 5, "refine": True}})
        h.set_refinements({"Sample Parameters": ["Shift"]})
        proj.data["Controls"]["data"]["max cyc"] = maxcyc
        return proj, h

    def _scales_on(proj, names):
        for i in range(len(names)):
            proj.phase(i).set_HAP_refinements({"Scale": True})

    def _wR(proj, h, tag2):
        """run one refine + return (wR, lst_text, conv, bad); None wR on fail"""
        proj.refine(makeBack=False)
        lstp = work / f"{tag2}.lst"
        if not lstp.exists():
            lstp = Path(str(work / tag2).replace(".gpx", ".lst") + ".lst")
        txt = lstp.read_text(errors="ignore") if lstp.exists() else ""
        wr = (re.findall(r"Final refinement wR =\s*([\d.]+)", txt)
              or [None])[-1]
        wr = float(wr) if wr is not None else None
        conv = bool(("Refinement successful" in txt)
                    or ("Final refinement" in txt))
        bad = any(mm in txt for mm in BAD_LST)
        return wr, txt, conv, bad

    stage_log = []
    # ---- Stage A: stand-alone model selection (with profile refine) ----
    solo = []
    for p in phases:
        try:
            proj, h = _mkproj([p], f"{tag}_solo_{p['name']}", None)
            _scales_on(proj, [p["name"]])
            h.set_refinements({"Instrument Parameters": ["U", "V", "W",
                                                         "X", "Y"]})
            wr, txt, conv, bad = _wR(proj, h, f"{tag}_solo_{p['name']}")
            solo.append({"phase": p, "wR": wr, "conv": conv, "bad": bad,
                         "proj": proj, "h": h})
            print(f"    solo {p['name']} ({p['canon']}): wR={wr} "
                  f"conv={conv} bad={bad}", flush=True)
        except Exception as e:
            solo.append({"phase": p, "wR": None, "error": str(e)})
            print(f"    solo {p['name']}: ERROR {e}", flush=True)
    ok = [x for x in solo if x["wR"] is not None]
    ok.sort(key=lambda x: x["wR"])
    if not ok:
        return {"wt": [], "wR": None, "rwp": None, "converged": False,
                "bad": True, "stage_log": stage_log, "model_select": solo}
    winner = ok[0]
    # freeze the winner's refined instrument profile (values [1], flags [2])
    frozen = deepcopy(winner["h"].data['Instrument Parameters'][0])
    for k, v in frozen.items():
        if isinstance(v, list) and len(v) > 2 and isinstance(
                v[2], (bool, np.bool_)):
            v[2] = False
    stage_log.append({"stage": "A_winner",
                      "phase": winner["phase"]["name"],
                      "canon": winner["phase"]["canon"],
                      "wR": winner["wR"]})

    # ---- Stage B: forward selection, frozen profile ----
    selected = [winner["phase"]]
    base_wr = winner["wR"]

    def _sel_lines(phases_sub) -> set:
        """strong lines (I>=0.3) of a phase list, in Angstrom."""
        out = set()
        for s in phases_sub:
            try:
                du, ii = cif_calc_lines(str(s["cif"]), dmin=1.0, dmax=22.0)
                for d, i in zip(du, ii):
                    if i >= 0.3:
                        out.add(round(float(d) / D_UNIT, 3))
            except Exception:
                pass
        return out

    for rnd in range(len(phases) - 1):
        remaining = [x for x in ok if x["phase"]["name"] not in
                     [s["name"] for s in selected]
                     and x["phase"].get("canon") not in
                     [s.get("canon") for s in selected]]
        if not remaining:
            break
        sel_d = _sel_lines(selected)
        best, best_d = None, 0.0
        for x in remaining:
            cand = x["phase"]
            # isostructural pre-filter vs the selected phases' strong lines:
            # silicon / mngse / naybf / sphalerite class never improves the
            # frozen-profile fit (empirical delta-wR < 0.5) -> skip the fit
            try:
                du, ii = cif_calc_lines(str(cand["cif"]), dmin=1.0,
                                        dmax=22.0)
                strong = [float(d) / D_UNIT for d, i in zip(du, ii)
                          if i >= 0.3]
            except Exception:
                strong = []
            if strong and sel_d:
                n_ov = sum(any(abs(d - s) <= 0.02 for s in sel_d)
                           for d in strong)
                if n_ov / len(strong) >= 0.5:
                    print(f"    fwd {cand['name']} ({cand['canon']}): "
                          f"isostructural skip ({n_ov}/{len(strong)})",
                          flush=True)
                    continue
            try:
                proj, h = _mkproj(selected + [cand], f"{tag}_fwd", frozen)
                _scales_on(proj, [s["name"] for s in selected] +
                           [cand["name"]])
                wr, txt, conv, bad = _wR(proj, h, f"{tag}_fwd")
            except Exception as e:
                wr, conv, bad = None, False, True
                print(f"    fwd {cand['name']}: ERROR {e}", flush=True)
            d = (base_wr - wr) if (wr is not None and conv and not bad) \
                else 0.0
            print(f"    fwd {cand['name']} ({cand['canon']}): wR={wr} "
                  f"d={d:+.3f}", flush=True)
            if wr is not None and d > best_d:
                best, best_d = cand, d
        if best is None or best_d <= dwr_gate:
            break
        selected.append(best)
        stage_log.append({"stage": f"B_accept_r{rnd + 1}",
                          "phase": best["name"], "canon": best["canon"],
                          "d_wR": round(best_d, 3)})
        # refit the accepted base set with the frozen profile
        proj, h = _mkproj(selected, f"{tag}_base", frozen)
        _scales_on(proj, [s["name"] for s in selected])
        base_wr, _, conv, bad = _wR(proj, h, f"{tag}_base")
        print(f"    base {[s['name'] for s in selected]}: wR={base_wr}",
              flush=True)
        if base_wr is None or bad:
            break
    stage_log.append({"stage": "B_done",
                      "selected": [s["name"] for s in selected],
                      "base_wR": base_wr})
    # always refit the accepted base set (frozen) so proj/h/base_wr reflect
    # the PURE selected-set model for the Stage C residual screening (the
    # forward-trial fits are contaminated by their candidate phases)
    proj, h = _mkproj(selected, f"{tag}_base", frozen)
    _scales_on(proj, [s["name"] for s in selected])
    base_wr, _, conv_b, bad_b = _wR(proj, h, f"{tag}_base")
    print(f"    base {[s['name'] for s in selected]}: wR={base_wr}",
          flush=True)

    # ---- Stage C v4: residual screening for weak true minors ----
    # Principle: rerank against the full COD index for UNMODELED intensity
    # in the clean residual (base fit, frozen profile) with an EXACT-
    # DISTINCTIVE pass: a candidate's line counts only if it sits at
    # d >= 2.0 A, outside every strong line of the selected phases
    # (0.02 A) and within 0.004 A of a residual fingerprint peak.  The
    # 0.004 A tolerance is what separates real phases from coincidental
    # riders (zincite/corundum index lines sit within ~0.001-0.003 of
    # their residual peaks; a 0.02 A pass admits ~10580 riders, 0.004 A
    # only 2771).  Candidates are queued by matched INTENSITY (sum of
    # the heights of their distinctive peaks) with duplicate line-sets
    # collapsed to one representative (the zincite-variant cluster is one
    # phase).  Each round fits the top-6 queued candidates against the
    # frozen profile and accepts ONLY the best delta-wR (junk riders that
    # coincide with a single peak collapse: nantokite-class +0.000;
    # zincite +2.9, corundum +0.19); the accepted phase is refit into
    # the base and the residual recomputed, so the next round's passers
    # automatically lose the accepted phase's peaks and the next true
    # minor surfaces.  Residual peaks within 0.025 A of a selected
    # phase's own strong lines are masked (frozen-profile misfit ridges
    # such as fluorite's 3.173 -- a whole unrelated-entry family rides
    # them and would crowd true minors out of the fit queue).
    resid_new = []
    if cod_ids is not None and d_units is not None and entry_of is not None \
            and metas_by_id is not None and base_wr is not None and not bad_b:
        try:
            id_to_i = {int(c): i for i, c in enumerate(cod_ids)}
            order = np.argsort(entry_of, kind="stable")
            se = entry_of[order]
            sd = d_units[order] / D_UNIT
            bounds = np.searchsorted(se, np.arange(len(cod_ids) + 1),
                                     side="left")
            arr = h.data['data'][1]  # 6 x N masked: 0 2th, 1 yobs, 3 ycalc
            tth_f = np.asarray(arr[0], dtype=float)
            yobs_f = np.asarray(arr[1], dtype=float)
            ycalc_f = np.asarray(arr[3], dtype=float)

            have_canon = set(s.get("canon") for s in selected
                             if s.get("canon"))
            have_cod = set(str(s["cif"].stem) for s in selected)
            tested = set(have_cod) | set(have_canon)
            for ph in phases:          # hypothesis already tested in A/B
                tested.add(str(ph["cif"].stem))
                if ph.get("canon"):
                    tested.add(ph["canon"])

            # strong-line catalog of the selected phases (base set only;
            # accepted residual phases are removed by the residual refit)
            sel_d = set()
            for s in selected:
                try:
                    du, ii = cif_calc_lines(str(s["cif"]), dmin=1.0,
                                            dmax=22.0)
                    sel_d.update(round(float(d) / D_UNIT, 3)
                                 for d, i in zip(du, ii) if i >= 0.3)
                except Exception:
                    pass

            resid_rnd = 0
            while len(resid_new) < n_resid and resid_rnd < n_resid_rounds:
                resid_rnd += 1
                resid = np.clip(yobs_f - ycalc_f, 0.0, None)
                rpat = PowderPattern(sample_name=f"{tag}_resid",
                                     source=pat.source, tth=tth_f,
                                     intensity=resid, instrument=pat.instrument)
                fp = sample_fingerprint(rpat, prominence=0.002)
                # deep pool (one screen; the pool is fixed across rounds,
                # the passers are re-derived each round from the new
                # residual peak set)
                if resid_rnd == 1:
                    hits = screen_fingerprint(fp, cod_ids, d_units, entry_of,
                                              metas_by_id, top_k=15000,
                                              pool_k=30000)
                rpk_d = np.array([p.d for p in fp.peaks])
                rpk_h = np.array([p.height for p in fp.peaks])
                # MASK residual peaks produced by the selected phases' own
                # lines (frozen-profile misfit: e.g. fluorite 111 at 3.154
                # leaves a 3.173 ridge that a whole family of unrelated
                # entries' lines coincide with; those are NOT evidence for
                # a new phase and would crowd true minors out of the fit
                # queue).  A residual peak within 0.025 A of a strong
                # selected line is dropped.
                if sel_d:
                    keep = np.array([
                        not any(abs(float(pd) - s) <= 0.025 for s in sel_d)
                        for pd in rpk_d])
                    rpk_d = rpk_d[keep]
                    rpk_h = rpk_h[keep]
                rd_sorted = np.sort(rpk_d)

                def _rpk_h(d):
                    # exact line match within RESID_D_TOL.  The COD line
                    # index stores d rounded to 0.001 A, so the clipping
                    # tolerance carries +0.001 A headroom (0.004 nominal +
                    # quantization): zincite's own 101 line (index d 2.480)
                    # sits 0.0041 A from the residual peak at 2.4759.
                    j = np.searchsorted(rd_sorted, d)
                    best = None
                    for k in (j - 1, j):
                        if 0 <= k < len(rd_sorted) \
                                and abs(rd_sorted[k] - d) <= resid_d_tol:
                            best = rpk_h[np.flatnonzero(rpk_d == rd_sorted[k])[0]]
                    return best

                # exact-distinctive passers, deduped by matched line-set
                best_by_set = {}
                for hh in hits:
                    cid = str(hh["cod_id"])
                    if cid in tested:
                        continue
                    canon = hh.get("canon") or hh.get("mineral") or cid
                    if canon in tested:
                        continue
                    if hh.get("significance", 0.0) < resid_sig:
                        continue
                    i = id_to_i.get(int(cid))
                    if i is None:
                        continue
                    ds = sd[bounds[i]:bounds[i + 1]]
                    if ds.size == 0:
                        continue
                    n_ov = sum(any(abs(float(d) - s) <= 0.02 for s in sel_d)
                               for d in ds)
                    if n_ov / ds.size >= 0.5:      # isostructural decoy
                        continue
                    dm = [(float(d), _rpk_h(float(d))) for d in ds
                          if float(d) >= 2.0
                          and not any(abs(float(d) - s) <= 0.02
                                      for s in sel_d)
                          and _rpk_h(float(d)) is not None]
                    if not dm:
                        continue
                    mi = sum(h for _, h in dm)
                    # named = a REAL canonical name (metas minname), not the
                    # cod-id fallback: nameless variants must never displace
                    # the named canonical phase in a line-set tie
                    named = bool(canon != cid and
                                 (hh.get("canon") or hh.get("mineral")))
                    key = (mi, named, hh.get("significance", 0.0))
                    mset = frozenset(round(d, 2) for d, _ in dm)
                    prev = best_by_set.get(mset)
                    if prev is None or key > prev[0]:
                        best_by_set[mset] = (key, cid, canon, dm, mi)
                queue = sorted(best_by_set.values(),
                               key=lambda t: (-t[0][0], -t[0][1], -t[0][2]))
                print(f"    C r{resid_rnd}: passers={len(queue)} "
                      f"rpk={len(rpk_d)} peaks", flush=True)
                for q in queue[:6]:
                    print(f"    C r{resid_rnd} cand {q[1]} ({q[2][:20]}) "
                          f"mi={q[4]} lines="
                          f"{[round(d, 3) for d, _ in q[3]]}", flush=True)
                best_cand, best_d = None, 0.0
                n_fit = 0
                for (_key, cid, canon, _dm, _mi) in queue:
                    if n_fit >= 6:
                        break
                    cif = download_cif(int(cid), delay=0.05)
                    if cif is None or not cif.exists():
                        continue
                    cand = {"name": f"p{cid}", "cif": cif, "canon": canon}
                    try:
                        proj2, h2 = _mkproj(selected + [cand], f"{tag}_resf",
                                            frozen)
                        _scales_on(proj2, [s["name"] for s in selected] +
                                   [cand["name"]])
                        wr2, _, conv2, bad2 = _wR(proj2, h2, f"{tag}_resf")
                    except Exception as e:
                        wr2, conv2, bad2 = None, False, True
                        print(f"    C r{resid_rnd} {cid} ({canon}): "
                              f"fit ERROR {e}", flush=True)
                    n_fit += 1
                    tested.add(cid)
                    tested.add(canon)
                    if wr2 is not None and base_wr is not None and conv2 \
                            and not bad2:
                        d2 = base_wr - wr2
                        print(f"    C r{resid_rnd} {cid} ({canon}) "
                              f"d_wR={d2:+.3f} wR={wr2}", flush=True)
                        if d2 > best_d:
                            best_cand, best_d = cand, d2
                    else:
                        print(f"    C r{resid_rnd} {cid} ({canon}) wR={wr2} "
                              f"conv={conv2} bad={bad2}", flush=True)
                if best_cand is None or best_d <= resid_dwr_gate:
                    break
                resid_new.append(best_cand)
                stage_log.append({"stage": f"C_accept_r{resid_rnd}",
                                  "phase": best_cand["name"],
                                  "canon": best_cand["canon"],
                                  "d_wR": round(best_d, 3)})
                # the accepted phase's strong lines join the mask catalog
                try:
                    du, ii = cif_calc_lines(str(best_cand["cif"]), dmin=1.0,
                                            dmax=22.0)
                    sel_d.update(round(float(d) / D_UNIT, 3)
                                 for d, i in zip(du, ii) if i >= 0.3)
                except Exception:
                    pass
                # refit the accepted set -> new base for the next round
                proj, h = _mkproj(selected + resid_new, f"{tag}_base",
                                  frozen)
                _scales_on(proj, [s["name"] for s in selected + resid_new])
                base_wr, _, conv_b, bad_b = _wR(proj, h, f"{tag}_base")
                print(f"    base {[s['name'] for s in selected + resid_new]}: "
                      f"wR={base_wr}", flush=True)
                if base_wr is None or bad_b:
                    break
                arr = h.data['data'][1]
                tth_f = np.asarray(arr[0], dtype=float)
                yobs_f = np.asarray(arr[1], dtype=float)
                ycalc_f = np.asarray(arr[3], dtype=float)
        except Exception as e:
            print(f"    residual screen ERROR {e}", flush=True)
    stage_log.append({"stage": "C_resid",
                      "added": [r["name"] for r in resid_new]})

    # ---- Stage D: final joint fit (frozen profile, scales) ----
    final_phases = selected + resid_new
    proj, h = _mkproj(final_phases, tag, frozen)
    _scales_on(proj, [p["name"] for p in final_phases])
    wr, lst_txt, conv, bad = _wR(proj, h, tag)
    stage_log.append({"stage": "D_final", "wR": wr, "converged": conv,
                      "bad": bad})
    # ---- Hill-Howard ----
    phs = proj.data["Phases"]
    per_phase = []
    for p in final_phases:
        pd = phs[p["name"]]
        mass = float(pd["General"]["Mass"])
        try:
            scale = float(pd["Histograms"][h.name]["Scale"][0])
        except (KeyError, TypeError):
            scale = 1.0
        per_phase.append({"name": p["name"], "cod": str(p["cif"].stem),
                          "canon": p["canon"], "scale": scale, "mass": mass})
    # GSAS-II mass fraction: W_i = Mass_i * Scale_i / sum(Mass_j*Scale_j)
    # (calcMassFracs; NO cell-volume factor)
    smv = {q["name"]: q["scale"] * q["mass"] for q in per_phase}
    tot = sum(smv.values())
    for q in per_phase:
        q["wt"] = round(100.0 * smv[q["name"]] / tot, 2) if tot else 0.0
    per_phase.sort(key=lambda q: -q["wt"])
    return {"wt": per_phase, "wR": wr, "rwp": wr, "converged": conv,
            "bad": bad, "stage_log": stage_log, "model_select": [
                {"name": x["phase"]["name"], "canon": x["phase"]["canon"],
                 "wR": x["wR"]} for x in ok]}

# ----------------------------------------------------------------------------
# gate
# ----------------------------------------------------------------------------
def gate(sample: str, truth: dict, reported) -> dict:
    """reported: list of {canon, wt}.  Returns result dict / gate outcome."""
    # aggregate by canon; unknown canons recorded separately
    agg = {}
    for r in reported:
        c = r.get("canon")
        if c:
            agg[c] = agg.get(c, 0.0) + r["wt"]
    res = {"sample": sample, "truth": truth, "refined": agg, "gaps": [],
           "ok": True}
    # crystalline normalization when sample has known amorphous content
    norm = 1.0
    if "glass" in truth:
        norm = 100.0 / (100.0 - truth["glass"])
    for canon, tgt in truth.items():
        if canon == "glass":
            continue
        tgt_n = tgt * norm
        got = agg.get(canon)
        if tgt_n >= 5.0:
            if got is None:
                res["gaps"].append(f"{canon}: MISSING (truth {tgt_n:.1f})")
                res["ok"] = False
            elif abs(got - tgt_n) > 3.0:
                res["gaps"].append(
                    f"{canon}: {got:.1f} vs {tgt_n:.1f} (Δ={got-tgt_n:+.1f})")
                res["ok"] = False
        else:  # minor phase: must be present if >=1
            if tgt_n >= 1.0 and got is None:
                res["gaps"].append(f"{canon}: minor MISSING (truth {tgt_n:.1f})")
                res["ok"] = False
    return res


# ----------------------------------------------------------------------------
# run
# ----------------------------------------------------------------------------
def build_manifest():
    for name in ("1a", "1b", "1c", "1d", "1e", "1f", "1g", "1h"):
        MANIFEST.append((f"qarr_{name}", QARR / f"cpd-{name}.dat", CU,
                         QARR_TRUTH[f"qarr_{name}"], (5.0, 150.0)))
    for name in ("2", "3", "4"):
        MANIFEST.append((f"qarr_{name}", QARR / f"cpd-{name}.dat", CU,
                         QARR_TRUTH[f"qarr_{name}"], (5.0, 150.0)))
    MANIFEST.append(("bauxite", QARR / "bauxite.dat", CU,
                     QARR_TRUTH["bauxite"], (5.0, 150.0)))
    for name, t in IRON_TRUTH.items():
        f = IRON / ("XRD_2024_04_05_%s.dat" %
                    name.replace("iron_", "").replace("_", "-"))
        if name == "iron_mexican":
            f = IRON / "XRD_2021_06_07_MexicanMagnetite.dat"
        MANIFEST.append((name, f, FE, t, (5.0, 90.0)))
    for name, t in SRM_TRUTH.items():
        f = SRM_IN / f"{name}.dat" if name == "Clinker_Synchrotron" \
            else SRM_IN / f"{name}.xrdml"
        instr = SYNC if "Synch" in name else CU
        lo = 3.0 if "Synch" in name else 4.0
        MANIFEST.append((f"srm_{name}", f, instr, t, (lo, 70.0)))
    return MANIFEST


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--step", choices=["screen", "qpa", "gate"], default=None)
    ap.add_argument("--pool", type=int, default=4000)
    ap.add_argument("--top", type=int, default=20)
    args = ap.parse_args()

    build_manifest()
    selected = [m for m in MANIFEST
                if args.only is None or m[0] in args.only]
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)

    # ---- COD full index (built once; cached to disk) ----
    if not CSV.exists():
        raise SystemExit("cod_metadata_full.csv missing — run fetch first")
    entries = parse_cod_csv(str(CSV))
    print(f"COD entries: {len(entries)}", flush=True)
    if IDX_NPZ.exists():
        (cod_ids, d_units, entry_of, dmin_eff, metas_by_id) = load_index(
            str(IDX_NPZ), str(IDX_JSON))
    else:
        t0 = time.time()
        cod_ids, d_units, entry_of, dmin_eff = build_index(
            entries, workers=N_WORKERS)
        save_index(cod_ids, d_units, entry_of, dmin_eff, entries,
                   str(IDX_NPZ), str(IDX_JSON))
        print(f"index built in {time.time()-t0:.0f}s", flush=True)

    results = {}
    for sid, fpath, instr, truth, (lo, hi) in selected:
        resj = OUT / f"{sid}.json"
        if resj.exists() and args.step != "gate":
            print(f"{sid}: cached", flush=True)
            results[sid] = json.loads(resj.read_text())
            continue
        t0 = time.time()
        print(f"--- {sid} ---", flush=True)
        pat = load_any(fpath, instr)
        if pat is None:
            print(f"{sid}: FAIL load", flush=True)
            results[sid] = {"error": "load"}
            continue
        step1 = screen_and_rank(pat, entries, cod_ids, d_units, entry_of,
                                metas_by_id, pool_k=args.pool,
                                top_n=args.top, sample=sid)
        print(f"{sid}: screen {len(step1)} hits in {time.time()-t0:.0f}s",
              flush=True)
        # hypothesis = stripping winners + runners (isostructural decoys
        # included; GSAS-II drives wrong phases to ~zero scale)
        hyp = step1
        if not hyp:
            print(f"{sid}: FAIL no hypothesis", flush=True)
            results[sid] = {"error": "no hypothesis", "screen": step1}
            continue
        # ensure CIFs on disk & parseable
        ok_phases = []
        for h in hyp:
            p = download_cif(h["cod_id"], delay=0.05)
            if p is None or not p.exists():
                continue
            ok_phases.append({"name": f"p{h['cod_id']}", "cif": p,
                              "canon": h["canon"]})
        if len(ok_phases) < 2:
            print(f"{sid}: FAIL <2 phases", flush=True)
            results[sid] = {"error": "<2 phases", "hyp": hyp}
            continue
        # GSAS-II staged QPA
        t1 = time.time()
        a1, a2, ratio = _wl(instr)
        q = gsas_qpa(pat, ok_phases, WORK / sid, sid, lo, hi,
                     a1, a2, ratio, sync=(a2 == 0.0),
                     cod_ids=cod_ids, d_units=d_units,
                     entry_of=entry_of, metas_by_id=metas_by_id)
        # report phases by canon+labe
        report = q["wt"]
        verdict = gate(sid, truth, report)
        results[sid] = {
            "hypothesis": [{"cod_id": h["cod_id"], "canon": h["canon"],
                            "score": h["score"], "mineral": h["mineral"]}
                           for h in hyp],
            "phases": [{"name": w["name"], "canon": w["canon"],
                        "wt": round(w["wt"], 2)} for w in report],
            "rwp": q["rwp"], "gate": verdict["ok"],
            "gaps": verdict["gaps"], "time_s": round(time.time() - t0, 1),
        }
        resj.write_text(json.dumps(results[sid], indent=2))
        print(f"{sid}: {[(w['canon'], round(w['wt'],1)) for w in report]} "
              f"rwp={q['rwp']} → "
              f"{'PASS' if verdict['ok'] else 'FAIL ' + str(verdict['gaps'])}",
              flush=True)

    # ---- gate summary ----
    if args.step in (None, "gate"):
        npass = nfail = 0
        for sid, _f, _i, t, _r in MANIFEST:
            r = results.get(sid)
            if r is None and (OUT / f"{sid}.json").exists():
                r = json.loads((OUT / f"{sid}.json").read_text())
            if r is None:
                print(f"{sid}: NO RESULT", flush=True)
                nfail += 1
                continue
            if "error" in r and not any(k == "gate" for k in r):
                print(f"{sid}: ERROR {r['error']}", flush=True)
                nfail += 1
            elif r.get("gate"):
                print(f"{sid}: PASS", flush=True)
                npass += 1
            else:
                print(f"{sid}: FAIL {r.get('gaps')}", flush=True)
                nfail += 1
        print(f"\n=== {npass} pass / {nfail} fail / {len(MANIFEST)} ===",
              flush=True)
        return 0 if nfail == 0 else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
