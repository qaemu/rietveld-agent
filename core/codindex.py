"""COD full-database d-line index and position-based screening.

Builds a searchable index from the COMPLETE Crystallography Open Database
metadata export (COD REST API, ``format=csv&id=%``): for every entry with a
valid cell and space group, the systematically allowed hkl reflections are
generated (vectorized general-position extinction per space group via gemmi
symops: a reflection hkl is allowed iff the phase sum over symmetry
operations that fix hkl is non-zero).

Index layout (memory-conscious, bucket-scatter build, no sort):
  * d_units: int16 in 0.01 A units, NON-DECREASING across the whole array,
    so ``np.searchsorted`` yields the full [d0-tol, d0+tol] window in O(log N)
    per query peak;
  * entry_of: int32, entry ordinal aligned with d_units;
  * dmin_eff: int16 per entry -- the ACTUAL d-min used for that entry when the
    cell axis forced the hkl cap up (big-organics honest note);
  * metadata json: full identity (mineral, formula, SG, year, DOI, journal...).

Screening: for every fingerprint peak (d-space, prominence-filtered), find the
window and histogram the matched entry ids weighted by peak height; score =
fraction of fingerprint intensity matched. Top-K reported.

Honest limitations:
  * POSITION-ONLY index (no intensities predicted from the cell): it is a
    candidate FILTER, not an identification;
  * index d-range [1.1, 22] A with an hkl cap of 28 (some huge-cell organics
    effectively start at dmin_eff > 1.1 A, recorded per entry);
  * the verification loop re-scores top candidates against FULL CIF patterns
    (GSAS-II simulation, all structure factors) on the shared d-grid;
  * the CSV export is the COD's curated selection (retracted / duplicate /
    erroneous / theoretical entries excluded by the database itself).

Sources: https://www.crystallography.net/cod/ (CC0); REST API:
https://wiki.crystallography.net/RESTful_API/  (id=% -> full metadata CSV).
"""
from __future__ import annotations

import json
import math
import os
import sys
from multiprocessing import Pool
from pathlib import Path
from typing import List

import numpy as np

from core.ingest import SampleFingerprint
from core.ingest.fingerprint import PEAK_D_TOL

#: index d-range [A]; 0.01 A units -> int16
D_MIN, D_MAX = 1.10, 22.0
D_UNIT = 100
#: hkl index cap (keeps huge-cell organics tractable; per-entry note in index)
HKL_CAP = 28
MAX_LINES_PER_ENTRY = 120_000
N_WORKERS = 8


def _dunit(d: float) -> int:
    return int(round(d * D_UNIT))


class CodEntry:
    """One COD record: identity + cell + space group (as parsed from CSV)."""

    __slots__ = ("cod_id", "mineral", "chemname", "formula", "sg", "sg_number",
                 "a", "b", "c", "alpha", "beta", "gamma",
                 "year", "doi", "journal", "authors", "title")

    def __init__(self, row: dict):
        self.cod_id = int(row["file"])
        self.mineral = (row.get("mineral") or "").strip()
        self.chemname = (row.get("chemname") or "").strip()
        self.formula = (row.get("formula") or "").strip("- ").strip()
        self.sg = (row.get("sg") or "").strip()
        self.year = (row.get("year") or "").strip()
        self.doi = (row.get("doi") or "").strip()
        self.journal = (row.get("journal") or "").strip()
        self.authors = (row.get("authors") or "").strip()
        self.title = (row.get("title") or "").strip()
        try:
            self.sg_number = int(row["sgNumber"]) if row.get("sgNumber") else 0
        except ValueError:
            self.sg_number = 0
        try:
            self.a, self.b, self.c = (float(row[k]) for k in ("a", "b", "c"))
            self.alpha, self.beta, self.gamma = (float(row[k])
                                                 for k in ("alpha", "beta", "gamma"))
        except (ValueError, TypeError):
            self.a = self.b = self.c = 0.0

    @property
    def valid(self) -> bool:
        return (self.sg_number > 0 and self.a > 0 and self.b > 0 and self.c > 0)

    def meta_dict(self) -> dict:
        return {"cod_id": self.cod_id, "mineral": self.mineral,
                "chemname": self.chemname, "formula": self.formula,
                "sg": self.sg, "sg_number": self.sg_number,
                "year": self.year, "doi": self.doi,
                "journal": self.journal, "authors": self.authors,
                "title": (self.title or "")[:240],   # trim: the meta JSON is
                "url": f"https://www.crystallography.net/cod/{self.cod_id}.html"}


def parse_cod_csv(path: str) -> List[CodEntry]:
    """Parse the COD metadata CSV (skip '#' comment lines) into entries."""
    import csv
    rows = [l for l in Path(path).read_text(errors="replace").splitlines()
            if l and not l.startswith("#")]
    out = []
    for row in csv.DictReader(rows):
        e = CodEntry(row)
        if e.valid:
            out.append(e)
    return out


class ReflectionGenerator:
    """hkl -> allowed reflections with exact general-position extinction."""

    def __init__(self, cod_entry: CodEntry):
        self.e = cod_entry
        import gemmi
        self.sg = gemmi.find_spacegroup_by_number(cod_entry.sg_number)
        # General-position structure factor test: F(h) = sum_i exp(-2pi i h . x_i)
        # over the orbit of an irrational general position x0. F = 0  <=>
        # the reflection is systematically absent (lattice centring, screw axes,
        # glide planes and all origins are captured exactly).
        self._x0 = np.array([0.3172103, 0.2819407, 0.1037702])
        self._orbit = []
        for op in self.sg.operations():
            rot = np.asarray(op.rot, dtype=np.int64)
            tran = np.asarray(op.tran, dtype=np.float64) / 24.0
            self._orbit.append((rot @ self._x0) / 24.0 + tran)

    def _allowed_mask(self, h: np.ndarray, k: np.ndarray,
                      l: np.ndarray) -> np.ndarray:
        """Vectorized general-position extinction: keep reflections whose
        single-atom structure factor does not vanish at a general site."""
        hm = np.stack([h, k, l], axis=-1).astype(np.int64)       # (N,3)
        acc = np.zeros(hm.shape[0], dtype=np.complex128)
        for xi in self._orbit:
            acc += np.exp(-2j * np.pi * (hm @ xi))
        return np.abs(acc) > 1e-6

    def line_list(self, dmin: float = D_MIN, dmax: float = D_MAX):
        """Allowed d-units (sorted, unique, int16) + effective dmin used.

        The hkl cap (HKL_CAP) keeps huge-cell entries tractable; when a cell
        axis would need a higher index, the effective dmin for that entry is
        raised accordingly (recorded honestly in the index dmin_eff).
        """
        e = self.e
        raw = [int(np.ceil(2.0 * ax / dmin)) for ax in (e.a, e.b, e.c)]
        capped = [min(HKL_CAP, max(1, r)) for r in raw]
        hmax, kmax, lmax = capped
        effs = [dmin] * 3
        for i, (r, c) in enumerate(zip(raw, capped)):
            if r > c:                       # axis index was capped
                effs[i] = 2.0 * (e.a, e.b, e.c)[i] / c
        dmin_eff = max(dmin, *effs)

        H, K, L = np.meshgrid(np.arange(-hmax, hmax + 1),
                              np.arange(-kmax, kmax + 1),
                              np.arange(-lmax, lmax + 1), indexing="ij")
        H, K, L = H.ravel(), K.ravel(), L.ravel()
        ca, cb, cg = (np.cos(np.radians(v)) for v in
                      (e.alpha, e.beta, e.gamma))
        sa, sb, sg = (np.sin(np.radians(v)) for v in
                      (e.alpha, e.beta, e.gamma))
        vcell = e.a * e.b * e.c * np.sqrt(
            1 - ca ** 2 - cb ** 2 - cg ** 2 + 2 * ca * cb * cg)
        s11 = (e.b * e.c * sa) ** 2
        s22 = (e.a * e.c * sb) ** 2
        s33 = (e.a * e.b * sg) ** 2
        s12 = e.a * e.b * e.c ** 2 * (ca * cb - cg)
        s23 = e.a ** 2 * e.b * e.c * (cb * cg - ca)
        s13 = e.a * e.b ** 2 * e.c * (ca * cg - cb)
        with np.errstate(divide="ignore", invalid="ignore"):
            d2 = vcell ** 2 / (s11 * H ** 2 + s22 * K ** 2 + s33 * L ** 2
                               + 2 * s12 * H * K + 2 * s23 * K * L
                               + 2 * s13 * H * L)
            d = np.sqrt(d2)
        mask = (d >= dmin_eff) & (d <= dmax) & np.isfinite(d)
        H, K, L, d = H[mask], K[mask], L[mask], d[mask]
        keep = self._allowed_mask(H, K, L)
        d = d[keep]
        du = np.sort(np.unique(np.round(d * D_UNIT).astype(np.int64)))
        if du.size > MAX_LINES_PER_ENTRY:
            du = du[:MAX_LINES_PER_ENTRY]
        return du.astype(np.int16), dmin_eff


def _entry_lines(entry: CodEntry):
    return ReflectionGenerator(entry).line_list()


def build_index(entries: List[CodEntry], dmin: float = D_MIN,
                dmax: float = D_MAX, workers: int = N_WORKERS) -> tuple:
    """Build the global sorted d-unit index (bucket-scatter, no sort).

    Returns (cod_ids, d_units int16 non-decreasing, entry_of int32,
             dmin_eff int16 per entry).
    """
    import concurrent.futures as cf
    from core.codsearch import _pool_ctx

    lists = []
    with cf.ProcessPoolExecutor(max_workers=workers,
                                mp_context=_pool_ctx()) as ex:
        for i, (du, dmin_e) in enumerate(
                ex.map(_entry_lines, entries, chunksize=32)):
            lists.append((du, dmin_e))
    n = sum(du.size for du, _ in lists)
    d_units = np.empty(n, dtype=np.int16)
    entry_of = np.empty(n, dtype=np.int32)
    dmin_eff = np.zeros(len(entries), dtype=np.int16)
    # bucket scatter: units are ints in [dmin_u, dmax_u] -> histogram
    dmin_u = _dunit(dmin)
    dmax_u = _dunit(dmax)
    width = dmax_u - dmin_u + 1
    counts = np.zeros(width, dtype=np.int64)
    for du, _ in lists:
        if du.size:
            np.add.at(counts, du.astype(np.int64) - dmin_u, 1)
    starts = np.concatenate([[0], np.cumsum(counts)])
    cursor = starts[:-1].copy()
    for i, (du, dmin_e) in enumerate(lists):
        dmin_eff[i] = int(round(dmin_e * D_UNIT))
        if du.size == 0:
            continue
        idx = du.astype(np.int64) - dmin_u
        pos = cursor[idx]
        d_units[pos] = du
        entry_of[pos] = i
        cursor[idx] += 1
    cod_ids = np.array([e.cod_id for e in entries], dtype=np.int64)
    return cod_ids, d_units, entry_of, dmin_eff


def save_index(cod_ids, d_units, entry_of, dmin_eff, entries: List[CodEntry],
               path: str, meta_path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, cod_ids=cod_ids, d_units=d_units,
                        entry_of=entry_of, dmin_eff=dmin_eff)
    Path(meta_path).write_text(json.dumps([e.meta_dict() for e in entries]))


def load_index(path: str, meta_path: str) -> tuple:
    z = np.load(path)
    metas = json.loads(Path(meta_path).read_text())
    metas_by_id = {m["cod_id"]: m for m in metas}
    return (z["cod_ids"], z["d_units"], z["entry_of"],
            z["dmin_eff"], metas_by_id)


def screen_fingerprint(fp: SampleFingerprint, cod_ids: np.ndarray,
                       d_units: np.ndarray, entry_of: np.ndarray,
                       metas_by_id: dict, top_k: int = 50,
                       tol: float = PEAK_D_TOL,
                       pool_k: int = 60) -> List[dict]:
    """Screen the fingerprint against the whole COD line index.

    Score = intensity-weighted SIGNIFICANCE of the coincidence: an entry's
    observed matched-peak intensity is compared against what its own line
    density would produce by chance (sigma units over the peaks). Dense-line
    organics that trivially coincide with every peak therefore score ~1-3
    sigma, while phases whose *strong* peaks line up score far higher.

    Returns a candidate POOL (union of the top-``pool_k`` by significance,
    by intensity-weighted match ``matched_intensity`` and by ``screen_score``
    coverage; deduped, sig-ordered) -- a pre-filter, NOT the final ranking.
    The final ranking is the caller's job (offline intensity rerank, see
    ``screen_cod``); significance alone overrates sparse-line high-symmetry
    entries whose few d-dense lines happen to sit inside many weak peaks.
    """
    peaks = [(p.d, p.height) for p in fp.peaks]
    if not peaks:
        return []
    n = len(peaks)
    total_h = sum(h for _, h in peaks) or 1e-30
    d_tol_units = int(np.ceil(tol * D_UNIT)) + 1
    # one COUNT per entry per peak (an entry may have several lines in
    # the window); never inflate beyond the peak's own intensity
    counts = np.zeros(len(cod_ids), dtype=np.int32)
    hit_h = np.zeros(len(cod_ids), dtype=np.float64)
    for d0, h0 in peaks:
        d0u = _dunit(d0)
        lo = np.searchsorted(d_units, d0u - d_tol_units, side="left")
        hi = np.searchsorted(d_units, d0u + d_tol_units, side="right")
        if lo >= hi:
            continue
        ids = np.unique(entry_of[lo:hi])
        hit_h[ids] += h0
        counts[ids] += 1
    # random-coincidence model: entry line density rho_e over [D_MIN, D_MAX];
    # per-peak chance p_i = 1 - (1 - w/range)^N_e ~= min(1, w*N_e/range)
    w_ang = 2 * d_tol_units / D_UNIT                    # full window width
    d_range = D_MAX - D_MIN
    n_lines = np.bincount(entry_of, minlength=len(cod_ids)).astype(np.float64)
    mu = np.zeros(len(cod_ids))
    var = np.zeros(len(cod_ids))
    for d0, h0 in peaks:
        p = np.minimum(1.0, w_ang * n_lines / d_range)
        mu += p * h0
        var += p * (1 - p) * h0 * h0
    sig = np.zeros(len(cod_ids))
    nz = var > 0
    sig[nz] = (hit_h[nz] - mu[nz]) / np.sqrt(var[nz])
    m_frac = hit_h / total_h
    cov = counts / n

    def _row(i: int, sig_i: float) -> dict:
        mid = metas_by_id.get(int(cod_ids[i]), {})
        return {
            "cod_id": int(cod_ids[i]),
            "mineral": mid.get("mineral", ""), "chemname": mid.get("chemname", ""),
            "formula": mid.get("formula", ""), "sg": mid.get("sg")
            or str(mid.get("sg_number", "")),
            "year": mid.get("year", ""), "doi": mid.get("doi", ""),
            "journal": mid.get("journal", ""), "title": mid.get("title", ""),
            "url": f"https://www.crystallography.net/cod/{int(cod_ids[i])}.html",
            "matched_peaks": int(counts[i]),
            "matched_intensity": round(float(m_frac[i]), 4),
            "screen_score": round(float(cov[i]), 4),
            "significance": round(float(sig_i), 2),
        }

    # candidate pool: union of the top-k by each criterion (deduped).
    # Deterministic: each criterion's tie-break is the composite key below.
    def comp(i):
        return -(sig[i] * 1e6 + hit_h[i] / total_h)

    picked = {}
    for criterion_idx in (np.argsort(-sig, kind="stable"),
                          np.argsort(-m_frac, kind="stable"),
                          np.argsort(-cov, kind="stable")):
        for i in criterion_idx[:pool_k]:
            if counts[i] == 0:
                continue
            picked[i] = comp(i)
    order = sorted(picked, key=lambda i: picked[i])
    return [_row(i, sig[i]) for i in order[:top_k]]

# --------------------------------------------------------------------------- #
# fast kinematic intensity stage (no GSAS-II): crude form factors, full site
# sum -> automatically includes systematic absences. Used to rank the screen
# candidates before the expensive GSAS-II confirmation sims.
# --------------------------------------------------------------------------- #
# gaussian-atom width per element [A]: rough f_j(s) ~ Z_j * exp(-2 pi^2 s^2 w^2)
_GAUSS_W = {"H": 0.28, "C": 0.30, "N": 0.28, "O": 0.26, "F": 0.26, "Na": 0.33,
            "Mg": 0.33, "Al": 0.33, "Si": 0.33, "P": 0.32, "S": 0.32, "Cl": 0.34,
            "K": 0.38, "Ca": 0.38, "Fe": 0.36, "Zn": 0.34, "Zr": 0.36,
            "Ag": 0.38, "Au": 0.40, "Pb": 0.42, "U": 0.40, "Ba": 0.42,
            "Sr": 0.40, "Ti": 0.36, "Mn": 0.36, "Cu": 0.34, "Ni": 0.34,
            "Co": 0.34, "Cr": 0.34, "B": 0.30, "Li": 0.30, "Be": 0.28,
            "Mo": 0.36, "W": 0.38, "V": 0.34, "Y": 0.38, "La": 0.40,
            "Nd": 0.40, "Pr": 0.40, "Pm": 0.40, "Sm": 0.40, "Eu": 0.40,
            "Gd": 0.40, "Tb": 0.40, "Dy": 0.40, "Ho": 0.40, "Er": 0.40,
            "Tm": 0.40, "Yb": 0.40, "Lu": 0.40, "Sc": 0.36, "Rb": 0.40,
            "Cs": 0.42, "Fr": 0.42, "Ra": 0.44, "Ac": 0.44, "Th": 0.42,
            "Pa": 0.42, "Np": 0.42, "Pu": 0.42, "Am": 0.42, "Cm": 0.42,
            "Bk": 0.42, "Cf": 0.42, "Es": 0.42, "Fm": 0.42, "Md": 0.42,
            "No": 0.42, "Lr": 0.42, "Rf": 0.42, "Db": 0.42, "Sg": 0.42,
            "Bh": 0.42, "Hs": 0.42, "Mt": 0.42, "Ds": 0.42, "Rg": 0.42,
            "Cn": 0.42, "Nh": 0.42, "Fl": 0.42, "Mc": 0.42, "Lv": 0.42,
            "Ts": 0.42, "Og": 0.42, "At": 0.36, "Se": 0.34, "Br": 0.34,
            "Kr": 0.34, "Te": 0.36, "I": 0.36, "Xe": 0.36, "Ga": 0.34,
            "Ge": 0.34, "As": 0.34, "Sb": 0.36, "Bi": 0.38, "Sn": 0.38,
            "In": 0.36, "Cd": 0.36, "Hg": 0.36, "Tl": 0.38, "Ta": 0.38,
            "Nb": 0.36, "Hf": 0.38, "Re": 0.36, "Os": 0.36, "Ir": 0.36,
            "Pt": 0.38, "Pd": 0.36, "Ru": 0.36, "Rh": 0.36, "Tc": 0.36}
ELEM_Z = {"H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
          "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
          "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22,
          "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28,
          "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34,
          "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40,
          "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46,
          "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52,
          "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58,
          "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64,
          "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70,
          "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76,
          "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82,
          "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88,
          "Ac": 89, "Th": 90, "Pa": 91, "U": 92, "Np": 93, "Pu": 94,
          "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100,
          "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105, "Sg": 106,
          "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111, "Cn": 112}


def _f(s) -> float:
    """float() tolerant of COD uncertainty notation '3.5786(5)'."""
    if s is None:
        return 0.0
    s = str(s).strip().strip("'\"")
    if "(" in s:
        s = s.split("(", 1)[0]
    try:
        return float(s)
    except ValueError:
        return 0.0


def cif_calc_lines(cif_path: str, dmin: float = 1.0, dmax: float = 22.0,
                   hkl_cap: int = HKL_CAP) -> tuple:
    """Kinematic d-spacing + relative intensity list from a CIF (no GSAS-II).

    Uses gaussian-atom form factors f_j ~ Z_j exp(-2 pi^2 s^2 w_j^2) and the
    full site sum F(h) = sum_j f_j exp(2 pi i h.x_j) over ALL sites, which
    reproduces systematic absences exactly. Returns (d_units int16, I float64)
    normalized to max I = 1.
    """
    text = Path(cif_path).read_text(errors="replace")
    lines = text.splitlines()

    def tagval(tag: str):
        for ln in lines:
            s = ln.strip()
            if s.startswith(tag):
                v = s[len(tag):].strip().strip("'\"")
                if v and not v.startswith(";"):
                    return v
        return None

    a, b, c = (_f(tagval("_cell_length_a")), _f(tagval("_cell_length_b")),
               _f(tagval("_cell_length_c")))
    al = _f(tagval("_cell_angle_alpha"))
    be = _f(tagval("_cell_angle_beta"))
    ga = _f(tagval("_cell_angle_gamma"))
    if a <= 0 or b <= 0 or c <= 0:
        raise ValueError("no cell in CIF")
    # _atom_site loop: collect column indices then the data rows
    def loop_spec():
        for i, ln in enumerate(lines):
            if ln.strip().startswith("_atom_site_fract_x"):
                j = i
                while j > 0 and (lines[j - 1].strip().startswith("_atom_site")
                                 or lines[j - 1].strip() == "loop_"):
                    j -= 1
                if lines[j].strip() == "loop_":
                    j += 1
                tags = []
                jj = j
                while jj < len(lines) and lines[jj].strip().startswith("_atom_site"):
                    tags.append(lines[jj].strip())
                    jj += 1
                return jj, tags
        return None, None
    row0, tags = loop_spec()
    if row0 is None:
        raise ValueError("no atom sites in CIF")
    idx = {t: tags.index(t) for t in tags if t in (
        "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z",
        "_atom_site_type_symbol", "_atom_site_label",
        "_atom_site_occupancy")}
    try:
        fx = idx["_atom_site_fract_x"]
    except KeyError:
        raise ValueError("no atom sites in CIF")
    fy = idx.get("_atom_site_fract_y", fx)
    fz = idx.get("_atom_site_fract_z", fx)
    elcol = idx.get("_atom_site_type_symbol", idx.get("_atom_site_label"))
    occ_col = idx.get("_atom_site_occupancy", -1)
    sites = []
    for ln in lines[row0:]:
        s = ln.strip()
        if not s or s.startswith("_") or s.startswith(";") or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) <= max(fx, fy, fz):
            continue
        parts = s.split()
        if len(parts) <= max(fx, fy, fz):
            continue
        if occ_col >= len(parts):
            continue
        if elcol >= len(parts):
            continue
        name = parts[elcol]
        sym = "".join(ch for ch in name if ch.isalpha())
        if sym not in ELEM_Z:
            continue
        x, y, z = _f(parts[fx]), _f(parts[fy]), _f(parts[fz])
        if "?" in (parts[fx] + parts[fy] + parts[fz]):
            continue                              # unknown coordinate
        if max(abs(x), abs(y), abs(z)) > 2.0:
            continue
        o = _f(parts[occ_col]) if occ_col >= 0 and occ_col < len(parts) else 1.0
        sites.append((o * ELEM_Z[sym], _GAUSS_W.get(sym, 0.36),
                      (x, y, z)))
    if not sites:
        raise ValueError("no atom sites in CIF")
    # space group: expand the asymmetric unit over the orbit (systematic
    # absences + centring come out exactly this way)
    import gemmi
    itn = tagval("_space_group_IT_number")
    sg = None
    if itn:
        try:
            sg = gemmi.find_spacegroup_by_number(int(float(itn)))
        except Exception:
            sg = None
    if sg is None:
        hm_name = tagval("_symmetry_space_group_name_H-M")
        if hm_name:
            try:
                sg = gemmi.find_spacegroup_by_name(hm_name)
            except Exception:
                sg = None
    orbit_ops = []
    if sg is not None:
        for op in sg.operations():
            rot = np.asarray(op.rot, dtype=np.float64) / 24.0
            tran = np.asarray(op.tran, dtype=np.float64) / 24.0
            orbit_ops.append((rot, tran))
    else:
        orbit_ops = [(np.eye(3), np.zeros(3))]   # no SG info: P1 fallback
    site_list = []
    for zj, wj, (x, y, z) in sites:
        for rot, tran in orbit_ops:
            xp = rot @ np.array([x, y, z]) + tran
            site_list.append((zj, wj, xp))
    n_terms = len(site_list)
    if n_terms > 40000:
        # keep memory bounded: subsample is NOT allowed, so raise honestly
        raise ValueError(f"too many site terms ({n_terms}) for one CIF")
    hmax = min(hkl_cap, max(1, int(np.ceil(2 * a / dmin))))
    kmax = min(hkl_cap, max(1, int(np.ceil(2 * b / dmin))))
    lmax = min(hkl_cap, max(1, int(np.ceil(2 * c / dmin))))
    H, K, L = np.meshgrid(np.arange(-hmax, hmax + 1),
                          np.arange(-kmax, kmax + 1),
                          np.arange(-lmax, lmax + 1), indexing="ij")
    H, K, L = H.ravel(), K.ravel(), L.ravel()
    ca, cb, cg = (np.cos(np.radians(v)) for v in (al, be, ga))
    sa, sb, sg = (np.sin(np.radians(v)) for v in (al, be, ga))
    vc = a * b * c * np.sqrt(1 - ca ** 2 - cb ** 2 - cg ** 2
                             + 2 * ca * cb * cg)
    s11 = (b * c * sa) ** 2
    s22 = (a * c * sb) ** 2
    s33 = (a * b * sg) ** 2
    s12 = a * b * c ** 2 * (ca * cb - cg)
    s23 = a ** 2 * b * c * (cb * cg - ca)
    s13 = a * b ** 2 * c * (ca * cg - cb)
    with np.errstate(divide="ignore", invalid="ignore"):
        d2 = vc ** 2 / (s11 * H ** 2 + s22 * K ** 2 + s33 * L ** 2
                        + 2 * s12 * H * K + 2 * s23 * K * L + 2 * s13 * H * L)
        d = np.sqrt(d2)
    mask = (d >= dmin) & (d <= dmax) & np.isfinite(d)
    H, K, L, d = H[mask], K[mask], L[mask], d[mask]
    hm = np.stack([H, K, L], axis=-1).astype(np.float64)
    s2 = 1.0 / (4.0 * d * d)                    # s^2 = (1/2d)^2
    F = np.zeros(hm.shape[0], dtype=np.complex128)
    for zj, wj, xp in site_list:
        phase = np.exp(2j * np.pi * (hm @ xp))
        f = zj * np.exp(-2 * np.pi ** 2 * s2 * wj * wj)
        F += f * phase
    I = np.abs(F) ** 2
    order = np.argsort(d)
    d = d[order]
    I = I[order]
    du = np.round(d * D_UNIT).astype(np.int64)
    # merge identical d-units (average I)
    u, inv = np.unique(du, return_inverse=True)
    Isum = np.zeros(len(u))
    np.add.at(Isum, inv, I)
    du = u.astype(np.int16)
    I = Isum
    m = I.max()
    if m > 0:
        I = I / m
    return du, I


def intensity_match(measured_fp, du: np.ndarray, I: np.ndarray,
                    tol: float = PEAK_D_TOL) -> dict:
    """Compare measured fingerprint peaks against calculated CIF lines.

    Returns dict: coverage (fraction of measured peaks matched), corr
    (cosine similarity of the measured heights vs calculated intensities at
    the measured peak positions), matched list.
    """
    d = du / D_UNIT
    tol_units = int(np.ceil(tol * D_UNIT)) + 1
    v_m = np.array([p.height for p in measured_fp.peaks], dtype=np.float64)
    v_c = np.zeros(len(measured_fp.peaks))
    n_match = 0
    for i, p in enumerate(measured_fp.peaks):
        d0u = _dunit(p.d)
        lo = np.searchsorted(du, d0u - tol_units, side="left")
        hi = np.searchsorted(du, d0u + tol_units, side="right")
        if lo < hi:
            seg = I[lo:hi]
            v_c[i] = seg.max()
            n_match += 1
    norm = np.linalg.norm(v_m) * np.linalg.norm(v_c) or 1e-30
    corr = float(v_m @ v_c / norm) if v_c.any() else 0.0
    n = len(measured_fp.peaks) or 1
    return {"coverage": round(n_match / n, 4), "corr": round(corr, 4),
            "matched": int(n_match), "n_peaks": n}
