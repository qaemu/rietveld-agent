"""COD wide-open search: rank a measured fingerprint against the COMPLETE COD.

The COD's own metadata CSV export endpoint is size-capped (it silently
truncates around ~100-220k rows), so a *complete* index cannot be built from
it. Instead this module builds the index from the full rsync'd CIF tree
(``data/cod_index/cifs/``, the COD's sanctioned bulk-download method:
``rsync -a rsync://www.crystallography.net/cif/ data/cod_index/cifs/``),
parsing cell + space group + identity metadata directly out of every CIF.
The tree is gitignored (``data/cod_index/``) and rebuilt deterministically
via ``make cod-index`` (``python -m core.codsearch build-index``).

Screening reuses ``core.codindex``: a lattice-geometry line index
(d-spacings of the systematically-allowed hkl reflections, d in [1.1, 22] A,
0.01 A units) with significance-ranked coincidence scoring. The top-K are
then re-scored offline against their LOCAL CIFs with kinematic intensities
(``cif_calc_lines``/``intensity_match``) - no per-request downloads.

Honest limitations (unchanged from unit 12):
  * the index is POSITION-ONLY (lattice geometry, no intensities): it is a
    candidate FILTER, not an identification;
  * on multi-phase patterns the screen dilutes (every phase competes with
    the full database); single-phase patterns with distinctive lattices
    rank cleanly;
  * significance is a coincidence score, not a calibrated p-value;
  * huge-cell entries have their index truncated at hkl cap 28, recorded
    per entry in ``dmin_eff``.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from multiprocessing import get_all_start_methods, get_context
from pathlib import Path
from typing import List, Optional

import numpy as np

from core._paths import repo_root
from core.codindex import (D_MAX, D_MIN, CodEntry, build_index,
                           cif_calc_lines, intensity_match,
                           screen_fingerprint)

ROOT = repo_root()
CIF_ROOT = ROOT / "data" / "cod_index" / "cifs"
NPZ_PATH = ROOT / "data" / "cod_index" / "cod_entire_v1.npz"
META_PATH = ROOT / "data" / "cod_index" / "cod_entire_v1.meta.json"
PATHS_PATH = ROOT / "data" / "cod_index" / "cod_entire_v1.paths.npy"
COD_TOTAL_AT_BUILD = "534674"      # COD counter, 2026-08-10 (site front page)
RERANK_MAX_LINES = 2500            # cap for the offline intensity rerank:
                                   # giant-cell organics are skipped (ranked
                                   # by geometry only) to bound the cost

_CACHE = None                       # (cod_ids, d_units, entry_of, dmin_eff,
                                    #  metas_by_id, paths_by_id)


# --------------------------------------------------------------------------- #
# fast CIF scanner (cell + SG + identity only; no atom sites needed)         #
# --------------------------------------------------------------------------- #
def _first_value(lines: List[str], tag: str) -> str:
    """Value of the first ``_tag`` line (single-line form only)."""
    want = "_" + tag
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith(want):
            rest = s[len(want):].strip().strip("'\"")
            if rest and not rest.startswith(";"):
                return rest
            if rest.startswith(";"):        # quoted text block
                out = [rest[1:]]
                for nxt in lines[i + 1:]:
                    if nxt.strip().startswith(";"):
                        return " ".join(out).strip()
                    out.append(nxt.strip())
    return ""


def _loop_first_column(lines: List[str], tag: str) -> str:
    """First column of a ``loop_`` block whose header contains ``tag``."""
    for i, ln in enumerate(lines):
        if ln.strip().startswith("_" + tag):
            j = i
            while j > 0 and (lines[j - 1].strip().startswith("_")
                             or lines[j - 1].strip() == "loop_"):
                j -= 1
            if lines[j].strip() == "loop_":
                j += 1
            jj = j
            while jj < len(lines) and lines[jj].strip().startswith("_"):
                jj += 1
            vals = []
            for row in lines[jj:]:
                s = row.strip()
                if not s or s.startswith("_") or s.startswith(";"):
                    break
                vals.append(s.split()[0])
            return "; ".join(vals)
    return ""


def _strip_esd(value: str) -> str:
    """'3.5786(5)' -> '3.5786' (COD uncertainty notation)."""
    v = value.strip().strip("'\"")
    return v.split("(", 1)[0].strip() if "(" in v else v


def _row_from_block(block: List[str], cod_id: int) -> dict:
    """Row dict with the keys CodEntry expects."""
    def t(tag: str) -> str:
        return _first_value(block, tag)

    row = {
        "file": str(cod_id),
        "sgNumber": _strip_esd(t("space_group_IT_number")) or _strip_esd(
            t("space_group_IT_number_alt")) or "",
        "a": _strip_esd(t("cell_length_a")),
        "b": _strip_esd(t("cell_length_b")),
        "c": _strip_esd(t("cell_length_c")),
        "alpha": _strip_esd(t("cell_angle_alpha")),
        "beta": _strip_esd(t("cell_angle_beta")),
        "gamma": _strip_esd(t("cell_angle_gamma")),
        "mineral": t("chemical_name_mineral"),
        "chemname": (t("chemical_name_common") or t("chemical_name_systematic")),
        "formula": t("chemical_formula_sum"),
        "year": t("journal_year"),
        "doi": t("database_code_doi"),
        "journal": t("journal_name_full"),
        "authors": _loop_first_column(block, "publ_author_name"),
        "title": t("publ_section_title"),
    }
    # SG fallback: H-M name -> number (via gemmi) when the IT number is
    # missing (rare in COD; keeps coverage maximal).
    if not row["sgNumber"]:
        hm = (t("space_group_name_H-M_alt") or t("space_group_name_H-M")
              or t("symmetry_space_group_name_H-M") or t("space_group_name_Hall"))
        if hm:
            try:
                import gemmi
                sg = gemmi.find_spacegroup_by_name(hm)
                if sg is not None:
                    row["sgNumber"] = str(sg.number)
            except Exception:            # noqa: BLE001 -- best effort only
                pass
    return row


def _scan_one(path: Path) -> Optional[tuple]:
    """One CIF file -> (cod_id, row dict) for every parseable block, or None."""
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return None
    cod_id = int(path.stem)
    out = []
    blocks, cur = [], []
    for ln in text.splitlines():
        if ln.startswith("data_"):
            if cur:
                blocks.append(cur)
            cur = [ln]
        else:
            cur.append(ln)
    if cur:
        blocks.append(cur)
    if not blocks:
        return None
    for b in blocks:
        row = _row_from_block(b, cod_id)
        e = CodEntry(row)
        if e.valid:
            out.append((cod_id, row))
    return out or None


def _pool_ctx():
    """fork when available (macOS/Linux; fast and stdin-safe), spawn fallback
    (Windows, where fork does not exist)."""
    if "fork" in get_all_start_methods():
        return get_context("fork")
    return get_context()


def scan_cif_tree(cif_root: Path = CIF_ROOT, workers: int = 8) -> List[CodEntry]:
    """Walk the rsync'd CIF tree and parse every entry (cell + identity).

    Returns entries sorted by cod_id (deterministic build order). Entries
    without a usable cell + space group are skipped (counted in the caller).
    """
    files = sorted(cif_root.rglob("*.cif"))
    entries, files_no_entry, n_files = [], 0, len(files)
    with _pool_ctx().Pool(workers) as pool:
        for res in pool.imap_unordered(_scan_one, files, chunksize=64):
            if res is None:
                files_no_entry += 1
                continue
            for _cod_id, row in res:
                entries.append(CodEntry(row))
    entries.sort(key=lambda e: e.cod_id)
    return entries, files_no_entry, n_files


def build_entire_index(cif_root: Path = CIF_ROOT, workers: int = 8) -> dict:
    """Scan the full CIF tree and rebuild the complete COD line index.

    Outputs (all under data/cod_index/, gitignored):
      cod_entire_v1.npz        the d-unit line index (codindex layout)
      cod_entire_v1.meta.json  per-entry identity + build manifest (md5)
      cod_entire_v1.paths.npy  cod_id -> relative cif path (offline CIFs)

    The build is deterministic (entries sorted by cod_id; bucket-scatter
    behind the index is order-stable), and the manifest md5 lets an
    interchangeable host verify the payload (repository invariant 5).
    """
    t0 = time.time()
    print(f"[cod] scanning CIF tree {cif_root} ...", flush=True)
    entries, n_skipped, n_files = scan_cif_tree(cif_root, workers=workers)
    print(f"[cod] {len(entries)} entries with valid cell+SG "
          f"({n_files - len(entries)} skipped), {time.time() - t0:.0f}s",
          flush=True)
    if not entries:
        raise SystemExit("no entries parsed - is the CIF tree rsync'd?")
    print(f"[cod] building reflection index ({workers} workers) ...", flush=True)
    cod_ids, d_units, entry_of, dmin_eff = build_index(entries)
    out = ROOT / "data" / "cod_index"
    # uncompressed on purpose: npz decompression of ~35 MB on every screen
    # costs tens of seconds per fresh process; raw saves load in ~1 s.
    np.savez(NPZ_PATH, cod_ids=cod_ids, d_units=d_units,
             entry_of=entry_of, dmin_eff=dmin_eff)
    # deterministic path map
    paths = np.empty(len(entries), dtype="<U48")
    rel = {p.stem: p.relative_to(cif_root).as_posix()
           for p in cif_root.rglob("*.cif")}
    for i, e in enumerate(entries):
        paths[i] = rel.get(str(e.cod_id), "")
    np.save(PATHS_PATH, paths)
    metas = [e.meta_dict() for e in entries]
    manifest = {
        "manifest": "cod_entire_v1",
        "entries": len(entries),
        "cif_files": n_files,
        "skipped_no_cell_sg": n_files - len(entries),
        "cod_total_front_page": COD_TOTAL_AT_BUILD,
        "cif_tree": "data/cod_index/cifs (rsync of rsync://www.crystallography.net/cif/)",
        "d_range_A": [D_MIN, D_MAX],
        "hkl_cap": 28,
        "built_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "index_md5": hashlib.md5(
            (NPZ_PATH.read_bytes() + PATHS_PATH.read_bytes())).hexdigest(),
    }
    META_PATH.write_text(json.dumps({"manifest": manifest, "entries": metas}))
    print(f"[cod] wrote {NPZ_PATH.name} ({len(entries)} entries, "
          f"{d_units.size} lines, {time.time() - t0:.0f}s)", flush=True)
    return manifest


def load_entire_index():
    """Load the complete-COD index (cached in-process)."""
    global _CACHE
    if _CACHE is not None:
        return _CACHE
    if not NPZ_PATH.exists():
        raise FileNotFoundError(
            f"complete-COD index not built ({NPZ_PATH}); run "
            f"`make cod-index` (requires the rsync'd CIF tree)")
    z = np.load(NPZ_PATH)
    datas = json.loads(META_PATH.read_text())
    metas = {m["cod_id"]: m for m in datas["entries"]}
    paths = np.load(PATHS_PATH) if PATHS_PATH.exists() else None
    paths_by_id = {int(z["cod_ids"][i]): str(paths[i]) for i in
                   range(len(z["cod_ids"]))} if paths is not None else {}
    _CACHE = (z["cod_ids"], z["d_units"], z["entry_of"], z["dmin_eff"],
              metas, paths_by_id)
    return _CACHE


def cif_path_for(cod_id: int) -> Optional[Path]:
    _, _, _, _, _, paths_by_id = load_entire_index()
    rel = paths_by_id.get(cod_id)
    return CIF_ROOT / rel if rel else None


def screen_cod(fp, top_k: int = 8, rerank: bool = True,
               rerank_top: int = 200) -> dict:
    """Screen a measured fingerprint against the COMPLETE COD.

    Stage 1 (geometry): ``screen_fingerprint`` builds a candidate pool --
    union of the top entries by significance, by intensity-weighted match
    and by peak coverage. Stage 2 (rerank): the pool is cross-checked
    against the local CIF with offline kinematic intensities
    (``intensity_match``); the final order is by intensity agreement
    (coverage, then corr), which on real patterns beats pure geometry
    significance (dense-line high-symmetry entries can cover many weak
    peaks by chance while missing the strong ones).

    Returns bundle-ready ``cod_screen`` evidence: index size, top-K hits
    (intensity-ranked when rerank is on, geometry-ranked otherwise), and
    the per-hit ``intensity_rerank`` cross-check (null when the local CIF
    is unavailable).
    """
    cod_ids, d_units, entry_of, dmin_eff, metas, _ = load_entire_index()
    pool = screen_fingerprint(fp, cod_ids, d_units, entry_of, metas,
                              top_k=max(top_k, rerank_top))
    if rerank:
        for hit in pool:
            try:
                path = cif_path_for(hit["cod_id"])
                if path is None or not path.exists():
                    hit["intensity_rerank"] = None
                    continue
                du, I = cif_calc_lines(str(path))
                if len(du) > RERANK_MAX_LINES:
                    # giant-cell organics: too many lines to recompute
                    # tractably, and their intensity agreement is
                    # uninformative; rank them by geometry only.
                    hit["intensity_rerank"] = None
                    continue
                hit["intensity_rerank"] = intensity_match(fp, du, I)
            except Exception as exc:               # noqa: BLE001
                hit["intensity_rerank"] = {"error": str(exc)[:200]}
        # final order: entries with a usable rerank first (coverage, then
        # corr, then geometry significance), CIF-less entries fall back to
        # their geometry pool order.
        def key(h):
            rr = h.get("intensity_rerank") or {}
            if not rr or "error" in rr:
                return (0, 0.0, 0.0, -h["significance"])
            return (1, rr.get("coverage", 0.0), rr.get("corr", 0.0),
                    -h["significance"])
        pool.sort(key=key, reverse=True)
    manifest = json.loads(META_PATH.read_text())["manifest"]
    return {
        "index_entries": int(len(cod_ids)),
        "index_manifest": manifest.get("manifest", ""),
        "index_md5": manifest.get("index_md5", ""),
        "d_range_A": manifest.get("d_range_A", [D_MIN, D_MAX]),
        "top": pool[:top_k],
    }


def _demo_fingerprint(path: str):
    """Standalone screen helper: parse + fingerprint an XRDML measurement."""
    from core.ingest import parse_xrdml, sample_fingerprint
    return sample_fingerprint(parse_xrdml(path))


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(prog="core.codsearch")
    sub = ap.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build-index", help="rebuild the complete-COD index "
                                           "from the rsync'd CIF tree")
    b.add_argument("--cifs", default=str(CIF_ROOT))
    b.add_argument("--workers", type=int, default=8)
    s = sub.add_parser("screen", help="screen one measurement vs the full COD")
    s.add_argument("xrdml")
    s.add_argument("--top", type=int, default=8)
    s.add_argument("--no-rerank", action="store_true")
    args = ap.parse_args(argv)
    if args.command == "build-index":
        build_entire_index(Path(args.cifs), workers=args.workers)
        return 0
    res = screen_cod(_demo_fingerprint(args.xrdml), top_k=args.top,
                     rerank=not args.no_rerank)
    print(f"index entries : {res['index_entries']}")
    print(f"index md5     : {res['index_md5'][:16]}...")
    for i, hit in enumerate(res["top"], 1):
        rr = hit.get("intensity_rerank") or {}
        print(f"#{i} COD {hit['cod_id']:<8} {(hit.get('formula') or '')[:42]:42} "
              f"sig={hit['significance']:>6} "
              f"corr={rr.get('corr', '-'):>6}")
    return 0


if __name__ == "__main__":
    sys.exit(main())