# Spike 18: complete-COD screening (core.codsearch)

Goal (from the spike-12 retrofit): the position-only lattice-geometry
screen (peak window +-0.02 A) must run against the *entire* COD, not the
size-capped REST/CSV export (~91k of 534,674 rows on fresh requests).

## Approach

- Mirror the full CIF tree with rsync:
  `rsync -a --partial rsync://www.crystallography.net/cif/ data/cod_index/cifs/`
  (~26 GB, 534k files; resumable).
- `core/codsearch.scan_cif_tree` parses every CIF (gemmi), keeps entries
  with a cell + a resolvable space group (SG number or Hall -> number),
  and extracts formula / mineral / chemname / year / doi / journal /
  title from the CIF header block.
- `core/codindex.build_index` builds the bucket-scattered int16 d-unit
  line index (hkl cap 28, d in [0.8, 18] A) + per-entry dmin_eff; the
  sort-free bucket layout keeps the build O(n) and order-stable
  (deterministic across hosts; manifest md5).
- Screening: sorted-search `screen_fingerprint` (counts matched lines in
  each window) -> `screen_cod` (top-k + significance incl. stochastic
  congruent-group size term), then offline kinematic-intensity rerank
  against the local CIF (`intensity_match`) for the top candidates.

## Rough scale (measured on a partial tree, 8 workers)

- ~45k CIFs: scan 20 s; index build 153 s (0.0034 s/entry) -> the full
  index (534k) should scan in ~5 min and build in ~31 min (8 workers).
- Skip rate w/o cell/SG ~0.8% (mostly odd organic/file captures).
- ~122 lines/entry -> ~65M lines at full scale -> ~0.6-1 GB npz (int16
  + int32 (+journal text kept as JSON), gitignored).

## CLI / bundle wiring

- `cli analyze --full-cod` -> `cod_screen` block in the run bundle
  (index_entries, index_manifest, index_md5, d_range_A, top[] with
  significance + intensity_rerank). Schema: run_bundle.schema.json.
- `make cod-tree` / `make cod-index` (entry points for the mirror+build).
- Tests: tests/test_cod_full.py (skip when index absent): index >= 500k
  entries, known structures present (anglesite 1010950, alite M3
  9008366, belite 9012794, corundum 1000017, ferrite 1200009, aluminate
  8103596), anglesite rank 1 on the Cu and Fe fixtures, rerank agrees,
  --full-cod bundle validates against the schema.

## Notes / failure modes seen

- multiprocessing default on macOS is spawn: a scanner/index test invoked
  without a `__main__` guard re-imports the driver script in workers and
  recursively re-spawns pools (BrokenProcessPool). Fixed by using the
  fork start method when available (`core.codsearch._pool_ctx`, also
  passed as `mp_context` to the index ProcessPoolExecutor).
- An rsync retry loop from an earlier triple-invocation collided three
  receivers on the same tree (transfer truncated, code 12); the single-
  process restart with `--partial` is stable.
- Do not run the index build while the tree is still syncing (workers
  read half-written CIFs; wait for the mirror to finish).
- **Geometry significance is a pre-filter, not the final ranking.** On the
  partial index, dense-line high-symmetry entries (13-line FCC Sn-Tl
  alloy) out-ranked the true phase: their lines sit inside many of the
  pattern's crowded *weak* peaks (m=0.41) while missing the strong ones,
  and the coincidence model rewards that overlap. The true phase wins on
  intensity-weighted coverage (m=0.995, coverage=1.0). `screen_cod` now
  builds a candidate pool as the union of the top entries by
  significance / matched intensity / coverage, then ranks it by the
  offline kinematic-intensity rerank (coverage, then corr) -- matching
  the spike-12 architecture (significant filter -> intensity ranking ->
  (optional) GSAS-II confirmation). Giant-cell organics (> 2500 lines)
  are ranked by geometry only.
- Index is stored uncompressed (npz): compressed saves cost tens of
  seconds of decompression on every fresh screen process; raw ~40 MB
  loads in ~1 s (16 s cold, page-cache dependent).
- Screen + rerank cost on 53k entries, 8 workers: load ~2-16 s,
  geometry screen ~3 s, 120-CIF rerank ~4 s -> ~10-25 s per sample;
  linear-ish in pool size, not index size (npz searchsorted is O(log n)).

## Reproduce

```
make cod-tree    # ~2-3 h, resumable
make cod-index   # ~1 h
python3 -m pytest tests/test_cod_full.py -q
python3 -m cli analyze tests/fixtures/xrdml/cu_PbSO4.xrdml --full-cod
```