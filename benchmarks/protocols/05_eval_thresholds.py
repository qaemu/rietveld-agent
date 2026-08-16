"""Unit 05: threshold-robustness evaluation for the M1 fingerprint thresholds.

Question
--------
Do the M1 constants (min_top_similarity = 0.35, min_margin = 0.10,
core/verdict/verdict.py) hold for laboratory-grade data, and where do they
break?

Method
------
For each material in the M1 library (PbSO4 Cu+Fe, SiO2, NaCl) generate noisy
realizations of the clean fixture at 8 cumulative severity levels:
    L0 clean                          counts=None   s=0.00mm  bg=0.00
    L1 1e6 peak counts, 0.02 mm,  0.5% bg
    L2 3e5 peak counts, 0.05 mm,  1%   bg
    L3 1e5 peak counts, 0.10 mm,  2%   bg
    L4 3e4 peak counts, 0.20 mm,  5%   bg
    L5 1e4 peak counts, 0.30 mm, 10%   bg
    L6 3e3 peak counts, 0.50 mm, 15%   bg
    L7 1e3 peak counts, 0.80 mm, 20%   bg
(counts = Poisson peak-max counts; s = sample displacement; bg = background
drift fraction), 50 seeded replicates each. Every realization flows through
the real pipeline: ingest -> fingerprint -> rank vs the M1 library -> verdict.

Controls
--------
* An amorphous (featureless hump) synthetic pattern at every level: supported
  verdicts are false positives and must be 0 across all severity.
* L0 (clean) must reproduce the M1 e2e results exactly.

Output
------
data/unit05/results/unit05_report.{json,md} with per-material stats
(median/p10 similarity, margin, abstain rate, flip level) and a verdict on
whether 0.35 / 0.10 survive, plus the statistics gate the data supports.

Run:  python benchmarks/protocols/05_eval_thresholds.py
"""
from __future__ import annotations

import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES_DIR = os.path.join(ROOT, "data", "unit05", "results")
FIX_DIR = os.path.join(ROOT, "tests", "fixtures", "xrdml")
REGISTRY = os.path.join(ROOT, "data", "unit03", "results", "registry.json")
LIBRARY = os.path.join(ROOT, "data", "candidates", "library.json")

sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from benchmarks.eval.noise import amorphous_pattern, perturb  # noqa: E402
from core.calibration import CalibrationRegistry  # noqa: E402
from core.hypothesis import load_library, rank_candidates  # noqa: E402
from core.ingest import InstrumentParams, PowderPattern, parse_xrdml, sample_fingerprint  # noqa: E402
from core.verdict import decide  # noqa: E402

LEVELS = [
    {"id": "L0", "counts": None,      "s_mm": 0.00, "bg_fraction": 0.00},
    {"id": "L1", "counts": 1_000_000, "s_mm": 0.02, "bg_fraction": 0.005},
    {"id": "L2", "counts": 300_000,   "s_mm": 0.05, "bg_fraction": 0.01},
    {"id": "L3", "counts": 100_000,   "s_mm": 0.10, "bg_fraction": 0.02},
    {"id": "L4", "counts": 30_000,    "s_mm": 0.20, "bg_fraction": 0.05},
    {"id": "L5", "counts": 10_000,    "s_mm": 0.30, "bg_fraction": 0.10},
    {"id": "L6", "counts": 3_000,     "s_mm": 0.50, "bg_fraction": 0.15},
    {"id": "L7", "counts": 1_000,     "s_mm": 0.80, "bg_fraction": 0.20},
]
REPLICATES = 50
SEED_BASE = 1000


def _halite_pattern(work_dir: str) -> PowderPattern:
    """Clean real-halite (NaCl, COD 1000041) pattern via GSAS-II.

    Unit 06 chemistry validation showed the old 'halite_9011025.cif'
    (used by unit 04/05) is Ag0.5Bi0.5S, so the eval material is the
    genuine NaCl structure (1000041)."""
    from benchmarks.eval import sim as _sim
    assert _sim is not None
    vendor = os.path.join(ROOT, ".vendor", "GSAS-II")
    prm = _sim.ensure_gsasii(ROOT, vendor, "")
    cif = os.path.join(ROOT, "data", "unit04", "input", "halite_1000041.cif")
    return _sim.sim_cif_to_pattern(cif, work_dir, prm_path=prm)


def _stats(values) -> dict:
    a = np.asarray(values, dtype=float)
    return {"median": float(np.median(a)), "p10": float(np.percentile(a, 10)),
            "p90": float(np.percentile(a, 90)), "min": float(np.min(a))}


def _run_material(clean: PowderPattern, library, names, families, resolution,
                  label: str, min_peak_counts: float) -> dict:
    rows = []
    for li, level in enumerate(LEVELS):
        sims, margins, fam_margins, abstains, primaries = [], [], [], 0, {}
        for r in range(REPLICATES):
            noisy = perturb(clean, counts=level["counts"], s_mm=level["s_mm"],
                            bg_fraction=level["bg_fraction"],
                            seed=SEED_BASE + r)
            ranking = rank_candidates(sample_fingerprint(noisy), library,
                                      names=names, families=families)
            verdict = decide(ranking, resolution, peak_max=noisy.peak_max,
                             min_peak_counts=min_peak_counts)
            sims.append(ranking.top_similarity)
            margins.append(ranking.margin)
            fam_margins.append(ranking.family_margin)
            if verdict.status != "supported":
                abstains += 1
            elif verdict.primary:
                primaries[verdict.primary.material_id] = \
                    primaries.get(verdict.primary.material_id, 0) + 1
        row = {"level": level["id"], "counts": level["counts"],
               "s_mm": level["s_mm"], "bg_fraction": level["bg_fraction"],
               "sim": _stats(sims), "margin": _stats(margins),
               "family_margin": _stats(fam_margins),
               "n": len(sims), "abstain_rate": abstains / len(sims),
               "primary_mode": max(primaries, key=primaries.get)
               if primaries else None}
        rows.append(row)
    # flip level: first level where at least half the realizations abstain
    flip = next((r["level"] for r in rows if r["abstain_rate"] >= 0.5), None)
    return {"material": label, "rows": rows, "flip_level": flip,
            "clean_supported": rows[0]["abstain_rate"] == 0.0}


def _run_amorphous(tth, cu_params, library, names, families, reg,
                   min_peak_counts: float) -> list:
    out = []
    for level in LEVELS:
        abstains = 0
        for r in range(REPLICATES):
            pat = amorphous_pattern(tth, counts=level["counts"], s_mm=level["s_mm"],
                                    bg_fraction=level["bg_fraction"],
                                    seed=SEED_BASE + r, instrument=cu_params)
            res = reg.lookup(pat.instrument)
            ranking = rank_candidates(sample_fingerprint(pat), library,
                                      names=names, families=families)
            verdict = decide(ranking, res, peak_max=pat.peak_max,
                             min_peak_counts=min_peak_counts)
            if verdict.status != "supported":
                abstains += 1
        out.append({"level": level["id"], "abstain_rate": abstains / REPLICATES,
                    "fp_rate": 1.0 - abstains / REPLICATES})
    return out


def main() -> None:
    t0 = time.time()
    os.makedirs(RES_DIR, exist_ok=True)
    work = os.path.join(ROOT, "data", "unit05", "work")
    os.makedirs(work, exist_ok=True)

    lib_payload = json.load(open(LIBRARY))
    library = load_library(lib_payload["materials"])
    names = {m["id"]: m["name"] for m in lib_payload["materials"]}
    families = {m["id"]: m["phase_family"] for m in lib_payload["materials"]}
    reg = CalibrationRegistry.load(REGISTRY)

    fixtures = {
        "pbso4-cu":   (os.path.join(FIX_DIR, "cu_PbSO4.xrdml"), None),
        "sio2-cu":    (os.path.join(FIX_DIR, "cu_quartz.xrdml"), None),
        "pbso4-fe":   (os.path.join(FIX_DIR, "fe_PbSO4.xrdml"), None),
    }
    clean = {k: parse_xrdml(p) for k, (p, _) in fixtures.items()}
    print("[nacl] GSAS-II halite sim (once)...")
    clean["nacl-cu"] = _halite_pattern(work)
    print(f"[nacl] sim ok, npts={clean['nacl-cu'].npts}")
    cu_params = clean["pbso4-cu"].instrument
    tth = clean["pbso4-cu"].tth

    # --- pass 1: evidence-only (no statistics gate) ------------------------
    ev = {}
    for label, pattern in clean.items():
        resolution = reg.lookup(pattern.instrument)
        if resolution.status.value != "resolved":
            raise SystemExit(f"calibration not resolved for {label}: "
                             f"{resolution.reason}")
        ev[label] = _run_material(pattern, library, names, families,
                                  resolution, label, min_peak_counts=0.0)
    ev_amorph = _run_amorphous(tth, cu_params, library, names, families, reg,
                               min_peak_counts=0.0)
    fp_total = max(a["fp_rate"] for a in ev_amorph)
    print(f"[amorphous] max fp_rate={fp_total:.4f}  (must be 0)")

    # evidence flip (first level with >=50% abstain) per material
    ev_flip = {}
    for label, m in ev.items():
        ev_flip[label] = next((r["level"] for r in m["rows"]
                               if r["abstain_rate"] >= 0.5), None)
        print(f"[evidence {label}] flip={ev_flip[label]} "
              f"L0 abstain={m['rows'][0]['abstain_rate']}")

    # --- derive statistics gate ----------------------------------------------
    # gate = counts of the last severity level where EVERY material is still
    # fully supported (abstain_rate == 0); None when even L0 is not clean.
    last_clean = None
    for i in range(len(LEVELS) - 1, -1, -1):
        if all(ev[l]["rows"][i]["abstain_rate"] == 0.0 for l in ev):
            last_clean = i
            break
    if last_clean is None:
        raise SystemExit("no level is clean for every material; aborting gate derivation")
    gate_counts = LEVELS[last_clean]["counts"]   # None for L0 -> no gate
    print(f"[gate] last fully-supported level = {LEVELS[last_clean]['id']} "
          f"-> min_peak_counts = {gate_counts}")

    # --- pass 2: fused pipeline (evidence + statistics gate) ----------------
    fused = {}
    for label, pattern in clean.items():
        resolution = reg.lookup(pattern.instrument)
        fused[label] = _run_material(pattern, library, names, families,
                                     resolution, label,
                                     min_peak_counts=gate_counts or 0.0)
    fused_amorph = _run_amorphous(tth, cu_params, library, names, families, reg,
                                  min_peak_counts=gate_counts or 0.0)
    fused_fp = max(a["fp_rate"] for a in fused_amorph)
    fused_flip = {label: next((r["level"] for r in m["rows"]
                               if r["abstain_rate"] >= 0.5), None)
                  for label, m in fused.items()}

    # --- slope sanity: abstain rate should be non-decreasing with level ----
    monotone_ok = True
    for label, m in fused.items():
        rates = [r["abstain_rate"] for r in m["rows"]]
        if any(b < a for a, b in zip(rates, rates[1:])):
            monotone_ok = False
            print(f"[warn] {label} fused abstain rate not monotone: {rates}")

    # --- verdict ------------------------------------------------------------
    # thresholds judged on the EVIDENCE pass (the gate is a separate
    # conservative layer protecting decisions below the validated envelope)
    idx_map = {lv["id"]: i for i, lv in enumerate(LEVELS)}
    ev_first = min((idx_map[f] for f in ev_flip.values() if f), default=len(LEVELS))
    critical = [label for label, f in ev_flip.items()
                if f and idx_map[f] == ev_first]
    keep_ok = ev_first >= 4
    if fp_total > 0.0:
        verdict_text = f"FAIL: amorphous false positives at fp_rate={fp_total:.3f}"
        ok = False
    elif ev_first == 0:
        verdict_text = "FAIL: clean (L0) patterns abstain; pipeline bug."
        ok = False
    else:
        fused_first = min((idx_map[f] for f in fused_flip.values() if f),
                          default=len(LEVELS))
        gate_str = (f"Statistics gate min_peak_counts = {gate_counts} "
                    f"(derived at envelope level {LEVELS[last_clean]['id']}); "
                    f"the fused pipeline intentionally abstains from the gate "
                    f"boundary (L{fused_first}) onward."
                    if gate_counts else "No statistics gate needed.")
        verdict_text = (
            f"{'KEEP' if keep_ok else 'ADJUST'} thresholds 0.35/0.10 with "
            f"family-aware margin: evidence remains majority-supported through "
            f"{LEVELS[ev_first - 1]['id']} (first majority flip at "
            f"{LEVELS[ev_first]['id']}; critical: {' / '.join(critical)}). "
            f"{gate_str}")
        ok = True

    # --- report --------------------------------------------------------------
    report = {
        "unit": "unit_05",
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "method": {"levels": LEVELS, "replicates": REPLICATES,
                   "rng": "numpy default_rng, seed_base=1000 (deterministic)",
                   "pipeline": "perturb -> sample_fingerprint -> rank -> decide"},
        "library_manifest": lib_payload["manifest_sha256"][:16],
        "constants_under_test": {"min_top_similarity": 0.35, "min_margin": 0.10,
                                 "margin_kind": "family_margin (vs best other phase family)"},
        "evidence_pass": {"materials": ev, "amorphous_control": ev_amorph,
                          "flip_levels": ev_flip},
        "statistics_gate": {"min_peak_counts": gate_counts,
                            "derived_from_level": LEVELS[last_clean]["id"],
                            "rule": "counts at the last severity level at which "
                                    "every material is fully supported"},
        "fused_pass": {"materials": fused, "amorphous_control": fused_amorph,
                       "flip_levels": fused_flip},
        "monotone_ok": monotone_ok,
        "threshold_verdict": verdict_text,
        "passed": ok,
        "wall_clock_s": round(time.time() - t0, 1),
    }
    with open(os.path.join(RES_DIR, "unit05_report.json"), "w") as fh:
        json.dump(report, fh, indent=2)

    lines = [
        "# Unit 05: threshold-robustness evaluation (M1 fingerprint thresholds)",
        "",
        f"Constants under test: `min_top_similarity = 0.35`, "
        f"`min_margin = 0.10` with the **family-aware margin** "
        f"(top similarity minus the best different-phase-family candidate; "
        f"introduced by unit 05: statuses are per phase family, so "
        f"differently-realized entries of the SAME family must not compete). "
        f"Library manifest `{report['library_manifest']}...`, {REPLICATES} "
        f"seeded replicates per level (deterministic RNG).",
        "",
        "## Noise envelope (cumulative severity)",
        "",
        "| level | peak counts | displacement (mm) | background drift |",
        "|---|---|---|---|",
    ]
    for lv in LEVELS:
        lines.append(f"| {lv['id']} | {lv['counts']} | {lv['s_mm']} | {lv['bg_fraction']} |")
    lines += [
        "",
        "## Pass 1: evidence-only (no statistics gate)",
        "",
        "Abstain rate; median sim / median family margin.",
        "",
    ]
    lines.append("| material | " + " | ".join(l["id"] for l in LEVELS)
                 + " | evidence flip |")
    lines.append("|---|" + "---|" * len(LEVELS) + "---|")
    for label, m in ev.items():
        cells = [f"{r['abstain_rate']:.2f} "
                 f"({r['sim']['median']:.3f}/{r['family_margin']['median']:.3f})"
                 for r in m["rows"]]
        lines.append(f"| {label} | " + " | ".join(cells)
                     + f" | {m['flip_level']} |")
    lines += [
        "",
        "## Statistics gate",
        "",
        f"- Derived rule: counts at the last severity level at which EVERY "
        f"material is fully supported -> "
        f"**min_peak_counts = {gate_counts}** "
        f"(level {LEVELS[last_clean]['id']}).",
        f"- Implemented in `core/verdict` ({'decide(..., min_peak_counts='
        f'{gate_counts})' if gate_counts else 'off'}).",
        "",
        "## Pass 2: fused pipeline (evidence + statistics gate, as shipped)",
        "",
    ]
    lines.append("| material | " + " | ".join(l["id"] for l in LEVELS)
                 + " | fused flip |")
    lines.append("|---|" + "---|" * len(LEVELS) + "---|")
    for label, m in fused.items():
        cells = [f"{r['abstain_rate']:.2f}" for r in m["rows"]]
        lines.append(f"| {label} | " + " | ".join(cells)
                     + f" | {m['flip_level']} |")
    lines += [
        "",
        "## Controls",
        f"- Amorphous negative control (fused): max false-positive rate "
        f"**{fused_fp:.4f}** (must be 0).",
        "- L0 clean must reproduce M1 e2e (all supported): "
        f"{'yes' if all(m['rows'][0]['abstain_rate'] == 0.0 for m in fused.values()) else 'NO'}.",
        f"- Abstain-rate monotonicity (fused): "
        f"{'ok' if monotone_ok else 'VIOLATED'}.",
        "",
        "## Verdict",
        f"- **{verdict_text}**",
        f"- Fused flips: {json.dumps(fused_flip)}",
        f"- Wall clock: {report['wall_clock_s']}s",
        f"- **passed = {report['passed']}**",
    ]
    with open(os.path.join(RES_DIR, "unit05_report.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print()
    print(verdict_text)
    print("wrote", os.path.join(RES_DIR, "unit05_report.json"))
    sys.exit(0 if ok else 2)


if __name__ == "__main__":
    main()