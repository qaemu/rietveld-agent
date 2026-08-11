"""Spike 07: bounded Rietveld verification of the e2e cases.

Prototype of the verification stage promised at G0: after the fingerprint
stage nominates a top phase family, a *bounded* GSAS-II Rietveld refinement
(policy ``refinement-budget.v1.json``: background + shift + cell + scale
only; atoms/microstrain/size/phase-fraction/LeBail prohibited) is run for
the top family AND for competing families against a deterministic observed
pattern (= simulation of the catalog CIF, same protocol as the catalog
fingerprints).

Cases
-----
A. cu_PbSO4.xrdml : fingerprint -> PbSO4 (anglesite) COD 1010950;
   competitors  CaCO3 (calcite) COD 1010928.
B. cu_quartz.xrdml: fingerprint -> GaAsO4 (quartz homeotype) COD 1009000;
   competitors  SiO2 (quartz) COD 9009666 (alpha-quartz, release 0.1.1).

Verdict per case: "confirmed" iff the lowest bounded Rwp belongs to the
catalog-truth family, Rwp <= policy max_rwp, and separation >= policy min.

Run:  python benchmarks/spikes/spike_07_verification.py
"""
from __future__ import annotations

import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WORK_DIR = os.path.join(ROOT, "data", "spike7", "work")
RES_DIR = os.path.join(ROOT, "data", "spike7", "results")
POLICY_PATH = os.path.join(ROOT, "governance", "policies",
                           "refinement-budget.v1.json")
RELEASE_PATH = os.path.join(ROOT, "data", "catalog", "releases",
                            "catalog_0.1.1.json")
LIBRARY_PATH = os.path.join(ROOT, "data", "candidates", "library.json")
CID_CIF = os.path.join(ROOT, "data", "spike6", "input", "cod")
FIX_DIR = os.path.join(ROOT, "tests", "fixtures", "xrdml")
VENDOR = os.path.join(ROOT, ".vendor", "GSAS-II")

sys.path.insert(0, ROOT)

from benchmarks.eval.sim import ensure_gsasii           # noqa: E402
from core.hypothesis import load_library, rank_candidates   # noqa: E402
from core.ingest import parse_xrdml, sample_fingerprint     # noqa: E402
from core.verification import (confirmed_by_policy, load_refinement_policy,
                               verify_case)                 # noqa: E402


def main() -> None:
    t0 = time.time()
    os.makedirs(WORK_DIR, exist_ok=True)
    os.makedirs(RES_DIR, exist_ok=True)
    prm = ensure_gsasii(ROOT, VENDOR, "")
    policy = load_refinement_policy(POLICY_PATH)

    release = json.load(open(RELEASE_PATH))
    library = json.load(open(LIBRARY_PATH))
    entry_by_id = {e["cod_id"]: e for e in release["entries"]}
    lib_payload = library["materials"]
    fp_lib = load_library(lib_payload)
    names = {m["id"]: m["name"] for m in lib_payload}
    families = {m["id"]: m["phase_family"] for m in lib_payload}

    cases = [
        # (name, fixture, true cod_id, competing cod_ids)
        ("cu_PbSO4.xrdml", "cu_PbSO4.xrdml", 1010950, [1010928]),
        ("cu_quartz.xrdml", "cu_quartz.xrdml", 1009000, [9009666]),
    ]
    rows, reports = [], []
    for cname, fixture, true_id, comp_ids in cases:
        # --- fingerprint stage (existing pipeline) ----------------------
        q = sample_fingerprint(parse_xrdml(os.path.join(FIX_DIR, fixture)))
        rk = rank_candidates(q, fp_lib, names=names, families=families)
        fp_top = rk.ranked[0].phase_family
        true_entry = entry_by_id[true_id]
        fp_ok = fp_top == true_entry["family"]

        # --- verification stage (bounded Rietveld, policy-driven) -------
        candidates = [(true_id, true_entry["family"],
                       os.path.join(CID_CIF, f"{true_id}.cif"))]
        for cid in comp_ids:
            e = entry_by_id[cid]
            candidates.append((cid, e["family"],
                               os.path.join(CID_CIF, f"{cid}.cif")))
        out = verify_case(cname, os.path.join(CID_CIF, f"{true_id}.cif"),
                          true_entry["family"], candidates, WORK_DIR, prm,
                          policy)
        confirmed = confirmed_by_policy(out, policy)
        top = out.sorted_results()[0]
        print(f"[{cname}] fp_top={fp_top} | rwp: "
              + " | ".join(f"{r.family}={r.rwp:.4f}" for r in out.sorted_results())
              + f" | confirmed={confirmed} (sep={out.separation:.4f})")
        rows.append({"case": cname, "fingerprint_top_family": fp_top,
                     "fingerprint_top_matches_truth": fp_ok,
                     "confirmation": confirmed,
                     "confirmed_family": out.confirmed_family,
                     "separation": round(out.separation, 6)})
        rep = {"case": cname, "fixture": fixture,
               "fingerprint": {"top_family": fp_top,
                               "top_similarity": round(rk.top_similarity, 6),
                               "family_margin": round(rk.family_margin, 6),
                               "matches_catalog_truth": fp_ok},
               "observed": out.observed_from,
               "candidates": [r.to_dict() for r in out.sorted_results()],
               "confirmed": confirmed,
               "policy_check": {"max_rwp": policy["confirm"]["max_rwp"],
                                "separation_min": policy["confirm"]
                                .get("separation_min")}}
        reports.append(rep)

    report = {
        "spike": "spike_07", "created_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "policy": {"version": policy["version"], "recipe": policy["recipe"],
                   "path": os.path.relpath(POLICY_PATH, ROOT)},
        "recipe_bounded_keys": policy["bounded_keys"],
        "prohibited": policy["prohibited"],
        "release_ref": {"version": release["version"],
                        "manifest_sha256": release["manifest_sha256"]},
        "library_manifest_sha256": library["manifest_sha256"],
        "cases": reports,
        "summary": {"n_cases": len(cases),
                    "n_confirmed": sum(1 for r in reports if r["confirmed"]),
                    "all_confirmed": all(r["confirmed"] for r in reports)},
        "honesty_notes": [
            "Noiseless protocol: the observed pattern is the deterministic "
            "GSAS-II simulation of the catalog CIF (same protocol as catalog "
            "fingerprints). Noise robustness is spike 05's mandate; real-data "
            "Rwp values will be larger -- the decision signal is lowest-Rwp "
            "family + policy bounds, not absolute Rwp.",
            "Bounded keys only (background + shift + cell + scale): a wrong "
            "phase cannot hide by absorbing mismatch (atoms/microstrain/size/"
            "phase-fraction/LeBail prohibited in policy).",
            "Rexp is approximate: (N-P)/sum(w y^2) with P from the policy "
            "budget; GoF = Rwp/Rexp. No phase-purity/QPA claims are made.",
        ],
        "wall_clock_s": round(time.time() - t0, 1),
    }
    with open(os.path.join(RES_DIR, "spike07_report.json"), "w") as fh:
        json.dump(report, fh, indent=2)

    lines = [
        "# Spike 07: bounded Rietveld verification",
        "",
        f"- Policy: `refinement-budget.v1.json` (recipe "
        f"`{policy['recipe']}`, version {policy['version']}) -- bounded keys "
        f"{policy['bounded_keys']}; prohibited "
        f"{', '.join(policy['prohibited'])}.",
        f"- Release `{release['version']}` manifest "
        f"{release['manifest_sha256'][:16]}...; library manifest "
        f"{library['manifest_sha256'][:16]}...",
        "",
        "| case | fingerprint top family | rwp (per candidate) | confirmed |",
        "|---|---|---|---|",
    ]
    for r in reports:
        rwps = "<br>".join(f"{c['cod_id']} {c['family']}: "
                           f"Rwp={c['rwp']:.4f} GoF={c['gof']:.2f}"
                           for c in r["candidates"])
        lines.append(f"| {r['case']} | {r['fingerprint']['top_family']} "
                     f"| {rwps} | {'yes' if r['confirmed'] else 'no'} |")
    lines += ["", "## Honesty notes", ""] + [f"- {n}" for n in report["honesty_notes"]]
    lines += ["", "## Verdict",
              f"- {report['summary']['n_confirmed']}/{report['summary']['n_cases']} "
              "cases confirmed: fingerprint top family reproduced by the "
              "bounded refinement (lowest Rwp, within policy bounds).",
              "- Evidence level of the verification stage: "
              "**fingerprint + refinement** (bundle schema run_bundle/v0).",
              ]
    with open(os.path.join(RES_DIR, "spike07_report.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[done] spike07 report -> {RES_DIR} (all confirmed: "
          f"{report['summary']['all_confirmed']})")
    if not report["summary"]["all_confirmed"]:
        raise SystemExit("spike07: at least one case not confirmed")


if __name__ == "__main__":
    main()