"""Unit 04 (b): M1 vertical slice -- analyze everything, record bundles.

Cases:
  cu_PbSO4  -> supported PbSO4          (same instrument as library)
  cu_quartz -> supported SiO2           (same instrument)
  fe_PbSO4  -> supported PbSO4          (cross-anode; d-space library matching)
  unknown   -> abstain                  (MoKa: no released calibration - hard stop)

Each case produces a schema-validated RunBundle under data/unit04/results/
"""
from __future__ import annotations

import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES_DIR = os.path.join(ROOT, "data", "unit04", "results")
BUNDLES = os.path.join(RES_DIR, "bundles")
FIX_DIR = os.path.join(ROOT, "tests", "fixtures", "xrdml")

sys.path.insert(0, ROOT)

from cli.analyze import analyze  # noqa: E402
from core.calibration import CalibrationRegistry, ResolutionStatus  # noqa: E402
from core.ingest import InstrumentParams  # noqa: E402
from core.report import load_bundle_schema  # noqa: E402

REGISTRY = os.path.join(ROOT, "data", "unit03", "results", "registry.json")
LIBRARY = os.path.join(ROOT, "data", "candidates", "library.json")


def main() -> None:
    os.makedirs(BUNDLES, exist_ok=True)
    report: dict = {"started_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    schema = load_bundle_schema()

    cases = [
        ("cu_PbSO4.xrdml", "supported", "PbSO4 (anglesite)"),
        ("cu_quartz.xrdml", "supported", "GaAsO4 (quartz homeotype)"),
        ("fe_PbSO4.xrdml", "supported", "PbSO4 (anglesite)"),
    ]
    rows = []
    for fname, exp_status, exp_family in cases:
        path = os.path.join(FIX_DIR, fname)
        out = os.path.join(BUNDLES, f"{fname}.bundle.json")
        bundle = analyze(path, REGISTRY, LIBRARY, out)
        with open(out, "w") as fh:
            json.dump(bundle, fh, indent=2)
        json_schema = __import__("jsonschema")
        json_schema.validate(bundle, schema)
        verdict = bundle["verdicts"][0] if bundle["verdicts"] else None
        got_status = verdict["status"] if verdict else "abstain"
        got_family = verdict["phase_family"] if verdict else None
        ok = got_status == exp_status and (exp_family is None or got_family == exp_family)
        sim = verdict["evidence"]["top_similarity"] if verdict else 0.0
        print(f"[case:{fname}] status={got_status} primary={got_family} "
              f"sim={sim:.3f}  ok={ok}")
        rows.append({"fixture": fname, "expected_status": exp_status,
                     "expected_family": exp_family, "got_status": got_status,
                     "got_family": got_family, "sim": round(sim, 4), "ok": ok,
                     "bundle": out})
        if not ok:
            raise SystemExit(f"case failed: {fname}")

    # --- unknown instrument: hard stop before ranking -----------------------
    reg = CalibrationRegistry.load(REGISTRY)
    mo = InstrumentParams(anode="MoKa", wavelengths=(0.7093, 0.7136),
                          scan_axis="2Theta/Theta")
    res = reg.lookup(mo)
    ok = res.status == ResolutionStatus.UNKNOWN
    print(f"[case:unknown] calibration={res.status.value}  ok={ok}")
    rows.append({"fixture": "unknown_MoKa(no file, params only)",
                 "expected_status": "unknown_calibration", "got_status": res.status.value,
                 "ok": ok, "reason": res.reason})
    if not ok:
        raise SystemExit("unknown-anode case failed")

    report["cases"] = rows
    report["library_manifest"] = json.load(open(LIBRARY))["manifest_sha256"][:16]
    report["n_all_passed"] = all(r["ok"] for r in rows)

    with open(os.path.join(RES_DIR, "unit04_report.json"), "w") as fh:
        json.dump(report, fh, indent=2, default=str)

    md = [
        "# Unit 04: M1 vertical slice -- analyze",
        "",
        f"- library: {os.path.relpath(LIBRARY, ROOT)} "
        f"(manifest {report['library_manifest']}...)",
        f"- registry: {os.path.relpath(REGISTRY, ROOT)}",
        "",
        "## Cases",
        "",
        "| fixture | expected | got | top sim | bundle |",
        "|---|---|---|---|---|",
    ]
    for r in rows:
        bundle = os.path.relpath(r["bundle"], ROOT) if r.get("bundle") else "-"
        md.append(f"| {r['fixture']} | {r['expected_status']} | {r['got_status']} "
                  f"| {r.get('sim', '-')} | {bundle} |")
    md += [
        "",
        "## Findings",
        "- Fe measurement resolves against a Cu-only fingerprint library in "
        "d-space (mat-pbso4 cross-anode), so calibration (not library) carries "
        "the instrument coupling.",
        "- Abstention is total-precedence: calibration unknown/ambiguous stops "
        "before any hypothesis work (policy: released-only).",
        "- M1 evidence is fingerprint-only; thresholds (top sim >= 0.35, "
        "margin >= 0.10) are explicit constants in core/verdict pending policy "
        "promotion when refinement-backed evidence lands.",
        "- SUPERSEDED by unit 05: the 0.108 cross-anode margin shown in M1 was "
        "an exact-match artifact (query == its own library entry). The margin "
        "is now family-aware (against the best DIFFERENT phase family); "
        "see data/unit05/results/unit05_report.md. M1 verdicts unchanged, "
        "same-anode and cross-anode both hold.",
        "- Catalogue-backing (unit 06): the M1 sampling of quartz-family "
        "patterns was chemistry-validated as GaAsO4 (quartz homeotype) and "
        "the 'NaCl' spot as Ag0.5Bi0.5S; the library and these expectations "
        "now reflect the validated identities (data/catalog/releases/"
        "catalog_0.1.0.json).",
        "- Every bundle validates against run_bundle.schema.json at write time "
        "and in the test suite.",
        "",
        "## Verdict",
        "- [ ] cu_PbSO4 -> supported PbSO4",
        "- [ ] cu_quartz -> supported SiO2",
        "- [ ] fe_PbSO4 -> supported PbSO4 (cross-anode)",
        "- [ ] unknown MoKa -> hard-stop abstain (calibration unknown)",
        "- [ ] all bundles schema-valid",
    ]
    with open(os.path.join(RES_DIR, "unit04_report.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    print("wrote", os.path.join(RES_DIR, "unit04_report.json"))
    if not report["n_all_passed"]:
        raise SystemExit("unit04: at least one case failed")


if __name__ == "__main__":
    main()