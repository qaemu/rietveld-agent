"""CLI: rietveld-agent analyze -- the full M1 vertical slice.

    python -m cli analyze data.foo.xrdml [--output bundle.json]
                                       [--registry registry.json]
                                       [--library library.json]

Pipeline: parse XRDML -> instrument fingerprint -> calibration resolution
(released-only; hard stop on unknown/ambiguous) -> hypothesis ranking
(d-space similarity vs the candidate library) -> governed verdict -> [bounded
Rietveld verification of the top family + nearest competing family against
the MEASURED pattern (spike 08; skipped on abstain/held or when catalog
CIFs are unavailable; any measured grid is resampled onto the protocol
grid, and counting statistics are assessed, never gating the
identification)] -> schema-validated RunBundle with optional
``verification`` evidence.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

from core._paths import repo_root
from core.calibration import CalibrationRegistry, ResolutionStatus
from core.hypothesis import load_library, rank_candidates
from core.ingest import parse_xrdml, sample_fingerprint
from core.report import build_run_bundle
from core.verdict import decide

_here_root = repo_root()
ROOT = _here_root if (_here_root / "data").is_dir() else None
if ROOT is not None:
    DEFAULT_REGISTRY = ROOT / "data" / "spike3" / "results" / "registry.json"
    DEFAULT_LIBRARY = ROOT / "data" / "candidates" / "library.json"
    DEFAULT_POLICY = ROOT / "governance" / "policies" / "refinement-budget.v1.json"
    VERIFY_WORK = ROOT / "data" / "spike8" / "work"
    CID_CIF_DIR = ROOT / "data" / "spike6" / "input" / "cod"
    VENDOR = ROOT / ".vendor" / "GSAS-II"
else:
    DEFAULT_REGISTRY = DEFAULT_LIBRARY = DEFAULT_POLICY = None
    VERIFY_WORK = CID_CIF_DIR = VENDOR = None
#: simulated-protocol grid produced by add_simulated_powder_histogram for the
#: fixture range; measured patterns on ANY grid are accepted for verification
#: (resampled onto this protocol grid, recorded in the evidence; see
#: _verify_supported).
_GRID_POINTS = 6251


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def _cod_id(material_id: str, lib_entry: dict) -> int:
    """Material -> catalog cod_id (cod-* ids directly; legacy via catalog_ref)."""
    if material_id.startswith("cod-"):
        return int(material_id[len("cod-"):])
    ref = (lib_entry or {}).get("catalog_ref") or {}
    return ref.get("cod_id", 0)


def _verification_candidates(ranking, lib_payload: dict) -> list:
    """(cod_id, family, cif_path) for the top family + nearest competing
    family; only when both resolve to catalog CIFs on disk."""
    entries = {m["id"]: m for m in lib_payload["materials"]}
    out, seen_families = [], set()
    for c in ranking.ranked:
        if len(out) >= 2:
            break
        if c.phase_family in seen_families and out:
            continue                      # same-family twins don't compete
        lib = entries.get(c.material_id)
        cod = _cod_id(c.material_id, lib)
        cif = CID_CIF_DIR / f"{cod}.cif"
        if cod <= 0 or not cif.exists():
            continue
        seen_families.add(c.phase_family)
        out.append((cod, c.phase_family, str(cif)))
    return out


def _verify_supported(pattern, ranking, lib_payload: dict,
                      policy_path: str, work_dir: str) -> dict:
    """Bounded verification of top + nearest competing family against the
    measured pattern; returns bundle `verification` evidence (or {} if the
    grid does not resolve against the protocol or no catalog CIFs resolve).

    Measured patterns on ANY grid are accepted: the observed profile is
    resampled (linear interpolation in 2theta, zero outside the measured
    range) onto the simulated-protocol grid the candidate models are built
    on, and the resampling is recorded in the evidence.
    """
    candidates = _verification_candidates(ranking, lib_payload)
    if len(candidates) < 2:
        return {}
    from benchmarks.eval.sim import ensure_gsasii
    from core.verification import load_refinement_policy, verify_measured
    prm = ensure_gsasii(str(ROOT), str(VENDOR), "")
    policy = load_refinement_policy(str(policy_path))
    import numpy as np

    measured_tth = np.asarray(pattern.tth)
    measured_y = np.asarray(pattern.intensity)
    grid_note = None
    if measured_tth.size != _GRID_POINTS:
        # resample the measurement onto the protocol grid the candidate
        # models are simulated on (honest: recorded, zero-extended outside
        # the measured range; zero-count bins are dropped by the engine)
        proto = policy["protocol"]
        step = float(proto["step"])
        tmin, tmax = float(proto["tmin"]), float(proto["tmax"])
        proto_tth = np.linspace(tmin, tmax, _GRID_POINTS)
        measured_y = np.interp(proto_tth, measured_tth, measured_y,
                               left=0.0, right=0.0)
        grid_note = {
            "resampled": True,
            "method": "linear interpolation in 2theta onto the "
                      "simulated-protocol grid (zero outside measured range)",
            "measured_grid": {"n_points": int(measured_tth.size),
                              "tmin": float(measured_tth[0]),
                              "tmax": float(measured_tth[-1])},
            "protocol_grid": {"n_points": int(_GRID_POINTS),
                              "tmin": tmin, "tmax": tmax, "step": step},
        }
        measured_tth = proto_tth

    os.makedirs(work_dir, exist_ok=True)
    outcome = verify_measured(
        measured_tth, measured_y,
        case=Path(pattern.sample_name or "measurement").stem,
        candidates=candidates, work_dir=work_dir, prm_path=prm, policy=policy)
    evidence = outcome.to_bundle_evidence()
    evidence["consistent_with_fingerprint"] = (
        evidence["confirmed_family"] == (ranking.ranked[0].phase_family
                                         if ranking.ranked else None))
    if grid_note:
        evidence["grid"] = grid_note
    return evidence


def analyze(xrdml_path: str, registry_path: str = str(DEFAULT_REGISTRY),
            library_path: str = str(DEFAULT_LIBRARY),
            output_path: str = "",
            verification: bool = True,
            full_cod: bool = False,
            policy_path: str = str(DEFAULT_POLICY),
            work_dir: str = str(VERIFY_WORK)) -> dict:
    """Run the M1 analyze pipeline; returns the schema-valid bundle dict."""
    raw_sha = _sha256(xrdml_path)
    pattern = parse_xrdml(xrdml_path)

    registry = CalibrationRegistry.load(registry_path)
    lib_payload = json.loads(Path(library_path).read_text())

    # --- hard stop first: instrument must be uniquely resolved
    resolution = registry.lookup(pattern)
    empty_ranking = rank_candidates(sample_fingerprint(pattern), {})
    verdict = decide(empty_ranking, resolution)
    if resolution.status != ResolutionStatus.RESOLVED:
        bundle = build_run_bundle(
            input_path=xrdml_path, raw_sha256=raw_sha,
            sample_name=pattern.sample_name,
            resolution=resolution, ranking=empty_ranking, verdict=verdict,
            library_manifest=lib_payload.get("manifest_sha256", "")[:12],
            artifact_path=output_path)
        return bundle

    # --- hypothesis + verdict
    library = load_library(lib_payload["materials"])
    fingerprint = sample_fingerprint(pattern)
    ranking = rank_candidates(fingerprint, library,
                              names={m["id"]: m["name"] for m in
                                     lib_payload["materials"]},
                              families={m["id"]: m["phase_family"] for m in
                                        lib_payload["materials"]})
    verdict = decide(ranking, resolution)

    cod_screen = None
    if full_cod:
        try:
            from core.codsearch import screen_cod
            cod_screen = screen_cod(fingerprint)
        except FileNotFoundError as exc:
            print(f"warning: full-COD screen skipped ({exc})",
                  file=sys.stderr)
        except Exception as exc:               # noqa: BLE001 -- screen must
            print(f"warning: full-COD screen failed ({exc})", file=sys.stderr)

    evidence = {}
    if verification and verdict.status == "supported" and ranking.ranked:
        try:
            evidence = _verify_supported(pattern, ranking, lib_payload,
                                         policy_path, work_dir)
        except Exception as exc:          # noqa: BLE001 -- verification must
            print(f"warning: verification skipped ({exc})", file=sys.stderr)

    bundle = build_run_bundle(
        input_path=xrdml_path, raw_sha256=raw_sha,
        sample_name=pattern.sample_name,
        resolution=resolution, ranking=ranking, verdict=verdict,
        library_manifest=lib_payload.get("manifest_sha256", "")[:12],
        artifact_path=output_path,
        verification=evidence or None,
        cod_screen=cod_screen)
    return bundle


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="rietveld-agent",
                                 description="policy-governed XRD phase identification")
    if ROOT is None:
        ap.error(
            "repository data not found: run rietveld-agent from inside a "
            "clone of github.com/qaemu/rietveld-agent"
        )
    sub = ap.add_subparsers(dest="command", required=True)

    p = sub.add_parser("analyze", help="identify phases in an XRDML measurement")
    p.add_argument("xrdml", help="input XRDML file")
    p.add_argument("--output", default="", help="bundle JSON output path "
                                                "(default: <input>.bundle.json)")
    p.add_argument("--registry", default=str(DEFAULT_REGISTRY))
    p.add_argument("--library", default=str(DEFAULT_LIBRARY))
    p.add_argument("--no-verification", action="store_true",
                   help="skip the bounded Rietveld verification stage")
    p.add_argument("--full-cod", action="store_true",
                   help="also screen the fingerprint against the COMPLETE "
                        "COD line index (core.codsearch; requires "
                        "`make cod-index`)")
    p.add_argument("-v", "--verbose", action="store_true")

    args = ap.parse_args(argv)
    if args.command != "analyze":
        ap.error(f"unknown command: {args.command}")

    out = args.output or f"{args.xrdml}.bundle.json"
    bundle = analyze(args.xrdml, args.registry, args.library, out,
                     verification=not args.no_verification,
                     full_cod=args.full_cod)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps(bundle, indent=2))

    print(f"run_id      : {bundle['run_id']}")
    print(f"status      : {bundle['status']}")
    print(f"verdict     : {bundle['verdicts'][0]['status'] if bundle['verdicts'] else 'abstain'}")
    if bundle.get("cod_screen"):
        top = bundle["cod_screen"]["top"][0]
        print(f"full-COD    : {bundle['cod_screen']['index_entries']} entries "
              f"screened; top COD {top['cod_id']} "
              f"({(top.get('formula') or '')[:28]}) sig={top['significance']}")
    if bundle.get("verification"):
        v = bundle["verification"]
        print(f"verified    : {v['confirmed_family']} "
              f"(evidence={v['evidence_level']}, "
              f"consistent={v.get('consistent_with_fingerprint')})")
    if bundle["verdicts"]:
        v = bundle["verdicts"][0]
        print(f"primary     : {v['phase_family']}  "
              f"(sim={v['evidence']['top_similarity']:.3f}, "
              f"margin={v['evidence']['margin']:.3f})")
    else:
        ck = [c for c in bundle["checkpoints"] if c["stage"] == "verdict"][0]
        print(f"reasons     : {ck['metrics']['reasons']}")
    print(f"bundle      : {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())