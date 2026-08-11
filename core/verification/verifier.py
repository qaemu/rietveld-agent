"""Bounded Rietveld verification (spike 07).

Deterministic, policy-driven verification stage: the *observed* pattern is a
GSAS-II simulation of the catalog CIF (same protocol that builds the catalog
fingerprints), and each candidate phase is refined against that observation
using ONLY the bounded parameter budget from ``refinement-budget.v1.json``
(background + sample shift + cell + scale). Prohibited keys (atoms,
microstrain, size, phase fractions, LeBail, ...) can never be used, so a
wrong phase cannot hide by absorbing mismatch.

Decision: the phase family with the lowest bounded Rwp is the refinement-
supported family; ``confirm`` additionally requires Rwp <= policy max_rwp
and a minimum separation to the next-best candidate.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Optional

import jsonschema

from benchmarks.eval.sim import ensure_gsasii, sim_cif_to_pattern

#: refined background coeff count (protocol.background) + scale + shift
_FIXED_PARAMS = 5
#: lattice cell parameters refined when phase "Cell" is on, per crystal system
_NCELL = {"Cubic": 1, "Tetragonal": 2, "Hexagonal": 2, "Trigonal": 2,
          "Orthorhombic": 3, "Monoclinic": 4, "Triclinic": 6}
#: counting-statistics floor for MEASURED-path verification claims
#: (spike 05, envelope L3 = laboratory-grade counting statistics). The
#: fingerprint identification itself is position-based and counts-
#: independent (core.verdict); refinements against weak data are
#: documented as "statistics-below-gate" rather than claiming a
#: verified fit and never block the identification verdict.
VERIFY_MIN_PEAK_COUNTS = 100_000.0


@dataclass
class RefinementResult:
    family: str
    cod_id: int
    rwp: float
    rexp: float
    gof: float
    n_points: int
    n_params: int
    converged: bool
    cell: Optional[list] = None

    def to_dict(self) -> dict:
        return {"family": self.family, "cod_id": self.cod_id,
                "rwp": round(float(self.rwp), 6),
                "rexp": round(float(self.rexp), 6),
                "gof": round(float(self.gof), 6),
                "n_points": int(self.n_points), "n_params": int(self.n_params),
                "converged": bool(self.converged),
                "cell": [round(float(v), 4) for v in self.cell]
                if self.cell else None}


@dataclass
class VerificationOutcome:
    case: str
    observed_from: dict
    truthtable_family: str
    results: list = field(default_factory=list)
    confirmed_family: Optional[str] = None
    separation: Optional[float] = None
    policy: dict = field(default_factory=dict)
    #: measured-path statistics assessment ("verified" when the counting
    #: statistics support a refinement claim, else "statistics-below-gate")
    status: str = "verified"
    statistics: Optional[dict] = None
    reasons: list = field(default_factory=list)

    def sorted_results(self) -> list:
        return sorted(self.results, key=lambda r: r.rwp)

    def to_bundle_evidence(self) -> dict:
        ev = {
            "recipe": self.policy.get("recipe", "refinement-verify-v1"),
            "policy_version": self.policy.get("version", "?"),
            "evidence_level": "fingerprint + refinement",
            "status": self.status,
            "confirmed_family": self.confirmed_family,
            "separation": round(float(self.separation), 6)
            if self.separation is not None else None,
            "candidates": [r.to_dict()
                           for r in self.sorted_results()],
        }
        if self.statistics is not None:
            ev["statistics"] = dict(self.statistics)
        if self.reasons:
            ev["reasons"] = list(self.reasons)
        if self.observed_from is not None:
            ev["observed_from"] = self.observed_from
        return ev


def load_refinement_policy(path: str) -> dict:
    """Load and schema-validate the bounded-verification policy record."""
    with open(path) as fh:
        policy = json.load(fh)
    schema_path = os.path.join(os.path.dirname(os.path.abspath(path)),
                               "..", "schemas",
                               "refinement_policy.schema.json")
    with open(schema_path) as fh:
        schema = json.load(fh)
    jsonschema.validate(policy, schema)
    return policy


def _n_cell_params(proj) -> int:
    try:
        sys_name = proj.phase(0).data["General"]["SGData"]["SGSys"]
        for key, n in _NCELL.items():
            if key in sys_name:
                return n
    except Exception:                                   # noqa: BLE001
        pass
    return 3                                            # conservative fallback


def _refined_param_count(policy: dict, proj) -> int:
    coeffs = int(policy["protocol"]["background"]["no. coeffs"])
    return coeffs + 1 + 1 + _n_cell_params(proj)  # bkg+scale+shift+cell


def make_observed(true_cif: str, work_dir: str, prm_path: str,
                  policy: dict) -> tuple:
    """Observed pattern = deterministic simulation of the catalog CIF.

    Cached to keep repeated runs fast; the cache content is deterministic
    (same CIF + protocol -> same arrays).
    """
    import numpy as np
    os.makedirs(work_dir, exist_ok=True)
    cache = os.path.join(work_dir,
                         f"obs_{os.path.splitext(os.path.basename(true_cif))[0]}.npz")
    if os.path.exists(cache):
        with np.load(cache) as z:
            return z["tth"], z["yobs"], z["wy"]
    proto = policy["protocol"]
    pat = sim_cif_to_pattern(
        true_cif, work_dir, prm_path=prm_path,
        tmin=float(proto["tmin"]), tmax=float(proto["tmax"]),
        step=float(proto["step"]), scale=float(proto["scale"]))
    tth = np.asarray(pat.tth)
    yobs = np.asarray(pat.intensity)
    wy = 1.0 / np.maximum(yobs, 1e-9)
    np.savez(cache, tth=tth, yobs=yobs, wy=wy)
    return tth, yobs, wy


def _prepare_refinements(hist, proj, policy: dict) -> None:
    """Apply the bounded parameter budget to histogram + phase 0."""
    proto = policy["protocol"]
    hist.set_refinements({"Background": dict(proto["background"],
                                             refine=True)})
    hist.set_refinements({"Sample Parameters": ["Shift"]})
    proj.phase(0).set_refinements({"Cell": True})
    proj.phase(0).set_HAP_refinements({"Scale": True})
    proj.data["Controls"]["data"]["max cyc"] = int(policy["max_cycles"])


def refine_candidate(cod_id: int, family: str, cif_path: str,
                     tth, yobs, wy=None, work_dir: str = "",
                     prm_path: str = "", policy: dict = None) -> RefinementResult:
    """Bounded Rietveld refinement of one candidate against the observed
    pattern (spike-07 sim-observed path, deterministic GSAS-II).

    The observed pattern is a GSAS-II simulation of the catalog CIF; the
    candidate is refined against it with ONLY the bounded budget. Repeated
    runs are bit-identical: fresh project per call (stale artifacts purged)
    and the Dummy/simulate histogram re-materializes the same deterministic
    pattern. The measured pattern path uses refine_measured_candidate.
    """
    import numpy as np
    from GSASII.GSASIIscriptable import G2Project

    if wy is None:
        wy = 1.0 / np.maximum(yobs, 1e-9)
    tag = f"{cod_id}_{family.replace(' ', '_')}"
    base = os.path.join(work_dir, tag)
    # always refine from a FRESH project: G2Project(newgpx=...) re-opens an
    # existing file and accumulates phases/histograms, corrupting repeats
    for suffix in (".gpx", ".bak0.gpx", "_final.gpx", ".lst"):
        stale = base + suffix
        if os.path.exists(stale):
            os.remove(stale)
    gpx = base + ".gpx"
    proj = G2Project(newgpx=gpx)
    proj.add_phase(cif_path, phasename=tag, fmthint="CIF")

    proto = policy["protocol"]
    hist = proj.add_simulated_powder_histogram(
        "obs", prm_path, float(np.min(tth)), float(np.max(tth)),
        Tstep=float(np.median(np.diff(tth))),
        scale=float(proto["scale"]), phases=proj.phases())
    # overwrite the simulated observation with the true observed arrays
    d1 = hist.data["data"][1]                 # [x, yobs, wy, ycalc, ...]
    d1[1] = yobs.copy(); d1[2] = wy.copy(); d1[3] = np.zeros_like(yobs)
    d0 = hist.data["data"][0]                 # points dict (I/W aliases)
    d0["I"] = yobs.copy(); d0["W"] = wy.copy()

    _prepare_refinements(hist, proj, policy)
    proj.refine()

    yc = np.asarray(hist.getdata("ycalc"))
    n = int(yc.size)
    num = float(np.sum(wy * (yobs - yc) ** 2))
    den = float(np.sum(wy * yobs ** 2))
    rwp = float(np.sqrt(num / den)) if den > 0 else 1.0
    n_params = _refined_param_count(policy, proj)
    rexp = float(np.sqrt(max(n - n_params, 1) / den)) if den > 0 else 1.0
    gof = float(rwp / rexp) if rexp > 0 else 0.0

    cell = None
    try:
        cell = list(proj.phase(0).data["General"]["Cell"][1][:6])
    except Exception:                                   # noqa: BLE001
        pass
    lst = os.path.splitext(gpx)[0] + ".lst"
    converged = bool(os.path.exists(lst)
                     and "Refinement successful" in open(lst, errors="ignore").read())
    proj.save(os.path.splitext(gpx)[0] + "_final.gpx")
    return RefinementResult(family=family, cod_id=cod_id, rwp=rwp,
                            rexp=rexp, gof=gof, n_points=n, n_params=n_params,
                            converged=converged, cell=cell)


def verify_case(case: str, true_cif: str, true_family: str,
                candidates: list, work_dir: str, prm_path: str,
                policy: dict) -> VerificationOutcome:
    """Refine every candidate (cod_id, family, cif_path) against the
    observed pattern of the catalog CIF; return the sorted outcome."""
    tth, yobs, wy = make_observed(true_cif, work_dir, prm_path, policy)
    outcome = VerificationOutcome(
        case=case,
        observed_from={"cod_id": int(os.path.splitext(
            os.path.basename(true_cif))[0]),
                       "family": true_family,
                       "cif": os.path.basename(true_cif)},
        truthtable_family=true_family,
        policy={"version": policy["version"],
                "recipe": policy["recipe"]})
    for cod_id, family, cif_path in candidates:
        outcome.results.append(refine_candidate(
            cod_id, family, cif_path, tth, yobs, wy, work_dir, prm_path,
            policy))
    ranked = outcome.sorted_results()
    outcome.confirmed_family = ranked[0].family if ranked else None
    outcome.separation = (ranked[1].rwp - ranked[0].rwp) if len(ranked) > 1 \
        else None
    return outcome


def _sim_candidate_pattern(cod_id: int, family: str, cif_path: str,
                           work_dir: str, prm_path: str,
                           policy: dict) -> tuple:
    """Deterministic candidate model: the protocol simulation of the
    candidate CIF (same recipe/build as the catalog fingerprints). Returns
    (tth, ycand) with GSAS-II provenance in the pattern metadata. Cached in
    ``work_dir`` (deterministic: fixed CIF + protocol -> fixed pattern)."""
    import numpy as np
    proto = policy["protocol"]
    os.makedirs(work_dir, exist_ok=True)
    cache = os.path.join(work_dir, f"sim_{int(cod_id)}.npz")
    if os.path.exists(cache):
        with np.load(cache) as z:
            return z["tth"], z["yc"]
    pat = sim_cif_to_pattern(
        cif_path, work_dir, prm_path=prm_path,
        anode=proto["anode"], wavelengths=tuple(proto["wavelengths"]),
        tmin=float(proto["tmin"]), tmax=float(proto["tmax"]),
        step=float(proto["step"]), scale=float(proto["scale"]))
    tth = np.asarray(pat.tth)
    yc = np.asarray(pat.intensity)
    np.savez(cache, tth=tth, yc=yc)
    return tth, yc


def refine_measured_candidate(cod_id: int, family: str, cif_path: str,
                              tth, yobs, work_dir: str, prm_path: str,
                              policy: dict,
                              wy: Optional[list] = None) -> RefinementResult:
    """Bounded verification fit of one candidate against the MEASURED
    pattern (numpy engine, deterministic).

    The observed grid + counting weights come from the measurement; the
    model is the protocol simulation of the candidate CIF (GSAS-II
    forward model, same recipe as the catalog fingerprints). Only the
    bounded key budget of the policy is fitted: histogram background
    (chebyschev-1, ``no. coeffs``) + Sample Shift + phase Scale. Cell is
    fixed (the fingerprint stage already constrained the phase family);
    the budget mirrors the allowlist, so a wrong phase cannot absorb
    structural mismatch. scipy trust-region LM is deterministic for the
    fixed start point and bounds.
    """
    import numpy as np
    from scipy.optimize import least_squares

    tth = np.asarray(tth, dtype=float)
    yobs = np.asarray(yobs, dtype=float)
    if wy is None:
        # counting weights with a 1-count floor: yobs == 0 bins must not
        # blow up (1e9 weight) nor vanish (sigma floor of 1 count)
        wy = 1.0 / np.maximum(yobs, 1.0)
    else:
        wy = np.asarray(wy, dtype=float)
    # EXCLUDED REGIONS (Rietveld practice): zero-count bins are not
    # informative (e.g. sample-holder shadows) and cannot be fit by any
    # bounded model -- drop them from the residual sum, Rwp and n_points.
    keep = yobs >= 1.0
    tth_m, yobs_m, wy_m = tth[keep], yobs[keep], wy[keep]
    n = int(yobs_m.size)
    ctth, ycand = _sim_candidate_pattern(cod_id, family, cif_path,
                                         work_dir, prm_path, policy)
    assert ctth.size == tth.size, "candidate sim grid must match observed"
    tm = np.arange(n, dtype=float) / max(n - 1, 1)   # 0..1 over fitted grid

    n_bkg = int(policy["protocol"]["background"]["no. coeffs"])

    def model(par):
        ln_s, shift, *c = par
        s = float(np.exp(ln_s))
        xt = tth_m - shift
        yc = np.interp(xt, ctth, ycand, left=0.0, right=0.0)
        bg = np.zeros_like(yc)
        for i in range(n_bkg):
            bg = bg + c[i] * np.cos(i * np.arccos(2.0 * tm - 1.0))
        return s * yc + bg

    s0 = max(float(np.max(yobs_m)), 1e-9) / max(float(np.max(ycand)), 1e-9)
    p0 = [np.log(s0), 0.0] + [0.0] * n_bkg
    bounds = ([-np.inf, -0.75] + [-1e9] * n_bkg,
              [np.inf, 0.75] + [1e9] * n_bkg)
    sqw = np.sqrt(wy_m)
    resid = lambda p: sqw * (yobs_m - model(p))      # noqa: E731
    res = least_squares(resid, p0, bounds=bounds,
                        max_nfev=80, x_scale="jac", verbose=0)
    yc = model(res.x)
    num = float(np.sum(wy_m * (yobs_m - yc) ** 2))
    den = float(np.sum(wy_m * yobs_m ** 2))
    rwp = float(np.sqrt(num / den)) if den > 0 else 1.0
    n_params = 2 + n_bkg                     # scale + shift + bkg coeffs
    rexp = float(np.sqrt(max(n - n_params, 1) / den)) if den > 0 else 1.0
    gof = float(rwp / rexp) if rexp > 0 else 0.0
    converged = bool(res.status in (1, 2, 3))     # gtol/ftol/xtol reached
    return RefinementResult(family=family, cod_id=cod_id, rwp=rwp,
                            rexp=rexp, gof=gof, n_points=n, n_params=n_params,
                            converged=converged, cell=None)


def verify_measured(tth, yobs, case: str, candidates: list,
                    work_dir: str, prm_path: str, policy: dict,
                    wy: Optional[list] = None,
                    min_peak_counts: float = VERIFY_MIN_PEAK_COUNTS
                    ) -> VerificationOutcome:
    """Verification against a MEASURED pattern (no truth table).

    ``candidates``: (cod_id, family, cif_path) triples. The observed grid and
    counting weights come from the measurement itself (wy = 1/yobs when not
    given). outcome.truthtable_family stays None -> only lowest-Rwp + policy
    bounds decide; the fingerprint stage carries the identification.

    Counting statistics are ASSESSED, never gate the analysis: the fit still
    runs so the report documents what the weak data can and cannot support,
    and ``status`` records whether the statistics meet the calibrated floor
    (``min_peak_counts``, spike-05 envelope L3) for a refinement claim.
    """
    import numpy as np
    tth = np.asarray(tth, dtype=float)
    yobs = np.asarray(yobs, dtype=float)
    if wy is None:
        wy = 1.0 / np.maximum(yobs, 1.0)
    peak_max = float(np.max(yobs)) if yobs.size else 0.0
    statistics = {"peak_max": peak_max,
                  "min_peak_counts": float(min_peak_counts),
                  "satisfied": bool(peak_max >= min_peak_counts)}
    outcome = VerificationOutcome(
        case=case, observed_from=None, truthtable_family=None,
        policy={"version": policy["version"], "recipe": policy["recipe"]},
        status="verified" if statistics["satisfied"] else "statistics-below-gate",
        statistics=statistics)
    if not statistics["satisfied"]:
        outcome.reasons.append(
            f"counting statistics below verification gate "
            f"(peak_max={peak_max:.5g} < {min_peak_counts:.5g}): the "
            f"refinement cannot claim verification of the phase family")
    for cod_id, family, cif_path in candidates:
        outcome.results.append(refine_measured_candidate(
            cod_id, family, cif_path, tth, yobs, work_dir, prm_path,
            policy, wy=wy))
    ranked = outcome.sorted_results()
    outcome.confirmed_family = ranked[0].family if ranked else None
    outcome.separation = (ranked[1].rwp - ranked[0].rwp) if len(ranked) > 1 \
        else None
    return outcome


def confirmed_by_policy(outcome: VerificationOutcome, policy: dict) -> bool:
    """Confirm only if: lowest-Rwp family matches the expectation, its Rwp is
    within policy max, and separation >= policy minimum."""
    conf = policy["confirm"]
    ranked = outcome.sorted_results()
    if not ranked:
        return False
    top = ranked[0]
    if top.family != outcome.truthtable_family:
        return False
    if top.rwp > float(conf["max_rwp"]):
        return False
    if conf.get("require_separation", True) and outcome.separation is not None \
            and outcome.separation < float(conf.get("separation_min", 0.0)):
        return False
    return True