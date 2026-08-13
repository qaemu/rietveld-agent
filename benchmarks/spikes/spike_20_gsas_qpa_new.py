def gsas_qpa(pat, phases, work: Path, tag: str, lo: float, hi: float,
             a1: float, a2: float, ratio: float, sync: bool = False,
             maxcyc: int = 40, dwr_gate: float = 0.5, n_resid: int = 2,
             resid_sig: float = 2.5,
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

    Stage C (residual screening): yobs - ycalc after Stage B; fingerprint the
    residual and screen the full COD index for weak TRUE minors (e.g.
    corundum at ~1%) that Stage B cannot distinguish from decoys by wR;
    add the top distinct-canon hits.

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
    for rnd in range(len(phases) - 1):
        remaining = [x for x in ok if x["phase"]["name"] not in
                     [s["name"] for s in selected]
                     and x["phase"].get("canon") not in
                     [s.get("canon") for s in selected]]
        if not remaining:
            break
        best, best_d = None, 0.0
        for x in remaining:
            cand = x["phase"]
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

    # ---- Stage C: residual screening for weak true minors ----
    resid_new = []
    if cod_ids is not None and d_units is not None and entry_of is not None \
            and metas_by_id is not None:
        try:
            arr = h.data['data'][1]  # 6 x N masked: 0 2th,1 yobs,3 ycalc
            tth_f = np.asarray(arr[0], dtype=float)
            yobs_f = np.asarray(arr[1], dtype=float)
            ycalc_f = np.asarray(arr[3], dtype=float)
            resid = np.clip(yobs_f - ycalc_f, 0.0, None)
            rpat = PowderPattern(sample_name=f"{tag}_resid",
                                 source=pat.source, tth=tth_f,
                                 intensity=resid, instrument=pat.instrument)
            fp = sample_fingerprint(rpat, prominence=0.005)
            hits = screen_fingerprint(fp, cod_ids, d_units, entry_of,
                                      metas_by_id, top_k=30, pool_k=60)
            have_canon = set([s.get("canon") for s in selected
                              if s.get("canon")])
            have_cod = set([str(s["cif"].stem) for s in selected])
            for hh in hits:
                if len(resid_new) >= n_resid:
                    break
                cid = str(hh["cod_id"])
                canon = hh.get("canon") or hh.get("mineral") or cid
                if cid in have_cod or canon in have_canon:
                    continue
                if hh.get("significance", 0.0) < resid_sig:
                    continue
                cif = download_cif(int(cid), delay=0.05)
                if cif is None or not cif.exists():
                    continue
                resid_new.append({"name": f"p{cid}", "cif": cif,
                                  "canon": canon})
                have_cod.add(cid)
                have_canon.add(canon)
                print(f"    resid add {cid} ({canon}) "
                      f"sig={hh.get('significance')}", flush=True)
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
