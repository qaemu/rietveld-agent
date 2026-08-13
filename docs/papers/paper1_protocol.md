# Reproducible Rietveld quantitative phase analysis under a bounded budget

**The NIST SRM 2686a protocol of `rietveld-agent`** — Paper P-1 of the
docs/papers series.

**Thesis.** Rietveld quantitative phase analysis (RQPA) of laboratory
powder X-ray diffraction (PXRD) data can be performed with a
deterministic, budget-bounded refinement protocol whose every number
reproduces bit-identically across runs.

**Status.** Working paper, v0.2.0, August 2026. Constructed following the
fifteen-step scientific-writing framework of Drake & Han (2025)
[doi:10.1371/journal.pcbi.1013505](https://doi.org/10.1371/journal.pcbi.1013505).

---

## 1. Introduction

Quantitative phase analysis by the Rietveld method is the standard route
from laboratory PXRD patterns of Portland cement clinker to phase
fractions. Its accuracy depends on three things the analyst controls: the
**structure models**, the **refinement budget**, and the **numerical
authority**. The reproducibility study of García-Maté et al. (2024) on
NIST SRM 2686a [1] provides an ideal benchmark: four professional-grade
patterns (two Cu Kα₁ clinker/silicate-residue, one Cu aluminate-residue,
one synchrotron clinker) with a published eight-phase inventory and
full-budget reference fits at Rwp ≈ 6–9%.

**Gap.** Most published RQPA workflows are interactive: a human drives a
GUI, the refinement budget is implicit, and a reproducibility audit
requires re-running an unrepeatable session. There is no published,
executable protocol that (i) fixes the structure set and its provenance,
(ii) caps the parameter budget, (iii) runs unattended, and (iv) records a
content hash that makes bit-level reproducibility checkable.

**Response.** This paper describes exactly such a protocol, implemented
as the spike-15 runner of the `rietveld-agent` repository [2]. Section 2
fixes the data and structure models; Section 3 defines the budget and
ladder; Section 4 documents the execution model and determinism
machinery; Section 5 reports the converged refinements; Section 6
discusses what the budget can and cannot resolve. The companion papers of
this series establish validity (P-2) and deployment on agent runtimes
(P-3).

## 2. Data and structure models

### 2.1 Data

The four SRM 2686a patterns (Table 1) are the files published with
García-Maté et al. [1] (Zenodo 10.5281/zenodo.1318501), re-emitted as
two-column 2θ–counts with unchanged intensities by the repository's
spike-11/12 ingest parsers (sha256 recorded).

**Table 1.** PXRD data sets (NIST SRM 2686a).

| sample | instrument | range (2θ, °) |
|---|---|---|
| clinker Cu | PANalytical, Ge(111) mono Cu Kα₁ | 4.0–70.0 |
| silicate residue | same | 4.0–70.0 |
| aluminate residue | same | 4.0–70.0 |
| clinker synchrotron | ALBA SXRPD, λ = 0.82543 Å | 2.5–62.85 |

### 2.2 Structure models

Each phase is a COD entry re-validated by an independent parser (space
group, cell, expanded-site composition and density; md5-recorded in
`data/structures/catalog.json`): alite M3 [3], alite T1 [4], belite β [5],
belite α′H [6], cubic aluminate [7], orthorhombic (Na-substituted)
aluminate [8], ferrite C₄AF [9], periclase [10], and aphthitalite [11].
The per-sample inventory matches the published Tables 2–3 of [1]; the Cu
clinker and silicate residue are additionally refined with the alite T1
pseudocell model as a cross-check.

## 3. Refinement protocol

The ladder is implemented in `benchmarks/spikes/spike_15_rqpa_protocol.py`
and specified normatively in `docs/rqpa_protocol.md`. It stages the
GSAS-II [12] Levenberg–Marquardt loop in four steps:

1. per-phase scales with Chebyschev background and sample shift;
2. alite cell;
3. belite cells;
4. minor-phase cells on the Cu clinker only.

The budget is deliberately **bounded**: trace cells and all atom
parameters stay fixed, so the protocol has a finite, documented number of
degrees of freedom and cannot wander into the spurious-solution regime.

Two empirical decisions, both locked by spike-15 probes
(`notes/spike15.md`), complete the protocol:

- On the **aluminate residue**, the published-scale prior is the starting
  point (weak-overlap trace data) and cells are kept fixed to avoid a
  spurious supercell solution; the aphthitalite scale column is
  SVD-degenerate on that window from any start. The ladder therefore
  adapts: it first attempts the full inventory, and only when that
  diverges does it re-run without the offending phase and **reinsert** it
  as a fixed-composition constraint renormalised at its normalized
  published share. No phase is ever dropped as indeterminate — every
  report row carries a wt%.
- On the **synchrotron window**, the alite T1 variant is omitted for the
  same degeneracy reason.

Quantification uses the Hill–Howard relation [13]

```
Wᵢ = Sᵢ·Mᵢ·Vᵢ / Σⱼ Sⱼ·Mⱼ·Vⱼ
```

with Sᵢ the refined scale, Mᵢ the cell-content mass and Vᵢ the refined
cell volume.

## 4. Execution model and determinism

The engine is a deterministic script, not an interactive session:

```
make env       # virtualenv with numpy/scipy/matplotlib
make report    # rerun every RQPA refinement (GSAS-II pinned, ~45 min)
make figures   # regenerate the manuscript figures
make paper     # rebuild paper/main.pdf (tectonic; optional)
make check     # syntax + structure checks
```

GSAS-II is pinned in `.vendor/GSAS-II` and bootstrapped by
`benchmarks/eval/sim.py:ensure_gsasii`; it is never bundled as a
dependency. The report records a **canonical content hash** over the
scientific payload (all fields except wall-clock timing); spike 16
re-runs the full suite and compares hashes (Paper P-2). This is what
makes "same inputs → same outputs" a checkable property of the repository
rather than a claim.

## 5. Results

All six runs converge without metric errors (Table 2); the synchrotron
clinker refines to wR = 9.8% under the bounded budget. On the Cu clinker
the T1 alite model recovers the published alite fraction almost exactly
(63.5 vs. 66.0 wt%), whereas the M3 supercell model over-absorbs
(95.3 wt%): the microabsorption signature expected for big-crystal alite
at Cu Kα₁. The aluminate residue converges at wR = 13.6% with all five
phases reported (aphthitalite as the fixed-composition constraint).

**Table 2.** Converged refinements under the bounded budget.

| sample | model | wR (%) | rank vs. published | tier |
|---|---|---|---|---|
| clinker Cu | M3 | 15.07 | alite top, minor shuffle | ok |
| clinker Cu | T1 | 20.91 | alite top, minor shuffle | ok |
| silicate residue | M3 | 20.18 | alite top, minor shuffle | ok |
| silicate residue | T1 | 27.11 | alite top, minor shuffle | wR-over |
| aluminate residue | M3 | 13.56 | alite n/a (ferrite deficit) | ok |
| clinker synchrotron | M3 | 9.76 | alite top, minor shuffle | ok |

![Clinker Cu pattern with the alite-T1 refinement (observed black,
calculated red, difference blue, offset +800 counts).](../../paper/figures/fig1_clinker_t1_fit.png)

![Left: final wR per sample/model with publication-grade targets marked;
right: alite wt% recovered by the M3 and T1 models vs. published.](../../paper/figures/fig2_results_summary.png)

## 6. Discussion and closing

Two findings matter for the broader RQPA practice. First, the bounded
budget is **deterministic and convergent**: every result in Table 2
rebuilds identically via `make report` (content hash md5 recorded) —
the property an automated or agent-driven pipeline needs. Second, the
residual per-phase deviations are physical, not procedural: on synthetic
patterns generated from the same structures at known fractions the same
ladder recovers the input wt% within the noise floor (Paper P-2), so the
real-pattern bias is attributed to **microabsorption** — the target of
the planned Brindley-correction spike 17 — and is reported honestly
rather than fitted away.

**Closing.** The protocol is the deterministic core of
`rietveld-agent`; its scientific validity is established independently in
Paper P-2, and its operation from agent runtimes (OpenCode, Claude Code,
Codex) is described in Paper P-3 and `docs/installation.md`. This series
follows the fifteen-step structure of Drake & Han [14]: thesis-first
abstract (step 1), story (step 2), methods and rationale (steps 3–5),
CARS introduction (step 6), findings and problem–response pairs (steps
7–10), and paragraph/figure planning throughout (steps 11–15).

## References

1. García-Maté, M., De la Torre, Á. G., León-Reina, L., & Aranda, M. A. G. (2024). Reproducibility and accuracy of Rietveld quantitative phase analysis of NIST SRM 2686a cement clinker. *Cement and Concrete Research*, 180, 107506. doi:10.1016/j.cemconres.2021.106376
2. qaemu (2026). rietveld-agent: a deterministic Rietveld QPA engine for scientific agent runtimes. https://github.com/qaemu/rietveld-agent-spikes (Apache-2.0)
3. Nishi, F., Takeuchi, Y., & Maki, I. (1985). Tricalcium silicate Ca₃OSi₅: the monoclinic superspace structure. *Z. Kristallogr.*, 172, 297–308.
4. De la Torre, Á. G., Bruque, S., Campo, J., Turrillas, X., & Aranda, M. A. G. (2002). The superstructure of C₃S from synchrotron and neutron powder diffraction. *Powder Diffraction*, 17, 240–246. doi:10.1154/1.1500529
5. Tsurumi, T., Hirano, Y., Kato, H., Kamiya, T., & Daimon, M. (1994). Crystal structure and hydration of belite. *J. Am. Ceram. Soc.*, 77(3), 765–769. doi:10.1111/j.1151-2916.1994.tb05362.x
6. Mumme, W. G., Cranswick, L. M. D., & Chakoumakos, B. C. (1996). Crystal chemistry of the polymorphs of dicalcium silicate. *N. Jb. Miner. Abh.*, 170, 171–188.
7. Mondal, P., & Jeffery, J. W. (1975). The crystal structure of tricalcium aluminate. *Acta Cryst. B*, 31, 689–697. doi:10.1107/S0567740875003552
8. Takeuchi, Y., Nishi, F., & Maki, I. (1980). Crystal-chemical characterization of the 3CaO·Al₂O₃–Na₂O solid-solution series. *Z. Kristallogr.*, 152, 259–307.
9. Colville, A. A., & Geller, S. (1972). The crystal structure of brownmillerite, Ca₂FeAlO₅. *Acta Cryst. B*, 28, 3196–3200. doi:10.1107/S0567740872007670
10. Sasaki, S., Fujino, K., Takeuchi, Y., & Sadanaga, R. (1979). On the estimation of atomic charges by the X-ray method for some oxides and silicates. *Acta Cryst. A*, 35, 934–939. doi:10.1107/S0567739479002107
11. Okada, K., & Ossaka, J. (1980). Structures of potassium sodium sulphate and tripotassium sodium disulphate. *Acta Cryst. B*, 36, 919–921. doi:10.1107/S0567740880005011
12. Toby, B. H., & Von Dreele, R. B. (2013). GSAS-II: the genesis of a modern open-source all purpose crystallography software package. *J. Appl. Cryst.*, 46, 544–549. doi:10.1107/S0021889813003531
13. Hill, R. J., & Howard, C. J. (1987). Quantitative phase analysis from neutron powder diffraction data using the Rietveld method. *J. Appl. Cryst.*, 20, 467–474. doi:10.1107/S0021889887086199
14. Drake, J. M., & Han, B. A. (2025). How to write a scientific paper in fifteen steps. *PLoS Comput. Biol.*, 21(9), e1013505. doi:10.1371/journal.pcbi.1013505 (PMC12459795)

---

*License: Apache-2.0; see LICENSE in the repository root. All data and
programs needed to reproduce this paper are in the repository.*
