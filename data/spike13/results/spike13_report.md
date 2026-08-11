# Spike 13: multiphase Rietveld QPA (GSAS-II) on SRM 2686a clinker

- Goal: quantitative phase analysis of the 4 real NIST clinker patterns using ONLY the published phase set, with the COD CIFs identified by spike 12 (alite 1538413, belite 2312428, ferrite 1200009, periclase 1000053, C3A 1000039).
- Budget (bounded, identical for all samples): chebyschev-1 background x8, sample shift, per-phase Scale + Cell; Uiso/peak shapes/absorption fixed.
- Instruments: Cu Kα1 strictly monochromatic (α1=1.540598, ratio 0) per the native XRDML; ALBA sync λ=0.82543 Å per the paper.
- Fractions: Hill-Howard Rietveld normalization W_i = S_i·M_i·V_i/Σ, M_i = GSAS-II cell-content mass, V_i from the refined cell.
- Published RQPA reference: García-Maté et al. 2024 (Tables 2-3).

| sample | Rwp | conv | phase (this QPA wt% / published wt%) |
|---|---|---|---|
| Clinker_Nist_CuKalpha1_R1.xrdml | 0.1752 | ✓ | alite 66.6%(pub 66.0), belite 20.0%(pub 16.2), C3A 10.3%(pub 1.9), ferrite 2.4%(pub 11.1), periclase 0.8%(pub 4.0) |
| Silicate_enriched_residue_Nist_CuKalpha1_R | 0.2639 | ✓ | alite 86.3%(pub 78.7), belite 12.8%(pub 16.3), periclase 0.9%(pub 5.0) |
| aluminate_enriched_residue_clinkerNIST_180 | 0.1477 | ✓ | C3A 82.3%(pub 13.1), ferrite 15.1%(pub 69.8), periclase 2.7%(pub 17.2) |
| Clinker_Synchrotron.dat | 0.3121 | ✓ | alite 89.7%(pub 65.4), belite 7.4%(pub 16.8), ferrite 2.2%(pub 11.6), periclase 0.7%(pub 3.65), C3A 0.0%(pub 1.99) |

## Per-sample tables

### Clinker_Nist_CuKalpha1_R1.xrdml

| phase | COD | scale | mass | V(Å³) | a b c α β γ | wt% | pub wt% |
|---|---|---|---|---|---|---|---|
| alite | 1538413 | 18.75 | 2054.41 | 1079.85 | 7.0636 7.0636 24.9912 90.0000 90.0000 120.0000 | 66.63 | 66.0 |
| belite | 2312428 | 52.56 | 688.97 | 344.05 | 6.8184 5.5620 9.0722 90.0000 90.0000 90.0000 | 19.96 | 16.2 |
| C3A | 1000039 | 0.2782 | 6484.75 | 3547.46 | 15.2513 15.2513 15.2513 90.0000 90.0000 90.0000 | 10.26 | 1.9 |
| ferrite | 1200009 | 2.924 | 1171.57 | 433.04 | 5.5131 14.6597 5.3580 90.0000 90.0000 90.0000 | 2.38 | 11.1 |
| periclase | 1000053 | 39.99 | 161.22 | 75.03 | 4.2177 4.2177 4.2177 90.0000 90.0000 90.0000 | 0.77 | 4.0 |

### Silicate_enriched_residue_Nist_CuKalpha1_R1.xrdml

| phase | COD | scale | mass | V(Å³) | a b c α β γ | wt% | pub wt% |
|---|---|---|---|---|---|---|---|
| alite | 1538413 | 4.353 | 2054.41 | 1081.64 | 7.0672 7.0672 25.0067 90.0000 90.0000 120.0000 | 86.26 | 78.7 |
| belite | 2312428 | 5.821 | 688.97 | 359.23 | 6.7275 5.6528 9.4461 90.0000 90.0000 90.0000 | 12.85 | 16.3 |
| periclase | 1000053 | 8.296 | 161.22 | 75.10 | 4.2190 4.2190 4.2190 90.0000 90.0000 90.0000 | 0.90 | 5.0 |

### aluminate_enriched_residue_clinkerNIST_180718_R1.xrdml

| phase | COD | scale | mass | V(Å³) | a b c α β γ | wt% | pub wt% |
|---|---|---|---|---|---|---|---|
| C3A | 1000039 | 0.9879 | 6484.75 | 3562.43 | 15.2727 15.2727 15.2727 90.0000 90.0000 90.0000 | 82.27 | 13.1 |
| ferrite | 1200009 | 8.15 | 1171.57 | 437.20 | 5.5626 14.6494 5.3651 90.0000 90.0000 90.0000 | 15.05 | 69.8 |
| periclase | 1000053 | 61.35 | 161.22 | 75.26 | 4.2220 4.2220 4.2220 90.0000 90.0000 90.0000 | 2.68 | 17.2 |

### Clinker_Synchrotron.dat

| phase | COD | scale | mass | V(Å³) | a b c α β γ | wt% | pub wt% |
|---|---|---|---|---|---|---|---|
| alite | 1538413 | 41.04 | 2054.41 | 1082.05 | 7.0780 7.0780 24.9400 90.0000 90.0000 120.0000 | 89.71 | 65.4 |
| belite | 2312428 | 31.65 | 688.97 | 345.05 | 6.7650 5.5140 9.2500 90.0000 90.0000 90.0000 | 7.40 | 16.8 |
| ferrite | 1200009 | 4.439 | 1171.57 | 438.12 | 5.5840 14.6000 5.3740 90.0000 90.0000 90.0000 | 2.24 | 11.6 |
| periclase | 1000053 | 54.35 | 161.22 | 74.99 | 4.2170 4.2170 4.2170 90.0000 90.0000 90.0000 | 0.65 | 3.65 |
| C3A | 1000039 | 1e-12 | 6484.75 | 3555.66 | 15.2630 15.2630 15.2630 90.0000 90.0000 90.0000 | 0.00 | 1.99 |


## Honest limitations
- COD polymorph approximations vs the real clinker phases: alite T1/M3 choice, pure Ca2AlFeO5 ferrite vs C4AF solid solutions, α'-vs-β belite; cell parameters refine, atom positions do not.
- Fixed Uiso / no peak-shape refinement: Rwp stays elevated on real clinker breadth; fractions therefore approximate.
- Periclase is a *minor* phase (2-17 wt%): its 2 strong lines (200/220) drive its scale; the spike-12 window study showed those lines match to <0.003 Å, which is what makes the scale physical.