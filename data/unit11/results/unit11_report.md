# Unit 11 - real multiphase inorganics: NIST SRM 2686a clinker

Source: **10.5281/zenodo.1318501** (https://zenodo.org/records/1318501), CC-BY-4.0 -- García-Maté et al., 'Rietveld Quantitative Phase Analyses of SRM 2686a: a Standard Portland Clinker' (rev. ms, Cement and Concrete Research).

Four real multiphase inorganic patterns (Portland cement clinker and its selective-dissolution residues) through `cli analyze`.
Native XRDML 1.3 / ALBA .dat -> XRDML 1.0 container: format conversion only, intensities unchanged (sha256 + md5 verified).

| case | radiation | peak_max | verdict | top-ranked family | published top-5 |
|---|---|---|---|---|---|
| clinker_cu | CuKa1-LXRPD | 1.594e+04 | abstain | (abstain) CaCO3 (calcite) (0.2114) | alite (C3S) … |
| residue_silicate_cu | CuKa1-LXRPD | 1.473e+04 | abstain | (abstain) CaCO3 (calcite) (0.2472) | alite (C3S) … |
| residue_aluminate_cu | CuKa1-LXRPD | 1.135e+04 | abstain | (abstain) PbSO4 (anglesite) (0.2566) | ferrite (C4AF) … |
| clinker_sync | SXRPD λ=0.82543 Å | 5.309e+05 | abstain | (abstain) PbSO4 (anglesite) (0.2203) | alite (C3S) … |

## Verdicts vs published inventory
The candidate library (12 minerals) contains NO cement phases: alite, belite, aluminate and ferrite are not entries. The governed verdict therefore abstains / reports low confidence instead of forcing a wrong single-phase identity -- this is the designed anti-hallucination behaviour. The top-5 ranking is reported against the published Rietveld QPA inventory for each sample (tables in the JSON report).

## Sources
- Zenodo record: https://zenodo.org/records/1318501 (DOI 10.5281/zenodo.1318501)
- File: https://zenodo.org/records/1318501/files/Clinker_Nist_CuKalpha1_R1.xrdml
- File: https://zenodo.org/records/1318501/files/Silicate_enriched_residue_Nist_CuKalpha1_R1.xrdml
- File: https://zenodo.org/records/1318501/files/aluminate_enriched_residue_clinkerNIST_180718_R1.xrdml
- File: https://zenodo.org/records/1318501/files/Clinker_Synchrotron.dat
- Paper: García-Maté et al., rev. ms v3 submitted to Cement and Concrete Research (methods: instrument details, RQPA Tables 2-3).