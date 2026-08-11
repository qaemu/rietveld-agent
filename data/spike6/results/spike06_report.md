# Spike 06: COD catalog release 0.1.1 + catalog-backed library

- 12 validated entries (release `/Users/lucas/Documents/CNEA/rietveld-agent/data/catalog/releases/catalog_0.1.1.json`), 0 rejected, cif_validation_rate **1.0**.
- Newly fetched from COD: 9 (first-introduced by spike 06: all but the three pre-existing fixture CODs [1000041, 1009000, 9011025]).
- Runtime library rebuilt from the release: all `cod-<id>` materials catalog-backed; legacy M1 fixture entries provenance-corrected.

## Entries

| cod_id | family | formula | space group | flags |
|---|---|---|---|---|
| 1000017 | Al2O3 (corundum) | Al2 O3 | R -3 c :H | - |
| 1000041 | NaCl (halite) | Cl Na | F m -3 m | - |
| 1000043 | CaF2 (fluorite) | Ca F2 | F m -3 m | - |
| 1009000 | GaAsO4 (quartz homeotype) | As Ga O4 | P 31 2 1 | quartz-homeotype, ex-SiO2-mislabel |
| 1010928 | CaCO3 (calcite) | C Ca O3 | R -3 c RS | - |
| 1010950 | PbSO4 (anglesite) | O4 Pb S | P n m a | - |
| 1530150 | TiO2 (rutile) | O2 Ti | P 42/m n m | - |
| 1544891 | FeS2 (pyrite) | Fe S2 | P a -3 | - |
| 9000927 | magnetite-family | Fe2.75 O4 Ti0.25 | F d -3 m :2 | ti-substituted |
| 9006758 | MgO (periclase) | Mg O | F m -3 m | - |
| 9009666 | SiO2 (quartz) | O2 Si | P 31 2 1 | - |
| 9011025 | AgBiS2 (schapbachite) | Ag0.5 Bi0.5 S | F m -3 m | synthetic, ex-NaCl-mislabel |

## Rejected


## M1 library provenance audit

| m1 entry | source cif | claimed | actual | action |
|---|---|---|---|---|
| mat-nacl | halite_9011025.cif | NaCl (halite) | Ag0.5Bi0.5S (schapbachite, NaCl-type lattice) | replaced by cod-1000041 (real NaCl); 9011025 becomes its own schapbachite entry |
| mat-sio2 | quartz_1009000.cif | SiO2 (quartz-family) | GaAsO4 (quartz homeotype, P3121) | replaced by cod-1009000 (GaAsO4 homeotype) + real SiO2 quartz entry |
| mat-pbso4 / mat-pbso4-fe | PbSO4-Wyckoff.cif | PbSO4 | PbSO4 (anglesite) -- verified | kept; catalog reference = COD anglesite |

## Release 0.1.1 structural corrections

| entry | family | actual structure | evidence | action |
|---|---|---|---|---|
| cod-1010942 (release 0.1.0) | TiO2 (rutile) | ANATASE, I41/amd, a=3.73 c=9.37 | published rutile R040049 cell a=4.5955(1) c=2.9598(1) (RRUFF REFINE v3.0) | replaced by cod-1530150 (O2 Ti, P42/mnm, a=4.59 c=2.96; Khitrova et al. 1977, Kristallografiya 22) |
| cod-9012601 (release 0.1.0) | SiO2 (quartz) | compressed quartz variant, P3121, a=4.812 c=5.327 | real alpha-quartz cell a~4.913 c~5.405; RRUFF quartz R040031 refined a=4.9134 c=5.4042 | pool now prefers cod-9009666 (alpha-quartz, a=4.9158 c=5.4091); 9012601 no longer a library material |

- Library manifest: 11d12a30671362b8... (deterministic rebuild: True)
- Library materials: 14 (cod-1000017, cod-1000041, cod-1000043, cod-1009000, cod-1010928, cod-1010950, cod-1530150, cod-1544891, cod-9000927, cod-9006758, cod-9009666, cod-9011025, mat-pbso4, mat-pbso4-fe)
- Wall clock: 11.1s

## Verdict
- Release schema-valid, entries chemistry-validated, rejected list recorded.
- Release 0.1.1 structural audit: the rutile family now hosts the true rutile structure (1530150, P42/mnm, anterior anatase 1010942 removed) and the quartz family the true alpha-quartz (9009666, compressed 9012601 removed) -- both match the RRUFF refined cells (R040049 / R040031) within profile tolerance, so published real rutile/quartz patterns can now identify.
- M1 library corrected: no remaining mislabelled chemistry (PbSO4 verified; NaCl/SiO2 fixed).
- Spike 04/05 re-run on the catalog-backed library: e2e green (PbSO4 (anglesite) same-anode + cross-anode; cu_quartz resolves to GaAsO4 (quartz homeotype); NaCl eval now simulates the real halite COD 1000041; family-aware margins >= 0.9).
