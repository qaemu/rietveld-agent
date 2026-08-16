# SRM 2686a RQPA structure set

Canonical structure models for the Rietveld quantitative phase analysis of
the NIST SRM 2686a Portland cement clinker, matching the polymorph set used
by the published reference analysis (Garcia-Mate et al. 2024, *Cem. Concr.
Res.* 180, 107506, DOI 10.1016/j.cemconres.2021.106376, Table 1).

## Files

| COD id | phase | formula (cell) | space group | cell (A, deg) | rho (g/cm3) | source |
|---|---|---|---|---|---|---|
| 9008366 | alite M3 | Ca108 Si36 O180 | C 1 m 1 (#8) | 33.083 7.027 18.499, b=94.12 | 3.182 | Nishi, Takeuchi & Maki 1985, Z. Kristallogr. 172, 297 (AMCSD 0010938) |
| 1538413 | alite T1 (cross-check) | Ca36 Si12 O60 | P -1 (#2) | 11.67 14.24 13.72, a=105.3 b=94.5 g=90.0 | 3.130 | de la Torre et al. 2002, Powder Diffr. 17, 240; triclinic T1 pseudocell (unit-14 T1 variant, unit-15 cross-check model) |
| 9012794 | belite beta (larnite) | Ca8 Si4 O16 | P 1 21/n 1 (#14) | 5.5075 6.7508 9.3054, b=94.59 | 3.317 | Tsurumi et al. 1994 (larnite, neutron) |
| 1546027 | belite alpha'H | Ca8 Si4 O16 | P n m a (#62) | 6.8709 5.6010 9.5563 | 3.111 | Mumme et al. 1996 (alpha'H model of Mumme 1995) |
| 1000039 | aluminate cubic (C3A) | Ca72 Al48 O144 | P a -3 (#205) | 15.263 15.263 15.263 | 3.028 | Mondal & Jeffery 1975, Acta Cryst. B31, 689 |
| 8103596 | aluminate orthorhombic (C3A-Na) | Ca33.57 Al20.70 Fe1.8 Na3.43 Si1.50 O72 | P b c a (#61) | 10.879 10.845 15.106 | 3.054 | Takeuchi, Nishi & Maki 1980, Z. Kristallogr. 152, 259 (NIST C3A1 model) |
| 1200009 | ferrite C4AF | Ca8 Al4 Fe4 O20 | I b m 2 (#46) | 5.584 14.60 5.374 | 3.684 | Colville & Geller 1972, Acta Cryst. B28, 3196 |
| 1000053 | periclase | Mg4 O4 | F m -3 m (#225) | 4.217 4.217 4.217 | 3.57 | Sasaki et al. 1979 (electron density) |
| 9007639 | aphthitalite | K3 Na1 S2 O8 | P -3 m 1 (#164) | 5.6801 5.6801 7.309 | 2.703 | Okada & Ossaka 1980, Acta Cryst. B36, 919 |

All files are CC0 from the Crystallography Open Database
(https://www.crystallography.net/cod/). `catalog.json` records the md5 of
every file and the full per-entry validation (cell, space group, operator
set used, cell composition, density).

## Origin and substitutions

The reference analysis used the structure files bundled with the NIST GSAS
tutorial (`Cements_Data.zip`, concrete.nist.gov/~bullard/, referenced in
Garcia-Mate et al. 2024). That host is dead and not archived (checked 2026);
each structure is therefore re-sourced from COD as the closest verified
publication match:

* **alite M3** -- the tutorial's ALITE_M3 model corresponds to the
  monoclinic M3 supercell of Nishi/Takeuchi/Maki (1985), the same cell
  (33.08/7.03/18.50 A, beta=94.12) refined later by de la Torre et al.
  (2002). rho = 3.182 g/cm3 agrees with the CIF.
* **beta-belite** -- Mumme et al. (1995) (ICSD 81096) is not in COD; the
  Tsurumi et al. (1994) larnite structure (same P 1 21/n 1 polymorph,
  neutron data) is used instead.
* **alpha'H-belite** -- Mumme et al. (1996) Pnma entry, the alpha'H model
  of Mumme (1995) republished for the neutron polymorph study.
* **orthorhombic aluminate** -- Takeuchi/Nishi/Maki (1980) Pbca
  Na-substituted C3A, the same model behind the NIST C3A1 file.
* **ferrite** -- Colville & Geller (1972) Ibm2, the standard C4AF model.
* **aphthitalite** -- Okada & Ossaka (1980) P-3m1. The CIF's operator loop
  is actually a complete International-Tables P-3m1 set (verified op-by-op
  in unit 15); the unusual bits are the Hall symbol string (`-P 3 2"`)
  and the non-canonical op ordering, which is why the unit-14 parser
  validates it with the International Tables set (`FORCE_HARD`). GSAS-II
  imports it correctly (SGData P-3m1, trigonal, 3m1 Laue), but its Scale
  column is numerically degenerate on the aluminate-residue data window
  (see unit 15 notes); clinker/sync patterns refine it normally.

## Validation

`benchmarks/protocols/unit_14_structures.py` re-validates every file:

1. space group number matches the literature polymorph;
2. cell constants match the published values;
3. the site list expanded over the space-group operators (special
   positions deduped, occupancies applied) reproduces the expected per-cell
   composition, cross-checked via the crystallographic density (rho within
   1% of the CIF/literature value);
4. compositions: alite Ca:Si:O = 3:1:5, belite 2:1:4, C3A cub/ort
   Ca:Al:O = 3:2:6, ferrite Ca:Al:Fe:O = 2:1:1:5, periclase 1:1,
   aphthitalite K:Na:S:O = 3:1:2:8.

Run: `.venv/bin/python benchmarks/protocols/unit_14_structures.py`
