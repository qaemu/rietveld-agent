# Unit 12: whole-COD matching (with any file given)

- Goal: match arbitrary powder files against the **entire** COD 
  (216,262 entries, 38,943,380 allowed d-lines, d ∈ [1.1, 22.0] Å).
- Index: space-group extinction (vectorized via gemmi symops), position-only (no cell intensities).
- Screening window: ±0.02 Å per fingerprint peak; significance score vs the entry's own line density.
- Two-stage verification: (1) fast kinematic intensity calc from the COD CIF (gaussian-atom form factors, exact orbit/absorption via the site sum) ranks all screen survivors by coverage + intensity correlation; (2) the top 8 per sample plus the top 5 by screen significance are confirmed with full GSAS-II CIF simulation (unit06 protocol, Cu Kα1/α2) on the shared d-space grid.
- Curated ground-truth anchor phases (periclase, brownmillerite, belite/larnite, C3A, lime) are force-included in the GSAS stage and tagged [anchor] in the tables, so the verifier is also tested on the real minor clinker phases even if screening or coverage ranking dropped them.
- Sources: Crystallography Open Database (CC0); metadata via REST `id=%` export.

## Clinker_Nist_CuKalpha1_R1.xrdml
- fingerprint: 29 peaks (λ=1.54060 Å); published phases: alite C3S / beta-belite C2S / ferrite C4AF / periclase MgO / alpha'H-belite / cub+ortho-aluminate C3A / aphthitalite

| # | COD | mineral | formula | SG | yr | fast cv/cr | simII | match |
|---|-----|---------|---------|----|----|------------|-------|-------|
| 1 | 1530012 | Te2 O (O H)2 F4 | F4 H2 O3 Te2 | F 2 d d | 1982 | 0.931/0.5538 | 0.2093 | 18/29 |
| 2 | 1538413 | Ca3 (Si O4) O | Ca3 O5 Si | R 3 m :H | 1950 | 0.8966/0.5472 | 0.5459 | 17/29 |
| 3 | 2310872 | Ca3 (Si O4) O | Ca3 O5 Si | R 3 m :H | 1952 | 0.8966/0.4716 | 0.2961 | 15/29 |
| 4 | 1525851 | Te2 V2 O9 | O9 Te2 V2 | F d d 2 | 1973 | 0.8621/0.3627 | 0.1906 | 14/29 |
| 5 | 1520006 | COD 1520006 | Bi9 O7.5 S6 | R -3 m :H | 2015 | 0.7931/0.6716 | 0.1967 | 13/29 |
| 6 | 4104493 | catena-(dimethylammonium)-(tris(μ~2~-Formato-O,O')-iron(ii) | C5 H15 Fe N O6 | R -3 c :H | 2009 | 0.8621/0.1185 | 0.194 | 13/29 |
| 7 | 2312428 | COD 2312428 **[anchor]** | Ca8 O16 Si4 | P n m a | 2024 | 0/0 | 0.3286 | 12/29 |
| 8 | 1539095 | Ca.80 Bi3.80 O4 Cl5 | Bi3.8 Ca0.8 Cl5 O4 | I 4/m m m | 1941 | 0.8276/0.5786 | 0.2213 | 12/29 |
| 9 | 1539099 | Cd0.5 Bi4 O4 Cl5 | Bi4 Cd0.5 Cl5 O4 | I 4/m m m | 1941 | 0.8966/0.6401 | 0.2749 | 11/29 |
| 10 | 1525036 | COD 1525036 | Co7 Er2 | R -3 m :H | 1967 | 0.7931/0.351 | 0.2642 | 11/29 |
| 11 | 1200009 | Brownmillerite **[anchor]** | Al Ca2 Fe O5 | I b m 2 | 1971 | 0/0 | 0.2272 | 11/29 |
| 12 | 1522366 | COD 1522366 | Gd2 Ni7 | R -3 m :H | 1969 | 0.7586/0.2821 | 0.1996 | 11/29 |
| 13 | 2002081 | COD 2002081 | Cr0.5 K0.5 S2 Sn0.5 | R -3 m :H | 1993 | 0.6897/0.5404 | 0.3063 | 10/29 |
| 14 | 1518102 | COD 1518102 | Na2 O11 Ta4 | R -3 c :H | 1985 | 0.8966/0.6814 | 0.2233 | 10/29 |
| 15 | 1522301 | Bi2 Sr2 Nb2 Mn O11.57 | Bi2 Mn Nb2 O11.57 Sr2 | F m m m | 1999 | 0.8966/0.4895 | 0.1731 | 10/29 |
| 16 | 2106554 | (Yb0.5 Eu0.5) Fe2 O4 | Eu0.5 Fe2 O4 Yb0.5 | R -3 m :H | 1975 | 0.7586/0.5807 | 0.164 | 10/29 |
| 17 | 2002673 | Tricopper magnesium antimony oxide (3/1/1.4/6) | Cu3 Mg O6 Sb1.4 | R 3 c :H | 1995 | 0.7931/0.2016 | 0.1633 | 10/29 |
| 18 | 1526132 | Er (Fe Mn) O4 | Er Fe Mn O4 | R -3 m :H | 2000 | 0.7586/0.4578 | 0.1594 | 10/29 |
| 19 | 2234027 | Didysprosium heptanickel | Dy2 Ni7 | R -3 m :H | 2012 | 0.7586/0.3158 | 0.2648 | 9/29 |
| 20 | 1521199 | Tl1.772 (Ba0.80 Sr0.20)2 Ca2 Cu3 O9.752 | Ba1.6 Ca2 Cu3 O9.752 Sr0.4 Tl1.772 | I 4/m m m | 1999 | 0.8966/0.6009 | 0.2438 | 9/29 |
| 21 | 1532574 | (Tl0.889 Cu0.067)2 Ba2 Ca2 Cu3 O9.808 | Ba2 Ca2 Cu3.134 O9.808 Tl1.778 | I 4/m m m | 2003 | 0.8966/0.5385 | 0.2011 | 9/29 |
| 22 | 1561710 | COD 1561710 | F1.5 Na4.5 O6 S1.5 | R -3 m :H | 2015 | 0.8621/0.4397 | 0.1584 | 9/29 |
| 23 | 1521200 | Tl1.810 (Ba0.75 Sr0.25)2 Ca2 Cu3 O9.784 | Ba1.5 Ca2 Cu3 O9.784 Sr0.5 Tl1.81 | I 4/m m m | 1999 | 0.8966/0.598 | 0.2547 | 8/29 |
| 24 | 1521198 | Tl1.778 (Ba0.85 Sr0.15)2 Ca2 Cu3 O9.728 | Ba1.7 Ca2 Cu3 O9.728 Sr0.3 Tl1.778 | I 4/m m m | 1999 | 0.8966/0.5531 | 0.2356 | 8/29 |
| 25 | 1573840 | Di(furan-2-yl)ethanedione | C10 H6 O4 | F d d 2 | 2025 | 0.931/0.5864 | 0.1485 | 7/29 |
| 26 | 1524402 | In Cu Al O4 | Al Cu In O4 | R -3 m :H | 1980 | 0.6207/0.338 | 0.1267 | 7/29 |
| 27 | 1000039 | Calcium cyclo-hexaaluminate **[anchor]** | Al6 Ca9 O18 | P a -3 | 1975 | 0/0 | 0.1593 | 6/29 |
| 28 | 1540838 | (Y3 Zr2)0.4 | Y1.2 Zr0.8 | P 63/m m c | 1972 | 0.4828/0.468 | 0.1526 | 6/29 |
| 29 | 1538059 | Ti Se | Se Ti | P 63/m m c | 1959 | 0.5172/0.356 | 0.2671 | 5/29 |
| 30 | 1538290 | V Se | Se V | P 63/m m c | 1939 | 0.5172/0.5072 | 0.2452 | 5/29 |
| 31 | 2300518 | COD 2300518 | Tl | P 63/m m c | 1994 | 0.4828/0.445 | 0.1412 | 5/29 |
| 32 | 1528102 | (Sn0.017 Tl0.983) | Sn0.017 Tl0.983 | P 63/m m c | 1960 | 0.4828/0.445 | 0.1406 | 5/29 |
| 33 | 2002599 | Tribarium tetradysprosium oxide | Ba3 Dy4 O9 | R 3 :H | 1993 | 0.8276/0.3433 | 0.1238 | 5/29 |
| 34 | 1000053 | Periclase **[anchor]** | Mg O | F m -3 m | 1979 | 0.2069/0.0571 | 0.0901 | 3/29 |
| 35 | 1563205 | COD 1563205 **[anchor]** | Ca O | F m -3 m | 2011 | 0/0 | 0.0752 | 3/29 |
| 36 | 1522937 | (Fe0.7 Zn0.3) | Fe0.7 Zn0.3 | I m -3 m | 1964 | 0.1379/0.0718 | 0.0327 | 2/29 |
| 37 | 1522976 | COD 1522976 | Fe0.8 Ti0.2 | I m -3 m | 1986 | 0.1379/0.0718 | 0.0326 | 2/29 |
| 38 | 1524987 | (Cr0.75 Ru0.25) | Cr0.75 Ru0.25 | I m -3 m | 1986 | 0.1379/0.0718 | 0.0393 | 1/29 |
| 39 | 1524276 | COD 1524276 | Cr0.7 Mo0.3 | I m -3 m | 1982 | 0.1379/0.0718 | 0.0343 | 1/29 |
| 40 | 1525300 | (Cr0.333 Mo0.333 V0.334) | Cr0.333 Mo0.333 V0.334 | I m -3 m | 1976 | 0.2414/0.0927 | 0.0324 | 1/29 |
| 41 | 1521913 | Na0.288 Hf N Cl | Cl Hf N Na0.288 | R -3 m :H | 1999 | 0.8621/0.6468 | 0.0 | 0/29 |
| 42 | 1001366 | Thallium(I,III) barium calcium copper(III) oxide (1.82/2/1.9/3/10.94) | Ba2 Ca1.9 Cu3 O10.94 Tl1.82 | I 4/m m m | 1988 | 0.8966/0.5372 | 0.0 | 0/29 |

CIF records: [2002081](https://www.crystallography.net/cod/2002081.html)  , [2106554](https://www.crystallography.net/cod/2106554.html)  , [1524276](https://www.crystallography.net/cod/1524276.html)  , [1524987](https://www.crystallography.net/cod/1524987.html)  , [1522976](https://www.crystallography.net/cod/1522976.html)  , [1522937](https://www.crystallography.net/cod/1522937.html)  , [1525300](https://www.crystallography.net/cod/1525300.html)  , [1524402](https://www.crystallography.net/cod/1524402.html)  , [1540838](https://www.crystallography.net/cod/1540838.html)  , [1525036](https://www.crystallography.net/cod/1525036.html)  , [2234027](https://www.crystallography.net/cod/2234027.html)  , [1522366](https://www.crystallography.net/cod/1522366.html)  , [1520006](https://www.crystallography.net/cod/1520006.html)  , [1521199](https://www.crystallography.net/cod/1521199.html)  , [2002673](https://www.crystallography.net/cod/2002673.html)  , [1538290](https://www.crystallography.net/cod/1538290.html)  , [1521913](https://www.crystallography.net/cod/1521913.html)  , [1521200](https://www.crystallography.net/cod/1521200.html)  , [1526132](https://www.crystallography.net/cod/1526132.html)  , [2310872](https://www.crystallography.net/cod/2310872.html)  , [1528102](https://www.crystallography.net/cod/1528102.html)  , [2300518](https://www.crystallography.net/cod/2300518.html)  , [1518102](https://www.crystallography.net/cod/1518102.html)  , [1538059](https://www.crystallography.net/cod/1538059.html)  , [1539099](https://www.crystallography.net/cod/1539099.html)  , [1539095](https://www.crystallography.net/cod/1539095.html)  , [1538413](https://www.crystallography.net/cod/1538413.html)  , [1522301](https://www.crystallography.net/cod/1522301.html)  , [1525851](https://www.crystallography.net/cod/1525851.html)  , [1561710](https://www.crystallography.net/cod/1561710.html)  , [1530012](https://www.crystallography.net/cod/1530012.html)  , [1001366](https://www.crystallography.net/cod/1001366.html)  , [4104493](https://www.crystallography.net/cod/4104493.html)  , [1573840](https://www.crystallography.net/cod/1573840.html)  , [2002599](https://www.crystallography.net/cod/2002599.html)  , [1521198](https://www.crystallography.net/cod/1521198.html)  , [1532574](https://www.crystallography.net/cod/1532574.html)  , [1000053](https://www.crystallography.net/cod/1000053.html)  , [1200009](https://www.crystallography.net/cod/1200009.html)  , [2312428](https://www.crystallography.net/cod/2312428.html)  , [1000039](https://www.crystallography.net/cod/1000039.html)  , [1563205](https://www.crystallography.net/cod/1563205.html)  

## Silicate_enriched_residue_Nist_CuKalpha1_R1.xrdml
- fingerprint: 28 peaks (λ=1.54060 Å); published phases: alite C3S / beta-belite C2S / alpha'H-belite / periclase MgO

| # | COD | mineral | formula | SG | yr | fast cv/cr | simII | match |
|---|-----|---------|---------|----|----|------------|-------|-------|
| 1 | 1530012 | Te2 O (O H)2 F4 | F4 H2 O3 Te2 | F 2 d d | 1982 | 1.0/0.5259 | 0.1921 | 19/28 |
| 2 | 1538413 | Ca3 (Si O4) O | Ca3 O5 Si | R 3 m :H | 1950 | 1.0/0.5784 | 0.6695 | 17/28 |
| 3 | 2310872 | Ca3 (Si O4) O | Ca3 O5 Si | R 3 m :H | 1952 | 0.9643/0.4805 | 0.2557 | 16/28 |
| 4 | 1525851 | Te2 V2 O9 | O9 Te2 V2 | F d d 2 | 1973 | 0.9286/0.3998 | 0.1513 | 14/28 |
| 5 | 2312428 | COD 2312428 **[anchor]** | Ca8 O16 Si4 | P n m a | 2024 | 0/0 | 0.3394 | 13/28 |
| 6 | 1520006 | COD 1520006 | Bi9 O7.5 S6 | R -3 m :H | 2015 | 0.8929/0.7027 | 0.1677 | 13/28 |
| 7 | 4104493 | catena-(dimethylammonium)-(tris(μ~2~-Formato-O,O')-iron(ii) | C5 H15 Fe N O6 | R -3 c :H | 2009 | 0.8929/0.1443 | 0.1463 | 13/28 |
| 8 | 1539099 | Cd0.5 Bi4 O4 Cl5 | Bi4 Cd0.5 Cl5 O4 | I 4/m m m | 1941 | 0.9286/0.6142 | 0.2823 | 11/28 |
| 9 | 1525036 | COD 1525036 | Co7 Er2 | R -3 m :H | 1967 | 0.8571/0.377 | 0.2669 | 11/28 |
| 10 | 1539095 | Ca.80 Bi3.80 O4 Cl5 | Bi3.8 Ca0.8 Cl5 O4 | I 4/m m m | 1941 | 0.8571/0.546 | 0.2067 | 11/28 |
| 11 | 1522366 | COD 1522366 | Gd2 Ni7 | R -3 m :H | 1969 | 0.8214/0.3039 | 0.2037 | 11/28 |
| 12 | 1518102 | COD 1518102 | Na2 O11 Ta4 | R -3 c :H | 1985 | 0.8929/0.6482 | 0.2016 | 11/28 |
| 13 | 1200009 | Brownmillerite **[anchor]** | Al Ca2 Fe O5 | I b m 2 | 1971 | 0/0 | 0.1602 | 11/28 |
| 14 | 2106554 | (Yb0.5 Eu0.5) Fe2 O4 | Eu0.5 Fe2 O4 Yb0.5 | R -3 m :H | 1975 | 0.8214/0.6088 | 0.114 | 11/28 |
| 15 | 2002081 | COD 2002081 | Cr0.5 K0.5 S2 Sn0.5 | R -3 m :H | 1993 | 0.75/0.5523 | 0.2714 | 10/28 |
| 16 | 1522301 | Bi2 Sr2 Nb2 Mn O11.57 | Bi2 Mn Nb2 O11.57 Sr2 | F m m m | 1999 | 0.9286/0.4656 | 0.1515 | 10/28 |
| 17 | 2002673 | Tricopper magnesium antimony oxide (3/1/1.4/6) | Cu3 Mg O6 Sb1.4 | R 3 c :H | 1995 | 0.8571/0.2187 | 0.1302 | 10/28 |
| 18 | 1526132 | Er (Fe Mn) O4 | Er Fe Mn O4 | R -3 m :H | 2000 | 0.8571/0.4839 | 0.1207 | 10/28 |
| 19 | 2234027 | Didysprosium heptanickel | Dy2 Ni7 | R -3 m :H | 2012 | 0.8214/0.3361 | 0.2702 | 9/28 |
| 20 | 1521200 | Tl1.810 (Ba0.75 Sr0.25)2 Ca2 Cu3 O9.784 | Ba1.5 Ca2 Cu3 O9.784 Sr0.5 Tl1.81 | I 4/m m m | 1999 | 0.9643/0.5242 | 0.2498 | 9/28 |
| 21 | 1521199 | Tl1.772 (Ba0.80 Sr0.20)2 Ca2 Cu3 O9.752 | Ba1.6 Ca2 Cu3 O9.752 Sr0.4 Tl1.772 | I 4/m m m | 1999 | 0.9643/0.5282 | 0.2349 | 9/28 |
| 22 | 1532574 | (Tl0.889 Cu0.067)2 Ba2 Ca2 Cu3 O9.808 | Ba2 Ca2 Cu3.134 O9.808 Tl1.778 | I 4/m m m | 2003 | 0.9643/0.4632 | 0.1583 | 9/28 |
| 23 | 1521198 | Tl1.778 (Ba0.85 Sr0.15)2 Ca2 Cu3 O9.728 | Ba1.7 Ca2 Cu3 O9.728 Sr0.3 Tl1.778 | I 4/m m m | 1999 | 0.9643/0.4713 | 0.2132 | 8/28 |
| 24 | 1561710 | COD 1561710 | F1.5 Na4.5 O6 S1.5 | R -3 m :H | 2015 | 0.9286/0.4729 | 0.1308 | 8/28 |
| 25 | 1573840 | Di(furan-2-yl)ethanedione | C10 H6 O4 | F d d 2 | 2025 | 0.9286/0.5309 | 0.1263 | 8/28 |
| 26 | 1524402 | In Cu Al O4 | Al Cu In O4 | R -3 m :H | 1980 | 0.75/0.3639 | 0.0824 | 7/28 |
| 27 | 1000039 | Calcium cyclo-hexaaluminate **[anchor]** | Al6 Ca9 O18 | P a -3 | 1975 | 0/0 | 0.1132 | 6/28 |
| 28 | 1538059 | Ti Se | Se Ti | P 63/m m c | 1959 | 0.5/0.3771 | 0.316 | 5/28 |
| 29 | 1538290 | V Se | Se V | P 63/m m c | 1939 | 0.5/0.4927 | 0.2804 | 5/28 |
| 30 | 1540838 | (Y3 Zr2)0.4 | Y1.2 Zr0.8 | P 63/m m c | 1972 | 0.5/0.5629 | 0.1318 | 5/28 |
| 31 | 2002599 | Tribarium tetradysprosium oxide | Ba3 Dy4 O9 | R 3 :H | 1993 | 0.8929/0.3753 | 0.0964 | 5/28 |
| 32 | 2300518 | COD 2300518 | Tl | P 63/m m c | 1994 | 0.5/0.5417 | 0.0867 | 4/28 |
| 33 | 1528102 | (Sn0.017 Tl0.983) | Sn0.017 Tl0.983 | P 63/m m c | 1960 | 0.5/0.5416 | 0.083 | 4/28 |
| 34 | 1000053 | Periclase **[anchor]** | Mg O | F m -3 m | 1979 | 0.2143/0.096 | 0.1019 | 3/28 |
| 35 | 1563205 | COD 1563205 **[anchor]** | Ca O | F m -3 m | 2011 | 0/0 | 0.0557 | 3/28 |
| 36 | 1524987 | (Cr0.75 Ru0.25) | Cr0.75 Ru0.25 | I m -3 m | 1986 | 0.1429/0.1023 | 0.0258 | 1/28 |
| 37 | 1525300 | (Cr0.333 Mo0.333 V0.334) | Cr0.333 Mo0.333 V0.334 | I m -3 m | 1976 | 0.25/0.1276 | 0.024 | 1/28 |
| 38 | 1524276 | COD 1524276 | Cr0.7 Mo0.3 | I m -3 m | 1982 | 0.1429/0.1023 | 0.0224 | 1/28 |
| 39 | 1522976 | COD 1522976 | Fe0.8 Ti0.2 | I m -3 m | 1986 | 0.1429/0.1021 | 0.0209 | 1/28 |
| 40 | 1522937 | (Fe0.7 Zn0.3) | Fe0.7 Zn0.3 | I m -3 m | 1964 | 0.1429/0.1022 | 0.0209 | 1/28 |
| 41 | 1521913 | Na0.288 Hf N Cl | Cl Hf N Na0.288 | R -3 m :H | 1999 | 0.9286/0.6167 | 0.0 | 0/28 |
| 42 | 1001366 | Thallium(I,III) barium calcium copper(III) oxide (1.82/2/1.9/3/10.94) | Ba2 Ca1.9 Cu3 O10.94 Tl1.82 | I 4/m m m | 1988 | 0.9643/0.4637 | 0.0 | 0/28 |

CIF records: [2002081](https://www.crystallography.net/cod/2002081.html)  , [2106554](https://www.crystallography.net/cod/2106554.html)  , [1524276](https://www.crystallography.net/cod/1524276.html)  , [1524987](https://www.crystallography.net/cod/1524987.html)  , [1522976](https://www.crystallography.net/cod/1522976.html)  , [1522937](https://www.crystallography.net/cod/1522937.html)  , [1525300](https://www.crystallography.net/cod/1525300.html)  , [1524402](https://www.crystallography.net/cod/1524402.html)  , [1540838](https://www.crystallography.net/cod/1540838.html)  , [1525036](https://www.crystallography.net/cod/1525036.html)  , [2234027](https://www.crystallography.net/cod/2234027.html)  , [1522366](https://www.crystallography.net/cod/1522366.html)  , [1520006](https://www.crystallography.net/cod/1520006.html)  , [1521199](https://www.crystallography.net/cod/1521199.html)  , [2002673](https://www.crystallography.net/cod/2002673.html)  , [1538290](https://www.crystallography.net/cod/1538290.html)  , [1521913](https://www.crystallography.net/cod/1521913.html)  , [1521200](https://www.crystallography.net/cod/1521200.html)  , [1526132](https://www.crystallography.net/cod/1526132.html)  , [2310872](https://www.crystallography.net/cod/2310872.html)  , [1528102](https://www.crystallography.net/cod/1528102.html)  , [2300518](https://www.crystallography.net/cod/2300518.html)  , [1518102](https://www.crystallography.net/cod/1518102.html)  , [1538059](https://www.crystallography.net/cod/1538059.html)  , [1539099](https://www.crystallography.net/cod/1539099.html)  , [1539095](https://www.crystallography.net/cod/1539095.html)  , [1538413](https://www.crystallography.net/cod/1538413.html)  , [1522301](https://www.crystallography.net/cod/1522301.html)  , [1525851](https://www.crystallography.net/cod/1525851.html)  , [1561710](https://www.crystallography.net/cod/1561710.html)  , [1530012](https://www.crystallography.net/cod/1530012.html)  , [1001366](https://www.crystallography.net/cod/1001366.html)  , [4104493](https://www.crystallography.net/cod/4104493.html)  , [1573840](https://www.crystallography.net/cod/1573840.html)  , [2002599](https://www.crystallography.net/cod/2002599.html)  , [1521198](https://www.crystallography.net/cod/1521198.html)  , [1532574](https://www.crystallography.net/cod/1532574.html)  , [1000053](https://www.crystallography.net/cod/1000053.html)  , [1200009](https://www.crystallography.net/cod/1200009.html)  , [2312428](https://www.crystallography.net/cod/2312428.html)  , [1000039](https://www.crystallography.net/cod/1000039.html)  , [1563205](https://www.crystallography.net/cod/1563205.html)  

## aluminate_enriched_residue_clinkerNIST_180718_R1.xrdml
- fingerprint: 62 peaks (λ=1.54060 Å); published phases: ferrite C4AF / periclase MgO / ortho+cub-aluminate C3A / aphthitalite

| # | COD | mineral | formula | SG | yr | fast cv/cr | simII | match |
|---|-----|---------|---------|----|----|------------|-------|-------|
| 1 | 1525851 | Te2 V2 O9 | O9 Te2 V2 | F d d 2 | 1973 | 0.871/0.339 | 0.2295 | 27/62 |
| 2 | 2312428 | COD 2312428 **[anchor]** | Ca8 O16 Si4 | P n m a | 2024 | 0/0 | 0.2254 | 27/62 |
| 3 | 1200009 | Brownmillerite **[anchor]** | Al Ca2 Fe O5 | I b m 2 | 1971 | 0/0 | 0.2846 | 26/62 |
| 4 | 1532574 | (Tl0.889 Cu0.067)2 Ba2 Ca2 Cu3 O9.808 | Ba2 Ca2 Cu3.134 O9.808 Tl1.778 | I 4/m m m | 2003 | 0.8065/0.3241 | 0.2126 | 26/62 |
| 5 | 1521199 | Tl1.772 (Ba0.80 Sr0.20)2 Ca2 Cu3 O9.752 | Ba1.6 Ca2 Cu3 O9.752 Sr0.4 Tl1.772 | I 4/m m m | 1999 | 0.7903/0.2979 | 0.1943 | 26/62 |
| 6 | 1521198 | Tl1.778 (Ba0.85 Sr0.15)2 Ca2 Cu3 O9.728 | Ba1.7 Ca2 Cu3 O9.728 Sr0.3 Tl1.778 | I 4/m m m | 1999 | 0.7903/0.3005 | 0.2068 | 25/62 |
| 7 | 1521200 | Tl1.810 (Ba0.75 Sr0.25)2 Ca2 Cu3 O9.784 | Ba1.5 Ca2 Cu3 O9.784 Sr0.5 Tl1.81 | I 4/m m m | 1999 | 0.7903/0.2911 | 0.1933 | 25/62 |
| 8 | 4104493 | catena-(dimethylammonium)-(tris(μ~2~-Formato-O,O')-iron(ii) | C5 H15 Fe N O6 | R -3 c :H | 2009 | 0.8226/0.3286 | 0.223 | 24/62 |
| 9 | 1530012 | Te2 O (O H)2 F4 | F4 H2 O3 Te2 | F 2 d d | 1982 | 0.8871/0.4799 | 0.2178 | 24/62 |
| 10 | 1539099 | Cd0.5 Bi4 O4 Cl5 | Bi4 Cd0.5 Cl5 O4 | I 4/m m m | 1941 | 0.7903/0.2577 | 0.1779 | 23/62 |
| 11 | 1520006 | COD 1520006 | Bi9 O7.5 S6 | R -3 m :H | 2015 | 0.6935/0.308 | 0.2072 | 22/62 |
| 12 | 2310872 | Ca3 (Si O4) O | Ca3 O5 Si | R 3 m :H | 1952 | 0.8226/0.263 | 0.2158 | 21/62 |
| 13 | 1522366 | COD 1522366 | Gd2 Ni7 | R -3 m :H | 1969 | 0.7097/0.5249 | 0.1967 | 21/62 |
| 14 | 1539095 | Ca.80 Bi3.80 O4 Cl5 | Bi3.8 Ca0.8 Cl5 O4 | I 4/m m m | 1941 | 0.8548/0.2548 | 0.1897 | 21/62 |
| 15 | 1525036 | COD 1525036 | Co7 Er2 | R -3 m :H | 1967 | 0.7419/0.5367 | 0.2169 | 20/62 |
| 16 | 2002673 | Tricopper magnesium antimony oxide (3/1/1.4/6) | Cu3 Mg O6 Sb1.4 | R 3 c :H | 1995 | 0.7903/0.2327 | 0.1873 | 19/62 |
| 17 | 1538413 | Ca3 (Si O4) O | Ca3 O5 Si | R 3 m :H | 1950 | 0.7581/0.3177 | 0.201 | 18/62 |
| 18 | 2234027 | Didysprosium heptanickel | Dy2 Ni7 | R -3 m :H | 2012 | 0.6935/0.5211 | 0.2134 | 17/62 |
| 19 | 2106554 | (Yb0.5 Eu0.5) Fe2 O4 | Eu0.5 Fe2 O4 Yb0.5 | R -3 m :H | 1975 | 0.6452/0.3437 | 0.1908 | 17/62 |
| 20 | 1524402 | In Cu Al O4 | Al Cu In O4 | R -3 m :H | 1980 | 0.6129/0.3553 | 0.1708 | 17/62 |
| 21 | 1518102 | COD 1518102 | Na2 O11 Ta4 | R -3 c :H | 1985 | 0.9032/0.2873 | 0.1915 | 16/62 |
| 22 | 1526132 | Er (Fe Mn) O4 | Er Fe Mn O4 | R -3 m :H | 2000 | 0.6613/0.2586 | 0.1657 | 15/62 |
| 23 | 1522301 | Bi2 Sr2 Nb2 Mn O11.57 | Bi2 Mn Nb2 O11.57 Sr2 | F m m m | 1999 | 0.7581/0.2082 | 0.15 | 15/62 |
| 24 | 2002081 | COD 2002081 | Cr0.5 K0.5 S2 Sn0.5 | R -3 m :H | 1993 | 0.5806/0.3517 | 0.1895 | 14/62 |
| 25 | 2002599 | Tribarium tetradysprosium oxide | Ba3 Dy4 O9 | R 3 :H | 1993 | 0.8226/0.3576 | 0.1369 | 13/62 |
| 26 | 1000039 | Calcium cyclo-hexaaluminate **[anchor]** | Al6 Ca9 O18 | P a -3 | 1975 | 0/0 | 0.1986 | 12/62 |
| 27 | 1561710 | COD 1561710 | F1.5 Na4.5 O6 S1.5 | R -3 m :H | 2015 | 0.8226/0.3109 | 0.1817 | 12/62 |
| 28 | 1573840 | Di(furan-2-yl)ethanedione | C10 H6 O4 | F d d 2 | 2025 | 0.9032/0.3271 | 0.1608 | 11/62 |
| 29 | 1528102 | (Sn0.017 Tl0.983) | Sn0.017 Tl0.983 | P 63/m m c | 1960 | 0.371/0.3196 | 0.1946 | 9/62 |
| 30 | 2300518 | COD 2300518 | Tl | P 63/m m c | 1994 | 0.371/0.3196 | 0.1854 | 9/62 |
| 31 | 1540838 | (Y3 Zr2)0.4 | Y1.2 Zr0.8 | P 63/m m c | 1972 | 0.3548/0.3184 | 0.1541 | 8/62 |
| 32 | 1538290 | V Se | Se V | P 63/m m c | 1939 | 0.3065/0.2339 | 0.1077 | 8/62 |
| 33 | 1538059 | Ti Se | Se Ti | P 63/m m c | 1959 | 0.3065/0.1801 | 0.1022 | 6/62 |
| 34 | 1563205 | COD 1563205 **[anchor]** | Ca O | F m -3 m | 2011 | 0/0 | 0.0902 | 5/62 |
| 35 | 1000053 | Periclase **[anchor]** | Mg O | F m -3 m | 1979 | 0.2419/0.4892 | 0.1269 | 3/62 |
| 36 | 1525300 | (Cr0.333 Mo0.333 V0.334) | Cr0.333 Mo0.333 V0.334 | I m -3 m | 1976 | 0.1613/0.6103 | 0.0531 | 3/62 |
| 37 | 1524987 | (Cr0.75 Ru0.25) | Cr0.75 Ru0.25 | I m -3 m | 1986 | 0.1452/0.5043 | 0.0591 | 2/62 |
| 38 | 1522937 | (Fe0.7 Zn0.3) | Fe0.7 Zn0.3 | I m -3 m | 1964 | 0.1452/0.5041 | 0.0522 | 2/62 |
| 39 | 1522976 | COD 1522976 | Fe0.8 Ti0.2 | I m -3 m | 1986 | 0.1452/0.5038 | 0.0521 | 2/62 |
| 40 | 1524276 | COD 1524276 | Cr0.7 Mo0.3 | I m -3 m | 1982 | 0.1452/0.5043 | 0.0513 | 2/62 |
| 41 | 1521913 | Na0.288 Hf N Cl | Cl Hf N Na0.288 | R -3 m :H | 1999 | 0.6774/0.2738 | 0.0 | 0/62 |
| 42 | 1001366 | Thallium(I,III) barium calcium copper(III) oxide (1.82/2/1.9/3/10.94) | Ba2 Ca1.9 Cu3 O10.94 Tl1.82 | I 4/m m m | 1988 | 0.8065/0.3159 | 0.0 | 0/62 |

CIF records: [2002081](https://www.crystallography.net/cod/2002081.html)  , [2106554](https://www.crystallography.net/cod/2106554.html)  , [1524276](https://www.crystallography.net/cod/1524276.html)  , [1524987](https://www.crystallography.net/cod/1524987.html)  , [1522976](https://www.crystallography.net/cod/1522976.html)  , [1522937](https://www.crystallography.net/cod/1522937.html)  , [1525300](https://www.crystallography.net/cod/1525300.html)  , [1524402](https://www.crystallography.net/cod/1524402.html)  , [1540838](https://www.crystallography.net/cod/1540838.html)  , [1525036](https://www.crystallography.net/cod/1525036.html)  , [2234027](https://www.crystallography.net/cod/2234027.html)  , [1522366](https://www.crystallography.net/cod/1522366.html)  , [1520006](https://www.crystallography.net/cod/1520006.html)  , [1521199](https://www.crystallography.net/cod/1521199.html)  , [2002673](https://www.crystallography.net/cod/2002673.html)  , [1538290](https://www.crystallography.net/cod/1538290.html)  , [1521913](https://www.crystallography.net/cod/1521913.html)  , [1521200](https://www.crystallography.net/cod/1521200.html)  , [1526132](https://www.crystallography.net/cod/1526132.html)  , [2310872](https://www.crystallography.net/cod/2310872.html)  , [1528102](https://www.crystallography.net/cod/1528102.html)  , [2300518](https://www.crystallography.net/cod/2300518.html)  , [1518102](https://www.crystallography.net/cod/1518102.html)  , [1538059](https://www.crystallography.net/cod/1538059.html)  , [1539099](https://www.crystallography.net/cod/1539099.html)  , [1539095](https://www.crystallography.net/cod/1539095.html)  , [1538413](https://www.crystallography.net/cod/1538413.html)  , [1522301](https://www.crystallography.net/cod/1522301.html)  , [1525851](https://www.crystallography.net/cod/1525851.html)  , [1561710](https://www.crystallography.net/cod/1561710.html)  , [1530012](https://www.crystallography.net/cod/1530012.html)  , [1001366](https://www.crystallography.net/cod/1001366.html)  , [4104493](https://www.crystallography.net/cod/4104493.html)  , [1573840](https://www.crystallography.net/cod/1573840.html)  , [2002599](https://www.crystallography.net/cod/2002599.html)  , [1521198](https://www.crystallography.net/cod/1521198.html)  , [1532574](https://www.crystallography.net/cod/1532574.html)  , [1000053](https://www.crystallography.net/cod/1000053.html)  , [1200009](https://www.crystallography.net/cod/1200009.html)  , [2312428](https://www.crystallography.net/cod/2312428.html)  , [1000039](https://www.crystallography.net/cod/1000039.html)  , [1563205](https://www.crystallography.net/cod/1563205.html)  

## Clinker_Synchrotron.dat
- fingerprint: 31 peaks (λ=0.82543 Å); published phases: alite C3S / beta-belite C2S / ferrite C4AF / periclase MgO / aluminate

| # | COD | mineral | formula | SG | yr | fast cv/cr | simII | match |
|---|-----|---------|---------|----|----|------------|-------|-------|
| 1 | 1538413 | Ca3 (Si O4) O | Ca3 O5 Si | R 3 m :H | 1950 | 0.9677/0.5375 | 0.4531 | 22/31 |
| 2 | 2310872 | Ca3 (Si O4) O | Ca3 O5 Si | R 3 m :H | 1952 | 0.9032/0.4763 | 0.2302 | 19/31 |
| 3 | 1530012 | Te2 O (O H)2 F4 | F4 H2 O3 Te2 | F 2 d d | 1982 | 0.9032/0.532 | 0.2279 | 18/31 |
| 4 | 1525851 | Te2 V2 O9 | O9 Te2 V2 | F d d 2 | 1973 | 0.9355/0.3583 | 0.2094 | 17/31 |
| 5 | 1520006 | COD 1520006 | Bi9 O7.5 S6 | R -3 m :H | 2015 | 0.871/0.6366 | 0.2098 | 15/31 |
| 6 | 4104493 | catena-(dimethylammonium)-(tris(μ~2~-Formato-O,O')-iron(ii) | C5 H15 Fe N O6 | R -3 c :H | 2009 | 0.871/0.141 | 0.2065 | 14/31 |
| 7 | 2312428 | COD 2312428 **[anchor]** | Ca8 O16 Si4 | P n m a | 2024 | 0/0 | 0.2637 | 13/31 |
| 8 | 2002081 | COD 2002081 | Cr0.5 K0.5 S2 Sn0.5 | R -3 m :H | 1993 | 0.7742/0.4266 | 0.2332 | 13/31 |
| 9 | 1539099 | Cd0.5 Bi4 O4 Cl5 | Bi4 Cd0.5 Cl5 O4 | I 4/m m m | 1941 | 0.9032/0.6038 | 0.228 | 13/31 |
| 10 | 1518102 | COD 1518102 | Na2 O11 Ta4 | R -3 c :H | 1985 | 0.871/0.6629 | 0.209 | 13/31 |
| 11 | 1522301 | Bi2 Sr2 Nb2 Mn O11.57 | Bi2 Mn Nb2 O11.57 Sr2 | F m m m | 1999 | 0.9355/0.5108 | 0.1531 | 13/31 |
| 12 | 2106554 | (Yb0.5 Eu0.5) Fe2 O4 | Eu0.5 Fe2 O4 Yb0.5 | R -3 m :H | 1975 | 0.8387/0.5634 | 0.1524 | 13/31 |
| 13 | 1525036 | COD 1525036 | Co7 Er2 | R -3 m :H | 1967 | 0.9032/0.3437 | 0.2233 | 12/31 |
| 14 | 1522366 | COD 1522366 | Gd2 Ni7 | R -3 m :H | 1969 | 0.871/0.2738 | 0.1677 | 12/31 |
| 15 | 1526132 | Er (Fe Mn) O4 | Er Fe Mn O4 | R -3 m :H | 2000 | 0.8065/0.4594 | 0.1511 | 12/31 |
| 16 | 1539095 | Ca.80 Bi3.80 O4 Cl5 | Bi3.8 Ca0.8 Cl5 O4 | I 4/m m m | 1941 | 0.871/0.54 | 0.1969 | 11/31 |
| 17 | 1200009 | Brownmillerite **[anchor]** | Al Ca2 Fe O5 | I b m 2 | 1971 | 0/0 | 0.2347 | 10/31 |
| 18 | 2234027 | Didysprosium heptanickel | Dy2 Ni7 | R -3 m :H | 2012 | 0.8387/0.3049 | 0.2213 | 10/31 |
| 19 | 2002673 | Tricopper magnesium antimony oxide (3/1/1.4/6) | Cu3 Mg O6 Sb1.4 | R 3 c :H | 1995 | 0.8065/0.2026 | 0.1799 | 10/31 |
| 20 | 1532574 | (Tl0.889 Cu0.067)2 Ba2 Ca2 Cu3 O9.808 | Ba2 Ca2 Cu3.134 O9.808 Tl1.778 | I 4/m m m | 2003 | 0.9032/0.5314 | 0.175 | 9/31 |
| 21 | 1521199 | Tl1.772 (Ba0.80 Sr0.20)2 Ca2 Cu3 O9.752 | Ba1.6 Ca2 Cu3 O9.752 Sr0.4 Tl1.772 | I 4/m m m | 1999 | 0.9355/0.6153 | 0.2062 | 8/31 |
| 22 | 1521198 | Tl1.778 (Ba0.85 Sr0.15)2 Ca2 Cu3 O9.728 | Ba1.7 Ca2 Cu3 O9.728 Sr0.3 Tl1.778 | I 4/m m m | 1999 | 0.9032/0.5613 | 0.1996 | 8/31 |
| 23 | 1524402 | In Cu Al O4 | Al Cu In O4 | R -3 m :H | 1980 | 0.6774/0.3135 | 0.124 | 8/31 |
| 24 | 1521200 | Tl1.810 (Ba0.75 Sr0.25)2 Ca2 Cu3 O9.784 | Ba1.5 Ca2 Cu3 O9.784 Sr0.5 Tl1.81 | I 4/m m m | 1999 | 0.9355/0.6022 | 0.2155 | 7/31 |
| 25 | 1538290 | V Se | Se V | P 63/m m c | 1939 | 0.5161/0.4759 | 0.1848 | 7/31 |
| 26 | 1561710 | COD 1561710 | F1.5 Na4.5 O6 S1.5 | R -3 m :H | 2015 | 0.871/0.4353 | 0.1742 | 7/31 |
| 27 | 1573840 | Di(furan-2-yl)ethanedione | C10 H6 O4 | F d d 2 | 2025 | 0.9355/0.5294 | 0.1699 | 7/31 |
| 28 | 2002599 | Tribarium tetradysprosium oxide | Ba3 Dy4 O9 | R 3 :H | 1993 | 0.8387/0.3947 | 0.1165 | 7/31 |
| 29 | 1538059 | Ti Se | Se Ti | P 63/m m c | 1959 | 0.5161/0.3583 | 0.2042 | 6/31 |
| 30 | 1000039 | Calcium cyclo-hexaaluminate **[anchor]** | Al6 Ca9 O18 | P a -3 | 1975 | 0/0 | 0.1475 | 6/31 |
| 31 | 1540838 | (Y3 Zr2)0.4 | Y1.2 Zr0.8 | P 63/m m c | 1972 | 0.5806/0.4825 | 0.1267 | 6/31 |
| 32 | 2300518 | COD 2300518 | Tl | P 63/m m c | 1994 | 0.5161/0.4496 | 0.1158 | 5/31 |
| 33 | 1528102 | (Sn0.017 Tl0.983) | Sn0.017 Tl0.983 | P 63/m m c | 1960 | 0.5161/0.4494 | 0.1153 | 5/31 |
| 34 | 1000053 | Periclase **[anchor]** | Mg O | F m -3 m | 1979 | 0.2581/0.0715 | 0.0732 | 4/31 |
| 35 | 1563205 | COD 1563205 **[anchor]** | Ca O | F m -3 m | 2011 | 0/0 | 0.0623 | 4/31 |
| 36 | 1524987 | (Cr0.75 Ru0.25) | Cr0.75 Ru0.25 | I m -3 m | 1986 | 0.1935/0.0769 | 0.0328 | 2/31 |
| 37 | 1522976 | COD 1522976 | Fe0.8 Ti0.2 | I m -3 m | 1986 | 0.1935/0.0767 | 0.0294 | 2/31 |
| 38 | 1522937 | (Fe0.7 Zn0.3) | Fe0.7 Zn0.3 | I m -3 m | 1964 | 0.1935/0.0768 | 0.0294 | 2/31 |
| 39 | 1524276 | COD 1524276 | Cr0.7 Mo0.3 | I m -3 m | 1982 | 0.1935/0.0769 | 0.0284 | 2/31 |
| 40 | 1525300 | (Cr0.333 Mo0.333 V0.334) | Cr0.333 Mo0.333 V0.334 | I m -3 m | 1976 | 0.3226/0.1053 | 0.0297 | 1/31 |
| 41 | 1521913 | Na0.288 Hf N Cl | Cl Hf N Na0.288 | R -3 m :H | 1999 | 0.8387/0.6173 | 0.0 | 0/31 |
| 42 | 1001366 | Thallium(I,III) barium calcium copper(III) oxide (1.82/2/1.9/3/10.94) | Ba2 Ca1.9 Cu3 O10.94 Tl1.82 | I 4/m m m | 1988 | 0.9032/0.5534 | 0.0 | 0/31 |

CIF records: [2002081](https://www.crystallography.net/cod/2002081.html)  , [2106554](https://www.crystallography.net/cod/2106554.html)  , [1524276](https://www.crystallography.net/cod/1524276.html)  , [1524987](https://www.crystallography.net/cod/1524987.html)  , [1522976](https://www.crystallography.net/cod/1522976.html)  , [1522937](https://www.crystallography.net/cod/1522937.html)  , [1525300](https://www.crystallography.net/cod/1525300.html)  , [1524402](https://www.crystallography.net/cod/1524402.html)  , [1540838](https://www.crystallography.net/cod/1540838.html)  , [1525036](https://www.crystallography.net/cod/1525036.html)  , [2234027](https://www.crystallography.net/cod/2234027.html)  , [1522366](https://www.crystallography.net/cod/1522366.html)  , [1520006](https://www.crystallography.net/cod/1520006.html)  , [1521199](https://www.crystallography.net/cod/1521199.html)  , [2002673](https://www.crystallography.net/cod/2002673.html)  , [1538290](https://www.crystallography.net/cod/1538290.html)  , [1521913](https://www.crystallography.net/cod/1521913.html)  , [1521200](https://www.crystallography.net/cod/1521200.html)  , [1526132](https://www.crystallography.net/cod/1526132.html)  , [2310872](https://www.crystallography.net/cod/2310872.html)  , [1528102](https://www.crystallography.net/cod/1528102.html)  , [2300518](https://www.crystallography.net/cod/2300518.html)  , [1518102](https://www.crystallography.net/cod/1518102.html)  , [1538059](https://www.crystallography.net/cod/1538059.html)  , [1539099](https://www.crystallography.net/cod/1539099.html)  , [1539095](https://www.crystallography.net/cod/1539095.html)  , [1538413](https://www.crystallography.net/cod/1538413.html)  , [1522301](https://www.crystallography.net/cod/1522301.html)  , [1525851](https://www.crystallography.net/cod/1525851.html)  , [1561710](https://www.crystallography.net/cod/1561710.html)  , [1530012](https://www.crystallography.net/cod/1530012.html)  , [1001366](https://www.crystallography.net/cod/1001366.html)  , [4104493](https://www.crystallography.net/cod/4104493.html)  , [1573840](https://www.crystallography.net/cod/1573840.html)  , [2002599](https://www.crystallography.net/cod/2002599.html)  , [1521198](https://www.crystallography.net/cod/1521198.html)  , [1532574](https://www.crystallography.net/cod/1532574.html)  , [1000053](https://www.crystallography.net/cod/1000053.html)  , [1200009](https://www.crystallography.net/cod/1200009.html)  , [2312428](https://www.crystallography.net/cod/2312428.html)  , [1000039](https://www.crystallography.net/cod/1000039.html)  , [1563205](https://www.crystallography.net/cod/1563205.html)  

## Periclase (MgO) deep-dive — why only 3-4 matched peaks?
- Diagnosis (`benchmarks/protocols/periclase_diag.py` → `data/unit12/work/periclase_diag.json`): the strong MgO lines 200 (d=a/2) and 220 (d=a/√8) land within 0.003 Å of a picked measured peak in **all four** samples (exact offsets +0.0004..−0.0013 Å). No systematic shift, no zero-error, no wavelength mis-set.
- Every out-of-window miss is a chemically weak rocksalt difference-reflection (222, 311) or a line buried under a stronger C3S/C2S/ferrite reflection (±0.04..0.17 Å) — expected even for a perfectly correct minor-phase match.
- **Conclusion:** periclase's 3-4 matched peaks / sim 0.07-0.13 is the *honest* score for a minor phase whose strong lines are already hit to <0.003 Å; no tolerance or peak-position fix is warranted. Caveat: the 111 (d=a/√3) window hits at +0.004..+0.017 Å are borderline and probably owned by neighbouring alite/ferrite lines, so the counting periclase's own fingerprint contribution as the 200+220 pair (~2 peaks/sample) is the safer reading.

## Honest limitations
- Index is position-only; intensities come from the GSAS-II CIF simulation step only.
- Index d-range is [1.1, 22] Å with an hkl cap of 28; huge-cell organics report a per-entry effective dmin.
- The COD metadata export itself excludes retracted / duplicate / erroneous entries.
- Clinker phases often appear as several COD records (mineral names): collapsed to the best record per name here.