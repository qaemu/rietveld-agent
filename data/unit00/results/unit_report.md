# Unit 01: GSAS-II scriptable validation

- Python 3.13.7 on Darwin x86_64
- GSAS-II v5.6.3 (vendor commit 013259a0defab16f8414e1fcfa6b0274f1c61ecf)
- macOS-15.6.1-x86_64-i386-64bit-Mach-O

## Tutorial reproduction (PbSO4, Cu Ka):

| stage | label | wR |
|---|---|---|
| background | background 3x chebyshev-1 | 40.8842525760248 |
| shift_scale | sample shift + HAP scale | 26.584410372762573 |
| cell | cell | 25.727293228673776 |
| instrument | instrument U,V,W,X,Y | 13.605982005811496 |
| mustrain_size | HAP mustrain(iso) + size(iso) flags | 13.607504527069162 |
| lebail | LeBail on -> refine -> off -> refine | 13.748274413305401 |
| phase_fraction | HAP PhaseFraction flag on/off | 13.765836678026371 |
| limits | Limits trim 15..100 -> restore 15..140 | 14.247419236836462 |
| atoms | Atoms xyz+Uiso (all) | 11.996066279972306 |

Final wR = 11.9961, cell = {'length_a': 8.48008, 'length_b': 5.39852, 'length_c': 6.95978, 'angle_alpha': 90.0, 'angle_beta': 90.0, 'angle_gamma': 90.0, 'volume': 318.61788}

## Simulator
- synthetic pattern: 6251 pts from 2th=15.0..140.0, peak at 2th=29.68 I=47462238, wR=0.11723020247314325

## Verdict
- [ ] tutorial Rwp trajectory matches published tutorial ballpark
- [ ] all allowlist refinement keys applied successfully
- [ ] simulator produced a physical pattern and exported XYE
- [ ] per-stage .gpx checkpoints saved
