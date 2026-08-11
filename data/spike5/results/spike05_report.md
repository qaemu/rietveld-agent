# Spike 05: threshold-robustness evaluation (M1 fingerprint thresholds)

Constants under test: `min_top_similarity = 0.35`, `min_margin = 0.10` with the **family-aware margin** (top similarity minus the best different-phase-family candidate; introduced by spike 05: statuses are per phase family, so differently-realized entries of the SAME family must not compete). Library manifest `340127a8343ea8f0...`, 50 seeded replicates per level (deterministic RNG).

## Noise envelope (cumulative severity)

| level | peak counts | displacement (mm) | background drift |
|---|---|---|---|
| L0 | None | 0.0 | 0.0 |
| L1 | 1000000 | 0.02 | 0.005 |
| L2 | 300000 | 0.05 | 0.01 |
| L3 | 100000 | 0.1 | 0.02 |
| L4 | 30000 | 0.2 | 0.05 |
| L5 | 10000 | 0.3 | 0.1 |
| L6 | 3000 | 0.5 | 0.15 |
| L7 | 1000 | 0.8 | 0.2 |

## Pass 1: evidence-only (no statistics gate)

Abstain rate; median sim / median family margin.

| material | L0 | L1 | L2 | L3 | L4 | L5 | L6 | L7 | evidence flip |
|---|---|---|---|---|---|---|---|---|---|
| pbso4-cu | 0.00 (1.000/0.914) | 0.00 (0.982/0.900) | 0.00 (0.903/0.825) | 0.00 (0.721/0.613) | 0.00 (0.570/0.426) | 0.48 (0.355/0.275) | 1.00 (0.130/0.016) | 0.00 (0.421/0.370) | L6 |
| sio2-cu | 0.00 (1.000/0.926) | 0.00 (0.977/0.915) | 0.00 (0.880/0.834) | 0.00 (0.671/0.640) | 1.00 (0.262/0.228) | 1.00 (0.065/0.037) | 1.00 (0.119/0.103) | 1.00 (0.074/0.029) | L4 |
| pbso4-fe | 0.00 (1.000/0.920) | 0.00 (0.969/0.892) | 0.00 (0.892/0.821) | 0.00 (0.791/0.730) | 0.00 (0.516/0.396) | 0.00 (0.476/0.354) | 1.00 (0.150/0.100) | 1.00 (0.098/0.025) | L6 |
| nacl-cu | 0.00 (1.000/0.900) | 0.00 (0.985/0.902) | 0.00 (0.893/0.832) | 0.00 (0.652/0.614) | 1.00 (0.332/0.298) | 1.00 (0.119/0.076) | 1.00 (0.047/0.017) | 1.00 (0.079/0.070) | L4 |

## Statistics gate

- Derived rule: counts at the last severity level at which EVERY material is fully supported -> **min_peak_counts = 100000** (level L3).
- Implemented in `core/verdict` (decide(..., min_peak_counts=100000)).

## Pass 2: fused pipeline (evidence + statistics gate, as shipped)

| material | L0 | L1 | L2 | L3 | L4 | L5 | L6 | L7 | fused flip |
|---|---|---|---|---|---|---|---|---|---|
| pbso4-cu | 0.00 | 0.00 | 0.00 | 0.52 | 1.00 | 1.00 | 1.00 | 1.00 | L3 |
| sio2-cu | 0.00 | 0.00 | 0.00 | 0.58 | 1.00 | 1.00 | 1.00 | 1.00 | L3 |
| pbso4-fe | 0.00 | 0.00 | 0.00 | 0.42 | 1.00 | 1.00 | 1.00 | 1.00 | L4 |
| nacl-cu | 0.00 | 0.00 | 0.00 | 0.58 | 1.00 | 1.00 | 1.00 | 1.00 | L3 |

## Controls
- Amorphous negative control (fused): max false-positive rate **0.0000** (must be 0).
- L0 clean must reproduce M1 e2e (all supported): yes.
- Abstain-rate monotonicity (fused): ok.

## Verdict
- **KEEP thresholds 0.35/0.10 with family-aware margin: evidence remains majority-supported through L3 (first majority flip at L4; critical: sio2-cu / nacl-cu). Statistics gate min_peak_counts = 100000 (derived at envelope level L3); the fused pipeline intentionally abstains from the gate boundary (L3) onward.**
- Fused flips: {"pbso4-cu": "L3", "sio2-cu": "L3", "pbso4-fe": "L4", "nacl-cu": "L3"}
- Wall clock: 144.7s
- **passed = True**
