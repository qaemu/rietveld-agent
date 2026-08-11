"""Instrument-aware synthetic noise model (Spike 05 threshold evaluation).

Every perturbation maps a clean PowderPattern to a realistic noisy one and is
deterministic given a seed (so the eval and its regression tests reproduce):

  * Counting statistics: Poisson sampling at a target peak-max count level.
  * Background drift: smooth quadratic + low-frequency cosine bump, as a
    fraction of the clean peak maximum (lab background wander).
  * Sample displacement: physically-motivated, angle-dependent 2theta shift
    for a Bragg-Brentano goniometer,  d(2th) = -2*s*cos(theta)/R  [rad].

The perturbed pattern still flows through the full ingest path
(find_peaks -> sample_fingerprint -> rank -> verdict), which is exactly the
point: we measure how decision metrics degrade, not how the model behaves.
"""
from __future__ import annotations

import numpy as np

from core.ingest.models import PowderPattern

#: typical Bragg-Brentano goniometer radius [mm]
R_GONIO_MM = 200.0


def add_poisson(y: np.ndarray, peak_max_counts: float, rng: np.random.Generator) -> np.ndarray:
    """Scale so the peak maximum equals ``peak_max_counts`` and sample Poisson."""
    scale = peak_max_counts / max(float(np.nanmax(y)), 1e-9)
    return rng.poisson(y * scale).astype(float)


def add_background_drift(y: np.ndarray, fraction: float,
                         rng: np.random.Generator) -> np.ndarray:
    """Smooth background wander; ``fraction`` of the clean peak maximum."""
    if fraction <= 0.0:
        return y
    n = y.size
    x = np.linspace(0.0, 1.0, n)
    r1 = rng.uniform(0.3, 1.0)
    r2 = rng.uniform(0.0, 2.0 * np.pi)
    drift = (r1 * x ** 2 + 0.5 * np.cos(3.0 * x + r2)) * fraction * float(np.nanmax(y))
    return y + drift


def sample_displacement(tth: np.ndarray, y: np.ndarray, s_mm: float) -> np.ndarray:
    """Apply a uniform sample displacement s [mm].

    In Bragg-Brentano geometry d(2theta) = -2*s*cos(theta)/R [radians], so the
    shift is largest at low angle. Edges keep their boundary values.
    """
    if s_mm == 0.0:
        return y
    delta_deg = np.degrees(-2.0 * s_mm * np.cos(np.radians(tth) / 2.0) / R_GONIO_MM)
    return np.interp(tth + delta_deg, tth, y, left=float(y[0]), right=float(y[-1]))


def perturb(pattern: PowderPattern, *, counts: float | None, s_mm: float,
            bg_fraction: float, seed: int) -> PowderPattern:
    """Perturbed copy of ``pattern`` (geometry first, counting noise last)."""
    rng = np.random.default_rng(seed)
    y = np.asarray(pattern.intensity, dtype=float)
    y = sample_displacement(pattern.tth, y, s_mm)
    y = add_background_drift(y, bg_fraction, rng)
    y = np.clip(y, 0.0, None)   # counts cannot be negative
    if counts is not None and counts > 0:
        y = add_poisson(y, counts, rng)
    return PowderPattern(sample_name=pattern.sample_name, source=pattern.source,
                         tth=pattern.tth, intensity=y, instrument=pattern.instrument,
                         metadata=dict(pattern.metadata, perturbed={"seed": seed,
                                                                    "counts": counts,
                                                                    "s_mm": s_mm,
                                                                    "bg_fraction": bg_fraction}))


def amorphous_pattern(tth: np.ndarray, *, counts: float | None, s_mm: float,
                      bg_fraction: float, seed: int,
                      instrument) -> PowderPattern:
    """Synthetic featureless (amorphous-like) pattern: broad humps only.

    Negative control: must never produce a supported verdict, at any noise
    level, because no library candidate genuinely explains it.
    """
    rng = np.random.default_rng(seed + 10 ** 6)
    y = np.zeros_like(tth)
    for center, width, rel in ((22.0, 7.0, 1.0), (40.0, 10.0, 0.6),
                               (65.0, 14.0, 0.35)):
        y += rel * np.exp(-0.5 * ((tth - center) / width) ** 2)
    y *= 1e6
    y = sample_displacement(tth, y, s_mm)
    y = add_background_drift(y, bg_fraction, rng)
    y = np.clip(y, 0.0, None)
    if counts is not None and counts > 0:
        y = add_poisson(y, counts, rng)
    return PowderPattern(sample_name="amorphous-control", source="synthetic",
                         tth=tth, intensity=y, instrument=instrument)