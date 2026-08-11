"""Instrument and sample fingerprints for powder patterns.

Design goals
------------
* Instrument fingerprint: a stable, hashable description of the measuring
  setup (anode, wavelengths, goniometer geometry, scan grid) used to map a
  data file onto the correct instrument parameter (PRM) in the calibration
  registry, and to group "same instrument" observations.
* Sample fingerprint: a material signature in *d-spacing space* so that the
  same sample measured on different anodes (Cu vs Fe ...) still matches;
  supplied as both a human-readable peak line list and a grid profile, with
  a similarity metric for the "same sample as previous batch?" hypothesis.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np

from .models import InstrumentParams, PowderPattern

#: d-spacing tolerance [Angstrom] for peak matching
PEAK_D_TOL = 0.02
#: d-spacing grid for profile similarity (Angstrom)
D_GRID_MIN, D_GRID_STEP = 0.5, 0.005
#: minimum prominence (fraction of max intensity) for a peak to be reported
PEAK_PROMINENCE = 0.02


def tth_to_d(tth_deg: np.ndarray | float, wavelength: float) -> np.ndarray | float:
    """Bragg 2theta [deg] -> d-spacing [Angstrom] for a given wavelength."""
    theta = np.radians(tth_deg) / 2.0
    return wavelength / (2.0 * np.sin(theta))


@dataclass
class Peak:
    tth: float          # degrees
    d: float            # Angstrom (from kalpha1)
    height: float       # intensity above local baseline
    fwhm: float         # approximate full width at half max, degrees

    def as_dict(self) -> dict:
        return {"tth": round(float(self.tth), 4), "d": round(float(self.d), 4),
                "height": float(self.height), "fwhm": round(float(self.fwhm), 4)}


@dataclass
class InstrumentFingerprint:
    params: InstrumentParams
    id: str

    @classmethod
    def from_pattern(cls, pattern: PowderPattern) -> "InstrumentFingerprint":
        return cls(params=pattern.instrument, id=pattern.instrument.fingerprint_id)


@dataclass
class SampleFingerprint:
    source: str
    peaks: List[Peak] = field(default_factory=list)
    #: normalized intensity profile on a uniform d grid (d -> intensity)
    d_grid: np.ndarray = field(default_factory=lambda: np.array([]))
    d_profile: np.ndarray = field(default_factory=lambda: np.array([]))

    def peak_line_list(self, max_n: int = 12) -> List[dict]:
        return [p.as_dict() for p in self.peaks[:max_n]]

    def to_dict(self) -> dict:
        """Serialize for a candidate library / registry (plain JSON types)."""
        return {
            "source": self.source,
            "peaks": [p.as_dict() for p in self.peaks],
            "d_grid": [float(v) for v in self.d_grid],
            "d_profile": [float(v) for v in self.d_profile],
        }

    @classmethod
    def from_dict(cls, data: dict) -> "SampleFingerprint":
        peaks = [Peak(tth=p["tth"], d=p["d"], height=p["height"], fwhm=p["fwhm"])
                 for p in data.get("peaks", [])]
        return cls(
            source=data.get("source", ""),
            peaks=peaks,
            d_grid=np.asarray(data.get("d_grid", []), dtype=float),
            d_profile=np.asarray(data.get("d_profile", []), dtype=float),
        )


def find_peaks(tth: np.ndarray, y: np.ndarray, prominence: float = PEAK_PROMINENCE,
               window: int = 5) -> List[Peak]:
    """Simple local-maximum peak finder with parabolic refinement.

    A candidate is kept when it is the strict maximum inside a window of
    ``window`` bins on each side and its prominence (rise above the lower of
    the left/right window minima) exceeds ``prominence * y.max()``.
    """
    ymax = float(np.nanmax(y))
    if ymax <= 0:
        return []
    thr = prominence * ymax
    peaks: List[Peak] = []
    t = np.asarray(tth)
    for i in range(window, len(y) - window):
        seg = y[i - window: i + window + 1]
        if y[i] != seg.max():
            continue
        left_min = float(seg[:window].min())
        right_min = float(seg[window + 1:].min())
        base = min(left_min, right_min)
        if y[i] - base < thr:
            continue
        # parabolic refinement of the apex using neighbors
        denom = y[i - 1] - 2.0 * y[i] + y[i + 1]
        off = 0.0 if abs(denom) < 1e-12 else 0.5 * (y[i - 1] - y[i + 1]) / denom
        off = max(-0.5, min(0.5, off))
        tth_peak = float(t[i]) + off * (t[1] - t[0])
        height = float(y[i]) - base
        # rough FWHM by scanning for half-height crossings
        half = base + 0.5 * (y[i] - base)
        il = i
        while il > 0 and y[il] > half:
            il -= 1
        ir = i
        while ir < len(y) - 1 and y[ir] > half:
            ir += 1
        fwhm = max(float(t[ir] - t[il]), float(t[1] - t[0]))
        peaks.append(Peak(tth_peak, 0.0, height, fwhm))
    peaks.sort(key=lambda p: p.height, reverse=True)
    return peaks[:64]


def sample_fingerprint(pattern: PowderPattern,
                       prominence: float = PEAK_PROMINENCE) -> SampleFingerprint:
    """Build the sample fingerprint (peak list + d-space profile)."""
    wl = pattern.instrument.kalpha1 or 1.5406
    peaks = find_peaks(pattern.tth, pattern.intensity, prominence=prominence)
    for p in peaks:
        p.d = float(tth_to_d(p.tth, wl))
    # d-space profile on a fixed grid, index from d_min (high angle) upwards
    d_min = D_GRID_MIN
    d_max = max(2.0, float(tth_to_d(pattern.tth[0], wl)))
    grid = np.arange(d_min, d_max, D_GRID_STEP)
    # intensity is roughly linear in sin(theta) -> resample via sin(theta) axis
    sin_theta = np.sin(np.radians(pattern.tth) / 2.0)
    sin_grid = wl / (2.0 * grid)
    prof = np.interp(sin_grid, sin_theta, pattern.intensity, left=0.0, right=0.0)
    m = float(np.nanmax(prof))
    if m > 0:
        prof = prof / m
    return SampleFingerprint(source=pattern.source, peaks=peaks,
                             d_grid=grid, d_profile=prof)


def profile_similarity(a: SampleFingerprint, b: SampleFingerprint) -> float:
    """Cosine similarity of the normalized d-space profiles (0..1)."""
    if a.d_grid.size == 0 or b.d_grid.size == 0:
        return 0.0
    lo = max(a.d_grid[0], b.d_grid[0])
    hi = min(a.d_grid[-1], b.d_grid[-1])
    if hi <= lo:
        return 0.0
    # resample both profiles onto b's grid, restrict to shared d range
    g = b.d_grid
    pa = np.interp(g, a.d_grid, a.d_profile, left=0.0, right=0.0)
    pb = b.d_profile
    mask = (g >= lo + 1e-9) & (g <= hi - 1e-9)
    x, y = pa[mask], pb[mask]
    nx = math.sqrt(float(x @ x)) + 1e-30
    ny = math.sqrt(float(y @ y)) + 1e-30
    if nx <= 0 or ny <= 0:
        return 0.0
    return float((x @ y) / (nx * ny))


def match_peaks(a: SampleFingerprint, b: SampleFingerprint,
                tol: float = PEAK_D_TOL) -> Tuple[int, List[Tuple[float, float]]]:
    """Count d-spacing matches between two peak lists (greedy nearest)."""
    pb = [(p.d, p.height) for p in b.peaks]
    used = [False] * len(pb)
    n = 0
    pairs: List[Tuple[float, float]] = []
    for p in a.peaks:
        best, bi = None, -1
        for i, (d, _h) in enumerate(pb):
            if used[i]:
                continue
            if best is None or abs(d - p.d) < abs(best - p.d):
                best, bi = d, i
        if best is not None and abs(best - p.d) <= tol:
            used[bi] = True
            n += 1
            pairs.append((p.d, best))
    return n, pairs


def instrument_match(a: InstrumentFingerprint, b: InstrumentFingerprint) -> bool:
    """True when two observations come from an equivalent instrument setup."""
    return a.id == b.id


def instrument_equivalent(pa: InstrumentParams, pb: InstrumentParams,
                          wl_tol: float = 1e-4, range_tol: float = 0.05) -> bool:
    """Slightly looser comparison (same anode, wavelengths, axis; grid close)."""
    if pa.anode.strip().upper() != pb.anode.strip().upper():
        return False
    if len(pa.wavelengths) != len(pb.wavelengths):
        return False
    if any(abs(x - y) > wl_tol for x, y in zip(pa.wavelengths, pb.wavelengths)):
        return False
    if (abs(pa.tmin - pb.tmin) > range_tol or abs(pa.tmax - pb.tmax) > range_tol
            or pa.npts != pb.npts):
        return False
    return True