"""Data models for ingested powder diffraction data (XRDML and beyond).

Kept dependency-light (numpy only) so downstream modules (calibration,
hypothesis, verdict) never need to import GSAS-II just for the core types.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Tuple

import numpy as np


@dataclass(frozen=True)
class InstrumentParams:
    """Normalized instrument parameters extracted from a data file."""

    anode: str
    wavelengths: Tuple[float, ...]          # K alpha1, K alpha2, ...
    scan_type: str = "continuous"
    scan_axis: str = "2Theta/Theta"
    tmin: float = 0.0
    tmax: float = 0.0
    step: float = 0.0
    npts: int = 0

    def canonical(self) -> str:
        """Canonical string used to build the instrument fingerprint id."""
        wl = ",".join(f"{w:.6f}" for w in self.wavelengths)
        return "|".join(
            [
                self.anode.strip().upper(),
                wl,
                self.scan_axis.strip().upper(),
                f"{self.tmin:.4f}",
                f"{self.tmax:.4f}",
                f"{self.step:.6f}",
                str(self.npts),
            ]
        )

    @property
    def fingerprint_id(self) -> str:
        """Stable id: equal iff the canonical instrument description is equal."""
        return hashlib.sha1(self.canonical().encode("utf-8")).hexdigest()[:16]

    @property
    def kalpha1(self) -> float:
        return self.wavelengths[0] if self.wavelengths else 0.0


@dataclass
class PowderPattern:
    """A parsed powder pattern: 2-theta grid + measured intensities + instrument."""

    sample_name: str
    source: str
    tth: np.ndarray                    # 2theta in degrees
    intensity: np.ndarray              # raw counts
    instrument: InstrumentParams
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        self.tth = np.asarray(self.tth, dtype=float)
        self.intensity = np.asarray(self.intensity, dtype=float)
        if self.tth.shape != self.intensity.shape:
            raise ValueError("tth and intensity must have identical shape")
        if self.tth.size == 0:
            raise ValueError("pattern is empty")

    @property
    def npts(self) -> int:
        return int(self.tth.size)

    @property
    def peak_max(self) -> float:
        return float(np.nanmax(self.intensity))