"""GSAS-II instrument parameter file (PRM) parser.

Reads the fixed-field text format used by GSAS-II (and GSAS) instruments and
produces structured data compatible with :class:`core.ingest.InstrumentParams`,
plus a sha256 so registry content references the exact bytes.
"""
from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass
from typing import Optional

from core.ingest import InstrumentParams

#: GSAS-II 'INS  1 IRAD' codes -> anode material
IRAD_ANODE = {0: "?", 1: "CrKa", 2: "FeKa", 3: "CuKa", 4: "MoKa",
              5: "AgKa", 6: "TiKa", 7: "CoKa"}

#: GSAS-II 'INS   HTYPE ' codes -> instrument type (first 3 chars)
HTYPE_NAMES = {"PXC": "Bragg-Brentano", "PNC": "Neutron-PositionSensitive",
               "TXC": "TOF-Xray", "TNC": "TOF-Neutron", "TTW": "TOF-EngDiff"}


class PrmError(Exception):
    """Raised for malformed or unsupported PRM content."""


def _field(line: str, key: str) -> Optional[str]:
    """Return the remainder of a line starting with the fixed-width key."""
    if not line.startswith(key):
        return None
    return line[len(key):].strip()


def _anode_from_irad(irad: str) -> str:
    try:
        return IRAD_ANODE[int(irad)]
    except (ValueError, KeyError):
        return "?"


@dataclass(frozen=True)
class PrmContent:
    """Structured representation of a GSAS-II PRM file."""

    source_path: str
    sha256: str
    instrument_type: str                 # e.g. 'PXCR'
    anode: str                           # e.g. 'CuKa'
    wavelengths: tuple[float, ...]       # (kalpha1, kalpha2, ...)
    profile_function: str = ""
    profile_coefficients: tuple[float, ...] = ()

    def as_instrument_params(self, tmin: float = 0.0, tmax: float = 0.0,
                             step: float = 0.0, npts: int = 0) -> InstrumentParams:
        return InstrumentParams(
            anode=self.anode, wavelengths=self.wavelengths,
            tmin=tmin, tmax=tmax, step=step, npts=npts,
        )


def parse_prm(path: str) -> PrmContent:
    """Parse a GSAS-II instrument parameter file.

    :param path: path to the .prm file
    :returns: :class:`PrmContent`
    """
    if not os.path.exists(path):
        raise PrmError(f"file not found: {path}")
    agg = hashlib.sha256()
    lines: list[str] = []
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            agg.update(chunk)
    with open(path, "r", errors="replace") as fh:
        lines = fh.readlines()

    htype = "?"
    irad = "0"
    icons: list[str] = []
    prcf_type = ""
    prcf_coeffs: list[str] = []
    for line in lines:
        s = _field(line, "INS   HTYPE ")
        if s:
            htype = s.split()[0][:3].upper()
            continue
        s = _field(line, "INS  1 IRAD")
        if s is not None:
            irad = s.split()[0]
            continue
        s = _field(line, "INS  1 ICONS")
        if s:
            icons = s.split()
            continue
        # PRCF1 has numeric-suffixed siblings PRCF11/PRCF12 (profile
        # coefficient blocks) -- check the specific keys first.
        if line.startswith("INS  1PRCF12"):
            rest = line[len("INS  1PRCF12"):]
            if rest and rest[0].isspace():
                prcf_coeffs.extend(rest.split())        # row index was part of key
            continue
        if line.startswith("INS  1PRCF11"):
            rest = line[len("INS  1PRCF11"):]
            if rest and rest[0].isspace():
                prcf_coeffs.extend(rest.split())
            continue
        if line.startswith("INS  1PRCF1"):
            rest = line[len("INS  1PRCF1"):]
            if (not rest or rest[0].isspace()) and not prcf_type:
                prcf_type = rest.split()[0]
            continue

    if not icons:
        raise PrmError(f"{path}: no 'INS  1 ICONS' line (wavelengths required)")
    try:
        wl = tuple(float(x) for x in icons[:2])
    except ValueError as exc:
        raise PrmError(f"{path}: non-numeric wavelengths in ICONS: {icons[:2]}") from exc
    if not wl[0] > 0:
        raise PrmError(f"{path}: invalid wavelength {wl[0]}")

    try:
        profile_coeffs = tuple(float(x) for x in prcf_coeffs[:8])
    except ValueError:
        profile_coeffs = ()

    return PrmContent(
        source_path=os.path.abspath(path),
        sha256=agg.hexdigest(),
        instrument_type=htype,
        anode=_anode_from_irad(irad),
        wavelengths=wl,
        profile_function=prcf_type,
        profile_coefficients=profile_coeffs,
    )


def instrument_type_name(htype: str) -> str:
    return HTYPE_NAMES.get(htype[:3], htype)