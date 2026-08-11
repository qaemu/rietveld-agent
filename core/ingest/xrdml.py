"""XRDML (Malvern Panalytical X'Pert XML) reader.

Parses the common XRDML 1.x structure while being tolerant of namespace
declarations and element ordering:

  xrdMeasurements
    xrdMeasurement (measurementType, status)
      sample > name
      scan
        xrdElement
          commonParameters
            scanType, tubeAnode
            wavelengths > kAlpha1, kAlpha2, kBeta1
            goniometer > thetaMin, thetaMax, stepSize, scanAxis
          intensities   (space separated ASCII, or base64)

Only stdlib + numpy are used; GSAS-II is not required to read data files.
"""
from __future__ import annotations

import base64
import io
import xml.etree.ElementTree as ET
from typing import Optional

import numpy as np

from .models import InstrumentParams, PowderPattern


class XRDMLError(Exception):
    """Raised for malformed or unsupported XRDML content."""


def _local(tag: str) -> str:
    """Strip an XML namespace from a tag name."""
    return tag.rsplit("}", 1)[-1]


def _first_text(elem: Optional[ET.Element], name: str) -> Optional[str]:
    """Return the text of the first descendant with local-name ``name``."""
    if elem is None:
        return None
    for el in elem.iter():
        if _local(el.tag) == name and el.text and el.text.strip():
            return el.text.strip()
    return None


def _first_float(elem: Optional[ET.Element], name: str) -> Optional[float]:
    txt = _first_text(elem, name)
    if txt is None:
        return None
    try:
        return float(txt)
    except ValueError as exc:
        raise XRDMLError(f"non-numeric value for <{name}>: {txt!r}") from exc


def _descendant_float_map(elem: ET.Element, parent: str) -> dict:
    """Read all numeric children of the first <parent> block (e.g. wavelengths)."""
    out: dict[str, float] = {}
    for el in elem.iter():
        children = list(el)
        names = [_local(c.tag) for c in children]
        if _local(el.tag) == parent or parent in names:
            continue
    # simpler: iterate all elements, key floats by local name
    for el in elem.iter():
        nm = _local(el.tag)
        if el.text and el.text.strip():
            try:
                out[nm] = float(el.text.strip())
            except ValueError:
                pass
    return out


def _parse_intensities(elem: ET.Element) -> np.ndarray:
    for el in elem.iter():
        if _local(el.tag) != "intensities":
            continue
        text = (el.text or "").strip()
        if not text:
            raise XRDMLError("empty <intensities> element")
        enc = (el.get("encoding") or "").lower()
        if enc in ("base64", "base-64"):
            try:
                text = base64.b64decode(text).decode("ascii", errors="strict")
            except Exception as exc:  # binascii.Error, UnicodeDecodeError
                raise XRDMLError(f"invalid base64 intensities: {exc}") from exc
        try:
            vals = np.fromstring(text, sep=" ")
        except Exception as exc:
            raise XRDMLError(f"cannot parse intensities: {exc}") from exc
        if vals.size == 0:
            raise XRDMLError("intensities contain no numeric values")
        return vals
    raise XRDMLError("no <intensities> element found")


def parse_xrdml(source: str | bytes | io.IOBase, name: Optional[str] = None) -> PowderPattern:
    """Parse XRDML content into a :class:`PowderPattern`.

    :param source: file path (str), raw bytes, or an open file-like object.
    :param name: display name override (defaults to the <sample><name> or path).
    """
    if isinstance(source, (str, io.IOBase)):
        data = source
        label = name or getattr(source, "name", "<xrdml>")
        if isinstance(source, str) and not source.lstrip().startswith("<"):
            try:
                with open(source, "r", encoding="utf-8") as fh:
                    tree = ET.parse(fh)
            except FileNotFoundError:
                raise XRDMLError(f"file not found: {source}")
            except ET.ParseError as exc:
                raise XRDMLError(f"malformed XML in {source}: {exc}") from exc
            label = name or source
        else:
            try:
                tree = ET.parse(io.StringIO(source if isinstance(source, str) else source.read()))
            except ET.ParseError as exc:
                raise XRDMLError(f"malformed XML: {exc}") from exc
    else:  # bytes
        try:
            tree = ET.parse(io.BytesIO(source))
        except ET.ParseError as exc:
            raise XRDMLError(f"malformed XML: {exc}") from exc
        label = name or "<xrdml>"

    root = tree.getroot()
    if _local(root.tag) != "xrdMeasurements":
        raise XRDMLError(f"root element is <{_local(root.tag)}>, expected <xrdMeasurements>")

    measurement = next(
        (el for el in root.iter() if _local(el.tag) == "xrdMeasurement"), root
    )

    sample_name = _first_text(measurement, "name") or label
    anode = _first_text(measurement, "tubeAnode") or "?"
    scan_type = _first_text(measurement, "scanType") or "continuous"
    scan_axis = _first_text(measurement, "scanAxis") or "2Theta/Theta"

    tmin = _first_float(measurement, "thetaMin")
    tmax = _first_float(measurement, "thetaMax")
    step = _first_float(measurement, "stepSize")
    tube_voltage = _first_float(measurement, "tubeVoltage")
    tube_current = _first_float(measurement, "tubeCurrent")

    wl = _descendant_float_map(measurement, "wavelengths")
    wavelengths = tuple(
        v for k, v in sorted(wl.items())
        if k.lower() in ("kalpha1", "kalpha2", "kbeta1")
    )

    intensities = _parse_intensities(measurement)
    npts = len(intensities)

    if step is not None and tmin is not None:
        tth = tmin + step * np.arange(npts)
    elif tmin is not None and tmax is not None and npts > 1:
        tth = np.linspace(tmin, tmax, npts)
    else:
        raise XRDMLError(
            "cannot construct 2theta axis: need (thetaMin+stepSize) or (thetaMin+thetaMax)"
        )

    instrument = InstrumentParams(
        anode=anode,
        wavelengths=wavelengths,
        scan_type=scan_type,
        scan_axis=scan_axis,
        tmin=float(tmin) if tmin is not None else float(tth[0]),
        tmax=float(tmax) if tmax is not None else float(tth[-1]),
        step=float(step) if step is not None else float(tth[1] - tth[0]),
        npts=npts,
    )
    return PowderPattern(
        sample_name=sample_name,
        source=label,
        tth=tth,
        intensity=intensities,
        instrument=instrument,
        metadata={"xrdml.measurement_type": _first_text(measurement, "xrdMeasurement") or None,
                  "xrdml.tube_voltage_kv": tube_voltage,
                  "xrdml.tube_current_ma": tube_current},
    )