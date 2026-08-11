"""Shared COD catalog helpers (Spike 06): CIF parsing, chemistry validation,
formula normalization, cell extraction. numpy-free, importable by tests
without GSAS-II.
"""
from __future__ import annotations

import hashlib
import re
from typing import Dict, List, Optional, Tuple

_ELEMENT_GROUP = re.compile(r"([A-Z][a-z]?)(\d*\.?\d*)")


def formula_dict(formula: str) -> Dict[str, float]:
    """'Ag0.5 Bi0.5 S' / 'O4 Pb S' -> {'Ag': 0.5, 'Bi': 0.5, 'S': 1.0, ...}."""
    out: Dict[str, float] = {}
    for el, num in _ELEMENT_GROUP.findall(formula.replace(" ", "")):
        out[el] = out.get(el, 0.0) + float(num or 1.0)
    return out


def ratios_close(actual: Dict[str, float], expected: Dict[str, float],
                 tol: float = 0.02) -> bool:
    """Element set identical and all ratios within tol after scaling."""
    if set(actual) != set(expected):
        return False
    key = max(expected, key=expected.get)
    scale = expected[key] / actual[key]
    for el in expected:
        if abs(expected[el] - actual[el] * scale) > tol * max(
                expected[el], actual[el] * scale):
            return False
    return True


def chemistry_matches(formula: str, expected: Dict[str, float],
                      tol: float = 0.02) -> bool:
    return ratios_close(formula_dict(formula), expected, tol)


def cif_field(text: str, pattern: str) -> Optional[str]:
    m = re.search(pattern + r"\s+(.+)", text)
    if not m:
        return None
    return m.group(1).strip().strip('"\'') or None


def _to_float(v: Optional[str]):
    if v is None:
        return None
    try:
        return float(re.sub(r"\([^)]*\)", "", v).strip())
    except ValueError:
        return None


def parse_cif(text: str) -> dict:
    """Extract identity fields from a CIF (COD-flavoured)."""
    cell = {}
    for k in "abc":
        cell[k] = _to_float(cif_field(text, rf"_cell_length_{k}"))
    for k in ("alpha", "beta", "gamma"):
        cell[k] = _to_float(cif_field(text, rf"_cell_angle_{k}"))
    authors = [a for a in re.findall(r"_publ_author_name\s+(.+)", text)
               if a.strip()]
    title = cif_field(text, r"_publ_section_title")
    journal = cif_field(text, r"_journal_name_full")
    year = cif_field(text, r"_journal_year")
    formula = (cif_field(text, r"_chemical_formula_sum")
               or cif_field(text, r"_chemical_formula_structural")
               or "?")
    sg = (cif_field(text, r"_symmetry_space_group_name_H-M")
          or cif_field(text, r"_space_group_name_H-M_alt") or "?")
    return {
        "formula": formula,
        "space_group": sg,
        "cell": cell,
        "citation": _citation(authors, title, journal, year),
        "title": title,
        "authors": authors,
        "year": year,
    }


def _citation(authors: List[str], title: Optional[str],
              journal: Optional[str], year: Optional[str]) -> str:
    parts = []
    if authors:
        parts.append(", ".join(authors))
    if title:
        parts.append(title)
    bits = [b for b in (journal, year) if b]
    if bits:
        parts.append(", ".join(bits))
    return "; ".join(parts) if parts else "no publication metadata in CIF"


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()