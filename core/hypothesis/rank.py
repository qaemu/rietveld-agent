"""Hypothesis stage: rank candidate materials against an observed pattern.

Candidate matching uses the d-space profile cosine (validated in unit 02:
0.908 same material across anodes vs <=0.05 for distinct materials) plus a
peak-line-list cross-check for reporting. Deterministic, numpy only, no
GSAS-II at runtime.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

from core.ingest import SampleFingerprint, match_peaks, profile_similarity


@dataclass
class CandidateMatch:
    material_id: str
    name: str
    similarity: float
    matched_peaks: int
    phase_family: str = ""
    top_peaks: List[dict] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {"material_id": self.material_id, "name": self.name,
                "phase_family": self.phase_family,
                "similarity": round(float(self.similarity), 6),
                "matched_peaks": int(self.matched_peaks),
                "top_peaks": self.top_peaks[:5]}


@dataclass
class HypothesisRanking:
    query_source: str
    n_candidates: int
    ranked: List[CandidateMatch]
    top_similarity: float = 0.0
    family_margin: float = 0.0      # top vs best OTHER phase family
    margin: float = 0.0             # top vs second best (any candidate)

    def to_dict(self) -> dict:
        return {"query_source": self.query_source,
                "n_candidates": self.n_candidates,
                "top_similarity": round(float(self.top_similarity), 6),
                "family_margin": round(float(self.family_margin), 6),
                "margin": round(float(self.margin), 6),
                "ranked": [c.to_dict() for c in self.ranked]}


def load_library(entries: list[dict],
                 families: Optional[dict[str, str]] = None) -> dict[str, SampleFingerprint]:
    """Turn library entries (SampleFingerprint.to_dict payloads) into
    fingerprints keyed by material id."""
    out: dict[str, SampleFingerprint] = {}
    for e in entries:
        out[e["id"]] = SampleFingerprint.from_dict(e["fingerprint"])
    return out


def rank_candidates(query: SampleFingerprint,
                    library: dict[str, SampleFingerprint],
                    names: Optional[dict[str, str]] = None,
                    families: Optional[dict[str, str]] = None,
                    top_k: int = 5) -> HypothesisRanking:
    """Rank all library materials by d-space profile similarity against the
    query fingerprint. Deterministic; ties broken by id for stability.

    ``families`` maps material id -> phase family. When provided, the ranking
    also computes ``family_margin`` (top similarity minus the best similarity
    among *different* phase families, which is the decision-relevant margin:
    statuses are per phase family, and differently-realized entries of the
    same family must not compete against each other). Without ``families``,
    ``family_margin`` falls back to the plain top-vs-second margin.
    """
    scored: List[CandidateMatch] = []
    for mid, ref in library.items():
        sim = profile_similarity(query, ref)
        n_matched, _pairs = match_peaks(query, ref)
        scored.append(CandidateMatch(
            material_id=mid,
            name=(names or {}).get(mid, mid),
            phase_family=(families or {}).get(mid, ""),
            similarity=float(sim),
            matched_peaks=int(n_matched),
            top_peaks=ref.peak_line_list(5),
        ))
    scored.sort(key=lambda c: (-c.similarity, c.material_id))
    top = scored[:top_k]
    top_sim = top[0].similarity if top else 0.0
    margin = (top[0].similarity - top[1].similarity) if len(top) > 1 else top_sim
    family_margin = margin
    if top and families:
        top_fam = top[0].phase_family
        others = [c.similarity for c in top[1:] if c.phase_family != top_fam]
        family_margin = (top_sim - max(others)) if others else top_sim
    return HypothesisRanking(
        query_source=query.source,
        n_candidates=len(scored),
        ranked=top,
        top_similarity=top_sim,
        family_margin=float(family_margin),
        margin=float(margin),
    )