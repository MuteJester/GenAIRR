"""Explicit unfitted defaults profile.

Records, as data, the placeholder behaviour a simulation uses when no
fitted model is attached: uniform allele usage, uniform NP lengths and
bases, and no trimming or P-nucleotide addition. Nothing here invents
new modelling; it only makes today's implicit defaults visible and
hashable.
"""
from __future__ import annotations

from typing import Any, Dict, List

from GenAIRR._artifact_io import typed_digest

PROFILE = "genairr.unfitted_defaults.v1"
DIGEST_DOMAIN = "genairr_reference_models_v1"

_NP_SUPPORT = [0, 1, 2, 3, 4, 5, 6]
_BASES = ["A", "C", "G", "T"]


def _has_d(chain_type: str) -> bool:
    if chain_type not in ("vdj", "vj"):
        raise ValueError(f"chain_type must be 'vdj' or 'vj', got {chain_type!r}")
    return chain_type == "vdj"


def _label(allele_id: str) -> str:
    return allele_id.split(":", 1)[1] if ":" in allele_id else allele_id


def unfitted_defaults_v1(chain_type: str, view: Dict[str, List[str]]) -> dict:
    has_d = _has_d(chain_type)
    segments = ("V", "D", "J") if has_d else ("V", "J")
    nps = ("NP1", "NP2") if has_d else ("NP1",)
    ends = ("V_3", "D_5", "D_3", "J_5") if has_d else ("V_3", "J_5")
    usage: Dict[str, Any] = {"unit": "allele", "kind": "uniform"}
    for seg in segments:
        usage[seg] = [_label(a) for a in view[seg]]
    return {
        "profile": PROFILE,
        "allele_usage": usage,
        "np_lengths": {k: {"kind": "uniform", "support": list(_NP_SUPPORT)} for k in nps},
        "np_bases": {k: {"kind": "uniform", "support": list(_BASES)} for k in nps},
        "trims": {k: {"kind": "disabled"} for k in ends},
        "p_nucleotides": {k: {"kind": "disabled"} for k in ends},
    }


def models_hash(record: dict) -> str:
    return typed_digest(record, DIGEST_DOMAIN)


def recombine_defaults(record: dict, *, chain_type: str) -> dict:
    """Lower ``record`` to the dict :meth:`Experiment.recombine` consumes."""
    has_d = _has_d(chain_type)
    usage = record["allele_usage"]

    def weights(seg):
        if seg not in usage:
            return None
        return {label: 1.0 for label in usage[seg]}

    def np_len(key):
        spec = record["np_lengths"].get(key)
        return None if spec is None else [(v, 1.0) for v in spec["support"]]

    return {
        "np1": np_len("NP1"),
        "np2": np_len("NP2") if has_d else None,
        "np1_bases": None,
        "np2_bases": None,
        "np1_markov_transitions": None,
        "np2_markov_transitions": None,
        "p_v_3_lengths": None,
        "p_d_5_lengths": None,
        "p_d_3_lengths": None,
        "p_j_5_lengths": None,
        "trim_v_3": None,
        "trim_d_5": None,
        "trim_d_3": None,
        "trim_j_5": None,
        "allele_usage_v": weights("V"),
        "allele_usage_d": weights("D") if has_d else None,
        "allele_usage_j": weights("J"),
    }
