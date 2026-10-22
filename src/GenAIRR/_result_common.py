"""Shared truth-column and reference-binding helpers for SimulationResult
and its validation mixin. Kept here (imported by both) so there is one
definition and no result <-> _result_validation import cycle.
"""
from __future__ import annotations

from typing import Any, Dict, Optional


def _allele_name_or_empty(refdata: Any, segment: str, allele_id: Optional[int]) -> str:
    """Look up an allele name from refdata by id; return ``""`` when
    the id is None or the lookup fails (defensive).
    """
    if allele_id is None:
        return ""
    try:
        if segment == "V":
            return refdata.v_allele(int(allele_id)).name
        if segment == "D":
            return refdata.d_allele(int(allele_id)).name
        if segment == "J":
            return refdata.j_allele(int(allele_id)).name
    except Exception:
        return ""
    return ""


def _inject_truth_columns(outcome: Any, refdata: Any, record: Dict[str, Any]) -> None:
    """Append `truth_v_call` / `truth_d_call` / `truth_j_call`
    columns to ``record`` from the originally-sampled allele ids
    stored in the simulation's `assignments`. These compatibility
    aliases equal the standard `v_call` / `d_call` / `j_call` fields.
    Evidence-based calls are available as `recalled_*_call`.
    """
    sim = outcome.final_simulation()
    record["truth_v_call"] = _allele_name_or_empty(refdata, "V", sim.v_allele_id())
    record["truth_d_call"] = _allele_name_or_empty(refdata, "D", sim.d_allele_id())
    record["truth_j_call"] = _allele_name_or_empty(refdata, "J", sim.j_allele_id())


def _as_native_refdata(refdata: Any) -> Any:
    """Return the native configuration of a portable ``Reference``; any other
    value is returned unchanged."""
    from ._reference.artifact import Reference

    if isinstance(refdata, Reference):
        return refdata.native_refdata()
    return refdata


def _require_same_reference(bound_hash: Optional[str], refdata: Any) -> None:
    """Raise ``ValueError`` when a result bound to ``bound_hash`` is checked
    against a native reference with a different content hash.

    Unbound results (built from records alone) and values without a
    content hash are left to the caller's own checks.
    """
    if bound_hash is None:
        return
    content_hash = getattr(refdata, "content_hash", None)
    if content_hash is None:
        return
    given = content_hash()
    if given != bound_hash:
        raise ValueError(
            f"records were generated against a different reference (content hash {bound_hash}); "
            f"the reference given has content hash {given}. Validate against the reference "
            f"the records were generated with."
        )
