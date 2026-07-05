"""Shared truth-column helpers for SimulationResult and its validation
mixin. Kept here (imported by both) so there is one definition and no
result <-> _result_validation import cycle.
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
    stored in the simulation's `assignments`. Distinct from
    `v_call` / `d_call` / `j_call`, which are evidence-driven and
    can change under heavy SHM.
    """
    sim = outcome.final_simulation()
    record["truth_v_call"] = _allele_name_or_empty(refdata, "V", sim.v_allele_id())
    record["truth_d_call"] = _allele_name_or_empty(refdata, "D", sim.d_allele_id())
    record["truth_j_call"] = _allele_name_or_empty(refdata, "J", sim.j_allele_id())
