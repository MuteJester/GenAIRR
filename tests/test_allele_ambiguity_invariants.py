"""
Allele-call ambiguity invariants (Tier 1).

These tests verify two reactive properties that the engine's metadata
must satisfy at every pipeline stage:

    Invariant 1 — Allele set ⇄ evidence.
        For each segment (V/D/J), the comma-separated call list emitted
        by the engine must equal the *complete* set of reference alleles
        tied at the best god-aligner score on the observable bases inside
        ``[seg_sequence_start, seg_sequence_end)`` against germline
        coordinates ``[seg_germline_start, seg_germline_end)``.

        ``test_god_aligner.py::TestAlleleCallGodAligner`` already verifies
        that *every called allele* ties at max score (no false positives).
        This file additionally verifies *no allele tied at max is missing*
        from the call list (no false negatives) — i.e. completeness.

    Invariant 2 — Coordinates ⇄ evidence.
        The segment span ``[seg_start, seg_end)`` should be *maximal*:
        the engine's boundary-extension logic must not leave matchable
        adjacent bases unclaimed. We assert that extending the span by
        one base on either side would no longer be consistent with any
        allele in the called set.

The two invariants are tested under five event combinations so a
failure points cleanly at the offending pipeline stage:

    * pure recombination (no SHM, no corruption)
    * SHM only
    * SHM + indels
    * SHM + Ns
    * full pipeline (SHM + 5'-loss + indels + Ns)

If a stage breaks consistency, only the corresponding test fails.
"""

import pytest

from GenAIRR import Experiment
from GenAIRR.ops import rate, with_5prime_loss, with_indels, with_ns

from tests.test_god_aligner import (
    CONFIGS,
    _flatten_alleles,
    _has_d_segment,
    _score_allele,
)
from GenAIRR.protocol import _resolve_config


# ====================================================================
# Oracle
# ====================================================================

def _oracle_consistent_set(seq, allele_pool, seq_start, seq_end,
                            germ_start, germ_end):
    """Return the set of allele names tied at the best god-aligner score.

    Mirrors ``_score_allele`` semantics from ``test_god_aligner`` so the
    completeness check is consistent with the existing best-match check.

    Returns ``None`` (caller should skip) when:
      * the segment range is empty / invalid, or
      * ``(seq_end - seq_start) != (germ_end - germ_start)``: an indel
        sits inside the segment, breaking the index-aligned comparison
        below. A correct oracle for indel-affected segments would need
        the per-node ``germline_pos`` map (not exposed in the AIRR row)
        — out of scope for Tier 1, which targets non-indel cases.
    """
    if seq_end <= seq_start or germ_start < 0 or germ_end <= germ_start:
        return None
    if (seq_end - seq_start) != (germ_end - germ_start):
        return None
    scores = {}
    for name, allele_seq in allele_pool.items():
        s = _score_allele(seq, allele_seq, seq_start, seq_end,
                          germ_start, germ_end)
        if s is not None:
            scores[name] = s
    if not scores:
        return None
    best = max(scores.values())
    return {name for name, s in scores.items() if s == best}


# ====================================================================
# Fixtures — one record set per scenario
# ====================================================================

_SCENARIO_PARAMS = {
    "pure_recomb": dict(
        builder=lambda exp: exp,                       # no SHM, no corruption
        n=120, seed=42),
    "shm_only": dict(
        builder=lambda exp: exp.mutate(rate(0.02, 0.08)),
        n=120, seed=42),
    "shm_indels": dict(
        builder=lambda exp: (exp.mutate(rate(0.02, 0.08))
                                 .observe(with_indels(0.01))),
        n=120, seed=42),
    "shm_ns": dict(
        builder=lambda exp: (exp.mutate(rate(0.02, 0.08))
                                 .observe(with_ns(0.005))),
        n=120, seed=42),
    "full_pipeline": dict(
        builder=lambda exp: (exp.mutate(rate(0.02, 0.08))
                                 .sequence(with_5prime_loss())
                                 .observe(with_indels(0.01),
                                          with_ns(0.005))),
        n=120, seed=42),
}


@pytest.fixture(scope="module")
def all_pools():
    """Reference allele pools per config: {config: {seg: {name: seq}}}."""
    out = {}
    for name in CONFIGS:
        dc = _resolve_config(name)
        out[name] = {
            "v": _flatten_alleles(dc.v_alleles),
            "d": _flatten_alleles(dc.d_alleles),
            "j": _flatten_alleles(dc.j_alleles),
            "has_d": _has_d_segment(name),
        }
    return out


@pytest.fixture(scope="module", params=list(_SCENARIO_PARAMS.keys()))
def scenario_records(request):
    """Simulated records for one event-combination scenario across all
    configs. Keyed by config name, value is the SimulationResult."""
    params = _SCENARIO_PARAMS[request.param]
    builder = params["builder"]
    out = {}
    for name in CONFIGS:
        exp = builder(Experiment.on(name))
        out[name] = exp.run(n=params["n"], seed=params["seed"])
    out["__scenario__"] = request.param
    return out


# ====================================================================
# Invariant 1: call set is complete
# ====================================================================

def _assert_call_set_complete(records, pools, scenario, segment):
    """Across configs, assert engine's segment call set == oracle set."""
    seg_seq_start = f"{segment}_sequence_start"
    seg_seq_end = f"{segment}_sequence_end"
    seg_germ_start = f"{segment}_germline_start"
    seg_germ_end = f"{segment}_germline_end"
    seg_call = f"{segment}_call"

    n_checked = 0
    for config_name in CONFIGS:
        if segment == "d" and not pools[config_name]["has_d"]:
            continue
        pool = pools[config_name][segment]
        for i, rec in enumerate(records[config_name]):
            call = rec.get(seg_call) or ""
            # Short-D is a sentinel — the engine reports it when there are
            # ≤ 2 D bases left; no per-allele claim is being made, so the
            # oracle isn't applicable.
            if call == "Short-D" or not call:
                continue
            # Skip records with any indels: the AIRR row exposes only
            # totals (``n_insertions`` / ``n_deletions``), not their
            # positions, so the index-aligned oracle below cannot map
            # ``seq[ss + i]`` to the right germline coordinate when an
            # indel sits inside the segment — even when seg/germ net
            # lengths happen to coincide. Indel-aware completeness needs
            # per-node introspection (out of scope for Tier 1).
            if (rec.get("n_insertions", 0) or 0) > 0 \
               or (rec.get("n_deletions", 0) or 0) > 0:
                continue
            ss = rec[seg_seq_start]
            se = rec[seg_seq_end]
            gs = rec[seg_germ_start]
            ge = rec[seg_germ_end]
            oracle = _oracle_consistent_set(
                rec["sequence"], pool, ss, se, gs, ge)
            if oracle is None:
                continue
            engine = set(call.split(","))
            # Drop any names the engine reports that aren't in our pool
            # (defensive — should never happen for well-formed configs).
            engine = {n for n in engine if n in pool}
            if not engine:
                continue
            n_checked += 1
            missing = oracle - engine
            extra = engine - oracle
            assert not missing and not extra, (
                f"[{scenario}/{config_name} seq {i}] "
                f"{segment}_call set mismatch:\n"
                f"  engine: {sorted(engine)}\n"
                f"  oracle: {sorted(oracle)}\n"
                f"  missing from engine: {sorted(missing)}\n"
                f"  extra in engine: {sorted(extra)}\n"
                f"  range seq[{ss}:{se}] germ[{gs}:{ge}]")
    assert n_checked > 0, (
        f"[{scenario}] no {segment} records were validated — "
        f"fixture may be empty")


class TestCallSetCompleteness:
    """Invariant 1: call set equals the oracle's consistent set."""

    def test_v_call_set_complete(self, scenario_records, all_pools):
        scenario = scenario_records["__scenario__"]
        _assert_call_set_complete(scenario_records, all_pools,
                                   scenario, "v")

    def test_d_call_set_complete(self, scenario_records, all_pools):
        scenario = scenario_records["__scenario__"]
        _assert_call_set_complete(scenario_records, all_pools,
                                   scenario, "d")

    def test_j_call_set_complete(self, scenario_records, all_pools):
        scenario = scenario_records["__scenario__"]
        _assert_call_set_complete(scenario_records, all_pools,
                                   scenario, "j")


# ====================================================================
# Invariant 2 (coordinates are maximal): TODO
# ====================================================================
#
# The "no false negative" direction of Invariant 2 — that the engine's
# segment span cannot be extended one base further into the adjacent
# NP region while still being consistent with the called allele set —
# is intentionally NOT tested in this file.
#
# Reason: the engine's boundary-extension logic stops the walk when it
# encounters a P-nucleotide node (``NUC_FLAG_P_NUCLEOTIDE``) on the
# rationale that P-nucs have a deterministic palindromic origin and
# are "not aligner-ambiguous" with the trimmed germline tail. From an
# external-aligner perspective, however, a P-nuc is just a base in the
# flat sequence, indistinguishable from an N-nuc — so a true god-
# aligner would treat it the same.
#
# A test that fires when "you could extend by one more base" therefore
# generates a false alarm whenever the adjacent base is a P-nuc, even
# though the engine is behaving as designed. To distinguish the two
# cases the test would need per-node flag access (via the snapshot
# API), which is out of scope for the AIRR-row-based oracle here.
#
# The forward direction of Invariant 2 — every base inside
# ``[seg_start, seg_end)`` is consistent with at least one called
# allele at the corresponding germline position — is already covered
# by ``test_god_aligner.TestSegmentBoundaryAccuracy``.
#
# Open design question for the engine: should the boundary-extension
# walk stop at P-nucs (current behaviour) or treat them as ambiguous
# bases (god-aligner semantics)? This decision affects both the
# coordinate display (``airr.c::check_*_ambiguity``) and the call list
# (``allele_bitmap.c::allele_call_derive``) — they should remain
# consistent with each other.
