"""Tests for the redesigned GenAIRR MCP server.

Layered per the spec at docs/superpowers/specs/2026-05-18-mcp-redesign-v2-design.md:
- Tier 1: happy path per tool (14 tests)
- Tier 2: error envelope per error code (~6 tests)
- Tier 3: one end-to-end agent-style chain (1 test)

Plus a foundation block here for the envelope decorator itself.
"""
from __future__ import annotations

import pytest

from GenAIRR.mcp_errors import (
    CONFIG_NOT_FOUND,
    ALLELE_NOT_FOUND,
    INVALID_PRESET,
    INVALID_PARAMETER,
    MALFORMED_RECORD,
    SEED_REPLAY_MISMATCH,
    MCPError,
    envelope,
)


# -- Foundation: envelope decorator ---------------------------------


def test_envelope_wraps_success_with_ok_true_and_elapsed_ms():
    @envelope("test_tool")
    def my_tool() -> dict:
        return {"x": 1}

    resp = my_tool()
    assert resp["ok"] is True
    assert resp["tool"] == "test_tool"
    assert isinstance(resp["elapsed_ms"], int) and resp["elapsed_ms"] >= 0
    assert resp["result"] == {"x": 1}
    assert "error" not in resp


def test_envelope_catches_mcp_error_into_failure_shape():
    @envelope("test_tool")
    def my_tool() -> dict:
        raise MCPError(CONFIG_NOT_FOUND, "Config 'foo' missing.", hint="Try list_configs.")

    resp = my_tool()
    assert resp["ok"] is False
    assert resp["tool"] == "test_tool"
    assert isinstance(resp["elapsed_ms"], int) and resp["elapsed_ms"] >= 0
    assert resp["error"]["code"] == CONFIG_NOT_FOUND
    assert resp["error"]["message"] == "Config 'foo' missing."
    assert resp["error"]["hint"] == "Try list_configs."
    assert "result" not in resp


def test_envelope_catches_unexpected_exception_into_internal_error():
    @envelope("test_tool")
    def my_tool() -> dict:
        raise ValueError("something went wrong")

    resp = my_tool()
    assert resp["ok"] is False
    assert resp["error"]["code"] == "internal_error"
    assert "ValueError" in resp["error"]["message"]
    assert "something went wrong" in resp["error"]["message"]


def test_envelope_preserves_function_metadata():
    @envelope("test_tool")
    def my_tool(x: int) -> dict:
        """My docstring."""
        return {"x": x}

    assert my_tool.__name__ == "my_tool"
    assert my_tool.__doc__ == "My docstring."


def test_error_code_constants_are_stable_strings():
    # The agent branches on these -- they must be plain strings, not enums.
    assert CONFIG_NOT_FOUND == "config_not_found"
    assert ALLELE_NOT_FOUND == "allele_not_found"
    assert INVALID_PRESET == "invalid_preset"
    assert INVALID_PARAMETER == "invalid_parameter"
    assert MALFORMED_RECORD == "malformed_record"
    assert SEED_REPLAY_MISMATCH == "seed_replay_mismatch"


# -- Foundation: _mcp_summary helpers --------------------------------

from GenAIRR._mcp_summary import compute_repertoire_summary


def _fake_record(
    *,
    v_call="IGHV1*01",
    d_call="IGHD1*01",
    j_call="IGHJ1*01",
    productive=True,
    junction_aa="CARDW",
    n_mutations=0,
    n_pcr_errors=0,
    n_indels=0,
    n_quality_errors=0,
    sequence_length=300,
):
    return {
        "v_call": v_call, "d_call": d_call, "j_call": j_call,
        "productive": productive,
        "junction_aa": junction_aa,
        "n_mutations": n_mutations, "n_pcr_errors": n_pcr_errors,
        "n_indels": n_indels, "n_quality_errors": n_quality_errors,
        "sequence_length": sequence_length,
    }


def test_summary_v_usage_top_caps_at_top_n_and_sorts_by_count():
    records = (
        [_fake_record(v_call="IGHV1*01") for _ in range(10)]
        + [_fake_record(v_call="IGHV2*01") for _ in range(5)]
        + [_fake_record(v_call="IGHV3*01") for _ in range(3)]
        + [_fake_record(v_call="IGHV4*01") for _ in range(1)]
    )
    summary = compute_repertoire_summary(records, summary_top_n=2)
    # top 2 only, in descending count order
    assert list(summary["v_usage_top"].items()) == [("IGHV1*01", 10), ("IGHV2*01", 5)]
    assert summary["n_unique_v"] == 4


def test_summary_cdr3_histogram_keys_are_stringified_ints():
    records = [_fake_record(junction_aa="CARDW") for _ in range(3)]  # length 5
    records += [_fake_record(junction_aa="CARDIW") for _ in range(2)]  # length 6
    summary = compute_repertoire_summary(records, summary_top_n=10)
    assert summary["cdr3_length_histogram"] == {"5": 3, "6": 2}
    assert all(isinstance(k, str) for k in summary["cdr3_length_histogram"])


def test_summary_productive_rate_counts_only_true():
    records = (
        [_fake_record(productive=True) for _ in range(7)]
        + [_fake_record(productive=False) for _ in range(2)]
        + [_fake_record(productive=None) for _ in range(1)]
    )
    summary = compute_repertoire_summary(records, summary_top_n=10)
    assert summary["productive_count"] == 7
    assert summary["productive_rate"] == 0.7


def test_summary_mutation_quartiles_present_when_any_record_has_mutations():
    records = [_fake_record(n_mutations=k) for k in range(0, 21)]
    summary = compute_repertoire_summary(records, summary_top_n=10)
    # quartiles[0,1,2,3,4] = min, q1, median, q3, max
    assert summary["mutation_count_quartiles"][0] == 0
    assert summary["mutation_count_quartiles"][4] == 20
    assert summary["mutation_count_quartiles"][2] == 10  # median


def test_summary_omits_corruption_quartiles_when_all_zero():
    records = [_fake_record() for _ in range(5)]  # all zeros
    summary = compute_repertoire_summary(records, summary_top_n=10)
    assert "n_pcr_errors_quartiles" not in summary
    assert "n_indels_quartiles" not in summary


def test_summary_handles_empty_records():
    summary = compute_repertoire_summary([], summary_top_n=10)
    assert summary["productive_count"] == 0
    assert summary["productive_rate"] == 0.0
    assert summary["v_usage_top"] == {}
    assert summary["n_unique_v"] == 0
    assert summary["cdr3_length_histogram"] == {}


# -- Foundation: _mcp_validators ------------------------------------

from GenAIRR._mcp_validators import validate_one_record


def test_validator_required_fields_present_passes():
    rec = _fake_record()
    issues = validate_one_record(rec, refdata=None)
    assert issues == []


def test_validator_missing_sequence_field_flags():
    rec = _fake_record()
    del rec["v_call"]
    issues = validate_one_record(rec, refdata=None)
    assert any("missing required field" in i and "v_call" in i for i in issues)


def test_validator_junction_bounds_out_of_range_flags():
    rec = _fake_record()
    rec["sequence_length"] = 200
    rec["junction_start"] = 180
    rec["junction_end"] = 250  # > sequence_length
    issues = validate_one_record(rec, refdata=None)
    assert any("junction_end" in i and "out of range" in i for i in issues)


def test_validator_productive_with_stop_codon_in_junction_flags():
    rec = _fake_record(productive=True, junction_aa="CAR*W")
    issues = validate_one_record(rec, refdata=None)
    assert any("stop" in i.lower() for i in issues)


def test_validator_cigar_with_soft_clip_op_flags():
    rec = _fake_record()
    rec["v_cigar"] = "5S280M"
    issues = validate_one_record(rec, refdata=None)
    assert any("v_cigar" in i and "S" in i for i in issues)


def test_validator_refdata_required_checks_skipped_without_config():
    # Tests that checks 6-9 are skipped (not flagged) when refdata=None.
    rec = _fake_record(v_call="IGHV_DOES_NOT_EXIST*01")
    issues = validate_one_record(rec, refdata=None)
    assert not any("allele not found" in i.lower() for i in issues)
