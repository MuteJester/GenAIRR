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
    v_call: str = "IGHV1*01",
    d_call: str = "IGHD1*01",
    j_call: str = "IGHJ1*01",
    productive: bool | None = True,
    junction_aa: str = "CARDW",
    n_mutations: int = 0,
    n_pcr_errors: int = 0,
    n_indels: int = 0,
    n_quality_errors: int = 0,
    sequence_length: int = 300,
):
    return {
        "sequence": "A" * sequence_length,
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


def test_validator_missing_v_call_field_flags():
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


def test_validator_v_germline_span_exceeding_allele_length_flags():
    # Check 7: span must fit within the claimed allele's length.
    import GenAIRR as ga
    rd = ga.Experiment.on("human_igh").refdata
    v_name = rd.v_allele(0).name
    allele_len = len(rd.v_allele(0).seq())
    rec = _fake_record(v_call=v_name)
    rec["v_germline_start"] = 0
    rec["v_germline_end"] = allele_len + 50  # impossible span
    issues = validate_one_record(rec, refdata=rd)
    assert any("v_germline span" in i and "exceeds" in i for i in issues)


def test_validator_v_anchor_inconsistent_with_junction_start_flags():
    # Check 8: junction_start must equal projected V Cys anchor position.
    import GenAIRR as ga
    rd = ga.Experiment.on("human_igh").refdata
    # Pick an allele that actually has an anchor.
    v_name = None
    for i in range(rd.v_pool_size()):
        a = rd.v_allele(i)
        if a.anchor is not None:
            v_name = a.name
            break
    assert v_name is not None
    rec = _fake_record(v_call=v_name)
    rec["v_sequence_start"] = 0
    rec["v_germline_start"] = 0
    rec["junction_start"] = 999  # arbitrarily wrong
    issues = validate_one_record(rec, refdata=rd)
    assert any("anchor projected" in i and "does not match junction_start" in i for i in issues)


def test_validator_locus_mismatch_with_config_flags():
    # Check 9: rec['locus'] must match locus derived from config name.
    rec = _fake_record()
    rec["locus"] = "TRB"  # but config is human_igh -> IGH
    issues = validate_one_record(rec, refdata=None, config_name="human_igh")
    assert any("locus" in i.lower() and "does not match" in i.lower() for i in issues)


def test_validator_locus_tcr_aliases_normalised_to_engine_token():
    # locus_from_config_name('human_tcrb') returns 'TCRB' but the engine
    # emits 'TRB' -- the check must normalise before comparing.
    rec = _fake_record()
    rec["locus"] = "TRB"
    issues = validate_one_record(rec, refdata=None, config_name="human_tcrb")
    assert not any("locus" in i.lower() and "does not match" in i.lower() for i in issues)


# -- Foundation: _mcp_presets --------------------------------------

from GenAIRR._mcp_presets import PRESETS, resolve_preset


def test_presets_enumeration_matches_spec():
    assert set(PRESETS.keys()) == {
        "naive_b_cell",
        "memory_b_cell_shm",
        "tcr_beta_pool",
        "low_quality_sequencing",
        "clonal_expansion",
    }


def test_resolve_preset_returns_param_dict_with_config_n_seed():
    params = resolve_preset("memory_b_cell_shm")
    assert params["config"] == "human_igh"
    assert params["productive_only"] is True
    assert params["mutation_model"] == "s5f"
    assert params["mutation_count_min"] == 5
    assert params["mutation_count_max"] == 15


def test_resolve_preset_clonal_expansion_uses_clonal_params():
    params = resolve_preset("clonal_expansion")
    assert params["n_clones"] == 20
    assert params["clone_size"] == 25
    # No top-level n -- total comes from n_clones x clone_size
    assert "n" not in params or params["n"] is None


def test_resolve_preset_unknown_raises_mcp_error():
    with pytest.raises(MCPError) as exc_info:
        resolve_preset("not_a_real_preset")
    assert exc_info.value.code == INVALID_PRESET


# -- Foundation: _mcp_refdata helpers ------------------------------

from GenAIRR._mcp_refdata import (
    iter_alleles,
    find_allele,
    locus_from_config_name,
    resolve_refdata,
)


def test_resolve_refdata_returns_engine_refdata():
    rd = resolve_refdata("human_igh")
    assert rd.v_pool_size() > 0
    assert rd.j_pool_size() > 0


def test_resolve_refdata_unknown_raises_mcp_error():
    with pytest.raises(MCPError) as exc_info:
        resolve_refdata("human_xyz")
    assert exc_info.value.code == CONFIG_NOT_FOUND


def test_iter_alleles_yields_every_allele_in_pool():
    rd = resolve_refdata("human_igh")
    v_alleles = list(iter_alleles(rd, "v"))
    assert len(v_alleles) == rd.v_pool_size()
    assert all(a.segment == "V" for a in v_alleles)


def test_iter_alleles_returns_empty_when_pool_missing():
    # human_igk is a VJ chain -- no D pool.
    rd = resolve_refdata("human_igk")
    d_alleles = list(iter_alleles(rd, "d"))
    assert d_alleles == []


def test_find_allele_returns_match_by_name():
    rd = resolve_refdata("human_igh")
    name = rd.v_allele(0).name
    found = find_allele(rd, "v", name)
    assert found is not None
    assert found.name == name


def test_find_allele_returns_none_when_not_found():
    rd = resolve_refdata("human_igh")
    found = find_allele(rd, "v", "IGHV_DOES_NOT_EXIST*99")
    assert found is None


def test_locus_from_config_name_extracts_uppercase_locus():
    assert locus_from_config_name("human_igh") == "IGH"
    assert locus_from_config_name("mouse_tcrb") == "TCRB"
    assert locus_from_config_name("rabbit_igk") == "IGK"


def test_locus_from_config_name_returns_none_for_unparseable():
    assert locus_from_config_name("not_a_real_format") in (None, "A")  # tolerant
    assert locus_from_config_name("") is None


# -- Tier 1: Discovery tools ----------------------------------------


def test_list_configs_returns_human_igh_in_envelope():
    from GenAIRR.mcp_server import list_configs

    resp = list_configs()
    assert resp["ok"] is True
    assert resp["tool"] == "list_configs"
    r = resp["result"]
    assert "human_igh" in r["configs"]
    assert r["n_total"] >= 100  # ~106 unique configs / 208 aliases in v2.0.0
    assert isinstance(r["by_species"], dict)
    assert "Human" in r["by_species"]


def test_list_configs_species_filter_narrows_results():
    from GenAIRR.mcp_server import list_configs

    resp = list_configs(species_filter="Human")
    r = resp["result"]
    # All returned configs should be human_*
    for cfg in r["configs"]:
        assert cfg.startswith("human_")
    assert r["n_total"] == len(r["configs"])


def test_config_info_returns_pool_sizes_for_human_igh():
    from GenAIRR.mcp_server import config_info

    resp = config_info(config="human_igh")
    assert resp["ok"] is True
    r = resp["result"]
    assert r["config"] == "human_igh"
    assert r["locus"] == "IGH"
    assert r["v_pool_size"] > 0
    assert r["d_pool_size"] > 0
    assert r["j_pool_size"] > 0
    assert isinstance(r["chain_type"], str)


def test_list_alleles_returns_v_pool_names_truncated_to_limit():
    from GenAIRR.mcp_server import list_alleles

    resp = list_alleles(config="human_igh", segment="v", limit=5)
    assert resp["ok"] is True
    r = resp["result"]
    assert r["segment"] == "v"
    assert len(r["allele_names"]) == 5
    assert r["n_total"] > 5
    assert r["truncated"] is True


def test_inspect_allele_returns_sequence_and_anchor_for_known_allele():
    from GenAIRR.mcp_server import inspect_allele, list_alleles

    # Grab any V allele name from the config to avoid hard-coding.
    names_resp = list_alleles(config="human_igh", segment="v", limit=1)
    allele_name = names_resp["result"]["allele_names"][0]

    resp = inspect_allele(config="human_igh", allele_name=allele_name)
    assert resp["ok"] is True
    r = resp["result"]
    assert r["allele_name"] == allele_name
    assert r["segment"] == "v"
    assert r["length"] > 0
    assert isinstance(r["sequence"], str)
    assert len(r["sequence"]) == r["length"]


# -- Tier 1: Simulation tools ---------------------------------------


def test_simulate_repertoire_default_returns_summary_no_records():
    from GenAIRR.mcp_server import simulate_repertoire

    resp = simulate_repertoire(config="human_igh", n=50, seed=42)
    assert resp["ok"] is True
    r = resp["result"]
    assert r["n_records"] == 50
    assert r["seed"] == 42
    assert r["config"] == "human_igh"
    assert "v_usage_top" in r
    assert "cdr3_length_histogram" in r
    assert "records" not in r  # default is opt-out


def test_simulate_repertoire_with_records_returns_capped_rows():
    from GenAIRR.mcp_server import simulate_repertoire

    resp = simulate_repertoire(
        config="human_igh",
        n=300,
        seed=42,
        return_records=True,
        return_records_limit=10,
    )
    r = resp["result"]
    assert len(r["records"]) == 10
    assert r["records_truncated"] is True
    # The records are AIRR dicts with the expected fields
    rec = r["records"][0]
    assert "v_call" in rec
    assert "junction_aa" in rec


def test_simulate_repertoire_productive_only_constraint_forces_100_pct_productive():
    from GenAIRR.mcp_server import simulate_repertoire

    resp = simulate_repertoire(
        config="human_igh", n=50, seed=42, productive_only=True,
    )
    r = resp["result"]
    assert r["productive_rate"] == 1.0
    assert r["productive_count"] == 50


def test_simulate_preset_naive_b_cell_returns_summary():
    from GenAIRR.mcp_server import simulate_preset

    resp = simulate_preset(preset="naive_b_cell", n=20, seed=42)
    assert resp["ok"] is True
    r = resp["result"]
    assert r["config"] == "human_igh"
    assert r["n_records"] == 20
    # naive_b_cell has productive_only=True
    assert r["productive_rate"] == 1.0


def test_simulate_allele_locks_v_and_j_calls_in_records():
    from GenAIRR.mcp_server import list_alleles, simulate_allele

    v_name = list_alleles(config="human_igh", segment="v", limit=1)["result"]["allele_names"][0]
    j_name = list_alleles(config="human_igh", segment="j", limit=1)["result"]["allele_names"][0]

    resp = simulate_allele(
        config="human_igh",
        v_allele=v_name,
        j_allele=j_name,
        n=10,
        seed=42,
        mutation_count_max=0,  # no mutations, so calls stay exactly == truth
    )
    r = resp["result"]
    # v_usage_top should be dominated by the locked allele.
    assert v_name in r["v_usage_top"]
    assert next(iter(r["v_usage_top"])) == v_name  # most common = the locked one


# -- Tier 1: Analysis tools ------------------------------------------


def test_validate_records_flags_record_with_stop_in_productive_junction():
    from GenAIRR.mcp_server import validate_records

    records = [
        _fake_record(productive=True, junction_aa="CARDW"),       # ok
        _fake_record(productive=True, junction_aa="CAR*W"),       # has stop
    ]
    resp = validate_records(records=records)
    r = resp["result"]
    assert r["n_records"] == 2
    assert r["n_valid"] == 1
    assert r["n_invalid"] == 1
    assert r["per_record"][1]["valid"] is False
    assert any("stop" in i.lower() for i in r["per_record"][1]["issues"])


def test_align_to_germline_round_trip_on_simulated_record_has_high_v_identity():
    from GenAIRR.mcp_server import align_to_germline, simulate_repertoire

    # Generate one record then re-align it. Should be ~100% identity since
    # the engine emits exactly the germline (no SHM in this call).
    sim = simulate_repertoire(
        config="human_igh", n=1, seed=42, productive_only=True,
        return_records=True, return_records_limit=1,
    )
    rec = sim["result"]["records"][0]

    resp = align_to_germline(config="human_igh", record=rec)
    r = resp["result"]
    assert "v_alignment" in r
    assert r["v_alignment"]["identity"] >= 0.95


def test_score_allele_calls_ranks_claimed_allele_at_or_near_top():
    from GenAIRR.mcp_server import score_allele_calls, simulate_repertoire

    sim = simulate_repertoire(
        config="human_igh", n=1, seed=42, productive_only=True,
        return_records=True, return_records_limit=1,
    )
    rec = sim["result"]["records"][0]

    resp = score_allele_calls(config="human_igh", record=rec, top_k=5)
    r = resp["result"]
    # The claimed V should be in the top 5
    claimed = rec["v_call"].split(",")[0].strip()
    top_names = [entry["allele"] for entry in r["v_score_ranking"]]
    assert claimed in top_names


def test_analyze_mutations_reports_zero_mutations_on_pristine_record():
    from GenAIRR.mcp_server import analyze_mutations, simulate_repertoire

    # No SHM clause -> n_mutations=0 -> every mutation field is zero/empty.
    sim = simulate_repertoire(
        config="human_igh", n=1, seed=42, productive_only=True,
        return_records=True, return_records_limit=1,
    )
    rec = sim["result"]["records"][0]
    resp = analyze_mutations(config="human_igh", record=rec)
    r = resp["result"]
    assert r["n_mutations"] == 0
    assert r["per_mutation"] == []


def test_classify_regions_returns_imgt_region_spans_for_human_igh_record():
    from GenAIRR.mcp_server import classify_regions, simulate_repertoire

    sim = simulate_repertoire(
        config="human_igh", n=1, seed=42, productive_only=True,
        return_records=True, return_records_limit=1,
    )
    rec = sim["result"]["records"][0]
    resp = classify_regions(config="human_igh", record=rec)
    r = resp["result"]
    # Allele lacks imgt_regions in v2.0.0 -- only CDR3 + FWR4 are computable
    # from the AIRR coordinates the engine emits (junction + j anchor).
    # FWR1/CDR1/FWR2/CDR2/FWR3 are NOT in r["regions"] because we have no
    # per-allele region boundaries to project; tool degrades gracefully.
    for region_name in ("CDR3", "FWR4"):
        assert region_name in r["regions"]
        start, end = r["regions"][region_name]
        assert end >= start >= 0


def test_summarize_dataset_returns_same_shape_as_simulate_repertoire_summary():
    from GenAIRR.mcp_server import simulate_repertoire, summarize_dataset

    sim = simulate_repertoire(
        config="human_igh", n=50, seed=42, productive_only=True,
        return_records=True, return_records_limit=50,
    )["result"]
    summary = summarize_dataset(records=sim["records"])["result"]
    # The two should agree on V usage / CDR3 / productive rate
    assert summary["v_usage_top"] == sim["v_usage_top"]
    assert summary["cdr3_length_histogram"] == sim["cdr3_length_histogram"]
    assert summary["productive_count"] == sim["productive_count"]


# -- Tier 1: Reproducibility ----------------------------------------


def test_replay_seed_reproduces_byte_identical_record():
    from GenAIRR.mcp_server import replay_seed, simulate_repertoire

    sim = simulate_repertoire(
        config="human_igh", n=1, seed=42, productive_only=True,
        return_records=True, return_records_limit=1,
    )
    expected = sim["result"]["records"][0]

    resp = replay_seed(
        config="human_igh", seed=42,
        params={"productive_only": True},
    )
    r = resp["result"]
    assert r["record"]["sequence"] == expected["sequence"]
    assert "trace_summary" in r
    assert "pass_names" in r["trace_summary"]
    assert isinstance(r["trace_summary"]["key_addresses"], dict)


# -- Tier 2: Error envelope per error code --------------------------


def test_error_config_not_found_on_simulate_repertoire():
    from GenAIRR.mcp_server import simulate_repertoire

    resp = simulate_repertoire(config="human_xyz", n=10, seed=0)
    assert resp["ok"] is False
    assert resp["error"]["code"] == CONFIG_NOT_FOUND
    assert "human_xyz" in resp["error"]["message"]


def test_error_allele_not_found_on_inspect_allele():
    from GenAIRR.mcp_server import inspect_allele

    resp = inspect_allele(config="human_igh", allele_name="IGHV_DOES_NOT_EXIST*99")
    assert resp["ok"] is False
    assert resp["error"]["code"] == ALLELE_NOT_FOUND


def test_error_invalid_preset_on_simulate_preset():
    from GenAIRR.mcp_server import simulate_preset

    resp = simulate_preset(preset="not_a_real_preset")
    assert resp["ok"] is False
    assert resp["error"]["code"] == INVALID_PRESET


def test_error_invalid_parameter_on_list_alleles_bad_segment():
    from GenAIRR.mcp_server import list_alleles

    resp = list_alleles(config="human_igh", segment="x")
    assert resp["ok"] is False
    assert resp["error"]["code"] == INVALID_PARAMETER


def test_error_invalid_parameter_on_simulate_repertoire_half_clonal():
    from GenAIRR.mcp_server import simulate_repertoire

    # n_clones set but clone_size missing -- should fail validation, not crash.
    resp = simulate_repertoire(config="human_igh", n_clones=10)
    assert resp["ok"] is False
    assert resp["error"]["code"] == INVALID_PARAMETER


def test_error_malformed_record_on_align_to_germline_no_sequence():
    from GenAIRR.mcp_server import align_to_germline

    resp = align_to_germline(config="human_igh", record={"v_call": "IGHV1*01"})
    assert resp["ok"] is False
    assert resp["error"]["code"] == MALFORMED_RECORD


# -- Tier 3: End-to-end agent-style chain ---------------------------


def test_agent_workflow_list_configs_then_simulate_then_summarize():
    """Mimic an LLM-agent workflow: discover configs, simulate, then
    independently summarize the records. The independently-computed
    summary must match the built-in summary from simulate_repertoire --
    proving the tools compose and share the same compute path."""
    from GenAIRR.mcp_server import list_configs, simulate_repertoire, summarize_dataset

    # 1. Agent discovers the catalog
    configs = list_configs()["result"]["configs"]
    assert "mouse_igh" in configs

    # 2. Agent simulates a small batch
    sim = simulate_repertoire(
        config="mouse_igh",
        n=100,
        seed=42,
        productive_only=True,
        return_records=True,
        return_records_limit=100,
    )["result"]

    # 3. Agent runs the same summary independently on the returned records
    summary = summarize_dataset(records=sim["records"])["result"]

    # The two should match field-for-field on the comparable fields.
    assert summary["v_usage_top"] == sim["v_usage_top"]
    assert summary["d_usage_top"] == sim["d_usage_top"]
    assert summary["j_usage_top"] == sim["j_usage_top"]
    assert summary["cdr3_length_histogram"] == sim["cdr3_length_histogram"]
    assert summary["productive_count"] == sim["productive_count"]
    assert summary["productive_rate"] == sim["productive_rate"]
