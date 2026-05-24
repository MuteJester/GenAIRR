"""PyO3 wrapper smoke tests for the engine.

Exercises every read-only Python wrapper exposed by the Rust kernel
(`Outcome`, `Simulation`, `Region`, `Trace`, `ChoiceRecord`) by
running the built-in `run_smoke_plan(seed)` and reading every
accessor on the returned `Outcome`.

The wrappers are pure read-only views.
"""
from __future__ import annotations

from GenAIRR import _engine as genairr_engine
import pytest


# ──────────────────────────────────────────────────────────────────
# Module surface
# ──────────────────────────────────────────────────────────────────


def test_module_exports_expected_classes():
    for name in (
        "Outcome",
        "Simulation",
        "Region",
        "Trace",
        "ChoiceRecord",
        "CompiledSimulator",
    ):
        assert hasattr(genairr_engine, name), f"missing {name!r} export"


def test_module_exports_smoke_runner():
    assert callable(genairr_engine.run_smoke_plan)


# ──────────────────────────────────────────────────────────────────
# PyOutcome
# ──────────────────────────────────────────────────────────────────


def test_smoke_outcome_has_expected_revision_count():
    # Plan = 8 EchoPass + 8 SampleBasePass = 16 passes → 17 revisions.
    o = genairr_engine.run_smoke_plan(42)
    assert o.revision_count() == 17


def test_smoke_outcome_pass_names_alternate_echo_sample_base():
    o = genairr_engine.run_smoke_plan(42)
    names = o.pass_names()
    assert len(names) == 16
    # Plan alternates echo, sample_base, echo, sample_base, ...
    for i, name in enumerate(names):
        expected = "echo" if i % 2 == 0 else "sample_base"
        assert name == expected, f"position {i}: {name!r} != {expected!r}"


def test_smoke_outcome_revision_indexing():
    o = genairr_engine.run_smoke_plan(42)
    initial = o.revision(0)
    final = o.revision(o.revision_count() - 1)
    assert len(initial) == 0  # initial sim is empty
    assert len(final) == 16  # 8 echo nucs + 8 sample-base nucs


def test_smoke_outcome_revision_out_of_range_raises_index_error():
    o = genairr_engine.run_smoke_plan(42)
    with pytest.raises(IndexError):
        o.revision(o.revision_count())


def test_smoke_outcome_revision_after_returns_first_match():
    # revision_after("echo") should return the first revision after
    # the FIRST echo pass — that is, revision(1).
    o = genairr_engine.run_smoke_plan(42)
    after_first_echo = o.revision_after("echo")
    assert after_first_echo is not None
    assert len(after_first_echo) == 1


def test_smoke_outcome_revision_after_unknown_pass_returns_none():
    o = genairr_engine.run_smoke_plan(42)
    assert o.revision_after("does_not_exist") is None


def test_smoke_outcome_repr_is_informative():
    o = genairr_engine.run_smoke_plan(42)
    r = repr(o)
    assert "Outcome" in r
    assert "revisions=17" in r
    assert "passes=16" in r


# ──────────────────────────────────────────────────────────────────
# PySimulation
# ──────────────────────────────────────────────────────────────────


def test_smoke_final_simulation_pool_is_sixteen_bases():
    o = genairr_engine.run_smoke_plan(42)
    sim = o.final_simulation()
    assert len(sim) == 16
    assert not sim.is_empty()


def test_smoke_final_simulation_bases_are_bytes_and_canonical():
    o = genairr_engine.run_smoke_plan(42)
    sim = o.final_simulation()
    bases = sim.bases()
    assert isinstance(bases, bytes)
    assert len(bases) == 16
    # Every base is one of A, C, G, T (upper).
    assert all(b in b"ACGT" for b in bases)


def test_smoke_final_simulation_germline_position_for_germline_nucs():
    # The Echo passes push germline nucleotides at germline_pos = 0..7.
    # The SampleBase passes push synthetic nucs (no germline pos).
    # Pool indexing alternates echo, sample-base, ... → even indices
    # are germline.
    o = genairr_engine.run_smoke_plan(42)
    sim = o.final_simulation()
    for i in range(16):
        gp = sim.germline_position(i)
        if i % 2 == 0:
            assert gp == i // 2, f"echo at pool {i} expected gp={i // 2}, got {gp}"
        else:
            assert gp is None, f"sample-base at pool {i} should have no germline pos"


def test_smoke_final_simulation_repr_is_informative():
    o = genairr_engine.run_smoke_plan(42)
    sim = o.final_simulation()
    r = repr(sim)
    assert "Simulation" in r
    assert "pool_len=16" in r


def test_smoke_simulation_has_no_assigned_alleles():
    # The smoke plan doesn't include SampleAllelePass, so all
    # assignment slots stay None.
    o = genairr_engine.run_smoke_plan(42)
    sim = o.final_simulation()
    assert sim.v_allele_id() is None
    assert sim.d_allele_id() is None
    assert sim.j_allele_id() is None


# ──────────────────────────────────────────────────────────────────
# PyTrace + PyChoiceRecord
# ──────────────────────────────────────────────────────────────────


def test_smoke_trace_has_one_record_per_sample_base_pass():
    o = genairr_engine.run_smoke_plan(42)
    trace = o.trace()
    assert len(trace) == 8
    assert not trace.is_empty()


def test_smoke_trace_choices_are_choice_records_in_order():
    o = genairr_engine.run_smoke_plan(42)
    trace = o.trace()
    choices = trace.choices()
    assert len(choices) == 8
    for i, c in enumerate(choices):
        assert isinstance(c, genairr_engine.ChoiceRecord)
        assert c.address == f"np.np1.bases[{i}]"


def test_smoke_trace_find_returns_record_at_address():
    o = genairr_engine.run_smoke_plan(42)
    trace = o.trace()
    rec = trace.find("np.np1.bases[3]")
    assert rec is not None
    assert rec.address == "np.np1.bases[3]"
    # A Base value comes back as a 1-byte bytes object.
    assert isinstance(rec.value, bytes)
    assert len(rec.value) == 1
    assert rec.value in (b"A", b"C", b"G", b"T")


def test_smoke_trace_find_returns_none_for_unknown_address():
    o = genairr_engine.run_smoke_plan(42)
    trace = o.trace()
    assert trace.find("not.an.address") is None


def test_smoke_trace_prefix_query_returns_matching_subset():
    o = genairr_engine.run_smoke_plan(42)
    trace = o.trace()
    matches = trace.prefix_query("np.np1.bases")
    assert len(matches) == 8
    assert trace.prefix_count("np.np1.bases") == 8


def test_smoke_trace_prefix_query_with_unmatched_prefix_is_empty():
    o = genairr_engine.run_smoke_plan(42)
    trace = o.trace()
    assert trace.prefix_query("nope.") == []
    assert trace.prefix_count("nope.") == 0


def test_smoke_trace_repr_includes_length():
    o = genairr_engine.run_smoke_plan(42)
    r = repr(o.trace())
    assert "Trace" in r
    assert "len=8" in r


# ──────────────────────────────────────────────────────────────────
# Determinism end-to-end
# ──────────────────────────────────────────────────────────────────


def test_smoke_run_is_deterministic_under_same_seed():
    o1 = genairr_engine.run_smoke_plan(0xC0FFEE)
    o2 = genairr_engine.run_smoke_plan(0xC0FFEE)
    assert o1.final_simulation().bases() == o2.final_simulation().bases()
    # Trace addresses + values match.
    a, b = o1.trace().choices(), o2.trace().choices()
    assert [c.address for c in a] == [c.address for c in b]
    assert [c.value for c in a] == [c.value for c in b]


def test_smoke_run_with_different_seeds_produces_different_traces():
    o1 = genairr_engine.run_smoke_plan(1)
    o2 = genairr_engine.run_smoke_plan(2)
    # 8 independent base draws — overwhelmingly unlikely to match.
    assert o1.final_simulation().bases() != o2.final_simulation().bases()


# ──────────────────────────────────────────────────────────────────
# PyRegion (smoke plan has zero regions, but the type must be usable)
# ──────────────────────────────────────────────────────────────────


def test_region_class_is_exposed():
    # Smoke plan creates no regions, but the class itself must be
    # importable + introspectable. F.2's recombination runner will
    # produce instances we can fully exercise.
    cls = genairr_engine.Region
    assert cls.__name__ == "Region"
    # Frozen pyclass: no constructor exposed.
    with pytest.raises(TypeError):
        cls()


# ──────────────────────────────────────────────────────────────────
# F.2 — full VJ recombination through the bridge
# ──────────────────────────────────────────────────────────────────


def test_vj_recombination_runs_five_passes():
    o = genairr_engine.run_smoke_vj_recombination(42)
    assert o.revision_count() == 6  # initial + 5 passes
    assert o.pass_names() == [
        "sample_allele.v",
        "sample_allele.j",
        "assemble.v",
        "generate_np.np1",
        "assemble.j",
    ]


def test_vj_recombination_produces_three_regions_in_order():
    o = genairr_engine.run_smoke_vj_recombination(42)
    sim = o.final_simulation()
    regions = sim.regions()
    assert [r.segment for r in regions] == ["V", "NP1", "J"]


def test_vj_recombination_v_and_j_regions_have_canonical_lengths():
    # V is 9bp, J is 6bp in the smoke refdata. NP1 length varies
    # by seed.
    o = genairr_engine.run_smoke_vj_recombination(42)
    regions = o.final_simulation().regions()
    v, np1, j = regions
    assert len(v) == 9
    assert len(j) == 6
    assert 0 <= len(np1) <= 6  # NP1 length distribution is 0..=6


def test_vj_recombination_pool_length_matches_region_chain():
    # pool_len = V_len + NP1_len + J_len.
    o = genairr_engine.run_smoke_vj_recombination(42)
    sim = o.final_simulation()
    total = sum(len(r) for r in sim.regions())
    assert len(sim) == total


def test_vj_recombination_assigns_v_and_j_allele_ids():
    o = genairr_engine.run_smoke_vj_recombination(42)
    sim = o.final_simulation()
    # Single-allele pool → always id 0.
    assert sim.v_allele_id() == 0
    assert sim.j_allele_id() == 0
    assert sim.d_allele_id() is None  # VJ chain — no D


def test_vj_recombination_v_region_carries_germline_provenance():
    # V is assembled from `AAACCCGGG`; every V-region nucleotide
    # should have a germline_position equal to its position in the
    # source allele.
    o = genairr_engine.run_smoke_vj_recombination(42)
    sim = o.final_simulation()
    v = sim.regions()[0]
    for i in range(v.start, v.end):
        gp = sim.germline_position(i)
        assert gp == i - v.start, (
            f"pool {i} has germline_pos {gp}, expected {i - v.start}"
        )


def test_vj_recombination_np1_region_has_no_germline_provenance():
    # NP1 is synthetic — every nucleotide must report germline_position = None.
    o = genairr_engine.run_smoke_vj_recombination(42)
    sim = o.final_simulation()
    np1 = sim.regions()[1]
    for i in range(np1.start, np1.end):
        assert sim.germline_position(i) is None


def test_vj_recombination_v_region_codon_rail_is_KPG():
    # V allele AAACCCGGG → codons AAA, CCC, GGG → K, P, G.
    o = genairr_engine.run_smoke_vj_recombination(42)
    v = o.final_simulation().regions()[0]
    assert v.amino_acids() == b"KPG"


def test_vj_recombination_trace_records_canonical_addresses():
    o = genairr_engine.run_smoke_vj_recombination(42)
    trace = o.trace()
    # Always present.
    for addr in ("sample_allele.v", "sample_allele.j", "np.np1.length"):
        assert trace.find(addr) is not None, f"missing {addr}"


def test_vj_recombination_np1_length_in_trace_matches_region_length():
    o = genairr_engine.run_smoke_vj_recombination(42)
    sim = o.final_simulation()
    np1 = sim.regions()[1]
    rec = o.trace().find("np.np1.length")
    assert rec is not None
    assert rec.value == len(np1)


def test_vj_recombination_np_base_records_match_pool_bases():
    # For each np.np1.bases[i] in the trace, the i-th NP1 nucleotide
    # in the pool must carry that base.
    o = genairr_engine.run_smoke_vj_recombination(42)
    sim = o.final_simulation()
    np1 = sim.regions()[1]
    np_records = o.trace().prefix_query("np.np1.bases")
    assert len(np_records) == len(np1)
    pool_bases = sim.bases()
    for i, rec in enumerate(np_records):
        assert rec.address == f"np.np1.bases[{i}]"
        assert rec.value == bytes([pool_bases[np1.start + i]])


def test_vj_recombination_revision_after_each_pass_is_visible():
    o = genairr_engine.run_smoke_vj_recombination(42)
    after_v = o.revision_after("assemble.v")
    after_np1 = o.revision_after("generate_np.np1")
    after_j = o.revision_after("assemble.j")

    # After V assembly: 9 nucs, 1 region.
    assert len(after_v) == 9
    assert after_v.region_count() == 1
    # After NP1: 9 + np1_len nucs, 2 regions.
    assert after_np1.region_count() == 2
    # After J assembly: pool grew by 6, 3 regions total.
    assert after_j.region_count() == 3
    assert len(after_j) == len(after_np1) + 6


def test_vj_recombination_is_deterministic_under_same_seed():
    o1 = genairr_engine.run_smoke_vj_recombination(0xCAFE)
    o2 = genairr_engine.run_smoke_vj_recombination(0xCAFE)
    assert o1.final_simulation().bases() == o2.final_simulation().bases()
    # Per-record equality on the trace.
    a, b = o1.trace().choices(), o2.trace().choices()
    assert [c.address for c in a] == [c.address for c in b]
    assert [c.value for c in a] == [c.value for c in b]


def test_vj_recombination_region_frame_phases_chain_correctly():
    # V starts at frame_phase 0. NP1 starts at (9 % 3) = 0. J starts
    # at (9 + np1_len) % 3.
    o = genairr_engine.run_smoke_vj_recombination(42)
    regions = o.final_simulation().regions()
    v, np1, j = regions

    assert v.frame_phase == 0
    assert np1.frame_phase == 0  # 9 % 3 == 0
    expected_j_phase = (9 + len(np1)) % 3
    assert j.frame_phase == expected_j_phase


# ──────────────────────────────────────────────────────────────────
# F.3 — RefDataConfig + Allele Python bridge
# ──────────────────────────────────────────────────────────────────


def test_refdata_config_vj_factory_starts_empty():
    cfg = genairr_engine.RefDataConfig.vj()
    assert cfg.chain_type == "vj"
    assert cfg.has_d() is False
    assert cfg.v_pool_size() == 0
    assert cfg.d_pool_size() == 0
    assert cfg.j_pool_size() == 0


def test_refdata_config_vdj_factory_starts_empty():
    cfg = genairr_engine.RefDataConfig.vdj()
    assert cfg.chain_type == "vdj"
    assert cfg.has_d() is True
    assert cfg.v_pool_size() == 0
    assert cfg.d_pool_size() == 0
    assert cfg.j_pool_size() == 0


def test_refdata_config_constructor_accepts_chain_type_string():
    a = genairr_engine.RefDataConfig("vj")
    b = genairr_engine.RefDataConfig("VJ")  # case-insensitive
    c = genairr_engine.RefDataConfig("Vdj")
    assert a.chain_type == "vj"
    assert b.chain_type == "vj"
    assert c.chain_type == "vdj"


def test_refdata_config_constructor_rejects_unknown_chain_type():
    with pytest.raises(ValueError, match="unknown chain type"):
        genairr_engine.RefDataConfig("igh")


def test_refdata_config_repr_includes_pool_sizes():
    cfg = genairr_engine.RefDataConfig.vdj()
    cfg.add_v_allele("v1*01", "v1", b"AAA")
    cfg.add_v_allele("v2*01", "v2", b"CCC")
    cfg.add_d_allele("d1*01", "d1", b"GGG")
    r = repr(cfg)
    assert "vdj" in r
    assert "V=2" in r
    assert "D=1" in r
    assert "J=0" in r


def test_refdata_config_add_v_allele_returns_id():
    cfg = genairr_engine.RefDataConfig.vj()
    a_id = cfg.add_v_allele("v1*01", "v1", b"AAA", anchor=0)
    b_id = cfg.add_v_allele("v2*01", "v2", b"CCCGGG", anchor=3)
    assert a_id == 0
    assert b_id == 1
    assert cfg.v_pool_size() == 2


def test_refdata_config_add_j_allele_returns_id():
    cfg = genairr_engine.RefDataConfig.vj()
    j_id = cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    assert j_id == 0
    assert cfg.j_pool_size() == 1


def test_refdata_config_add_d_allele_works_on_vdj():
    cfg = genairr_engine.RefDataConfig.vdj()
    d_id = cfg.add_d_allele("d1*01", "d1", b"GGGAAA")
    assert d_id == 0
    assert cfg.d_pool_size() == 1


def test_refdata_config_add_d_allele_rejects_vj_chain():
    cfg = genairr_engine.RefDataConfig.vj()
    with pytest.raises(ValueError, match="cannot add D allele"):
        cfg.add_d_allele("d1*01", "d1", b"GGG")


def test_refdata_config_add_allele_anchor_defaults_to_none():
    cfg = genairr_engine.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAA")  # no anchor kwarg
    a = cfg.v_allele(0)
    assert a.anchor is None


def test_refdata_config_v_allele_round_trips_metadata():
    cfg = genairr_engine.RefDataConfig.vj()
    cfg.add_v_allele("IGHV1-2*01", "IGHV1-2", b"ACGTACGT", anchor=3)
    a = cfg.v_allele(0)
    assert a.name == "IGHV1-2*01"
    assert a.gene == "IGHV1-2"
    assert a.seq() == b"ACGTACGT"
    assert len(a) == 8
    assert a.segment == "V"
    assert a.anchor == 3


def test_refdata_config_v_allele_out_of_range_raises_index_error():
    cfg = genairr_engine.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAA")
    with pytest.raises(IndexError, match="V allele id 5"):
        cfg.v_allele(5)


def test_refdata_config_d_allele_out_of_range_raises_index_error():
    cfg = genairr_engine.RefDataConfig.vdj()
    with pytest.raises(IndexError, match="D allele id 0"):
        cfg.d_allele(0)


def test_allele_repr_includes_segment_and_anchor():
    cfg = genairr_engine.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
    r = repr(cfg.v_allele(0))
    assert "Allele" in r
    assert "V" in r
    assert "v1*01" in r
    assert "9bp" in r
    assert "anchor=6" in r


def test_allele_repr_for_anchorless():
    cfg = genairr_engine.RefDataConfig.vdj()
    cfg.add_d_allele("d1*01", "d1", b"GGGAAA")
    r = repr(cfg.d_allele(0))
    assert "anchor=None" in r


# ── Recombination via Python-supplied refdata ─────────────────────


def _build_smoke_vj_refdata():
    cfg = genairr_engine.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    return cfg


def test_run_vj_recombination_with_supplied_refdata_matches_smoke_runner():
    # The Python-built fixture is byte-for-byte identical to the
    # baked-in fixture, so the resulting traces and pools must
    # match the existing smoke runner under the same seed.
    cfg = _build_smoke_vj_refdata()
    seed = 0xCAFE
    o_supplied = genairr_engine.run_vj_recombination(cfg, seed)
    o_baked = genairr_engine.run_smoke_vj_recombination(seed)

    assert o_supplied.final_simulation().bases() == o_baked.final_simulation().bases()
    a, b = o_supplied.trace().choices(), o_baked.trace().choices()
    assert [c.address for c in a] == [c.address for c in b]
    assert [c.value for c in a] == [c.value for c in b]


def test_run_vj_recombination_uses_supplied_v_allele_sequence():
    # A different V allele → different V region bases in the output.
    cfg = genairr_engine.RefDataConfig.vj()
    cfg.add_v_allele("v_custom*01", "v_custom", b"GGGAAACCC", anchor=6)
    cfg.add_j_allele("j_custom*01", "j_custom", b"AAATTT", anchor=0)

    o = genairr_engine.run_vj_recombination(cfg, seed=42)
    sim = o.final_simulation()
    bases = sim.bases()
    # First 9 bases come from V.
    assert bases[:9] == b"GGGAAACCC"
    # Last 6 bases come from J.
    assert bases[-6:] == b"AAATTT"


def test_run_vj_recombination_samples_different_v_alleles_across_seeds():
    cfg = genairr_engine.RefDataConfig.vj()
    cfg.add_v_allele("vA*01", "vA", b"AAAAAA", anchor=0)
    cfg.add_v_allele("vB*01", "vB", b"CCCCCC", anchor=0)
    cfg.add_v_allele("vC*01", "vC", b"GGGGGG", anchor=0)
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)

    seen_ids = set()
    for seed in range(50):
        o = genairr_engine.run_vj_recombination(cfg, seed)
        seen_ids.add(o.final_simulation().v_allele_id())
        if len(seen_ids) == 3:
            break

    # 50 seeds × 3 alleles → overwhelmingly likely to see all three.
    assert seen_ids == {0, 1, 2}, f"sampled only {seen_ids}"


def test_run_vj_recombination_rejects_vdj_refdata():
    cfg = genairr_engine.RefDataConfig.vdj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    with pytest.raises(ValueError, match="requires a VJ-chain"):
        genairr_engine.run_vj_recombination(cfg, seed=0)


def test_run_vj_recombination_rejects_empty_v_pool():
    cfg = genairr_engine.RefDataConfig.vj()
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    with pytest.raises(ValueError, match="V pool is empty"):
        genairr_engine.run_vj_recombination(cfg, seed=0)


def test_run_vj_recombination_rejects_empty_j_pool():
    cfg = genairr_engine.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
    with pytest.raises(ValueError, match="J pool is empty"):
        genairr_engine.run_vj_recombination(cfg, seed=0)


def test_run_vj_recombination_rejects_invalid_np1_max_length():
    cfg = _build_smoke_vj_refdata()
    with pytest.raises(ValueError, match="np1_max_length"):
        genairr_engine.run_vj_recombination(cfg, seed=0, np1_max_length=0)


def test_run_vj_recombination_respects_np1_max_length():
    # max_length=4 → NP1 length distribution {0, 1, 2, 3}. Sample many
    # seeds and verify every NP1 length stays in that range.
    cfg = _build_smoke_vj_refdata()
    for seed in range(50):
        o = genairr_engine.run_vj_recombination(cfg, seed, np1_max_length=4)
        np1 = o.final_simulation().regions()[1]
        assert 0 <= len(np1) <= 3


# ──────────────────────────────────────────────────────────────────
# F.4 — PyPassPlan builder + generic run()
# ──────────────────────────────────────────────────────────────────


def _build_smoke_vj_plan(cfg):
    plan = genairr_engine.PassPlan()
    plan.push_sample_allele("V", cfg)
    plan.push_sample_allele("J", cfg)
    plan.push_assemble("V")
    plan.push_generate_np("NP1", [(i, 1.0) for i in range(7)])
    plan.push_assemble("J")
    return plan


def test_pass_plan_starts_empty():
    plan = genairr_engine.PassPlan()
    assert len(plan) == 0
    assert plan.is_empty() is True
    assert "len=0" in repr(plan)


def test_pass_plan_grows_with_each_push():
    cfg = _build_smoke_vj_refdata()
    plan = genairr_engine.PassPlan()
    plan.push_sample_allele("V", cfg)
    assert len(plan) == 1
    plan.push_sample_allele("J", cfg)
    plan.push_assemble("V")
    plan.push_generate_np("NP1", [(0, 1.0), (3, 1.0), (6, 1.0)])
    plan.push_assemble("J")
    assert len(plan) == 5
    assert plan.is_empty() is False


def test_pass_plan_run_matches_run_vj_recombination_under_same_seed():
    # Build a plan from Python that's structurally identical to
    # run_vj_recombination's hardcoded recipe. Same seed → same outcome.
    cfg = _build_smoke_vj_refdata()
    seed = 0xC0FFEE

    plan = _build_smoke_vj_plan(cfg)

    o_built = genairr_engine.run(plan, seed=seed, refdata=cfg)
    o_recipe = genairr_engine.run_vj_recombination(cfg, seed)
    assert o_built.final_simulation().bases() == o_recipe.final_simulation().bases()
    assert [c.address for c in o_built.trace().choices()] == [
        c.address for c in o_recipe.trace().choices()
    ]
    assert [c.value for c in o_built.trace().choices()] == [
        c.value for c in o_recipe.trace().choices()
    ]


def test_pass_plan_exposes_last_pushed_node_and_explicit_edges():
    # Build a plan with two independent allele samples (V and D) and
    # force D to run before V via an explicit edge. The compile step
    # should topo-sort accordingly.
    cfg = _build_smoke_vj_refdata()
    plan = genairr_engine.PassPlan()

    plan.push_sample_allele("V", cfg)
    v_idx = plan.last_pushed_node()
    plan.push_sample_allele("J", cfg)
    j_idx = plan.last_pushed_node()
    assert v_idx == 0
    assert j_idx == 1

    # `.after(j_idx, v_idx)` means "j runs after v". This is already
    # the natural order; the API should accept it without complaint.
    plan.after(j_idx, v_idx)

    # `.add_edge` with an out-of-range index errors at the Python boundary.
    with pytest.raises(ValueError):
        plan.add_edge(99, 0)
    with pytest.raises(ValueError):
        plan.add_edge(0, 0)  # self-edge


def test_pass_plan_explicit_edge_reorders_independent_passes():
    # Push V then D, then declare D-before-V. The compiled report's
    # pass_names should reflect the topo-sorted order with D first.
    cfg = _build_smoke_vj_refdata()
    plan = genairr_engine.PassPlan()

    plan.push_sample_allele("V", cfg)
    v_idx = plan.last_pushed_node()
    plan.push_sample_allele("J", cfg)
    j_idx = plan.last_pushed_node()
    plan.push_assemble("V")
    plan.push_generate_np("NP1", [(0, 1.0), (3, 1.0)])
    plan.push_assemble("J")

    # Force J's allele draw to come before V's.
    plan.before(j_idx, v_idx)

    compiled = plan.compile(refdata=cfg)
    names = compiled.pass_names()
    assert names.index("sample_allele.j") < names.index("sample_allele.v"), (
        f"explicit edge `j before v` must reorder allele draws: {names}"
    )


def test_pass_plan_run_without_refdata_works_for_pure_np_plan():
    # NP-only plan: no allele sampling, no assembly. Doesn't need refdata.
    plan = genairr_engine.PassPlan()
    plan.push_generate_np("NP1", [(5, 1.0)])
    o = genairr_engine.run(plan, seed=42)
    assert o.revision_count() == 2
    assert o.pass_names() == ["generate_np.np1"]
    sim = o.final_simulation()
    assert len(sim) == 5
    assert sim.region_count() == 1
    assert sim.regions()[0].segment == "NP1"


def test_pass_plan_push_trim_reduces_assembled_v_length():
    # V is 9bp (AAACCCGGG). Trim 5'=2, 3'=3 → 4 bases assembled.
    cfg = _build_smoke_vj_refdata()
    plan = genairr_engine.PassPlan()
    plan.push_sample_allele("V", cfg)
    plan.push_trim("V", "5", [(2, 1.0)])
    plan.push_trim("V", "3", [(3, 1.0)])
    plan.push_assemble("V")

    o = genairr_engine.run(plan, seed=0, refdata=cfg)
    sim = o.final_simulation()
    assert len(sim) == 4
    assert sim.bases() == b"ACCC"  # AAACCCGGG[2..6)


def test_pass_plan_push_trim_records_trim_address():
    cfg = _build_smoke_vj_refdata()
    plan = genairr_engine.PassPlan()
    plan.push_sample_allele("V", cfg)
    plan.push_trim("V", "3", [(2, 1.0)])

    o = genairr_engine.run(plan, seed=0, refdata=cfg)
    rec = o.trace().find("trim.v_3")
    assert rec is not None
    assert rec.value == 2


def test_pass_plan_push_sample_allele_supports_all_three_segments():
    cfg = genairr_engine.RefDataConfig.vdj()
    cfg.add_v_allele("v*01", "v", b"AAA", anchor=0)
    cfg.add_d_allele("d*01", "d", b"CCC")
    cfg.add_j_allele("j*01", "j", b"GGG", anchor=0)

    plan = genairr_engine.PassPlan()
    plan.push_sample_allele("V", cfg)
    plan.push_sample_allele("D", cfg)
    plan.push_sample_allele("J", cfg)
    o = genairr_engine.run(plan, seed=0, refdata=cfg)

    assert o.pass_names() == [
        "sample_allele.v",
        "sample_allele.d",
        "sample_allele.j",
    ]


def test_pass_plan_push_generate_np_supports_np1_and_np2():
    plan = genairr_engine.PassPlan()
    plan.push_generate_np("NP1", [(2, 1.0)])
    plan.push_generate_np("NP2", [(3, 1.0)])
    o = genairr_engine.run(plan, seed=0)
    assert o.pass_names() == ["generate_np.np1", "generate_np.np2"]


def test_pass_plan_push_segment_is_case_insensitive():
    cfg = _build_smoke_vj_refdata()
    plan = genairr_engine.PassPlan()
    plan.push_sample_allele("v", cfg)  # lowercase
    plan.push_assemble("V")  # uppercase
    plan.push_generate_np("np1", [(0, 1.0)])
    assert len(plan) == 3


def test_pass_plan_push_sample_allele_rejects_np_segment():
    cfg = _build_smoke_vj_refdata()
    plan = genairr_engine.PassPlan()
    with pytest.raises(ValueError, match="V', 'D', or 'J'"):
        plan.push_sample_allele("NP1", cfg)


def test_pass_plan_push_sample_allele_rejects_empty_pool():
    # Build a refdata that has J alleles but no V alleles, then try
    # to push_sample_allele("V").
    cfg = genairr_engine.RefDataConfig.vj()
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    plan = genairr_engine.PassPlan()
    with pytest.raises(ValueError, match="V pool is empty"):
        plan.push_sample_allele("V", cfg)


def test_pass_plan_push_assemble_rejects_np_segment():
    plan = genairr_engine.PassPlan()
    with pytest.raises(ValueError, match="V', 'D', or 'J'"):
        plan.push_assemble("NP1")


def test_pass_plan_push_generate_np_rejects_v_segment():
    plan = genairr_engine.PassPlan()
    with pytest.raises(ValueError, match="NP1' or 'NP2'"):
        plan.push_generate_np("V", [(1, 1.0)])


def test_pass_plan_push_generate_np_rejects_empty_pairs():
    plan = genairr_engine.PassPlan()
    with pytest.raises(ValueError, match="length_pairs"):
        plan.push_generate_np("NP1", [])


def test_pass_plan_push_trim_rejects_bad_end():
    plan = genairr_engine.PassPlan()
    with pytest.raises(ValueError, match="end must be"):
        plan.push_trim("V", "7", [(1, 1.0)])


def test_pass_plan_push_trim_accepts_either_string_form_for_end():
    plan = genairr_engine.PassPlan()
    plan.push_trim("V", "5", [(1, 1.0)])
    plan.push_trim("V", "3", [(1, 1.0)])
    plan.push_trim("V", "5'", [(1, 1.0)])
    plan.push_trim("V", "three", [(1, 1.0)])
    assert len(plan) == 4


def test_pass_plan_push_trim_rejects_empty_pairs():
    plan = genairr_engine.PassPlan()
    with pytest.raises(ValueError, match="length_pairs"):
        plan.push_trim("V", "3", [])


def test_pass_plan_run_is_deterministic_under_same_seed():
    cfg = _build_smoke_vj_refdata()
    plan = _build_smoke_vj_plan(cfg)

    o1 = genairr_engine.run(plan, seed=0xCAFE, refdata=cfg)
    o2 = genairr_engine.run(plan, seed=0xCAFE, refdata=cfg)
    assert o1.final_simulation().bases() == o2.final_simulation().bases()


def test_pass_plan_can_be_re_run_with_different_seeds():
    # The plan is borrowed by the runner, not consumed — so we can
    # call run() multiple times against the same plan.
    cfg = _build_smoke_vj_refdata()
    plan = _build_smoke_vj_plan(cfg)

    runs = [genairr_engine.run(plan, seed=s, refdata=cfg) for s in range(5)]
    bases = [r.final_simulation().bases() for r in runs]
    # Not all 5 should be identical (NP1 varies with seed).
    assert len(set(bases)) > 1


# ──────────────────────────────────────────────────────────────────
# F.8 — owning CompiledSimulator
# ──────────────────────────────────────────────────────────────────


def test_pass_plan_compile_returns_owning_compiled_simulator():
    cfg = _build_smoke_vj_refdata()
    plan = _build_smoke_vj_plan(cfg)

    compiled = plan.compile(refdata=cfg)

    assert isinstance(compiled, genairr_engine.CompiledSimulator)
    assert compiled.policy() == "permissive"
    assert compiled.pass_names() == [
        "sample_allele.v",
        "sample_allele.j",
        "assemble.v",
        "generate_np.np1",
        "assemble.j",
    ]
    assert "consumed=true" in repr(plan)
    with pytest.raises(ValueError, match="consumed"):
        len(plan)
    with pytest.raises(ValueError, match="consumed"):
        plan.push_generate_np("NP1", [(0, 1.0)])


def test_compiled_simulator_runs_repeatedly_without_recompiling():
    cfg = _build_smoke_vj_refdata()
    compiled = _build_smoke_vj_plan(cfg).compile(refdata=cfg)

    a = compiled.run(0xCAFE)
    b = compiled.run(0xCAFE)

    assert a.final_simulation().bases() == b.final_simulation().bases()
    assert a.pass_names() == b.pass_names() == compiled.pass_names()


def test_compiled_simulator_run_batch_uses_seed_stitching():
    cfg = _build_smoke_vj_refdata()
    compiled = _build_smoke_vj_plan(cfg).compile(refdata=cfg)

    batch = compiled.run_batch(3, 100)
    expected = [compiled.run(seed) for seed in range(100, 103)]

    assert len(batch) == 3
    assert [o.final_simulation().bases() for o in batch] == [
        o.final_simulation().bases() for o in expected
    ]


def test_pass_plan_compile_failure_does_not_consume_builder():
    plan = genairr_engine.PassPlan()
    plan.push_assemble("V")

    with pytest.raises(ValueError, match="requires reference data"):
        plan.compile()

    assert len(plan) == 1
    assert "len=1" in repr(plan)


def test_pass_plan_productive_compile_rejects_anchorless_v_support():
    cfg = genairr_engine.RefDataConfig.vj()
    cfg.add_v_allele("v_anchorless*01", "v_anchorless", b"AAACCCGGG")
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    plan = _build_smoke_vj_plan(cfg)

    with pytest.raises(ValueError, match="V sample support has no alleles"):
        plan.compile(refdata=cfg, respect=genairr_engine.productive())

    assert len(plan) == 5
    assert "len=5" in repr(plan)


def test_pass_plan_productive_compile_rejects_no_in_frame_np_mass():
    cfg = _build_smoke_vj_refdata()
    plan = genairr_engine.PassPlan()
    plan.push_sample_allele("V", cfg)
    plan.push_sample_allele("J", cfg)
    plan.push_assemble("V")
    plan.push_generate_np("NP1", [(1, 1.0), (2, 1.0)])
    plan.push_assemble("J")

    with pytest.raises(ValueError, match="NP1 length support has no in-frame mass"):
        plan.compile(refdata=cfg, respect=genairr_engine.productive())

    assert len(plan) == 5


def test_compiled_simulator_captures_contracts_and_strict_policy():
    cfg = _build_smoke_vj_refdata()
    compiled = _build_smoke_vj_plan(cfg).compile(
        refdata=cfg,
        respect=genairr_engine.productive(),
        strict=True,
    )

    assert compiled.policy() == "strict"
    assert compiled.active_contracts() == [
        "productive_junction_frame",
        "no_stop_codon_in_junction",
        "anchor_preserved.v",
        "anchor_preserved.j",
    ]
    out = compiled.run(0)
    assert out.final_simulation().regions()[1].frame_phase == 0
