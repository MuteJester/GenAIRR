"""Public continuation oracles for certified duplex landmarks and coding phase."""


import pytest

from GenAIRR import _engine as ge
from GenAIRR.experimental import project_duplex_simulation
from tests.support.fixtures.duplex_observation import _compile, _fixture

pytestmark = pytest.mark.behavior


COMPLEMENT = str.maketrans("ACGT", "TGCA")
STOPS = {"TAA", "TAG", "TGA"}








def _project(outcome, refs):
    return project_duplex_simulation(outcome.final_simulation(), refs)


def _independent_productive_checks(projection, *, phase):
    """Use explicit positions, literal residue codons and the three stop codons."""
    rec = projection["record"]
    sequence = rec["sequence"].upper()
    if rec["rev_comp"]:
        sequence = sequence.translate(COMPLEMENT)[::-1]
    start, end = rec["junction_start"], rec["junction_end"]
    assert sequence[start:start + 3] in {"TGT", "TGC"}
    assert sequence[end - 3:end] in {"TTT", "TTC", "TGG"}
    assert start % 3 == phase
    assert (end - start) % 3 == 0
    assert rec["junction"].upper() == sequence[start:end]
    coding_start = rec["v_sequence_start"]
    coding_start += (phase - coding_start) % 3
    coding_end = rec["j_sequence_end"]
    assert all(sequence[p:p + 3] not in STOPS for p in range(coding_start, coding_end - 2, 3))
    assert rec["productive"] is True
    assert rec["vj_in_frame"] is True
    assert rec["stop_codon"] is False
    assert rec["productivity_policy"] == "duplex_observed_vdj_v1"
    assert rec["productivity_scope"] == "observed_coding_span"


def _find_outcome(compiled, seed, predicate):
    # Locate a specific ordinary sampled edit for a literal geometry check.
    # This is test-fixture selection, never a production rejection wrapper.
    for rng_seed in range(256):
        outcome = compiled.run_from(seed, rng_seed)
        if predicate(outcome):
            return outcome
    pytest.fail("the deterministic fixture seeds did not contain the requested edit")


@pytest.mark.parametrize("order", [("vd", "dj"), ("dj", "vd")])
@pytest.mark.parametrize("inverted", [False, True])
def test_restored_anchor_remains_productive_through_public_seed_continuation(order, inverted):
    refs, repair, seed = _fixture(order=order, inverted=inverted, d="ACG")
    outcome = _compile(ge.PassPlan(), refs, productive=True).run_from(seed, 3)
    p = _project(outcome, refs)
    rec = p["record"]
    assert rec["sequence"] == ("AAATGTCGTTTTAAA" if inverted else "AAATGTACGTTTAAA")
    assert rec["junction_aa"] == ("CRF" if inverted else "CTF")
    assert (rec["junction_start"], rec["junction_end"]) == (3, 12)
    assert rec["d_inverted"] is inverted
    assert [point["reference_position"] for point in p["paths"]["d"]["body"]] == (
        [2, 1, 0] if inverted else [0, 1, 2]
    )
    founder = rec["duplex_repair"]["product"]
    assert founder == repair.run()
    assert founder["columns"][5]["top"]["kind"] == "templated_synthesis"
    assert founder["columns"][5]["top"]["original"] is None
    assert founder["columns"][5]["top"]["template_id"] == 11  # Original V-bottom position 5.
    assert rec["p_v_3_length"] is None
    _independent_productive_checks(p, phase=0)


@pytest.mark.parametrize("strict", [False, True])
def test_productive_mutation_loss_and_reverse_output_keep_founder_and_event_burden(strict):
    refs, repair, seed = _fixture()
    plan = ge.PassPlan()
    plan.push_mutate_uniform([(3, 1.0)])
    plan.push_corrupt_5prime_loss([(1, 1.0)])
    plan.push_corrupt_rev_comp(1.0)
    compiled = _compile(plan, refs, productive=True, strict=strict)
    founder = repair.run()
    for rng_seed in range(12):
        outcome = compiled.run_from(seed, rng_seed)
        p = _project(outcome, refs)
        rec = p["record"]
        coding_sequence = outcome.final_simulation().bases().decode().upper()
        assert rec["sequence"].upper() == coding_sequence.translate(COMPLEMENT)[::-1]
        assert rec["rev_comp"] is True
        assert rec["v_cigar"] == "5M"
        assert rec["v_trim_5"] == 0 and rec["v_germline_start"] == 1
        assert rec["n_mutations"] == 3
        assert rec["n_v_mutations"] + rec["n_d_mutations"] + rec["n_j_mutations"] == 3
        assert rec["mutation_rate"] == pytest.approx(3 / 15)
        assert rec["end_loss_5_length"] == 1
        assert rec["duplex_repair"]["product"] == founder
        observation = rec["duplex_repair"]["observation"]
        assert [removed["cause"] for removed in observation["removed"]] == ["observation_loss"]
        assert len(observation["links"]) == 14
        _independent_productive_checks(p, phase=2)
    # Every run branches from the unchanged persistent founder.
    assert seed.bases() == b"AAATGTCCCTTTAAA"


@pytest.mark.parametrize("strict", [False, True])
def test_exact_productive_indels_preserve_restored_landmarks_and_report_real_events(strict):
    refs, repair, seed = _fixture()
    plan = ge.PassPlan()
    plan.push_corrupt_indel([(2, 1.0)], rejection_attempts=0, exact_node_budget=100_000)
    compiled = _compile(plan, refs, productive=True, strict=strict)
    for rng_seed in (1, 7, 17):
        outcome = compiled.run_from(seed, rng_seed)
        p = _project(outcome, refs)
        rec = p["record"]
        assert rec["n_indels"] == 2
        assert rec["duplex_repair"]["product"] == repair.run()
        assert outcome.trace().diagnostics()[0]["applied_count"] == 2
        _independent_productive_checks(p, phase=0)


@pytest.mark.parametrize("end,site,expected_aa", [("five", 0, "NVPL"), ("three", 14, "KCPF")])
def test_structural_terminal_deletion_and_read_clipping_have_distinct_phase_truth(end, site, expected_aa):
    refs, _, seed = _fixture()
    deletion = ge.PassPlan()
    deletion.push_corrupt_indel([(1, 1.0)], insertion_prob=0.0)
    compiled = _compile(deletion, refs)
    deleted = _find_outcome(
        compiled, seed, lambda out: out.trace().find("corrupt.indel.site[0]").value == site,
    )
    loss = ge.PassPlan()
    if end == "five":
        loss.push_corrupt_5prime_loss([(1, 1.0)])
    else:
        loss.push_corrupt_3prime_loss([(1, 1.0)])
    clipped = _compile(loss, refs, productive=True).run_from(seed, 0)
    deleted_rec = _project(deleted, refs)["record"]
    clipped_rec = _project(clipped, refs)["record"]
    assert deleted_rec["sequence"] == clipped_rec["sequence"]
    assert deleted_rec["productive"] is False
    assert clipped_rec["productive"] is True
    assert deleted_rec["sequence_aa"] == expected_aa
    assert deleted_rec["stop_codon"] is False and clipped_rec["stop_codon"] is False
    assert deleted_rec["duplex_repair"]["observation"]["removed"][0]["cause"] == "structural_deletion"
    assert clipped_rec["duplex_repair"]["observation"]["removed"][0]["cause"] == "observation_loss"
    if end == "five":
        assert deleted_rec["v_cigar"] == "1D5M"
        assert clipped_rec["v_cigar"] == "5M"
        assert deleted_rec["vj_in_frame"] is False and clipped_rec["vj_in_frame"] is True
        assert clipped_rec["sequence_aa"] == "CPFK"
    else:
        assert deleted_rec["j_cigar"] == "5M1D"
        assert clipped_rec["j_cigar"] == "5M"
        assert deleted_rec["vj_in_frame"] is True and clipped_rec["vj_in_frame"] is True


def test_clipping_through_part_of_a_restored_anchor_is_unknown_without_losing_founder_evidence():
    refs, repair, seed = _fixture()
    loss = ge.PassPlan()
    loss.push_corrupt_5prime_loss([(4, 1.0)])
    p = _project(_compile(loss, refs).run_from(seed, 0), refs)
    rec = p["record"]
    assert rec["sequence"] == "GTCCCTTTAAA"
    assert rec["productive"] is None and rec["junction"] is None and rec["junction_aa"] is None
    assert "anchor_not_observed" in rec["productivity_reasons"]
    assert rec["v_trim_5"] == 0 and rec["v_germline_start"] == 4
    assert rec["duplex_repair"]["product"] == repair.run()


@pytest.mark.parametrize("strict", [False, True])
@pytest.mark.parametrize("failed_predicate", ["frameshift", "unknown_anchor"])
def test_final_productive_assertion_rejects_failed_or_unknown_duplex_continuations(strict, failed_predicate):
    refs, _, seed = _fixture()
    corrupt = ge.PassPlan()
    if failed_predicate == "frameshift":
        corrupt.push_corrupt_indel([(1, 1.0)], insertion_prob=0.0)
        predicate = lambda out: out.trace().find("corrupt.indel.site[0]").value == 0
    else:
        corrupt.push_corrupt_ns([(1, 1.0)])
        predicate = lambda out: out.final_simulation().bases()[5] == ord("N")
    changed = _find_outcome(_compile(corrupt, refs), seed, predicate)
    rec = _project(changed, refs)["record"]
    if failed_predicate == "frameshift":
        assert rec["productive"] is False
    else:
        assert rec["productive"] is None and rec["stop_codon"] is None
        assert rec["junction"] == "TGNCCCTTT" and rec["junction_aa"] == "XPF"
    constrained = _compile(ge.PassPlan(), refs, productive=True, strict=strict)
    with pytest.raises(ge.StrictSamplingError, match="productive.final"):
        constrained.run_from(changed.final_simulation(), 0)


def test_verified_noncanonical_restoration_retains_junction_without_inventing_phase():
    refs, _, seed = _fixture(v="AAAGGG", verified=True)
    p = _project(_compile(ge.PassPlan(), refs).run_from(seed, 0), refs)
    rec = p["record"]
    assert rec["productive"] is False
    assert rec["vj_in_frame"] is None and rec["stop_codon"] is None
    assert rec["junction"] == "GGGCCCTTT" and rec["junction_aa"] == "GPF"
    assert (rec["junction_start"], rec["junction_end"]) == (3, 12)
    assert "anchor_residue_rejected" in rec["productivity_reasons"]
    assert "coding_phase_unresolved" in rec["productivity_reasons"]
    assert rec["productivity_policy"] == "duplex_observed_vdj_v1"


@pytest.mark.parametrize("strict", [False, True])
@pytest.mark.parametrize("prior_edit", ["inserted_suffix", "deleted_suffix"])
def test_j_tail_uses_observed_reference_endpoint_after_clipping_a_prior_edit(strict, prior_edit):
    refs, repair, seed = _fixture()
    edit = ge.PassPlan()
    insertion = prior_edit == "inserted_suffix"
    edit.push_corrupt_indel([(1, 1.0)], insertion_prob=1.0 if insertion else 0.0)
    edited = _find_outcome(
        _compile(edit, refs), seed,
        lambda out: out.trace().find("corrupt.indel.site[0]").value == (15 if insertion else 14),
    )
    assert _project(edited, refs)["record"]["productive"] is False
    clipping = ge.PassPlan()
    clipping.push_corrupt_3prime_loss([(1, 1.0)])
    # The allowed final observation loses the whole earlier terminal defect.
    # The productive contract must admit that actual event in both modes.
    clipped = _compile(clipping, refs, productive=True, strict=strict).run_from(
        edited.final_simulation(), 0,
    )
    p = _project(clipped, refs)
    rec = p["record"]
    assert rec["sequence"] == ("AAATGTCCCTTTAAA" if insertion else "AAATGTCCCTTTA")
    assert rec["j_cigar"] == ("6M" if insertion else "4M")
    assert rec["j_germline_end"] == (6 if insertion else 4)
    assert rec["j_trim_3"] == 0
    assert rec["end_loss_3_length"] == 1
    assert rec["duplex_repair"]["product"] == repair.run()
    _independent_productive_checks(p, phase=0)
