"""Independent literal oracles for persistent duplex observation projection."""
import pytest

from GenAIRR import _engine as ge
from GenAIRR.experimental import CoupledDuplexRepair, DuplexJoin, project_duplex_simulation

pytestmark = pytest.mark.behavior


def fixture(*, inverted=False, np1=""):
    v, d, j = ("AAATGTA" if np1 else "AAATGT"), "ACT", "TTTAAA"
    ref = ge.RefDataConfig.vdj()
    ref.add_v_allele("V*01", "V", v.encode(), anchor=3)
    ref.add_d_allele("D*01", "D", d.encode())
    ref.add_j_allele("J*01", "J", j.encode(), anchor=0)
    vd = DuplexJoin("vd").open("left", offset=0).open("right", offset=0)
    if np1:
        vd = vd.add_n("left", np1)
    vd = vd.synapse(right_start=len(v) + len(np1))
    if np1:
        vd = vd.fill("bottom")
    script = CoupledDuplexRepair(v, d, j, d_inverted=inverted).join(vd.ligate())
    script = script.join(
        DuplexJoin("dj").open("left", offset=0).open("right", offset=0)
        .synapse(right_start=len(d)).ligate()
    )
    calls = dict(v_call="V*01", d_call="D*01", j_call="J*01", reference_ends_are_coding_ends=True)
    return ref, script.seed(ref, **calls), script, calls


def compiled(ref, configure):
    plan = ge.PassPlan()
    configure(plan)
    return plan.compile(refdata=ref, strict=True)


def add_loss(plan, end, count):
    getattr(plan, f"push_corrupt_{end}prime_loss")([(count, 1.0)])


@pytest.mark.parametrize("inverted", [False, True])
@pytest.mark.parametrize("end,count", [(5, 6), (5, 9), (5, 15), (3, 6), (3, 9), (3, 15)])
def test_complete_source_clipping_keeps_generating_identity_without_false_deletions(inverted, end, count):
    ref, seed, _, _ = fixture(inverted=inverted)
    source = seed.bases().decode()
    outcome = compiled(ref, lambda p: add_loss(p, end, count)).run_from(seed, 13)
    projection = project_duplex_simulation(outcome.final_simulation(), ref)
    record = projection["record"]
    low, high = (count, 15) if end == 5 else (0, 15-count)
    assert record["sequence"] == source[low:high]
    assert record[f"end_loss_{end}_length"] == count
    for name, (start, stop) in zip(("v", "d", "j"), ((0, 6), (6, 9), (9, 15))):
        assert record[f"{name}_call"] == f"{name.upper()}*01"
        assert record[f"{name}_trim_5"] == record[f"{name}_trim_3"] == 0
        observed = max(low, start) < min(high, stop)
        assert projection["paths"][name]["available"] is observed
        if not observed:
            assert record[f"{name}_cigar"] is None
            assert record[f"{name}_germline_start"] is None
            assert record[f"{name}_germline_end"] is None
    assert record["p_d_5_length"] is None
    assert record["duplex_repair"]["product"]["sequence"] == source


def test_insertions_at_all_source_boundaries_and_terminal_edge_remain_in_alignment():
    ref, seed, _, _ = fixture()
    simulator = compiled(ref, lambda p: p.push_corrupt_indel([(1, 1.0)], insertion_prob=1.0))
    required = {0: ("v", "1I6M"), 6: ("d", "1I3M"), 9: ("j", "1I6M"), 15: ("j", "6M1I")}
    seen = set()
    for random_seed in range(128):
        outcome = simulator.run_from(seed, random_seed)
        edits = project_duplex_simulation(outcome.final_simulation(), ref)
        insertion = edits["record"]["duplex_repair"]["observation"]["edits"][0]
        at = insertion["action"]["at"] if "action" in insertion else insertion["at"]
        # Edit serialization nests the tagged operation under action.
        if at not in required:
            continue
        record = edits["record"]
        segment, cigar = required[at]
        assert record[f"{segment}_cigar"] == cigar
        assert record["sequence_alignment"].replace("-", "") == record["sequence"]
        assert record["n_indels"] == 1
        assert len(record["sequence"]) == 16
        seen.add(at)
        if seen == set(required):
            break
    assert seen == set(required)
    assert len(seed) == 15


def test_inserted_only_observation_after_complete_clipping_preserves_empty_source_order():
    ref, seed, _, _ = fixture()
    def configure(plan):
        add_loss(plan, 5, 15)
        plan.push_corrupt_indel([(1, 1.0)], insertion_prob=1.0)
    outcome = compiled(ref, configure).run_from(seed, 7)
    simulation = outcome.final_simulation()
    record = project_duplex_simulation(simulation, ref)["record"]
    assert len(record["sequence"]) == 1
    assert [(r.segment, r.start, r.end) for r in simulation.regions()] == [
        ("V", 0, 1), ("D", 1, 1), ("J", 1, 1),
    ]
    assert record["sequence_alignment"] is None
    assert record["np1"] is None and record["np2"] is None
    assert record["n_indels"] == 1 and record["end_loss_5_length"] == 15
    assert record["v_call"] == "V*01" and record["p_v_3_length"] is None


def test_np_amino_acids_follow_complete_coding_codons_in_founder_and_descendant():
    ref, seed, script, calls = fixture(np1="CCC")
    standalone = script.project(ref, **calls)
    descendant = project_duplex_simulation(seed, ref)
    for projection in (standalone, descendant):
        record = projection["record"]
        assert (record["v_sequence_end"], record["d_sequence_start"]) == (7, 10)
        assert record["np1"] == "CCC"
        assert record["np1_aa"] == ""  # q9 codon ends at12, outside NP [7,10)


def test_reverse_presentation_keeps_clipped_inverted_alignment_in_coding_orientation():
    ref, seed, _, _ = fixture(inverted=True)
    def configure(plan):
        add_loss(plan, 5, 1)
        plan.push_corrupt_rev_comp(1.0)
    outcome = compiled(ref, configure).run_from(seed, 19)
    projection = project_duplex_simulation(outcome.final_simulation(), ref)
    record = projection["record"]
    coding = seed.bases().decode()[1:]
    assert record["sequence"] == coding.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    assert record["sequence_alignment"] == coding
    assert record["rev_comp"] is True and record["d_inverted"] is True
    assert (record["d_sequence_start"], record["d_sequence_end"]) == (5, 8)
    assert [p["reference_position"] for p in projection["paths"]["d"]["body"]] == [2, 1, 0]
    assert record["duplex_repair"]["observation"]["coordinate_system"] == "zero_based_coding_oriented"


def test_attempt_counters_survive_empty_observation_and_noop_continuation():
    ref, seed, _, _ = fixture()
    def configure(plan):
        add_loss(plan, 5, 15)
        plan.push_corrupt_pcr([(2, 1.0)])
        plan.push_corrupt_quality([(3, 1.0)])
    outcome = compiled(ref, configure).run_from(seed, 2)
    child = compiled(ref, lambda p: None).run_from(outcome.final_simulation(), 3)
    for simulation in (outcome.final_simulation(), child.final_simulation()):
        record = project_duplex_simulation(simulation, ref)["record"]
        assert record["sequence"] == ""
        assert record["n_pcr_errors"] == 2 and record["n_quality_errors"] == 3
        assert record["n_mutations"] == 0 and record["n_indels"] == 0


def test_clipping_only_an_inserted_j_suffix_restores_founder_productivity():
    ref, seed, _, _ = fixture()
    insertion = compiled(ref, lambda p: p.push_corrupt_indel([(1, 1.0)], insertion_prob=1.0))
    inserted = insertion.run_from(seed, 0).final_simulation()
    before = project_duplex_simulation(inserted, ref)
    assert before["record"]["j_cigar"] == "6M1I"
    restored = compiled(ref, lambda p: add_loss(p, 3, 1)).run_from(inserted, 1).final_simulation()
    parent = project_duplex_simulation(seed, ref)["record"]
    child = project_duplex_simulation(restored, ref)["record"]
    for field in ("sequence", "junction", "junction_aa", "sequence_aa", "j_cigar", "vj_in_frame", "stop_codon", "productive"):
        assert child[field] == parent[field], field
    assert child["productive"] is True
    assert child["end_loss_3_length"] == 1 and child["n_j_indels"] == 1
