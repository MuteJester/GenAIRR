"""Independent current-query/immutable-founder checks for native continuations."""
import csv
import json
import re

import pytest

from GenAIRR import _engine as ge
from GenAIRR.experimental import (
    CoupledDuplexRepair, DuplexJoin, export_duplex_projection,
    project_duplex_simulation,
)

pytestmark = pytest.mark.behavior


def fixture(*, inverted=False, order=("vd", "dj")):
    r = ge.RefDataConfig.vdj()
    refs = ("AAATGT", "CCA", "TTTAAA")
    r.add_v_allele("V*01", "V", refs[0].encode(), anchor=3)
    r.add_d_allele("D*01", "D", refs[1].encode())
    r.add_j_allele("J*01", "J", refs[2].encode(), anchor=0)
    script = CoupledDuplexRepair(*refs, d_inverted=inverted)
    for name in order:
        script = script.join(DuplexJoin(name).open("left", offset=0)
            .open("right", offset=0).synapse(right_start=6 if name == "vd" else 3).ligate())
    seed = script.seed(r, v_call="V*01", d_call="D*01", j_call="J*01",
                      reference_ends_are_coding_ends=True)
    return r, script, seed


def continue_with(r, seed, method, *args, random_seed=5, **kwargs):
    plan = ge.PassPlan()
    getattr(plan, method)(*args, **kwargs)
    return plan.compile(refdata=r, strict=True).run_from(seed, random_seed).final_simulation()


def observation(projection):
    return projection["record"]["duplex_repair"]["observation"]


def assert_paths_consume_current_query(p):
    rec = p["record"]
    if rec["sequence_alignment"] is not None:
        expected = rec["sequence"]
        if rec["rev_comp"]:
            expected = expected.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]
        assert rec["sequence_alignment"].replace("-", "") == expected
    for name, path in p["paths"].items():
        if not path["available"]:
            continue
        ops = [(int(n), op) for n, op in re.findall(r"(\d+)([MID])", path["cigar"])]
        qlen = sum(n for n, op in ops if op in "MI")
        glen = sum(n for n, op in ops if op in "MD")
        assert glen == path["reference_end"] - path["reference_start"]
        if qlen:
            assert qlen == path["query_end"] - path["query_start"]
        else:
            assert path["query_start"] is path["query_end"] is None
        assert len(path["germline_alignment"]) == len(path["sequence_alignment"])


@pytest.mark.parametrize("inverted", [False, True])
@pytest.mark.parametrize("order", [("vd", "dj"), ("dj", "vd")])
def test_seed_roundtrip_preserves_bound_projection_and_canonical_coordinates(inverted, order):
    r, script, seed = fixture(inverted=inverted, order=order)
    initial = script.project(r, v_call="V*01", d_call="D*01", j_call="J*01",
                             reference_ends_are_coding_ends=True)
    p = project_duplex_simulation(seed, r)
    for field, value in initial["record"].items():
        if field not in {"sequence_id", "duplex_repair", "productivity_reasons"}:
            assert p["record"][field] == value, field
    assert observation(p)["founder_id"] == initial["record"]["sequence_id"]
    assert p["record"]["duplex_repair"]["product"] == script.run()
    assert [x["origin"] for x in observation(p)["links"]] == list(range(15))
    assert [seed.germline_position(i) for i in range(6, 9)] == ([2, 1, 0] if inverted else [0, 1, 2])
    assert_paths_consume_current_query(p)


@pytest.mark.parametrize("method,args,kwargs", [
    ("push_mutate_uniform", ([(3, 1.0)],), {}),
    ("push_mutate_s5f", ([(3, 1.0)], [1.0]*1024, [0.25]*4096), {}),
    ("push_corrupt_pcr", ([(3, 1.0)],), {}),
    ("push_corrupt_quality", ([(3, 1.0)],), {}),
    ("push_corrupt_ns", ([(3, 1.0)],), {}),
    ("push_corrupt_indel", ([(3, 1.0)],), {"insertion_prob": 1.0}),
    ("push_corrupt_indel", ([(3, 1.0)],), {"insertion_prob": 0.0}),
    ("push_corrupt_5prime_loss", ([(8, 1.0)],), {}),
    ("push_corrupt_3prime_loss", ([(8, 1.0)],), {}),
    ("push_corrupt_rev_comp", (1.0,), {}),
])
def test_audited_native_continuations_preserve_founder_and_query_ancestry(method, args, kwargs):
    r, script, seed = fixture(inverted=True)
    initial = project_duplex_simulation(seed, r)
    child = continue_with(r, seed, method, *args, **kwargs)
    p = project_duplex_simulation(child, r)
    assert project_duplex_simulation(seed, r) == initial
    assert p["record"]["duplex_repair"]["product"] == script.run()
    assert p["record"]["duplex_repair"]["counts"] == initial["record"]["duplex_repair"]["counts"]
    assert len(observation(p)["links"]) == len(child)
    assert observation(p)["founder_id"] == observation(initial)["founder_id"]
    assert_paths_consume_current_query(p)
    assert p["record"]["v_call"] == "V*01"
    assert p["record"]["d_call"] == "D*01"
    assert p["record"]["j_call"] == "J*01"


def test_mutation_burden_uses_stage_molecule_not_later_short_read():
    r, _, seed = fixture()
    child = continue_with(r, seed, "push_mutate_uniform", [(3, 1.0)])
    child = continue_with(r, child, "push_corrupt_3prime_loss", [(6, 1.0)])
    p = project_duplex_simulation(child, r)
    rec = p["record"]
    # The subject is the denominator -- the stage molecule (15), not the
    # 9-base read left after the 3' loss. Stated against the record's own
    # count, since `count=3` is three attempts.
    assert 0 < rec["n_mutations"] <= 3
    assert rec["mutation_rate"] == pytest.approx(rec["n_mutations"] / 15)
    assert rec["sequence_length"] == 9
    assert rec["end_loss_3_length"] == 6
    assert rec["n_indels"] == 0
    assert rec["j_germline_end"] is None
    assert rec["j_trim_3"] == 0
    assert all(x["cause"] == "observation_loss" for x in observation(p)["removed"])


def test_complete_deletion_exports_reference_path_without_fake_query_interval(tmp_path):
    r, _, seed = fixture()
    child = continue_with(r, seed, "push_corrupt_indel", [(15, 1.0)], insertion_prob=0.0)
    p = project_duplex_simulation(child, r)
    assert p["record"]["sequence"] == ""
    assert p["record"]["n_indels"] == 15
    assert p["record"]["d_cigar"] == "3D"
    assert p["record"]["d_sequence_start"] is None
    assert p["record"]["d_germline_start"] == 0
    assert p["record"]["d_germline_end"] == 3
    assert_paths_consume_current_query(p)
    path = tmp_path / "deleted.tsv"
    export_duplex_projection(p, path)
    export_duplex_projection(p, path, format="tsv")
    row = next(csv.DictReader(path.open(), delimiter="\t"))
    assert row["d_cigar"] == "3D"
    assert row["d_sequence_start"] == ""
    assert row["d_germline_start"] == "1"
    assert json.loads(row["duplex_repair"])["counts"]["founder_query"]["columns"] == 15


def test_forks_keep_immutable_founder_and_separate_observation_links():
    r, _, seed = fixture()
    a = continue_with(r, seed, "push_corrupt_indel", [(2, 1.0)], insertion_prob=1.0, random_seed=2)
    b = continue_with(r, seed, "push_corrupt_5prime_loss", [(2, 1.0)])
    pa, pb = project_duplex_simulation(a, r), project_duplex_simulation(b, r)
    assert len(observation(pa)["links"]) == 17
    assert len(observation(pb)["links"]) == 13
    assert observation(pa)["founder_id"] == observation(pb)["founder_id"]
    assert pa["record"]["sequence_id"] != pb["record"]["sequence_id"]
    assert observation(pa)["insertion_identity_scope"] == "observation_branch_local"
    assert project_duplex_simulation(seed, r)["record"]["sequence_length"] == 15


def test_sampled_error_attempts_survive_empty_pool_and_accumulate_across_continuations():
    r, _, seed = fixture()
    child = continue_with(r, seed, "push_corrupt_5prime_loss", [(15, 1.0)])
    child = continue_with(r, child, "push_corrupt_pcr", [(2, 1.0)])
    first = project_duplex_simulation(child, r)
    child = continue_with(r, child, "push_corrupt_pcr", [(3, 1.0)])
    child = continue_with(r, child, "push_corrupt_quality", [(4, 1.0)])
    final = project_duplex_simulation(child, r)
    assert first["record"]["n_pcr_errors"] == 2
    assert final["record"]["n_pcr_errors"] == 5
    assert final["record"]["n_quality_errors"] == 4
    assert not [x for x in observation(final)["edits"] if x["action"]["kind"] == "substitution"]
    assert final["record"]["sequence_id"] != first["record"]["sequence_id"]


@pytest.mark.parametrize("phase,expected", [(0, ""), (1, "P")])
def test_np_amino_acids_use_whole_query_coding_phase(phase, expected):
    import hashlib
    r = ge.RefDataConfig.vdj()
    for name, seq, anchor in [("v", "AAATGTA", 3 + phase), ("d", "CCA", None), ("j", "TTTAAA", 0)]:
        getattr(r, f"add_{name}_allele")(f"{name.upper()}*01", name.upper(), seq.encode(), **({"anchor": anchor} if name != "d" else {}))
        if name != "d":
            r.set_allele_evidence(name, 0, json.dumps({
                "schema_version": 1, "sequence_sha256": "sha256:" + hashlib.sha256(seq.encode()).hexdigest(),
                "anchor_status": "verified", "coding_phase": phase if name == "v" else 0,
                "boundary_5": "complete", "boundary_3": "complete",
                "source": "literal control", "method": "authored_coordinate", "notes": [],
            }))
    script = CoupledDuplexRepair("AAATGTA", "CCA", "TTTAAA").join(
        DuplexJoin("vd").open("left", offset=0).open("right", offset=0)
        .add_n("left", "CCC").add_n("right", "GGG").synapse(right_start=10).ligate()
    ).join(DuplexJoin("dj").open("left", offset=0).open("right", offset=0).synapse(right_start=3).ligate())
    args = dict(v_call="V*01", d_call="D*01", j_call="J*01", reference_ends_are_coding_ends=True)
    standalone = script.project(r, **args)
    seeded = project_duplex_simulation(script.seed(r, **args), r)
    assert standalone["record"]["np1"] == seeded["record"]["np1"] == "CCC"
    assert standalone["record"]["np1_aa"] == seeded["record"]["np1_aa"] == expected


def test_recalled_alleles_never_replace_generating_calls_or_alignment_reference():
    r, script, _ = fixture()
    r.add_v_allele("V*02", "V", b"CCCTGT", anchor=3)
    seed = script.seed(r, v_call="V*01", d_call="D*01", j_call="J*01", reference_ends_are_coding_ends=True)
    plan = ge.PassPlan()
    plan.push_mutate_uniform([(5, 1.0)], segment_rates=[1.0, 0.0, 0.0, 0.0])
    compiled = plan.compile(refdata=r, strict=True)
    for random_seed in range(128):
        child = compiled.run_from(seed, random_seed).final_simulation()
        p = project_duplex_simulation(child, r)
        recalled = p["record"]["recalled_v_call"]
        if recalled and recalled.split(",")[0] == "V*02":
            assert p["record"]["v_call"] == "V*01"
            assert p["paths"]["v"]["germline_alignment"] == "AAATGT"
            assert p["record"]["v_germline_start"] == 0
            assert p["record"]["v_germline_end"] == 6
            assert p["record"]["v_cigar"] == "6M"
            break
    else:
        pytest.fail("fixed native seeds did not yield the authored V*02 re-call")


@pytest.mark.parametrize("inverted", [False, True])
def test_structural_edits_keep_live_recall_cache_equal_to_fresh_walk(inverted):
    r, _, seed = fixture(inverted=inverted)
    plan = ge.PassPlan()
    plan.push_corrupt_indel([(2, 1.0)], insertion_prob=0.5)
    compiled = plan.compile(refdata=r, strict=True)
    for random_seed in range(64):
        out = compiled.run_from(seed, random_seed)
        for parity in out.check_live_call_cache_parity(r):
            assert parity["tie_set_matches"], (random_seed, parity)
            assert parity["hypothesis_bounds_match"] is not False, (random_seed, parity)


def test_partial_terminal_codon_after_two_sided_clipping_has_empty_translation():
    r, _, seed = fixture()
    plan = ge.PassPlan()
    plan.push_corrupt_5prime_loss([(1, 1.0)])
    plan.push_corrupt_3prime_loss([(13, 1.0)])
    child = plan.compile(refdata=r, strict=True).run_from(seed, 4).final_simulation()
    p = project_duplex_simulation(child, r)
    assert p["record"]["sequence"] == "A"
    assert p["record"]["sequence_aa"] == ""
    assert p["record"]["stop_codon"] is None
