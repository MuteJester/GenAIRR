"""End-outcome ancestry across genotype, observation, continuation and export."""
import csv
import json

import pytest

import GenAIRR as ga

from tests.support.cases.productive import reverse_complement, synthetic_cartridge
from tests.support.fixtures.phased_anchors import _phased_fixture

pytestmark = pytest.mark.composition


def configured(*, d=None, v="AAACCCTGT", **ends):
    cfg = synthetic_cartridge({}, v=v, d=d)
    lengths = {"np1_lengths": [(0, 1)]}
    if d is not None:
        lengths["np2_lengths"] = [(0, 1)]
    return ga.Experiment.on(cfg).recombine(**lengths).process_coding_ends(
        reference_ends_are_coding_ends=True, **ends
    )


@pytest.mark.parametrize("inverted", [False, True])
def test_original_p_can_survive_complete_deletion_from_opposite_d_end(inverted):
    exp = configured(d="ACG", d_5=[(3, 1)], d_3=[(-2, 1)]).invert_d(prob=float(inverted))
    result = exp.run_records(n=4, seed=19)
    for rec in result:
        p = "CG" if not inverted else "CG"  # RC of the original terminal CG.
        assert rec["sequence"].upper() == "AAACCCTGT" + p + "TGGAAA"
        assert rec["d_germline_start"] is None
        assert rec["d_germline_end"] is None
        meta = rec["coding_end_processing"]
        origins = [o for o in meta["p_origins"] if o["source_segment"] == "D"]
        assert len(origins) == 2
        assert {o["template_position"] for o in origins} == {1, 2}
        assert {o["canonical_end"] for o in origins} == {"3"}
    assert result.validate_records(exp.refdata).ok


@pytest.mark.parametrize("strict", [False, True])
def test_productive_outcomes_honor_linked_genotype_and_true_replay(strict):
    cfg, genotype, expected_v, expected_j = _phased_fixture()
    exp = (ga.Experiment.on(cfg).with_genotype(genotype)
           .recombine(np1_lengths=[(0, 1)])
           .process_coding_ends(v_3=[(-1, 1)], reference_ends_are_coding_ends=True)
           .productive_only())
    result = exp.run_records(n=24, seed=43, strict=strict)
    for rec in result:
        assert (rec["haplotype"], rec["v_call"], rec["j_call"]) == (1, expected_v, expected_j)
        assert rec["sequence"].upper() == "TGTACCTTTGGG"
    assert result.validate_records(exp.refdata).ok
    compiled = exp.compile()
    fresh = compiled.simulator.run(seed=43, strict=strict)
    replay = compiled.simulator.replay_from_trace_file(compiled.simulator.trace_file_from(fresh, seed=43))
    assert fresh.airr_record(compiled.refdata) == replay.airr_record(compiled.refdata)


@pytest.mark.parametrize("feature", ["substitutions", "indels", "loss", "mask", "reverse", "revision", "clones"])
def test_template_identity_survives_edits_and_continuation(feature):
    if feature == "revision":
        exp = ga.Experiment.on("human_igh").recombine().process_coding_ends(
            v_3=[(-2, 1)], d_5=[(-1, 1)], d_3=[(-2, 1)], j_5=[(-1, 1)],
            reference_ends_are_coding_ends=True,
        ).receptor_revision(prob=1.0)
    else:
        exp = configured(d="ACGTAC", v_3=[(-2, 1)], d_5=[(-1, 1)], d_3=[(-2, 1)], j_5=[(-1, 1)])
        if feature == "substitutions": exp.mutate(model="uniform", count=10)
        elif feature == "indels": exp.polymerase_indels(count=8)
        elif feature == "loss": exp.end_loss_5prime(length=12).end_loss_3prime(length=10)
        elif feature == "mask": exp.ambiguous_base_calls(count=10)
        elif feature == "reverse": exp.random_strand_orientation(prob=1.0)
        elif feature == "clones": exp.expand_clones(n_clones=2, per_clone=3).polymerase_indels(count=4)
    result = exp.run_records(seed=29, **({} if feature == "clones" else {"n": 16}))
    for rec in result:
        meta = rec["coding_end_processing"]
        assert sum(meta["retained_p_at_join"].values()) == 6
        assert len(meta["p_origins"]) == 6
        coding = reverse_complement(rec["sequence"]) if rec["rev_comp"] else rec["sequence"]
        for origin in meta["p_origins"]:
            pos = origin["query_position"]
            if pos is None:
                assert origin["current_base"] is None
            else:
                assert coding[pos].upper() == origin["current_base"].upper()
        if feature == "revision":
            v_choice = next(c for c in meta["choices"] if c["segment"] == "V")
            assert v_choice["original_allele"] == rec["original_v_call"]
            assert v_choice["original_allele"] != rec["v_call"]
    assert result.validate_records(exp.refdata).ok


@pytest.mark.parametrize("strict", [False, True])
def test_tabular_export_retains_json_ancestry_and_its_explicit_coordinates(tmp_path, strict):
    exp = configured(v_3=[(-1, 1)], j_5=[(-2, 1)])
    result = exp.run_records(n=1, seed=2)
    output = tmp_path / "processing.tsv"
    result.to_tsv(str(output), airr_strict=strict)
    with output.open(newline="") as source:
        exported = next(csv.DictReader(source, delimiter="\t"))
    assert json.loads(exported["coding_end_processing"]) == result[0]["coding_end_processing"]
    assert json.loads(exported["coding_end_processing"])["coordinate_system"] == "zero_based_coding_oriented"
    assert int(exported["j_sequence_start"]) == result[0]["j_sequence_start"] + int(strict)


def test_public_validator_detects_nested_ancestry_tampering_and_metadata_is_reserved():
    exp = configured(v_3=[(-1, 1)])
    with pytest.raises(ValueError, match="reserved"):
        ga.Experiment.on("human_igk").with_metadata(coding_end_processing={})
    result = exp.run_records(n=1, seed=0)
    result[0]["coding_end_processing"]["p_origins"][0]["template_position"] -= 1
    report = result.validate_records(exp.refdata)
    assert not report.ok
    assert any(i.get("field") == "coding_end_processing" for failure in report.failures for i in failure["issues"])
