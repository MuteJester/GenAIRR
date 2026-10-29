"""Independent finite-history law, compile-preflight and DSL integration checks."""
from collections import Counter
from dataclasses import FrozenInstanceError, replace
import json
import warnings

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge
from GenAIRR.experimental import DuplexHistory, DuplexRepairModel, CoupledDuplexRepair, DuplexJoin
from tests.support.cases.duplex import references, script, history, model

pytestmark = pytest.mark.behavior


@pytest.mark.parametrize("weight", [-1, float("nan"), float("inf"), True, "invalid"])
def test_authored_history_weights_require_finite_nonnegative_values(weight):
    with pytest.raises(ValueError, match="weight"):
        history(weight=weight)


def test_model_is_immutable_roundtrippable_and_preserves_duplicate_products():
    r = references()
    values = [history("blunt", weight=1), history("restored", restore=True, weight=3)]
    m = model(r, values)
    values.clear()
    assert len(m.rows) == 2
    assert m.rows[0].script.run()["sequence"] == m.rows[1].script.run()["sequence"]
    assert DuplexRepairModel.from_dict(json.loads(m.to_json())) == m
    with pytest.raises(FrozenInstanceError):
        m.rows = ()
    with pytest.raises(ValueError, match="unique"):
        model(r, [history("same"), history("same", restore=True)])
    with pytest.raises(ValueError, match="positive"):
        model(r, [history(weight=0)])
    changed = m.to_dict()
    changed["unexpected"] = 0
    with pytest.raises(ValueError, match="schema"):
        DuplexRepairModel.from_dict(changed)


def test_duplex_dsl_has_one_owner_and_no_legacy_fallback_warnings():
    r = references()
    with warnings.catch_warnings(record=True) as seen:
        exp = ga.Experiment.on(r).recombine(duplex_repair=model(r))
        compiled = exp.compile()
    assert not seen
    assert compiled.pass_names == ("sample_allele.v", "sample_allele.d", "sample_allele.j", "duplex.recombine")
    assert "coupled duplex repair" in exp.describe()
    assert "NP1" not in exp.describe()
    rec = compiled.run_records(seed=3, validate_records=True)[0]
    assert rec["sequence"] == "AAATGTCCCTTTAAA"
    assert rec["productivity_policy"] == "duplex_observed_vdj_v1"
    assert rec["productive"] is True
    assert rec["duplex_repair"]["selection"]["history_id"] == "blunt"
    assert rec["p_v_3_length"] is None


@pytest.mark.parametrize("strict", [False, True])
@pytest.mark.parametrize("inverted", [False, True])
@pytest.mark.parametrize("order", [("vd", "dj"), ("dj", "vd")])
def test_generated_experiment_preserves_restored_landmarks_through_native_edits(strict, inverted, order):
    r = references()
    m = model(r, [history("restored", restore=True, inverted=inverted, order=order)])
    exp = ga.Experiment.on(r).recombine(duplex_repair=m)
    if inverted:
        exp.invert_d(prob=1)
    compiled = exp.mutate(model="uniform", count=2).end_loss_5prime(length=1).random_strand_orientation(prob=1).productive_only().compile()
    out = compiled.run(n=1, seed=11, strict=strict)[0]
    rec = out.airr_record(r)
    assert rec["productive"] is True
    assert rec["d_inverted"] is inverted
    assert rec["v_germline_start"] == 1 and rec["v_trim_5"] == 0
    assert rec["n_mutations"] == 2 and rec["mutation_rate"] == pytest.approx(2 / 15)
    assert rec["n_indels"] == 0 and rec["end_loss_5_length"] == 1
    assert rec["rev_comp"] is True
    founder = rec["duplex_repair"]["product"]
    assert founder == m.rows[0].script.run()
    assert founder["columns"][5]["top"]["kind"] == "templated_synthesis"
    assert rec["duplex_repair"]["selection"]["history_id"] == "restored"


def test_productivity_conditions_within_context_without_success_mass_reweighting():
    r = references(v_count=2)
    m = model(r, [history("v1_good", weight=1), history("v1_bad", weight=99, added="A"),
                  history("v2_good", weight=1, v_call="V2*01"),
                  history("v2_bad", weight=1, v_call="V2*01", added="A")])
    compiled = ga.Experiment.on(r).recombine(duplex_repair=m, v_allele_weights={"V2*01": 3}).productive_only().compile()
    sampled = Counter()
    for out in compiled.stream(n=4096, seed=17, strict=True):
        idx = out.trace().find("duplex_repair.v1.history").value
        sampled[idx] += 1
        assert idx in {0, 2}
    assert sampled[2] / sum(sampled.values()) == pytest.approx(0.75, abs=0.035)


def test_different_histories_with_identical_products_keep_their_natural_mass():
    r = references()
    m = model(r, [history("blunt", weight=1), history("restored", weight=3, restore=True)])
    compiled = ga.Experiment.on(r).recombine(duplex_repair=m).productive_only().compile()
    counts = Counter(out.trace().find("duplex_repair.v1.history").value
                     for out in compiled.stream(n=2048, seed=29, strict=True))
    assert set(counts) == {0, 1}
    assert counts[1] / sum(counts.values()) == pytest.approx(0.75, abs=0.045)


def test_missing_natural_context_fails_before_sampling_and_explicit_lock_is_respected():
    r = references(v_count=2)
    m = model(r)
    with pytest.raises(ValueError, match="(?i)(context|history)"):
        ga.Experiment.on(r).recombine(duplex_repair=m).compile()
    compiled = ga.Experiment.on(r).recombine(duplex_repair=m).restrict_alleles(v="V1*01").compile()
    assert compiled.run_records(seed=1)[0]["v_call"] == "V1*01"


def test_every_naturally_selectable_orientation_requires_a_history_table():
    r = references()
    m = model(r)
    with pytest.raises(ValueError, match="(?i)(context|orientation|inverted)"):
        ga.Experiment.on(r).recombine(duplex_repair=m).invert_d(prob=0.5).compile()
    both = model(r, [history("forward"), history("reverse", inverted=True)])
    compiled = ga.Experiment.on(r).recombine(duplex_repair=both).invert_d(prob=0.25).compile()
    choices = Counter(out.trace().find("sample_allele.d.inverted").value
                      for out in compiled.stream(n=16384, seed=3, strict=True))
    assert choices[True] / sum(choices.values()) == pytest.approx(0.25, abs=0.02)


@pytest.mark.parametrize("method,kwargs", [
    ("trim", {"enabled":False}), ("trim", {"v_3":[(0,1)]}),
    ("process_coding_ends", {"reference_ends_are_coding_ends":True}),
])
def test_explicit_legacy_processing_is_rejected_even_when_configured_as_zero(method, kwargs):
    r = references()
    exp = ga.Experiment.on(r).recombine(duplex_repair=model(r))
    with pytest.raises(ValueError, match="duplex"):
        getattr(exp, method)(**kwargs)
    with pytest.raises(ValueError, match="NP"):
        ga.Experiment.on(r).recombine(duplex_repair=model(r), np1_lengths=[(0,1)])


def test_reference_change_and_preflight_budget_fail_without_productive_fallback():
    r = references()
    m = model(r)
    r.add_v_allele("V2*01", "V2", b"AAATGT", anchor=3)
    with pytest.raises(ValueError, match="(?i)reference"):
        ga.Experiment.on(r).recombine(duplex_repair=m).compile()
    r = references()
    for budget in ({"max_preflight_operations":1}, {"max_preflight_bases":1}):
        with pytest.raises(ValueError, match="(?i)(resource|budget|limit)"):
            ga.Experiment.on(r).recombine(duplex_repair=model(r, **budget)).productive_only().compile()


def test_positive_unsupported_layout_is_not_hidden_by_productivity_conditioning():
    r = ge.RefDataConfig.vdj()
    for seg, seq in (("v", "ACG"), ("d", "GTA"), ("j", "CCC")):
        getattr(r, f"add_{seg}_allele")(f"{seg.upper()}*01", seg.upper(), seq.encode())
    overlap = CoupledDuplexRepair("ACG", "GTA", "CCC").join(
        DuplexJoin("vd").open("left", offset=0).open("right", offset=0)
        .resect("left", "bottom", count=1).resect("right", "top", count=1)
        .synapse(right_start=2).ligate()
    ).join(DuplexJoin("dj").open("left", offset=0).open("right", offset=0).synapse(right_start=3).ligate())
    m = model(r, [DuplexHistory("overlap", overlap, "V*01", "D*01", "J*01")])
    with pytest.raises(ValueError, match="(?i)(overlap|layout|capability)"):
        ga.Experiment.on(r).allow_curatable_refdata().recombine(duplex_repair=m).productive_only().compile()


@pytest.mark.parametrize("productive", [False, True])
@pytest.mark.parametrize("weights", [(5e-324, 1e308), (1.0, 1e-20)])
def test_floating_point_support_loss_is_an_error_before_sampling(productive, weights):
    r = references()
    m = model(r, [history("first", weight=weights[0]), history("second", weight=weights[1], restore=True)])
    exp = ga.Experiment.on(r).recombine(duplex_repair=m)
    if productive:
        exp.productive_only()
    with pytest.raises(ValueError, match="(?i)(weight|numeric|underflow|support)"):
        exp.compile()


def test_zero_weight_history_is_never_sampled_or_used_as_required_coverage():
    r = references()
    m = model(r, [history("chosen"), history("disabled", weight=0, added="A")])
    compiled = ga.Experiment.on(r).recombine(duplex_repair=m).compile()
    assert {out.trace().find("duplex_repair.v1.history").value for out in compiled.stream(n=16)} == {0}
    r = references(v_count=2)
    m = model(r, [history("chosen"), history("disabled", weight=0, v_call="V2*01")])
    with pytest.raises(ValueError, match="(?i)(context|positive|support|history)"):
        ga.Experiment.on(r).recombine(duplex_repair=m).compile()


@pytest.mark.parametrize("configure", [
    lambda exp: exp.contaminate(prob=0),
    lambda exp: exp.receptor_revision(prob=0),
    lambda exp: exp.paired_end(r1_length=10, insert_size=15),
])
def test_unsupported_suffix_is_rejected_at_compile_even_with_zero_effect_probability(configure):
    r = references()
    exp = configure(ga.Experiment.on(r).recombine(duplex_repair=model(r)))
    with pytest.raises(ValueError, match="(?i)(duplex|unsupported)"):
        exp.compile()


def _phased_fixture(chromosome_weights=(0.7, 0.3)):
    from GenAIRR.alleles import VAllele, DAllele, JAllele
    from GenAIRR.genotype import Genotype
    from GenAIRR.reference_models import ReferenceEmpiricalModels
    vs = [VAllele(f"IGHV-DUPLEX{i}*01", "AAATGT", 6, anchor_override=3) for i in (1, 2)]
    ds = [DAllele("IGHD-DUPLEX*01", "CCC", 3)]
    js = [JAllele(f"IGHJ-DUPLEX{i}*01", "TGGAAA", 6, anchor_override=0) for i in (1, 2)]
    cfg = ga.DataConfig(name="duplex-phased-control", metadata=ga.HUMAN_IGH_OGRDB.metadata,
                        v_alleles={v.gene:[v] for v in vs}, d_alleles={d.gene:[d] for d in ds},
                        j_alleles={j.gene:[j] for j in js}, reference_models=ReferenceEmpiricalModels())
    genotype = Genotype.from_dataconfig(cfg).complete_from_reference().with_subject("duplex-subject")
    for seg, alleles in (("V", vs), ("J", js)):
        genotype.delete_gene(alleles[1].gene, haplotype=0, segment=seg)
        genotype.delete_gene(alleles[0].gene, haplotype=1, segment=seg)
    genotype.chromosome_weights(*chromosome_weights)
    r = ga.dataconfig_to_refdata(cfg)
    rows = [DuplexHistory("h0_good", replace(script(), j="TGGAAA"), vs[0].name, ds[0].name, js[0].name),
            DuplexHistory("h0_bad", replace(script(added="A"), j="TGGAAA"), vs[0].name, ds[0].name, js[0].name, 99)]
    if chromosome_weights[1] > 0:
        rows.append(DuplexHistory("h1_good", replace(script(restore=True), j="TGGAAA"), vs[1].name, ds[0].name, js[1].name, 100))
    return cfg, genotype, model(r, rows), vs, js


@pytest.mark.parametrize("strict", [False, True])
@pytest.mark.parametrize("chromosome_weights", [(0.7, 0.3), (1.0, 0.0)])
def test_genotype_coverage_uses_expressible_joint_haplotypes_and_keeps_chromosome_weights(strict, chromosome_weights):
    cfg, genotype, m, vs, js = _phased_fixture(chromosome_weights)
    compiled = ga.Experiment.on(cfg).with_genotype(genotype).recombine(duplex_repair=m).productive_only().compile()
    counts = Counter()
    for out in compiled.stream(n=1024, seed=59, strict=strict):
        hap = out.trace().find("sample_haplotype").value
        counts[hap] += 1
        sim = out.final_simulation()
        assert compiled.refdata.v_allele(sim.v_allele_id()).name == vs[hap].name
        assert compiled.refdata.j_allele(sim.j_allele_id()).name == js[hap].name
    assert counts[1] / sum(counts.values()) == pytest.approx(chromosome_weights[1], abs=0.045)
    result = compiled.run_records(n=2, seed=3, strict=strict, validate_records=True)
    assert all(rec["subject_id"] == "duplex-subject" for rec in result)
    assert result.genotypes[0].subject_id == "duplex-subject"
