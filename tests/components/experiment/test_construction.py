"""Building an Experiment and compiling it.

What on() accepts, how steps chain, and what compile() produces.
"""
from __future__ import annotations
from GenAIRR import _engine as ge
import pytest
from GenAIRR import CompiledExperiment, Experiment
from tests.support.fixtures.experiment_dsl import _vdj_refdata, _vj_refdata

def test_experiment_on_returns_experiment_instance():
    exp = Experiment.on(_vj_refdata()).allow_curatable_refdata()
    assert isinstance(exp, Experiment)
    assert exp.chain_type == "vj"
    assert exp.step_count == 0


def test_experiment_on_passes_through_vdj_chain_type():
    assert Experiment.on(_vdj_refdata()).allow_curatable_refdata().chain_type == "vdj"


def test_experiment_on_rejects_unknown_config_name_string():
    # Strings are valid (config-name lookup), but unknown names get a
    # clear ValueError that lists alternatives.
    with pytest.raises(ValueError, match="Unknown config name"):
        Experiment.on("not_a_real_config_name")


def test_experiment_on_rejects_non_refdata_object():
    # Anything that is not a string / DataConfig / RefDataConfig is
    # rejected by the type guard with TypeError.
    with pytest.raises(TypeError, match="config-name string"):
        Experiment.on(42)


def test_experiment_on_rejects_none():
    with pytest.raises(TypeError, match="config-name string"):
        Experiment.on(None)


def test_experiment_on_resolves_known_config_name():
    # "human_igh" routes through the alias table to HUMAN_IGH_OGRDB
    # and produces a working VDJ Experiment with a non-empty refdata.
    exp = Experiment.on("human_igh")
    assert exp.chain_type == "vdj"
    assert exp.refdata.v_pool_size() > 0
    assert exp.refdata.j_pool_size() > 0


def test_experiment_on_accepts_dataconfig_object():
    # A loaded DataConfig is also a valid input.
    import GenAIRR

    cfg = GenAIRR.HUMAN_IGK_OGRDB
    exp = Experiment.on(cfg).allow_curatable_refdata()
    assert exp.chain_type == "vj"
    assert exp.refdata.v_pool_size() > 0
    assert exp.refdata.j_pool_size() > 0


def test_experiment_repr_includes_chain_type_and_step_count():
    exp = Experiment.on(_vj_refdata()).allow_curatable_refdata().recombine()
    r = repr(exp)
    assert "Experiment" in r
    assert "chain=vj" in r
    assert "steps=1" in r


def test_recombine_returns_same_experiment_for_chaining():
    exp = Experiment.on(_vj_refdata()).allow_curatable_refdata()
    same = exp.recombine()
    assert same is exp
    assert exp.step_count == 1


def test_recombine_can_be_called_multiple_times():
    # Two recombine() calls = two recombination shapes in the plan.
    # Not biologically meaningful, but the builder permits it.
    exp = Experiment.on(_vj_refdata()).allow_curatable_refdata().recombine().recombine()
    assert exp.step_count == 2


def test_compile_produces_compiled_experiment():
    compiled = Experiment.on(_vj_refdata()).allow_curatable_refdata().recombine().compile(allow_curatable_refdata=True)
    assert isinstance(compiled, CompiledExperiment)


def test_compile_vj_recombine_produces_five_pass_plan():
    compiled = Experiment.on(_vj_refdata()).allow_curatable_refdata().recombine().compile(allow_curatable_refdata=True)
    assert len(compiled.pass_plan) == 5


def test_compile_vdj_recombine_produces_eight_pass_plan():
    compiled = Experiment.on(_vdj_refdata()).allow_curatable_refdata().recombine().compile(allow_curatable_refdata=True)
    assert len(compiled.pass_plan) == 8


def test_compile_is_idempotent_and_returns_distinct_instances():
    exp = Experiment.on(_vj_refdata()).allow_curatable_refdata().recombine()
    a = exp.compile(allow_curatable_refdata=True)
    b = exp.compile(allow_curatable_refdata=True)
    assert a is not b
    assert a.simulator is not b.simulator
    assert len(a.pass_plan) == len(b.pass_plan) == 5


def test_compiled_experiment_repr_includes_plan_len_and_chain():
    compiled = Experiment.on(_vj_refdata()).allow_curatable_refdata().recombine().compile(allow_curatable_refdata=True)
    r = repr(compiled)
    assert "CompiledExperiment" in r
    assert "plan_len=5" in r
    assert "chain=vj" in r


def test_compiled_experiment_exposes_simulator_pass_summary_and_refdata():
    cfg = _vj_refdata()
    compiled = Experiment.on(cfg).allow_curatable_refdata().recombine().compile(allow_curatable_refdata=True)
    assert compiled.refdata is cfg
    assert isinstance(compiled.simulator, ge.CompiledSimulator)
    assert compiled.pass_plan == (
        "sample_allele.v",
        "sample_allele.j",
        "assemble.v",
        "generate_np.np1",
        "assemble.j",
    )


def test_vj_recombine_emits_canonical_pass_names():
    outcomes = Experiment.on(_vj_refdata()).allow_curatable_refdata().recombine().run(n=1, seed=0)
    assert outcomes[0].pass_names() == [
        "sample_allele.v",
        "sample_allele.j",
        "assemble.v",
        "generate_np.np1",
        "assemble.j",
    ]


def test_vdj_recombine_emits_canonical_pass_names():
    outcomes = Experiment.on(_vdj_refdata()).allow_curatable_refdata().recombine().run(n=1, seed=0)
    assert outcomes[0].pass_names() == [
        "sample_allele.v",
        "sample_allele.d",
        "sample_allele.j",
        "assemble.v",
        "generate_np.np1",
        "assemble.d",
        "generate_np.np2",
        "assemble.j",
    ]
