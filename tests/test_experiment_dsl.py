"""
Tests for the Experiment DSL — fluent protocol builder (clause-based API).

Covers construction, clause storage, repr formatting,
compile/run terminals, and error handling.
"""

import pytest


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def igh_config():
    from GenAIRR.data import HUMAN_IGH_OGRDB
    return HUMAN_IGH_OGRDB


# ---------------------------------------------------------------------------
# Construction and chaining
# ---------------------------------------------------------------------------

class TestExperimentConstruction:

    def test_on_creates_experiment(self):
        from GenAIRR import Experiment
        exp = Experiment.on("human_igh")
        assert exp._config == "human_igh"
        assert exp._recombine_clauses == []
        assert exp._mutate_clauses == []

    def test_chaining_returns_self(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import rate
        exp = Experiment.on("human_igh")
        same = exp.mutate(rate(0.01, 0.05))
        assert same is exp

    def test_recombine_stores_clauses(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import using, with_d_inversion
        exp = (Experiment.on("human_igh")
               .recombine(using(v="IGHV1-2*01"), with_d_inversion(0.2)))
        assert len(exp._recombine_clauses) == 2

    def test_mutate_stores_clauses(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import rate, model, with_antigen_selection
        exp = (Experiment.on("human_igh")
               .mutate(rate(0.02, 0.08), model("s5f"),
                       with_antigen_selection(0.7)))
        assert len(exp._mutate_clauses) == 3

    def test_prepare_stores_clauses(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import with_primer_mask, with_umi, with_pcr
        exp = (Experiment.on("human_igh")
               .prepare(with_primer_mask(), with_umi(10), with_pcr()))
        assert len(exp._prepare_clauses) == 3

    def test_sequence_stores_clauses(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import (with_5prime_loss, with_3prime_loss,
                                 with_quality_profile, with_reverse_complement)
        exp = (Experiment.on("human_igh")
               .sequence(with_5prime_loss(), with_3prime_loss(),
                         with_quality_profile(), with_reverse_complement()))
        assert len(exp._sequence_clauses) == 4

    def test_observe_stores_clauses(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import with_contaminants, with_indels, with_ns
        exp = (Experiment.on("human_igh")
               .observe(with_contaminants(0.05), with_indels(0.02),
                        with_ns(0.005)))
        assert len(exp._observe_clauses) == 3

    def test_cross_phase_rejection(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import rate, with_ns
        with pytest.raises(TypeError, match="RecombineClause"):
            Experiment.on("human_igh").recombine(rate(0.01, 0.05))
        with pytest.raises(TypeError, match="MutateClause"):
            Experiment.on("human_igh").mutate(with_ns())


# ---------------------------------------------------------------------------
# Repr / display
# ---------------------------------------------------------------------------

class TestExperimentRepr:

    def test_repr_contains_config(self):
        from GenAIRR import Experiment
        r = repr(Experiment.on("human_igh"))
        assert "human_igh" in r

    def test_repr_contains_clause_phases(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import (rate, model, with_d_inversion,
                                 with_5prime_loss, with_ns)
        exp = (Experiment.on("human_igh")
               .recombine(with_d_inversion(0.15))
               .mutate(rate(0.02, 0.08), model("s5f"))
               .sequence(with_5prime_loss())
               .observe(with_ns(0.01)))
        r = repr(exp)
        assert "Recombine" in r
        assert "Mutate" in r
        assert "Sequence" in r
        assert "Observe" in r

    def test_repr_shows_vdj_for_empty_recombine(self):
        from GenAIRR import Experiment
        r = repr(Experiment.on("human_igh"))
        assert "V(D)J rearrangement" in r

    def test_repr_hides_empty_phases(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import rate
        exp = Experiment.on("human_igh").mutate(rate(0.01, 0.05))
        r = repr(exp)
        assert "Prepare" not in r
        assert "Observe" not in r

    def test_str_equals_repr(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import rate
        exp = Experiment.on("human_igh").mutate(rate(0.01, 0.05))
        assert str(exp) == repr(exp)

    def test_repr_with_dataconfig_object(self, igh_config):
        from GenAIRR import Experiment
        r = repr(Experiment.on(igh_config))
        assert "Experiment" in r


# ---------------------------------------------------------------------------
# Build steps
# ---------------------------------------------------------------------------

class TestBuildSteps:

    def test_empty_experiment_builds_rearrange_only(self):
        from GenAIRR import Experiment
        from GenAIRR.steps import Rearrange
        steps = Experiment.on("human_igh")._build_steps()
        assert len(steps) == 1
        assert isinstance(steps[0], Rearrange)

    def test_implicit_rearrange_in_steps(self):
        from GenAIRR import Experiment
        from GenAIRR.steps import Rearrange
        from GenAIRR.ops import rate
        exp = Experiment.on("human_igh").mutate(rate(0.01, 0.05))
        steps = exp._build_steps()
        assert isinstance(steps[0], Rearrange)

    def test_step_types_correct(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import rate, with_antigen_selection, with_5prime_loss
        from GenAIRR.steps import (
            Rearrange, Mutate, SelectionPressure, Corrupt5Prime
        )
        exp = (Experiment.on("human_igh")
               .mutate(rate(0.01, 0.05), with_antigen_selection())
               .sequence(with_5prime_loss()))
        steps = exp._build_steps()
        assert isinstance(steps[0], Rearrange)
        assert isinstance(steps[1], Mutate)
        assert isinstance(steps[2], SelectionPressure)
        assert isinstance(steps[3], Corrupt5Prime)

    def test_full_pipeline_step_count(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import (
            with_d_inversion, rate, with_antigen_selection,
            with_primer_mask, with_umi, with_pcr,
            with_5prime_loss, with_quality_profile, with_reverse_complement,
            with_contaminants, with_indels, with_ns,
        )
        exp = (Experiment.on("human_igh")
               .recombine(with_d_inversion())
               .mutate(rate(0.01, 0.05), with_antigen_selection())
               .prepare(with_primer_mask(), with_umi(), with_pcr())
               .sequence(with_5prime_loss(), with_quality_profile(),
                         with_reverse_complement())
               .observe(with_contaminants(), with_indels(), with_ns()))
        steps = exp._build_steps()
        # Rearrange + DInversion + Mutate + Selection + 3 prepare + 3 sequence + 3 observe = 13
        assert len(steps) == 13


# ---------------------------------------------------------------------------
# Compile and run (integration)
# ---------------------------------------------------------------------------

class TestExperimentExecution:

    def test_run_basic(self):
        from GenAIRR import Experiment
        result = Experiment.on("human_igh").run(n=5, seed=42)
        assert len(result) == 5

    def test_run_with_mutation(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import rate
        result = (Experiment.on("human_igh")
                  .mutate(rate(0.02, 0.06))
                  .run(n=10, seed=42))
        assert len(result) == 10
        n_mutated = sum(1 for r in result if r["n_mutations"] > 0)
        assert n_mutated > 0

    def test_run_productive(self):
        from GenAIRR import Experiment
        result = Experiment.on("human_igh").run(n=20, seed=42, productive=True)
        for rec in result:
            assert rec["productive"] is True

    def test_compile_and_reuse(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import rate
        sim = (Experiment.on("human_igh")
               .mutate(rate(0.01, 0.05))
               .compile(seed=42))
        r1 = sim.simulate(n=5)
        r2 = sim.simulate(n=5)
        assert len(r1) == 5
        assert len(r2) == 5

    def test_run_with_corruption(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import with_5prime_loss
        result = (Experiment.on("human_igh")
                  .sequence(with_5prime_loss())
                  .run(n=10, seed=42))
        assert len(result) == 10

    def test_run_full_pipeline(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import (
            with_d_inversion, rate,
            with_umi, with_pcr,
            with_5prime_loss, with_reverse_complement,
            with_indels, with_ns,
        )
        result = (Experiment.on("human_igh")
                  .recombine(with_d_inversion())
                  .mutate(rate(0.01, 0.05))
                  .prepare(with_umi(12), with_pcr())
                  .sequence(with_5prime_loss(), with_reverse_complement())
                  .observe(with_indels(), with_ns())
                  .run(n=20, seed=42))
        assert len(result) == 20
        df = result.to_dataframe()
        assert df.shape[0] == 20
        assert "sequence" in df.columns

    def test_run_with_dataconfig_object(self, igh_config):
        from GenAIRR import Experiment
        result = Experiment.on(igh_config).run(n=5, seed=42)
        assert len(result) == 5

    def test_run_multi_species(self):
        from GenAIRR import Experiment
        for config_name in ["human_igh", "human_igk", "human_tcrb", "mouse_igh"]:
            result = Experiment.on(config_name).run(n=3, seed=42)
            assert len(result) == 3

    def test_seed_reproducibility(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import rate
        r1 = (Experiment.on("human_igh")
              .mutate(rate(0.02, 0.06))
              .run(n=10, seed=42))
        r2 = (Experiment.on("human_igh")
              .mutate(rate(0.02, 0.06))
              .run(n=10, seed=42))
        for a, b in zip(r1, r2):
            assert a["sequence"] == b["sequence"]
            assert a["v_call"] == b["v_call"]

    def test_run_with_uniform_model(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import rate, model
        result = (Experiment.on("human_igh")
                  .mutate(rate(0.02, 0.06), model("uniform"))
                  .run(n=10, seed=42))
        assert len(result) == 10
        n_mutated = sum(1 for r in result if r["n_mutations"] > 0)
        assert n_mutated > 0

    def test_run_with_isotype_rates(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import rate, with_isotype_rates
        result = (Experiment.on("human_igh")
                  .mutate(rate(0.01, 0.05), with_isotype_rates())
                  .run(n=10, seed=42))
        assert len(result) == 10
