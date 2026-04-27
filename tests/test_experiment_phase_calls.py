"""
Tests for the new clause-based phase-call methods on Experiment.

Phase 3 of the DSL redesign: verifies that .recombine(*clauses),
.mutate(*clauses), .prepare(*clauses), .sequence(*clauses), and
.observe(*clauses) work correctly.
"""

import pytest

from GenAIRR import Experiment
from GenAIRR.ops import (
    # Recombine
    using, with_d_inversion, with_receptor_revision,
    # Mutate
    model, rate, with_antigen_selection, with_isotype_rates,
    # Prepare
    with_umi, with_pcr, with_primer_mask,
    # Sequence
    paired_end, long_read,
    with_5prime_loss, with_3prime_loss,
    with_quality_profile, with_reverse_complement,
    # Observe
    with_contaminants, with_indels, with_ns,
    # Base types (for cross-phase rejection tests)
    RecombineClause, MutateClause, PrepareClause,
    SequenceClause, ObserveClause,
)


# ── Helpers ─────────────────────────────────────────────────────


def _exp():
    """Fresh experiment on human_igh (string config, no C needed)."""
    return Experiment.on("human_igh")


# =====================================================================
# Phase method basics: clauses are stored
# =====================================================================


class TestRecombinePhaseCall:

    def test_recombine_stores_clauses(self):
        exp = _exp().recombine(using(v="IGHV1-2*01"), with_d_inversion(0.2))
        assert len(exp._recombine_clauses) == 2

    def test_recombine_no_args_is_noop(self):
        """No-arg .recombine() is a no-op (returns self)."""
        exp = _exp().recombine()
        assert len(exp._recombine_clauses) == 0

    def test_recombine_multiple_calls_accumulate(self):
        exp = (_exp()
               .recombine(using(v="IGHV1-2*01"))
               .recombine(with_d_inversion(0.1)))
        assert len(exp._recombine_clauses) == 2

    def test_recombine_rejects_mutate_clause(self):
        with pytest.raises(TypeError, match="RecombineClause"):
            _exp().recombine(rate(0.01, 0.05))

    def test_recombine_rejects_observe_clause(self):
        with pytest.raises(TypeError, match="RecombineClause"):
            _exp().recombine(with_ns())


class TestMutatePhaseCall:

    def test_mutate_stores_clauses(self):
        exp = _exp().mutate(rate(0.02, 0.08), model("s5f"))
        assert len(exp._mutate_clauses) == 2

    def test_mutate_all_clause_types(self):
        exp = _exp().mutate(
            rate(0.01, 0.05),
            model("s5f"),
            with_isotype_rates(),
            with_antigen_selection(0.7),
        )
        assert len(exp._mutate_clauses) == 4

    def test_mutate_rejects_recombine_clause(self):
        with pytest.raises(TypeError, match="MutateClause"):
            _exp().mutate(using(v="IGHV1-2*01"))

    def test_mutate_rejects_sequence_clause(self):
        with pytest.raises(TypeError, match="MutateClause"):
            _exp().mutate(paired_end(300))


class TestPreparePhaseCall:

    def test_prepare_stores_clauses(self):
        exp = _exp().prepare(with_umi(12), with_pcr(), with_primer_mask())
        assert len(exp._prepare_clauses) == 3

    def test_prepare_no_args_is_noop(self):
        exp = _exp().prepare()
        assert len(exp._prepare_clauses) == 0

    def test_prepare_rejects_mutate_clause(self):
        with pytest.raises(TypeError, match="PrepareClause"):
            _exp().prepare(rate(0.01, 0.05))


class TestSequencePhaseCall:

    def test_sequence_stores_clauses(self):
        exp = _exp().sequence(
            with_5prime_loss(), with_3prime_loss(),
            paired_end(300), with_quality_profile(),
        )
        assert len(exp._sequence_clauses) == 4

    def test_sequence_no_args_is_noop(self):
        exp = _exp().sequence()
        assert len(exp._sequence_clauses) == 0

    def test_sequence_rejects_prepare_clause(self):
        with pytest.raises(TypeError, match="SequenceClause"):
            _exp().sequence(with_umi(12))


class TestObservePhaseCall:

    def test_observe_stores_clauses(self):
        exp = _exp().observe(with_contaminants(0.01), with_indels(), with_ns())
        assert len(exp._observe_clauses) == 3

    def test_observe_rejects_recombine_clause(self):
        with pytest.raises(TypeError, match="ObserveClause"):
            _exp().observe(with_d_inversion(0.1))

    def test_observe_no_args_is_noop(self):
        """observe() with no args is a no-op (returns self)."""
        exp = _exp()
        result = exp.observe()
        assert result is exp
        assert len(exp._observe_clauses) == 0


# =====================================================================
# Chaining across phases
# =====================================================================


class TestPhaseChainingNewDSL:

    def test_full_pipeline_chaining(self):
        exp = (
            _exp()
            .recombine(using(v="IGHV1-2*01"), with_d_inversion(0.15))
            .mutate(rate(0.02, 0.08), model("s5f"))
            .prepare(with_primer_mask(), with_umi(12))
            .sequence(with_5prime_loss(), paired_end(300))
            .observe(with_contaminants(0.01))
        )
        assert len(exp._recombine_clauses) == 2
        assert len(exp._mutate_clauses) == 2
        assert len(exp._prepare_clauses) == 2
        assert len(exp._sequence_clauses) == 2
        assert len(exp._observe_clauses) == 1

    def test_returns_self_for_chaining(self):
        exp = _exp()
        result = exp.mutate(rate(0.01, 0.05))
        assert result is exp

    def test_clauses_populated_after_phase_call(self):
        exp = _exp()
        assert len(exp._mutate_clauses) == 0
        exp.mutate(rate(0.01, 0.05))
        assert len(exp._mutate_clauses) == 1


# =====================================================================
# _build_steps() uses clause path
# =====================================================================


class TestBuildStepsClausePath:

    def test_build_steps_clause_path_basic(self):
        exp = _exp().mutate(rate(0.02, 0.08))
        steps = exp._build_steps()
        # Should have Rearrange + Mutate
        assert len(steps) == 2
        assert steps[0].__class__.__name__ == "Rearrange"
        assert steps[1].__class__.__name__ == "Mutate"
        assert steps[1].min_rate == 0.02
        assert steps[1].max_rate == 0.08

    def test_build_steps_clause_path_full(self):
        exp = (
            _exp()
            .recombine(with_d_inversion(0.15))
            .mutate(rate(0.01, 0.05))
            .prepare(with_umi(10))
            .sequence(with_5prime_loss())
            .observe(with_indels(0.02))
        )
        steps = exp._build_steps()
        step_names = [s.__class__.__name__ for s in steps]
        assert "Rearrange" in step_names
        assert "SimulateDGeneInversion" in step_names
        assert "Mutate" in step_names
        assert "SimulateUMI" in step_names
        assert "Corrupt5Prime" in step_names
        assert "InsertIndels" in step_names

    def test_build_steps_empty_experiment(self):
        """Empty experiment produces just Rearrange."""
        exp = _exp()
        steps = exp._build_steps()
        step_names = [s.__class__.__name__ for s in steps]
        assert step_names == ["Rearrange"]

    def test_build_steps_csr_via_clauses(self):
        exp = _exp().mutate(rate(0.02, 0.08), with_isotype_rates())
        steps = exp._build_steps()
        step_names = [s.__class__.__name__ for s in steps]
        assert "SimulateCSR" in step_names
        assert "Mutate" not in step_names

    def test_build_steps_selection_via_clauses(self):
        exp = _exp().mutate(
            rate(0.01, 0.05), with_antigen_selection(strength=0.7)
        )
        steps = exp._build_steps()
        step_names = [s.__class__.__name__ for s in steps]
        assert "Mutate" in step_names
        assert "SelectionPressure" in step_names


# =====================================================================
# __repr__ for clause-based experiments
# =====================================================================


class TestClauseRepr:

    def test_repr_shows_clause_phases(self):
        exp = (
            _exp()
            .mutate(rate(0.02, 0.08), model("s5f"))
            .observe(with_ns(0.01))
        )
        r = repr(exp)
        assert "human_igh" in r
        assert "Mutate" in r
        assert "Observe" in r
        assert "Rate:" in r
        assert "Model: S5F" in r
        assert "N bases" in r

    def test_repr_always_shows_recombine(self):
        exp = _exp().mutate(rate(0.01, 0.05))
        r = repr(exp)
        assert "Recombine" in r
        assert "V(D)J rearrangement" in r

    def test_repr_hides_empty_phases(self):
        exp = _exp().mutate(rate(0.01, 0.05))
        r = repr(exp)
        # Prepare, Sequence, Observe should not appear
        assert "Prepare" not in r
        assert "Sequence" not in r
        assert "Observe" not in r

    def test_repr_mutate_shows_rate(self):
        exp = _exp().mutate(rate(0.01, 0.05))
        r = repr(exp)
        assert "Mutate" in r
        assert "Rate:" in r

    def test_str_equals_repr(self):
        exp = _exp().mutate(rate(0.01, 0.05))
        assert str(exp) == repr(exp)


# =====================================================================
# Cross-phase rejection comprehensive
# =====================================================================


class TestCrossPhaseRejection:
    """Every phase method rejects clauses from other phases."""

    @pytest.mark.parametrize("wrong_clause", [
        rate(0.01, 0.05), with_umi(12), paired_end(300), with_ns(),
    ])
    def test_recombine_rejects_non_recombine(self, wrong_clause):
        with pytest.raises(TypeError):
            _exp().recombine(wrong_clause)

    @pytest.mark.parametrize("wrong_clause", [
        using(v="X"), with_umi(12), paired_end(300), with_ns(),
    ])
    def test_mutate_rejects_non_mutate(self, wrong_clause):
        with pytest.raises(TypeError):
            _exp().mutate(wrong_clause)

    @pytest.mark.parametrize("wrong_clause", [
        using(v="X"), rate(0.01, 0.05), paired_end(300), with_ns(),
    ])
    def test_prepare_rejects_non_prepare(self, wrong_clause):
        with pytest.raises(TypeError):
            _exp().prepare(wrong_clause)

    @pytest.mark.parametrize("wrong_clause", [
        using(v="X"), rate(0.01, 0.05), with_umi(12), with_ns(),
    ])
    def test_sequence_rejects_non_sequence(self, wrong_clause):
        with pytest.raises(TypeError):
            _exp().sequence(wrong_clause)

    @pytest.mark.parametrize("wrong_clause", [
        using(v="X"), rate(0.01, 0.05), with_umi(12), paired_end(300),
    ])
    def test_observe_rejects_non_observe(self, wrong_clause):
        with pytest.raises(TypeError):
            _exp().observe(wrong_clause)


# =====================================================================
# T1-9: Step classes are not exported at the top level
# =====================================================================


class TestStepNotTopLevelExport:
    """Internal Step classes (Mutate, Rearrange, Corrupt5Prime, etc.)
    must NOT be importable from the top-level GenAIRR package.

    They were trap-bait — tab-completion surfaced them next to
    Experiment, but passing them to phase methods raised TypeError
    because the phases only accept *Clause from GenAIRR.ops."""

    def test_mutate_not_at_top_level(self):
        import GenAIRR
        assert not hasattr(GenAIRR, "Mutate"), \
            "GenAIRR.Mutate is internal — should not be top-level"

    def test_rearrange_not_at_top_level(self):
        import GenAIRR
        assert not hasattr(GenAIRR, "Rearrange")

    def test_step_base_not_at_top_level(self):
        import GenAIRR
        assert not hasattr(GenAIRR, "Step")

    def test_all_internal_steps_not_at_top_level(self):
        """No Step subclass leaked into the top-level namespace."""
        import GenAIRR
        leaked = []
        for name in [
            "Step", "Rearrange", "Mutate", "SimulateCSR", "SelectionPressure",
            "SimulateDGeneInversion", "SimulateReceptorRevision",
            "Corrupt5Prime", "Corrupt3Prime", "CorruptQuality",
            "SimulatePairedEnd", "PCRAmplification", "SimulateUMI",
            "PrimerMask", "ReverseComplement", "SpikeContaminants",
            "SkewBaseComposition", "InsertIndels", "InsertNs",
        ]:
            if hasattr(GenAIRR, name):
                leaked.append(name)
        assert not leaked, f"top-level Step exports leaked: {leaked}"

    def test_subpackage_access_still_works(self):
        """Power users / tests can still import via GenAIRR.steps.*"""
        from GenAIRR.steps import (
            Step, Rearrange, Mutate, Corrupt5Prime, InsertIndels
        )
        # And they can be instantiated
        assert Mutate(min_rate=0.01, max_rate=0.05) is not None
        assert Step is not None

    def test_passing_step_to_phase_gives_redirect_hint(self):
        """When a user does mistakenly pass a Step (via the subpackage
        import) to a phase, the TypeError must redirect them to
        GenAIRR.ops — that's the audit's signposting requirement."""
        from GenAIRR.steps import Mutate as _MutateStep
        with pytest.raises(TypeError) as excinfo:
            _exp().mutate(_MutateStep(min_rate=0.01, max_rate=0.05))
        msg = str(excinfo.value)
        assert "Step" in msg, f"hint missing Step keyword: {msg}"
        assert "GenAIRR.ops" in msg, f"hint missing GenAIRR.ops redirect: {msg}"
        assert "rate" in msg, f"hint missing example clause name: {msg}"

    def test_normal_clause_error_is_unchanged(self):
        """Passing a wrong-CLASS clause (not a Step) still gets the
        plain TypeError — no Step hint."""
        with pytest.raises(TypeError) as excinfo:
            _exp().mutate(with_umi(12))   # Prepare-clause into mutate
        msg = str(excinfo.value)
        assert "MutateClause" in msg
        assert "Step" not in msg, \
            f"wrong-class clause should not get Step hint: {msg}"
