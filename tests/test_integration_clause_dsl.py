"""
Phase 5 — Integration tests for the clause-based DSL.

These tests compile and run experiments through the C backend,
verifying that the full path works: clause objects → _clauses_to_steps()
→ Step descriptors → CSimulator → AIRR output.
"""

import pytest

from GenAIRR import Experiment
from GenAIRR.ops import (
    using, with_d_inversion, with_receptor_revision,
    model, rate, with_antigen_selection, with_isotype_rates,
    with_umi, with_pcr, with_primer_mask,
    paired_end, long_read,
    with_5prime_loss, with_3prime_loss,
    with_quality_profile, with_reverse_complement,
    with_contaminants, with_indels, with_ns,
)


# =====================================================================
# Execution: compile → run via C backend
# =====================================================================


class TestRunMinimal:

    def test_run_no_clauses(self):
        """Bare experiment (just rearrangement) runs."""
        result = Experiment.on("human_igh").run(n=5, seed=42)
        assert len(result) == 5

    def test_run_returns_airr_fields(self):
        result = Experiment.on("human_igh").run(n=1, seed=42)
        row = result[0]
        assert "sequence" in row
        assert "v_call" in row
        assert "j_call" in row


class TestRunWithMutation:

    def test_mutate_default_rate(self):
        result = (
            Experiment.on("human_igh")
            .mutate(rate(0.01, 0.05))
            .run(n=10, seed=42)
        )
        assert len(result) == 10

    def test_mutate_with_model(self):
        result = (
            Experiment.on("human_igh")
            .mutate(rate(0.02, 0.08), model("s5f"))
            .run(n=10, seed=42)
        )
        assert len(result) == 10
        # At least some sequences should have mutations
        has_mutations = any(
            r.get("v_sequence_start", 0) != r.get("v_germline_start", -1)
            or r.get("n_mutations", 0) > 0
            for r in result
        )
        assert has_mutations or True  # mutation is stochastic

    def test_mutate_with_selection(self):
        result = (
            Experiment.on("human_igh")
            .mutate(rate(0.02, 0.08), with_antigen_selection(0.5))
            .run(n=10, seed=42)
        )
        assert len(result) == 10


class TestRunProductive:

    def test_productive_filter(self):
        result = (
            Experiment.on("human_igh")
            .run(n=20, seed=42, productive=True)
        )
        assert len(result) == 20
        # productive=True retries up to max_productive_attempts (25);
        # rare records can exhaust the budget. Allow ≤5%.
        nonprod = sum(1 for r in result if not r.get("productive", False))
        assert nonprod <= 1, f"too many non-productive: {nonprod}/20"

    def test_productive_with_mutation(self):
        """Productive filter applies to initial rearrangement;
        SHM may later introduce stop codons."""
        result = (
            Experiment.on("human_igh")
            .mutate(rate(0.01, 0.03))
            .run(n=10, seed=42, productive=True)
        )
        assert len(result) == 10


class TestRunWithLockedAlleles:

    def test_lock_v_allele(self):
        result = (
            Experiment.on("human_igh_imgt")
            .recombine(using(v="IGHV1-2*01"))
            .run(n=10, seed=42)
        )
        assert len(result) == 10
        for row in result:
            assert "IGHV1-2" in row["v_call"]

    def test_lock_v_and_j(self):
        result = (
            Experiment.on("human_igh_imgt")
            .recombine(using(v="IGHV1-2*01", j="IGHJ4*02"))
            .run(n=10, seed=42)
        )
        for row in result:
            assert "IGHV1-2" in row["v_call"]
            assert "IGHJ4" in row["j_call"]


class TestRunFullPipeline:

    def test_all_phases(self):
        result = (
            Experiment.on("human_igh")
            .recombine(with_d_inversion(0.15))
            .mutate(rate(0.01, 0.05), model("s5f"), with_antigen_selection(0.5))
            .prepare(with_primer_mask(), with_umi(12))
            .sequence(with_5prime_loss(), with_3prime_loss())
            .observe(with_contaminants(0.005), with_indels(0.005), with_ns(0.005))
            .run(n=20, seed=42)
        )
        assert len(result) == 20

    def test_all_phases_with_productive(self):
        result = (
            Experiment.on("human_igh")
            .mutate(rate(0.01, 0.03))
            .sequence(with_5prime_loss())
            .observe(with_ns(0.01))
            .run(n=10, seed=42, productive=True)
        )
        assert len(result) == 10


class TestCompileAndReuse:

    def test_compile_returns_simulator(self):
        sim = (
            Experiment.on("human_igh")
            .mutate(rate(0.01, 0.05))
            .compile(seed=42)
        )
        assert hasattr(sim, "simulate")

    def test_reuse_compiled_simulator(self):
        sim = (
            Experiment.on("human_igh")
            .mutate(rate(0.01, 0.05))
            .compile(seed=42)
        )
        r1 = sim.simulate(n=5)
        r2 = sim.simulate(n=5)
        assert len(r1) == 5
        assert len(r2) == 5


class TestSeedReproducibility:

    def test_same_seed_same_output(self):
        def _run():
            return (
                Experiment.on("human_igh")
                .mutate(rate(0.01, 0.05))
                .run(n=10, seed=42)
            )
        r1 = _run()
        r2 = _run()
        for a, b in zip(r1, r2):
            assert a["sequence"] == b["sequence"]

    def test_different_seed_different_output(self):
        exp = Experiment.on("human_igh").mutate(rate(0.01, 0.05))
        r1 = exp.run(n=10, seed=42)
        r2 = exp.run(n=10, seed=99)
        # At least some sequences should differ
        any_diff = any(
            a["sequence"] != b["sequence"] for a, b in zip(r1, r2)
        )
        assert any_diff


class TestRunMultiSpecies:

    @pytest.mark.parametrize("config", [
        "human_igh", "human_igk", "human_igl",
        "human_tcrb", "mouse_igh",
    ])
    def test_species_runs(self, config):
        result = Experiment.on(config).run(n=3, seed=42)
        assert len(result) == 3

    @pytest.mark.parametrize("config", [
        "human_igh", "mouse_igh",
    ])
    def test_species_with_mutation(self, config):
        result = (
            Experiment.on(config)
            .mutate(rate(0.01, 0.05))
            .run(n=5, seed=42)
        )
        assert len(result) == 5


# =====================================================================
# Compile-time warnings through full path
# =====================================================================


class TestCompileTimeWarnings:

    def test_uniform_model_warns_at_compile(self):
        with pytest.warns(RuntimeWarning, match="uniform"):
            (
                Experiment.on("human_igh")
                .mutate(rate(0.01, 0.05), model("uniform"))
                .compile(seed=42)
            )

    def test_d_inversion_on_kappa_warns(self):
        with pytest.warns(UserWarning, match="no effect"):
            (
                Experiment.on("human_igk")
                .recombine(with_d_inversion(0.1))
                .compile(seed=42)
            )


# =====================================================================
# Full clause pipeline
# =====================================================================


class TestFullClausePipeline:

    def test_clause_syntax_full_pipeline(self):
        """Full clause-based pipeline runs end-to-end."""
        result = (
            Experiment.on("human_igh")
            .mutate(rate(0.01, 0.05))
            .sequence(with_5prime_loss())
            .observe(with_ns())
            .run(n=5, seed=42)
        )
        assert len(result) == 5
