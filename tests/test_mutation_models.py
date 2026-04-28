"""T2-4: end-to-end tests for the mutation-model selector.

Before T2-4, ``model("uniform")`` emitted a ``RuntimeWarning`` then
silently ran S5F anyway. The C uniform kernel existed (and was
unit-tested in `test_ops.c::test_uniform_mutate_*`) but was not
reachable from the Python API.

After T2-4:
  * ``CSimulator.set_mutation_model("s5f"|"uniform")`` is the new
    ctypes binding.
  * ``Mutate`` / ``SimulateCSR`` Step dataclasses gain a ``model``
    field; ``configure()`` propagates it to the C side.
  * ``model("uniform")`` is no longer a no-op — it actually runs the
    uniform kernel, which produces statistically different sequences
    than S5F at the same rate.
"""
from __future__ import annotations

import warnings

import pytest

from GenAIRR import Experiment
from GenAIRR.ops import model, rate, with_isotype_rates


# ── CSimulator.set_mutation_model contract ──────────────────────────


class TestSetMutationModelBinding:
    """Direct exercise of the ctypes binding, independent of the DSL."""

    def _make_sim(self):
        # Use compile() to build a CSimulator without running any sequences.
        return Experiment.on("human_igh").compile(seed=42)

    def test_default_is_s5f(self):
        """A freshly compiled simulator runs S5F by default — verified
        indirectly via simulation determinism vs an explicit s5f setter."""
        sim_a = self._make_sim()
        sim_b = self._make_sim()
        sim_b._sim.set_mutation_model("s5f")
        # Both should produce the same sequences for the same seed.
        recs_a = list(sim_a.simulate(n=5))
        recs_b = list(sim_b.simulate(n=5))
        for a, b in zip(recs_a, recs_b):
            assert a["sequence"] == b["sequence"]

    def test_unknown_name_raises(self):
        sim = self._make_sim()
        with pytest.raises(ValueError, match="Unknown mutation model"):
            sim._sim.set_mutation_model("nonsense")

    def test_set_uniform_then_back_to_s5f(self):
        """The setter is idempotent and reversible — flipping uniform
        and back to s5f should not corrupt simulator state."""
        sim = self._make_sim()
        sim._sim.set_mutation_model("uniform")
        sim._sim.set_mutation_model("s5f")
        # Should still produce valid sequences (non-empty).
        recs = list(sim.simulate(n=2))
        assert len(recs) == 2
        for r in recs:
            assert r["sequence"]


# ── DSL: model() now propagates to the C kernel ─────────────────────


class TestDslModelPropagation:

    def test_uniform_compiles_without_warning(self):
        """T2-4: the silent-fallback RuntimeWarning is gone."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            (Experiment.on("human_igh")
             .mutate(rate(0.01, 0.05), model("uniform"))
             .compile(seed=42))
        rt = [x for x in w if issubclass(x.category, RuntimeWarning)]
        assert rt == [], \
            f"Unexpected RuntimeWarning(s): {[str(x.message) for x in rt]}"

    def test_uniform_produces_mutations(self):
        """``model("uniform")`` is no longer a no-op — at a non-trivial
        rate, every sequence accumulates mutations."""
        recs = list(
            Experiment.on("human_igh")
            .mutate(rate(0.05, 0.10), model("uniform"))
            .run(n=20, seed=42)
        )
        n_with_mut = sum(1 for r in recs if r["n_mutations"] > 0)
        assert n_with_mut == 20, \
            f"Expected mutations on every record, got {n_with_mut}/20"

    def test_uniform_and_s5f_produce_different_sequences(self):
        """Same seed, same rate, different model → different sequences.
        If they were identical, uniform would be silently running S5F
        (the pre-T2-4 bug)."""
        seed = 4242
        n = 10
        s5f = list(Experiment.on("human_igh")
                   .mutate(rate(0.05, 0.10), model("s5f"))
                   .run(n=n, seed=seed))
        uni = list(Experiment.on("human_igh")
                   .mutate(rate(0.05, 0.10), model("uniform"))
                   .run(n=n, seed=seed))

        differ = sum(1 for a, b in zip(s5f, uni)
                     if a["sequence"] != b["sequence"])
        assert differ >= n - 1, (
            f"Only {differ}/{n} sequences differ between s5f and uniform — "
            f"the kernels are running the same code path (T2-4 regression)"
        )

    def test_default_model_is_s5f(self):
        """No ``model()`` clause → default is S5F. Verified by parity
        with explicit ``model('s5f')`` (same seed → same output)."""
        n, seed = 5, 7
        default = list(Experiment.on("human_igh")
                       .mutate(rate(0.02, 0.06))
                       .run(n=n, seed=seed))
        explicit = list(Experiment.on("human_igh")
                        .mutate(rate(0.02, 0.06), model("s5f"))
                        .run(n=n, seed=seed))
        for a, b in zip(default, explicit):
            assert a["sequence"] == b["sequence"]


# ── CSR + uniform: rate adjustment still applies ────────────────────


class TestCSRWithUniform:
    """T2-4: CSR rate-adjustment was previously baked into the S5F
    code path. Under uniform, rate adjustment must still apply (CSR
    tunes ``cfg->{min,max}_mutation_rate`` for the call, which the
    uniform kernel reads directly)."""

    def test_csr_plus_uniform_runs_and_mutates(self):
        recs = list(
            Experiment.on("human_igh")
            .mutate(rate(0.05, 0.10), with_isotype_rates(), model("uniform"))
            .run(n=10, seed=42)
        )
        assert len(recs) == 10
        # CSR + uniform must still produce some mutations.
        assert sum(r["n_mutations"] for r in recs) > 0


# ── Step dataclass: model field is plumbed end-to-end ───────────────


class TestStepDataclassModelField:
    def test_mutate_step_carries_model(self):
        from GenAIRR.steps import Mutate
        from GenAIRR.protocol import _clauses_to_steps
        steps = _clauses_to_steps(
            recombine_clauses=[],
            mutate_clauses=[model("uniform"), rate(0.01, 0.05)],
            prepare_clauses=[], sequence_clauses=[], observe_clauses=[],
        )
        m = [s for s in steps if isinstance(s, Mutate)][0]
        assert m.model == "uniform"

    def test_csr_step_carries_model(self):
        from GenAIRR.steps import SimulateCSR
        from GenAIRR.protocol import _clauses_to_steps
        steps = _clauses_to_steps(
            recombine_clauses=[],
            mutate_clauses=[model("uniform"), rate(0.01, 0.05),
                            with_isotype_rates()],
            prepare_clauses=[], sequence_clauses=[], observe_clauses=[],
        )
        csr = [s for s in steps if isinstance(s, SimulateCSR)][0]
        assert csr.model == "uniform"
