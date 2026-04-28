"""T2-5: lifecycle tests for ``CSimulator`` and ``CompiledSimulator``.

The audit's central concern is that ``CSimulator.__del__`` could crash
during interpreter shutdown when ``self._lib``'s ctypes function
pointers are torn down before the GC runs. Plus there was no
``close()`` or context-manager protocol on ``CompiledSimulator`` — so
long-running scripts that build many simulators leaked C handles
until the next GC pass.

After T2-5:
  * ``CSimulator.__del__`` swallows all exceptions (interpreter is
    exiting; nothing to do but get out of the way).
  * Both classes expose ``close()`` (idempotent), ``closed`` (bool
    property), and ``__enter__/__exit__``.
  * ``CompiledSimulator``'s public methods raise ``RuntimeError`` if
    used after ``close()`` — clean error vs. ``AttributeError``.
"""
from __future__ import annotations

import pytest

from GenAIRR import Experiment


# ── CSimulator (low-level, ctypes-backed) ───────────────────────────


class TestCSimulatorLifecycle:
    """Direct exercise of the underlying ctypes wrapper."""

    def _new(self):
        # CompiledSimulator wraps a CSimulator at ``._sim``.
        return Experiment.on("human_igh").compile(seed=42)._sim

    def test_close_is_idempotent(self):
        c = self._new()
        c.close()
        c.close()  # must not raise
        c.close()
        assert c.closed

    def test_destroy_alias_still_works(self):
        """``destroy()`` is kept as a back-compat alias for ``close()``."""
        c = self._new()
        c.destroy()
        assert c.closed

    def test_closed_property_reflects_state(self):
        c = self._new()
        assert not c.closed
        c.close()
        assert c.closed

    def test_simulate_after_close_raises(self):
        c = self._new()
        c.close()
        with pytest.raises(RuntimeError, match="destroyed"):
            c.simulate(1)

    def test_context_manager_closes_on_exit(self):
        c = self._new()
        with c as inner:
            assert inner is c
            assert not c.closed
            recs = c.simulate(2)
            assert len(recs) == 2
        assert c.closed

    def test_context_manager_closes_on_exception(self):
        """If the body of ``with`` raises, the simulator is still closed
        and the exception propagates."""
        c = self._new()
        with pytest.raises(ValueError):
            with c:
                raise ValueError("from the body")
        assert c.closed

    def test_del_safe_when_lib_torn_down(self):
        """Simulates interpreter shutdown: ``self._lib`` is gone.
        ``__del__`` must not raise — otherwise users see a Python-emitted
        traceback at process exit."""
        c = self._new()
        # Force the shutdown-like state: nuke the lib reference.
        c._lib = None
        # __del__ -> close() -> getattr(self, '_lib', None) is None,
        # so the genairr_destroy call is skipped. Must not raise.
        c.__del__()
        # And calling again is still fine.
        c.__del__()

    def test_del_safe_when_handle_already_freed(self):
        """If ``close()`` was called explicitly, the later GC-driven
        ``__del__`` must be a no-op."""
        c = self._new()
        c.close()
        c.__del__()  # must not raise

    def test_del_swallows_arbitrary_lib_errors(self):
        """Even if the lib's ``genairr_destroy`` raised something
        unexpected (e.g. during teardown), ``__del__`` must not
        propagate."""
        c = self._new()

        class _BoomLib:
            def genairr_destroy(self, _h):
                raise RuntimeError("simulated shutdown chaos")

        c._lib = _BoomLib()
        # Must not raise. (close() WOULD raise, but __del__ wraps it.)
        c.__del__()


# ── CompiledSimulator (high-level user API) ─────────────────────────


class TestCompiledSimulatorLifecycle:

    def test_close_is_idempotent(self):
        sim = Experiment.on("human_igh").compile(seed=42)
        sim.close()
        sim.close()
        sim.close()
        assert sim.closed

    def test_close_releases_underlying_csimulator(self):
        sim = Experiment.on("human_igh").compile(seed=42)
        inner = sim._sim
        assert not inner.closed
        sim.close()
        assert inner.closed
        assert sim._sim is None

    def test_simulate_after_close_raises(self):
        sim = Experiment.on("human_igh").compile(seed=42)
        sim.close()
        with pytest.raises(RuntimeError, match="closed"):
            sim.simulate(1)

    def test_set_seed_after_close_raises(self):
        sim = Experiment.on("human_igh").compile(seed=42)
        sim.close()
        with pytest.raises(RuntimeError, match="closed"):
            sim.set_seed(123)

    def test_simulate_to_file_after_close_raises(self, tmp_path):
        sim = Experiment.on("human_igh").compile(seed=42)
        sim.close()
        with pytest.raises(RuntimeError, match="closed"):
            sim.simulate_to_file(1, str(tmp_path / "out.tsv"))

    def test_stream_after_close_raises(self):
        sim = Experiment.on("human_igh").compile(seed=42)
        sim.close()
        with pytest.raises(RuntimeError, match="closed"):
            sim.stream()

    def test_iter_after_close_raises(self):
        sim = Experiment.on("human_igh").compile(seed=42)
        sim.close()
        with pytest.raises(RuntimeError, match="closed"):
            iter(sim)

    def test_context_manager_basic(self):
        with Experiment.on("human_igh").compile(seed=42) as sim:
            recs = sim.simulate(3)
            assert len(recs) == 3
            assert not sim.closed
        # After the block, sim is closed.
        assert sim.closed

    def test_context_manager_closes_on_exception(self):
        sim_ref = []
        with pytest.raises(ZeroDivisionError):
            with Experiment.on("human_igh").compile(seed=42) as sim:
                sim_ref.append(sim)
                _ = sim.simulate(1)
                _ = 1 / 0
        # Simulator still got closed despite the exception.
        assert sim_ref[0].closed

    def test_repr_reflects_closed_state(self):
        sim = Experiment.on("human_igh").compile(seed=42)
        assert "closed" not in repr(sim)
        sim.close()
        assert "closed" in repr(sim)

    def test_many_create_and_close_does_not_leak(self):
        """Simulates a parameter-sweep style loop: create N simulators
        with explicit close() — none should accumulate handles. (We
        can't directly assert non-leak, but we can confirm 50 cycles
        complete without OOM/crash and each close() finalizes
        deterministically.)"""
        for _ in range(50):
            sim = Experiment.on("human_igh").compile(seed=42)
            sim.close()
            assert sim.closed
