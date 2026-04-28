"""T2-6: partial-write detection across the C → Python boundary.

Before T2-6:
  * ``airr_write_tsv_row`` was ``void`` — write failures (disk full,
    broken pipe) were tracked internally via ``TsvWriter.error`` but
    never returned to the caller.
  * ``genairr_simulate`` / ``genairr_simulate_to_fd`` ignored the C
    writer entirely and unconditionally returned ``n`` after the loop.
  * Python ``simulate_to_file`` only checked ``rc < 0``, which could
    never trigger from the loop body.

After T2-6:
  * ``airr_write_tsv_row`` returns ``int`` (``0`` ok, ``-1`` fail).
  * The simulate loops break on the first row failure and return the
    actual count written.
  * ``fflush`` runs before ``fclose`` so block-buffered tail rows can't
    silently disappear.
  * Python raises ``RuntimeError`` with the partial count whenever
    ``rc != n``.

These tests probe both layers:
  * **Mock layer** — patch the C return to fake a partial. Cross-platform.
  * **/dev/full integration** — Linux-only; exercises the real C-side
    write-failure path end-to-end.
"""
from __future__ import annotations

import os
import sys
from unittest.mock import patch

import pytest

from GenAIRR import Experiment


# ── Mock layer: Python check fires when C reports rc != n ───────────


class TestMockedPartialWrite:
    """Patches ``_lib.genairr_simulate`` to return a short count and
    confirms the Python wrapper raises a clear ``RuntimeError``. Does
    not require any specific filesystem or platform."""

    def test_simulate_to_file_raises_on_short_return(self, tmp_path):
        sim = Experiment.on("human_igh").compile(seed=42)
        out = str(tmp_path / "out.tsv")
        with patch.object(sim._sim._lib, "genairr_simulate", return_value=3):
            with pytest.raises(RuntimeError) as excinfo:
                sim.simulate_to_file(10, out)
        msg = str(excinfo.value)
        assert "Partial" in msg
        assert "3" in msg and "10" in msg
        assert out in msg

    def test_simulate_to_file_raises_on_negative_return(self, tmp_path):
        """``rc < 0`` (header failed or fopen failed) → distinct
        message that does not claim a partial count."""
        sim = Experiment.on("human_igh").compile(seed=42)
        out = str(tmp_path / "out.tsv")
        with patch.object(sim._sim._lib, "genairr_simulate", return_value=-1):
            with pytest.raises(RuntimeError) as excinfo:
                sim.simulate_to_file(10, out)
        msg = str(excinfo.value)
        assert "C simulation failed" in msg
        assert "Partial" not in msg

    def test_simulate_to_file_succeeds_when_rc_equals_n(self, tmp_path):
        """Sanity: when rc == n, no exception is raised."""
        sim = Experiment.on("human_igh").compile(seed=42)
        out = str(tmp_path / "out.tsv")
        with patch.object(sim._sim._lib, "genairr_simulate", return_value=10):
            rc = sim.simulate_to_file(10, out)
        assert rc == 10

    def test_simulate_dict_path_raises_on_partial(self):
        """The dict-returning ``simulate(n)`` writes to a temp file via
        the same C path. Disk-full on /tmp would silently truncate the
        list. T2-6 closes that hole too."""
        sim = Experiment.on("human_igh").compile(seed=42)
        with patch.object(sim._sim._lib, "genairr_simulate", return_value=4):
            with pytest.raises(RuntimeError) as excinfo:
                sim._sim.simulate(10)
        msg = str(excinfo.value)
        assert "Partial" in msg
        assert "4" in msg and "10" in msg

    def test_run_to_file_propagates_partial_error(self, tmp_path):
        """End-to-end via the ``Experiment.run_to_file`` public API.
        Any caller using the user-facing entry point gets the same
        guard as direct ``simulate_to_file``."""
        out = str(tmp_path / "out.tsv")
        # Patch through the underlying CSimulator that compile() builds.
        from GenAIRR._native import _load_lib
        lib = _load_lib()
        with patch.object(lib, "genairr_simulate", return_value=2):
            with pytest.raises(RuntimeError, match=r"Partial"):
                Experiment.on("human_igh").run_to_file(
                    n=5, output_path=out, seed=42)


# ── /dev/full integration: real C-side error path ──────────────────


@pytest.mark.skipif(
    sys.platform != "linux" or not os.path.exists("/dev/full"),
    reason="/dev/full only exists on Linux",
)
class TestDevFullIntegration:
    """Exercises the full C → Python pipeline against a write target
    that always returns ENOSPC. This is the only test that proves the
    actual TsvWriter / fflush / fclose error propagation works — mocks
    cannot."""

    def test_simulate_to_file_dev_full_raises(self):
        sim = Experiment.on("human_igh").compile(seed=42)
        with pytest.raises(RuntimeError) as excinfo:
            sim.simulate_to_file(50, "/dev/full")
        msg = str(excinfo.value)
        # Either "Partial write" (some rows buffered + flush failed) or
        # "C simulation failed" (header itself failed). Either is
        # acceptable; what matters is we raise, not silently succeed.
        assert "Partial" in msg or "C simulation failed" in msg

    def test_run_to_file_dev_full_raises(self):
        with pytest.raises(RuntimeError):
            Experiment.on("human_igh").run_to_file(
                n=50, output_path="/dev/full", seed=42)

    def test_partial_count_is_less_than_n(self):
        """When stdio block buffering absorbs a few rows before the
        flush actually hits the device, rc should be 0 ≤ rc < n —
        never ``rc == n``."""
        sim = Experiment.on("human_igh").compile(seed=42)
        # Call the raw C entry point to inspect the integer return.
        rc = sim._sim._lib.genairr_simulate(
            sim._sim._handle, 50, b"/dev/full")
        assert rc < 50, (
            f"Expected partial count or -1 from /dev/full write, "
            f"got {rc} — C-side write-failure detection is broken")


# ── Happy path regression (rc == n must still work) ────────────────


class TestHappyPathStillWorks:
    """T2-6 changed C return semantics from "always n on success" to
    "actual count of rows written." For successful runs that count
    must equal n — otherwise we'd false-positive on every batch."""

    def test_run_to_file_normal_run_returns_n(self, tmp_path):
        out = str(tmp_path / "happy.tsv")
        rc = Experiment.on("human_igh").run_to_file(
            n=15, output_path=out, seed=42)
        assert rc == 15
        with open(out) as f:
            # Header + 15 rows.
            assert sum(1 for _ in f) == 16

    def test_simulate_dicts_still_returns_n_records(self):
        recs = list(Experiment.on("human_igh").run(n=12, seed=42))
        assert len(recs) == 12
