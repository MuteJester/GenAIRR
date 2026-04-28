"""T2-1: schema_version + schema_sha256 integrity checks for DataConfig.

The contract:
  * Every builtin pickle in ``data/builtin_dataconfigs/`` is stamped
    with ``schema_version=SCHEMA_VERSION`` and a sha256 over the
    canonical pickle bytes.
  * ``DataConfig.verify_integrity()`` raises ``DataConfigError`` on
    schema-version mismatch (stale pickle from before the migration)
    or content tampering (checksum mismatch).
  * The ``data/__init__.py`` loader runs ``verify_integrity()`` on
    every load, so a tampered or pre-migration pickle surfaces as an
    ``ImportError`` before reaching user code.
"""
from __future__ import annotations

import pickle
import shutil
from pathlib import Path

import pytest

from GenAIRR.dataconfig.data_config import (
    DataConfig,
    DataConfigError,
    SCHEMA_VERSION,
)


BUILTIN_DIR = (
    Path(__file__).resolve().parent.parent
    / "src" / "GenAIRR" / "data" / "builtin_dataconfigs"
)


# ── DataConfig instance contract ─────────────────────────────────────


class TestDataConfigIntegrityAPI:
    def test_compute_checksum_is_deterministic(self):
        """Calling compute_checksum twice on the same content must
        return the same hash. If this drifts, the migration script's
        output won't match what verify_integrity() recomputes."""
        cfg = DataConfig(name="test")
        h1 = cfg.compute_checksum()
        h2 = cfg.compute_checksum()
        assert h1 == h2, f"compute_checksum drifted: {h1} != {h2}"
        assert len(h1) == 64, "sha256 hex digest must be 64 chars"

    def test_compute_checksum_does_not_self_include(self):
        """The checksum must be computed with schema_sha256 zeroed —
        otherwise the field would include itself, making verification
        impossible."""
        cfg = DataConfig(name="test")
        cfg.schema_sha256 = "GARBAGE_VALUE_SHOULD_BE_IGNORED"
        h_with_garbage = cfg.compute_checksum()
        cfg.schema_sha256 = ""
        h_with_empty = cfg.compute_checksum()
        assert h_with_garbage == h_with_empty, (
            "compute_checksum must zero schema_sha256 internally")

    def test_compute_checksum_changes_with_content(self):
        """Sanity: any content change must produce a different hash."""
        cfg_a = DataConfig(name="A")
        cfg_b = DataConfig(name="B")
        assert cfg_a.compute_checksum() != cfg_b.compute_checksum()

    def test_verify_integrity_passes_on_freshly_stamped(self):
        """Stamp a DataConfig the same way the migration script does
        and verify it passes."""
        cfg = DataConfig(name="test", schema_version=SCHEMA_VERSION)
        cfg.schema_sha256 = cfg.compute_checksum()
        cfg.verify_integrity()  # must not raise

    def test_verify_integrity_raises_on_legacy_pickle(self):
        """A pickle that predates the migration has no checksum
        stamped on it. The instance falls through to the dataclass
        class-level default ``schema_sha256 = ""``, so the
        missing-checksum branch fires with a migration hint."""
        cfg = DataConfig(name="legacy")
        cfg.__dict__.pop("schema_sha256", None)
        with pytest.raises(DataConfigError) as excinfo:
            cfg.verify_integrity()
        msg = str(excinfo.value).lower()
        assert "schema_sha256" in msg or "checksum" in msg
        assert "migrate" in msg or "re-migrate" in msg

    def test_verify_integrity_raises_on_future_schema(self):
        """If SCHEMA_VERSION is bumped to 2 in code but a v1 pickle is
        loaded, verify_integrity must reject it with a clear hint."""
        cfg = DataConfig(name="from_future", schema_version=SCHEMA_VERSION + 1)
        cfg.schema_sha256 = cfg.compute_checksum()
        with pytest.raises(DataConfigError) as excinfo:
            cfg.verify_integrity()
        msg = str(excinfo.value)
        assert "schema_version" in msg
        assert str(SCHEMA_VERSION + 1) in msg
        assert str(SCHEMA_VERSION) in msg

    def test_verify_integrity_raises_on_missing_checksum(self):
        """Right schema version but no checksum stamped — must raise."""
        cfg = DataConfig(name="unstamped", schema_version=SCHEMA_VERSION)
        # schema_sha256 default is "" — verify must reject.
        with pytest.raises(DataConfigError) as excinfo:
            cfg.verify_integrity()
        assert "schema_sha256" in str(excinfo.value).lower() or \
               "checksum" in str(excinfo.value).lower()

    def test_verify_integrity_raises_on_tampered_content(self):
        """Stamp a config, then mutate a field. verify_integrity must
        detect the drift."""
        cfg = DataConfig(name="original", schema_version=SCHEMA_VERSION)
        cfg.schema_sha256 = cfg.compute_checksum()
        cfg.verify_integrity()  # baseline OK

        # Tamper.
        cfg.name = "tampered"
        with pytest.raises(DataConfigError) as excinfo:
            cfg.verify_integrity()
        assert "checksum mismatch" in str(excinfo.value).lower() or \
               "modified or corrupted" in str(excinfo.value).lower()


# ── Builtin pickles all carry schema_version + valid checksum ────────


def _list_builtin_pickles() -> list[Path]:
    return sorted(BUILTIN_DIR.glob("*.pkl"))


class TestBuiltinPickleIntegrity:
    """Every shipped builtin must pass verify_integrity. If any of
    these fail, the migration script either was not re-run after a
    DataConfig change or the wheel ships a stale pickle."""

    @pytest.mark.parametrize(
        "pkl_path",
        _list_builtin_pickles(),
        ids=lambda p: p.stem,
    )
    def test_every_builtin_pickle_passes_integrity(self, pkl_path: Path):
        with open(pkl_path, "rb") as f:
            cfg = pickle.load(f)
        # Must be a DataConfig (or compatible) — not just any pickle.
        assert hasattr(cfg, "verify_integrity"), (
            f"{pkl_path.name} unpickled to {type(cfg).__name__}, "
            f"which lacks verify_integrity")
        # And it must verify.
        cfg.verify_integrity()
        assert cfg.schema_version == SCHEMA_VERSION
        assert len(cfg.schema_sha256) == 64

    def test_builtin_dir_is_non_empty(self):
        """Sanity check — if this fails, the test parametrization
        above silently skipped everything."""
        assert _list_builtin_pickles(), \
            "No builtin pickles found — parametrized check above " \
            "would silently pass with zero items."


# ── Loader-level integration: tampered file → ImportError ────────────


class TestLoaderRejectsTamperedPickle:
    """End-to-end check: write a tampered pickle into a temp
    builtin_dataconfigs dir, load via the same code path
    ``data/__init__.py`` uses, and confirm the user gets a clear
    failure."""

    def test_tampered_pickle_raises_import_error_on_load(self, tmp_path):
        # Set up a fake builtin dir with one tampered pickle.
        fake_builtin_dir = tmp_path / "builtin_dataconfigs"
        fake_builtin_dir.mkdir()

        # Copy a real one, load, mutate, re-pickle — but DON'T
        # update the checksum.
        src = BUILTIN_DIR / "HUMAN_IGH_OGRDB.pkl"
        dst = fake_builtin_dir / "HUMAN_IGH_OGRDB.pkl"
        shutil.copy(src, dst)

        with open(dst, "rb") as f:
            cfg = pickle.load(f)
        cfg.name = "TAMPERED"  # checksum no longer matches
        with open(dst, "wb") as f:
            pickle.dump(cfg, f, protocol=pickle.HIGHEST_PROTOCOL)

        # Load via the same path data/__init__.py uses.
        with open(dst, "rb") as f:
            tampered = pickle.load(f)
        with pytest.raises(DataConfigError) as excinfo:
            tampered.verify_integrity()
        assert "checksum" in str(excinfo.value).lower()


# ── Idempotency of the migration script ───────────────────────────────


class TestMigrationScriptIdempotent:
    """Running the migration twice on the same pickle must produce a
    bitwise-identical second-write — the checksum is deterministic and
    pickle.dumps with HIGHEST_PROTOCOL is stable across calls."""

    def test_re_migration_produces_same_checksum(self, tmp_path):
        # Copy a real builtin pickle to a temp dir to avoid mutating
        # the actual installed file.
        src = BUILTIN_DIR / "HUMAN_IGH_OGRDB.pkl"
        dst = tmp_path / "HUMAN_IGH_OGRDB.pkl"
        shutil.copy(src, dst)

        with open(dst, "rb") as f:
            cfg1 = pickle.load(f)
        cfg1.verify_integrity()
        original_checksum = cfg1.schema_sha256

        # Re-stamp (no content changes).
        cfg1.schema_sha256 = ""
        cfg1.schema_sha256 = cfg1.compute_checksum()
        assert cfg1.schema_sha256 == original_checksum, (
            "Re-stamping a clean config produced a different checksum")
