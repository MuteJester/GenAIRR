"""Spec tests for the SHM-metadata extension to
``DataConfig.cartridge_manifest()``.

The slice adds a ``models.shm`` block surfacing the SHM model
inventory: available models, bundled S5F kernel names, the DSL
default kernel, an optional digest of the default kernel's
bundled bytes, and the explicit ``in_content_hash=False`` v1
boundary the audit pinned.

Spec coverage (from the user brief):

1. Manifest contains ``models["shm"]``.
2. It lists ``uniform`` and ``s5f``.
3. It lists the bundled S5F kernels discovered by
   ``_s5f_loader.py``.
4. It explicitly says ``in_content_hash is False``.
5. ``json.dumps(manifest)`` still works.
6. Manifest call does not mutate ``compute_checksum``.
7. Companion contract pins flip (covered in the SHM model
   contract file).
8. Content-hash absence pin stays as-is.

Out of scope: no Rust changes, no simulation behaviour change,
no hash semantics change.
"""
from __future__ import annotations

import json

import GenAIRR as ga
from GenAIRR._s5f_loader import (
    DEFAULT_S5F_KERNEL,
    _BUILTIN_S5F_MODELS,
    available_s5f_kernels,
    builtin_s5f_kernel_digest,
)


def _cfg():
    return ga.HUMAN_IGH_OGRDB


# ──────────────────────────────────────────────────────────────────
# Spec 1 — manifest contains models["shm"]
# ──────────────────────────────────────────────────────────────────


def test_manifest_models_carries_shm_block() -> None:
    """``cartridge_manifest()["models"]["shm"]`` is present."""
    m = _cfg().cartridge_manifest()
    assert "shm" in m["models"]
    assert isinstance(m["models"]["shm"], dict)


def test_shm_block_has_documented_keys() -> None:
    """The ``shm`` block carries exactly the nine documented keys
    (five from the manifest slice + ``segment_rate_support`` from
    the per-segment SHM rates slice + ``v_subregion_support`` from
    the V-subregion cartridge annotation slice +
    ``v_subregion_rate_support`` from the V-subregion SHM rate
    slice + ``v_subregion_counter_support`` from the V-subregion
    mutation counters slice)."""
    shm = _cfg().cartridge_manifest()["models"]["shm"]
    assert set(shm.keys()) == {
        "available_models",
        "s5f_kernels_available",
        "default_s5f_kernel",
        "s5f_kernel_digest",
        "in_content_hash",
        "segment_rate_support",
        "v_subregion_support",
        "v_subregion_rate_support",
        "v_subregion_counter_support",
    }


# ──────────────────────────────────────────────────────────────────
# Spec 2 — lists uniform AND s5f
# ──────────────────────────────────────────────────────────────────


def test_shm_block_lists_uniform_and_s5f() -> None:
    """``available_models`` enumerates the two SHM model families
    the engine ships."""
    shm = _cfg().cartridge_manifest()["models"]["shm"]
    assert shm["available_models"] == ["uniform", "s5f"]


# ──────────────────────────────────────────────────────────────────
# Spec 3 — lists bundled S5F kernels
# ──────────────────────────────────────────────────────────────────


def test_shm_block_lists_bundled_s5f_kernels() -> None:
    """``s5f_kernels_available`` is the sorted list of bundled
    kernel short names, exactly matching the loader's inventory.
    Two surfaces with the same source of truth: the loader's
    ``_BUILTIN_S5F_MODELS`` and the manifest's exposed list."""
    shm = _cfg().cartridge_manifest()["models"]["shm"]
    assert shm["s5f_kernels_available"] == sorted(_BUILTIN_S5F_MODELS)


def test_shm_block_default_s5f_kernel() -> None:
    """``default_s5f_kernel`` mirrors the loader's
    ``DEFAULT_S5F_KERNEL`` constant (which the DSL's
    ``Experiment.mutate(model="s5f", s5f_model=...)`` default
    follows)."""
    shm = _cfg().cartridge_manifest()["models"]["shm"]
    assert shm["default_s5f_kernel"] == DEFAULT_S5F_KERNEL
    # And the default name appears in the kernels list.
    assert DEFAULT_S5F_KERNEL in shm["s5f_kernels_available"]


def test_shm_block_default_kernel_digest_is_stable_sha256() -> None:
    """``s5f_kernel_digest`` is a SHA-256 of the default kernel's
    bundled bytes. Two manifest calls produce the same digest;
    the digest matches the loader's independent helper."""
    shm = _cfg().cartridge_manifest()["models"]["shm"]
    digest = shm["s5f_kernel_digest"]
    assert digest is not None, (
        "default S5F kernel digest is None; the bundled .pkl bytes "
        "couldn't be read."
    )
    assert isinstance(digest, str)
    assert digest.startswith("sha256:")
    # 64 hex chars after the prefix.
    assert len(digest.split(":", 1)[1]) == 64
    # Matches the loader's helper directly.
    assert digest == builtin_s5f_kernel_digest(DEFAULT_S5F_KERNEL)


def test_shm_block_digest_stable_across_calls() -> None:
    """Idempotent: two ``cartridge_manifest()`` calls produce the
    same digest. The loader caches via ``lru_cache`` and the
    bytes are static, so this is a sanity pin against accidental
    randomness."""
    a = _cfg().cartridge_manifest()["models"]["shm"]["s5f_kernel_digest"]
    b = _cfg().cartridge_manifest()["models"]["shm"]["s5f_kernel_digest"]
    assert a == b


def test_s5f_kernel_digest_differs_across_kernels() -> None:
    """Different bundled kernels produce different digests. Pins
    that the digest is doing real work — a refactor that
    accidentally hashed the wrong file would surface here."""
    digests = {
        name: builtin_s5f_kernel_digest(name)
        for name in available_s5f_kernels()
    }
    # All non-None.
    assert all(d is not None for d in digests.values())
    # All distinct.
    assert len(set(digests.values())) == len(digests)


def test_s5f_kernel_digest_unknown_name_returns_none() -> None:
    """``builtin_s5f_kernel_digest`` returns ``None`` for an
    unknown kernel name rather than raising. Pinned so the
    manifest path stays safe under typos."""
    assert builtin_s5f_kernel_digest("totally_made_up") is None


# ──────────────────────────────────────────────────────────────────
# Spec 4 — explicitly says in_content_hash is False
# ──────────────────────────────────────────────────────────────────


def test_shm_block_in_content_hash_is_false() -> None:
    """The v1 boundary the SHM audit pinned: kernel choice does
    NOT enter ``content_hash``. The manifest documents this
    explicitly so users don't infer otherwise."""
    shm = _cfg().cartridge_manifest()["models"]["shm"]
    assert shm["in_content_hash"] is False


# ──────────────────────────────────────────────────────────────────
# Spec 5 — JSON-serialisable
# ──────────────────────────────────────────────────────────────────


def test_manifest_with_shm_block_is_json_serialisable() -> None:
    """The whole manifest round-trips through ``json.dumps`` /
    ``json.loads`` after the SHM block is added."""
    m = _cfg().cartridge_manifest()
    blob = json.dumps(m)
    parsed = json.loads(blob)
    assert parsed == m
    # And the shm sub-block survives the round-trip unchanged.
    assert parsed["models"]["shm"] == m["models"]["shm"]


# ──────────────────────────────────────────────────────────────────
# Spec 6 — manifest doesn't mutate compute_checksum
# ──────────────────────────────────────────────────────────────────


def test_manifest_with_shm_block_does_not_mutate_checksum() -> None:
    """The Slice 1 idempotence pin extends to the SHM block: a
    manifest call doesn't mutate the cartridge, so
    ``compute_checksum`` stays stable across repeated calls that
    now include the SHM extension."""
    cfg = _cfg()
    before = cfg.compute_checksum()
    cfg.cartridge_manifest()
    cfg.cartridge_manifest()
    after = cfg.compute_checksum()
    assert before == after


def test_manifest_idempotent_with_shm_block() -> None:
    """Two manifest calls produce equal dicts — same SHM block,
    same digest, same kernel list."""
    a = _cfg().cartridge_manifest()
    b = _cfg().cartridge_manifest()
    assert a == b
    assert a["models"]["shm"] == b["models"]["shm"]


# ──────────────────────────────────────────────────────────────────
# Spec 7 / 8 — content_hash NOT affected by S5F kernel
# ──────────────────────────────────────────────────────────────────


def test_content_hash_does_not_depend_on_s5f_kernel() -> None:
    """Sanity pin alongside the manifest's
    ``in_content_hash=False``: an actual ``content_hash``
    comparison across two experiments using different S5F kernels
    on the SAME cartridge produces the same hash. The manifest
    correctly reports the v1 boundary."""
    h_default = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", count=5, s5f_model="hh_s5f")
        .compile()
        .refdata
        .content_hash()
    )
    h_other = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", count=5, s5f_model="hkl_s5f")
        .compile()
        .refdata
        .content_hash()
    )
    assert h_default == h_other
