"""Functional-status curation v1.

Pins the Python-facing slice that adds:

- a per-allele ``functional_status`` field on the engine-native
  ``RefDataConfig`` (``"functional"`` / ``"orf"`` / ``"pseudogene"`` /
  ``"unknown"`` strings, case-insensitive at the boundary, ``None`` when
  the cartridge did not annotate the allele);
- a new ``"functional_status"`` curation policy with ``allowed=`` and
  ``keep_unannotated=`` kwargs;
- preservation of existing bundled cartridges (status is uniformly
  ``None`` today; the default ``keep_unannotated=True`` keeps them
  intact).

No data migration in this slice. The bundled ``.pkl`` builders are NOT
rebuilt — those cartridges still carry zero status annotations.
"""
from __future__ import annotations

import pytest

from GenAIRR import _engine
from GenAIRR._refdata_resolver import _resolve_config_name


# ──────────────────────────────────────────────────────────────────
# 1. Engine-native: allele can store and read functional status
# ──────────────────────────────────────────────────────────────────


def test_engine_native_allele_stores_status_round_trip() -> None:
    cfg = _engine.RefDataConfig.vdj()
    cfg.add_v_allele(
        "MYV*01", "MYV", b"TGTAAACCC", anchor=0, functional_status="functional"
    )
    cfg.add_v_allele(
        "MYV*02", "MYV", b"TAGAAACCC", anchor=0, functional_status="orf"
    )
    cfg.add_v_allele(
        "MYV*03", "MYV", b"TAAAAACCC", anchor=0, functional_status="pseudogene"
    )
    cfg.add_v_allele(
        "MYV*04", "MYV", b"TGTAAACCC", anchor=0, functional_status="unknown"
    )
    cfg.add_v_allele("MYV*05", "MYV", b"TGTAAACCC", anchor=0)

    statuses = [cfg.v_allele(i).functional_status for i in range(5)]
    assert statuses == ["functional", "orf", "pseudogene", "unknown", None]


def test_engine_native_allele_status_accepts_case_and_alias() -> None:
    """``"F"``, ``"P"``, mixed case all normalise to canonical lowercase."""
    cfg = _engine.RefDataConfig.vj()
    cfg.add_v_allele("a*01", "a", b"TGTAAACCC", anchor=0, functional_status="F")
    cfg.add_v_allele("a*02", "a", b"TGTAAACCC", anchor=0, functional_status="ORF")
    cfg.add_v_allele("a*03", "a", b"TGTAAACCC", anchor=0, functional_status="P")
    cfg.add_v_allele(
        "a*04", "a", b"TGTAAACCC", anchor=0, functional_status="UnKnOwN"
    )
    assert cfg.v_allele(0).functional_status == "functional"
    assert cfg.v_allele(1).functional_status == "orf"
    assert cfg.v_allele(2).functional_status == "pseudogene"
    assert cfg.v_allele(3).functional_status == "unknown"


def test_engine_native_allele_status_absent_is_none() -> None:
    cfg = _engine.RefDataConfig.vj()
    cfg.add_v_allele("a*01", "a", b"TGTAAACCC", anchor=0)
    assert cfg.v_allele(0).functional_status is None


def test_engine_native_invalid_status_raises_value_error() -> None:
    cfg = _engine.RefDataConfig.vj()
    with pytest.raises(ValueError, match="unknown functional_status"):
        cfg.add_v_allele(
            "a*01", "a", b"TGTAAACCC", anchor=0, functional_status="not-a-status"
        )


# ──────────────────────────────────────────────────────────────────
# 2. Python bridge: DataConfig allele.functional_status → Rust
# ──────────────────────────────────────────────────────────────────


def test_bundled_human_igh_round_trips_status_as_none() -> None:
    """Bundled human_igh has no functional_status annotations
    today. After bridging, every Rust allele reports ``None``."""
    from GenAIRR._refdata_resolver import dataconfig_to_refdata

    cfg = _resolve_config_name("human_igh")
    refdata = dataconfig_to_refdata(cfg)
    for i in range(refdata.v_pool_size()):
        assert refdata.v_allele(i).functional_status is None
    for i in range(refdata.j_pool_size()):
        assert refdata.j_allele(i).functional_status is None


def test_bridge_translates_native_enum_via_name_attribute() -> None:
    """Python ``Allele`` may carry a native enum from ``_native._anchor``.
    The bridge reads its ``.name`` attribute and normalises it."""
    from GenAIRR._refdata_resolver import _normalise_functional_status

    class _FakeEnum:
        name = "FUNCTIONAL"

    assert _normalise_functional_status(_FakeEnum()) == "functional"
    assert _normalise_functional_status("ORF") == "orf"
    assert _normalise_functional_status("p") == "pseudogene"
    assert _normalise_functional_status(None) is None
    # Unknown strings silently collapse to None — bridge stays quiet
    # on legacy data with surprise labels.
    assert _normalise_functional_status("FNG") is None


# ──────────────────────────────────────────────────────────────────
# 3. Curation policy: functional_status filtering
# ──────────────────────────────────────────────────────────────────


def _build_mixed_cfg() -> "_engine.RefDataConfig":
    """Build a tiny VDJ catalogue with V/D/J pools holding one of every
    status (Functional, ORF, Pseudogene, Unknown, None).
    """
    cfg = _engine.RefDataConfig.vdj()
    seq = b"TGTAAACCC"
    cfg.add_v_allele("v-f*01", "v-f", seq, anchor=0, functional_status="functional")
    cfg.add_v_allele("v-o*01", "v-o", seq, anchor=0, functional_status="orf")
    cfg.add_v_allele("v-p*01", "v-p", seq, anchor=0, functional_status="pseudogene")
    cfg.add_v_allele("v-u*01", "v-u", seq, anchor=0, functional_status="unknown")
    cfg.add_v_allele("v-na*01", "v-na", seq, anchor=0)  # unannotated

    cfg.add_d_allele("d-f*01", "d-f", seq, functional_status="functional")
    cfg.add_d_allele("d-na*01", "d-na", seq)

    cfg.add_j_allele("j-f*01", "j-f", b"TGGAAACCC", anchor=0, functional_status="functional")
    cfg.add_j_allele("j-p*01", "j-p", b"TGGAAACCC", anchor=0, functional_status="pseudogene")
    return cfg


def test_curation_functional_only_drops_non_functional_and_unannotated() -> None:
    cfg = _build_mixed_cfg()
    curated = cfg.curated("functional_status", allowed=["functional"], keep_unannotated=False)
    v_names = [curated.v_allele(i).name for i in range(curated.v_pool_size())]
    d_names = [curated.d_allele(i).name for i in range(curated.d_pool_size())]
    j_names = [curated.j_allele(i).name for i in range(curated.j_pool_size())]
    assert v_names == ["v-f*01"]
    assert d_names == ["d-f*01"]
    assert j_names == ["j-f*01"]


def test_curation_functional_only_keeps_unannotated_when_flagged() -> None:
    cfg = _build_mixed_cfg()
    curated = cfg.curated("functional_status", allowed=["functional"], keep_unannotated=True)
    v_names = [curated.v_allele(i).name for i in range(curated.v_pool_size())]
    d_names = [curated.d_allele(i).name for i in range(curated.d_pool_size())]
    assert "v-f*01" in v_names
    assert "v-na*01" in v_names  # unannotated kept
    assert "v-o*01" not in v_names
    assert "v-p*01" not in v_names
    assert "v-u*01" not in v_names
    assert "d-f*01" in d_names
    assert "d-na*01" in d_names


def test_curation_functional_plus_orf_allows_both() -> None:
    cfg = _build_mixed_cfg()
    curated = cfg.curated(
        "functional_status",
        allowed=["functional", "orf"],
        keep_unannotated=False,
    )
    v_names = [curated.v_allele(i).name for i in range(curated.v_pool_size())]
    assert v_names == ["v-f*01", "v-o*01"]


def test_curation_default_allowed_is_functional() -> None:
    """Calling without ``allowed=`` keeps ``"functional"`` only."""
    cfg = _build_mixed_cfg()
    curated = cfg.curated("functional_status", keep_unannotated=False)
    v_names = [curated.v_allele(i).name for i in range(curated.v_pool_size())]
    assert v_names == ["v-f*01"]


def test_curation_default_keep_unannotated_is_true() -> None:
    """Calling without ``keep_unannotated=`` defaults to True so bundled
    cartridges (status=None everywhere) survive the filter."""
    cfg = _build_mixed_cfg()
    curated = cfg.curated("functional_status", allowed=["functional"])
    v_names = [curated.v_allele(i).name for i in range(curated.v_pool_size())]
    assert "v-na*01" in v_names


def test_curation_status_tags_identity_source() -> None:
    cfg = _build_mixed_cfg()
    curated = cfg.curated(
        "functional_status",
        allowed=["functional", "orf"],
        keep_unannotated=False,
    )
    src = curated.identity()["source"]
    # The tag carries the canonical allowed list (sorted) and the flag.
    assert src is not None
    assert "curated:functional_status:functional,orf|keep_unannotated=false" in src


def test_curation_status_does_not_touch_raw_or_anchors_policies() -> None:
    """Raw + functional_anchors_only must reject the ``allowed=`` kwarg —
    it would be silently ignored otherwise, masking a user bug."""
    cfg = _build_mixed_cfg()
    with pytest.raises(ValueError, match="does not accept an 'allowed'"):
        cfg.curated("raw", allowed=["functional"])
    with pytest.raises(ValueError, match="does not accept an 'allowed'"):
        cfg.curated("functional_anchors_only", allowed=["functional"])


def test_curation_status_invalid_string_raises() -> None:
    cfg = _build_mixed_cfg()
    with pytest.raises(ValueError, match="unknown functional_status"):
        cfg.curated("functional_status", allowed=["bogus"])


# ──────────────────────────────────────────────────────────────────
# 4. Content hash includes functional status
# ──────────────────────────────────────────────────────────────────


def test_content_hash_changes_when_only_status_changes() -> None:
    """Two cartridges identical except for ``functional_status`` produce
    different content hashes."""
    a = _engine.RefDataConfig.vj()
    a.add_v_allele("v*01", "v", b"TGTAAACCC", anchor=0, functional_status="functional")
    a.add_j_allele("j*01", "j", b"TGGAAACCC", anchor=0, functional_status="functional")

    b = _engine.RefDataConfig.vj()
    b.add_v_allele("v*01", "v", b"TGTAAACCC", anchor=0, functional_status="pseudogene")
    b.add_j_allele("j*01", "j", b"TGGAAACCC", anchor=0, functional_status="functional")

    assert a.content_hash() != b.content_hash()


def test_content_hash_stable_when_status_absent_on_both() -> None:
    a = _engine.RefDataConfig.vj()
    a.add_v_allele("v*01", "v", b"TGTAAACCC", anchor=0)
    a.add_j_allele("j*01", "j", b"TGGAAACCC", anchor=0)

    b = _engine.RefDataConfig.vj()
    b.add_v_allele("v*01", "v", b"TGTAAACCC", anchor=0)
    b.add_j_allele("j*01", "j", b"TGGAAACCC", anchor=0)
    assert a.content_hash() == b.content_hash()


def test_none_distinct_from_explicit_unknown_in_hash() -> None:
    """``None`` (no annotation) and ``Some(Unknown)`` are different
    cartridge states — the hash reflects the distinction."""
    a = _engine.RefDataConfig.vj()
    a.add_v_allele("v*01", "v", b"TGTAAACCC", anchor=0)
    a.add_j_allele("j*01", "j", b"TGGAAACCC", anchor=0)

    b = _engine.RefDataConfig.vj()
    b.add_v_allele("v*01", "v", b"TGTAAACCC", anchor=0, functional_status="unknown")
    b.add_j_allele("j*01", "j", b"TGGAAACCC", anchor=0)
    assert a.content_hash() != b.content_hash()


# ──────────────────────────────────────────────────────────────────
# 5. Empty-pool compile failure
# ──────────────────────────────────────────────────────────────────


def test_status_curation_empty_v_pool_surfaces_empty_required_pool() -> None:
    """If status curation drops every V allele, ``validate()`` reports
    ``EmptyRequiredPool`` — the same diagnostic compile() will trip on."""
    cfg = _engine.RefDataConfig.vj()
    cfg.add_v_allele("v*01", "v", b"TGTAAACCC", anchor=0, functional_status="pseudogene")
    cfg.add_j_allele("j*01", "j", b"TGGAAACCC", anchor=0, functional_status="functional")
    curated = cfg.curated(
        "functional_status",
        allowed=["functional"],
        keep_unannotated=False,
    )
    issues = curated.validate()
    kinds = [(i["kind"], i.get("segment")) for i in issues]
    assert ("EmptyRequiredPool", "V") in kinds


# ──────────────────────────────────────────────────────────────────
# 6. Bundled cartridge invariance under default policy
# ──────────────────────────────────────────────────────────────────


def test_bundled_cartridge_unchanged_under_default_status_policy() -> None:
    """Bundled human_igh leaves status unannotated. With the default
    ``keep_unannotated=True``, curating with status policy is a no-op
    on the pool sizes."""
    from GenAIRR._refdata_resolver import dataconfig_to_refdata

    cfg = _resolve_config_name("human_igh")
    refdata = dataconfig_to_refdata(cfg)
    sizes_before = (
        refdata.v_pool_size(),
        refdata.d_pool_size(),
        refdata.j_pool_size(),
    )
    curated = refdata.curated("functional_status")
    sizes_after = (
        curated.v_pool_size(),
        curated.d_pool_size(),
        curated.j_pool_size(),
    )
    assert sizes_before == sizes_after


def test_bundled_cartridge_emptied_when_keep_unannotated_false() -> None:
    """Conversely, the same bundled cartridge becomes empty when
    ``keep_unannotated=False``, because every allele is unannotated.
    This is the diagnostic the user would see if they ran a
    'strict-status' policy without first rebuilding the bundled data."""
    from GenAIRR._refdata_resolver import dataconfig_to_refdata

    cfg = _resolve_config_name("human_igh")
    refdata = dataconfig_to_refdata(cfg)
    curated = refdata.curated(
        "functional_status",
        allowed=["functional"],
        keep_unannotated=False,
    )
    assert curated.v_pool_size() == 0
    assert curated.j_pool_size() == 0


# ──────────────────────────────────────────────────────────────────
# 7. Experiment.curate_refdata("functional_status", ...) integration
# ──────────────────────────────────────────────────────────────────


def test_experiment_curate_refdata_accepts_status_kwargs() -> None:
    from GenAIRR.experiment import Experiment

    refdata = _build_mixed_cfg()
    exp = Experiment.on(refdata)
    exp.curate_refdata(
        "functional_status",
        allowed=["functional"],
        keep_unannotated=False,
    )
    # After curation V pool drops to one functional allele.
    curated = exp._refdata
    assert curated.v_pool_size() == 1
    assert curated.v_allele(0).name == "v-f*01"


def test_experiment_curate_refdata_chainable() -> None:
    """``curate_refdata`` returns self for fluent chaining."""
    from GenAIRR.experiment import Experiment

    refdata = _build_mixed_cfg()
    out = Experiment.on(refdata).curate_refdata(
        "functional_status",
        allowed=["functional", "orf"],
    )
    assert isinstance(out, Experiment)
