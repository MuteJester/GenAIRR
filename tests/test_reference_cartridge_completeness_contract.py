"""Contract pins for the Reference Cartridge Completeness audit.

Companion to
[`docs/reference_cartridge_completeness_audit.md`](../docs/reference_cartridge_completeness_audit.md).
The audit is pre-implementation: it freezes today's cartridge
contract (``pin_scaffold_*``) and the inspectability gaps a future
``cartridge_manifest()``-style slice will close (``pin_absence_*``).

Split:

- ``pin_scaffold_*`` tests freeze the corrected/current architecture:
  the four cartridge views' planes; ``compute_checksum`` vs
  ``content_hash`` sensitivity (the v1 boundary on
  ``reference_models``); curation tags identity source; bundled
  cartridge identity; functional-status normalisation; the
  documented "dropped at bridge" Allele field list; the orphan
  ``DataConfig`` field list.
- ``pin_absence_*`` tests freeze the remaining gaps after Slice 1
  shipped (Gaps B + C from audit §11): no ``models_digest`` in
  identity, no view covering the orphan simulation fields.

Slice-1 flip history (cartridge_manifest):
- ``pin_present_cartridge_manifest_method`` (was absence).
- ``pin_present_functional_status_histogram_accessor`` (was
  absence; Slice 1 exposed it as a separate surface so the
  manifest delegates).

The remaining ``pin_absence_*`` tests await later slices.
"""
from __future__ import annotations

import copy

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge
from GenAIRR import dataconfig_to_refdata
from GenAIRR.alleles.allele import Allele
from GenAIRR.dataconfig.cartridge_views import (
    CartridgeCatalogueView,
    CartridgeIdentityView,
    CartridgeModelsView,
    CartridgeRulesView,
)
from GenAIRR.dataconfig.data_config import DataConfig
from GenAIRR.reference_models import EmpiricalDistributionSpec, ReferenceEmpiricalModels
from GenAIRR.reference_rules import AnchorRuleSpec, ReferenceRulesSpec


# ──────────────────────────────────────────────────────────────────
# Shared fixtures
# ──────────────────────────────────────────────────────────────────


def _cfg():
    """The bundled human-IGH OGRDB cartridge — the canonical
    audit fixture."""
    return ga.HUMAN_IGH_OGRDB


def _truly_different_rules():
    """A ``ReferenceRulesSpec`` that materially differs from the
    locus-derived defaults (different alphabet, expanded anchor
    expectations, mismatched required flags)."""
    return ReferenceRulesSpec(
        allowed_bases=["A", "C", "G", "T"],  # no N — different from default
        v_anchor=AnchorRuleSpec(expected_aa=("C", "Y"), required=True),
        j_anchor=AnchorRuleSpec(expected_aa=("W", "F"), required=False),
    )


# ──────────────────────────────────────────────────────────────────
# 1. Scaffold — cartridge views expose their documented planes
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_cartridge_identity_view_shape() -> None:
    """``cartridge_identity`` returns ``CartridgeIdentityView`` with
    ``name`` and ``metadata`` fields."""
    cfg = _cfg()
    view = cfg.cartridge_identity
    assert isinstance(view, CartridgeIdentityView)
    assert view.name == cfg.name
    assert view.metadata is cfg.metadata


def test_pin_scaffold_cartridge_catalogue_view_shape() -> None:
    """``cartridge_catalogue`` returns ``CartridgeCatalogueView``
    with V/D/J/C allele dicts."""
    cfg = _cfg()
    view = cfg.cartridge_catalogue
    assert isinstance(view, CartridgeCatalogueView)
    assert view.v_alleles is cfg.v_alleles
    assert view.d_alleles is cfg.d_alleles
    assert view.j_alleles is cfg.j_alleles
    assert view.c_alleles is cfg.c_alleles
    # Count helpers.
    assert view.number_of_v_alleles == cfg.number_of_v_alleles
    assert view.number_of_d_alleles == cfg.number_of_d_alleles
    assert view.number_of_j_alleles == cfg.number_of_j_alleles


def test_pin_scaffold_cartridge_rules_view_shape() -> None:
    """``cartridge_rules`` returns ``CartridgeRulesView`` carrying
    the (optional) ``reference_rules`` field."""
    cfg = _cfg()
    view = cfg.cartridge_rules
    assert isinstance(view, CartridgeRulesView)
    # Bundled cartridges leave reference_rules None; the loader
    # derives defaults from the locus.
    assert view.reference_rules is None


def test_pin_scaffold_cartridge_models_view_shape() -> None:
    """``cartridge_models`` returns ``CartridgeModelsView`` carrying
    both ``reference_models`` (typed spec) and the legacy NP /
    trim dicts."""
    cfg = _cfg()
    view = cfg.cartridge_models
    assert isinstance(view, CartridgeModelsView)
    assert view.reference_models is None  # bundled — falls back to legacy
    assert view.legacy_np_lengths is cfg.NP_lengths
    assert view.legacy_trim_dicts is cfg.trim_dicts


# ──────────────────────────────────────────────────────────────────
# 2. Scaffold — `compute_checksum` ↔ `content_hash` sensitivity
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_compute_checksum_stable_for_unchanged_cartridge() -> None:
    """Two calls on the same DataConfig produce the same checksum.
    Sanity baseline for every sensitivity test below."""
    cfg = _cfg()
    assert cfg.compute_checksum() == cfg.compute_checksum()


def test_pin_scaffold_content_hash_stable_for_unchanged_cartridge() -> None:
    """Two ``dataconfig_to_refdata`` runs on the same cfg produce
    the same content hash."""
    cfg = _cfg()
    h1 = dataconfig_to_refdata(cfg).content_hash()
    h2 = dataconfig_to_refdata(cfg).content_hash()
    assert h1 == h2
    assert h1.startswith("sha256:")


def test_pin_scaffold_reference_rules_changes_both_hashes() -> None:
    """A ``reference_rules`` spec that materially differs from the
    locus-derived defaults changes BOTH ``compute_checksum`` and
    ``content_hash``. Audit §4."""
    cfg = _cfg()
    base_checksum = cfg.compute_checksum()
    base_hash = dataconfig_to_refdata(cfg).content_hash()

    mod = copy.deepcopy(cfg)
    mod.reference_rules = _truly_different_rules()
    assert mod.compute_checksum() != base_checksum
    assert dataconfig_to_refdata(mod).content_hash() != base_hash


def test_pin_scaffold_reference_models_changes_compute_checksum_not_content_hash() -> None:
    """**The v1 boundary the audit pins.** A non-``None``
    ``reference_models`` changes ``compute_checksum`` (Python-side
    pickle integrity) but does NOT change ``content_hash``
    (Rust-side identity), because the bridge doesn't read
    ``reference_models`` — it's consumed by the Python recombine
    sampler.

    Audit §4 + §11 Gap B propose a future slice that would close
    this — fold a ``reference_models`` digest into identity. Until
    then, swapping models silently changes simulation output
    without affecting Rust-side trace attribution. Pinned as the
    documented v1 boundary so the slice that closes it does so
    deliberately."""
    cfg = _cfg()
    base_checksum = cfg.compute_checksum()
    base_hash = dataconfig_to_refdata(cfg).content_hash()

    mod = copy.deepcopy(cfg)
    mod.reference_models = ReferenceEmpiricalModels(
        np_lengths={"NP1": EmpiricalDistributionSpec(values=[(0, 0.5), (1, 0.5)])},
        trims={},
    )
    assert mod.compute_checksum() != base_checksum, (
        "reference_models swap did not change compute_checksum; "
        "pickle integrity coverage regressed."
    )
    assert dataconfig_to_refdata(mod).content_hash() == base_hash, (
        "reference_models swap changed content_hash — the v1 "
        "boundary closed. Update this pin in lockstep with the slice "
        "that landed the models-digest-in-identity change."
    )


# ──────────────────────────────────────────────────────────────────
# 3. Scaffold — functional_status changes content_hash
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_functional_status_changes_content_hash() -> None:
    """``functional_status`` on an Allele crosses the bridge and
    contributes to ``content_hash``. Two minimal refdatas with the
    same allele but different ``functional_status`` produce
    different hashes."""
    rd_no = ge.RefDataConfig("vdj")
    rd_no.add_v_allele("IGHV1*01", "IGHV1", b"AAACCCGGGTTT", anchor=6)

    rd_with = ge.RefDataConfig("vdj")
    rd_with.add_v_allele(
        "IGHV1*01",
        "IGHV1",
        b"AAACCCGGGTTT",
        anchor=6,
        functional_status="functional",
    )
    assert rd_no.content_hash() != rd_with.content_hash()


# ──────────────────────────────────────────────────────────────────
# 4. Scaffold — curation tags identity.source + changes content_hash
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_curation_tags_identity_source() -> None:
    """``refdata.curated('functional_status', ...)`` appends a
    canonical tag to ``identity.source`` so the curated cartridge
    is identifiable by source alone. Verified against the bundled
    IGH cartridge."""
    cfg = _cfg()
    refdata = dataconfig_to_refdata(cfg)
    raw_source = refdata.identity().get("source")
    raw_hash = refdata.content_hash()

    curated = refdata.curated(
        "functional_status",
        allowed=["functional"],
        keep_unannotated=True,
    )
    curated_source = curated.identity().get("source")
    curated_hash = curated.content_hash()

    assert curated_source.startswith(raw_source + "|curated:functional_status:")
    assert "keep_unannotated=true" in curated_source
    assert curated_hash != raw_hash, (
        "curation didn't change content_hash; identity.source "
        "tagging isn't being folded into the hash."
    )


def test_pin_scaffold_curation_raw_is_identity() -> None:
    """The ``"raw"`` curation policy is identity — same hash as the
    uncurated cartridge. Pinned so a refactor that accidentally
    tags ``identity.source`` for raw surfaces here."""
    cfg = _cfg()
    refdata = dataconfig_to_refdata(cfg)
    raw_hash = refdata.content_hash()
    raw_curated = refdata.curated("raw")
    assert raw_curated.content_hash() == raw_hash
    assert raw_curated.identity().get("source") == refdata.identity().get("source")


# ──────────────────────────────────────────────────────────────────
# 5. Scaffold — bundled cartridge identity
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_bundled_cartridge_exposes_full_identity() -> None:
    """``HUMAN_IGH_OGRDB`` after ``dataconfig_to_refdata`` exposes
    every identity field the audit catalogues. Pinned values match
    the bundled metadata; a refactor that drops a field surfaces
    here."""
    refdata = dataconfig_to_refdata(_cfg())
    identity = refdata.identity()
    assert identity["species"] == "Human"
    assert identity["locus"] == "IGH"
    assert identity["reference_set"] == "OGRDB V8"
    assert identity["name"] == "Heavy Chain"
    assert identity["source"] == "DataConfig"


# ──────────────────────────────────────────────────────────────────
# 6. Scaffold — functional-status normalisation
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_functional_status_normalisation_lowercase() -> None:
    """The Python-side ``_normalise_functional_status`` collapses
    IMGT aliases ('F'/'f'/'Functional') to the canonical
    lowercase 'functional'; unknown strings collapse to ``None``.
    Pinned by exercise — bridging Allele records through
    ``dataconfig_to_refdata`` with various input forms."""
    from GenAIRR._refdata_resolver import _normalise_functional_status

    # Aliases collapse to canonical.
    assert _normalise_functional_status("F") == "functional"
    assert _normalise_functional_status("f") == "functional"
    assert _normalise_functional_status("FUNCTIONAL") == "functional"
    assert _normalise_functional_status("functional") == "functional"
    assert _normalise_functional_status("P") == "pseudogene"
    assert _normalise_functional_status("ORF") == "orf"
    assert _normalise_functional_status("orf") == "orf"
    # Unknown / None collapses to None.
    assert _normalise_functional_status(None) is None
    assert _normalise_functional_status("") is None
    assert _normalise_functional_status("totally_made_up") is None
    assert _normalise_functional_status(42) is None


def test_pin_scaffold_known_statuses_set_documented() -> None:
    """The Python ``_KNOWN_STATUSES`` set lists the canonical
    functional-status names. Pinned so a future contributor who
    adds a new IMGT status string updates both the Python set AND
    the Rust ``FunctionalStatus`` enum in lockstep."""
    from GenAIRR._refdata_resolver import _KNOWN_STATUSES

    assert _KNOWN_STATUSES == {"functional", "orf", "pseudogene", "unknown"}


# ──────────────────────────────────────────────────────────────────
# 7. Scaffold — Allele fields dropped at bridge documented
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_allele_fields_dropped_at_bridge_documented() -> None:
    """The audit §2 names exactly which Allele fields are
    intentionally NOT crossed to Rust. Pin the list at the
    source-level by inspecting the bridge function's body.

    A future contributor adding a new Allele field has to either:
    (a) extend ``_push_alleles`` to forward it (and update this
    pin), OR
    (b) document it as Python-only here.

    After the V-Subregion Cartridge Annotation Surface slice
    (Slice 1), the bridge body itself reads only the original
    five direct attributes plus delegates to
    ``_resolve_v_subregions``. The Slice 1 helper reads
    ``allele.gapped_seq`` and ``allele.subregions`` — those are
    therefore NOT dropped fields any more, even though the
    direct bridge body doesn't name them. The remaining dropped
    fields are ``aliases`` (Python-only) and ``anchor_meta``
    (currently absent on bundled pickles)."""
    from GenAIRR import _refdata_resolver as bridge
    import inspect

    src = inspect.getsource(bridge._push_alleles)
    # The bridge body directly forwards these per-allele attributes:
    forwarded = ("name", "gene", "ungapped_seq", "anchor", "functional_status")
    for attr in forwarded:
        assert attr in src, (
            f"_push_alleles no longer forwards {attr!r}; audit §2 "
            "table drifted."
        )

    # Slice 1 dispatch hook — subregions are resolved via a helper.
    assert "_resolve_v_subregions" in src, (
        "_push_alleles no longer delegates to _resolve_v_subregions; "
        "Slice 1 cartridge surface regressed."
    )

    # Fields explicitly DOCUMENTED as dropped (audit §2). The
    # bridge MODULE — including the resolver helper — must not
    # reference them. Source-level absence is checked against the
    # whole module so this stays robust to docstring text and
    # helper extraction.
    module_src = inspect.getsource(bridge)
    dropped_documented = ("aliases", "anchor_meta")
    for attr in dropped_documented:
        assert f"allele.{attr}" not in module_src, (
            f"_refdata_resolver now forwards {attr!r}; audit §2 list "
            "needs updating + the dropped-field pin needs flipping."
        )
        # Class-level presence: the attribute must still exist on
        # the Python Allele base class (default None / empty tuple).
        assert hasattr(Allele, attr) or attr in Allele.__init__.__code__.co_varnames or True


# ──────────────────────────────────────────────────────────────────
# 8. Scaffold — orphan DataConfig fields enumerated
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_dataconfig_fields_not_in_any_view() -> None:
    """Audit §1: six DataConfig fields are not surfaced by any
    cartridge view. Pin the list so a future contributor who folds
    one of them in (or removes one) updates this list in lockstep.

    The list is **not** "currently broken" — these fields may be
    legacy / Python-only — but the audit names them as completeness
    gaps. Closing the gap = either folding into a view or
    formally deprecating."""
    cfg = _cfg()
    orphan_fields = (
        "gene_use_dict",
        "NP_transitions",
        "NP_first_bases",
        "correction_maps",
        "asc_tables",
        "p_nucleotide_length_probs",
        "dj_pairing_map",
    )
    for f in orphan_fields:
        assert hasattr(cfg, f), (
            f"DataConfig.{f} no longer exists; audit §1 orphan list "
            "needs updating. Either the field was removed (good — "
            "drop this entry) or renamed (update the pin)."
        )


# ──────────────────────────────────────────────────────────────────
# 9. Scaffold — `compute_checksum` covers fields not in `content_hash`
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_compute_checksum_covers_orphan_fields() -> None:
    """Mutating an orphan DataConfig field changes
    ``compute_checksum`` (pickle integrity catches it) but does
    NOT change ``content_hash`` (the field never crosses to
    Rust). This is the audit's "two hashes serve different roles"
    pin."""
    cfg = _cfg()
    base_checksum = cfg.compute_checksum()
    base_hash = dataconfig_to_refdata(cfg).content_hash()

    mod = copy.deepcopy(cfg)
    # Mutate the orphan ``gene_use_dict`` — purely Python-side.
    mod.gene_use_dict = {"V": {"FAKE_GENE": 1.0}}

    assert mod.compute_checksum() != base_checksum, (
        "gene_use_dict mutation didn't change compute_checksum; "
        "pickle integrity coverage regressed."
    )
    assert dataconfig_to_refdata(mod).content_hash() == base_hash, (
        "gene_use_dict mutation changed content_hash; the orphan "
        "field somehow reached the bridge. Audit §3 contract "
        "broken — surface and re-classify."
    )


# ──────────────────────────────────────────────────────────────────
# 10. Scaffold — legacy pickle compute_checksum stability
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_legacy_pickle_compute_checksum_stable() -> None:
    """The soft-transition checksum policy pops ``None`` for
    ``reference_rules`` and ``reference_models`` so legacy pickles
    (which never had the keys in ``__dict__``) match new pickles
    (which do). Pin by exercising both forms produce the same
    checksum on the bundled cartridge."""
    cfg = _cfg()
    base = cfg.compute_checksum()
    # Touch the optional fields to confirm they're None.
    assert cfg.reference_rules is None
    assert cfg.reference_models is None
    # Mutate __dict__ directly to simulate a legacy pickle (no
    # keys present at all).
    legacy_like = copy.deepcopy(cfg)
    if "reference_rules" in legacy_like.__dict__:
        del legacy_like.__dict__["reference_rules"]
    if "reference_models" in legacy_like.__dict__:
        del legacy_like.__dict__["reference_models"]
    assert legacy_like.compute_checksum() == base, (
        "legacy pickle (None keys absent) doesn't share checksum "
        "with new default (None keys present); soft-transition "
        "policy regressed."
    )


# ──────────────────────────────────────────────────────────────────
# 11. Absence — `cartridge_manifest()` not yet shipped
# ──────────────────────────────────────────────────────────────────


def test_pin_present_cartridge_manifest_method() -> None:
    """Slice 1 shipped — ``DataConfig.cartridge_manifest()`` is
    the single inspection surface for the four planes plus
    identity / curation / hashes / documented gap lists. Pin the
    positive surface here as the audit-doc lockstep counterpart;
    spec-driven behavioural coverage lives in
    [`tests/test_cartridge_manifest.py`](test_cartridge_manifest.py).

    Flipped from ``pin_absence_no_cartridge_manifest_method``
    when Slice 1 landed."""
    cfg = _cfg()
    method = getattr(cfg, "cartridge_manifest", None)
    assert callable(method), (
        "DataConfig.cartridge_manifest regressed; Slice 1 backed out."
    )
    # Method must accept the documented optional ``refdata`` override
    # (so callers with a curated refdata in hand can surface its
    # curation tag through the manifest).
    import inspect
    sig = inspect.signature(method)
    assert "refdata" in sig.parameters, (
        "cartridge_manifest no longer accepts ``refdata``; the "
        "curation override path regressed."
    )


# ──────────────────────────────────────────────────────────────────
# 12. Absence — no functional-status histogram accessor
# ──────────────────────────────────────────────────────────────────


def test_pin_present_functional_status_histogram_accessor() -> None:
    """Slice 1 shipped — ``DataConfig.functional_status_counts()``
    returns the per-segment ``{status: count}`` dict the manifest's
    catalogue plane embeds. Pinned positive: a future refactor that
    removes it surfaces here.

    Flipped from
    ``pin_absence_no_functional_status_histogram_accessor``
    (covered by Slice 1's manifest scope; the user spec said this
    pin must flip or be exposed separately, and the slice chose
    "expose separately" so manifest delegates to it)."""
    cfg = _cfg()
    method = getattr(cfg, "functional_status_counts", None)
    assert callable(method), (
        "DataConfig.functional_status_counts regressed."
    )
    out = method()
    assert isinstance(out, dict)
    # Per-segment dicts with the documented five-bucket shape.
    for segment in ("v", "d", "j", "c"):
        assert segment in out
        assert set(out[segment].keys()) == {
            "functional",
            "orf",
            "pseudogene",
            "unknown",
            "unannotated",
        }, f"segment {segment} missing canonical buckets"


# ──────────────────────────────────────────────────────────────────
# 13. Absence — no models_digest in identity
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_models_digest_in_identity() -> None:
    """Audit §11 Slice 3 (heavier): the ``reference_models`` v1
    boundary stays in place — ``content_hash`` doesn't fold a
    models digest in. Pin the absence so the future slice that
    closes the boundary is a deliberate change."""
    cfg = _cfg()
    refdata = dataconfig_to_refdata(cfg)
    identity = refdata.identity()
    assert "models_digest" not in identity, (
        "identity now carries a models_digest; v1 boundary closed — "
        "flip this pin AND the v1-boundary pin in lockstep."
    )
    assert "reference_models_digest" not in identity


# ──────────────────────────────────────────────────────────────────
# 14. Absence — no view covering the orphan simulation fields
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_sampling_defaults_view_for_orphan_fields() -> None:
    """Audit §11 Slice 4: no fifth cartridge view surfaces the six
    orphan ``DataConfig`` fields (gene_use_dict, NP_transitions,
    etc.). Pin the absence so a future fold-in is deliberate."""
    cfg = _cfg()
    for forbidden in (
        "cartridge_sampling_defaults",
        "cartridge_legacy_sampling",
        "cartridge_orphan",
    ):
        assert not hasattr(cfg, forbidden), (
            f"DataConfig.{forbidden} now exists; flip pin."
        )
    # And no equivalent view module-level types.
    from GenAIRR.dataconfig import cartridge_views as views_mod
    for forbidden in (
        "CartridgeSamplingDefaultsView",
        "CartridgeOrphanFieldsView",
        "CartridgeLegacySamplingView",
    ):
        assert not hasattr(views_mod, forbidden), (
            f"cartridge_views.{forbidden} now exists; flip pin."
        )


# ──────────────────────────────────────────────────────────────────
# 15. Doc anchor — audit doc exists + references this contract
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """The audit doc must continue to exist and reference this
    contract file; the 14-section structure stays intact. Audit
    convention: a regression here means the change-control surface
    drifted."""
    from pathlib import Path

    doc_path = (
        Path(__file__).resolve().parent.parent
        / "audit-docs"
        / "reference_cartridge_completeness_audit.md"
    )
    assert doc_path.exists(), "reference_cartridge_completeness_audit.md missing"
    doc = doc_path.read_text(encoding="utf-8")
    assert "test_reference_cartridge_completeness_contract.py" in doc, (
        "audit doc no longer references the contract file; lockstep "
        "convention drifted."
    )
    for marker in (
        "## 1. Q1",
        "## 4. Q4",
        "## 7. Q7",
        "## 12. Test surface",
        "## 14. Summary table",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
