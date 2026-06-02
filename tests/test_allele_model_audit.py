"""Allele model audit — pin current behaviour before any biology
or DSL work runs on top of the cartridge.

Companion to ``docs/allele_model_audit.md``. These tests
characterize the legacy Allele layer at the moment the cartridge
architecture (identity / rules / catalogue / curation / empirical
models) was completed. They are deliberately read-only: no
refactor, no behaviour change, just guarantees that the next
refactor must preserve (or knowingly break).

Key claims pinned here:

1. ``dataconfig_to_refdata`` transfers exactly four fields per
   allele into Rust: ``name``, ``gene``, ``ungapped_seq``, ``anchor``.
2. ``functional_status``, ``anchor_meta``, ``aliases``, ``locus``,
   ``species``, ``source``, and ``gapped_seq`` do not cross.
3. ``reference_models`` (and the legacy nested-dict extractor)
   drive recombine trim defaults; the Allele-resident
   ``get_trimmed`` / ``_get_trim_length`` methods are orphaned.
4. ``anchor_override`` sets ``anchor`` but leaves
   ``anchor_meta is None`` — a known provenance gap.
5. ``gene`` / ``family`` parsing by string split has predictable
   edge cases (single ``*``, no ``-``, multi-hyphen names).
"""
from __future__ import annotations

import inspect
import pickle
from pathlib import Path
from typing import Any, Dict, List

import pytest

import GenAIRR as ga
from GenAIRR._dataconfig_extract import extract_recombine_defaults
from GenAIRR._refdata_resolver import _resolve_config_name, dataconfig_to_refdata
from GenAIRR.alleles.allele import Allele, CAllele, DAllele, JAllele, VAllele
from GenAIRR.dataconfig.data_config import DataConfig
from GenAIRR.reference_models import (
    EmpiricalDistributionSpec,
    ReferenceEmpiricalModels,
)


# ──────────────────────────────────────────────────────────────────
# 1. dataconfig_to_refdata: which fields cross into Rust?
# ──────────────────────────────────────────────────────────────────


def _first_allele(cfg: DataConfig, segment: str):
    """Pull one allele from the bundled IGH catalogue, indexed by
    ``cfg.X_alleles[gene][0]``. Used to inspect what the Python side
    holds before the bridge."""
    pool = getattr(cfg, f"{segment}_alleles")
    assert pool, f"bundled IGH cfg has empty {segment}_alleles"
    first_gene = next(iter(pool.keys()))
    return pool[first_gene][0]


def test_bridge_transfers_expected_fields_per_allele() -> None:
    """The Rust ``Allele`` carries the cross-bridge field set pinned
    by ``docs/allele_model_audit.md``. After the V-Subregion
    Cartridge Annotation Surface slice (Slice 1) the bridge reads
    five direct fields — ``name`` + ``gene`` + ``ungapped_seq`` +
    ``anchor`` + ``functional_status`` — and derives a sixth
    payload, ``subregions``, via the resolver helper
    ``_resolve_v_subregions`` (which itself reads
    ``allele.gapped_seq`` / ``allele.subregions``). Other
    Python-only fields stay dropped; promoting one of them is a
    documented event that must update both the audit doc and
    this test."""
    from GenAIRR import _refdata_resolver as resolver

    push_src = inspect.getsource(resolver._push_alleles)
    # Positive — directly forwarded fields.
    for needed in (
        "allele.name",
        "allele.gene",
        "allele.ungapped_seq",
        "allele.anchor",
        "allele.functional_status",
    ):
        assert needed in push_src, (
            f"_push_alleles must read {needed} when lowering into Rust; "
            f"if it stops reading this, the bridge contract changed"
        )
    # Slice 1 — subregion payload is resolved through a dedicated
    # helper rather than read directly in the bridge body.
    assert "_resolve_v_subregions" in push_src, (
        "_push_alleles no longer delegates to _resolve_v_subregions; "
        "Slice 1 cartridge surface regressed."
    )
    resolver_src = inspect.getsource(resolver._resolve_v_subregions)
    for needed in ("allele.gapped_seq", "allele.subregions", "compute_v_region_boundaries"):
        assert (
            needed in resolver_src
            or needed.split(".")[-1] in resolver_src
        ), (
            f"_resolve_v_subregions must touch {needed} to deliver the "
            "Slice 1 derivation; cartridge surface regressed."
        )
    # Negative: dropped fields must NOT be referenced anywhere in
    # the bridge body — otherwise the audit's claim about which
    # fields cross is stale and the bridge has silently expanded.
    for dropped in (
        "allele.anchor_meta",
        "allele.aliases",
        "allele.locus",
        "allele.species",
        "allele.source",
    ):
        assert dropped not in push_src, (
            f"_push_alleles reads {dropped} but the audit claims that field "
            f"is dropped; either the bridge changed or the audit is stale"
        )


def test_rust_refdata_exposes_only_documented_fields() -> None:
    """The Rust-side ``Allele`` PyO3 class exposes the documented
    cartridge fields. As of the functional-status slice this is
    name / gene / seq / segment / anchor / functional_status. If a
    future slice adds another field (anchor_meta, aliases, etc.),
    this test surfaces it so the audit doc can be updated in
    lockstep."""
    rd = ga.Experiment.on("human_igh").refdata
    py_allele = rd.v_allele(0)
    surface = {a for a in dir(py_allele) if not a.startswith("_")}
    # Currently exposed (pinned).
    for expected in (
        "name",
        "gene",
        "seq",
        "anchor",
        "segment",
        "functional_status",
        "subregions",  # Slice 1: V-region substructure annotations
    ):
        assert expected in surface, f"Rust Allele lost expected attr {expected!r}"
    # Confirmed-dropped — if any of these appear, update the audit.
    for dropped in ("anchor_meta", "aliases", "locus", "source"):
        assert dropped not in surface, (
            f"Rust Allele unexpectedly exposes {dropped!r}; update "
            f"docs/allele_model_audit.md if this is intentional"
        )


def test_python_to_rust_round_trip_preserves_only_four_fields() -> None:
    """End-to-end: load human_igh DataConfig, bridge to RefDataConfig,
    read back the first V allele. The values are byte-equal to the
    Python source for the four expected fields."""
    cfg = _resolve_config_name("human_igh")
    py_v = _first_allele(cfg, "v")
    rd = dataconfig_to_refdata(cfg)
    rust_v = rd.v_allele(0)
    assert rust_v.name == py_v.name
    assert rust_v.gene == py_v.gene
    assert rust_v.seq() == py_v.ungapped_seq.encode("ascii")
    assert rust_v.anchor == py_v.anchor


# ──────────────────────────────────────────────────────────────────
# 2. Allele "extra" fields — what's still on the Python side?
# ──────────────────────────────────────────────────────────────────


def test_allele_extra_fields_have_class_level_defaults() -> None:
    """``Allele`` defines class-level defaults so legacy pickles
    built before T2-8 (when these fields were added) continue to
    surface them as ``None``/``()``. Pin those defaults."""
    assert Allele.anchor_meta is None
    assert Allele.functional_status is None
    assert Allele.locus is None
    assert Allele.aliases == ()
    assert Allele.species is None
    assert Allele.source is None


def test_bundled_alleles_do_not_currently_populate_extras() -> None:
    """The bundled builders populate ``name`` / ``gene`` / ``family``
    / ``gapped_seq`` / ``ungapped_seq`` / ``length`` / ``anchor`` /
    ``anchor_meta`` (the last only when ``_find_anchor`` ran), and
    leave the other extras at their class defaults.

    Pseudogenes for which the C resolver rejects the anchor produce
    ``anchor=None`` with ``anchor_meta`` carrying the rejection
    reason — useful diagnostic, but not transferred to Rust.
    """
    cfg = _resolve_config_name("human_igh")
    v0 = _first_allele(cfg, "v")
    j0 = _first_allele(cfg, "j")
    # Core fields the engine consumes:
    for a in (v0, j0):
        assert isinstance(a.name, str) and a.name
        assert isinstance(a.gene, str) and a.gene
        assert isinstance(a.ungapped_seq, str) and a.ungapped_seq
    # T2-8 anchor_meta IS populated when the C resolver ran.
    # (May still be None on legacy pickles built before T2-8 — the
    # check is "if present, it's an AnchorResult"; we just confirm
    # the attribute is accessible without error.)
    _ = v0.anchor_meta
    _ = j0.anchor_meta
    # Extras the cartridge could carry but currently doesn't.
    for a in (v0, j0):
        assert a.functional_status is None
        assert a.aliases == ()
        assert a.locus is None
        # ``species`` and ``source`` may be class-default None.
        assert a.species in (None, "Human")
        assert a.source is None or isinstance(a.source, str)


# ──────────────────────────────────────────────────────────────────
# 3. Trim methods on Allele are orphaned (not called in production)
# ──────────────────────────────────────────────────────────────────


def test_get_trimmed_has_no_production_caller_outside_allele_module() -> None:
    """The trim methods (intentionally not naming them literally in
    this docstring, so the audit test doesn't grep itself) are
    defined on every Allele subclass but never invoked from outside
    ``alleles/allele.py`` — the engine's trim path goes through the
    typed/legacy extractor in ``_dataconfig_extract``.

    The grep below is the audit's actual contract: if any other
    source file starts calling the trim helpers, this test surfaces
    it so the audit doc can be re-evaluated before that call site
    lands.
    """
    import subprocess
    repo = subprocess.run(
        ["git", "rev-parse", "--show-toplevel"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    hits = subprocess.run(
        [
            "grep", "-rn", "-E",
            # The two trim entry points the audit pins. The pattern
            # itself does not appear in this test file body — only
            # in the assembled regex below — so we never grep
            # ourselves.
            r"\.get_trimmed\(|\._get_trim_length\(",
            f"{repo}/src/GenAIRR",
            f"{repo}/tests",
        ],
        capture_output=True, text=True,
    ).stdout.splitlines()
    # Strip self-references inside the allele module (legitimate
    # internal dispatch) and this audit file (it documents the
    # method names in code, not as caller).
    external = [
        h for h in hits
        if "src/GenAIRR/alleles/allele.py" not in h
        and __file__ not in h
    ]
    assert external == [], (
        f"unexpected production caller(s) of trim methods: {external!r}; "
        f"update docs/allele_model_audit.md and decide whether trim should "
        f"stay on Allele or fully move into the cartridge."
    )


def test_recombine_trim_defaults_come_from_extractor_not_alleles() -> None:
    """The extractor's output drives recombine defaults — not the
    Allele methods. Setting ``reference_models.trims["V_3"]`` shapes
    the recombine step; the underlying allele's ``_get_trim_length``
    is never called for it."""
    cfg = _resolve_config_name("human_igh").copy()
    cfg.reference_models = ReferenceEmpiricalModels(
        trims={"V_3": EmpiricalDistributionSpec([(2, 1.0), (5, 1.0)])},
    )
    defaults = extract_recombine_defaults(cfg)
    pairs = defaults["trim_v_3"]
    # The extractor lowered the explicit model directly — no
    # Allele-side trim sampling happened.
    assert pairs is not None
    values = {v for v, _ in pairs}
    assert {2, 5}.issubset(values)


# ──────────────────────────────────────────────────────────────────
# 4. anchor_override sets anchor but leaves anchor_meta None
# ──────────────────────────────────────────────────────────────────


class _SilentVAllele(VAllele):
    """V allele that quietly does nothing in ``_find_anchor``.
    Lets the parsing tests exercise ``__init__`` without depending
    on the C resolver (``_native._anchor`` may not be built in
    every dev environment).
    """

    def _find_anchor(self):
        # Quietly leave self.anchor as the value __init__ set
        # (None by default; the override path bypasses this method
        # entirely).
        return


class _ExplosiveVAllele(VAllele):
    """V allele whose ``_find_anchor`` raises — used to verify
    that ``anchor_override`` bypasses the resolver entirely.
    """

    def _find_anchor(self):  # pragma: no cover — should not run
        raise AssertionError("anchor_override should suppress _find_anchor")


def test_anchor_override_sets_anchor_and_leaves_meta_none() -> None:
    """Known provenance gap: ``anchor_override`` is opaque. The
    audit doc lists this as a remediation candidate (record a
    synthetic ``AnchorResult`` so ``anchor_meta`` is never
    silently None on the override path). This test pins the
    current behaviour so any future change is intentional."""
    a = _ExplosiveVAllele(
        name="MyV1*01", gapped_sequence="A" * 300, length=300,
        anchor_override=108,
    )
    assert a.anchor == 108
    assert a.anchor_meta is None, (
        "Audit pins: anchor_override leaves anchor_meta=None. "
        "If you fix this by stamping a synthetic AnchorResult, "
        "update the audit doc + this test."
    )


def test_anchor_override_bypasses_find_anchor() -> None:
    """Override should *replace* resolution. _ExplosiveVAllele's
    `_find_anchor` raises — if `anchor_override` were to bypass
    correctly, no exception."""
    _ExplosiveVAllele(
        name="MyV1*01", gapped_sequence="A" * 300, length=300,
        anchor_override=42,
    )


# ──────────────────────────────────────────────────────────────────
# 5. gene / family parsing by string split — edge cases
# ──────────────────────────────────────────────────────────────────


def _make_v(name: str) -> VAllele:
    """Construct a VAllele just to exercise the name-parsing
    code inside ``__init__``. The gapped sequence is dummy and the
    silent stub skips the C resolver."""
    return _SilentVAllele(name=name, gapped_sequence="A" * 30, length=30)


@pytest.mark.parametrize(
    "name,expected_family,expected_gene",
    [
        # Canonical IMGT names (the common case).
        ("IGHV1-2*01",      "IGHV1",   "IGHV1-2"),
        ("IGHJ4*02",        "IGHJ4",   "IGHJ4"),
        ("IGKV3-20*01",     "IGKV3",   "IGKV3-20"),
        ("TRBV20-1*01",     "TRBV20",  "TRBV20-1"),
        # OGRDB-style F-numbered families.
        ("IGHVF1-G3*01",    "IGHVF1",  "IGHVF1-G3"),
        # Multi-hyphen gene names — split('-')[0] stops at the first
        # hyphen, so family loses the back half (documented).
        ("IGHV1-2-3*01",    "IGHV1",   "IGHV1-2-3"),
        # No hyphen in the name — family falls back to gene root.
        ("IGHV*01",         "IGHV",    "IGHV"),
        # No allele suffix at all (rare).
        ("IGHV1-2",         "IGHV1",   "IGHV1-2"),
    ],
)
def test_gene_family_split_pins_documented_cases(
    name: str, expected_family: str, expected_gene: str,
) -> None:
    """The split logic in ``Allele.__init__``:
        family = name.split('-')[0] if '-' in name else name.split('*')[0]
        gene   = name.split('*')[0]
    Pin the cases the audit doc lists so any future refactor (e.g.
    locus-aware parsing) is a deliberate decision, not silent drift.
    """
    a = _make_v(name)
    assert a.family == expected_family, f"name={name!r}: family"
    assert a.gene == expected_gene, f"name={name!r}: gene"


# ──────────────────────────────────────────────────────────────────
# 6. The audit documents itself
# ──────────────────────────────────────────────────────────────────


def test_audit_doc_exists_and_lists_every_pinned_topic() -> None:
    """The audit answers seven questions. If the doc loses any
    section header, the test surfaces it so the doc and the
    behaviour-pinning tests don't drift apart."""

    doc = Path(__file__).resolve().parent.parent / "audit-docs" / "allele_model_audit.md"
    assert doc.is_file(), f"docs/allele_model_audit.md is missing at {doc}"
    text = doc.read_text(encoding="utf-8")
    for header in (
        "Production paths that still use legacy methods",
        "Fields that cross into Rust",
        "Fields that are dropped",
        "Candidates for the cartridge catalogue",
        "Is gene / family parsing reliable?",
        "Is `anchor_override` provenance acceptable?",
        "Should trimming still belong on allele objects?",
    ):
        assert header in text, (
            f"docs/allele_model_audit.md is missing section {header!r}; "
            f"audit doc and test corpus must stay in sync"
        )
