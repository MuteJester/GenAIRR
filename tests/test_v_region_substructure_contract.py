"""Contract pins for the V-Region Substructure / CDR-FR Targeting
audit.

Companion to
[`docs/v_region_substructure_audit.md`](../docs/v_region_substructure_audit.md).
The audit's pre-implementation pins have been partially
discharged by Slice 1 (V-Subregion Cartridge Annotation
Surface). The remaining absence pins gate Slice 2
(``v_subregion_rates``) and Slice 3 (CDR/FR mutation counters +
``*SubregionMismatch`` validator kinds).

Split:

- ``pin_scaffold_*`` tests freeze the bundled metadata that
  the cartridge slices build on: ``gapped_seq`` coverage on all
  V alleles, IMGT-gapped boundary constants, the derivation
  helper producing clean intervals on every bundled allele, the
  Python ``Allele.anchor_meta`` field's class-level presence,
  the bridge's existing forwarding contract, the audit doc
  anchor.
- ``pin_present_*`` tests pin Slice 1's surface so a regression
  surfaces here: the Rust ``Allele.subregions`` field, the Python
  ``Allele.subregions`` attribute, the ``VSubregionLabel`` /
  ``VSubregion`` Rust types, the manifest's
  ``v_subregion_support`` block, the ``_refdata_resolver``
  consumer of ``imgt_regions``.
- ``pin_absence_*`` tests still freeze the gaps Slice 2+ closes:
  no ``v_subregion_rates`` kwarg on ``Experiment.mutate``, no
  CDR/FR mutation counters on AIRR records, no
  ``*SubregionMismatch`` validator issue kinds, and (unrelated)
  the ``anchor_meta`` absence on bundled pickles.
"""
from __future__ import annotations

import inspect
from pathlib import Path

import pytest

import GenAIRR as ga


_REPO_ROOT = Path(__file__).resolve().parent.parent


# ──────────────────────────────────────────────────────────────────
# 1. Scaffold — bundled V alleles carry gapped_seq
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_bundled_v_alleles_carry_gapped_seq() -> None:
    """Audit §1 pre-flight: every V allele in
    ``HUMAN_IGH_OGRDB`` has a populated ``gapped_seq``. The
    derivation slice depends on this; pin it as a hard
    requirement on bundled cartridges."""
    cfg = ga.HUMAN_IGH_OGRDB
    total = 0
    with_gapped = 0
    for _gene, alleles in cfg.v_alleles.items():
        for a in alleles:
            total += 1
            if getattr(a, "gapped_seq", None):
                with_gapped += 1
    assert total > 0, "bundled cartridge has no V alleles"
    assert with_gapped == total, (
        f"only {with_gapped}/{total} V alleles carry gapped_seq; "
        "subregion derivation cannot run on the missing ones."
    )


def test_pin_scaffold_bundled_v_gapped_seq_uses_imgt_dot_convention() -> None:
    """Audit §1: bundled V alleles' ``gapped_seq`` uses the IMGT
    dot convention (``.`` characters for gap positions). Pin
    presence of at least one gap on a sample allele so a refactor
    that switches conventions surfaces here."""
    cfg = ga.HUMAN_IGH_OGRDB
    sample = list(cfg.v_alleles.values())[0][0]
    gapped = sample.gapped_seq
    assert "." in gapped, (
        f"V allele {sample.name} gapped_seq has no '.' gaps; the "
        "IMGT dot convention may have been replaced."
    )
    # gapped_seq must be at least as long as ungapped_seq.
    assert len(gapped) >= len(sample.ungapped_seq)
    # Number of gaps = len_gapped - len_ungapped.
    assert gapped.count(".") == len(gapped) - len(sample.ungapped_seq)


# ──────────────────────────────────────────────────────────────────
# 2. Scaffold — Python Allele.anchor_meta field exists
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_anchor_meta_field_exists_on_python_allele() -> None:
    """Audit §1: ``Allele.anchor_meta`` is declared on the Python
    class (T2-8 anchor provenance). Pin the class-level presence
    so a refactor that removes the field surfaces here. Whether
    bundled pickles populate it is a separate pin (§ absence)."""
    from GenAIRR.alleles.allele import Allele

    assert hasattr(Allele, "anchor_meta"), (
        "Allele.anchor_meta class attribute removed; T2-8 anchor "
        "provenance surface lost."
    )


# ──────────────────────────────────────────────────────────────────
# 3. Scaffold — IMGT gapped boundaries constants
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_imgt_gapped_boundaries_constants_pinned() -> None:
    """Audit §3: the IMGT-canonical gapped nucleotide boundaries
    are pinned in ``utilities/imgt_regions.IMGT_GAPPED_BOUNDARIES``.
    Five regions, half-open intervals over gapped positions.
    A refactor that changes the constants forces the slice
    author to update both sides in lockstep."""
    from GenAIRR.utilities.imgt_regions import IMGT_GAPPED_BOUNDARIES

    assert IMGT_GAPPED_BOUNDARIES == {
        "FWR1": (0, 78),     # IMGT aa 1-26
        "CDR1": (78, 114),   # IMGT aa 27-38
        "FWR2": (114, 165),  # IMGT aa 39-55
        "CDR2": (165, 195),  # IMGT aa 56-65
        "FWR3": (195, 312),  # IMGT aa 66-104
    }


def test_pin_scaffold_imgt_gapped_boundaries_match_imgt_aa_positions() -> None:
    """The gapped nucleotide ranges align with the documented IMGT
    AA codon positions: FWR1 = codons 1-26 (78 nt), CDR1 = codons
    27-38 (36 nt), FWR2 = codons 39-55 (51 nt), CDR2 = codons
    56-65 (30 nt), FWR3 = codons 66-104 (117 nt). Each region
    must be a multiple of 3 (codon-aligned)."""
    from GenAIRR.utilities.imgt_regions import IMGT_GAPPED_BOUNDARIES

    for region, (start, end) in IMGT_GAPPED_BOUNDARIES.items():
        size = end - start
        assert size % 3 == 0, (
            f"{region}: gapped boundary size {size} is not codon-"
            f"aligned (multiple of 3)."
        )


# ──────────────────────────────────────────────────────────────────
# 4. Scaffold — derivation helper works on bundled cartridge
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_compute_v_region_boundaries_works_on_all_bundled_v_alleles() -> None:
    """Audit §3: ``compute_v_region_boundaries`` produces a clean
    five-region dict for every bundled V allele in
    ``HUMAN_IGH_OGRDB``. Coverage = 100% (audit-probed). Pin both
    the helper's interface AND the coverage claim so a future
    bundled-data refresh that breaks derivation surfaces here."""
    from GenAIRR.utilities.imgt_regions import compute_v_region_boundaries

    cfg = ga.HUMAN_IGH_OGRDB
    total = 0
    fully_derived = 0
    for _gene, alleles in cfg.v_alleles.items():
        for a in alleles:
            total += 1
            try:
                boundaries = compute_v_region_boundaries(a)
            except Exception as exc:  # pragma: no cover — should not fire
                pytest.fail(
                    f"compute_v_region_boundaries failed on "
                    f"{a.name}: {exc}"
                )
            # Five expected regions.
            assert set(boundaries.keys()) == {
                "FWR1", "CDR1", "FWR2", "CDR2", "FWR3"
            }
            # Each region's interval must be half-open with
            # start <= end.
            for region, (s, e) in boundaries.items():
                assert s <= e, (
                    f"{a.name} {region}: start={s} > end={e}"
                )
            # FWR3 must have non-zero length (the conserved Cys
            # codon is in FWR3 for every functional V allele).
            if boundaries["FWR3"][1] > boundaries["FWR3"][0]:
                fully_derived += 1
    assert total > 0
    assert fully_derived == total, (
        f"only {fully_derived}/{total} V alleles produced a non-"
        "empty FWR3; bundled cartridge gapped_seq integrity "
        "regressed."
    )


def test_pin_scaffold_v_region_boundaries_intervals_are_monotonic() -> None:
    """For a canonical bundled allele, the five regions appear in
    biological order (FWR1 < CDR1 < FWR2 < CDR2 < FWR3) and form
    a contiguous partition of the V allele up to FWR3 end."""
    from GenAIRR.utilities.imgt_regions import compute_v_region_boundaries

    sample = list(ga.HUMAN_IGH_OGRDB.v_alleles.values())[0][0]
    b = compute_v_region_boundaries(sample)
    # Biological order.
    assert b["FWR1"][1] == b["CDR1"][0]
    assert b["CDR1"][1] == b["FWR2"][0]
    assert b["FWR2"][1] == b["CDR2"][0]
    assert b["CDR2"][1] == b["FWR3"][0]
    # FWR1 starts at the V allele's 5' end.
    assert b["FWR1"][0] == 0


# ──────────────────────────────────────────────────────────────────
# 5. Scaffold — bridge field-forwarding contract
# ──────────────────────────────────────────────────────────────────


def test_pin_present_bridge_forwards_documented_fields_plus_subregions() -> None:
    """Audit §2 / Slice 1: ``_push_alleles`` forwards the five
    original per-allele attributes (`name`, `gene`,
    `ungapped_seq`, `anchor`, `functional_status`) and dispatches
    V-subregion derivation through ``_resolve_v_subregions`` for
    V alleles only (D / J skip subregion plumbing). Pin via
    source-level inspection: a future contributor who removes
    either the original forwarding contract OR the subregion
    derivation hook surfaces here.

    ``allele.gapped_seq`` and ``allele.anchor_meta`` are still
    NOT forwarded directly from the bridge body — they're either
    read at derivation time inside the resolver helper, or
    deliberately stripped (anchor_meta is held back on the
    Python side; see audit §1)."""
    from GenAIRR import _refdata_resolver as bridge

    src = inspect.getsource(bridge._push_alleles)
    forwarded = ("name", "gene", "ungapped_seq", "anchor", "functional_status")
    for attr in forwarded:
        assert attr in src, (
            f"_push_alleles no longer forwards {attr!r}; audit §2 "
            "contract drifted."
        )
    # Slice 1 dispatch hook — the helper does the derivation +
    # validation; the bridge body just delegates.
    assert "_resolve_v_subregions" in src, (
        "_push_alleles no longer delegates to _resolve_v_subregions; "
        "Slice 1 surface regressed."
    )
    # ``anchor_meta`` must still NOT cross the bridge — it stays
    # on the Python side per audit §1 (the pickles never populated
    # it, and the bridge never queries it).
    bridge_module_src = inspect.getsource(bridge)
    assert "anchor_meta" not in bridge_module_src, (
        "bridge module now references anchor_meta; cartridge "
        "refresh surface drifted."
    )


def test_pin_present_v_subregion_field_on_rust_allele_struct() -> None:
    """Audit §2 / Slice 1: the Rust ``RefDataConfig::Allele``
    struct now carries a ``subregions: Vec<VSubregion>`` field
    alongside the original five. Pinned at source so a later
    refactor that drops or renames the field surfaces here."""
    refdata_src = (
        _REPO_ROOT / "engine_rs" / "src" / "refdata.rs"
    ).read_text(encoding="utf-8")
    # Find the Allele struct body.
    import re

    match = re.search(
        r"pub struct Allele\s*\{(.*?)\n\}",
        refdata_src,
        re.DOTALL,
    )
    assert match is not None, (
        "RefDataConfig::Allele struct not found; refdata.rs shape "
        "drifted."
    )
    body = match.group(1)
    # The five original fields must still be there...
    for field in (
        "pub name:",
        "pub gene:",
        "pub seq:",
        "pub segment:",
        "pub anchor:",
        "pub functional_status:",
    ):
        assert field in body, (
            f"Allele struct missing {field!r}; bridge contract "
            "drifted."
        )
    # ...plus the Slice 1 subregions field.
    assert "pub subregions:" in body, (
        "Allele struct missing pub subregions field; Slice 1 surface "
        "regressed."
    )


# ──────────────────────────────────────────────────────────────────
# 6. Absence — anchor_meta is None on bundled cartridges
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_anchor_meta_unpopulated_on_bundled_v_alleles() -> None:
    """Audit §1: the bundled ``HUMAN_IGH_OGRDB`` pickle was built
    before the ``anchor_meta`` field was added, so the field is
    always ``None`` on bundled V alleles even though the class
    declares it. Pin the absence — a future cartridge-refresh
    task that re-pickles with the resolver-populated metadata
    would flip this."""
    cfg = ga.HUMAN_IGH_OGRDB
    populated = 0
    total = 0
    for _gene, alleles in cfg.v_alleles.items():
        for a in alleles:
            total += 1
            if a.anchor_meta is not None:
                populated += 1
    assert total > 0
    assert populated == 0, (
        f"{populated}/{total} bundled V alleles now carry "
        "anchor_meta; cartridge refresh has landed — flip pin."
    )


# ──────────────────────────────────────────────────────────────────
# 7. Absence — no subregions field on Python VAllele
# ──────────────────────────────────────────────────────────────────


def test_pin_present_subregions_field_on_python_v_allele() -> None:
    """Audit §4 / Slice 1: the ``subregions`` field is now a
    first-class attribute on the Python ``Allele`` base class
    (and therefore on bundled ``VAllele`` instances). The other
    names listed here remain forbidden — only ``subregions`` is
    the canonical surface."""
    sample = list(ga.HUMAN_IGH_OGRDB.v_alleles.values())[0][0]
    assert hasattr(sample, "subregions"), (
        "Allele.subregions attribute missing; Slice 1 surface "
        "regressed."
    )
    for forbidden in (
        "imgt_subregions",
        "v_subregions",
        "fwr_cdr_boundaries",
    ):
        assert not hasattr(sample, forbidden), (
            f"VAllele.{forbidden} now exists; only ``subregions`` is "
            "the canonical Slice 1 surface — additional aliases would "
            "require a manifest update."
        )


# ──────────────────────────────────────────────────────────────────
# 8. Absence — no subregion label enum / fields in Rust
# ──────────────────────────────────────────────────────────────────


def test_pin_present_subregion_label_enum_in_rust() -> None:
    """Audit §4 / Slice 1: the Rust crate now declares a
    ``VSubregionLabel`` enum (FWR1/CDR1/FWR2/CDR2/FWR3) plus a
    ``VSubregion`` carrier struct. Source-level pin so a future
    rename has to update the contract."""
    import subprocess

    result = subprocess.run(
        ["grep", "-rn", "enum VSubregionLabel", "engine_rs/src/"],
        capture_output=True,
        text=True,
        cwd=str(_REPO_ROOT),
    )
    assert result.stdout.strip(), (
        "Rust no longer declares VSubregionLabel; Slice 1 surface "
        "regressed."
    )
    # And the companion carrier struct.
    result_struct = subprocess.run(
        ["grep", "-rn", "struct VSubregion", "engine_rs/src/"],
        capture_output=True,
        text=True,
        cwd=str(_REPO_ROOT),
    )
    assert result_struct.stdout.strip(), (
        "Rust no longer declares VSubregion struct; Slice 1 surface "
        "regressed."
    )


# ──────────────────────────────────────────────────────────────────
# 9. Absence — no v_subregion_support in manifest
# ──────────────────────────────────────────────────────────────────


def test_pin_present_v_subregion_support_in_manifest() -> None:
    """Audit §6 / Slice 1: ``cartridge_manifest()["models"]["shm"]``
    now exposes a ``v_subregion_support`` block reporting the
    canonical labels, coverage statistics, derivation source, and
    hash-participation flag. The block must carry the six
    documented keys, list the five canonical labels in order, and
    advertise that subregions enter ``refdata_content_hash``."""
    sup = (
        ga.HUMAN_IGH_OGRDB.cartridge_manifest()["models"]["shm"][
            "v_subregion_support"
        ]
    )
    assert set(sup.keys()) == {
        "available",
        "labels",
        "annotated_v_count",
        "total_v_count",
        "derivation",
        "in_content_hash",
    }, "v_subregion_support shape regressed; Slice 1 surface drifted."
    assert sup["labels"] == ["FWR1", "CDR1", "FWR2", "CDR2", "FWR3"]
    assert sup["in_content_hash"] is True
    assert sup["derivation"] == "bridge_imgt_gapped_seq"
    # Bundled HUMAN_IGH_OGRDB has full IMGT coverage; both counters
    # must be > 0 and equal (audited 198/198).
    assert sup["total_v_count"] > 0
    assert sup["annotated_v_count"] == sup["total_v_count"]


# ──────────────────────────────────────────────────────────────────
# 10. Absence — no v_subregion_rates kwarg on mutate
# ──────────────────────────────────────────────────────────────────


def test_pin_present_v_subregion_rates_kwarg_on_mutate() -> None:
    """Audit §5 / §11 Slice 2 — landed. ``Experiment.mutate``
    now carries the ``v_subregion_rates`` kwarg (Slice B). The
    other listed names remain forbidden — only the canonical
    spelling is the user-facing surface."""
    sig = inspect.signature(ga.Experiment.mutate)
    assert "v_subregion_rates" in sig.parameters, (
        "Experiment.mutate is missing the v_subregion_rates kwarg; "
        "Slice B surface regressed."
    )
    for forbidden in (
        "subregion_rates",
        "cdr_rates",
        "fr_rates",
        "fwr_rates",
    ):
        assert forbidden not in sig.parameters, (
            f"Experiment.mutate now accepts {forbidden!r}; only "
            "``v_subregion_rates`` is the canonical Slice B surface."
        )


# ──────────────────────────────────────────────────────────────────
# 11. Absence — no CDR/FR mutation counter AIRR fields
# ──────────────────────────────────────────────────────────────────


def test_pin_present_v_subregion_mutation_counter_fields() -> None:
    """Audit §11 Slice 3 — landed. Five canonical V-subregion
    counter fields plus ``n_v_unannotated_mutations`` partition
    ``n_v_mutations`` on every AIRR record. The two-bucket
    aggregates (``n_cdr_mutations`` / ``n_fwr_mutations``) remain
    deliberately absent — the counters audit §1 recommended
    five-label fields only, downstream aggregation is trivial."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=5)
        .run_records(n=1, seed=0)
    )
    rec = result.records[0]
    for required in (
        "n_cdr1_mutations",
        "n_cdr2_mutations",
        "n_fwr1_mutations",
        "n_fwr2_mutations",
        "n_fwr3_mutations",
        "n_v_unannotated_mutations",
    ):
        assert required in rec, (
            f"AIRR record missing {required!r}; V-subregion counters "
            "slice surface regressed"
        )
    # Two-bucket aggregates explicitly not in v1.
    for forbidden in ("n_cdr_mutations", "n_fwr_mutations"):
        assert forbidden not in rec, (
            f"AIRR record now carries the two-bucket aggregate "
            f"{forbidden!r}; the counters audit §1 recommended "
            "five-label fields only"
        )


# ──────────────────────────────────────────────────────────────────
# 12. Absence — no subregion mismatch validator issue kinds
# ──────────────────────────────────────────────────────────────────


def test_pin_present_v_subregion_mutation_count_mismatch_validator_kinds() -> None:
    """Audit §6 — the counters slice landed. The Rust validator
    now carries the six per-field
    ``N{Fwr1,Cdr1,Fwr2,Cdr2,Fwr3,VUnannotated}MutationsMismatch``
    issue kinds plus the cross-field
    ``VSubregionMutationCountSumMismatch`` invariant. Two of the
    placeholder names from the pre-implementation audit
    (``SubregionsBoundariesMismatch`` and ``SubregionSumMismatch``)
    were renamed in the implementation: the former is an
    AIRR-record validation kind that didn't ship (the slice's
    validation is event-time, not coordinate-boundary), and the
    latter became ``VSubregionMutationCountSumMismatch`` per the
    counters audit §5."""
    validate_src = (
        _REPO_ROOT / "engine_rs" / "src" / "airr_record" / "validate.rs"
    ).read_text(encoding="utf-8")
    # The six per-field mismatch kinds are now present.
    for required in (
        "NCdr1MutationsMismatch",
        "NCdr2MutationsMismatch",
        "NFwr1MutationsMismatch",
        "NFwr2MutationsMismatch",
        "NFwr3MutationsMismatch",
        "NVUnannotatedMutationsMismatch",
        "VSubregionMutationCountSumMismatch",
    ):
        assert required in validate_src, (
            f"validate.rs is missing {required!r}; V-subregion "
            "counter validator surface regressed"
        )
    # Placeholder names from the pre-implementation audit that
    # didn't ship as-is must still be absent — the slice used the
    # canonical names listed above.
    for placeholder in (
        "SubregionsBoundariesMismatch",
        "SubregionSumMismatch",
    ):
        assert placeholder not in validate_src, (
            f"validate.rs carries placeholder name {placeholder!r}; "
            "the canonical names are listed in this test's docstring"
        )


# ──────────────────────────────────────────────────────────────────
# 13. Absence — imgt_regions helper not consumed by simulation pipeline
# ──────────────────────────────────────────────────────────────────


def test_pin_present_imgt_regions_consumed_by_bridge_and_mcp() -> None:
    """Audit §1 / Slice 1: the ``imgt_regions`` derivation helper
    is now consumed by both the MCP analysis tooling (its
    original consumer) AND the cartridge bridge
    (``_refdata_resolver`` derives subregions at bridge time).
    Any other Python module pulling it in means a new surface
    has appeared and needs a pin update."""
    import subprocess

    result = subprocess.run(
        ["grep", "-rln", "from .*imgt_regions\\|import imgt_regions",
         "src/GenAIRR/"],
        capture_output=True,
        text=True,
        cwd=str(_REPO_ROOT),
    )
    consumers = sorted(set(result.stdout.strip().splitlines()))
    expected = {
        "src/GenAIRR/_refdata_resolver.py",
        "src/GenAIRR/utilities/mcp_helpers.py",
        # `ReferenceCartridgeBuilder.infer_v_subregions` uses the
        # same derivation helper at authoring time — same
        # boundary the bridge resolver uses at load time. See
        # `docs/reference_cartridge_authoring_audit.md` §2.
        "src/GenAIRR/cartridge_builder.py",
    }
    assert set(consumers) == expected, (
        f"imgt_regions consumers changed: got {consumers}, expected "
        f"{sorted(expected)}. A new consumer means an additional "
        "surface has appeared — flip pin and document the surface."
    )


# ──────────────────────────────────────────────────────────────────
# 14. Doc anchor — audit doc references contract file
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """The audit doc must continue to exist and reference this
    contract file; the 14-section structure stays intact."""
    docs_dir = Path(__file__).resolve().parent.parent / "docs"
    if not docs_dir.is_dir():
        import pytest
        pytest.skip("docs/ is contributor-only; not present in this checkout")
    doc_path = _REPO_ROOT / "docs" / "v_region_substructure_audit.md"
    assert doc_path.exists(), "v_region_substructure_audit.md missing"
    doc = doc_path.read_text(encoding="utf-8")
    assert "test_v_region_substructure_contract.py" in doc, (
        "audit doc no longer references the contract file; lockstep "
        "convention drifted."
    )
    for marker in (
        "## 1. Q1",
        "## 4. Q4",
        "## 11. Implementation order",
        "## 14. Summary table",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
