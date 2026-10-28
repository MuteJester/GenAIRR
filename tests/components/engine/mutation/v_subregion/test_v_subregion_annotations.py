"""Spec tests for **Slice 1: V-Subregion Cartridge Annotation
Surface**.

This slice makes V-region substructure annotations a first-class
cartridge property — derived (or user-supplied), validated,
hashed, and exposed via the cartridge manifest. It does NOT
introduce CDR/FR mutation counters, ``v_subregion_rates`` on
``Experiment.mutate``, or any SHM behaviour change; those gaps
remain explicitly pinned by
``tests/test_v_region_substructure_contract.py``.

Spec coverage (from the slice brief):

1. Bundled human IGH/IGK/IGL all derive subregions for every
   V allele (full coverage), via the bridge.
2. Intervals are monotonic, in-bounds, non-overlapping, and use
   the five canonical labels FWR1/CDR1/FWR2/CDR2/FWR3.
3. Editing a single subregion interval changes
   ``refdata_content_hash``.
4. The cartridge manifest's ``v_subregion_support`` block
   reports coverage accurately.
5. A legacy allele with no ``gapped_seq`` and no explicit
   ``subregions`` still loads with zero coverage (absence
   valid, not fatal).
6. A user-supplied explicit ``subregions`` dict survives the
   bridge round-trip.
7. Malformed user subregions (overlapping intervals, out-of-bounds,
   start>=end, unknown label, duplicate label) fail bridge
   validation with a clear error.

Out of scope: no Rust simulation behaviour change, no SHM-rate
plumbing, no CDR/FR counters.
"""
from __future__ import annotations

import copy

import pytest

import GenAIRR as ga
from GenAIRR._refdata_resolver import (
    _resolve_v_subregions,
    dataconfig_to_refdata,
)

CANONICAL_LABELS = ("FWR1", "CDR1", "FWR2", "CDR2", "FWR3")


# ──────────────────────────────────────────────────────────────────
# Spec 1 — bundled human cartridges derive subregions for every V
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "name",
    ["HUMAN_IGH_OGRDB", "HUMAN_IGK_OGRDB", "HUMAN_IGL_OGRDB"],
)
def test_bundled_human_full_v_subregion_coverage(name: str) -> None:
    """Every V allele in the three bundled human cartridges gets a
    derived subregion list via the bridge (IMGT ``gapped_seq``
    boundaries → ungapped allele coordinates). Coverage is
    ``100%``; the slice's bridge wiring is uniform across the
    three loci."""
    cfg = getattr(ga, name)
    sup = cfg.cartridge_manifest()["models"]["shm"]["v_subregion_support"]
    assert sup["available"] is True
    assert sup["annotated_v_count"] == sup["total_v_count"]
    assert sup["annotated_v_count"] > 0


# ──────────────────────────────────────────────────────────────────
# Spec 2 — intervals are monotonic, in-bounds, non-overlapping,
#          labelled with the canonical IMGT strings
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "name",
    ["HUMAN_IGH_OGRDB", "HUMAN_IGK_OGRDB", "HUMAN_IGL_OGRDB"],
)
def test_subregions_are_monotonic_and_in_bounds(name: str) -> None:
    """Every V allele's subregion list, walked in ``label``
    iteration order, is sorted ascending by ``start``, every
    interval is strictly half-open (``start < end``), and the
    whole list stays inside ``[0, len(ungapped_seq))``. Labels
    are exactly the five canonical IMGT strings."""
    cfg = getattr(ga, name)
    rd = dataconfig_to_refdata(cfg)
    for v_id in range(rd.v_pool_size()):
        allele = rd.v_allele(v_id)
        subs = allele.subregions
        if not subs:
            continue
        seq_len = len(allele.seq())
        prev_end = 0
        seen_labels = set()
        for label, start, end in subs:
            assert label in CANONICAL_LABELS, (
                f"{allele.name}: unknown label {label}"
            )
            assert label not in seen_labels, (
                f"{allele.name}: duplicate label {label}"
            )
            seen_labels.add(label)
            assert start < end, (
                f"{allele.name}: empty interval {label}={start}-{end}"
            )
            assert end <= seq_len, (
                f"{allele.name}: out-of-bounds {label}={start}-{end} "
                f"(len={seq_len})"
            )
            assert start >= prev_end, (
                f"{allele.name}: non-monotonic {label}={start}-{end} "
                f"after end={prev_end}"
            )
            prev_end = end


# ──────────────────────────────────────────────────────────────────
# Spec 3 — interval edits change refdata_content_hash
# ──────────────────────────────────────────────────────────────────


def test_subregion_interval_change_flips_content_hash() -> None:
    """Mutating a single subregion interval on a single V allele
    flips ``refdata_content_hash`` — the slice promised that
    subregion intervals enter the hash (manifest reports
    ``in_content_hash=True``). Same allele set, same anchors,
    same functional status, only the interval differs."""
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    base_hash = dataconfig_to_refdata(cfg).content_hash()

    # Override one V allele's subregions; pick the first one.
    first_gene_alleles = next(iter(cfg.v_alleles.values()))
    target = first_gene_alleles[0]
    seq_len = len(target.ungapped_seq)
    # Build a subregion list that's structurally valid but
    # definitely different from whatever derivation produced.
    target.subregions = {
        "FWR1": (0, min(60, seq_len // 4)),
        "CDR1": (min(60, seq_len // 4), min(78, seq_len // 3)),
        "FWR2": (min(78, seq_len // 3), min(120, seq_len // 2)),
    }
    edited_hash = dataconfig_to_refdata(cfg).content_hash()
    assert edited_hash != base_hash


# ──────────────────────────────────────────────────────────────────
# Spec 4 — manifest reports coverage accurately
# ──────────────────────────────────────────────────────────────────


def test_manifest_v_subregion_support_shape() -> None:
    """The ``v_subregion_support`` block carries exactly the six
    documented keys (``available``, ``labels``,
    ``annotated_v_count``, ``total_v_count``, ``derivation``,
    ``in_content_hash``) — locked here so additions go through
    a documented manifest update."""
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
    }
    assert sup["labels"] == list(CANONICAL_LABELS)
    assert sup["derivation"] == "bridge_imgt_gapped_seq"
    assert sup["in_content_hash"] is True


def test_manifest_total_v_count_matches_dataconfig_count() -> None:
    """``v_subregion_support.total_v_count`` equals the bridged
    V pool size and ``DataConfig.number_of_v_alleles``."""
    cfg = ga.HUMAN_IGH_OGRDB
    sup = cfg.cartridge_manifest()["models"]["shm"]["v_subregion_support"]
    assert sup["total_v_count"] == cfg.number_of_v_alleles


# ──────────────────────────────────────────────────────────────────
# Spec 5 — legacy allele with no gapped_seq → zero coverage but
#          NOT a load failure
# ──────────────────────────────────────────────────────────────────


def test_legacy_allele_without_gapped_seq_loads_with_zero_coverage() -> None:
    """An allele with empty ``gapped_seq`` and no user-supplied
    ``subregions`` loads cleanly — the bridge resolver returns
    an empty list, the engine accepts it, and the manifest
    counts it as unannotated. Absence of metadata is valid."""
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    # Strip gapped_seq from one V allele and clear any subregion
    # override.
    target = next(iter(cfg.v_alleles.values()))[0]
    target.gapped_seq = ""
    target.subregions = None

    # Bridge must succeed.
    rd = dataconfig_to_refdata(cfg)
    # The targeted allele's subregions are now empty…
    for v_id in range(rd.v_pool_size()):
        a = rd.v_allele(v_id)
        if a.name == target.name:
            assert a.subregions == []
            break
    else:
        pytest.fail(f"target allele {target.name} not found in bridged pool")

    # …and the manifest reports the coverage drop.
    sup = cfg.cartridge_manifest()["models"]["shm"]["v_subregion_support"]
    assert sup["annotated_v_count"] == sup["total_v_count"] - 1
    assert sup["annotated_v_count"] >= 0


# ──────────────────────────────────────────────────────────────────
# Spec 6 — user-supplied subregions survive the bridge round-trip
# ──────────────────────────────────────────────────────────────────


def test_explicit_user_subregions_round_trip() -> None:
    """When the user sets ``allele.subregions`` directly with a
    valid dict ``{label: (start, end)}``, the bridge forwards
    those intervals verbatim — the Rust ``PyAllele.subregions``
    accessor surfaces them unchanged."""
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    target = next(iter(cfg.v_alleles.values()))[0]
    seq_len = len(target.ungapped_seq)
    expected = {
        "FWR1": (0, 12),
        "CDR1": (12, 24),
        "FWR2": (24, min(36, seq_len)),
    }
    target.subregions = expected

    rd = dataconfig_to_refdata(cfg)
    found = None
    for v_id in range(rd.v_pool_size()):
        a = rd.v_allele(v_id)
        if a.name == target.name:
            found = a.subregions
            break
    assert found is not None
    # Order is monotonic by ``start`` regardless of dict iteration.
    assert found == [
        ("FWR1", 0, 12),
        ("CDR1", 12, 24),
        ("FWR2", 24, min(36, seq_len)),
    ]


# ──────────────────────────────────────────────────────────────────
# Spec 7 — malformed user subregions reject at the bridge
# ──────────────────────────────────────────────────────────────────


def _cfg_with_v_subregions(subregions: dict):
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    target = next(iter(cfg.v_alleles.values()))[0]
    target.subregions = subregions
    return cfg


def test_malformed_subregions_overlap_rejected() -> None:
    """Overlapping intervals (CDR1 end > FWR2 start) raise at the
    bridge. The Rust validator (``parse_subregions``) walks the
    sorted list and rejects the pair."""
    bad = {
        "FWR1": (0, 30),
        "CDR1": (30, 60),
        "FWR2": (50, 90),  # overlaps CDR1
    }
    with pytest.raises(Exception):
        dataconfig_to_refdata(_cfg_with_v_subregions(bad))


def test_malformed_subregions_out_of_bounds_rejected() -> None:
    """An ``end`` past the ungapped sequence length raises."""
    # 100000 is unrealistically far past any V allele length.
    bad = {"FWR1": (0, 100000)}
    with pytest.raises(Exception):
        dataconfig_to_refdata(_cfg_with_v_subregions(bad))


def test_malformed_subregions_start_ge_end_rejected() -> None:
    """``start >= end`` raises — half-open intervals must be
    non-empty."""
    bad = {"FWR1": (10, 10)}
    with pytest.raises(Exception):
        dataconfig_to_refdata(_cfg_with_v_subregions(bad))


def test_malformed_subregions_unknown_label_rejected() -> None:
    """A label that isn't one of the five canonical IMGT strings
    raises. Case-sensitive."""
    bad = {"fwr1": (0, 30)}  # lowercase
    with pytest.raises(Exception):
        dataconfig_to_refdata(_cfg_with_v_subregions(bad))


def test_malformed_subregions_duplicate_label_rejected() -> None:
    """Duplicate label entries — even with disjoint intervals —
    raise at the bridge: each label appears at most once."""
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    target = next(iter(cfg.v_alleles.values()))[0]
    # Use a list-of-tuples override to bypass the dict
    # uniqueness; this is the wire shape the Rust side validates.
    target.subregions = None
    # We can't actually express duplicates through the dict-style
    # ``subregions`` attribute, so we verify the wire-format
    # validator directly via the resolver helper.
    from GenAIRR import _engine

    rd = _engine.RefDataConfig("vdj")
    with pytest.raises(Exception):
        rd.add_v_allele(
            "X*01",
            "X",
            b"A" * 300,
            anchor=None,
            functional_status=None,
            subregions=[
                ("FWR1", 0, 30),
                ("FWR1", 30, 60),  # duplicate label
            ],
        )


# ──────────────────────────────────────────────────────────────────
# Spec extras — resolver helper returns [] in safe-fail cases
# ──────────────────────────────────────────────────────────────────


def test_resolver_returns_empty_on_missing_gapped_seq() -> None:
    """Calling the resolver helper on an allele with no
    ``gapped_seq`` and no override returns ``[]`` — never raises."""
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    target = next(iter(cfg.v_alleles.values()))[0]
    target.gapped_seq = ""
    target.subregions = None
    assert _resolve_v_subregions(target) == []
