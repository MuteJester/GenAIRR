"""T2-8 closure tests — verify the full custom-config user journey:

  1. Build a DataConfig from custom references via `make_from_reference`.
  2. Every Allele carries the metadata that came in via the loader
     (functional_status, locus, source, aliases, species).
  3. The DataConfig has `schema_sha256` populated; `verify_integrity()`
     passes.
  4. `build_report` is attached to the DataConfig and survives pickle
     round-trip.
  5. Mutating `build_report` post-load does NOT invalidate the
     integrity check (build_report is excluded from the checksum).
  6. Every shipped builtin still verifies cleanly (regression).
"""
from __future__ import annotations

import json
import pickle
from pathlib import Path

import pytest

from GenAIRR._native._anchor import (
    AnchorConfidence, FunctionalStatus, Locus, Segment,
)
from GenAIRR.dataconfig.enums import ChainType, Species
from GenAIRR.dataconfig.make.random import BuildReport, RandomDataConfigBuilder


# ── Fixtures ───────────────────────────────────────────────────────


@pytest.fixture
def imgt_v_fasta(tmp_path):
    f = tmp_path / "v.fasta"
    f.write_text(
        ">X|IGHV1-1*01|Homo sapiens|F|V-REGION|...\n"
        + "a" * 309 + "tgt" + "a" * 30 + "\n"
    )
    return str(f)


@pytest.fixture
def imgt_j_fasta(tmp_path):
    f = tmp_path / "j.fasta"
    f.write_text(
        ">X|IGHJ1*01|Homo sapiens|F|J-REGION|...\n"
        "tggggcaaaggcacc\n"
    )
    return str(f)


@pytest.fixture
def imgt_d_fasta(tmp_path):
    f = tmp_path / "d.fasta"
    f.write_text(">X|IGHD1*01|Homo sapiens|F|D-REGION|...\nacgtacgt\n")
    return str(f)


@pytest.fixture
def custom_cfg_with_report(imgt_v_fasta, imgt_j_fasta, imgt_d_fasta):
    builder = RandomDataConfigBuilder(
        species=Species.HUMAN,
        chain_type=ChainType.BCR_HEAVY,
        reference_set="custom-test",
    )
    cfg, report = builder.make_from_reference(
        v_reference=imgt_v_fasta,
        j_reference=imgt_j_fasta,
        d_reference=imgt_d_fasta,
        v_format="imgt", j_format="imgt",
    )
    return cfg, report


# ── Gap 1 — Allele metadata preservation ───────────────────────────


class TestAlleleMetadataPreserved:
    """Per-allele fields that came in via the LoadedAlleleRecord must
    survive into the constructed Allele instance — previously dropped."""

    def test_v_allele_carries_functional_status(self, custom_cfg_with_report):
        cfg, _ = custom_cfg_with_report
        a = list(cfg.v_alleles.values())[0][0]
        assert a.functional_status is FunctionalStatus.F

    def test_v_allele_carries_locus(self, custom_cfg_with_report):
        cfg, _ = custom_cfg_with_report
        a = list(cfg.v_alleles.values())[0][0]
        assert a.locus is Locus.IGH

    def test_v_allele_carries_source(self, custom_cfg_with_report):
        cfg, _ = custom_cfg_with_report
        a = list(cfg.v_alleles.values())[0][0]
        assert a.source == "imgt-vquest"

    def test_v_allele_carries_species(self, custom_cfg_with_report):
        cfg, _ = custom_cfg_with_report
        a = list(cfg.v_alleles.values())[0][0]
        assert a.species == "Homo sapiens"

    def test_v_allele_aliases_default_to_empty_tuple(self, custom_cfg_with_report):
        # IMGT FASTAs don't carry aliases — should be empty tuple.
        cfg, _ = custom_cfg_with_report
        a = list(cfg.v_alleles.values())[0][0]
        assert a.aliases == ()

    def test_d_allele_also_carries_metadata(self, custom_cfg_with_report):
        cfg, _ = custom_cfg_with_report
        d = list(cfg.d_alleles.values())[0][0]
        assert d.locus is Locus.IGH
        assert d.source == "imgt-vquest"
        assert d.functional_status is FunctionalStatus.F

    def test_j_allele_also_carries_metadata(self, custom_cfg_with_report):
        cfg, _ = custom_cfg_with_report
        j = list(cfg.j_alleles.values())[0][0]
        assert j.locus is Locus.IGH
        assert j.source == "imgt-vquest"


# ── Gap 1+ — Aliases from OGRDB sidecar are surfaced ───────────────


class TestOgrdbAliasesPreserved:

    def test_aliases_round_trip_from_ogrdb(self, tmp_path):
        fasta = tmp_path / "v.fasta"
        fasta.write_text(
            ">IGHV1-1*01\n" + "a" * 309 + "tgt" + "a" * 30 + "\n"
        )
        sidecar = tmp_path / "v.json"
        sidecar.write_text(json.dumps({
            "GermlineSet": [{
                "allele_descriptions": [{
                    "label": "IGHV1-1*01",
                    "functional": True,
                    "v_gene_delineations": [
                        {"delineation_scheme": "IMGT", "cdr3_start": 310},
                    ],
                    "sequence_alias": ["IGHV1-1*01_alt", "VH1-1*01"],
                }],
            }],
        }))
        j = tmp_path / "j.fasta"
        j.write_text(">IGHJ1*01\ntggggcaaaggcacc\n")
        builder = RandomDataConfigBuilder(
            species=Species.HUMAN,
            chain_type=ChainType.BCR_HEAVY,
        )
        cfg, _ = builder.make_from_reference(
            v_reference=str(fasta),
            j_reference=str(j),
            v_sidecar=str(sidecar),
            v_format="ogrdb", j_format="plain",
            j_locus_hint=Locus.IGH,
        )
        a = list(cfg.v_alleles.values())[0][0]
        # Aliases captured from OGRDB sequence_alias array.
        assert "IGHV1-1*01_alt" in a.aliases
        assert "VH1-1*01" in a.aliases


# ── Gap 2 — schema_sha256 auto-stamped ─────────────────────────────


class TestChecksumAutoStamp:

    def test_schema_sha256_populated(self, custom_cfg_with_report):
        cfg, _ = custom_cfg_with_report
        assert cfg.schema_sha256 != ""
        assert len(cfg.schema_sha256) == 64    # sha256 hex

    def test_verify_integrity_passes(self, custom_cfg_with_report):
        cfg, _ = custom_cfg_with_report
        cfg.verify_integrity()  # must not raise

    def test_round_trip_pickle_verifies(self, custom_cfg_with_report):
        cfg, _ = custom_cfg_with_report
        restored = pickle.loads(pickle.dumps(cfg))
        restored.verify_integrity()


# ── Gap 3 — BuildReport attached to DataConfig ─────────────────────


class TestBuildReportAttached:

    def test_report_attribute_set(self, custom_cfg_with_report):
        cfg, returned_report = custom_cfg_with_report
        assert cfg.build_report is returned_report

    def test_report_survives_pickle(self, custom_cfg_with_report):
        cfg, _ = custom_cfg_with_report
        restored = pickle.loads(pickle.dumps(cfg))
        assert restored.build_report is not None
        assert isinstance(restored.build_report, BuildReport)
        assert restored.build_report.accepted >= 3   # 1 V + 1 D + 1 J

    def test_report_summary_works(self, custom_cfg_with_report):
        cfg, _ = custom_cfg_with_report
        s = cfg.build_report.summary()
        assert "accepted" in s.lower()


# ── Build report does NOT invalidate the checksum ──────────────────


class TestBuildReportOrthogonalToChecksum:
    """The build_report carries diagnostic data that has no semantic
    effect on simulation. It must be excluded from the integrity hash
    so users can inspect/mutate it without breaking verification."""

    def test_mutating_build_report_keeps_integrity(self, custom_cfg_with_report):
        cfg, _ = custom_cfg_with_report
        cfg.verify_integrity()
        # Mutate a report field.
        cfg.build_report.notes.append("user added a note")
        # Integrity still passes — report is excluded from checksum.
        cfg.verify_integrity()

    def test_clearing_build_report_keeps_integrity(self, custom_cfg_with_report):
        cfg, _ = custom_cfg_with_report
        cfg.build_report = None
        cfg.verify_integrity()


# ── Regression: shipped builtins still verify cleanly ──────────────


class TestShippedBuiltinsRegression:
    """Adding `build_report` to DataConfig must not change the
    pre-computed checksums of the 106 shipped builtin pickles. Their
    __dict__ doesn't contain `build_report`; `compute_checksum`
    pops the key from __dict__ during hashing so legacy and new
    pickles produce identical hash inputs."""

    def test_sample_builtin_still_verifies(self):
        from GenAIRR.data import HUMAN_IGH_OGRDB
        HUMAN_IGH_OGRDB.verify_integrity()
        assert HUMAN_IGH_OGRDB.build_report is None  # legacy: no report

    def test_recompute_matches_stored_for_builtin(self):
        from GenAIRR.data import HUMAN_IGH_OGRDB
        # The stored sha256 must match what compute_checksum produces
        # against the same instance — even though instance lacks
        # `build_report` in __dict__.
        recomputed = HUMAN_IGH_OGRDB.compute_checksum()
        assert recomputed == HUMAN_IGH_OGRDB.schema_sha256


# ── Legacy pickle (pre-build_report field) loads cleanly ───────────


class TestLegacyAlleleClassDefaults:
    """Allele class-level defaults guarantee that legacy pickled
    instances (built before the new fields existed) still expose
    those fields without AttributeError on access."""

    def test_legacy_pickled_v_allele_has_all_new_fields(self):
        from GenAIRR.data import HUMAN_IGH_OGRDB
        a = list(HUMAN_IGH_OGRDB.v_alleles.values())[0][0]
        # All five new fields readable; values are class-level defaults
        # (None / empty tuple) since the legacy pickle didn't set them.
        assert a.functional_status is None
        assert a.locus is None
        assert a.aliases == ()
        assert a.species is None
        assert a.source is None
