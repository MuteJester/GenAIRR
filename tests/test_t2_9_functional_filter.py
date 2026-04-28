"""T2-9: F/ORF/P/partial filter policy for custom DataConfigs.

Audit ask:
  1. Add gene_status to Allele                       — done in T2-8 closure
  2. Default productive sims to F-only               — addressed via flags
  3. Expose include_orf / include_pseudo flags       — landed here

Covers:
  * Loaders surface every record with functional_status (no
    loader-side filtering anymore)
  * Default policy keeps F + ORF; drops pseudogenes + partials
  * include_orf=False → F only
  * include_pseudo=True → keeps pseudogenes (with keep_anchorless=True
    since their anchors typically resolve to REJECTED)
  * include_partial=True → keeps partial sequences
  * BuildReport tracks filter exclusions separately from anchor
    rejections (different semantic categories)
"""
from __future__ import annotations

import pytest

from GenAIRR._native._anchor import (
    FunctionalStatus, Locus, ReferenceLoader, Segment,
)
from GenAIRR.dataconfig.enums import ChainType, Species
from GenAIRR.dataconfig.make.random import RandomDataConfigBuilder


# ── Synthetic IMGT FASTA covering all four statuses ────────────────


@pytest.fixture
def mixed_status_v_fasta(tmp_path):
    """V FASTA with 4 alleles: 1 F, 1 ORF, 1 pseudogene, 1 partial.
    All have valid Cys at IMGT pos 309 so the anchor resolver
    succeeds for any record that passes the filter."""
    f = tmp_path / "v.fasta"
    body = (
        ">X|IGHV1-1*01|Homo sapiens|F|V-REGION|...\n"
        + "a" * 309 + "tgt" + "a" * 30 + "\n"
        ">X|IGHV1-2*01|Homo sapiens|ORF|V-REGION|...\n"
        + "a" * 309 + "tgt" + "a" * 30 + "\n"
        ">X|IGHV1-3*01|Homo sapiens|P|V-REGION|...\n"
        + "a" * 309 + "tgt" + "a" * 30 + "\n"
        ">X|IGHV1-4*01|Homo sapiens|F|V-REGION|...|partial\n"
        + "a" * 309 + "tgt" + "a" * 30 + "\n"
    )
    f.write_text(body)
    return str(f)


@pytest.fixture
def simple_j_fasta(tmp_path):
    f = tmp_path / "j.fasta"
    f.write_text(
        ">X|IGHJ1*01|Homo sapiens|F|J-REGION|...\n"
        "tggggcaaaggcacc\n"
    )
    return str(f)


def _v_count(cfg) -> int:
    return sum(len(lst) for lst in cfg.v_alleles.values())


# ── Loader pass-through: every record surfaces with status tag ─────


class TestLoaderPassThrough:
    """T2-9 architectural change: loaders are now dumb pass-through.
    Filtering happens in Python, not C."""

    def test_imgt_loader_surfaces_partial_and_pseudo(self, mixed_status_v_fasta):
        with ReferenceLoader.open_imgt(mixed_status_v_fasta) as L:
            recs = list(L)
        # All 4 records surfaced; statuses preserved.
        statuses = [r.functional_status for r in recs]
        assert statuses.count(FunctionalStatus.F) == 1
        assert statuses.count(FunctionalStatus.ORF) == 1
        assert statuses.count(FunctionalStatus.PSEUDO) == 1
        assert statuses.count(FunctionalStatus.PARTIAL) == 1


# ── Default policy: F + ORF, no pseudo/partial ─────────────────────


class TestDefaultFilterPolicy:

    def test_default_keeps_f_and_orf_drops_pseudo_and_partial(
            self, mixed_status_v_fasta, simple_j_fasta):
        builder = RandomDataConfigBuilder(
            species=Species.HUMAN, chain_type=ChainType.BCR_HEAVY)
        cfg, report = builder.make_from_reference(
            v_reference=mixed_status_v_fasta,
            j_reference=simple_j_fasta,
            v_format="imgt", j_format="imgt",
        )
        # F + ORF kept; P + partial dropped.
        assert _v_count(cfg) == 2
        # Filtered names recorded.
        filtered_names = {n for n, _ in report.filtered_out}
        assert "IGHV1-3*01" in filtered_names    # pseudo
        assert "IGHV1-4*01" in filtered_names    # partial

    def test_filter_status_buckets_populated(
            self, mixed_status_v_fasta, simple_j_fasta):
        builder = RandomDataConfigBuilder(
            species=Species.HUMAN, chain_type=ChainType.BCR_HEAVY)
        _, report = builder.make_from_reference(
            v_reference=mixed_status_v_fasta,
            j_reference=simple_j_fasta,
            v_format="imgt", j_format="imgt",
        )
        assert report.by_filter_status["PSEUDO"] == 1
        assert report.by_filter_status["PARTIAL"] == 1


# ── F-only policy: include_orf=False ───────────────────────────────


class TestFOnlyPolicy:

    def test_f_only_keeps_one_v(self, mixed_status_v_fasta, simple_j_fasta):
        builder = RandomDataConfigBuilder(
            species=Species.HUMAN, chain_type=ChainType.BCR_HEAVY)
        cfg, report = builder.make_from_reference(
            v_reference=mixed_status_v_fasta,
            j_reference=simple_j_fasta,
            v_format="imgt", j_format="imgt",
            include_orf=False,
        )
        assert _v_count(cfg) == 1
        # ORF + P + partial all filtered out.
        assert report.by_filter_status["ORF"] == 1
        assert report.by_filter_status["PSEUDO"] == 1
        assert report.by_filter_status["PARTIAL"] == 1


# ── Pseudogene-keeping policy ──────────────────────────────────────


class TestPseudogenePolicy:

    def test_include_pseudo_keeps_them(
            self, mixed_status_v_fasta, simple_j_fasta):
        builder = RandomDataConfigBuilder(
            species=Species.HUMAN, chain_type=ChainType.BCR_HEAVY)
        cfg, report = builder.make_from_reference(
            v_reference=mixed_status_v_fasta,
            j_reference=simple_j_fasta,
            v_format="imgt", j_format="imgt",
            include_pseudo=True,
        )
        # F + ORF + P kept; partial still excluded.
        assert _v_count(cfg) == 3
        # Only partial filtered.
        assert "IGHV1-4*01" in {n for n, _ in report.filtered_out}
        assert "IGHV1-3*01" not in {n for n, _ in report.filtered_out}

    def test_pseudo_alleles_carry_pseudo_status(
            self, mixed_status_v_fasta, simple_j_fasta):
        """Verifies the T2-8 closure (Allele.functional_status) plays
        cleanly with T2-9: the pseudogene allele is in the config AND
        carries its status, so downstream code can identify it."""
        builder = RandomDataConfigBuilder(
            species=Species.HUMAN, chain_type=ChainType.BCR_HEAVY)
        cfg, _ = builder.make_from_reference(
            v_reference=mixed_status_v_fasta,
            j_reference=simple_j_fasta,
            v_format="imgt", j_format="imgt",
            include_pseudo=True,
        )
        # Find the pseudogene allele.
        all_v = [a for lst in cfg.v_alleles.values() for a in lst]
        pseudos = [a for a in all_v
                   if a.functional_status == FunctionalStatus.PSEUDO]
        assert len(pseudos) == 1
        assert pseudos[0].name == "IGHV1-3*01"


# ── Partial-keeping policy ─────────────────────────────────────────


class TestPartialPolicy:

    def test_include_partial_keeps_them(
            self, mixed_status_v_fasta, simple_j_fasta):
        builder = RandomDataConfigBuilder(
            species=Species.HUMAN, chain_type=ChainType.BCR_HEAVY)
        cfg, report = builder.make_from_reference(
            v_reference=mixed_status_v_fasta,
            j_reference=simple_j_fasta,
            v_format="imgt", j_format="imgt",
            include_partial=True,
        )
        # F + ORF + partial kept; P still excluded.
        assert _v_count(cfg) == 3
        assert "IGHV1-3*01" in {n for n, _ in report.filtered_out}


# ── Permissive policy: include everything ──────────────────────────


class TestAllInclusive:

    def test_include_all_keeps_every_record(
            self, mixed_status_v_fasta, simple_j_fasta):
        builder = RandomDataConfigBuilder(
            species=Species.HUMAN, chain_type=ChainType.BCR_HEAVY)
        cfg, report = builder.make_from_reference(
            v_reference=mixed_status_v_fasta,
            j_reference=simple_j_fasta,
            v_format="imgt", j_format="imgt",
            include_orf=True,
            include_pseudo=True,
            include_partial=True,
        )
        assert _v_count(cfg) == 4
        assert report.filtered_out == []


# ── BuildReport summary line includes filter info ──────────────────


class TestBuildReportSummary:

    def test_summary_mentions_filtered_count(
            self, mixed_status_v_fasta, simple_j_fasta):
        builder = RandomDataConfigBuilder(
            species=Species.HUMAN, chain_type=ChainType.BCR_HEAVY)
        _, report = builder.make_from_reference(
            v_reference=mixed_status_v_fasta,
            j_reference=simple_j_fasta,
            v_format="imgt", j_format="imgt",
        )
        s = report.summary()
        assert "filtered" in s.lower()
        assert "PSEUDO" in s
        assert "PARTIAL" in s


# ── Regression: shipped builtins still load ────────────────────────


class TestShippedBuiltinsRegression:

    def test_shipped_builtin_still_loads(self):
        """Removing loader-side filtering shouldn't affect the
        existing builtin pickles — they were pre-built and don't go
        through the loader path."""
        from GenAIRR.data import HUMAN_IGH_OGRDB
        HUMAN_IGH_OGRDB.verify_integrity()
        assert HUMAN_IGH_OGRDB.v_alleles is not None

    def test_simulation_still_runs(self):
        from GenAIRR import Experiment
        results = list(Experiment.on("human_igh").run(n=2, seed=42))
        assert len(results) == 2
