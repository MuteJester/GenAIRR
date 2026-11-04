"""`to_genotype_tsv` on a single-subject result and on a cohort.

A records file carries `subject_id` and nothing that says what that subject's
genome was. These tests pin the sidecar that closes the join, and they pin the
refusal harder than the happy path: a zero-row truth file written beside a
million records is the failure mode that looks like success, so a result with
no genotype must raise by name and leave no file behind.
"""
import csv

import pytest

import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR._genotype_truth import COLUMNS, DELETED_ALLELE, genotype_truth_rows
from GenAIRR.cohort import CohortResult, CohortSubjectResult
from GenAIRR.genotype import Genotype

pytestmark = pytest.mark.behavior


def _cfg():
    return gdata.HUMAN_IGH_OGRDB


def _read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _run(genotype, n=6, seed=1):
    return ga.Experiment.on(_cfg()).with_genotype(genotype).recombine().run_records(n=n, seed=seed)


class TestASingleSubjectSidecar:
    def test_its_subject_id_matches_every_record(self, tmp_path):
        g = Genotype.sample(_cfg(), seed=0, subject_id="S0")
        result = _run(g)
        out = tmp_path / "genotype.tsv"

        result.to_genotype_tsv(str(out))

        rows = _read(out)
        assert rows, "premise: the sidecar must not be empty"
        record_subjects = {r["subject_id"] for r in result}
        assert record_subjects == {"S0"}, "premise: the records must carry the subject"
        assert {r["subject_id"] for r in rows} == record_subjects

    def test_the_header_is_the_long_form_columns(self, tmp_path):
        g = Genotype.sample(_cfg(), seed=0, subject_id="S0")
        out = tmp_path / "genotype.tsv"

        _run(g).to_genotype_tsv(str(out))

        with open(out, newline="", encoding="utf-8") as fh:
            header = next(csv.reader(fh, delimiter="\t"))
        assert header == list(COLUMNS)

    def test_the_file_is_the_projection_and_nothing_else(self, tmp_path):
        g = Genotype.sample(_cfg(), seed=0, subject_id="S0")
        out = tmp_path / "genotype.tsv"

        _run(g).to_genotype_tsv(str(out))

        expected = genotype_truth_rows(g)
        rows = _read(out)
        assert len(rows) == len(expected)
        assert [r["allele"] for r in rows] == [r["allele"] for r in expected]
        assert [int(r["chromosome"]) for r in rows] == [r["chromosome"] for r in expected]

    def test_a_deleted_gene_reaches_the_file_as_the_sentinel(self, tmp_path):
        """End to end for the decision that a deleted gene stays in the table."""
        cfg = _cfg()
        gene = next(iter(cfg.v_alleles))
        allele = cfg.v_alleles[gene][0].name
        g = (Genotype.permissive(cfg)
             .homozygous(gene, allele)
             .delete_gene(gene, haplotype="both")
             .with_subject("S0"))
        out = tmp_path / "genotype.tsv"

        _run(g).to_genotype_tsv(str(out))

        deleted = [r for r in _read(out) if r["gene"] == gene]
        assert len(deleted) == 2, deleted
        assert {r["allele"] for r in deleted} == {DELETED_ALLELE}
        assert {r["gene_state"] for r in deleted} == {"deleted"}


class TestACohortSidecar:
    def test_it_carries_every_subject_and_nothing_else(self, tmp_path):
        cfg = _cfg()
        g0 = Genotype.sample(cfg, seed=0, subject_id="A")
        g1 = Genotype.sample(cfg, seed=1, subject_id="B")
        cohort = CohortResult([
            CohortSubjectResult(subject_id="A", genotype=g0, result=_run(g0, n=3, seed=1),
                                refdata=object(), seed=11, count=3),
            CohortSubjectResult(subject_id="B", genotype=g1, result=_run(g1, n=3, seed=2),
                                refdata=object(), seed=22, count=3),
        ])
        out = tmp_path / "cohort_genotype.tsv"

        cohort.to_genotype_tsv(str(out))

        rows = _read(out)
        assert {r["subject_id"] for r in rows} == {"A", "B"}
        assert len(rows) == len(genotype_truth_rows(g0)) + len(genotype_truth_rows(g1))

    def test_two_subjects_genotypes_are_not_collapsed(self, tmp_path):
        """Premise and point in one: the two sampled genotypes must actually
        differ, or the test could not tell a per-subject table from a
        duplicated one."""
        cfg = _cfg()
        g0 = Genotype.sample(cfg, seed=0, subject_id="A")
        g1 = Genotype.sample(cfg, seed=1, subject_id="B")
        a_truth = {(r["gene"], r["chromosome"], r["allele"]) for r in genotype_truth_rows(g0)}
        b_truth = {(r["gene"], r["chromosome"], r["allele"]) for r in genotype_truth_rows(g1)}
        assert a_truth != b_truth, "premise: two subjects must differ somewhere"

        cohort = CohortResult([
            CohortSubjectResult(subject_id="A", genotype=g0, result=_run(g0, n=3, seed=1),
                                refdata=object(), seed=11, count=3),
            CohortSubjectResult(subject_id="B", genotype=g1, result=_run(g1, n=3, seed=2),
                                refdata=object(), seed=22, count=3),
        ])
        out = tmp_path / "cohort_genotype.tsv"
        cohort.to_genotype_tsv(str(out))

        rows = _read(out)
        got_a = {(r["gene"], int(r["chromosome"]), r["allele"]) for r in rows
                 if r["subject_id"] == "A"}
        got_b = {(r["gene"], int(r["chromosome"]), r["allele"]) for r in rows
                 if r["subject_id"] == "B"}
        assert got_a != got_b


class TestARunWithNoGenotypeRefuses:
    """A zero-row truth file beside a million records is the failure that
    looks like success. It must raise, by name, and write nothing."""

    def test_a_plain_result_raises_naming_the_cause(self, tmp_path):
        result = ga.Experiment.on(_cfg()).recombine().run_records(n=4, seed=1)
        out = tmp_path / "genotype.tsv"

        with pytest.raises(ValueError, match="no genotype"):
            result.to_genotype_tsv(str(out))

    def test_the_refusal_leaves_no_file_behind(self, tmp_path):
        result = ga.Experiment.on(_cfg()).recombine().run_records(n=4, seed=1)
        out = tmp_path / "genotype.tsv"

        with pytest.raises(ValueError):
            result.to_genotype_tsv(str(out))

        assert not out.exists(), "a refused export must not leave a partial file"

    def test_a_genotype_with_no_subject_id_refuses(self, tmp_path):
        """A sidecar whose subject_id is blank cannot be joined to anything,
        so it is the same failure as an empty one: plausible and useless."""
        cfg = _cfg()
        gene = next(iter(cfg.v_alleles))
        nameless = Genotype.permissive(cfg).homozygous(gene, cfg.v_alleles[gene][0].name)
        assert not nameless.subject_id, "premise: the genotype must carry no subject_id"
        out = tmp_path / "genotype.tsv"

        with pytest.raises(ValueError, match="no subject_id"):
            _run(nameless).to_genotype_tsv(str(out))

        assert not out.exists()

    def test_a_cohort_with_no_genotypes_raises(self, tmp_path):
        plain = ga.Experiment.on(_cfg()).recombine().run_records(n=2, seed=1)
        cohort = CohortResult([
            CohortSubjectResult(subject_id="A", genotype=None, result=plain,
                                refdata=object(), seed=11, count=2),
        ])
        out = tmp_path / "cohort_genotype.tsv"

        with pytest.raises(ValueError, match="no genotype"):
            cohort.to_genotype_tsv(str(out))
        assert not out.exists()
