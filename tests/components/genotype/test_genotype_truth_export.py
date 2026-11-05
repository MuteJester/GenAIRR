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
        """End to end for the decision that a deleted gene stays in the table.

        Built from a complete genotype, not a permissive one. A permissive
        genotype reaches the writer's refusal now, and it was the wrong
        premise before that: the deletion was real, but every unspecified
        gene was missing from the file the test then read, so it asserted
        the sentinel inside a table that was already a subset of the truth.
        """
        cfg = _cfg()
        gene = next(iter(cfg.v_alleles))
        g = (Genotype.from_dataconfig(cfg).complete_from_reference()
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


class TestAnIncompleteGenotypeRefuses:
    """A permissive genotype is not the set of alleles the records came from.

    `Genotype.permissive` leaves unspecified genes to be sampled over the
    whole reference, unphased -- its own docstring says never to mistake it
    for a complete genotype. Projecting one yields a row only for the genes
    somebody happened to constrain, so the file is silently a subset of the
    truth: measured on this cartridge, a permissive genotype with one gene
    pinned writes two rows while the records draw 119 alleles that are not in
    it. A benchmark scoring a caller against that file marks correct calls
    wrong, which is worse than having no file, and it is the same failure the
    empty-result refusal already guards.
    """

    def test_a_permissive_genotype_with_unspecified_genes_refuses(self, tmp_path):
        g = Genotype.permissive(_cfg()).with_subject("P0")
        result = _run(g)
        out = tmp_path / "genotype.tsv"
        with pytest.raises(ValueError, match="unspecified"):
            result.to_genotype_tsv(str(out))

    def test_the_refusal_says_how_to_complete_the_genotype(self, tmp_path):
        g = Genotype.permissive(_cfg()).with_subject("P0")
        result = _run(g)
        with pytest.raises(ValueError, match="complete_from_reference"):
            result.to_genotype_tsv(str(tmp_path / "genotype.tsv"))

    def test_the_refusal_leaves_no_file_behind(self, tmp_path):
        g = Genotype.permissive(_cfg()).with_subject("P0")
        result = _run(g)
        out = tmp_path / "genotype.tsv"
        with pytest.raises(ValueError):
            result.to_genotype_tsv(str(out))
        assert not out.exists()

    def test_one_pinned_gene_does_not_buy_a_sidecar(self, tmp_path):
        """The dangerous shape: enough specified to look deliberate."""
        cfg = _cfg()
        gene = next(iter(cfg.v_alleles))
        g = (Genotype.permissive(cfg)
             .homozygous(gene, cfg.v_alleles[gene][0].name)
             .with_subject("P1"))
        result = _run(g, n=40, seed=2)
        with pytest.raises(ValueError, match="unspecified"):
            result.to_genotype_tsv(str(tmp_path / "genotype.tsv"))

    def test_a_completed_permissive_genotype_is_accepted(self, tmp_path):
        """Completing it is what makes the table the carried set.

        `permissive` stays set, so the guard must test what is specified
        rather than which constructor was used -- refusing every permissive
        genotype would refuse this one, whose table is complete.
        """
        cfg = _cfg()
        g = (Genotype.permissive(cfg).complete_from_reference()
             .with_subject("P2"))
        result = _run(g, n=40, seed=3)
        out = tmp_path / "genotype.tsv"
        result.to_genotype_tsv(str(out))
        carried = {row["allele"] for row in _read(str(out))}
        assert g.is_permissive, "premise: this genotype is still permissive"
        drawn = {rec[f] for rec in result
                 for f in ("v_call", "d_call", "j_call") if rec.get(f)}
        assert drawn, "premise: the run drew some alleles"
        assert drawn <= carried, f"drawn but absent: {sorted(drawn - carried)}"


class TestTheFileIsJoinableAsAFile:
    """Byte-level properties a consumer joins on, not just parsed values.

    The projection tests assert Python objects; these assert what lands on
    disk. Four of the eight columns had no file-level assertion, including
    `novel` -- whose serialized spelling the documented join compares against
    the literal string "True". A writer that emitted 1/0 or true would make
    the published snippet silently return False for every read while every
    other test stayed green.
    """

    def test_the_novel_column_is_spelled_the_way_the_documented_join_reads_it(
            self, tmp_path):
        cfg = _cfg()
        gene = next(iter(cfg.v_alleles))
        base = cfg.v_alleles[gene][0].name
        g = (Genotype.from_dataconfig(cfg).complete_from_reference()
             .add_novel_allele(f"{gene}*i01", base=base, mutations=[(10, "A")],
                               allow_nonfunctional=True)
             .heterozygous(gene, base, f"{gene}*i01")
             .with_subject("N0"))
        out = tmp_path / "genotype.tsv"
        _run(g).to_genotype_tsv(str(out))

        rows = _read(str(out))
        spellings = {row["novel"] for row in rows}
        assert spellings <= {"True", "False"}, spellings
        novel_rows = [r for r in rows if r["novel"] == "True"]
        assert [r["allele"] for r in novel_rows] == [f"{gene}*i01"], novel_rows

    def test_every_column_carries_a_value_on_every_row(self, tmp_path):
        """`locus`, `segment` and `gene` had no file-level assertion either."""
        g = Genotype.sample(_cfg(), seed=0, subject_id="S0")
        out = tmp_path / "genotype.tsv"
        _run(g).to_genotype_tsv(str(out))
        rows = _read(str(out))
        assert rows, "premise: the file has rows"
        for column in COLUMNS:
            blank = [r for r in rows if r[column] == ""]
            assert blank == [], f"{column} is blank on {len(blank)} row(s)"
        assert {r["locus"] for r in rows} == {"IGH"}
        assert {r["segment"] for r in rows} == {"V", "D", "J"}

    def test_its_line_endings_match_the_records_file_it_joins_to(self, tmp_path):
        """`_write_delimited` pins LF deliberately. A CRLF sidecar beside an LF
        records file leaves a stray \\r on the last column, so `gene_state ==
        "deleted"` fails in any shell pipeline that joins the two."""
        g = Genotype.sample(_cfg(), seed=0, subject_id="S0")
        out = tmp_path / "genotype.tsv"
        _run(g).to_genotype_tsv(str(out))
        raw = out.read_bytes()
        assert b"\r\n" not in raw, "sidecar is CRLF; the records writer pins LF"


class TestAnUnjoinableCohortRefuses:
    def test_a_missing_genotype_is_refused_not_skipped(self, tmp_path):
        """Dropping a `None` writes a partial table that looks complete."""
        cfg = _cfg()
        g = Genotype.sample(cfg, seed=0, subject_id="S0")
        result = _run(g)
        cohort = CohortResult([
            CohortSubjectResult(subject_id="S0", genotype=g, result=result,
                                refdata=object(), seed=11, count=6),
            CohortSubjectResult(subject_id="S1", genotype=None, result=result,
                                refdata=object(), seed=22, count=6),
        ])
        with pytest.raises(ValueError, match="no genotype"):
            cohort.to_genotype_tsv(str(tmp_path / "genotype.tsv"))

    def test_two_subjects_sharing_one_id_are_refused(self, tmp_path):
        """Both subjects' rows would merge under one join key."""
        cfg = _cfg()
        g0 = Genotype.sample(cfg, seed=0, subject_id="DUP")
        g1 = Genotype.sample(cfg, seed=1, subject_id="DUP")
        result = _run(g0)
        cohort = CohortResult([
            CohortSubjectResult(subject_id="DUP", genotype=g0, result=result,
                                refdata=object(), seed=11, count=6),
            CohortSubjectResult(subject_id="DUP", genotype=g1, result=result,
                                refdata=object(), seed=22, count=6),
        ])
        with pytest.raises(ValueError, match="subject_id"):
            cohort.to_genotype_tsv(str(tmp_path / "genotype.tsv"))
