"""The test that makes the genotype sidecar trustworthy.

A sidecar is only ground truth if it accounts for every allele the records
were actually generated from. Measuring that once, by hand, proves it for one
run; this turns the measurement into a standing guarantee.

The contract: for every record, the allele named in `truth_v_call` /
`truth_d_call` / `truth_j_call` appears in that subject's sidecar rows as
carried on some chromosome. A sidecar that misses one is worse than absent --
a benchmark scoring a caller against it would count a correct call as wrong.

Every assertion here is paired with a guard on its own premise. A completeness
test runs green on a degenerate repertoire: one subject, one allele per
segment, no deletions. Those guards are what stop this from passing for the
wrong reason.
"""
import pytest

import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR._genotype_truth import DELETED_ALLELE, genotype_truth_rows
from GenAIRR.genotype import Genotype

pytestmark = pytest.mark.behavior

SUBJECTS = 3
PER_SUBJECT = 150


def _cfg():
    return gdata.HUMAN_IGH_OGRDB


@pytest.fixture(scope="module")
def cohort():
    """Three sampled subjects, one carrying a deliberately deleted gene.

    The deletion is forced rather than hoped for: whether `Genotype.sample`
    happens to delete a gene at these seeds is not something this test should
    depend on, and a sidecar that mishandles deletions is exactly the defect
    worth catching.
    """
    cfg = _cfg()
    genotypes = [Genotype.sample(cfg, seed=i, subject_id=f"S{i}") for i in range(SUBJECTS)]
    gene = next(iter(cfg.v_alleles))
    genotypes[-1] = genotypes[-1].delete_gene(gene, haplotype="both")
    return ga.Experiment.on(cfg).recombine().run_cohort(
        genotypes, n_per_subject=PER_SUBJECT, seed=7, expose_provenance=True)


def _carried_by_subject(cohort):
    carried = {}
    for subject in cohort.subjects:
        rows = genotype_truth_rows(subject.genotype)
        carried[subject.subject_id] = {
            r["allele"] for r in rows if r["allele"] != DELETED_ALLELE}
    return carried


class TestThePremisesThisRunMustSatisfy:
    """Guards first. If any of these fails, the completeness assertions below
    would be passing on a repertoire too degenerate to mean anything."""

    def test_the_run_produced_records_for_every_subject(self, cohort):
        counts = {s.subject_id: len(s.result) for s in cohort.subjects}
        assert len(counts) == SUBJECTS, counts
        assert all(n == PER_SUBJECT for n in counts.values()), counts

    def test_the_records_carry_the_truth_call_fields(self, cohort):
        rec = cohort.records[0]
        for field in ("truth_v_call", "truth_d_call", "truth_j_call", "subject_id"):
            assert field in rec, f"{field} missing; expose_provenance did not take effect"

    def test_more_than_one_allele_was_drawn_per_segment(self, cohort):
        """A run that drew a single V allele would satisfy completeness
        trivially."""
        for segment in ("v", "d", "j"):
            drawn = {r[f"truth_{segment}_call"] for r in cohort.records}
            drawn.discard("")
            assert len(drawn) > 1, f"only {drawn} drawn for {segment}; too degenerate"

    def test_at_least_one_subject_carries_a_fully_deleted_gene(self, cohort):
        """Fully deleted, not merely one chromosome short.

        `test_a_fully_deleted_gene_is_never_drawn_from` skips any subject
        with no gene deleted on BOTH chromosomes. A guard keyed on the
        `(none)` sentinel alone is satisfied by a hemizygous deletion, which
        would leave that test asserting nothing while this one still passed,
        so the guard has to match the state the test actually needs.
        """
        deleted = [
            (s.subject_id, r["gene"])
            for s in cohort.subjects
            for r in genotype_truth_rows(s.genotype)
            if r["allele"] == DELETED_ALLELE and r["gene_state"] == "deleted"
        ]
        assert deleted, "no fully deleted gene anywhere; that path is untested"

    def test_the_subjects_genotypes_are_not_identical(self, cohort):
        signatures = {
            s.subject_id: frozenset(
                (r["gene"], r["chromosome"], r["allele"])
                for r in genotype_truth_rows(s.genotype))
            for s in cohort.subjects
        }
        assert len(set(signatures.values())) == SUBJECTS, (
            "subjects share a genotype; a per-subject sidecar could not be "
            "told apart from one table repeated")


class TestEveryDrawnAlleleIsInThatSubjectsSidecar:
    def test_no_record_names_an_allele_its_subject_does_not_carry(self, cohort):
        carried = _carried_by_subject(cohort)
        missing = []
        for rec in cohort.records:
            sid = rec["subject_id"]
            for segment in ("v", "d", "j"):
                allele = rec[f"truth_{segment}_call"]
                if allele and allele not in carried[sid]:
                    missing.append((sid, segment, allele))
        assert not missing, (
            f"{len(missing)} record/allele pairs name an allele absent from "
            f"the subject's genotype truth: {sorted(set(missing))[:10]}")

    def test_a_subjects_alleles_are_not_satisfied_by_another_subjects_sidecar(self, cohort):
        """The join is per subject. If completeness only held when pooling
        every subject's alleles, the sidecar would not be per-subject truth.
        """
        carried = _carried_by_subject(cohort)
        pooled = set().union(*carried.values())
        assert any(carried[sid] != pooled for sid in carried), (
            "every subject carries the pooled allele set, so this test cannot "
            "distinguish per-subject truth from pooled truth")

    def test_a_deleted_gene_is_never_drawn_from(self, cohort):
        """The other direction: a gene the genotype deleted must not appear in
        any of that subject's records."""
        for subject in cohort.subjects:
            rows = genotype_truth_rows(subject.genotype)
            fully_deleted = {
                r["gene"] for r in rows if r["gene_state"] == "deleted"}
            if not fully_deleted:
                continue
            drawn_genes = set()
            for rec in subject.result:
                for segment in ("v", "d", "j"):
                    allele = rec[f"truth_{segment}_call"]
                    if allele:
                        drawn_genes.add(allele.split("*")[0])
            assert not (drawn_genes & fully_deleted), (
                f"{subject.subject_id} drew from deleted gene(s): "
                f"{sorted(drawn_genes & fully_deleted)}")


class TestNovelAllelesAreMarkedInTheSidecar:
    def test_a_carried_novel_allele_is_flagged_and_reaches_the_records(self):
        """Pins both halves at once: the sidecar marks the novel allele, and
        the records really did draw it, so the flag describes live truth.

        The novel gene is made the subject's *only* V gene, because every V
        gene in this cartridge carries the same usage weight -- drawing one
        specific gene is a ~0.3% event per record, so a test that merely hoped
        for the draw would need thousands of records and still flake. Deleting
        the rest makes the draw a certainty rather than a probability.
        """
        cfg = _cfg()
        genes = list(cfg.v_alleles)
        target = genes[0]
        novel = f"{target}*i01"
        genotype = (Genotype.from_dataconfig(cfg)
                    .complete_from_reference()
                    .add_novel_allele(novel, base=cfg.v_alleles[target][0].name,
                                      mutations=[(10, "A")], allow_nonfunctional=True)
                    .homozygous(target, novel))
        for other in genes[1:]:
            genotype = genotype.delete_gene(other, haplotype="both")
        genotype = genotype.with_subject("N0")

        result = (ga.Experiment.on(cfg).with_genotype(genotype).recombine()
                  .run_records(n=30, seed=3, expose_provenance=True))

        drawn = {r["truth_v_call"] for r in result}
        assert drawn == {novel}, (
            f"premise: the novel allele must be the only V drawn; got {sorted(drawn)[:5]}")

        rows = genotype_truth_rows(genotype)
        flagged = {r["allele"] for r in rows if r["novel"]}
        assert flagged == {novel}, flagged
        carried = {r["allele"] for r in rows if r["allele"] != DELETED_ALLELE}
        assert drawn <= carried, "the drawn allele is absent from the sidecar"

    def test_the_other_v_genes_are_present_as_deleted_not_absent(self):
        """The same genotype's 51 deleted genes must still be in the table.
        This is the completeness property that a long form could silently
        drop, checked on a genotype where deletions are the majority."""
        cfg = _cfg()
        genes = list(cfg.v_alleles)
        genotype = Genotype.from_dataconfig(cfg).complete_from_reference()
        for other in genes[1:]:
            genotype = genotype.delete_gene(other, haplotype="both")
        genotype = genotype.with_subject("N0")

        rows = genotype_truth_rows(genotype)
        deleted_genes = {r["gene"] for r in rows if r["allele"] == DELETED_ALLELE}

        assert deleted_genes == set(genes[1:]), (
            f"{len(set(genes[1:]) - deleted_genes)} deleted V genes vanished "
            f"from the truth table instead of being reported as deleted")
