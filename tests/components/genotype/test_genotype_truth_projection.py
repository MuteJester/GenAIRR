"""The long, per-carried-allele projection of a genotype's truth.

`Genotype.to_table()` is per gene, with the alleles as list-valued columns
(`haplotype_0`, `haplotype_1`). That shape is wrong for a TSV a benchmark
joins against: a list inside a flat cell makes every consumer write its own
splitter. These tests pin the long form -- one row per (gene, chromosome,
carried allele) -- and, above all, pin the two things a long form can lose:

* a **deleted** gene still gets rows. "Absent from the table" and "present
  and deleted" are different truths, and a genotype caller is scored on
  exactly that difference.
* the `locus` column agrees with the locus the engine puts on the records,
  because the sidecar is useless if the join key disagrees with the join.
"""
import pytest

import GenAIRR.data as gdata
from GenAIRR._genotype_truth import DELETED_ALLELE, genotype_truth_rows
from GenAIRR.genotype import Genotype

pytestmark = pytest.mark.behavior


def _cfg():
    return gdata.HUMAN_IGH_OGRDB


def _gene_with_at_least(n, cfg=None):
    """A V gene carrying at least `n` distinct alleles, or skip.

    Guards the premise: a heterozygous test built on a single-allele gene
    would construct a homozygote and pass for the wrong reason.
    """
    cfg = cfg or _cfg()
    for gene, alleles in cfg.v_alleles.items():
        if len({a.name for a in alleles}) >= n:
            return gene, sorted({a.name for a in alleles})
    pytest.skip(f"no V gene with >= {n} distinct alleles in this cartridge")


def _rows_for(genotype, gene):
    return [r for r in genotype_truth_rows(genotype) if r["gene"] == gene]


class TestOneRowPerCarriedAllele:
    def test_a_homozygous_gene_gives_the_same_allele_on_both_chromosomes(self):
        gene, alleles = _gene_with_at_least(1)
        g = Genotype.from_dataconfig(_cfg()).homozygous(gene, alleles[0]).with_subject("S1")

        rows = _rows_for(g, gene)

        assert [r["chromosome"] for r in rows] == [0, 1]
        assert {r["allele"] for r in rows} == {alleles[0]}
        assert {r["gene_state"] for r in rows} == {"homozygous"}
        assert {r["subject_id"] for r in rows} == {"S1"}

    def test_a_heterozygous_gene_gives_a_different_allele_per_chromosome(self):
        gene, alleles = _gene_with_at_least(2)
        a0, a1 = alleles[0], alleles[1]
        assert a0 != a1, "premise: the two alleles must differ"
        g = Genotype.from_dataconfig(_cfg()).heterozygous(gene, a0, a1).with_subject("S1")

        rows = _rows_for(g, gene)

        assert [(r["chromosome"], r["allele"]) for r in rows] == [(0, a0), (1, a1)]
        assert {r["gene_state"] for r in rows} == {"heterozygous"}

    def test_every_row_names_its_segment(self):
        gene, alleles = _gene_with_at_least(1)
        g = Genotype.from_dataconfig(_cfg()).homozygous(gene, alleles[0]).with_subject("S1")

        assert {r["segment"] for r in _rows_for(g, gene)} == {"V"}


class TestADeletedGeneStillGetsRows:
    """The one place a long table can lose what a wide one keeps."""

    def test_deleted_on_one_chromosome_keeps_one_allele_row_and_one_sentinel(self):
        gene, alleles = _gene_with_at_least(1)
        g = (Genotype.from_dataconfig(_cfg())
             .homozygous(gene, alleles[0])
             .delete_gene(gene, haplotype=1)
             .with_subject("S1"))

        rows = _rows_for(g, gene)

        assert [(r["chromosome"], r["allele"]) for r in rows] == [
            (0, alleles[0]), (1, DELETED_ALLELE)]
        assert {r["gene_state"] for r in rows} == {"hemizygous"}

    def test_deleted_on_both_chromosomes_gives_two_sentinel_rows(self):
        gene, alleles = _gene_with_at_least(1)
        g = (Genotype.from_dataconfig(_cfg())
             .homozygous(gene, alleles[0])
             .delete_gene(gene, haplotype="both")
             .with_subject("S1"))

        rows = _rows_for(g, gene)

        assert [(r["chromosome"], r["allele"]) for r in rows] == [
            (0, DELETED_ALLELE), (1, DELETED_ALLELE)]
        assert {r["gene_state"] for r in rows} == {"deleted"}

    def test_the_sentinel_is_a_literal_token_and_not_an_empty_cell(self):
        """An empty TSV cell cannot be told apart from a missing value, which
        is the one confusion a truth file must not create."""
        assert DELETED_ALLELE == "(none)"
        assert DELETED_ALLELE.strip() != ""


class TestNovelAlleles:
    """`allow_nonfunctional=True` throughout: a single substitution at a fixed
    offset lands in a stop codon for some genes, and the projection does not
    read functionality at all. Hunting for a position that stays functional
    would couple these tests to the cartridge's codon layout for no gain."""

    def test_a_placed_novel_allele_is_marked_novel(self):
        cfg = _cfg()
        gene, alleles = _gene_with_at_least(1, cfg)
        base = alleles[0]
        g = (Genotype.from_dataconfig(cfg)
             .add_novel_allele(f"{gene}*i01", base=base, mutations=[(10, "A")],
                               allow_nonfunctional=True)
             .homozygous(gene, f"{gene}*i01")
             .with_subject("S1"))

        rows = _rows_for(g, gene)

        assert rows, "premise: the placed novel allele must produce rows"
        assert all(r["novel"] is True for r in rows), rows
        assert {r["allele"] for r in rows} == {f"{gene}*i01"}

    def test_a_carried_reference_allele_is_not_marked_novel(self):
        gene, alleles = _gene_with_at_least(1)
        g = Genotype.from_dataconfig(_cfg()).homozygous(gene, alleles[0]).with_subject("S1")

        assert all(r["novel"] is False for r in _rows_for(g, gene))

    def test_an_unplaced_novel_allele_has_no_row(self):
        """Defined but never assigned to a chromosome. It is not carried, so
        it is not truth about this subject's genome."""
        cfg = _cfg()
        gene, alleles = _gene_with_at_least(1, cfg)
        g = (Genotype.from_dataconfig(cfg)
             .add_novel_allele(f"{gene}*i01", base=alleles[0], mutations=[(10, "A")],
                               allow_nonfunctional=True)
             .homozygous(gene, alleles[0])
             .with_subject("S1"))

        # Premise: the allele really was defined. Without this the test
        # passes just as well against an `add_novel_allele` that did nothing.
        assert g.has_novel(), "premise: a novel allele must be defined"
        assert f"{gene}*i01" in g.novel_allele_names()

        rows = genotype_truth_rows(g)

        assert rows, "premise: the gene must still produce rows"
        assert all(r["allele"] != f"{gene}*i01" for r in rows), (
            "an unplaced novel allele must not appear as carried")


class TestTheLocusColumn:
    """What the locus column is, and what it is not.

    The class was once named for the records and never looked at one, which
    is how it missed the case below: on a cartridge whose V pool mixes loci,
    a record's locus and its own allele's sidecar locus genuinely differ.
    The record carries one locus for the whole rearrangement, taken from the
    V call; a row carries the locus of its gene. The join key is
    `(subject_id, allele)`, so the difference is a readability property of
    the column rather than a broken join -- but it has to be pinned, because
    the obvious wrong move is to filter the sidecar by the record's locus.
    """

    def test_the_locus_is_the_prefix_the_engine_derives(self):
        gene, alleles = _gene_with_at_least(1)
        g = Genotype.from_dataconfig(_cfg()).homozygous(gene, alleles[0]).with_subject("S1")

        assert {r["locus"] for r in _rows_for(g, gene)} == {"IGH"}

    def test_a_name_with_no_airr_locus_prefix_yields_an_empty_locus(self):
        """`derive_locus` in the engine accepts only the seven AIRR loci and
        returns empty otherwise; the projection must not invent one."""
        from GenAIRR._genotype_truth import locus_of

        assert locus_of("IGHVF1-G2") == "IGH"
        assert locus_of("TRBV20-1") == "TRB"
        assert locus_of("igkv1-5") == "IGK"
        assert locus_of("SYNTHETIC1") == ""
        assert locus_of("XX") == ""

    def test_on_a_single_locus_cartridge_it_agrees_with_every_record(self):
        """The property the class name used to assert, now against records."""
        import GenAIRR as ga

        g = (Genotype.from_dataconfig(_cfg()).complete_from_reference()
             .with_subject("S1"))
        result = (ga.Experiment.on(_cfg()).with_genotype(g).recombine()
                  .run_records(n=20, seed=3))
        rows = genotype_truth_rows(g)
        assert rows, "premise: the genotype projects to rows"
        record_loci = {rec["locus"] for rec in result}
        assert record_loci == {"IGH"}, record_loci
        assert {r["locus"] for r in rows} == record_loci

    def test_on_a_mixed_locus_cartridge_a_rows_locus_is_its_genes_not_the_reads(self):
        """`HUMAN_TCRD_IMGT` carries both TRAV and TRDV genes, so a read with
        a TRAV V call is labelled TRA while the TRDJ allele it was generated
        from is a TRD row. Pinned so nobody 'fixes' the column to match the
        record and loses the D/J truth."""
        import GenAIRR as ga
        import GenAIRR.data as gdata

        cfg = gdata.HUMAN_TCRD_IMGT
        v_prefixes = {gene[:4] for gene in cfg.v_alleles}
        assert v_prefixes >= {"TRAV", "TRDV"}, (
            f"premise: this cartridge's V pool mixes loci, got {v_prefixes}")

        g = (Genotype.from_dataconfig(cfg).complete_from_reference()
             .with_subject("T0"))
        result = (ga.Experiment.on(cfg).with_genotype(g).recombine()
                  .run_records(n=40, seed=5, expose_provenance=True))
        locus_by_allele = {r["allele"]: r["locus"] for r in genotype_truth_rows(g)}

        disagreeing = [rec for rec in result
                       if locus_by_allele.get(rec["truth_j_call"]) != rec["locus"]]
        assert disagreeing, (
            "premise: this cartridge is supposed to produce reads whose J "
            "allele sits in a different locus than the read's own label")
        # The join is on the allele, and it holds for every record regardless.
        for rec in result:
            for field in ("truth_v_call", "truth_d_call", "truth_j_call"):
                allele = rec.get(field)
                if allele:
                    assert allele in locus_by_allele, allele
