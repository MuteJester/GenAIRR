"""The long projection of a genotype's carried-allele truth.

``Genotype.to_table()`` is one row per gene, with the carried alleles as
list-valued columns (``haplotype_0``, ``haplotype_1``). That is the right
shape for a Python caller and the wrong shape for a TSV a benchmark joins
against: a list inside a flat cell is a nested format inside a flat one, and
every consumer ends up writing its own splitter.

This module projects that table into the long form -- one row per (gene,
chromosome, carried allele) -- which is what a join wants. ``to_table()`` is
left exactly as it is; it is a public, tested surface with its own callers,
and this is a projection of it rather than a replacement.

Two properties matter more than the shape:

**A deleted gene still gets rows.** This is the one point where a long table
can lose information the wide one keeps. "Gene absent from the table" and
"gene present and deleted on both chromosomes" are different truths, and a
genotype caller is scored on exactly that difference. A gene deleted on one
chromosome gets one allele row and one ``DELETED_ALLELE`` row; deleted on
both gets two. ``gene_state`` carries the zygosity verbatim, so the wide form
stays recoverable from the long one.

**The locus names the gene's locus, which is not always the record's.**
``locus_of`` is a port of the engine's naming rule (``locus_prefix`` in
``engine_rs/src/airr_record/locus.rs``): the uppercased three-character
prefix, accepted only when it is one of the seven AIRR loci, and empty
otherwise. A synthetic reference whose allele names are not locus-prefixed
gets an empty locus here, as it does on the record.

The two values are *not* interchangeable, and the difference is real rather
than theoretical. A record carries one locus for the whole rearrangement,
taken from its V call first (``derive_locus``). A sidecar row carries the
locus of the gene in that row. On a cartridge whose V pool mixes loci --
``HUMAN_TCRD_IMGT`` has both ``TRAV`` and ``TRDV`` genes, as do the two
mouse TCRD cartridges -- a read with a ``TRAV`` V call and a ``TRDJ`` J call
is labelled ``TRA``, while the sidecar row for the very ``TRDJ`` allele it
was generated from says ``TRD``. Measured on that cartridge: 134 of 200
records. So **join on ``(subject_id, allele)`` and do not filter the sidecar
by the record's locus first** -- that filter silently drops the D and J
truth for most reads. The locus column is there to make a row readable on
its own, not to be matched against the record's.
"""
from typing import Dict, List, Tuple

from .genotype._common import _SEGMENTS, _alleles_by_gene

# The seven AIRR loci, as the engine spells them. TCR loci are TRA/TRB/TRG/TRD
# here, not the IMGT-style TCRA/TCRB -- the engine emits the former.
AIRR_LOCI = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")

#: What a deleted chromosome's allele cell holds. A literal token, not an
#: empty cell: an empty TSV cell cannot be told apart from a missing value,
#: which is the one confusion a ground-truth file must not create.
DELETED_ALLELE = "(none)"

#: The long form's columns, in order.
COLUMNS = (
    "subject_id",
    "locus",
    "segment",
    "gene",
    "allele",
    "chromosome",
    "novel",
    "gene_state",
)


def locus_of(name: str) -> str:
    """The AIRR locus a gene or allele name belongs to, or ``""``.

    A port of the engine's ``locus_prefix``, kept deliberately identical so
    the sidecar's ``locus`` column equals the ``locus`` on the records it is
    joined to. Returning ``""`` rather than guessing is the point: a name
    outside the seven loci names no locus, and inventing one would put a
    value in the truth file that the records disagree with.
    """
    if len(name) < 3:
        return ""
    prefix = name[:3].upper()
    return prefix if prefix in AIRR_LOCI else ""


def genotype_truth_rows(genotype) -> List[Dict]:
    """One row per (gene, chromosome, carried allele) for one subject.

    Deleted chromosomes contribute a ``DELETED_ALLELE`` row rather than no
    row. A novel allele that was defined but never placed on a chromosome
    contributes nothing: it is not carried, so it is not truth about this
    subject's genome.
    """
    rows: List[Dict] = []
    for wide in genotype.to_table():
        novel = set(wide["novel"])
        gene = wide["gene"]
        locus = locus_of(gene)
        for chromosome in (0, 1):
            carried = wide[f"haplotype_{chromosome}"]
            # A chromosome with no allele is a deletion, and it is reported.
            alleles = carried if carried else [DELETED_ALLELE]
            for allele in alleles:
                rows.append({
                    "subject_id": wide["subject_id"],
                    "locus": locus,
                    "segment": wide["segment"],
                    "gene": gene,
                    "allele": allele,
                    "chromosome": chromosome,
                    "novel": allele in novel,
                    "gene_state": wide["zygosity"],
                })
    return rows


def unspecified_genes(genotype) -> List[Tuple[str, str]]:
    """The reference genes this genotype does not pin, as ``(segment, gene)``.

    Empty for a strict or completed genotype. Non-empty only for a
    ``Genotype.permissive`` that was never completed, where the engine samples
    the unspecified genes over the whole reference: those alleles reach the
    records but no row of the projection, which is the one way this table can
    be a silent subset of the truth.

    Keyed on what is specified rather than on ``is_permissive``, because
    ``permissive().complete_from_reference()`` stays permissive while its
    table becomes complete.
    """
    cfg = getattr(genotype, "_cfg", None)
    if cfg is None:
        return []
    missing: List[Tuple[str, str]] = []
    for segment in _SEGMENTS:
        for gene in _alleles_by_gene(cfg, segment):
            if not genotype.is_specified(segment, gene):
                missing.append((segment, gene))
    return missing


def write_genotype_tsv(path: str, genotypes) -> None:
    """Write the long truth table for ``genotypes`` to ``path`` as TSV.

    Refuses, before opening the file, when there is nothing to write or when
    a genotype cannot be joined. Both refusals exist for the same reason: a
    truth file that is empty, or whose ``subject_id`` is blank, sits beside a
    records file looking like an export that worked. Raising leaves the
    caller with an error instead of a plausible artifact.
    """
    import csv

    supplied = list(genotypes or ())
    carried = [g for g in supplied if g is not None]
    if not carried:
        raise ValueError(
            "this result has no genotype, so there is no genotype truth to "
            "write; run with `.with_genotype(...)` or `run_cohort(...)` to "
            "produce a result that carries one")

    # Dropping the Nones would write a table that is short by whole subjects
    # and says nothing about it -- a partial file that looks complete.
    missing_genotype = [i for i, g in enumerate(supplied) if g is None]
    if missing_genotype:
        raise ValueError(
            f"subject(s) at position(s) {missing_genotype} have no genotype, so "
            f"the sidecar would be missing those subjects entirely while "
            f"looking complete; give every subject a genotype or export the "
            f"ones that have them separately")

    unjoinable = [i for i, g in enumerate(carried) if not getattr(g, "subject_id", None)]
    if unjoinable:
        raise ValueError(
            f"genotype(s) at position(s) {unjoinable} carry no subject_id, so "
            f"the sidecar could not be joined to the records; set one with "
            f"`.with_subject(...)`")

    seen: Dict[str, int] = {}
    for index, genotype in enumerate(carried):
        sid = genotype.subject_id
        if sid in seen:
            raise ValueError(
                f"subject_id {sid!r} appears on genotypes at positions "
                f"{seen[sid]} and {index}, so their rows would merge under one "
                f"join key and neither subject's truth would be recoverable; "
                f"give each subject a distinct id with `.with_subject(...)`")
        seen[sid] = index

    for genotype in carried:
        missing = unspecified_genes(genotype)
        if missing:
            shown = ", ".join(f"{seg}:{gene}" for seg, gene in missing[:3])
            raise ValueError(
                f"genotype {getattr(genotype, 'subject_id', None)!r} leaves "
                f"{len(missing)} reference gene(s) unspecified (e.g. {shown}), "
                f"so this table is not the set of alleles the records were "
                f"generated from: the engine samples those genes over the "
                f"whole reference, unphased, and they would be missing from "
                f"the file. Call `complete_from_reference()` to make the "
                f"carried set explicit, or export the partial picture with "
                f"`Genotype.to_tsv(...)`, which records that it is permissive")

    rows = [row for g in carried for row in genotype_truth_rows(g)]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        # LF, matching `_write_delimited`'s deliberate `lineterminator="\n"`.
        # The excel dialect's default CRLF would leave a stray \r on
        # `gene_state`, the last column, in any shell pipeline that joins the
        # sidecar to the records file.
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for row in rows:
            writer.writerow([row[c] for c in COLUMNS])
