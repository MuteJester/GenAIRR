"""The coordinate books must balance against the reference, not just against themselves.

The assembly truth test next door compares each recorded germline *base* to the
reference base at the coordinate the record claims. That establishes allele
identity, and deliberately claims nothing about the coordinates: production
lays the germline byte and its coordinate down from the same two numbers, so a
trim that removes one base too many moves the byte and the coordinate together
and a base comparison still agrees.

This module never compares a base. It checks that four independently produced
numbers agree -- the trim amounts the trimming pass sampled, the end-loss
lengths the assay recorded, the coordinates the record builder exported from
the nucleotide pool, and the allele's own length in the reference table. A
trimming pass that removes four bases while recording three satisfies the base
comparison and breaks this one.
"""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import coordinates

pytestmark = pytest.mark.behavior

#: Loci spanning both chain shapes, so a VJ locus (no D) is covered too.
LOCI = ("human_igh", "human_igk", "human_tcrb")


def _observed(exp, *, n: int, seed: int):
    """``(record, allele_lengths)`` per outcome, lengths read from the reference."""
    refdata = exp.refdata
    getters = {"v": refdata.v_allele, "d": refdata.d_allele, "j": refdata.j_allele}
    for outcome in exp.run(n=n, seed=seed):
        simulation = outcome.final_simulation()
        identifiers = {"v": simulation.v_allele_id(), "d": simulation.d_allele_id(),
                       "j": simulation.j_allele_id()}
        lengths = {}
        for segment, identifier in identifiers.items():
            if identifier is None:
                continue
            allele = getters[segment](identifier)
            sequence = allele.seq() if callable(allele.seq) else allele.seq
            lengths[segment] = len(sequence)
        yield outcome.airr_record(refdata), lengths


def _trim_fields_present(record):
    """Every trim field, named literally rather than built from a segment letter.

    The oracle reads them through f-strings, which the catalog's claimed-field
    check cannot see -- and a claim it cannot see is a claim nobody verifies.
    Naming them here also makes a record that stopped exporting one fail.
    """
    for field in ("v_trim_5", "v_trim_3", "d_trim_5", "d_trim_3",
                  "j_trim_5", "j_trim_3"):
        value = record.get(field)
        assert value is None or isinstance(value, int), f"{field} is not an integer"


def _coordinate_fields_present(record):
    """Every read and alignment coordinate, named literally. See above."""
    for field in ("v_sequence_start", "v_sequence_end", "d_sequence_start",
                  "d_sequence_end", "j_sequence_start", "j_sequence_end",
                  "v_alignment_start", "v_alignment_end", "d_alignment_start",
                  "d_alignment_end", "j_alignment_start", "j_alignment_end"):
        value = record.get(field)
        assert value is None or isinstance(value, int), f"{field} is not an integer"


@verifies("recombination.trimming",
          claim="every segment's recorded germline span starts at its own trim_5 "
                "and ends at the reference allele's length minus its trim_3, so "
                "the trim each coding end reports is the germline actually "
                "missing -- a trim that removed a different number of bases than "
                "it recorded breaks this even though the bases it did keep are "
                "the allele's own",
          fields=("v_trim_5", "v_trim_3", "d_trim_5", "d_trim_3",
                  "j_trim_5", "j_trim_3"))
@pytest.mark.parametrize("locus", LOCI)
def test_recorded_trims_balance_against_the_reference_allele_length(locus):
    exp = experiment_builder.recombined(locus=locus)
    segments_seen: set[str] = set()
    for record, lengths in _observed(exp, n=40, seed=19):
        found = coordinates.disagreements(record, lengths)
        assert found == (), f"{locus}: {len(found)} coordinate disagreement(s); first: {found[0]}"
        _trim_fields_present(record)
        segments_seen.update(lengths)
    expected = {"v", "d", "j"} if locus != "human_igk" else {"v", "j"}
    assert segments_seen == expected, (
        f"{locus}: checked {sorted(segments_seen)}, expected {sorted(expected)}")


@verifies("projection.segment_coordinates",
          claim="each segment's span in the read equals its span in the germline, "
                "and its alignment coordinates equal its read coordinates while "
                "no indel has displaced them -- so the coordinates a downstream "
                "aligner is scored on are the same span the germline accounts for",
          fields=("v_sequence_start", "v_sequence_end", "d_sequence_start",
                  "d_sequence_end", "j_sequence_start", "j_sequence_end",
                  "v_alignment_start", "v_alignment_end", "d_alignment_start",
                  "d_alignment_end", "j_alignment_start", "j_alignment_end"))
@pytest.mark.parametrize("locus", LOCI)
def test_read_spans_and_alignment_coordinates_agree_with_the_germline_span(locus):
    exp = experiment_builder.recombined(locus=locus)
    checked = 0
    for record, lengths in _observed(exp, n=40, seed=23):
        found = coordinates.disagreements(record, lengths)
        assert found == (), f"{locus}: {len(found)} coordinate disagreement(s); first: {found[0]}"
        _coordinate_fields_present(record)
        for segment in lengths:
            for suffix in ("sequence_start", "sequence_end", "alignment_start", "alignment_end"):
                assert record[f"{segment}_{suffix}"] is not None, (
                    f"{locus}: {segment}_{suffix} is absent on a segment that was assembled")
        checked += 1
    assert checked == 40, f"{locus}: only {checked} records were checked"


@verifies("assay.end_loss",
          claim="end loss is accounted for in the coordinate books: the V germline "
                "span begins end_loss_5_length further in, and the J germline span "
                "ends end_loss_3_length earlier. The identity alone does not "
                "separate read truncation from exonuclease loss -- charging the "
                "6 lost bases to v_trim_5 instead satisfies it -- so the "
                "separation is asserted here directly, by the two trim fields "
                "that no biology may write staying zero",
          fields=("end_loss_5_length", "end_loss_3_length"))
def test_end_loss_is_accounted_for_in_the_germline_coordinates():
    """Checked together with trimming on purpose, because both shorten the
    germline span and the arithmetic cannot tell them apart on its own.

    ``germline_start == trim_5 + end_loss_5_length`` is satisfied by any split
    of the loss between the two fields. What makes the split checkable is that
    one of the two is structurally zero: RAG cleaves only the RSS-adjacent
    coding end, so a non-zero v_trim_5 would mean read truncation had been
    recorded as germline that was never there.
    """
    exp = (experiment_builder.recombined()
           .end_loss_5prime(length=6)
           .end_loss_3prime(length=5))
    losses = set()
    for record, lengths in _observed(exp, n=40, seed=29):
        found = coordinates.disagreements(record, lengths)
        assert found == (), f"{len(found)} coordinate disagreement(s); first: {found[0]}"
        assert record["v_trim_5"] == 0, (
            "the 5' read truncation was recorded as V exonuclease loss")
        assert record["j_trim_3"] == 0, (
            "the 3' read truncation was recorded as J exonuclease loss")
        losses.add((record["end_loss_5_length"], record["end_loss_3_length"]))
    assert losses == {(6, 5)}, f"the fixture did not apply the configured end loss: {losses}"


@pytest.mark.parametrize("stack, build", [
    ("indels", lambda e: e.polymerase_indels(count=3)),
    ("reverse-complemented", lambda e: e.random_strand_orientation(prob=1.0)),
    ("indels and end loss and rc",
     lambda e: (e.polymerase_indels(count=2)
                 .end_loss_5prime(length=8).end_loss_3prime(length=7)
                 .random_strand_orientation(prob=0.5))),
])
def test_the_books_still_balance_once_edits_displace_the_coordinates(stack, build):
    """Indels are the only thing allowed to perturb the identities, and only by
    their own recorded count. Reverse complement must not perturb them at all:
    the coordinate fields are reported in germline orientation."""
    exp = build(experiment_builder.recombined())
    perturbed = 0
    for record, lengths in _observed(exp, n=60, seed=31):
        found = coordinates.disagreements(record, lengths)
        assert found == (), f"{stack}: {len(found)} disagreement(s); first: {found[0]}"
        if (record.get("n_indels") or 0) > 0:
            perturbed += 1
    if stack != "reverse-complemented":
        assert perturbed > 0, f"{stack}: no record carried an indel, so the bound was never exercised"
