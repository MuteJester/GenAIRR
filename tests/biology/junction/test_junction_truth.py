"""The junction window, against the anchors the reference declares.

Every CDR3 benchmark reads these five fields. Clustering thresholds,
repertoire overlap and every length distribution are computed from them, and
`junction_start` / `junction_end` are what a tool is told the window is.
Nothing registered checked any of them before this module: both mechanisms
were `unverified` with zero Python tests.
"""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support import alleles
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import junction as junction_oracle
from tests.support.oracles import reconstruction

pytestmark = pytest.mark.behavior

#: TCRB is here because the pipeline-space sweep ran the junction oracle on it
#: and nothing else ever had: both junction mechanisms were verified on the two
#: BCR loci only, while TCRB is the locus whose reference carries six alleles
#: with an unresolvable anchor and so the one where the null-junction contract
#: actually fires.
LOCI = ("human_igh", "human_igk", "human_tcrb")


def _check(exp, *, n, seed):
    """Every disagreement the junction oracle reports, and how much it saw.

    The anchor coordinate and the reference bases are read here, from the
    refdata, and handed to the oracle -- which never imports GenAIRR and so
    cannot look them up itself. That is what keeps the expected value on the
    other side of the thing being checked.
    """
    found, checked, records = [], 0, []
    for outcome in exp.run(n=n, seed=seed):
        record = outcome.airr_record(exp.refdata)
        records.append(record)
        simulation = outcome.final_simulation()
        v_allele = alleles.by_name(exp.refdata, "V", record["v_call"])
        j_allele = alleles.by_name(exp.refdata, "J", record["j_call"])
        found.extend(junction_oracle.disagreements(
            record, simulation,
            reconstruction.segments_by_position(simulation),
            v_anchor=v_allele.anchor, j_anchor=j_allele.anchor,
            v_sequence=alleles.sequence(v_allele),
            j_sequence=alleles.sequence(j_allele),
            # Whether a bundled allele's anchor codon is canonical is a claim
            # about the reference cartridge, and six TCRB alleles are not --
            # correctly, with the engine refusing to emit a junction for them.
            # test_anchor_reference_quality.py owns that claim; leaving it on
            # here would make every TCRB record carrying one of those alleles
            # a failure of a per-record test.
            reference_residue=False))
        checked += 1
    return tuple(found), checked, records


def _junction_bound_fields_present(record):
    """Both bounds, named literally.

    The claimed fields have to be read by the body that claims them, and the
    oracle reads them one module away -- so they are named here, as the
    coordinate truth tests do for the same reason. The type check is the
    weakest part of this file; the oracle comparison above it is what
    establishes the values.
    """
    if _has_no_junction(record):
        return False
    for field in ("junction_start", "junction_end"):
        assert isinstance(record[field], int), f"{field} is not an integer"
    assert 0 <= record["junction_start"] < record["junction_end"] <= len(record["sequence"]), (
        f"bounds {record['junction_start']}..{record['junction_end']} fall "
        f"outside a {len(record['sequence'])}-base molecule")
    return True


def _reference_resolves_both_anchors(v_allele, j_allele):
    """Whether the reference resolves this record's two anchors.

    The engine's own rule: an anchor is resolved when the codon at it is the
    conserved residue for that segment. Six bundled TCRB alleles fail it, and
    for those the engine computes no junction at all -- which is right, and
    means there is no window for a test to check.
    """
    for allele, canonical in ((v_allele, junction_oracle.CYSTEINE),
                              (j_allele, junction_oracle.TRYPTOPHAN_OR_PHENYLALANINE)):
        codon = junction_oracle.reference_anchor_codon(
            alleles.sequence(allele), allele.anchor)
        if codon is None or codon not in canonical:
            return False
    return True


def _has_no_junction(record):
    """Whether this record legitimately reports no junction, asserted as a group.

    A record whose anchor the reference does not resolve -- six of the bundled
    TCRB alleles -- or whose anchor codon trimming removed has no window to
    report. The presence helpers must not read that as a missing field, and
    must not let it pass unchecked either: the whole group goes null together
    or the record is lying about something.
    """
    if record["junction_start"] is not None:
        return False
    assert record["junction_end"] is None, (
        "junction_start is null while junction_end is not, so half a window "
        "is being reported")
    assert record["junction"] == "" and record["junction_aa"] == "", (
        f"no junction bounds but junction={record['junction']!r} "
        f"junction_aa={record['junction_aa']!r}")
    assert record["junction_length"] is None, (
        f"no junction bounds but junction_length={record['junction_length']!r}")
    assert record["vj_in_frame"] is None, (
        "no junction bounds but vj_in_frame is "
        f"{record['vj_in_frame']!r}, which a consumer reads as out-of-frame")
    assert record["productivity_reasons"], (
        "a record with no junction says nothing about why")
    return True


def _junction_string_fields_present(record):
    """The three reported fields, named literally. See above."""
    if _has_no_junction(record):
        return False
    assert isinstance(record["junction"], str) and record["junction"], (
        "junction is empty")
    assert record["junction_length"] == len(record["junction"]), (
        f"junction_length {record['junction_length']} for a "
        f"{len(record['junction'])}-base junction")
    assert len(record["junction_aa"]) == record["junction_length"] // 3, (
        f"junction_aa has {len(record['junction_aa'])} residues for "
        f"{record['junction_length']} nucleotides")
    return True


@verifies("junction.coding_end_processing",
          claim="junction_start is the read index of the first surviving base "
                "of the V anchor codon and junction_end one past the last "
                "surviving base of the J anchor codon, with the anchor "
                "coordinate and the conserved residue read from the reference "
                "rather than from the run. What this does not establish is "
                "that the bound names the CDR3 anchor specifically: 197 of 198 "
                "IGH V alleles carry a second in-frame cysteine upstream, and "
                "a bound displaced onto it satisfies the residue check",
          fields=("junction_start", "junction_end"),
          basis="reference")
@pytest.mark.parametrize("locus", LOCI)
def test_the_junction_bounds_are_the_surviving_anchor_codons(locus):
    """Trimming can eat part of a conserved codon, which is why the bound is
    the first *surviving* base rather than the codon start: the oracle walks
    the provenance for whichever of the three anchor positions is still there.

    Both loci run because the J anchor residue differs -- tryptophan on the
    heavy chain, phenylalanine on the light -- and a check that only ever saw
    TGG would pass against a residue set hardcoded to it.
    """
    exp = experiment_builder.recombined(locus=locus)
    found, checked, records = _check(exp, n=40, seed=5)
    assert found == (), f"{locus}: {found[:3]}"
    assert checked == 40, f"{locus}: only {checked} records were checked"
    with_window = sum(_junction_bound_fields_present(record) for record in records)
    assert with_window >= 35, (
        f"{locus}: only {with_window} of {len(records)} records reported a "
        "junction window; on TCRB a few draw an allele whose anchor the "
        "reference does not resolve, many more than that is a defect")


@verifies("junction.assembly",
          claim="junction and junction_length are the window the REFERENCE "
                "anchors fix: the expected bounds are walked out of the "
                "recorded provenance to the surviving bases of the V and J "
                "anchor codons the reference declares, and the exported "
                "junction_start / junction_end are not consulted. A window "
                "displaced as a whole satisfies the slice identity beside "
                "this one and fails here",
          fields=("junction", "junction_length"),
          basis="reference")
@pytest.mark.parametrize("locus", LOCI)
def test_the_junction_is_the_window_the_reference_anchors_fix(locus):
    """Why the slice identity next door is not enough.

    That one compares ``junction`` to the molecule sliced at the recorded
    bounds. Both bounds come from the record, so a junction and a pair of
    bounds that moved together satisfy it -- the oracle's circularity map
    records the identity as CIRCULAR once the bounds move with it, and
    ``junction_length == end - start`` is arithmetic over those same two
    exported numbers.

    Here the expected window comes from the anchor coordinates the reference
    declares, walked through the recorded provenance to whichever anchor bases
    survived trimming. Nothing on the expected side is the record builder's
    output. Measured over both loci, 40 records each: 40 of 40 agree, and 0 of
    40 would agree with the window shifted one base along -- so a displaced
    junction has somewhere to show up.
    """
    exp = experiment_builder.recombined(locus=locus)
    checked = unresolved = 0
    for index, outcome in enumerate(exp.run(n=40, seed=5)):
        record = outcome.airr_record(exp.refdata)
        simulation = outcome.final_simulation()
        v_allele = alleles.by_name(exp.refdata, "V", record["v_call"])
        j_allele = alleles.by_name(exp.refdata, "J", record["j_call"])
        segment_of = reconstruction.segments_by_position(simulation)
        if not _reference_resolves_both_anchors(v_allele, j_allele):
            # The reference does not resolve one of these anchors, so there is
            # no anchor-fixed window to compute. Gated on the REFERENCE rather
            # than on the record reporting no junction: the record's own
            # null-ness is what test_anchor_reference_quality.py establishes,
            # and using it here would skip exactly the records a wrongly-null
            # junction would hide in.
            unresolved += 1
            continue
        v_codon = junction_oracle.anchor_codon_positions(
            simulation, segment_of, "V", v_allele.anchor)
        j_codon = junction_oracle.anchor_codon_positions(
            simulation, segment_of, "J", j_allele.anchor)
        if not v_codon or not j_codon:
            # No surviving anchor codon means there is no window to compare:
            # either trimming removed it, or -- on TCRB -- the reference never
            # resolved it. The empty-junction group is its own identity,
            # checked by the oracle run in the test below, and the six
            # unresolvable TCRB alleles are pinned by name in
            # test_anchor_reference_quality.py.
            unresolved += 1
            continue
        start, end = v_codon[0], j_codon[-1] + 1
        expected = junction_oracle.subject(record)[start:end]
        assert record["junction"].upper() == expected.upper(), (
            f"{locus}: record {index} reports a junction that is not the "
            f"molecule between the surviving anchor bases at {start}..{end}")
        assert record["junction_length"] == end - start, (
            f"{locus}: record {index} reports junction_length "
            f"{record['junction_length']} for an anchor-fixed window of "
            f"{end - start}")
        checked += 1
    assert checked + unresolved == 40, (
        f"{locus}: {checked} compared and {unresolved} unresolvable of 40 "
        "records, which do not add up -- a record was dropped silently")
    assert checked >= 35, (
        f"{locus}: only {checked} of 40 records had a window to compare. On "
        "TCRB a few draw one of the six alleles whose anchor the reference "
        "does not resolve; many more than that means something else is "
        "removing windows")


@verifies("junction.assembly",
          claim="junction is the molecule sliced at the recorded bounds -- in "
                "germline orientation, which is probed rather than assumed -- "
                "junction_length is the distance between them, and junction_aa "
                "is the translation of junction under an independent codon "
                "table",
          fields=("junction_aa",),
          basis="reference")
@pytest.mark.parametrize("locus", LOCI)
def test_the_junction_string_is_the_molecule_between_its_bounds(locus):
    """The three reported fields and the window they claim to describe.

    Of these, the translation is the one registered here: the codon table is
    outside the engine, so ``junction_aa`` has an independent expected value.
    The slice and the length are comparisons among the record's own fields --
    the oracle's circularity map says so, and the test above is what gives
    those two an external side. They still run here, because the oracle checks
    every identity it has on one fixture.
    """
    exp = experiment_builder.recombined(locus=locus)
    found, checked, records = _check(exp, n=40, seed=5)
    assert found == (), f"{locus}: {found[:3]}"
    assert checked == 40
    with_window = sum(_junction_string_fields_present(record) for record in records)
    assert with_window >= 35, (
        f"{locus}: only {with_window} of {len(records)} records reported a "
        "junction string; see the bound test for why a few may not")


@verifies("junction.assembly",
          claim="the junction is reported in germline orientation: on a "
                "reverse-complemented record it is the slice of the "
                "reverse complement of sequence, not of sequence",
          fields=("junction",),
          basis="reference")
def test_the_junction_is_reported_in_germline_orientation():
    """The rule that makes every junction field interpretable on a flipped
    record, and the one an oracle is most likely to get wrong silently -- the
    slice simply comes out as the wrong string.

    Run at prob=1.0 so every record is flipped, with the flag asserted, so a
    build that stopped flipping would fail here rather than quietly passing
    against unflipped reads.
    """
    exp = experiment_builder.recombined().random_strand_orientation(prob=1.0)
    found, checked, records = _check(exp, n=25, seed=5)
    assert found == (), f"{found[:3]}"
    assert checked == 25
    for record in records:
        assert record["junction"], "a flipped record reported an empty junction"
    flags = [o.airr_record(exp.refdata)["rev_comp"]
             for o in exp.run(n=25, seed=5)]
    assert all(flags), "no record was reverse-complemented; the rule is untested"


@verifies("junction.coding_end_processing",
          claim="a mutation pass does not move the junction bounds: "
                "hypermutation substitutes bases without changing which "
                "germline positions survive, so the window is the same window "
                "the unmutated molecule reports under the same seed",
          fields=("junction_start", "junction_end"),
          basis="cross-run")
def test_hypermutation_does_not_move_the_junction_window():
    """A cross-run differential: the bounds are derived from provenance, and
    a substitution changes a base without changing its germline coordinate,
    so adding SHM to the same seed must leave both bounds alone.

    This is the one check here whose expected value comes from a second run
    rather than from the reference, which is why it can see something the
    oracle cannot: a bound that moved because the base under it changed.
    """
    plain = experiment_builder.recombined()
    mutated = experiment_builder.recombined().mutate(model="s5f", count=8)
    before = [o.airr_record(plain.refdata) for o in plain.run(n=40, seed=5)]
    after = [o.airr_record(mutated.refdata) for o in mutated.run(n=40, seed=5)]
    moved = [(index, a["junction_start"], b["junction_start"],
              a["junction_end"], b["junction_end"])
             for index, (a, b) in enumerate(zip(before, after))
             if (a["junction_start"], a["junction_end"])
             != (b["junction_start"], b["junction_end"])]
    assert moved == [], f"{len(moved)} records moved a bound under SHM: {moved[:3]}"
    assert any(b["n_mutations"] for b in after), "no mutation was applied"
