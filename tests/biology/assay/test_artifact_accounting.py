"""Each assay artifact is counted as itself and not as biology.

This is the distinction that matters most to a benchmark: a PCR
misincorporation and a somatic hypermutation both look like a difference from
germline, and a tool scored against GenAIRR truth is only fair if the record
says which was which.
"""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import reconstruction

pytestmark = pytest.mark.behavior


@verifies("assay.pcr_error",
          claim="amplification error is counted in n_pcr_errors and not folded "
                "into the somatic mutation count",
          fields=("n_pcr_errors",),
          basis="exercised")
def test_pcr_error_is_counted_separately_from_hypermutation():
    exp = experiment_builder.recombined().pcr_amplify(count=3)
    seen = 0
    for outcome in exp.run(n=40, seed=5):
        record = outcome.airr_record(exp.refdata)
        assert record["n_mutations"] == 0, (
            "no SHM was configured, so a non-zero mutation count means an "
            "amplification error was counted as somatic")
        seen += record["n_pcr_errors"]
    assert seen > 0, "no PCR error was produced; the fixture establishes nothing"


@verifies("assay.quality_error",
          claim="base-calling error is counted in n_quality_errors, separately "
                "from both amplification error and somatic mutation",
          fields=("n_quality_errors",),
          basis="exercised")
def test_quality_error_is_counted_separately():
    exp = experiment_builder.recombined().sequencing_errors(count=3)
    seen = 0
    for outcome in exp.run(n=40, seed=5):
        record = outcome.airr_record(exp.refdata)
        assert record["n_mutations"] == 0 and record["n_pcr_errors"] == 0
        seen += record["n_quality_errors"]
    assert seen > 0, "no quality error was produced; the fixture establishes nothing"


@verifies("assay.ambiguous_calls",
          claim="an uncalled position appears as N and is not counted as a "
                "mutation -- the claim the mechanism's truth_gap says has no "
                "field of its own to carry it",
          fields=("sequence",),
          basis="exercised")
def test_ambiguous_positions_are_not_counted_as_mutations():
    """The mechanism declares a truth gap because no N-count field exists. The
    claim can still be established against the sequence itself."""
    exp = experiment_builder.recombined().ambiguous_base_calls(count=8)
    saw_n = False
    for outcome in exp.run(n=40, seed=5):
        record = outcome.airr_record(exp.refdata)
        ns = record["sequence"].upper().count("N")
        saw_n |= ns > 0
        assert record["n_mutations"] == 0, (
            f"{ns} ambiguous positions were counted as {record['n_mutations']} "
            "mutations; an N is an absence of information, not a substitution")
    assert saw_n, "no ambiguous base was produced; the fixture establishes nothing"


@verifies("assay.orientation",
          claim="rev_comp discriminates: it is False when no strand was "
                "flipped, True when one was, and both occur in a mixed run. "
                "When True the recorded sequence is the reverse complement of "
                "the sequence the same seed produces unflipped, so the flag "
                "names a real transformation rather than being set",
          fields=("rev_comp",),
          basis="cross-run")
def test_orientation_is_recorded_and_discriminates():
    """AIRR puts positions with respect to the reverse complement of `sequence`
    when rev_comp is true, so the flag is what makes every coordinate in the
    record interpretable.

    Checked on both sides. Run only at prob=1.0, this test passes against a
    field hardcoded to True -- which is no check on a boolean at all. The
    reverse-complement comparison against the same seed unflipped is what makes
    the True case say something.
    """
    forward = experiment_builder.recombined().random_strand_orientation(prob=0.0)
    flipped = experiment_builder.recombined().random_strand_orientation(prob=1.0)
    unflipped = [o.airr_record(forward.refdata) for o in forward.run(n=20, seed=5)]
    reversed_ = [o.airr_record(flipped.refdata) for o in flipped.run(n=20, seed=5)]

    for record in unflipped:
        assert record["rev_comp"] is False, "no strand was flipped, yet rev_comp is set"
    for plain, flip in zip(unflipped, reversed_):
        assert flip["rev_comp"] is True
        assert flip["sequence"].upper() == reconstruction.reverse_complement(plain["sequence"]), (
            "rev_comp is True but the sequence is not the reverse complement of "
            "what the same seed produces unflipped")
        assert 0 <= flip["v_sequence_start"] < len(flip["sequence"])

    mixed = experiment_builder.recombined().random_strand_orientation(prob=0.5)
    flags = [o.airr_record(mixed.refdata)["rev_comp"] for o in mixed.run(n=40, seed=5)]
    assert set(flags) == {True, False}, (
        f"a half-and-half run produced only {set(flags)}; the flag does not vary")


@verifies("assay.contamination",
          claim="is_contaminant discriminates: False on an uncontaminated run "
                "and on one with the pass configured off, True when every read "
                "is contaminated, and both occur in a mixed run. Which read "
                "carries the flag in a mixed run is not established here -- the "
                "body reads the boolean, not the read's provenance",
          fields=("is_contaminant",),
          basis="exercised")
def test_contamination_is_flagged_and_discriminates():
    """Both sides, for the same reason as orientation above: at prob=1.0 alone,
    a field hardcoded to True passes."""
    always = experiment_builder.recombined().contaminate(prob=1.0)
    assert {o.airr_record(always.refdata)["is_contaminant"]
            for o in always.run(n=20, seed=5)} == {True}, (
        "contaminated reads were not all flagged")

    never = experiment_builder.recombined().contaminate(prob=0.0)
    assert {o.airr_record(never.refdata)["is_contaminant"]
            for o in never.run(n=20, seed=5)} == {False}, (
        "reads were flagged as contaminants with the pass configured off")

    absent = experiment_builder.recombined()
    assert {o.airr_record(absent.refdata)["is_contaminant"]
            for o in absent.run(n=20, seed=5)} == {False}, (
        "reads were flagged as contaminants with no contamination pass at all")

    mixed = experiment_builder.recombined().contaminate(prob=0.5)
    flags = [o.airr_record(mixed.refdata)["is_contaminant"] for o in mixed.run(n=40, seed=5)]
    assert set(flags) == {True, False}, (
        f"a half-and-half run produced only {set(flags)}; the flag does not vary")


@verifies("assay.paired_end",
          claim="each read's recorded span is the length of the read and lies "
                "within the molecule, and R2 is positioned by the fragment "
                "length actually drawn for that record -- r2_end == insert_size "
                "and r2_start == insert_size - r2_length -- so a tool "
                "reconstructing the fragment from R1 and R2 is scored against "
                "coordinates that hold",
          fields=("r1_sequence", "r1_start", "r1_end", "r2_sequence",
                  "r2_start", "r2_end", "insert_size"),
          basis="consistent")
def test_each_read_is_the_stretch_its_offsets_claim():
    """Fields are read by name rather than through an f-string, so a static
    check can confirm the claim matches what the body actually inspects.

    insert_size is drawn from a distribution, not pinned. Fixed at a constant
    the relation below is satisfied by a field echoing the configuration back,
    which is how insert_size came to be filed as an echo in the first place --
    so the fixture draws over a range and the test requires it to have varied.
    """
    exp = (experiment_builder.recombined()
           .paired_end(r1_length=120, insert_size=(200, 320)))
    checked = 0
    drawn = set()
    for outcome in exp.run(n=25, seed=7):
        record = outcome.airr_record(exp.refdata)
        molecule = len(record["sequence"])
        for sequence, start, end, which in (
                (record["r1_sequence"], record["r1_start"], record["r1_end"], "r1"),
                (record["r2_sequence"], record["r2_start"], record["r2_end"], "r2")):
            if not sequence or start is None or end is None:
                continue
            assert end - start == len(sequence), (
                f"{which}: offsets span {end - start} but the read is "
                f"{len(sequence)} bases")
            assert 0 <= start <= end <= molecule, (
                f"{which}: offsets {start}..{end} fall outside a "
                f"{molecule}-base molecule")
            checked += 1
        # The relation the mechanism's own summary names: the insert size is
        # what places R2. Every assertion above holds with R2 sampled from a
        # different fragment length entirely.
        size = record["insert_size"]
        drawn.add(size)
        assert record["r2_end"] == size, (
            f"r2 ends at {record['r2_end']} but the fragment drawn for this "
            f"record is {size} bases")
        assert record["r2_start"] == size - len(record["r2_sequence"]), (
            f"r2 starts at {record['r2_start']}, not {size} minus its own "
            f"{len(record['r2_sequence'])} bases")
    assert checked, "no paired read was produced; the fixture establishes nothing"
    assert len(drawn) > 1, (
        f"insert_size took one value ({drawn}); the relation above would be "
        "satisfied by a field echoing the configuration back")


PAIRED_END_CONFIGS = (
    ("fixed_insert", dict(r1_length=120, r2_length=120, insert_size=250)),
    ("drawn_insert", dict(r1_length=120, r2_length=120, insert_size=(200, 320))),
    ("asymmetric", dict(r1_length=150, r2_length=75, insert_size=(200, 300))),
)


@verifies("assay.paired_end",
          claim="each read IS the stretch of the molecule its offsets name: R1 "
                "the forward slice and R2 the reverse complement of its slice, "
                "which is the orientation a paired-end consumer assembles "
                "from. Read against the molecule rather than against its own "
                "length, so a read whose content came from somewhere else is "
                "caught",
          fields=("r1_sequence", "r2_sequence", "r1_start", "r1_end",
                  "r2_start", "r2_end"),
          basis="reference")
@pytest.mark.parametrize("label,layout", PAIRED_END_CONFIGS,
                         ids=[name for name, _ in PAIRED_END_CONFIGS])
def test_each_read_is_the_molecule_at_its_own_offsets(label, layout):
    """What the existing span test could not see.

    That one checks each read's offsets against its own length and the
    molecule's bounds -- true of any read of the right size from anywhere,
    including one assembled from the wrong fragment. The content was never
    compared to the molecule at all.

    R2's orientation is the half that discriminates: comparing it forward
    rather than reverse-complemented matches on 0 of 40 records, so a build
    that stopped reverse-complementing R2 would fail here and pass every
    existing check.
    """
    exp = experiment_builder.recombined().paired_end(**layout)
    checked = 0
    for index, outcome in enumerate(exp.run(n=40, seed=7)):
        record = outcome.airr_record(exp.refdata)
        molecule = record["sequence"]
        assert record["r1_sequence"].upper() == \
            molecule[record["r1_start"]:record["r1_end"]].upper(), (
            f"{label}: record {index} R1 is not the molecule at "
            f"{record['r1_start']}..{record['r1_end']}")
        assert record["r2_sequence"].upper() == reconstruction.reverse_complement(
            molecule[record["r2_start"]:record["r2_end"]]).upper(), (
            f"{label}: record {index} R2 is not the reverse complement of the "
            f"molecule at {record['r2_start']}..{record['r2_end']}")
        assert record["r2_sequence"].upper() != \
            molecule[record["r2_start"]:record["r2_end"]].upper(), (
            f"{label}: record {index} R2 reads the same forward and reverse, "
            f"so the orientation half of this test establishes nothing")
        checked += 1
    assert checked == 40


@verifies("assay.end_loss",
          claim="5' end loss is recorded in its own field, so sequence missing "
                "because the read did not reach it is distinguishable from "
                "germline trimming",
          fields=("end_loss_5_length",),
          basis="exercised")
def test_end_loss_is_recorded_separately_from_trimming():
    exp = experiment_builder.recombined().end_loss_5prime(length=12)
    lost = 0
    for outcome in exp.run(n=30, seed=7):
        record = outcome.airr_record(exp.refdata)
        assert record["v_trim_5"] == 0, (
            "end loss was folded into the trim fields; a read that did not "
            "reach the 5' end is not exonuclease loss")
        lost += record["end_loss_5_length"]
    assert lost > 0, "no end loss was produced; the fixture establishes nothing"


def _positions_the_pass_changed(add_pass, *, n, seed, count):
    """Per record, the positions a corruption pass actually changed.

    A cross-run differential, not an oracle: it runs the library twice with
    one seed, once with the pass and once without, and diffs the sequences.
    Adding a corruption pass does not disturb the recombination draws that
    precede it, so every position that differs between the two runs is a base
    this pass wrote and changed. Measured on this fixture: 3 to 5 differing
    positions of ~380, the rest byte-identical.

    This cannot be an oracle under ``tests/support/oracles`` because it has to
    call GenAIRR to produce the second dataset, which an oracle may never do.
    What makes it independent of the thing under test is that it never reads
    the event ledger or the trace the counters are derived from -- only two
    emitted sequences.
    """
    plain = [o.airr_record(experiment_builder.recombined().refdata)["sequence"]
             for o in experiment_builder.recombined().run(n=n, seed=seed)]
    corrupt_exp = add_pass(experiment_builder.recombined(), count)
    corrupted = [o.airr_record(corrupt_exp.refdata)["sequence"]
                 for o in corrupt_exp.run(n=n, seed=seed)]
    changed = []
    for before, after in zip(plain, corrupted):
        assert len(before) == len(after), (
            "a substitution pass must not change the molecule's length")
        changed.append(sum(1 for a, b in zip(before, after)
                           if a.upper() != b.upper()))
    return changed, corrupt_exp


def _realized_and_reported(add_pass, field, *, n=200, seed=11, count=1):
    """``(changed positions per record, counter per record)`` for one pass."""
    changed, exp = _positions_the_pass_changed(add_pass, n=n, seed=seed, count=count)
    reported = [o.airr_record(exp.refdata)[field] for o in exp.run(n=n, seed=seed)]
    return changed, reported


def _assert_counts_realized_not_drawn(field, changed, reported, *, draws):
    assert sum(changed) > 0, f"{field}: the fixture produced no error at all"
    for index, (got, expected) in enumerate(zip(reported, changed)):
        assert got == expected, (
            f"{field}: record {index} reports {got} where the pass changed "
            f"{expected} positions")
    # The gap is the point: if every draw changed a base the two numbers would
    # coincide and this would establish nothing about which one is counted.
    assert sum(changed) < draws, (
        f"{field}: every one of {draws} draws changed a base, so a draw count "
        f"and a realized count are indistinguishable on this fixture")


@verifies("assay.pcr_error",
          claim="n_pcr_errors counts the amplification errors that changed a "
                "base, not the number of error draws. A replacement is drawn "
                "from A/C/G/T including the base already present, so a draw "
                "landing on it writes a byte and introduces no error",
          fields=("n_pcr_errors",),
          basis="cross-run")
def test_pcr_error_counter_counts_realized_errors_not_draws():
    """The same defect `n_mutations` carried: a counter read off the drawn
    count reports work that did not happen, and a tool scored against it
    carries a built-in false-negative rate equal to the discrepancy.

    At one error per record a repeat hit is impossible, so the realized count
    and the number of changed positions are the same number and can be
    compared exactly.
    """
    changed, reported = _realized_and_reported(
        lambda e, c: e.pcr_amplify(count=c), "n_pcr_errors")
    _assert_counts_realized_not_drawn(
        "n_pcr_errors", changed, reported, draws=200)


@verifies("assay.quality_error",
          claim="n_quality_errors counts the base-calling errors that changed "
                "a base, not the number of error draws. It is populated by "
                "sequencing_errors(); ambiguous_base_calls() masks positions "
                "with N and reports on neither error counter",
          fields=("n_quality_errors",),
          basis="cross-run")
def test_quality_error_counter_counts_realized_errors_not_draws():
    """Same claim for the base-calling counter, which shares the write path
    and shared the defect."""
    changed, reported = _realized_and_reported(
        lambda e, c: e.sequencing_errors(count=c), "n_quality_errors")
    _assert_counts_realized_not_drawn(
        "n_quality_errors", changed, reported, draws=200)
