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


@verifies("assay.ambiguous_calls",
          claim="the pass writes N and only N, at positions that carried a "
                "base before it ran, and changes nothing else: against the "
                "same seed without the pass the molecule is the same length "
                "and every differing position is an N. No N appears without "
                "the pass at all, so every one in the sequence came from it",
          fields=("sequence",),
          basis="cross-run")
def test_ambiguous_calls_write_n_and_change_nothing_else():
    """The existing test counts Ns and checks none was called a mutation,
    which an N written anywhere -- or a molecule mangled alongside -- would
    also satisfy.

    This diffs against the same seed without the pass, so the expected value
    comes from a second run rather than from the record. Measured: 40 of 40
    records identical in length, every differing position an N, and 0 of 40
    plain records containing an N anywhere.

    The count requested is an upper bound rather than an equality -- 5 or 6 of
    6 here -- because a draw can land on a position already masked, the same
    repeat-hit shape the mutation and error counters have.
    """
    requested = 6
    plain = experiment_builder.recombined()
    masked = experiment_builder.recombined().ambiguous_base_calls(count=requested)
    before = [outcome.airr_record(plain.refdata)["sequence"]
              for outcome in plain.run(n=40, seed=5)]
    after = [outcome.airr_record(masked.refdata)["sequence"]
             for outcome in masked.run(n=40, seed=5)]
    total_masked = 0
    for index, (clean, ambiguous) in enumerate(zip(before, after)):
        assert "N" not in clean.upper(), (
            f"record {index} carries an N with no ambiguous-call pass "
            f"configured, so the Ns below are not all this pass's work")
        assert len(clean) == len(ambiguous), (
            f"record {index}: masking changed the molecule's length from "
            f"{len(clean)} to {len(ambiguous)}")
        changed = [position for position in range(len(clean))
                   if clean[position].upper() != ambiguous[position].upper()]
        assert changed, f"record {index}: the pass masked nothing"
        for position in changed:
            assert ambiguous[position].upper() == "N", (
                f"record {index} position {position} changed to "
                f"{ambiguous[position]!r}, not to N")
        masked_count = ambiguous.upper().count("N")
        assert masked_count == len(changed), (
            f"record {index} carries {masked_count} Ns but only "
            f"{len(changed)} positions changed")
        assert masked_count <= requested, (
            f"record {index} masked {masked_count} positions for "
            f"{requested} requested")
        total_masked += masked_count
    assert total_masked < requested * 40, (
        "every draw masked a fresh position, so this fixture cannot tell a "
        "requested count from a realised one")


@verifies("assay.contamination",
          claim="is_contaminant says whether the READ is this rearrangement. "
                "On a flagged record 70 to 80 per cent of the templated "
                "positions differ from the germline the record names, because "
                "the molecule is background DNA; on an unflagged one exactly "
                "none do. So the flag names a real substitution of the "
                "molecule rather than being set",
          fields=("is_contaminant",),
          basis="provenance")
def test_a_contaminant_read_is_not_the_rearrangement_it_names():
    """The existing test reads the boolean and checks that it varies, which
    its own claim admitted establishes nothing about which read carries it.

    The separation is wide enough to assert per record rather than
    distributionally: measured over 80 records at prob=0.5, the flagged ones
    run 0.698 to 0.804 and the unflagged are all exactly 0.000. The threshold
    below sits in the empty space between.

    Note what this does NOT use: the recorded germline provenance survives
    contamination untouched, so comparing it to the reference allele still
    agrees and reports nothing. What moves is the observed base.
    """
    exp = experiment_builder.recombined().contaminate(prob=0.5)
    flagged = unflagged = 0
    for index, outcome in enumerate(exp.run(n=80, seed=5)):
        record = outcome.airr_record(exp.refdata)
        simulation = outcome.final_simulation()
        bases, germline = simulation.bases(), simulation.germline_bases()
        templated = [position for position in range(len(bases))
                     if simulation.germline_position(position) is not None]
        assert templated, f"record {index} has no templated base to compare"
        divergence = sum(
            1 for position in templated
            if chr(bases[position]).upper() != chr(germline[position]).upper()
        ) / len(templated)
        if record["is_contaminant"]:
            flagged += 1
            assert divergence > 0.5, (
                f"record {index} is flagged as a contaminant but only "
                f"{divergence:.1%} of its templated bases differ from the "
                f"germline it names, so it is that rearrangement after all")
        else:
            unflagged += 1
            assert divergence == 0, (
                f"record {index} is not flagged yet {divergence:.1%} of its "
                f"templated bases differ from its own germline")
    assert flagged and unflagged, (
        f"the run produced {flagged} flagged and {unflagged} unflagged "
        f"records; both are needed for the comparison to mean anything")


PAIRED_END_CONFIGS = (
    ("fixed_insert", dict(r1_length=120, r2_length=120, insert_size=250)),
    ("drawn_insert", dict(r1_length=120, r2_length=120, insert_size=(200, 320))),
    ("asymmetric", dict(r1_length=150, r2_length=75, insert_size=(200, 300))),
)


#: A layout with no drawn component, so every window the engine reports is
#: determined by the configuration alone. The drawn layouts above cannot do
#: this: where the insert size is sampled per record, the only thing that
#: knows which value a record drew is the record.
FIXED_LAYOUT = dict(r1_length=120, r2_length=90, insert_size=250)


@verifies("assay.paired_end",
          claim="the read windows are where the declared layout puts them: R1 "
                "at 0..r1_length and R2 ending at the declared insert size, "
                "with the expected offsets computed from the configuration "
                "and not read back from the record. A layout that cannot fit "
                "the molecule tags read_layout and leaves the windows empty "
                "rather than looking like a record with no layout at all",
          fields=("r1_start", "r1_end", "r2_start", "r2_end"),
          basis="reference")
def test_the_read_windows_are_where_the_declared_layout_puts_them():
    """The offsets against the configuration, which is the only external side.

    The breadth test below compares each read to the molecule at the offsets
    the record itself reports, so a projection that placed both windows
    consistently wrong would satisfy it. Here the expected offsets come from
    the three numbers the caller declared -- 0..120 for R1, 160..250 for R2 --
    and the record has no say in them. Measured: 40 of 40 records agree.

    The second half is the geometry the engine cannot honour. An insert size
    larger than the molecule is rejected by the projection kernel, and the
    builder deliberately tags `read_layout` without populating the windows so
    the record stays distinguishable from one that never asked for reads. That
    path had no test on this side of the boundary.
    """
    layout = FIXED_LAYOUT
    exp = experiment_builder.recombined().paired_end(**layout)
    checked = 0
    for index, outcome in enumerate(exp.run(n=40, seed=7)):
        record = outcome.airr_record(exp.refdata)
        assert (record["r1_start"], record["r1_end"]) == (0, layout["r1_length"]), (
            f"record {index}: R1 is at {record['r1_start']}..{record['r1_end']}, "
            f"not at 0..{layout['r1_length']} as declared")
        expected_r2 = (layout["insert_size"] - layout["r2_length"],
                       layout["insert_size"])
        assert (record["r2_start"], record["r2_end"]) == expected_r2, (
            f"record {index}: R2 is at {record['r2_start']}..{record['r2_end']}, "
            f"not at {expected_r2[0]}..{expected_r2[1]} as declared")
        checked += 1
    assert checked == 40

    # A declared insert size no molecule can satisfy. 900 is well past the
    # longest IGH rearrangement this builder produces (404 bases on the
    # longest of the records above).
    impossible = experiment_builder.recombined().paired_end(
        r1_length=120, r2_length=90, insert_size=900)
    rejected = 0
    for index, outcome in enumerate(impossible.run(n=5, seed=7)):
        record = outcome.airr_record(impossible.refdata)
        assert record["read_layout"] == "paired_end", (
            f"record {index} asked for paired-end reads and the geometry was "
            f"rejected, but read_layout is {record['read_layout']!r} -- the "
            "record is now indistinguishable from one with no layout")
        assert record["r1_start"] is None and record["r2_end"] is None, (
            f"record {index}: the rejected geometry still produced windows "
            f"{record['r1_start']}..{record['r1_end']} and "
            f"{record['r2_start']}..{record['r2_end']}")
        rejected += 1
    assert rejected == 5


@verifies("assay.paired_end",
          claim="each read is the recorded molecule at the window the declared "
                "layout fixes: R1 the forward slice and R2 the reverse "
                "complement of its slice. Compared against the engine's own "
                "per-base record rather than against the exported sequence "
                "field, so a read sliced from a mis-exported molecule is "
                "caught",
          fields=("r1_sequence", "r2_sequence"),
          basis="provenance")
def test_each_read_is_the_recorded_molecule_at_the_declared_window():
    """Read content against the nucleotide record, not against a sibling field.

    The breadth test below slices ``record["sequence"]``, which is the same
    projection pass's output as the reads themselves -- so the two agree
    whether or not that molecule is what the engine built. This slices the
    recorded per-base molecule at offsets computed from the declared layout,
    which leaves nothing on the expected side that the paired-end projection
    produced.

    No orientation step: with a reverse-complement pass in the pipeline the
    exported molecule is the flipped read and the recorded pool is not, so the
    two sides would legitimately differ. Flipping is `assay.orientation`'s
    own claim.

    Measured: 40 of 40 records agree, and comparing R2 forward rather than
    reverse-complemented agrees on 0 of 40.
    """
    layout = FIXED_LAYOUT
    exp = experiment_builder.recombined().paired_end(**layout)
    r2_window = slice(layout["insert_size"] - layout["r2_length"],
                      layout["insert_size"])
    checked = 0
    for index, outcome in enumerate(exp.run(n=40, seed=7)):
        record = outcome.airr_record(exp.refdata)
        recorded = outcome.final_simulation().bases().decode()
        assert record["r1_sequence"].upper() == \
            recorded[0:layout["r1_length"]].upper(), (
            f"record {index}: R1 is not the recorded molecule's first "
            f"{layout['r1_length']} bases")
        assert record["r2_sequence"].upper() == reconstruction.reverse_complement(
            recorded[r2_window]).upper(), (
            f"record {index}: R2 is not the reverse complement of the recorded "
            f"molecule at {r2_window.start}..{r2_window.stop}")
        assert record["r2_sequence"].upper() != recorded[r2_window].upper(), (
            f"record {index}: R2 reads the same forward and reverse, so the "
            "orientation half of this test establishes nothing")
        checked += 1
    assert checked == 40


@verifies("assay.paired_end",
          claim="each read is the stretch of the molecule its own offsets "
                "name, across a fixed insert size, a drawn one and an "
                "asymmetric pair. Both sides of this comparison are exported "
                "by the same projection, so it establishes that the reads and "
                "the molecule agree and not where either came from -- the two "
                "claims above do that for one declared layout",
          fields=("r1_sequence", "r2_sequence", "r1_start", "r1_end",
                  "r2_start", "r2_end"),
          basis="consistent")
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

    What this does NOT establish, and what it was registered as `reference`
    for three weeks claiming: where the reads or the molecule came from. Both
    sides are the same projection's output, so a window placed consistently
    wrong satisfies it. That is why the two tests above exist, and why this
    one is now registered on the basis it has.
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


@verifies("assay.paired_end",
          claim="insert_size is the fragment length the caller declared: "
                "inside a declared range on every record and reaching both of "
                "its endpoints, or constant where a constant was declared. "
                "The expected values come from the configuration, which is "
                "the one part of a read layout that is not the projection's "
                "own output",
          fields=("insert_size",),
          basis="reference")
@pytest.mark.parametrize("label,layout", PAIRED_END_CONFIGS,
                         ids=[name for name, _ in PAIRED_END_CONFIGS])
def test_insert_sizes_follow_the_declared_layout(label, layout):
    """Registered on its own rather than tacked onto the test above.

    It was the last assertion in the content test, which meant one field with
    an external expectation was carrying the `reference` label for six fields
    that had none. Splitting it is most of why that row's basis was wrong.
    """
    exp = experiment_builder.recombined().paired_end(**layout)
    _insert_sizes_follow_the_declared_layout(exp, layout)


def _insert_sizes_follow_the_declared_layout(exp, layout):
    """``insert_size`` against the fragment length the caller declared.

    The expected values come from the configuration rather than from the
    record, which is what makes this more than an echo -- and an echo is
    exactly what this field was once filed as. A declared range must be
    respected on every record and both of its endpoints must be reachable; a
    declared constant must stay constant.
    """
    declared = layout["insert_size"]
    sizes = [outcome.airr_record(exp.refdata)["insert_size"]
             for outcome in exp.run(n=200, seed=7)]
    assert sizes, "no record carried an insert size"
    if isinstance(declared, tuple):
        low, high = declared
        outside = [size for size in sizes if not low <= size <= high]
        assert not outside, (
            f"{len(outside)} insert sizes fall outside the declared "
            f"{low}..{high}: {sorted(set(outside))[:5]}")
        assert min(sizes) == low and max(sizes) == high, (
            f"the declared range {low}..{high} was never reached: observed "
            f"{min(sizes)}..{max(sizes)}, so the bound is untested")
    else:
        assert set(sizes) == {declared}, (
            f"a declared constant insert size of {declared} varied: "
            f"{sorted(set(sizes))[:5]}")


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


def _indel_counter_fields_present(record):
    """The four indel counters, named literally."""
    for field in ("n_indels", "n_v_indels", "n_d_indels", "n_j_indels"):
        assert isinstance(record[field], int), f"{field} is not an integer"


@verifies("assay.polymerase_indels",
          claim="n_indels counts the indels the pass APPLIED, not the number "
                "drawn: it equals the net change in the molecule's length "
                "against the same seed without the pass, 120 of 120 in each "
                "direction. The three per-segment counters sum to at most it, "
                "and to exactly it when no NP region exists -- an indel in an "
                "NP region belongs to no segment counter, which the "
                "mechanism's truth gap records",
          fields=("n_indels", "n_v_indels", "n_d_indels", "n_j_indels"),
          basis="cross-run")
@pytest.mark.parametrize("label,insertion_prob", [("insertions", 1.0),
                                                  ("deletions", 0.0)])
def test_the_indel_count_is_the_length_the_molecule_changed_by(label, insertion_prob):
    """A cross-run differential, which is the only way to see this: the record
    cannot say how long the molecule was before the pass ran.

    Worth checking at all because the neighbouring error counters turned out
    to report the number of errors REQUESTED rather than applied. This one
    does not -- it was already read off the event ledger -- and now something
    says so.
    """
    plain = experiment_builder.recombined()
    with_indels = experiment_builder.recombined().polymerase_indels(
        count=(1, 3), insertion_prob=insertion_prob)
    before = [outcome.airr_record(plain.refdata)["sequence"]
              for outcome in plain.run(n=120, seed=5)]
    after = [outcome.airr_record(with_indels.refdata)
             for outcome in with_indels.run(n=120, seed=5)]
    moved = 0
    for index, (clean, record) in enumerate(zip(before, after)):
        _indel_counter_fields_present(record)
        delta = len(record["sequence"]) - len(clean)
        expected = delta if insertion_prob == 1.0 else -delta
        assert record["n_indels"] == expected, (
            f"{label}: record {index} reports {record['n_indels']} indels "
            f"where the molecule's length moved by {delta}")
        per_segment = (record["n_v_indels"] + record["n_d_indels"]
                       + record["n_j_indels"])
        assert per_segment <= record["n_indels"], (
            f"{label}: record {index} attributes {per_segment} indels to "
            f"segments out of {record['n_indels']}")
        moved += record["n_indels"] > 0
    assert moved > 100, f"{label}: only {moved} records had an indel applied"


@verifies("assay.polymerase_indels",
          claim="the per-segment counters account for every indel once no NP "
                "region exists: with np1 and np2 pinned to zero the three sum "
                "to n_indels on 200 of 200 records, where at the default NP "
                "lengths 18 of 200 carry an indel that belongs to no segment",
          fields=("n_v_indels", "n_d_indels", "n_j_indels"),
          basis="provenance")
def test_segment_indels_account_for_everything_when_there_is_no_np_region():
    """Identifies the causes of the shortfall rather than tolerating it.

    Comparing the two configurations is what distinguishes "the per-segment
    attribution loses indels" from the two places an indel legitimately
    belongs to no segment: an NP region, which has no counter of its own, and
    a site past the last region, where an insertion lands outside every
    segment and no alignment has a column for it.

    Both are documented gaps rather than defects, and the second one is the
    narrower: it is reachable only at the molecule's 3' end, on about 1% of
    insertion-bearing records. See
    tests/biology/recombination/test_boundary_insertions.py for why the base
    is left there instead of being drawn into the J.
    """
    without_np = experiment_builder.recombined(
        np1_lengths=[(0, 1.0)], np2_lengths=[(0, 1.0)]
    ).polymerase_indels(count=(1, 3), insertion_prob=0.5)
    exact = outside = 0
    for index, outcome in enumerate(without_np.run(n=200, seed=5)):
        record = outcome.airr_record(without_np.refdata)
        per_segment = (record["n_v_indels"] + record["n_d_indels"]
                       + record["n_j_indels"])
        trailing = record["sequence_length"] - record["j_sequence_end"]
        if trailing:
            # An insertion past the J: counted in n_indels, in no per-segment
            # counter, and in no cigar -- the three agreeing is the contract.
            outside += 1
            assert per_segment == record["n_indels"] - trailing, (
                f"record {index}: {per_segment} attributed of "
                f"{record['n_indels']} with {trailing} base(s) past the J")
            assert "I" not in (record["j_cigar"] or ""), (
                f"record {index}: the J cigar claims an insertion the J's "
                "span does not cover")
            continue
        assert per_segment == record["n_indels"], (
            f"record {index}: {per_segment} attributed of "
            f"{record['n_indels']} with no NP region to lose one to")
        exact += 1
    assert exact + outside == 200
    assert exact > 150, (
        f"only {exact} of 200 records attributed every indel to a segment, "
        "which is too few for the comparison below to be about NP regions")

    with_np = experiment_builder.recombined().polymerase_indels(
        count=(1, 3), insertion_prob=0.5)
    unattributed = sum(
        1 for outcome in with_np.run(n=200, seed=5)
        for record in [outcome.airr_record(with_np.refdata)]
        if (record["n_v_indels"] + record["n_d_indels"] + record["n_j_indels"])
        < record["n_indels"])
    assert unattributed, (
        "no record lost an indel to an NP region, so this comparison does not "
        "establish where the shortfall comes from")
