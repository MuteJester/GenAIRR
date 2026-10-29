"""Would the oracle catch a regression? Measured, not assumed.

Everything else here asserts that GenAIRR's recorded truth is correct today. It
does not establish that a test would notice if that stopped being true -- and a
verification that cannot fail is indistinguishable from no verification.

So: take a correct recorded result, corrupt it in a specific, biologically
meaningful way, and assert the oracle reports it. Each mutant corresponds to a
defect the engine could plausibly acquire -- a swapped junction pair, a
coordinate off by one, an allele assembled from the wrong reference, a mutation
counted that changed nothing.

These are record-level mutants rather than engine-level ones. An engine mutant
is stronger, because it exercises the real path, but it needs a rebuild per
mutant and so cannot run in ordinary CI. This runs on every commit.
"""
from __future__ import annotations

from dataclasses import dataclass

import pytest

from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import coordinates, projection, reconstruction

pytestmark = pytest.mark.pin


class _Corrupted:
    """A Simulation view with one aspect deliberately wrong.

    Wraps rather than mutates, because the real object is read-only and because
    a wrapper makes the single corrupted aspect obvious at the call site.
    """

    def __init__(self, real, *, bases=None, germline=None, shift=0, drop_positions=()):
        self._real = real
        self._bases = bases
        self._germline = germline
        self._shift = shift
        self._drop = set(drop_positions)

    def bases(self):
        return self._bases if self._bases is not None else self._real.bases()

    def germline_bases(self):
        return self._germline if self._germline is not None else self._real.germline_bases()

    def germline_position(self, index):
        if index in self._drop:
            return None
        position = self._real.germline_position(index)
        if position is None:
            return None
        return position + self._shift

    def regions(self):
        return self._real.regions()

    def v_allele_id(self):
        return self._real.v_allele_id()

    def d_allele_id(self):
        return self._real.d_allele_id()

    def j_allele_id(self):
        return self._real.j_allele_id()


@dataclass(frozen=True)
class Mutant:
    name: str
    what_it_models: str


def _correct():
    exp = experiment_builder.recombined()
    return exp, exp.run(n=1, seed=0)[0]


def test_the_oracle_passes_on_an_uncorrupted_record():
    """The control. Without it, a mutant test that always fails would look like
    a working sensitivity check."""
    exp, outcome = _correct()
    check = reconstruction.templated_check(outcome.final_simulation(), exp.refdata)
    assert check.mismatches == ()
    assert check.clean, "nothing was compared, so agreeing means nothing"


def test_a_verdict_reached_without_comparing_anything_is_distinguishable():
    """Models the most dangerous shape of a passing test: one whose oracle was
    handed a record it could not read.

    An empty mismatch tuple describes both "every position agreed" and "no
    position was examined". Dropping every germline coordinate produces the
    second, and the result must not be mistakable for the first.
    """
    exp, outcome = _correct()
    simulation = outcome.final_simulation()
    blinded = _Corrupted(simulation, drop_positions=range(len(simulation.bases())))
    check = reconstruction.templated_check(blinded, exp.refdata)
    assert check.mismatches == (), "the blinded record has nothing to disagree about"
    assert check.checked == 0, "a record with no germline coordinates compared something"
    assert not check.clean, "an unexamined record reported itself clean"


def test_a_germline_coordinate_shifted_by_one_is_caught():
    """Models an off-by-one in a coordinate the engine records -- the defect
    class a benchmark would silently inherit as a systematic one-base offset."""
    exp, outcome = _correct()
    corrupted = _Corrupted(outcome.final_simulation(), shift=1)
    check = reconstruction.templated_check(corrupted, exp.refdata)
    assert check.mismatches, "a one-base coordinate shift went undetected"
    assert check.checked > 50, "too few positions compared for this to mean much"


def test_an_assembly_from_the_wrong_allele_is_caught():
    """Models the engine recording one allele while having copied another.

    The engine's own provenance would be internally consistent, so only a check
    against the reference tables can see it -- which is why the oracle reads
    the reference rather than the record alone.
    """
    exp, outcome = _correct()
    simulation = outcome.final_simulation()
    germline = bytearray(simulation.germline_bases())
    # rewrite the V region's germline as a different allele's sequence
    pool = exp.refdata.v_pool_size
    pool = pool() if callable(pool) else pool
    other_id = (simulation.v_allele_id() + 1) % max(1, pool)
    other = exp.refdata.v_allele(other_id)
    other_sequence = other.seq() if callable(other.seq) else other.seq
    changed = 0
    for index in range(len(germline)):
        position = simulation.germline_position(index)
        if position is None or position >= len(other_sequence):
            continue
        if germline[index] != other_sequence[position]:
            germline[index] = other_sequence[position]
            changed += 1
    if changed == 0:
        pytest.skip("the neighbouring allele is identical over the assembled span")
    corrupted = _Corrupted(simulation, germline=bytes(germline))
    assert reconstruction.templated_check(corrupted, exp.refdata).mismatches, (
        "germline rewritten from the wrong allele went undetected")


def test_a_swapped_np_pair_is_caught():
    """Models np1 and np2 being wired the wrong way round -- a defect that
    leaves both fields populated and plausible.

    Run against the same comparison the junction truth test uses, not a
    restatement of it: a mutant checked by reimplemented logic proves only
    that two pieces of test code agree with each other.
    """
    exp, outcome = _correct()
    simulation = outcome.final_simulation()
    runs = reconstruction.untemplated_runs(simulation)
    record = outcome.airr_record(exp.refdata)
    if reconstruction.np_fields_disagreement(record, runs) is not None:
        pytest.fail("the control record already disagrees; the fixture is wrong")
    if len(runs) != 2 or runs[0].bases == runs[1].bases:
        pytest.skip("this record has no distinguishable NP pair to swap")
    swapped = {**record, "np1": record["np2"], "np2": record["np1"],
               "np1_length": record["np2_length"], "np2_length": record["np1_length"]}
    assert reconstruction.np_fields_disagreement(swapped, runs) is not None, (
        "a swapped NP pair matched the untemplated runs, so the junction claim "
        "cannot tell the two apart")


def test_a_mutation_counted_that_changed_nothing_is_caught():
    """Models exactly the defect this branch found and pinned: a recorded
    substitution that left the base unchanged."""
    exp = experiment_builder.recombined().mutate(model="s5f", count=2)
    outcome = exp.run(n=1, seed=0)[0]
    simulation = outcome.final_simulation()
    # make the final bases identical to germline: every recorded mutation
    # becomes a no-op, which is what the uniform model does 27% of the time
    record = outcome.airr_record(exp.refdata)
    assert record["n_mutations"] > 0, "the fixture produced no mutations to contradict"
    # The control first: without it, an oracle that always returned () would
    # pass the mutant below and look like a working check.
    assert reconstruction.observable_substitutions(simulation), (
        "the uncorrupted record shows no observable substitution, so finding "
        "none in the mutant would prove nothing")
    corrupted = _Corrupted(simulation, bases=simulation.germline_bases())
    assert reconstruction.observable_substitutions(corrupted) == (), (
        f"{record['n_mutations']} mutations are recorded but none of them "
        "changed a base, and the oracle did not notice")


def test_a_p_nucleotide_that_is_not_a_palindrome_is_caught():
    """Models P addition emitting arbitrary bases instead of the reverse
    complement of the segment tip.

    The provenance oracle cannot see this on its own: it records P with no
    germline coordinate, exactly as it records N. Run against a real record
    carrying a real P addition, with the P run overwritten in place, so what
    is measured is the check's sensitivity and not two string literals.
    """
    exp = experiment_builder.with_p_nucleotides()
    for outcome in exp.run(n=25, seed=1):
        record = outcome.airr_record(exp.refdata)
        boundary, p_length = record["v_sequence_end"], record["p_v_3_length"]
        if not p_length or boundary is None:
            continue
        sequence = outcome.final_simulation().bases().decode()
        assert reconstruction.p_tip_disagreement(sequence, boundary, p_length) is None, (
            "the control record's P run is not a tip palindrome")
        # Overwrite the P run with bases that cannot be its tip's reverse
        # complement, leaving everything else as recorded.
        tip = sequence[boundary - p_length:boundary]
        forged = next(base * p_length for base in "ACGT"
                      if base * p_length != reconstruction.reverse_complement(tip))
        mutant = sequence[:boundary] + forged + sequence[boundary + p_length:]
        assert len(mutant) == len(sequence)
        assert reconstruction.p_tip_disagreement(mutant, boundary, p_length) is not None, (
            f"a P run forged as {forged!r} over tip {tip!r} went undetected")
        return
    pytest.fail("no record carried a V 3' P nucleotide; the fixture is wrong")


# ── the coordinate oracle ────────────────────────────────────────────────
#
# Mutants are applied to the exported record rather than the simulation,
# because the coordinate oracle reads the record. Each one models a defect in
# a different producer: the trimming pass, the record builder's germline
# export, its read export, and its alignment export.


def _coordinate_case():
    """A record with no indels, plus the allele lengths read from the reference."""
    exp = experiment_builder.recombined()
    outcome = exp.run(n=1, seed=0)[0]
    simulation = outcome.final_simulation()
    getters = {"v": exp.refdata.v_allele, "d": exp.refdata.d_allele, "j": exp.refdata.j_allele}
    identifiers = {"v": simulation.v_allele_id(), "d": simulation.d_allele_id(),
                   "j": simulation.j_allele_id()}
    lengths = {}
    for segment, identifier in identifiers.items():
        if identifier is None:
            continue
        allele = getters[segment](identifier)
        sequence = allele.seq() if callable(allele.seq) else allele.seq
        lengths[segment] = len(sequence)
    record = dict(outcome.airr_record(exp.refdata))
    assert (record.get("n_indels") or 0) == 0, "the fixture must be indel-free"
    return record, lengths


def test_the_coordinate_oracle_passes_on_an_uncorrupted_record():
    """The control, and the coverage check: a record whose coordinates were
    never examined would also report no disagreements."""
    record, lengths = _coordinate_case()
    assert coordinates.disagreements(record, lengths) == ()
    assert set(lengths) == {"v", "d", "j"}, (
        "the fixture must present all three segments, or the mutants below "
        "would be checking fewer identities than they claim")


@pytest.mark.parametrize("field, delta, identity", [
    ("v_trim_3", -1, "germline_end == allele_length - trim_3 - end_loss_3_length"),
    ("d_trim_5", 1, "germline_start == trim_5 + end_loss_5_length"),
    ("j_germline_start", 1, "germline_start == trim_5 + end_loss_5_length"),
    ("v_germline_end", -2, "germline_end == allele_length - trim_3 - end_loss_3_length"),
    ("j_sequence_end", 3, "germline span == read span"),
    ("d_alignment_start", 1, "alignment coordinates == read coordinates"),
])
def test_each_coordinate_identity_is_sensitive_to_its_own_producer(field, delta, identity):
    """Models a single producer disagreeing with the others by a few bases.

    This is the defect class the base-level oracle cannot see: shifting a trim
    amount or a germline coordinate leaves every recorded base equal to the
    reference base at the coordinate it claims, because production derives the
    byte and the coordinate from the same slice.
    """
    record, lengths = _coordinate_case()
    corrupted = {**record, field: record[field] + delta}
    found = coordinates.disagreements(corrupted, lengths)
    assert found, f"a {delta:+d} shift in {field} went undetected"
    assert any(disagreement.identity == identity for disagreement in found), (
        f"{field} was detected, but not by {identity!r}; "
        f"got {sorted({d.identity for d in found})}")


def test_a_reference_table_the_engine_disagrees_with_is_caught():
    """Models the engine assembling from one allele while the reference table
    publishes a different length for it -- a packaging defect that would make
    every downstream coordinate wrong while the record stayed self-consistent."""
    record, lengths = _coordinate_case()
    stretched = {**lengths, "v": lengths["v"] + 3}
    found = coordinates.disagreements(record, stretched)
    assert found, "a reference allele three bases longer than the engine used went undetected"
    assert {disagreement.segment for disagreement in found} == {"v"}


def test_an_indel_count_does_not_excuse_a_germline_bound_with_no_end_loss():
    """The first version of this test pinned the slack it warned against.

    It asserted that "two recorded indels should tolerate a two-base
    displacement" of v_germline_end -- granting exactly the licence its own
    docstring said must not exist. Measured, the germline bounds deviate only
    at an end where read truncation also applies; with no end loss they are
    exact however many indels the segment carries. So the tolerance is now
    conditioned on the end loss rather than on the indel count alone, and a
    one-base error here is a finding.
    """
    record, lengths = _coordinate_case()
    assert (record.get("end_loss_3_length") or 0) == 0, "the fixture must have no 3' end loss"
    excused = {**record, "n_v_indels": 3, "n_indels": 3,
               "v_germline_end": record["v_germline_end"] - 1}
    found = coordinates.disagreements(excused, lengths)
    assert found, (
        "a one-base error in v_germline_end was excused by three indels, on a "
        "record with no end loss at that end")
    assert any("germline_end ==" in d.identity for d in found)


def test_the_read_span_tolerance_is_bounded_by_the_recorded_indel_count():
    """The read span is the identity an indel genuinely displaces, so it does
    carry a tolerance -- and that tolerance must stay bounded."""
    record, lengths = _coordinate_case()
    within = {**record, "n_v_indels": 2, "n_indels": 2,
              "v_sequence_end": record["v_sequence_end"] + 2,
              "v_alignment_end": record["v_alignment_end"] + 2}
    assert any("span" in d.identity for d in coordinates.disagreements(within, lengths)) is False, (
        "two recorded indels should tolerate a two-base read-span displacement")
    beyond = {**within, "v_sequence_end": record["v_sequence_end"] + 3,
              "v_alignment_end": record["v_alignment_end"] + 3}
    assert coordinates.disagreements(beyond, lengths), (
        "a three-base read-span displacement with two recorded indels went undetected")


@pytest.mark.parametrize("segment, delta", [("v", 7), ("d", -4), ("j", 3)])
def test_a_segment_displaced_as_a_whole_is_caught(segment, delta):
    """The defect the difference-based identities cannot see.

    Shift a segment's read and alignment coordinates by the same amount and
    every span, every gap-free difference and every alignment-equals-read pair
    is preserved. Nothing anchored an absolute position, so this returned no
    disagreements -- while the twelve coordinates an aligner benchmark is
    scored on were all wrong by the same offset.

    The six per-field mutants above cannot reach this shape, because each of
    them moves exactly one field and is therefore caught by a difference.
    """
    record, lengths = _coordinate_case()
    displaced = dict(record)
    for suffix in ("sequence_start", "sequence_end", "alignment_start", "alignment_end"):
        displaced[f"{segment}_{suffix}"] = record[f"{segment}_{suffix}"] + delta
    found = coordinates.disagreements(displaced, lengths)
    assert found, (
        f"{segment.upper()} displaced by {delta:+d} as a whole went undetected")


def test_the_whole_rearrangement_displaced_is_caught():
    """Every segment moved together, which preserves even the inter-segment
    gaps. Only the two outer anchors can see this one."""
    record, lengths = _coordinate_case()
    displaced = dict(record)
    for segment in ("v", "d", "j"):
        for suffix in ("sequence_start", "sequence_end", "alignment_start", "alignment_end"):
            key = f"{segment}_{suffix}"
            if record.get(key) is not None:
                displaced[key] = record[key] + 5
    found = coordinates.disagreements(displaced, lengths)
    assert found, "the whole rearrangement displaced by five bases went undetected"
    identities = {d.identity for d in found}
    assert "the first segment starts at read position 0" in identities, (
        f"caught, but not by the leading anchor; got {sorted(identities)}")


def test_an_unaccounted_gap_between_segments_is_caught():
    """Models np1 being recorded shorter than the untemplated stretch actually
    present -- read positions that belong to nothing."""
    record, lengths = _coordinate_case()
    assert record["np1_length"] > 0, "the fixture has no np1 to shorten"
    shortened = {**record, "np1_length": record["np1_length"] - 1}
    found = coordinates.disagreements(shortened, lengths)
    assert found, "an np1 length one short of the gap it describes went undetected"


def test_a_reference_too_short_for_the_coordinates_is_not_reported_clean():
    """The third silent route to a passing verdict.

    `reference is None or not 0 <= position < len(reference)` was a bare
    `continue`: not compared, not counted, not reported. A reference truncated
    by thirty bases gave 0 mismatches and `clean`, and the comparison count
    added alongside it did not close this -- nothing asserted on the count's
    absolute value, only that it was non-zero.

    Models a packaging change shipping a truncated or renumbered allele.
    """
    exp = experiment_builder.recombined()
    outcome = exp.run(n=1, seed=0)[0]
    simulation = outcome.final_simulation()
    control = reconstruction.templated_check(simulation, exp.refdata)
    assert control.clean and control.out_of_range == 0

    class _Truncated:
        """The reference, with the V allele short of what was assembled from it."""

        def __init__(self, real, missing):
            self._real, self._missing = real, missing

        def v_allele(self, identifier):
            allele = self._real.v_allele(identifier)
            sequence = allele.seq() if callable(allele.seq) else allele.seq
            return type("_Allele", (), {"seq": sequence[:-self._missing]})

        def __getattr__(self, name):
            return getattr(self._real, name)

    for missing in (3, 30):
        check = reconstruction.templated_check(simulation, _Truncated(exp.refdata, missing))
        assert check.mismatches == (), "a short reference has no base to disagree with"
        assert check.out_of_range == missing, (
            f"{missing} positions fall outside the truncated allele but "
            f"{check.out_of_range} were reported")
        assert not check.clean, (
            f"a reference {missing} bases short of the assembly reported clean")


# ── the projection oracle ────────────────────────────────────────────────────
#
# The cigars and alignment strings are what a benchmarked annotator is compared
# to, and the existing cigar tests read the reference zero times -- their
# strongest member re-derives the cigar from the production helpers and says
# "equivalence is guaranteed by construction". These mutants are the measure of
# whether crossing to the reference actually catches anything.


def _projection_case(*, invert=False):
    """``(record, allele_sequences, inverted)`` for one indel-bearing record."""
    exp = experiment_builder.recombined().polymerase_indels(count=3)
    if invert:
        # Compose the two: the first version built the indel experiment and
        # then discarded it, so the inverted mutant only ever ran on a
        # gap-free record and the reflection-with-gaps path had no mutant.
        exp = (experiment_builder.recombined()
               .invert_d(prob=1.0)
               .polymerase_indels(count=2))
    refdata = exp.refdata
    getters = {"v": refdata.v_allele, "d": refdata.d_allele, "j": refdata.j_allele}
    for outcome in exp.run(n=30, seed=53):
        record = dict(outcome.airr_record(refdata))
        if not invert and (record.get("n_indels") or 0) == 0:
            continue
        simulation = outcome.final_simulation()
        sequences = {}
        for segment, getter in getters.items():
            identifier = {"v": simulation.v_allele_id(), "d": simulation.d_allele_id(),
                          "j": simulation.j_allele_id()}[segment]
            if identifier is None:
                continue
            allele = getter(identifier)
            raw = allele.seq() if callable(allele.seq) else allele.seq
            sequences[segment] = raw.decode().upper()
        inverted = {"d"} if record.get("d_inverted") else set()
        return record, sequences, inverted
    pytest.fail("no suitable record was produced; the fixture is wrong")


def test_the_projection_oracle_passes_on_an_uncorrupted_record():
    """The control, with the coverage check beside it: a record whose columns
    were never compared would also report no disagreements."""
    record, sequences, inverted = _projection_case()
    assert projection.disagreements(record, sequences, inverted=inverted) == ()
    columns = projection.columns_compared(record, sequences, inverted=inverted)
    assert sum(columns.values()) > 200, (
        f"only {sum(columns.values())} alignment columns were compared; the "
        "mutants below would not mean much")


@pytest.mark.parametrize("name, corrupt, identity", [
    ("germline alignment truncated",
     lambda r: {**r, "germline_alignment": r["germline_alignment"][:-1]},
     "len(sequence_alignment) == len(germline_alignment)"),
    ("a read coordinate moved",
     lambda r: {**r, "v_sequence_end": r["v_sequence_end"] + 2},
     "cigar M+I == read span"),
    ("a germline coordinate moved",
     lambda r: {**r, "j_germline_end": r["j_germline_end"] + 2},
     "cigar M+D == germline span"),
    ("the alignment shifted against the read",
     lambda r: {**r, "sequence_alignment": "A" + r["sequence_alignment"][1:]},
     "sequence_alignment without gaps is the read"),
    ("the D mask left a base unmasked",
     lambda r: {**r, "germline_alignment_d_mask": r["germline_alignment"]},
     "the D mask is the germline alignment with the D span's bases replaced by N"),
])
def test_each_projection_identity_is_sensitive_to_its_own_producer(name, corrupt, identity):
    record, sequences, inverted = _projection_case()
    found = projection.disagreements(corrupt(record), sequences, inverted=inverted)
    assert found, f"{name} went undetected"
    assert any(identity in d.identity for d in found), (
        f"{name} was detected, but not by {identity!r}; "
        f"got {sorted({d.identity for d in found})}")


def test_a_cigar_with_the_right_op_lengths_in_the_wrong_place_is_caught():
    """The defect no existing cigar test can see.

    Every structural property survives moving a segment's alignment start: the
    ops are still M/I/D, the run-lengths still collapse, the op counts still
    sum to the same spans. Only the comparison against the reference notices
    that the M columns now land on the wrong germline bases.
    """
    record, sequences, inverted = _projection_case()
    assert projection.disagreements(record, sequences, inverted=inverted) == ()
    shifted = {**record, "v_alignment_start": record["v_alignment_start"] + 3}
    found = projection.disagreements(shifted, sequences, inverted=inverted)
    assert found, "a V alignment shifted three columns went undetected"
    assert any("carries the reference base" in d.identity for d in found), (
        f"caught, but not by the reference comparison; "
        f"got {sorted({d.identity for d in found})}")


def test_an_assembly_from_the_wrong_allele_is_caught_in_the_projection_too():
    """Models the record naming one allele while the alignment was built from
    another -- the projection equivalent of the IR-level allele-identity check."""
    record, sequences, inverted = _projection_case()
    wrong = dict(sequences)
    reference = wrong["v"]
    # Rotate the reference so it is a real sequence of the right length that is
    # not the one the alignment was built from.
    wrong["v"] = reference[7:] + reference[:7]
    found = projection.disagreements(record, wrong, inverted=inverted)
    assert found, "a rotated V reference went undetected"
    assert {d.segment for d in found} == {"v"}


def test_an_inverted_d_compared_as_if_forward_is_caught():
    """The orientation has to be told to the oracle, so the oracle must fail
    when it is told the wrong thing -- otherwise passing the right thing proves
    nothing.

    This is also the shape of the defect found in the oracle next door: its
    `skip_inverted` read an attribute a Simulation does not have, so it skipped
    nothing and claimed to skip inverted D by default.
    """
    record, sequences, inverted = _projection_case(invert=True)
    assert inverted == {"d"}, "the fixture did not invert D"
    assert any("D" in (record.get(f"{s}_cigar") or "") or "I" in (record.get(f"{s}_cigar") or "")
               for s in "vdj"), (
        "the inverted fixture carries no gap, so the reflection-with-gaps path "
        "is still untested")
    assert projection.disagreements(record, sequences, inverted=inverted) == ()
    found = projection.disagreements(record, sequences, inverted=())
    assert found, "an inverted D compared against the forward reference agreed"
    assert {d.segment for d in found} == {"d"}


def test_an_inverted_d_is_compared_by_the_base_level_oracle_rather_than_skipped():
    """`skip_inverted` was dead: it read `getattr(simulation, "d_inverted")`,
    a Simulation has no such attribute, so it always found None. The docstring
    said inverted D was skipped by default; nothing was skipped, and no
    registered test ran an inverting stack to notice.

    The IR stores the complement of the reference base at the original allele
    coordinate, so the comparison applies with the expectation complemented.
    """
    exp = experiment_builder.recombined().invert_d(prob=1.0)
    outcome = exp.run(n=1, seed=43)[0]
    simulation = outcome.final_simulation()
    record = outcome.airr_record(exp.refdata)
    assert record["d_inverted"] is True, "the fixture did not invert D"

    told = reconstruction.templated_check(simulation, exp.refdata, inverted={"D"})
    assert told.clean, f"an inverted D disagreed when declared: {told.mismatches[:2]}"
    assert told.inverted_compared > 0, (
        "no inverted position was compared, so this passed by skipping")

    not_told = reconstruction.templated_check(simulation, exp.refdata)
    assert not_told.mismatches, (
        "an inverted D compared against the forward reference agreed, which "
        "means the comparison is not reading the reference")
    assert {m.segment for m in not_told.mismatches} == {"D"}
