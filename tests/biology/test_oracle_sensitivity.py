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
from tests.support.oracles import coordinates, reconstruction

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


def test_an_indel_count_cannot_excuse_an_arbitrary_coordinate_shift():
    """The indel tolerance is per recorded edit, so it must not become a blanket
    licence: a shift larger than the recorded indel count is still a defect."""
    record, lengths = _coordinate_case()
    excused = {**record, "n_v_indels": 2, "n_indels": 2,
               "v_germline_end": record["v_germline_end"] - 2}
    assert coordinates.disagreements(excused, lengths) == (), (
        "two recorded indels should tolerate a two-base displacement")
    unexcused = {**excused, "v_germline_end": record["v_germline_end"] - 3}
    assert coordinates.disagreements(unexcused, lengths), (
        "a three-base displacement with two recorded indels went undetected")
