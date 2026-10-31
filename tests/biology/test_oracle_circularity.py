"""The circularity probes, as code rather than as prose.

Every reporting oracle carries a ``CIRCULARITY`` mapping: per identity, what a
both-sides perturbation produces, measured. That mapping is the honest part of
this suite's claim -- it is where an identity that cannot see a defect says so
-- and until this module existed, every number in it was a sentence. The
measurements were real when they were taken; nothing re-took them, so a
rewritten oracle would keep its old notes and nobody would know.

A probe is not a mutant. A mutant corrupts one side and asserts the oracle
fires, which shows the check works. A probe corrupts BOTH sides consistently
and records what the oracle says, which shows what the check cannot see. Some
probes here therefore assert SILENCE, and each of those says in its docstring
why silence is the identity holding rather than a hole.

``PROSE_ONLY`` below names the identities whose probe is still a sentence, so
the gap is a number that can only go down.
"""
from __future__ import annotations

import importlib

import pytest

from tests.support import alleles
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import expression as expression_oracle
from tests.support.oracles import junction as junction_oracle
from tests.support.oracles import reconstruction

pytestmark = pytest.mark.behavior

#: ``(oracle, identity) -> [test names]``, filled at import by the decorator.
_PROBED: dict[tuple[str, str], list[str]] = {}


def probes(oracle: str, identity: str):
    """Mark a test as the executable probe for one ``CIRCULARITY`` identity.

    Validated at decoration time against that oracle's mapping, so a probe
    cannot claim an identity the oracle does not record -- the same rule the
    oracles' own ``IDENTITIES`` tuples enforce at construction, for the same
    reason.
    """
    module = importlib.import_module(f"tests.support.oracles.{oracle}")
    if identity not in module.CIRCULARITY:
        raise ValueError(
            f"{oracle} records no circularity note for {identity!r}; the "
            f"notes it has are {sorted(module.CIRCULARITY)}")

    def decorate(function):
        _PROBED.setdefault((oracle, identity), []).append(function.__name__)
        function.genairr_probes = (oracle, identity)
        return function
    return decorate


#: Identities whose probe is still only a sentence in the oracle's mapping,
#: with what it would take to execute it. Pinned so adding an identity without
#: a probe is a decision somebody records rather than a silent omission, and so
#: the list can only shrink.
#:
#: The three oracles whose strongest identities already have record-level
#: mutants in test_oracle_sensitivity.py are the reason this list is as long as
#: it is: a mutant is weaker than a probe, but it is not nothing, and the
#: identities with neither were the ones worth writing first.
PROSE_ONLY = {
    ("coordinates", "germline_start == trim_5 + end_loss_5_length"):
        "the five coordinate identities are probed by the parametrized "
        "producer mutants in test_oracle_sensitivity.py, which corrupt one "
        "side; the both-sides perturbation is the over-determined set the "
        "module's docstring describes and needs a second reference table",
    ("coordinates", "germline_end == allele_length - trim_3 - end_loss_3_length"):
        "as above",
    ("coordinates", "germline span == read span"): "as above",
    ("coordinates", "alignment coordinates == read coordinates"): "as above",
    ("coordinates", "alignment coordinates are at or after the read coordinates"):
        "as above",
    ("coordinates", "the first segment starts at read position 0"):
        "the three chain identities have displacement mutants in "
        "test_oracle_sensitivity.py",
    ("coordinates", "the gap between segments is the untemplated length recorded there"):
        "as above",
    ("coordinates", "the read ends where the last segment ends"): "as above",
    ("junction", "junction_start == first surviving base of the V anchor codon"):
        "an anchor-moved mutant exists in test_oracle_sensitivity.py; the "
        "both-sides probe needs a simulation whose provenance moves with the "
        "bound, which the record-level mutants cannot build",
    ("junction", "junction_end == one past the last surviving base of the J anchor codon"):
        "as above",
    ("junction", "junction == the subject sliced at the recorded bounds"):
        "recorded CIRCULAR: the bounds move with the slice. Executing it "
        "needs a consistent re-slice of both bounds and the string, which is "
        "the mutant test_a_junction_bound_moved_off_the_anchor_is_reported "
        "builds one side of",
    ("junction", "junction_length == junction_end - junction_start"): "as above",
    ("junction", "junction_aa == translation of junction"): "as above",
    ("junction", "a record with no junction reports it empty and says the frame is unresolved"):
        "the null-group contract is asserted directly by "
        "test_an_unresolvable_frame_is_null_rather_than_false; there is no "
        "second side to perturb",
    ("projection", "the alignment strings are present to be checked"):
        "the eleven projection identities are probed one-sidedly by the "
        "parametrized mutants in test_oracle_sensitivity.py",
    ("projection", "a segment the record names is checkable"): "as above",
    ("projection", "len(sequence_alignment) == len(germline_alignment)"): "as above",
    ("projection", "len(sequence_alignment) == len(sequence) + deletions"): "as above",
    ("projection", "the cigar uses only M, I and D"): "as above",
    ("projection", "every cigar column is inside the alignment"): "as above",
    ("projection", "a column's germline coordinate is inside the allele"): "as above",
    ("projection", "a column consuming a germline base carries the reference's base there"):
        "as above",
    ("projection", "a column declared M/D/I has that op's gap pattern"): "as above",
    ("projection", "cigar M+I == read span"): "as above",
    ("projection", "cigar M+D == germline span"): "as above",
    ("projection", "sequence_alignment without gaps is the read"): "as above",
    ("projection", "the D mask is the germline alignment with the D span's bases replaced by N"):
        "now checked against the provenance D span by "
        "test_the_alignment_and_the_d_mask_against_the_recorded_molecule, "
        "which is the independent side the note says this identity lacks",
    ("reconstruction", "recorded germline base == reference base at the claimed coordinate"):
        "recorded CIRCULAR by a coordinate shift over both the byte and the "
        "coordinate, which needs a corrupted Simulation view; the one-sided "
        "shift is test_a_germline_coordinate_shifted_by_one_is_caught",
    ("reconstruction", "observable_substitutions == the engine's n_observable_mutations"):
        "as above, plus a no-mutation pipeline with sequencing errors",
    ("reconstruction", "recorded germline base == base of the allele the CALL names"):
        "probed by test_a_call_naming_the_wrong_allele_is_reported, which "
        "perturbs the call rather than both sides",
    ("reconstruction", "an inverted segment's byte == complement of the reference base"):
        "probed by test_an_inverted_d_compared_as_if_forward_is_caught",
    ("reconstruction", "the np fields are the untemplated runs' bases and lengths"):
        "the both-sides rewrite needs a Simulation view whose untemplated "
        "bases are replaced; the one-sided probe is "
        "test_a_swapped_np_pair_is_caught",
    ("reconstruction", "a P run is the reverse complement of the segment's post-trim tip"):
        "the matched-pair overwrite is a string perturbation and could be "
        "executed; the one-sided probe is "
        "test_a_p_nucleotide_that_is_not_a_palindrome_is_caught",
}


def _mutated_records(n=60, seed=8):
    """Records from the fixture the expression notes were measured on.

    s5f at ten mutations: the fixture where anchors actually get destroyed and
    stops actually appear, which is why the notes were taken here rather than
    on plain records.
    """
    exp = experiment_builder.recombined().mutate(model="s5f", count=10)
    out = [outcome.airr_record(exp.refdata) for outcome in exp.run(n=n, seed=seed)]
    return [record for record in out if record["junction_length"] is not None]


def _fires(records, identity, *, mutate):
    """How many of ``records`` report ``identity`` once ``mutate`` is applied."""
    count = 0
    for record in records:
        copy = dict(record)
        mutate(copy)
        if any(d.identity == identity for d in expression_oracle.disagreements(copy)):
            count += 1
    return count


IN_FRAME = "vj_in_frame == the junction length being a whole number of codons"
STOP = "stop_codon == a stop in the translated coding span"
PRODUCTIVE = "productive == in frame, no stop, and both anchors intact"


# ── expression ───────────────────────────────────────────────────────────────

@probes("expression", IN_FRAME)
def test_flipping_the_frame_flag_alone_is_caught():
    """One side moved: the flag disagrees with the length it is derived from.

    The note records 199 of 199 on this fixture. Asserted as every record
    rather than as that count, because the count is a property of the sample
    size and the thing being established is that there is no record where a
    flipped flag goes unnoticed.
    """
    records = _mutated_records()
    assert len(records) > 40, f"only {len(records)} records kept a junction"

    def flip(record):
        record["vj_in_frame"] = not record["vj_in_frame"]

    assert _fires(records, IN_FRAME, mutate=flip) == len(records)
    assert _fires(records, IN_FRAME, mutate=lambda record: None) == 0, (
        "the unperturbed records already disagree, so the probe above "
        "establishes nothing")


@probes("expression", IN_FRAME)
def test_moving_the_length_with_the_flag_is_invisible_and_the_conjunction_is_not():
    """Both sides moved: this is the blind spot the note records.

    ``vj_in_frame`` and ``junction_length`` are two fields of one record
    written by one producer, so a consistent lie satisfies the identity. What
    catches it is the productive conjunction, which the note measured at 54 of
    199 -- the records where the flipped frame contradicts the recorded
    assessment. The junction oracle's bound identities are what establish the
    length underneath this one; see its own probe list.
    """
    records = _mutated_records()

    def flip_both(record):
        length = record["junction_length"]
        record["vj_in_frame"] = not record["vj_in_frame"]
        if record["vj_in_frame"]:
            # Must become a whole number of codons. Rounding up is what keeps
            # this right for both remainders; subtracting one only worked when
            # the length happened to be 1 mod 3, which is why the first
            # version of this probe reported 21 false catches.
            record["junction_length"] = length + (3 - length % 3) % 3
        else:
            record["junction_length"] = length + 1

    assert _fires(records, IN_FRAME, mutate=flip_both) == 0, (
        "the frame identity caught a consistent perturbation of both its "
        "sides, which would make the note wrong in the oracle's favour")
    caught = _fires(records, PRODUCTIVE, mutate=flip_both)
    assert caught > 0, (
        "no identity noticed a record whose frame flag and junction length "
        "were both rewritten, so that pair is unguarded")


@probes("expression", STOP)
def test_flipping_the_stop_flag_alone_is_caught():
    """One side moved. The note records 199 of 199."""
    records = _mutated_records()

    def flip(record):
        record["stop_codon"] = not record["stop_codon"]

    assert _fires(records, STOP, mutate=flip) == len(records)


@probes("expression", STOP)
def test_a_stop_written_into_the_molecule_with_the_flag_is_the_identity_holding():
    """Both sides moved -- and here silence is the identity working.

    The molecule is written by the assembly passes and the flag by the
    productivity assessment, so a stop codon placed in the sequence WITH the
    flag set is not a consistent lie: it is what a real stop looks like. The
    identity staying silent is it agreeing with a true state of affairs, which
    is why the note calls this one NOT CIRCULAR despite the silence.

    The productive conjunction fires on the same perturbation -- the note
    measured 27 of 199 -- because a record carrying a stop cannot express.
    """
    records = _mutated_records()

    def write_a_stop(record):
        # In the coding span and in frame: start from v_sequence_start, which
        # the oracle advances to the record's own reading frame.
        start = record["v_sequence_start"]
        offset = start + 3 * ((record["junction_start"] - start) // 3)
        sequence = record["sequence"]
        record["sequence"] = sequence[:offset] + "TAA" + sequence[offset + 3:]
        record["stop_codon"] = True

    assert _fires(records, STOP, mutate=write_a_stop) == 0, (
        "the stop identity reported a record whose sequence really does carry "
        "a stop in its coding span, which is the flag being right")
    assert _fires(records, PRODUCTIVE, mutate=write_a_stop) > 0, (
        "no record changed its productive assessment when a stop codon was "
        "written into its coding span")


@probes("expression", PRODUCTIVE)
def test_flipping_the_productive_flag_alone_is_caught():
    """One side moved. The note records 199 of 199, and this is the identity
    that fires on both joint probes above."""
    records = _mutated_records()

    def flip(record):
        record["productive"] = not record["productive"]

    assert _fires(records, PRODUCTIVE, mutate=flip) == len(records)


# ── junction: the reference-only identity ────────────────────────────────────

RESIDUE = "the reference anchor codon carries the conserved residue"


@probes("junction", RESIDUE)
def test_a_reference_anchor_codon_that_is_not_conserved_is_caught():
    """The identity with no side from the run to move, and no mutant until now.

    It is what makes the two bound identities mean something -- a bound at the
    first surviving base of a codon that is not a cysteine is a bound into the
    wrong place -- and it is the one that would catch a mis-curated anchor in a
    reference cartridge. The note sweeps the bundled loci and finds every
    anchor conserved, which is a measurement of the DATA; this is a measurement
    of the CHECK.
    """
    exp = experiment_builder.recombined()
    outcome = next(iter(exp.run(n=1, seed=5)))
    record = outcome.airr_record(exp.refdata)
    simulation = outcome.final_simulation()
    segment_of = reconstruction.segments_by_position(simulation)
    v_allele = alleles.by_name(exp.refdata, "V", record["v_call"])
    j_allele = alleles.by_name(exp.refdata, "J", record["j_call"])
    v_sequence = alleles.sequence(v_allele)
    j_sequence = alleles.sequence(j_allele)

    def run(v_bytes=v_sequence, j_bytes=j_sequence):
        return junction_oracle.disagreements(
            record, simulation, segment_of,
            v_anchor=v_allele.anchor, j_anchor=j_allele.anchor,
            v_sequence=v_bytes, j_sequence=j_bytes)

    assert run() == (), "the unperturbed record already disagrees"

    for label, anchor, sequence in (("V", v_allele.anchor, v_sequence),
                                    ("J", j_allele.anchor, j_sequence)):
        # Overwrite the anchor codon with one that is neither a cysteine nor a
        # tryptophan/phenylalanine. GCA is alanine, which is neither.
        broken = bytearray(sequence)
        broken[anchor:anchor + 3] = b"GCA"
        found = run(**{f"{label.lower()}_bytes": bytes(broken)})
        assert [d for d in found if d.identity == RESIDUE], (
            f"a {label} anchor codon rewritten to alanine went unreported, so "
            "a mis-curated anchor in a reference cartridge would too")

    # And the both-sides half: a DIFFERENT conserved codon is accepted, which
    # is the identity declining to care which of the synonymous codons the
    # reference uses. TGT and TGC are both cysteine; TTT and TTC both
    # phenylalanine.
    for label, anchor, sequence, codon in (("V", v_allele.anchor, v_sequence, b"TGC"),
                                           ("J", j_allele.anchor, j_sequence, b"TTC")):
        swapped = bytearray(sequence)
        swapped[anchor:anchor + 3] = codon
        found = run(**{f"{label.lower()}_bytes": bytes(swapped)})
        assert not [d for d in found if d.identity == RESIDUE], (
            f"a {label} anchor rewritten to {codon.decode()} was reported, "
            "but that codon carries the conserved residue")


def _region_counters(record):
    """The six per-subregion observable counters, named literally.

    The same six the registered truth test reads, spelled out here rather than
    imported from it: a test module may not import another test module, and an
    f-string loop would hide the field names from the static check.
    """
    return {
        "fwr1": record["n_fwr1_observable_mutations"],
        "cdr1": record["n_cdr1_observable_mutations"],
        "fwr2": record["n_fwr2_observable_mutations"],
        "cdr2": record["n_cdr2_observable_mutations"],
        "fwr3": record["n_fwr3_observable_mutations"],
        "v_unannotated": record["n_v_unannotated_observable_mutations"],
    }


# ── reconstruction: the bucketing helpers ────────────────────────────────────

BUCKETED = "a visible substitution is bucketed by the region or segment carrying it"


def _mutated_outcomes(*, count, n, seed):
    """``(experiment, outcomes)`` for a hypermutated fixture.

    The bucket comparison is an equality at one mutation per record -- an
    event and a visible difference are the same thing there -- and stays an
    equality against the OBSERVABLE counters at twenty, which is what these
    probes compare against.
    """
    exp = experiment_builder.recombined().mutate(model="s5f", count=count)
    return exp, list(exp.run(n=n, seed=seed))


@probes("reconstruction", BUCKETED)
def test_a_displaced_subregion_table_no_longer_matches_the_recorded_counters():
    """The subregion half: the bucketing is what the base comparison misses.

    The comparison is the registered truth test's -- the oracle's buckets
    against the record's own six counters -- run with the table perturbed. A
    table shifted by 15 and an emptied table both stop matching, so the
    reference's IMGT intervals are load-bearing in the result rather than
    decorative. Reversing the interval list changes nothing at all, because
    the intervals are disjoint and the same one still matches whichever end
    the search starts from; that measurement is recorded because it looks like
    a blind spot and is not one.

    Executing the note corrected it. It reported 58 of 60 for the shift and 60
    of 60 for the emptied table under a sentence that said "probed at one
    mutation per record"; at one mutation the shift stops matching on 11 of
    the 49 records whose mutation lands in V, and the emptied table on 48 of
    49. The 58 and the 60 are the twenty-mutation fixture, where every record
    carries V substitutions spread widely enough that a 15-base displacement
    moves at least one across a boundary. This probe runs that fixture, and
    the oracle's note now says which one the numbers came from.
    """
    exp, outcomes = _mutated_outcomes(count=20, n=60, seed=7)
    shifted = emptied = reversed_same = compared = 0
    for outcome in outcomes:
        simulation = outcome.final_simulation()
        record = outcome.airr_record(exp.refdata)
        table = list(alleles.v_subregions(exp.refdata, record["v_call"]))
        reported = _region_counters(record)
        truth = reconstruction.observable_substitutions_by_v_subregion(
            simulation, table)
        assert reported == truth, (
            f"the unperturbed record reports {reported} where the read shows "
            f"{truth}, so a perturbation of the table proves nothing")
        if sum(truth.values()) == 0:
            continue                     # the mutation did not land in V
        compared += 1
        moved = [(label, start + 15, end + 15) for label, start, end in table]
        shifted += reported != reconstruction.observable_substitutions_by_v_subregion(
            simulation, moved)
        emptied += reported != reconstruction.observable_substitutions_by_v_subregion(
            simulation, [])
        reversed_same += reported == reconstruction.observable_substitutions_by_v_subregion(
            simulation, list(reversed(table)))
    assert compared >= 20, (
        f"only {compared} records had a visible V substitution to bucket")
    assert shifted >= compared * 0.9, (
        f"shifting every interval by 15 stopped matching on only {shifted} of "
        f"{compared} records, so most of the table's positions are not "
        "distinguished by it")
    assert emptied == compared, (
        f"emptying the table still matched on {compared - emptied} records, "
        "so the table is not what decides the buckets")
    assert reversed_same == compared, (
        f"reversing the interval list changed the answer on "
        f"{compared - reversed_same} records; the intervals are disjoint, so "
        "the order must not matter")


@probes("reconstruction", BUCKETED)
def test_swapping_the_recorded_segment_counts_is_caught():
    """The segment half, and why the partition identity is not enough.

    Moving one mutation from D to J leaves the sum unchanged, so a partition
    check cannot see a misrouted event. Swapping the recorded V and D counts
    against the bucketed truth disagrees on most records -- the note measured
    55 of 60, the survivors being records whose two counts happen to be equal.
    """
    exp, outcomes = _mutated_outcomes(count=1, n=40, seed=4)
    caught = equal_counts = compared = 0
    for outcome in outcomes:
        simulation = outcome.final_simulation()
        record = outcome.airr_record(exp.refdata)
        buckets = reconstruction.observable_substitutions_by_segment(simulation)
        recorded = {segment: record[f"n_{segment}_mutations"]
                    for segment in ("v", "d", "j")
                    if record.get(f"n_{segment}_mutations") is not None}
        if set(recorded) < {"v", "d"}:
            continue
        compared += 1
        assert all(buckets.get(segment, 0) == count
                   for segment, count in recorded.items()), (
            f"the unperturbed buckets {buckets} disagree with the recorded "
            f"{recorded}, so the swap below would prove nothing")
        swapped = dict(recorded, v=recorded["d"], d=recorded["v"])
        if recorded["v"] == recorded["d"]:
            equal_counts += 1
            continue
        caught += any(buckets.get(segment, 0) != count
                      for segment, count in swapped.items())
    assert compared >= 20, f"only {compared} records carried V and D counts"
    assert caught == compared - equal_counts, (
        f"swapping V and D was caught on {caught} of "
        f"{compared - equal_counts} records where the two counts differ")


# ── the coverage pin ─────────────────────────────────────────────────────────

def test_every_circularity_identity_is_either_probed_or_pinned_as_prose():
    """The gap, as a number that can only go down.

    Before this module every one of the 36 identities was prose. An identity
    that is neither probed here nor named in PROSE_ONLY fails this test, so
    adding one is a decision somebody records.
    """
    from tests.support import oracles

    declared = {(name, identity)
                for name, module in oracles.reporting_modules().items()
                for identity in getattr(module, "CIRCULARITY", {})}
    probed = set(_PROBED)
    assert probed <= declared, (
        f"probes naming identities no oracle records: {sorted(probed - declared)}")
    unaccounted = sorted(declared - probed - set(PROSE_ONLY))
    assert unaccounted == [], (
        f"{len(unaccounted)} circularity identities are neither probed here "
        f"nor pinned as prose: {unaccounted}")
    stale = sorted(set(PROSE_ONLY) & probed)
    assert stale == [], (
        f"these identities are probed and still pinned as prose: {stale}")
    assert len(probed) >= 5, (
        f"only {len(probed)} identities have an executable probe; the "
        "decorator is not registering them")
