"""What the junction's untemplated content is, and what it is made of."""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import projection, reconstruction

pytestmark = pytest.mark.behavior


@verifies("junction.n_addition",
          claim="np1 and np2 are exactly the bases carrying no germline coordinate, "
                "in order, and their recorded lengths match",
          fields=("np1", "np2", "np1_length", "np2_length"),
          basis="provenance")
def test_np_fields_are_the_bases_with_no_germline_coordinate():
    """A zero-length NP region contributes no run, so the comparison is against
    the non-empty ones. Records with np1_length == 0 do occur and are the reason
    this is not a positional comparison."""
    exp = experiment_builder.recombined()
    saw_empty = False
    for outcome in exp.run(n=50, seed=0):
        record = outcome.airr_record(exp.refdata)
        runs = reconstruction.untemplated_runs(outcome.final_simulation())
        saw_empty |= record["np1_length"] == 0 or record["np2_length"] == 0
        disagreement = reconstruction.np_fields_disagreement(record, runs)
        assert disagreement is None, disagreement
    assert saw_empty, "no record had a zero-length NP region; widen the sample"


@verifies("junction.p_nucleotide",
          claim="the V 3' P nucleotides are the reverse complement of the V's "
                "own post-trim tip, and their recorded length matches the bases "
                "present",
          fields=("p_v_3_length",),
          basis="provenance")
def test_p_nucleotides_are_palindromes_of_the_segment_tip():
    """P is templated -- by the segment's own end, reversed and complemented.

    The provenance oracle records it with no germline coordinate, which is why
    it cannot verify this on its own: the distinction between P and N survives
    only in letter case.
    """
    exp = experiment_builder.with_p_nucleotides()
    checked = 0
    for outcome in exp.run(n=25, seed=1):
        simulation = outcome.final_simulation()
        record = outcome.airr_record(exp.refdata)
        bases = simulation.bases().decode()
        v_end = record["v_sequence_end"]
        p_length = record["p_v_3_length"]
        if not p_length or v_end is None:
            continue
        disagreement = reconstruction.p_tip_disagreement(bases, v_end, p_length)
        assert disagreement is None, (
            f"the V 3' P run is {disagreement[0]!r} but its tip's reverse "
            f"complement is {disagreement[1]!r}")
        checked += 1
    assert checked, "no record carried a V 3' P nucleotide; the fixture is wrong"


@verifies("junction.p_nucleotide",
          claim="P and N are distinguishable in the output, so a consumer can tell "
                "templated palindromic bases from untemplated ones",
          fields=("p_v_3_length", "p_d_5_length", "p_d_3_length", "p_j_5_length"),
          basis="consistent")
def test_p_and_n_are_distinguishable_in_the_untemplated_runs():
    exp = experiment_builder.with_p_nucleotides()
    seen = set()
    for outcome in exp.run(n=25, seed=1):
        record = outcome.airr_record(exp.refdata)
        total_p = sum(record[f] for f in ("p_v_3_length", "p_d_5_length",
                                          "p_d_3_length", "p_j_5_length"))
        runs = reconstruction.untemplated_runs(outcome.final_simulation())
        lowercase = sum(sum(1 for c in r.bases if c.islower()) for r in runs)
        assert lowercase == total_p, (
            f"{lowercase} lowercase untemplated bases but {total_p} P nucleotides "
            "recorded; the P/N channel does not agree with the P lengths")
        seen.update(r.kind for r in runs)
    assert {"P", "N"} & seen, f"no run was classified P or N: {seen}"


@verifies("junction.n_addition",
          claim="np1 and np2 are the stretch of the molecule between the "
                "segment bounds -- np1 == sequence[v_sequence_end:"
                "d_sequence_start], and its length -- which is what the "
                "fields promise and what holds when a polymerase insertion "
                "lands in the junction. What it does NOT establish is the "
                "TdT-added amount, which the interval also contains P "
                "nucleotides and any such insertion on top of; see the "
                "mechanism's truth_gap",
          fields=("np1", "np2", "np1_length", "np2_length"),
          basis="consistent")
@pytest.mark.parametrize("with_indels", [False, True], ids=["plain", "indels"])
def test_np_fields_are_the_interval_between_the_segments(with_indels):
    """The np fields' documented definition, checked on both sides of an indel
    pass because that is where it diverges from the mechanism's own summary.

    The run-splitting oracle cannot be used here: an inserted base carries no
    germline coordinate, so it forms an untemplated run of its own and the
    oracle reports a disagreement on 194 of 200 records whose np fields are
    correct. This compares the fields against the sequence at the coordinates
    the record itself publishes, which is why the basis is `consistent` -- two
    exported views of one molecule. It is enough to catch an np field drifting
    from its own interval, and deliberately not claimed as more.
    """
    exp = experiment_builder.recombined()
    if with_indels:
        exp = exp.polymerase_indels(count=(1, 3), insertion_prob=1.0)
    checked = nonempty = 0
    for index, outcome in enumerate(exp.run(n=100, seed=5)):
        record = outcome.airr_record(exp.refdata)
        sequence = record["sequence"]
        for np_field, length_field, start, end in (
                ("np1", "np1_length", "v_sequence_end", "d_sequence_start"),
                ("np2", "np2_length", "d_sequence_end", "j_sequence_start")):
            interval = sequence[record[start]:record[end]]
            assert record[np_field].upper() == interval.upper(), (
                f"record {index}: {np_field} is {record[np_field]!r} but the "
                f"molecule between {start} and {end} is {interval!r}")
            assert record[length_field] == len(record[np_field]), (
                f"record {index}: {length_field} is {record[length_field]} for "
                f"a {len(record[np_field])}-base {np_field}")
            checked += 1
            nonempty += bool(interval)
    assert checked == 200, f"only {checked} intervals were compared"
    assert nonempty > 100, (
        f"only {nonempty} intervals carried any base; the identity is "
        f"satisfied trivially by empty strings")


P_BOUNDARIES = (
    ("p_v_3_length", "v_sequence_end", "3prime"),
    ("p_d_5_length", "d_sequence_start", "5prime"),
    ("p_d_3_length", "d_sequence_end", "3prime"),
    ("p_j_5_length", "j_sequence_start", "5prime"),
)


def _p_length_fields_present(record):
    """The four P lengths, named literally.

    A claimed field has to be readable in the body that claims it, and a name
    arriving through a parametrize table is as invisible to the static check
    as an f-string is.
    """
    for field in ("p_v_3_length", "p_d_5_length", "p_d_3_length",
                  "p_j_5_length"):
        value = record.get(field)
        assert value is None or isinstance(value, int), f"{field} is not an integer"


@verifies("junction.p_nucleotide",
          claim="each recorded P length names a run that IS the palindrome of "
                "its segment's own post-trim tip, at all four coding ends. "
                "The 3' ends read the tip behind the run and the 5' ends the "
                "tip ahead of it, which is the asymmetry the engine's own "
                "derivation has. Measured 60 of 60 runs at every boundary",
          fields=("p_v_3_length", "p_d_5_length", "p_d_3_length",
                  "p_j_5_length"),
          basis="provenance")
@pytest.mark.parametrize("length_field,boundary_field,side", P_BOUNDARIES,
                         ids=[f for f, _, _ in P_BOUNDARIES])
def test_every_p_length_names_a_tip_palindrome(length_field, boundary_field, side):
    """Three of these four fields were previously checked only against the
    letter case the P pass itself writes -- two views of one producer. This
    reads the bases.

    It needs a fixture the old one could not provide: with P configured only
    at the V 3' and J 5' ends, the two D lengths have no run at all over 60
    records, so a palindrome check on them would have passed by never
    comparing anything.
    """
    exp = experiment_builder.with_p_nucleotides_at_every_end()
    checked = 0
    for index, outcome in enumerate(exp.run(n=60, seed=1)):
        record = outcome.airr_record(exp.refdata)
        bases = outcome.final_simulation().bases().decode()
        _p_length_fields_present(record)
        length = record.get(length_field) or 0
        boundary = record.get(boundary_field)
        if length <= 0 or boundary is None:
            continue
        if side == "3prime":
            run, tip = bases[boundary:boundary + length], bases[boundary - length:boundary]
        else:
            run, tip = bases[boundary - length:boundary], bases[boundary:boundary + length]
        assert run.upper() == reconstruction.reverse_complement(tip).upper(), (
            f"record {index}: the {length}-base run at {boundary_field} is "
            f"{run!r}, not the reverse complement of the tip {tip!r}")
        checked += 1
    assert checked == 60, (
        f"{length_field}: only {checked} of 60 records carried a run, so the "
        f"fixture does not exercise this boundary")


NP_TRANSLATION_FIXTURES = (
    ("plain", lambda: experiment_builder.recombined()),
    ("mutated", lambda: experiment_builder.recombined().mutate(model="s5f", count=8)),
    ("reverse_complemented", lambda: experiment_builder.recombined()
     .random_strand_orientation(prob=1.0)),
    ("light_chain", lambda: experiment_builder.recombined(locus="human_igk")),
)


def _np_translation_fields_present(record):
    """Both NP translations, named literally."""
    for field in ("np1_aa", "np2_aa"):
        value = record.get(field)
        assert value is None or isinstance(value, str), f"{field} is not a string"


@verifies("junction.n_addition",
          claim="np1_aa and np2_aa are the whole codons inside the NP "
                "interval, translated in the record's own reading frame with "
                "an independent codon table. A partial codon at either end is "
                "dropped rather than padded. Measured 60 of 60 on plain, "
                "hypermutated, reverse-complemented and light-chain records",
          fields=("np1_aa", "np2_aa"),
          basis="reference")
@pytest.mark.parametrize("label,make", NP_TRANSLATION_FIXTURES,
                         ids=[name for name, _ in NP_TRANSLATION_FIXTURES])
def test_the_np_translations_are_the_whole_codons_of_the_interval(label, make):
    """The last two fields of this mechanism that nothing read.

    The frame comes from the translation oracle rather than being
    reconstructed here, because getting it from `v_trim_5 + end_loss_5_length`
    is wrong whenever two 5' truncations compose -- a mistake that oracle
    records having made, with 30 of 30 false failures.

    The codon table is this package's own, so the expected residues do not
    come from the engine. What makes the interval independent is the
    coordinate pair: the NP span is bounded by the segment coordinates, which
    are themselves anchored to the span their bases occupy.
    """
    exp = make()
    checked = 0
    for index, outcome in enumerate(exp.run(n=60, seed=5)):
        record = outcome.airr_record(exp.refdata)
        _np_translation_fields_present(record)
        phase = projection.translation_phase(record, native_phase=0)
        subject = (reconstruction.reverse_complement(record["sequence"])
                   if record.get("rev_comp") else record["sequence"])
        # On a VJ chain there is no D, so np1 spans V to J directly and np2
        # does not exist. Keying the spans off d_sequence_* unconditionally
        # compared nothing at all on IGK -- 0 translations, which the
        # `checked` floor below is what caught.
        if record.get("d_sequence_start") is None:
            assert not record.get("np2_aa"), (
                f"{label}: record {index} has no D yet carries np2_aa "
                f"{record.get('np2_aa')!r}")
            spans = (("np1_aa", "v_sequence_end", "j_sequence_start"),)
        else:
            spans = (("np1_aa", "v_sequence_end", "d_sequence_start"),
                     ("np2_aa", "d_sequence_end", "j_sequence_start"))
        for field, start_field, end_field in spans:
            span_start = record.get(start_field)
            span_end = record.get(end_field)
            recorded = record.get(field)
            if span_start is None or span_end is None or recorded is None:
                continue
            first = span_start + ((phase - span_start) % 3)
            last = span_end - ((span_end - first) % 3)
            expected = projection.translate(subject[first:last]) if last > first else ""
            assert recorded.upper() == expected.upper(), (
                f"{label}: record {index} {field} is {recorded!r}, not the "
                f"whole codons of {span_start}..{span_end} in frame {phase}")
            checked += 1
    # A floor rather than a count, because a VJ chain contributes one
    # translation per record and a VDJ chain two. It exists because keying the
    # spans off the D coordinates silently compared nothing on IGK and the
    # test passed.
    assert checked >= 55, f"{label}: only {checked} translations were compared"
