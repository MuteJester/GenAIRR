"""The per-segment identity fields, against what the read actually shows.

`v_identity` is the fraction of the V alignment that matches the germline. Any
benchmark that reports "percent identity to germline" -- the standard summary
of how mutated a repertoire is -- reads these three fields, and nothing
registered checked any of them: the mechanism was reported no-truth-exported
because two thirds of its fields are structurally null.

Those two thirds, `*_score` and `*_support`, are an aligner's inference
outputs. GenAIRR assembles from a known allele, so there is no hypothesis to
score and no choice to be confident about; the catalog records them as
unclaimable rather than as missing truth, and the columns stay in the output as
null because the AIRR schema has them. What is left is identity, which the
engine does compute, and which this module establishes.
"""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import reconstruction

pytestmark = pytest.mark.behavior

#: Pass stacks where an alignment column and a read position are the same
#: thing, so a fractional identity is matches over templated positions. The
#: indel stack is deliberately absent and has its own test below.
STACKS = {
    "plain": lambda exp: exp,
    "mutated": lambda exp: exp.mutate(model="s5f", count=10),
    "reverse-complemented": lambda exp: exp.random_strand_orientation(prob=1.0),
    "inverted-d": lambda exp: exp.invert_d(prob=1.0),
}

SEGMENT_FIELDS = (("V", "v_identity"), ("D", "d_identity"), ("J", "j_identity"))


def _recorded(record):
    """The three fields, named literally so the static check can see them."""
    return {"V": record["v_identity"], "D": record["d_identity"],
            "J": record["j_identity"]}


@verifies("projection.alignment_scores",
          claim="each segment's identity is the fraction of its templated "
                "positions whose base still matches the germline base "
                "recorded there, counted off the engine's own per-base "
                "provenance. Exact across mutation, reverse complement and an "
                "inverted D; a segment the record has no bases for reports "
                "null rather than 1.0, which would read as a perfect match",
          fields=("v_identity", "d_identity", "j_identity"),
          basis="provenance")
@pytest.mark.parametrize("stack", sorted(STACKS))
def test_each_segment_identity_is_what_the_read_shows(stack):
    """Derived from the provenance rather than compared to another field.

    The engine divides matching alignment columns by all of that segment's
    columns. This counts matching templated positions over all of them, which
    is the same quantity wherever an alignment column is a read position --
    every stack here. Measured over 25 records per stack: exact on all of
    them, to the last bit.

    The mutated stack is the one that discriminates. On plain records every
    identity is 1.0, so an implementation returning the constant 1.0 would
    pass; under s5f at ten mutations, 53 of the 75 segment identities fall
    below 1.0, which the floor below requires.
    """
    exp = STACKS[stack](experiment_builder.recombined())
    checked = below_one = 0
    for index, outcome in enumerate(exp.run(n=25, seed=11)):
        record = outcome.airr_record(exp.refdata)
        counts = reconstruction.germline_match_counts(outcome.final_simulation())
        recorded = _recorded(record)
        for segment, field in SEGMENT_FIELDS:
            if segment not in counts:
                assert recorded[segment] is None, (
                    f"{stack}: record {index} reports {field}="
                    f"{recorded[segment]} for a segment with no templated "
                    "base in the read, where null is the only honest value")
                continue
            matching, total = counts[segment]
            assert recorded[segment] == pytest.approx(matching / total, abs=1e-12), (
                f"{stack}: record {index} reports {field}={recorded[segment]} "
                f"where the read matches the germline at {matching} of "
                f"{total} {segment} positions")
            below_one += recorded[segment] < 1.0
            checked += 1
    assert checked >= 50, f"{stack}: only {checked} segment identities were checked"
    if stack == "mutated":
        assert below_one >= 40, (
            f"only {below_one} identities fell below 1.0 under ten mutations, "
            "so a constant 1.0 would satisfy this test")


@verifies("projection.alignment_scores",
          claim="an indel lowers the identity it touches and never raises it: "
                "an inserted or deleted column has no matching germline base "
                "and the engine counts it in its own denominator, so the "
                "recorded value sits at or below the substitution-only "
                "fraction the provenance gives",
          fields=("v_identity", "d_identity", "j_identity"),
          basis="provenance")
def test_an_indel_only_lowers_the_identity_it_touches():
    """The one stack where an alignment column is not a read position.

    The engine's denominator is alignment columns, which include the I and D
    ones; the provenance count has no column for a base that is not there. So
    the two are not equal here, and asserting equality would be asserting the
    wrong thing. What holds is the inequality, and its direction is the claim:
    an indel is a mismatch to the germline, not a free pass.

    Measured over 120 segment identities under polymerase_indels(count=3):
    recorded lower on 51, equal on 69, higher on none.
    """
    exp = experiment_builder.recombined().polymerase_indels(count=3)
    lower = equal = 0
    for index, outcome in enumerate(exp.run(n=40, seed=11)):
        record = outcome.airr_record(exp.refdata)
        counts = reconstruction.germline_match_counts(outcome.final_simulation())
        recorded = _recorded(record)
        for segment, field in SEGMENT_FIELDS:
            if segment not in counts or recorded[segment] is None:
                continue
            matching, total = counts[segment]
            substitutions_only = matching / total
            assert recorded[segment] <= substitutions_only + 1e-12, (
                f"record {index} reports {field}={recorded[segment]} above the "
                f"{substitutions_only} the substitutions alone give, so an "
                "indel column raised the identity rather than lowering it")
            if recorded[segment] < substitutions_only - 1e-12:
                lower += 1
            else:
                equal += 1
    assert lower + equal >= 100, (
        f"only {lower + equal} segment identities were compared under indels")
    assert lower > 0, (
        "no record's identity was lowered by an indel, so the inequality held "
        "trivially and the fixture is not producing indels inside the segments")


def test_score_and_support_are_null_on_every_record():
    """The six fields the catalog calls unclaimable, pinned as such.

    Excluding a field from coverage raises its mechanism's ratio without
    anyone testing anything, so the exclusion has to stay true. If the engine
    ever starts populating a score, this fails and the catalog entry -- which
    says there is no hypothesis to score -- has to be revisited rather than
    quietly outlived.
    """
    exp = experiment_builder.recombined().mutate(model="s5f", count=10)
    checked = 0
    for index, outcome in enumerate(exp.run(n=20, seed=11)):
        record = outcome.airr_record(exp.refdata)
        populated = {field: record[field] for field in
                     ("v_score", "d_score", "j_score",
                      "v_support", "d_support", "j_support")
                     if record[field] is not None}
        assert populated == {}, (
            f"record {index} carries {populated}, but the catalog excludes "
            "those fields from coverage on the grounds that nothing could "
            "establish them")
        checked += 1
    assert checked == 20


def test_the_score_columns_are_still_exported():
    """Null is not the same as absent.

    The AIRR rearrangement schema has these columns, and a consumer reading by
    name should find them rather than raise a KeyError -- which is the reason
    the catalog entry calls them unclaimable rather than recommending their
    removal.
    """
    exp = experiment_builder.recombined()
    record = next(iter(exp.run(n=1, seed=11))).airr_record(exp.refdata)
    for field in ("v_score", "d_score", "j_score",
                  "v_support", "d_support", "j_support",
                  "v_identity", "d_identity", "j_identity"):
        assert field in record, f"{field} is not an exported column at all"
