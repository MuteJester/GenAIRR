"""Can a consumer reproduce the productivity flags from the record?

This is the question a benchmarking oracle has to be able to answer yes to.
GenAIRR exists so other people can score their tools against its ground truth,
and a flag nobody can recompute is a flag they have to take on trust -- so
when their aligner disagrees with ours they cannot tell which of the two is
wrong.

Before the fields this module reads existed, the answer was no on three
counts, all found by the pipeline-space sweep:

- `stop_codon` is "a stop in the coding span", and the span was not exported.
  A reconstruction from `v_sequence_start`, the frame and `j_sequence_end` is
  not the same span: with ambiguous bases, 32 of 325 IGH records had the N
  inside the engine's span and outside the reconstruction, or the reverse.
- `vj_in_frame` is not `junction_length % 3 == 0`. It also requires the
  junction to begin on a codon boundary of the coding frame, which is the
  clause an indel breaks -- 54 of 200 records disagreed with the naive
  computation -- and the frame's phase was not exported.
- `productive` is a five-predicate conjunction under three-valued logic, and
  three of the five were not exported: the J tail's frame and the two anchor
  predicates, whose tri-state depends on provenance a reader does not have
  (an ambiguous base in the anchor codon makes it unresolved rather than
  rejected; an insertion inside the codon disrupts it without changing the
  base a reader sees).

All three are reproducible now, and this module is what says so. It is
deliberately written as a consumer would write it -- from the exported fields
and an independent genetic code, with no access to the simulation -- because
that is the claim.
"""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import expression as expression_oracle

pytestmark = pytest.mark.behavior

#: Loci and pass stacks the claim is measured over. The stacks matter more
#: than the loci: each one broke a different part of the reproduction before
#: the fields existed.
LOCI = ("human_igh", "human_igk", "human_tcrb")
STACKS = {
    "plain": lambda experiment: experiment,
    "deletions": lambda experiment: experiment.polymerase_indels(
        count=1, insertion_prob=0.0),
    "insertions": lambda experiment: experiment.polymerase_indels(
        count=1, insertion_prob=1.0),
    "mutated": lambda experiment: experiment.mutate(model="s5f", count=10),
    "ambiguous": lambda experiment: experiment.ambiguous_base_calls(count=2),
    "trimmed": lambda experiment: experiment.trim(
        v_3=[(length, 1.0) for length in range(0, 9)]),
}
#: Somatic hypermutation is a B-cell phenomenon and the engine refuses it on a
#: TCR locus.
B_CELL_ONLY = {"mutated"}

RECORDS = 120
SEED = 0


def _in_frame(record) -> bool:
    """`vj_in_frame` from the exported fields.

    Both clauses: the junction is a whole number of codons AND it begins on a
    codon boundary of the coding frame. The frame's phase is
    `coding_frame_start % 3` -- the span's start is by construction the first
    position in the frame, so its residue mod three IS the phase.
    """
    return ((record["junction_end"] - record["junction_start"]) % 3 == 0
            and record["junction_start"] % 3 == record["coding_frame_start"] % 3)


def _productive(record):
    """`productive` as the five-predicate conjunction, three-valued.

    False if any predicate is false, True if all are true, null otherwise --
    the same rule the engine applies, from the five values the record now
    carries.
    """
    stop = expression_oracle.coding_stop_state(record)
    predicates = [record["v_anchor_accepted"], record["j_anchor_accepted"],
                  _in_frame(record),
                  None if stop is None else not stop,
                  record["j_tail_in_frame"]]
    if False in predicates:
        return False
    return True if all(one is True for one in predicates) else None


@verifies("expression.productivity",
          claim="all three productivity flags are reproducible from the "
                "exported record by a consumer with no access to the "
                "simulation: stop_codon by translating the recorded coding "
                "span with an independent codon table three-valued over "
                "ambiguity, vj_in_frame from the junction bounds and the "
                "frame's phase, and productive as the five-predicate "
                "conjunction the record now carries in full. Measured exactly "
                "on 2040 records across three loci and six pass stacks, "
                "including the indel and ambiguous-base cases where each of "
                "the three used to be unreproducible",
          # The three flags and the span: `reference`, because the arbiter is
          # this suite's own codon table and the junction bounds the junction
          # oracle establishes against the reference anchors. A wrong span
          # would not translate to the recorded stop state.
          fields=("productive", "vj_in_frame", "stop_codon",
                  "coding_frame_start", "coding_frame_end"),
          basis="reference")
@pytest.mark.parametrize("stack", sorted(STACKS))
@pytest.mark.parametrize("locus", LOCI)
def test_every_productivity_flag_is_reproducible_from_the_record(locus, stack):
    """The three formulas, over every record the stack produces.

    `reference` basis: the expected values come from the genetic code and from
    arithmetic over exported fields, not from a second view of the engine's
    own answer. The codon table is this suite's own, and the three-valued rule
    over ambiguity codes follows from the code's three stops rather than from
    anything the engine says.
    """
    if stack in B_CELL_ONLY and "tcrb" in locus:
        pytest.skip("somatic hypermutation does not occur in TCR sequences")
    experiment = STACKS[stack](experiment_builder.recombined(locus=locus))
    checked = unresolved = 0
    for index, outcome in enumerate(experiment.run(n=RECORDS, seed=SEED)):
        record = outcome.airr_record(experiment.refdata)
        if record["coding_frame_start"] is None:
            # No resolvable reading frame, so the stop cannot be assessed.
            # `productive` may still be False -- a false predicate decides it
            # whatever the unresolved ones do, which is what three-valued
            # logic means -- but it must not be True. Counted, so a stack that
            # produced nothing but these would not pass by having checked
            # nothing.
            unresolved += 1
            assert record["stop_codon"] is None, (
                f"{locus}/{stack}: record {index} has no coding frame but "
                f"reports stop_codon={record['stop_codon']!r}")
            assert record["productive"] is not True, (
                f"{locus}/{stack}: record {index} is productive with no "
                f"resolvable reading frame; reasons "
                f"{record['productivity_reasons']}")
            continue
        assert expression_oracle.coding_stop_state(record) == record["stop_codon"], (
            f"{locus}/{stack}: record {index} reports "
            f"stop_codon={record['stop_codon']!r} where its own coding span "
            f"{record['coding_frame_start']}..{record['coding_frame_end']} "
            "translates otherwise")
        if record["junction_start"] is not None:
            assert _in_frame(record) == bool(record["vj_in_frame"]), (
                f"{locus}/{stack}: record {index} reports "
                f"vj_in_frame={record['vj_in_frame']!r}, which the junction "
                "bounds and the frame phase do not give")
            assert _productive(record) == record["productive"], (
                f"{locus}/{stack}: record {index} reports "
                f"productive={record['productive']!r} where its five "
                f"predicates give {_productive(record)!r}; reasons "
                f"{record['productivity_reasons']}")
        checked += 1
    assert checked >= RECORDS * 0.5, (
        f"{locus}/{stack}: only {checked} of {RECORDS} records had a "
        f"resolvable frame ({unresolved} did not), so this ran on less than "
        "the fixture it claims")


@verifies("expression.productivity",
          claim="the three predicate fields are the values that make the "
                "conjunction come out right on every record: used as inputs, "
                "they reproduce `productive` exactly across three loci and "
                "six pass stacks. What this does NOT establish is their own "
                "values -- the anchor tri-state and the J tail's frame are "
                "germline-phase arithmetic over provenance, and deriving "
                "either independently is work nobody has done, so a wrong "
                "predicate that happened to keep the conjunction right would "
                "satisfy this",
          fields=("j_tail_in_frame", "v_anchor_accepted", "j_anchor_accepted"),
          basis="consistent")
def test_the_three_exported_predicates_are_consistent_with_the_verdict():
    """Registered separately and honestly.

    These three fields were added so `productive` could be reproduced, and
    this is what that buys: the conjunction is checkable end to end. It is a
    comparison among the record's own fields, so it is `consistent` and the
    mechanism reads `consistent` with it -- which is the status doing its job.
    Raising it needs an independent derivation of the anchor tri-state (which
    must distinguish a rejected anchor from an unresolved one, something
    reading the bases cannot do) and of the J tail's frame against the
    reference J end.
    """
    experiment = (experiment_builder.recombined(locus="human_igh")
                  .polymerase_indels(count=1, insertion_prob=0.0))
    checked = decisive = 0
    for index, outcome in enumerate(experiment.run(n=RECORDS, seed=SEED)):
        record = outcome.airr_record(experiment.refdata)
        if record["coding_frame_start"] is None or record["junction_start"] is None:
            continue
        checked += 1
        assert _productive(record) == record["productive"], (
            f"record {index}: the five predicates give {_productive(record)!r} "
            f"where the record reports {record['productive']!r}")
        # A record the other four conditions satisfy and the J tail does not:
        # the case that makes the fifth predicate load-bearing rather than
        # decorative.
        others = [record["v_anchor_accepted"], record["j_anchor_accepted"],
                  _in_frame(record),
                  expression_oracle.coding_stop_state(record) is False]
        if all(one is True for one in others) and record["j_tail_in_frame"] is False:
            decisive += 1
            assert record["productive"] is False
            assert "j_tail_out_of_frame" in record["productivity_reasons"]
    assert checked >= RECORDS * 0.5, f"only {checked} records were checked"
    assert decisive >= 1, (
        "no record was decided by the J tail alone, so this fixture does not "
        "exercise the predicate the three fields were added for")


def test_the_codon_rule_is_three_valued_over_ambiguity():
    """The genetic code's part of the claim, on its own.

    A codon containing an ambiguous base is a stop if every realization is,
    not a stop if none is, and unknown otherwise. The first version of this
    oracle called any ambiguous codon unknown and so disagreed with the engine
    on 70 of 150 records -- it was less careful than the thing it was checking.
    """
    assert expression_oracle.codon_stop("TAA") is True
    assert expression_oracle.codon_stop("TGG") is False
    # CTN cannot be a stop whichever base the N turns out to be.
    assert expression_oracle.codon_stop("CTN") is False
    # TAN is TAA or TAG (stops) and TAC or TAT (not), so it is unknown.
    assert expression_oracle.codon_stop("TAN") is None
    # TAR is exactly TAA or TAG, both stops.
    assert expression_oracle.codon_stop("TAR") is True
    # TRA is TAA (stop) or TGA (stop) -- also every realization.
    assert expression_oracle.codon_stop("TRA") is True
    assert expression_oracle.codon_stop("NNN") is None
    assert expression_oracle.codon_stop("TA") is None


def test_the_documentation_states_all_five_predicates():
    """The documentation defect the sweep exposed, pinned shut.

    The getting-started page called `productive` a "full triad" of four
    conditions, which made it a complete specification that was wrong: the
    engine applies a fifth, and `j_tail_out_of_frame` appeared nowhere in the
    documentation at all. A consumer implementing the documented rule got a
    different answer from ours on a few percent of indel-bearing records.

    Pinned on the content rather than the prose: each of the five predicates'
    field names has to appear on the page, and the phrase that made the old
    text a complete specification must not.
    """
    from tests.support import paths

    page = (paths.repo_root() / "site_docs" / "getting-started"
            / "first-airr-record.md").read_text(encoding="utf-8")
    for field in ("v_anchor_accepted", "j_anchor_accepted", "vj_in_frame",
                  "stop_codon", "j_tail_in_frame", "coding_frame_start",
                  "coding_frame_end"):
        assert field in page, (
            f"{field} is a predicate of `productive` or the span it is "
            "assessed over, and the field reference does not mention it")
    assert "j_tail_out_of_frame" in page, (
        "the reason a record can satisfy the four documented conditions and "
        "still be unproductive is not in the documentation")
    assert "full triad" not in page, (
        "the page still calls the conditions a 'full triad', which reads as a "
        "complete specification of a conjunction that has five terms")
