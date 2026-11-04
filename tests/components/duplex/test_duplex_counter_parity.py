"""A counter the native record exports, the duplex record must also export.

The duplex projection does not reuse the native builder: it walks its own
observation ledger and assembles the record field by field. Three separate
**hand-written lists of counter field names** decide which counters it
emits -- one for zero defaults, one for nulling them when the observation
history is unattributed, and one in the repair projection. A counter added
to the native record and not to those lists is not a crash and not an empty
column: it is a column missing from one of two record shapes that are
supposed to share a schema.

That has now happened twice. ``n_observable_mutations`` was absent, which
"wrote a blank column into every export" (the comment that fix left
behind). ``n_ambiguous_bases`` was absent, and worse, the ledger had an
explicit ``Some("corrupt.ns") => ()`` arm discarding the pass because no
field existed to put it in. Both were found by reading, not by a test.
This is the test, so the third time fails instead.

It is the same defect class PR #44 fixed for productivity reasons -- a
vocabulary spelled out as string literals in several places with nothing
tying them together -- and it is checked here the same way: by comparing
the two shapes against each other rather than against a list that would
itself need maintaining.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge
from GenAIRR.experimental import (
    CoupledDuplexRepair,
    DuplexJoin,
    project_duplex_simulation,
)

pytestmark = pytest.mark.behavior


def _native_record() -> dict:
    """A native record carrying every corruption counter non-zero.

    The counts are deliberately above 1: a substitution pass draws its
    replacement from all four bases including the base already there, so a
    requested count of 1 realizes as 0 about a quarter of the time. Field
    *presence* is what this module checks, so a zero would not actually
    break it -- but a record built with the passes switched on is the
    honest subject, and it keeps the parity check meaningful if a future
    field is only emitted when its pass ran.
    """
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=4)
        .pcr_amplify(count=4)
        .sequencing_errors(count=4)
        .ambiguous_base_calls(count=4)
        .polymerase_indels(count=2)
    )
    result = exp.run_records(n=1, seed=17)
    return result.outcomes[0].airr_record(exp.refdata)


def _duplex_record() -> dict:
    """The duplex projection's record, built from its own ledger walk."""
    refdata = ge.RefDataConfig.vdj()
    refs = ("AAATGT", "CCA", "TTTAAA")
    refdata.add_v_allele("V*01", "V", refs[0].encode(), anchor=3)
    refdata.add_d_allele("D*01", "D", refs[1].encode())
    refdata.add_j_allele("J*01", "J", refs[2].encode(), anchor=0)
    script = CoupledDuplexRepair(*refs)
    for name in ("vd", "dj"):
        script = (
            script.join(
                DuplexJoin(name)
                .open("left", offset=0)
                .open("right", offset=0)
                .synapse(right_start=6 if name == "vd" else 3)
                .ligate()
            )
        )
    seed = script.seed(
        refdata,
        v_call="V*01",
        d_call="D*01",
        j_call="J*01",
        reference_ends_are_coding_ends=True,
    )
    return project_duplex_simulation(seed, refdata)["record"]


def _counter_fields(record: dict) -> set[str]:
    """The record's count-valued provenance fields.

    Selected by name rather than from a list this test holds, so a new
    counter joins the check by being named like one. ``n_`` covers the
    event-ledger counts and ``_length`` the end-loss amounts, which are the
    two shapes the duplex lists already enumerate.
    """
    return {
        key
        for key in record
        if (key.startswith("n_") or key.endswith("_length"))
        and not key.startswith("np")
    }


def test_every_native_counter_field_exists_on_a_duplex_record() -> None:
    """Field-for-field parity, in the direction that has actually broken.

    Asserted one-way on purpose. The duplex record carries repair fields
    with no native counterpart, which is correct -- native records have no
    repair. What is never correct is a counter the native shape publishes
    and the duplex shape silently omits, because both flow into the same
    export with the same column order.
    """
    native = _counter_fields(_native_record())
    duplex = set(_duplex_record())

    missing = sorted(field for field in native if field not in duplex)
    assert not missing, (
        "the duplex projection omits counter field(s) the native record "
        f"exports: {missing}. These are emitted from hand-written field "
        "lists in engine_rs/src/duplex_repair/observation.rs and "
        "projection.rs -- add the field there, or the export writes a "
        "blank column for every duplex record."
    )


def test_the_parity_check_is_looking_at_real_counters() -> None:
    """Guard the guard.

    If ``_counter_fields`` ever selected nothing -- a rename of the
    counters, a change in the record's key spelling -- the test above would
    pass by vacuity and keep passing forever. Pin that it sees the counters
    whose omission caused the two known incidents, plus enough others to
    show the selector is not matching by accident.
    """
    native = _counter_fields(_native_record())
    for field in (
        "n_mutations",
        "n_observable_mutations",
        "n_ambiguous_bases",
        "n_pcr_errors",
        "n_quality_errors",
        "n_indels",
    ):
        assert field in native, (
            f"{field} is no longer selected as a counter field, so the "
            "parity check above has stopped covering it"
        )
