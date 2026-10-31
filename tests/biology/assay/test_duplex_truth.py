"""Double-strand break repair, against scripts authored outside the engine.

`assay.duplex_formation` was the last mechanism with no registered test. It is
the surface where GenAIRR models a junction as DNA chemistry rather than as a
length: a break is opened, ends are resected, bases are added or copied back
from the opposite strand, the ends synapse and ligate. The repair is authored
by the caller as a `CoupledDuplexRepair` script, which is what makes it
checkable here -- the expected molecule and the expected per-base accounting
are written by hand in this module, and the engine is asked to reproduce them.

The claim that matters most is the one a read cannot carry. A base resected
from one strand and restored by template-directed fill leaves a molecule
identical to the one a blunt ligation produces: `sequence` is the same string,
and only `duplex_repair` records that a base was cut and resynthesised. Any
analysis that treats the read as the whole truth about a junction is blind to
that, which is the reason the plane exists and the reason it needs a test.
"""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.cases import duplex as duplex_case
from tests.support.oracles import duplex as duplex_oracle

import GenAIRR as ga

pytestmark = pytest.mark.behavior

#: The three authored repairs, with the molecule each one must produce and the
#: physical accounting it must report. Written out rather than derived: these
#: are the expectation, and deriving them from the engine's own numbers is the
#: mistake this module exists not to make.
#:
#: The fixture's alleles are V=AAATGT, D=CCC, J=TTTAAA, joined at offset 0 with
#: no trimming, so a blunt repair spells them out in order. Thirty original
#: nucleotides is the two strands of those fifteen bases.
AUTHORED = (
    ("blunt", {}, "AAATGTCCCTTTAAA",
     {"original_nucleotides": 30, "resected_nucleotides": 0,
      "retained_nucleotides": 30, "templated_synthesis": 0,
      "untemplated_additions": 0, "hairpin_transfer_events": 0}),
    ("two-added", dict(added="GG"), "AAATGTGGCCCTTTAAA",
     {"original_nucleotides": 30, "resected_nucleotides": 0,
      # Two bases added to each strand: the mate is the reverse complement,
      # which is what makes the product a duplex rather than a nick.
      "retained_nucleotides": 34, "templated_synthesis": 0,
      "untemplated_additions": 4, "hairpin_transfer_events": 0}),
    ("four-added", dict(added="ACGT"), "AAATGTACGTCCCTTTAAA",
     {"original_nucleotides": 30, "resected_nucleotides": 0,
      "retained_nucleotides": 38, "templated_synthesis": 0,
      "untemplated_additions": 8, "hairpin_transfer_events": 0}),
    ("resect-and-fill", dict(restore=True), "AAATGTCCCTTTAAA",
     {"original_nucleotides": 30, "resected_nucleotides": 1,
      # One base cut and one resynthesised: the molecule is the blunt one and
      # the accounting is not.
      "retained_nucleotides": 30, "templated_synthesis": 1,
      "untemplated_additions": 0, "hairpin_transfer_events": 0}),
)


def _record(label, **script):
    """One record from the authored script named by ``script``."""
    reference = duplex_case.references()
    model = duplex_case.model(reference, rows=[duplex_case.history(**script)])
    experiment = ga.Experiment.on(reference).recombine(duplex_repair=model)
    outcome = experiment.run(n=1, seed=0)[0]
    return outcome.airr_record(reference)


@verifies("assay.duplex_formation",
          claim="the molecule is what the authored repair script ligates, "
                "base for base, and duplex_repair accounts for every "
                "nucleotide in it: how many were already there, how many were "
                "resected, how many were added untemplated and how many were "
                "resynthesised off the opposite strand. Both expectations are "
                "written by hand from the script, not read back from the "
                "engine",
          fields=("sequence", "duplex_repair"),
          basis="reference")
@pytest.mark.parametrize("label,script,expected_sequence,expected_counts", AUTHORED,
                         ids=[entry[0] for entry in AUTHORED])
def test_the_molecule_and_its_accounting_are_what_the_script_authored(
        label, script, expected_sequence, expected_counts):
    """The authored comparison. Nothing on the expected side came from the run.

    The scripts differ in exactly one respect each -- bases added at the V-D
    join, or a base resected and filled -- so a failure names which repair
    operation the engine stopped honouring.
    """
    record = _record(label, **script)
    assert record["sequence"] == expected_sequence, (
        f"{label}: the repair ligated {record['sequence']!r} where the script "
        f"says {expected_sequence!r}")
    physical = record["duplex_repair"]["counts"]["physical"]
    assert physical == expected_counts, (
        f"{label}: the physical accounting is {physical}, not the "
        f"{expected_counts} the script implies")
    found = duplex_oracle.disagreements(record, offered_histories=("blunt",))
    assert found == (), f"{label}: {len(found)} plane disagreement(s); first: {found[0]}"


@verifies("assay.duplex_formation",
          claim="a resected base restored by template-directed fill is "
                "invisible in the read and recorded in the plane: the "
                "molecule is byte-identical to the blunt repair's, while "
                "duplex_repair reports one resected nucleotide and one "
                "resynthesised, and the restored base is marked "
                "templated_synthesis rather than passed off as original",
          fields=("sequence", "duplex_repair"),
          basis="reference")
def test_a_resected_and_refilled_base_is_invisible_in_the_read():
    """Why this mechanism cannot be verified from `sequence` alone.

    This is the measurement that justifies the whole plane. Two repairs that
    differ in what happened to the DNA produce the same read, so a benchmark
    scoring junction inference against the read is scoring against a molecule
    that cannot distinguish them.
    """
    blunt = _record("blunt")
    restored = _record("resect-and-fill", restore=True)
    assert blunt["sequence"] == restored["sequence"], (
        "the fixture no longer produces the same molecule both ways, so it no "
        "longer shows that the read cannot tell them apart")

    blunt_counts = blunt["duplex_repair"]["counts"]["physical"]
    restored_counts = restored["duplex_repair"]["counts"]["physical"]
    assert blunt_counts["resected_nucleotides"] == 0
    assert blunt_counts["templated_synthesis"] == 0
    assert restored_counts["resected_nucleotides"] == 1, (
        "the resection is not reported, so the only record of it is gone")
    assert restored_counts["templated_synthesis"] == 1, (
        "the refilled base is not reported as resynthesised")

    kinds = duplex_oracle.kind_counts(restored)
    assert kinds.get("templated_synthesis") == 1, (
        f"the restored base is not marked templated_synthesis: {kinds}")
    assert kinds.get("original_strand") == 29, (
        f"the restored base is being counted as an original one: {kinds}")


@verifies("assay.duplex_formation",
          claim="a repair that leaves a gap cannot ligate. A resection with "
                "no fill is refused at plan compilation, naming the "
                "unresolved nick, rather than producing a molecule with a "
                "silent single-strand break in it",
          fields=("duplex_repair",),
          basis="reference")
def test_a_resection_with_no_fill_cannot_ligate():
    """The biochemical constraint, asserted where a caller would hit it.

    Ligase needs a 3'-OH and a 5'-phosphate meeting with no gap between them.
    A script that cuts a base and does not restore it leaves exactly that gap,
    and the engine refuses the plan instead of emitting a product whose
    accounting would not add up. Checked here because a silently accepted gap
    would be a record claiming a repair that chemistry does not allow.
    """
    from GenAIRR.experimental import CoupledDuplexRepair, DuplexJoin

    reference = duplex_case.references()
    script = CoupledDuplexRepair("AAATGT", "CCC", "TTTAAA")
    unfilled = (DuplexJoin("vd").open("left", offset=0).open("right", offset=0)
                .resect("left", "top", count=1).synapse(right_start=6).ligate())
    script = script.join(unfilled)
    script = script.join(DuplexJoin("dj").open("left", offset=0)
                         .open("right", offset=0).synapse(right_start=3).ligate())
    model = duplex_case.model(
        reference, rows=[duplex_case.history_from_script(script)])
    with pytest.raises(ValueError, match="nick|gap"):
        ga.Experiment.on(reference).recombine(duplex_repair=model).run(n=1, seed=0)
