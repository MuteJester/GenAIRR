"""The amino-acid projection must be the translation, in the frame it claims.

`sequence_aa` cannot be checked by inspecting it: whether it is right depends
on a reading frame that is not in the field. The existing tests assert a
property of it -- that no stop codon appears in a productive record -- and
compare it to hardcoded strings on synthetic micro-fixtures. Neither
establishes that the amino acids are the translation of the bases.

The frame was read from the engine rather than guessed, then confirmed against
records. See `tests/support/oracles/projection` for the derivation.
"""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import projection

pytestmark = pytest.mark.behavior

LOCI = ("human_igh", "human_igk", "human_tcrb")

STACKS = {
    "plain": lambda exp: exp,
    "indels": lambda exp: exp.polymerase_indels(count=2),
    "mutated": lambda exp: exp.mutate(model="s5f", count=5),
    "reverse-complemented": lambda exp: exp.random_strand_orientation(prob=1.0),
    "end-loss": lambda exp: exp.end_loss_5prime(length=4),
}
B_CELL_ONLY = frozenset({"mutated"})


def _native_phase(refdata, simulation):
    """The V allele's own coding phase, from the reference table.

    ``anchor % 3``. A curated reference can carry explicit ``coding_phase``
    evidence that takes precedence, which the packaged references used here do
    not; if one did, this would have to read it, and the oracle deliberately
    does not reimplement that precedence.
    """
    identifier = simulation.v_allele_id()
    if identifier is None:
        return None
    allele = refdata.v_allele(identifier)
    anchor = allele.anchor() if callable(allele.anchor) else allele.anchor
    return None if anchor is None else anchor % 3


def test_the_codon_table_is_complete_before_anything_relies_on_it():
    """A 63-codon table would silently produce an X and the comparison would
    fail somewhere unrelated."""
    table = projection.codon_table()
    assert len(table) == 64, f"the codon table has {len(table)} entries, not 64"
    assert sorted(k for k, v in table.items() if v == "*") == ["TAA", "TAG", "TGA"]
    assert set("".join(table)) == set("ACGT")
    assert projection.translate("ATGTGGTAA") == "MW*"
    assert projection.translate("ATGTG") == "M", "a partial trailing codon is dropped"
    assert projection.translate("ATGNNN") == "MX", "an ambiguous codon is X, not a stop"


@verifies("projection.translation",
          claim="sequence_aa is the translation of the rearrangement in the "
                "frame the V allele's coding phase and the bases removed from "
                "its 5' end imply -- taken in germline orientation, so the "
                "reverse complement of the read when the read is on the other "
                "strand. Checked against an independent codon table, not "
                "against the engine's",
          fields=("sequence_aa",))
@pytest.mark.parametrize("stack", sorted(STACKS))
@pytest.mark.parametrize("locus", LOCI)
def test_sequence_aa_is_the_translation_in_the_derived_frame(locus, stack):
    if "tcr" in locus and stack in B_CELL_ONLY:
        pytest.skip("somatic hypermutation does not occur in TCR sequences")
    exp = STACKS[stack](experiment_builder.recombined(locus=locus))
    checked = frameless = 0
    for outcome in exp.run(n=25, seed=61):
        record = outcome.airr_record(exp.refdata)
        native = _native_phase(exp.refdata, outcome.final_simulation())
        if native is None:
            pytest.fail(f"{locus}: the V allele carries no anchor, so no frame "
                        "can be derived; the fixture is wrong")
        if not record["sequence_aa"]:
            frameless += 1
            continue
        found = projection.translation_disagreement(record, native_phase=native)
        assert found is None, (
            f"{locus}/{stack}: sequence_aa is {found[0][:40]!r} but the "
            f"translation of its own bases is {found[1][:40]!r}")
        checked += 1
    assert checked, (
        f"{locus}/{stack}: every record was frameless ({frameless}), so "
        "nothing was translated and this establishes nothing")


def test_a_frameless_record_is_refused_rather_than_reported_as_agreeing():
    """An empty sequence_aa means the engine found no coding frame. Returning
    "no disagreement" for it would be the ambiguity this suite keeps removing:
    "agreed" and "there was nothing to compare" must not look alike."""
    exp = experiment_builder.recombined()
    record = dict(exp.run_records(n=1, seed=61)[0])
    frameless = {**record, "sequence_aa": ""}
    with pytest.raises(ValueError, match="no coding frame"):
        projection.translation_disagreement(frameless, native_phase=0)


@pytest.mark.parametrize("shift", [1, 2])
def test_a_translation_in_the_wrong_frame_is_caught(shift):
    """The defect the existing property test cannot see: a frame shifted by one
    or two still produces plausible amino acids, and a shifted translation of a
    productive sequence very often still has no stop codon in it."""
    exp = experiment_builder.recombined()
    outcome = exp.run(n=1, seed=61)[0]
    record = dict(outcome.airr_record(exp.refdata))
    native = _native_phase(exp.refdata, outcome.final_simulation())
    assert projection.translation_disagreement(record, native_phase=native) is None
    found = projection.translation_disagreement(record, native_phase=native + shift)
    assert found is not None, f"a frame shifted by {shift} went undetected"


def test_a_translation_of_the_wrong_strand_is_caught():
    """Models the projection being written in read orientation rather than
    germline orientation -- which is what the first draft of this check
    assumed, and it disagreed on 30 of 30 reverse-complemented records."""
    exp = experiment_builder.recombined().random_strand_orientation(prob=1.0)
    outcome = exp.run(n=1, seed=61)[0]
    record = dict(outcome.airr_record(exp.refdata))
    native = _native_phase(exp.refdata, outcome.final_simulation())
    assert record["rev_comp"] is True
    assert projection.translation_disagreement(record, native_phase=native) is None
    # Claim the read is forward, so the oracle translates the wrong strand.
    forward = {**record, "rev_comp": False}
    assert projection.translation_disagreement(forward, native_phase=native) is not None, (
        "translating the read strand instead of the germline strand agreed")
