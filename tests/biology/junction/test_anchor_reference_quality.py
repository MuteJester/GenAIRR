"""Which bundled alleles carry a usable anchor, and what happens when one does not.

The junction window is defined by two conserved residues: a cysteine in the V
and a tryptophan or phenylalanine in the J. The engine's rule for "this
reference resolves an anchor" is that the codon at the recorded anchor IS one
of those -- so an allele whose anchor codon is something else is an allele the
engine will refuse to compute a junction for, and every CDR3 metric over a
repertoire silently loses those sequences.

That makes the canonical-codon property a claim about the reference cartridge,
not about any record, and this is where it is tested. The pipeline-space sweep
is what forced the distinction: it ran the junction oracle on TCRB, which the
junction truth tests never did, and reported six alleles as findings. They are
not findings about a record -- the engine handles all six correctly -- but
they are facts about the shipped reference that nobody had written down.

What is pinned here:

- every V and J anchor codon in the bundled BCR references is canonical;
- exactly six TCRB alleles are not, by name;
- for those six, the engine emits no junction at all and says
  `reference_anchor_unresolved` rather than computing a window from an anchor
  it does not trust.

The third is the one that matters for truthfulness. A reference whose
annotation is incomplete is a data problem; a record that invents a junction
from an unusable anchor would be a lie.
"""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import junction as junction_oracle

pytestmark = pytest.mark.behavior

#: Loci whose packaged reference this module reads.
LOCI = ("human_igh", "human_igk", "human_igl", "human_tcrb")

#: The bundled alleles whose anchor codon is not the conserved residue, by
#: locus and segment. Pinned by name: each one is a sequence the engine will
#: not produce a junction for, so the set is a coverage fact about the shipped
#: reference and a change to it should be deliberate.
#:
#: TRBJ2-2P is a pseudogene, where no usable anchor is expected. The other five
#: are not marked as pseudogenes, and `functional_status` is None for all six,
#: so nothing on the public allele surface distinguishes them from an allele
#: with a usable anchor except reading the codon.
NON_CANONICAL_ANCHORS = {
    ("human_tcrb", "V"): {"TRBV17*01", "TRBV17*02", "TRBV7-3*02", "TRBV7-3*03"},
    ("human_tcrb", "J"): {"TRBJ2-2P*01", "TRBJ2-7*02"},
}

CANONICAL = {"V": junction_oracle.CYSTEINE,
             "J": junction_oracle.TRYPTOPHAN_OR_PHENYLALANINE}


def _alleles(refdata, segment):
    size, getter = {"V": (refdata.v_pool_size(), refdata.v_allele),
                    "J": (refdata.j_pool_size(), refdata.j_allele)}[segment]
    return [getter(index) for index in range(size)]


def _non_canonical(refdata, segment):
    """``{name: codon}`` for alleles whose anchor codon is not conserved."""
    from tests.support import alleles as allele_support

    out = {}
    for allele in _alleles(refdata, segment):
        codon = junction_oracle.reference_anchor_codon(
            allele_support.sequence(allele), allele.anchor)
        if codon is not None and codon not in CANONICAL[segment]:
            out[allele.name] = codon
    return out


@pytest.mark.parametrize("locus", LOCI)
@pytest.mark.parametrize("segment", ["V", "J"])
def test_the_bundled_anchor_codons_are_canonical_except_where_pinned(locus, segment):
    """The reference-data claim, per locus and segment.

    Measured over the bundled set: 198 IGH V and 7 IGH J, 168 IGK V and 8 IGK
    J, 181 IGL V and 10 IGL J anchor codons are all canonical; TCRB has 98 V of
    which 4 are not and 16 J of which 2 are not.
    """
    refdata = experiment_builder.recombined(locus=locus).refdata
    found = _non_canonical(refdata, segment)
    expected = NON_CANONICAL_ANCHORS.get((locus, segment), set())
    assert set(found) == expected, (
        f"{locus} {segment}: alleles whose anchor codon is not conserved "
        f"changed to {found}. Each one is a sequence the engine will not "
        "produce a junction for, so a change here is a change in what the "
        "shipped reference can be used for.")
    assert _alleles(refdata, segment), f"{locus} has no {segment} alleles at all"


@verifies("junction.coding_end_processing",
          claim="an allele whose anchor codon the reference does not resolve "
                "produces no junction at all -- the bounds, the string, the "
                "length and the frame flag all null or empty, with "
                "reference_anchor_unresolved in the reasons -- rather than a "
                "window computed from an anchor the engine does not trust. "
                "Asserted on all six such alleles in the bundled TCRB "
                "reference, drawn deliberately by weighting",
          fields=("junction_start", "junction_end"),
          basis="reference")
@pytest.mark.parametrize("segment,name", sorted(
    (segment, name)
    for (locus, segment), names in NON_CANONICAL_ANCHORS.items()
    for name in names))
def test_an_unresolvable_anchor_produces_no_junction_rather_than_a_wrong_one(
        segment, name):
    """The truthfulness claim. A data gap must not become a false metric.

    Each allele is forced by weighting rather than waited for: TRBV17*01 is one
    of 98 V alleles, so an unweighted run reaches it a few times in two
    hundred records and a test that relied on chance would pass by never
    drawing it.
    """
    weights = {f"{segment.lower()}_allele_weights": {name: 50_000}}
    experiment = experiment_builder.recombined(locus="human_tcrb", **weights)
    checked = 0
    for outcome in experiment.run(n=25, seed=2):
        record = outcome.airr_record(experiment.refdata)
        if record[f"{segment.lower()}_call"] != name:
            continue
        checked += 1
        assert record["junction_start"] is None and record["junction_end"] is None, (
            f"{name}: a junction window was computed from an anchor the "
            f"reference does not resolve: "
            f"{record['junction_start']}..{record['junction_end']}")
        assert record["junction"] == "" and record["junction_aa"] == ""
        assert record["junction_length"] is None
        assert record["vj_in_frame"] is None, (
            f"{name}: vj_in_frame is {record['vj_in_frame']!r}, which a "
            "consumer reads as out-of-frame rather than as unresolved")
        assert "reference_anchor_unresolved" in (record["productivity_reasons"] or ()), (
            f"{name}: nothing in the record says why there is no junction: "
            f"{record['productivity_reasons']}")
    assert checked >= 5, (
        f"only {checked} of 25 records drew {name}; the weighting is not "
        "forcing the allele and this test is passing on an empty loop")


def test_the_reason_counts_one_unresolved_anchor_per_segment():
    """How much a consumer can tell from the reasons, stated exactly.

    `reference_anchor_unresolved` appears once per unresolved anchor -- twice
    when both the V and the J are. Measured over 200 TCRB records: 172 with
    neither, 25 with the J alone, 3 with the V alone, each carrying exactly
    that many occurrences.

    What the record does NOT say is WHICH anchor. The count is information and
    the attribution is missing, so a consumer seeing one occurrence cannot
    tell a V problem from a J problem without looking the alleles up itself.
    Pinned here because it is a limit of the exported metadata, not a defect in
    it, and a reader of this suite should be able to find it stated.
    """
    experiment = experiment_builder.recombined(locus="human_tcrb")
    unresolved_v = NON_CANONICAL_ANCHORS[("human_tcrb", "V")]
    unresolved_j = NON_CANONICAL_ANCHORS[("human_tcrb", "J")]
    seen = {}
    for outcome in experiment.run(n=200, seed=4):
        record = outcome.airr_record(experiment.refdata)
        reasons = record["productivity_reasons"] or ()
        occurrences = sum(1 for reason in reasons
                          if reason == "reference_anchor_unresolved")
        expected = ((record["v_call"] in unresolved_v)
                    + (record["j_call"] in unresolved_j))
        assert occurrences == expected, (
            f"{record['v_call']}/{record['j_call']}: the reason appears "
            f"{occurrences} times where {expected} anchors are unresolved")
        seen[expected] = seen.get(expected, 0) + 1
    assert seen.get(1, 0) >= 10, (
        f"only {seen.get(1, 0)} of 200 records had exactly one unresolved "
        "anchor, so the per-segment count was barely exercised")
    assert seen.get(0, 0) >= 100, f"the run was mostly unresolved: {seen}"
