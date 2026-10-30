"""The junction window, against the anchors the reference declares.

Every CDR3 benchmark reads these five fields. Clustering thresholds,
repertoire overlap and every length distribution are computed from them, and
`junction_start` / `junction_end` are what a tool is told the window is.
Nothing registered checked any of them before this module: both mechanisms
were `unverified` with zero Python tests.
"""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import junction as junction_oracle
from tests.support.oracles import reconstruction

pytestmark = pytest.mark.behavior

LOCI = ("human_igh", "human_igk")


def _allele_by_name(refdata, segment, name):
    size, getter = {
        "V": (refdata.v_pool_size(), refdata.v_allele),
        "J": (refdata.j_pool_size(), refdata.j_allele),
    }[segment]
    for index in range(size):
        allele = getter(index)
        if allele.name == name:
            return allele
    return None


def _check(exp, *, n, seed):
    """Every disagreement the junction oracle reports, and how much it saw.

    The anchor coordinate and the reference bases are read here, from the
    refdata, and handed to the oracle -- which never imports GenAIRR and so
    cannot look them up itself. That is what keeps the expected value on the
    other side of the thing being checked.
    """
    found, checked, records = [], 0, []
    for outcome in exp.run(n=n, seed=seed):
        record = outcome.airr_record(exp.refdata)
        records.append(record)
        simulation = outcome.final_simulation()
        v_allele = _allele_by_name(exp.refdata, "V", record["v_call"])
        j_allele = _allele_by_name(exp.refdata, "J", record["j_call"])
        assert v_allele is not None and j_allele is not None, (
            f"{record['v_call']} / {record['j_call']} name no allele in the "
            "bound reference, so the record cannot be checked against it")
        found.extend(junction_oracle.disagreements(
            record, simulation,
            reconstruction.segments_by_position(simulation),
            v_anchor=v_allele.anchor, j_anchor=j_allele.anchor,
            v_sequence=bytes(v_allele.seq()), j_sequence=bytes(j_allele.seq())))
        checked += 1
    return tuple(found), checked, records


def _junction_bound_fields_present(record):
    """Both bounds, named literally.

    The claimed fields have to be read by the body that claims them, and the
    oracle reads them one module away -- so they are named here, as the
    coordinate truth tests do for the same reason. The type check is the
    weakest part of this file; the oracle comparison above it is what
    establishes the values.
    """
    for field in ("junction_start", "junction_end"):
        assert isinstance(record[field], int), f"{field} is not an integer"
    assert 0 <= record["junction_start"] < record["junction_end"] <= len(record["sequence"]), (
        f"bounds {record['junction_start']}..{record['junction_end']} fall "
        f"outside a {len(record['sequence'])}-base molecule")


def _junction_string_fields_present(record):
    """The three reported fields, named literally. See above."""
    assert isinstance(record["junction"], str) and record["junction"], (
        "junction is empty")
    assert record["junction_length"] == len(record["junction"]), (
        f"junction_length {record['junction_length']} for a "
        f"{len(record['junction'])}-base junction")
    assert len(record["junction_aa"]) == record["junction_length"] // 3, (
        f"junction_aa has {len(record['junction_aa'])} residues for "
        f"{record['junction_length']} nucleotides")


@verifies("junction.coding_end_processing",
          claim="junction_start is the read index of the first surviving base "
                "of the V anchor codon and junction_end one past the last "
                "surviving base of the J anchor codon, with the anchor "
                "coordinate and the conserved residue read from the reference "
                "rather than from the run. What this does not establish is "
                "that the bound names the CDR3 anchor specifically: 197 of 198 "
                "IGH V alleles carry a second in-frame cysteine upstream, and "
                "a bound displaced onto it satisfies the residue check",
          fields=("junction_start", "junction_end"),
          basis="reference")
@pytest.mark.parametrize("locus", LOCI)
def test_the_junction_bounds_are_the_surviving_anchor_codons(locus):
    """Trimming can eat part of a conserved codon, which is why the bound is
    the first *surviving* base rather than the codon start: the oracle walks
    the provenance for whichever of the three anchor positions is still there.

    Both loci run because the J anchor residue differs -- tryptophan on the
    heavy chain, phenylalanine on the light -- and a check that only ever saw
    TGG would pass against a residue set hardcoded to it.
    """
    exp = experiment_builder.recombined(locus=locus)
    found, checked, records = _check(exp, n=40, seed=5)
    assert found == (), f"{locus}: {found[:3]}"
    assert checked == 40, f"{locus}: only {checked} records were checked"
    for record in records:
        _junction_bound_fields_present(record)


@verifies("junction.assembly",
          claim="junction is the molecule sliced at the recorded bounds -- in "
                "germline orientation, which is probed rather than assumed -- "
                "junction_length is the distance between them, and junction_aa "
                "is the translation of junction under an independent codon "
                "table",
          fields=("junction", "junction_aa", "junction_length"),
          basis="reference")
@pytest.mark.parametrize("locus", LOCI)
def test_the_junction_string_is_the_molecule_between_its_bounds(locus):
    """The three reported fields and the window they claim to describe.

    Of these, the slice identity is the one with an independent side: the
    assembled molecule is written by the assembly passes and the junction by
    the record builder from the coding span. The length and the translation
    are consistency checks over one producer, recorded as such in the oracle's
    circularity map, and they ride along here because they share the fixture.
    """
    exp = experiment_builder.recombined(locus=locus)
    found, checked, records = _check(exp, n=40, seed=5)
    assert found == (), f"{locus}: {found[:3]}"
    assert checked == 40
    for record in records:
        _junction_string_fields_present(record)


@verifies("junction.assembly",
          claim="the junction is reported in germline orientation: on a "
                "reverse-complemented record it is the slice of the "
                "reverse complement of sequence, not of sequence",
          fields=("junction",),
          basis="reference")
def test_the_junction_is_reported_in_germline_orientation():
    """The rule that makes every junction field interpretable on a flipped
    record, and the one an oracle is most likely to get wrong silently -- the
    slice simply comes out as the wrong string.

    Run at prob=1.0 so every record is flipped, with the flag asserted, so a
    build that stopped flipping would fail here rather than quietly passing
    against unflipped reads.
    """
    exp = experiment_builder.recombined().random_strand_orientation(prob=1.0)
    found, checked, records = _check(exp, n=25, seed=5)
    assert found == (), f"{found[:3]}"
    assert checked == 25
    for record in records:
        assert record["junction"], "a flipped record reported an empty junction"
    flags = [o.airr_record(exp.refdata)["rev_comp"]
             for o in exp.run(n=25, seed=5)]
    assert all(flags), "no record was reverse-complemented; the rule is untested"


@verifies("junction.coding_end_processing",
          claim="a mutation pass does not move the junction bounds: "
                "hypermutation substitutes bases without changing which "
                "germline positions survive, so the window is the same window "
                "the unmutated molecule reports under the same seed",
          fields=("junction_start", "junction_end"),
          basis="cross-run")
def test_hypermutation_does_not_move_the_junction_window():
    """A cross-run differential: the bounds are derived from provenance, and
    a substitution changes a base without changing its germline coordinate,
    so adding SHM to the same seed must leave both bounds alone.

    This is the one check here whose expected value comes from a second run
    rather than from the reference, which is why it can see something the
    oracle cannot: a bound that moved because the base under it changed.
    """
    plain = experiment_builder.recombined()
    mutated = experiment_builder.recombined().mutate(model="s5f", count=8)
    before = [o.airr_record(plain.refdata) for o in plain.run(n=40, seed=5)]
    after = [o.airr_record(mutated.refdata) for o in mutated.run(n=40, seed=5)]
    moved = [(index, a["junction_start"], b["junction_start"],
              a["junction_end"], b["junction_end"])
             for index, (a, b) in enumerate(zip(before, after))
             if (a["junction_start"], a["junction_end"])
             != (b["junction_start"], b["junction_end"])]
    assert moved == [], f"{len(moved)} records moved a bound under SHM: {moved[:3]}"
    assert any(b["n_mutations"] for b in after), "no mutation was applied"
