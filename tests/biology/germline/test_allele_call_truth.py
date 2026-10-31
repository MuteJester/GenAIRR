"""The call names the allele whose bases were actually used.

`v_call` / `d_call` / `j_call` are the most-used output GenAIRR has: they are
what an allele caller is scored against. Nothing registered checked them, and
the one external check that touched `v_call` compared the *length* of the
allele it names -- which cannot separate 195 of the 198 IGH V alleles,
because they share an exact length with another.

The engine's own re-projection comparison does catch an edited call, but by
recomputing it from the same assignment: it agrees with itself whichever
allele the projection names. So the check here resolves the reference by the
exported NAME and compares the recorded germline bases to it.
"""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import reconstruction

pytestmark = pytest.mark.behavior

LOCI = ("human_igh", "human_igk")


def _alleles_by_name(refdata):
    pools = {
        "V": (refdata.v_pool_size(), refdata.v_allele),
        "D": (refdata.d_pool_size(), refdata.d_allele),
        "J": (refdata.j_pool_size(), refdata.j_allele),
    }
    return {segment: {get(index).name: bytes(get(index).seq())
                      for index in range(size)}
            for segment, (size, get) in pools.items()}


def _standard_calls(record):
    """The three standard calls, named literally.

    A claimed field has to be readable in the body that claims it, and an
    f-string is invisible to the static check -- which rejected the first
    version of this module for exactly that.
    """
    return {
        "V": record.get("v_call"),
        "D": record.get("d_call"),
        "J": record.get("j_call"),
    }


def _recalled_calls(record):
    """The three recalled calls, named literally. See above."""
    return {
        "V": record.get("recalled_v_call") or "",
        "D": record.get("recalled_d_call") or "",
        "J": record.get("recalled_j_call") or "",
    }


def _resolver(record, by_name):
    """Resolve each segment's reference by the name the record exports."""
    def resolve(segment):
        return by_name[segment].get(record.get(f"{segment.lower()}_call"))
    return resolve


@verifies("germline.allele_choice",
          claim="the recorded germline bases are the bases of the allele the "
                "exported call NAMES, not merely of the allele the engine "
                "remembers assigning. Resolving by name catches a call "
                "swapped to a same-length allele on 40 of 40 records where "
                "the assignment-resolved walk catches none, which matters "
                "because 195 of 198 IGH V alleles share a length with another",
          fields=("v_call", "d_call", "j_call"),
          basis="reference")
@pytest.mark.parametrize("locus", LOCI)
def test_each_call_names_the_allele_its_bases_came_from(locus):
    """Read `clean` rather than `mismatches`: a call naming a shorter allele
    puts recorded coordinates past its end, which the walk counts as
    out-of-range rather than as a mismatch. Probing `mismatches` alone made a
    shorter-allele swap look invisible.
    """
    exp = experiment_builder.recombined(locus=locus)
    by_name = _alleles_by_name(exp.refdata)
    checked = compared = 0
    for outcome in exp.run(n=40, seed=5):
        record = outcome.airr_record(exp.refdata)
        simulation = outcome.final_simulation()
        for segment, name in _standard_calls(record).items():
            if name:
                assert name in by_name[segment], (
                    f"the {segment} call names {name!r}, which is no allele in "
                    f"the bound reference")
        check = reconstruction.templated_check(
            simulation, exp.refdata,
            inverted=("D",) if record.get("d_inverted") else (),
            resolve=_resolver(record, by_name))
        assert check.clean, (
            f"{locus}: the bases do not come from the alleles the calls name: "
            f"{check.mismatches[:3]} "
            f"(out_of_range={check.out_of_range}, unresolved={check.unresolved_allele})")
        compared += check.checked
        checked += 1
    assert checked == 40
    assert compared > 1000, (
        f"{locus}: only {compared} positions were compared; the walk is "
        f"finding no templated bases to check")


@verifies("germline.allele_choice",
          claim="every name a recalled call lists is an allele of the bound "
                "reference, so a consumer can resolve it. Containing the "
                "generating allele is NOT promised and is not asserted: "
                "measured 200/200 for V and J but 196/200 for D, the four "
                "misses being same-gene or neighbouring alleles that fit the "
                "16-24 surviving D bases better under hypermutation",
          fields=("recalled_v_call", "recalled_d_call", "recalled_j_call"),
          basis="reference")
def test_every_recalled_name_resolves_in_the_reference():
    """The recalled calls are inference, not truth, and the documentation says
    so. What a consumer still needs is that the names resolve -- a tie set
    naming an allele the reference does not have cannot be scored at all.
    """
    exp = experiment_builder.recombined().mutate(model="s5f", count=8)
    by_name = _alleles_by_name(exp.refdata)
    listed = ties = 0
    for outcome in exp.run(n=120, seed=5):
        record = outcome.airr_record(exp.refdata)
        for segment, value in _recalled_calls(record).items():
            names = [part.strip() for part in value.split(",") if part.strip()]
            for name in names:
                assert name in by_name[segment], (
                    f"the recalled {segment} call lists {name!r}, which is no "
                    f"allele in the bound reference")
                listed += 1
            ties += len(names) > 1
    assert listed > 300, f"only {listed} recalled names were checked"
    assert ties, "no tie set was produced, so the comma-separated form is untested"


@verifies("germline.allele_choice",
          claim="on a revised receptor, original_v_call names the superseded "
                "allele and v_call the one now in the molecule: the recorded "
                "V bases match the allele v_call names on 60 of 60 revised "
                "records and the allele original_v_call names on 0 of 60. So "
                "the revision replaced the bases rather than relabelling them, "
                "and the pre-revision truth is recoverable",
          fields=("original_v_call",),
          basis="reference")
def test_a_revised_receptor_keeps_the_superseded_call_and_uses_the_new_bases():
    """`original_v_call` is empty unless `receptor_revision` runs, which is
    why the rest of this module cannot reach it.

    The discriminating half is the 0 of 60: a field that merely echoed the
    current call, or a revision that relabelled without re-assembling, would
    both pass the first assertion and fail this one.
    """
    exp = experiment_builder.recombined().receptor_revision(prob=1.0)
    by_name = _alleles_by_name(exp.refdata)
    revised = 0
    for outcome in exp.run(n=60, seed=5):
        record = outcome.airr_record(exp.refdata)
        if not record.get("receptor_revision_applied"):
            continue
        revised += 1
        original = record["original_v_call"]
        assert original in by_name["V"], (
            f"original_v_call names {original!r}, no allele of the reference")
        assert original != record["v_call"], (
            "the superseded call equals the current one, so the field carries "
            "nothing a consumer could use")

        def resolve_with(v_name):
            def resolve(segment):
                if segment == "V":
                    return by_name["V"].get(v_name)
                return by_name[segment].get(record.get(f"{segment.lower()}_call"))
            return resolve

        inverted = ("D",) if record.get("d_inverted") else ()
        current = reconstruction.templated_check(
            outcome.final_simulation(), exp.refdata, inverted=inverted,
            resolve=resolve_with(record["v_call"]))
        assert current.clean, (
            f"the V bases do not come from {record['v_call']}, which v_call "
            f"names: {current.mismatches[:3]}")
        superseded = reconstruction.templated_check(
            outcome.final_simulation(), exp.refdata, inverted=inverted,
            resolve=resolve_with(original))
        assert not superseded.clean, (
            f"the V bases also match {original}, the superseded allele, so "
            f"this record does not establish that the revision re-assembled")
    assert revised == 60, f"only {revised} of 60 records were revised"
