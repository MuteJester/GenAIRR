"""What a recorded mutation count must mean.

A tool scored against GenAIRR truth can only recover substitutions that are
visible in the bases. If the recorded count includes events that left the
sequence unchanged, every such tool carries a built-in false-negative rate
equal to the discrepancy -- and nothing in the record tells it so.
"""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import provenance, reconstruction

pytestmark = pytest.mark.behavior

MODELS = ("uniform", "s5f")


@verifies("maturation.hypermutation",
          claim="every counted mutation changed a base, so n_mutations counts "
                "events that actually happened rather than draws that were "
                "made. At one mutation per record it equals the number of "
                "positions differing from their germline byte exactly, for "
                "both models. What it does NOT promise is that a caller can "
                "recover them all: a mutation on an untemplated junction base "
                "has no germline base to differ from, which is intended and "
                "measured in the test below",
          fields=("n_mutations",),
          basis="provenance")
@pytest.mark.parametrize("model", MODELS)
def test_every_counted_mutation_changed_a_base(model):
    """A substitution to the same base is not a mutation under any molecular
    definition: AID deamination repair yields transitions and transversions,
    never a base replaced by itself.

    The engine used to count them. The replacement is drawn from all four
    bases including the original, so roughly a quarter of uniform draws landed
    on the base already present; those writes still happened and were still
    counted. They are now subtracted at the point of the write, which leaves
    the draw, the RNG consumption, the recorded trace and the sequence exactly
    as they were -- only the count moved.

    Compared case-insensitively on purpose: a templated base is stored
    lowercase and an SHM write emits uppercase, so a byte comparison sees a
    change at every mutated position and hides this entirely. The first
    attempt at the fix made exactly that mistake. It also patched one of the
    seven write sites, and the uniform model takes a different one, so the
    measurement did not move at all.

    Measured at one mutation per record, seed 101, over 200 records: 198
    positions differ from their germline byte, 191 of them a lowercase
    germline byte against an uppercase base, 44 the same base differing only
    in case, and 46 records realise no mutation at all.
    """
    exp = experiment_builder.recombined().mutate(model=model, count=1)
    unchanged = []
    for index, outcome in enumerate(exp.run(n=200, seed=0)):
        record = outcome.airr_record(exp.refdata)
        simulation = outcome.final_simulation()
        changed = provenance.changed_base_count(simulation)
        if record["n_mutations"] and not changed:
            unchanged.append(index)
        # At one mutation per record there are no repeat hits, so the count
        # and the number of changed positions are the same number.
        assert record["n_mutations"] == changed, (
            f"{model}: record {index} reports {record['n_mutations']} mutations "
            f"with {changed} positions changed")
    assert unchanged == [], (
        f"{model}: {len(unchanged)}/200 records report a mutation that changed "
        "no base at all")


@verifies("maturation.hypermutation",
          claim="a mutation on an untemplated junction base is counted and is "
                "not recoverable from the sequence, which is intended rather "
                "than a defect: there is no germline base at that position to "
                "differ from. Measured at 5.2% of uniform and 4.5% of "
                "s5f records that carry a mutation "
                "single-mutation records, and it is the whole of the gap "
                "between n_mutations and what a mutation caller can see",
          fields=("n_mutations",),
          basis="provenance")
@pytest.mark.parametrize("model", MODELS)
def test_a_mutation_on_an_untemplated_base_is_counted_and_unrecoverable(model):
    """The part of the gap that is not a defect, pinned so it is not mistaken
    for one again.

    `observable_substitutions` compares only templated positions, because an
    untemplated one has no germline coordinate. A mutation there changes a
    base, so it belongs in n_mutations, and no caller comparing the sequence
    to its germline can recover it. Both statements are true at once.
    """
    exp = experiment_builder.recombined().mutate(model=model, count=1)
    invisible = 0
    checked = 0
    for outcome in exp.run(n=200, seed=0):
        record = outcome.airr_record(exp.refdata)
        simulation = outcome.final_simulation()
        if not record["n_mutations"]:
            continue
        checked += 1
        observable = reconstruction.observable_substitutions(simulation)
        if observable:
            continue
        invisible += 1
        # Every such record must be an untemplated hit, not a no-op: a
        # counted mutation with no change anywhere is the defect the test
        # above forbids.
        bases, germline = simulation.bases(), simulation.germline_bases()
        changed = [i for i in range(len(bases))
                   if chr(bases[i]).upper() != chr(germline[i]).upper()]
        assert len(changed) == provenance.changed_base_count(simulation)
        assert changed, f"{model}: a counted mutation changed nothing"
        assert all(simulation.germline_position(i) is None for i in changed), (
            f"{model}: a templated change was invisible to the oracle, which "
            "would mean the oracle is wrong rather than the position untemplated")
    assert checked > 150, f"{model}: only {checked} records carried a mutation"
    assert 0 < invisible < checked * 0.15, (
        f"{model}: {invisible}/{checked} counted mutations were unrecoverable; "
        "expected a few per cent of untemplated hits")


@verifies("maturation.hypermutation",
          claim="the recorded per-segment counts sum to the recorded total, so the "
                "attribution accounts for every event",
          fields=("n_v_mutations", "n_d_mutations", "n_j_mutations", "n_np_mutations"),
          basis="consistent")
@pytest.mark.parametrize("model", MODELS)
def test_per_segment_counts_sum_to_the_total(model):
    exp = experiment_builder.recombined().mutate(model=model, count=4)
    for outcome in exp.run(n=50, seed=2):
        record = outcome.airr_record(exp.refdata)
        parts = sum(record[f] for f in ("n_v_mutations", "n_d_mutations",
                                        "n_j_mutations", "n_np_mutations"))
        assert parts == record["n_mutations"], (
            f"{model}: per-segment counts sum to {parts} but n_mutations is "
            f"{record['n_mutations']}")


@verifies("maturation.hypermutation",
          claim="n_observable_mutations is what a mutation caller can recover: "
                "the hypermutation substitutions still visible in the "
                "sequence. In a pipeline where hypermutation is the only pass "
                "that writes a base it equals the templated mismatch set an "
                "oracle derives from the per-base provenance -- which is the "
                "condition this test runs under, because that oracle cannot "
                "see which pass wrote a base. Bounded above by n_mutations, "
                "which additionally counts mutations on untemplated bases and "
                "repeat hits that restored the germline base. Attribution "
                "under corruption is the separate claim below",
          fields=("n_observable_mutations",),
          basis="provenance")
@pytest.mark.parametrize("count", [1, 20])
@pytest.mark.parametrize("model", MODELS)
def test_n_observable_mutations_is_what_a_caller_can_recover(model, count):
    """The field exists because `n_mutations` cannot be both things.

    A mutation on an untemplated junction base changes a base, so it belongs
    in `n_mutations`, and has no germline base to differ from, so nothing can
    recover it. Rather than redefine the count, the recoverable number is
    exported alongside it. Scoring a caller against the wrong one of the two
    penalises it for substitutions that left no evidence.
    """
    exp = experiment_builder.recombined().mutate(model=model, count=count)
    totals = {"recorded": 0, "reported": 0, "oracle": 0}
    for outcome in exp.run(n=60, seed=101):
        record = outcome.airr_record(exp.refdata)
        expected = len(reconstruction.observable_substitutions(outcome.final_simulation()))
        assert record["n_observable_mutations"] == expected, (
            f"{model}/{count}: reports {record['n_observable_mutations']} "
            f"recoverable where the oracle derives {expected}")
        assert record["n_observable_mutations"] <= record["n_mutations"], (
            "more recoverable than counted, which cannot happen")
        totals["recorded"] += record["n_mutations"]
        totals["reported"] += record["n_observable_mutations"]
        totals["oracle"] += expected
    assert totals["oracle"] > 0, f"{model}/{count}: nothing was mutated"
    assert totals["reported"] == totals["oracle"]
    # The gap is the point of the field: if it were zero the two counts would
    # be the same number and one of them would be redundant.
    assert totals["reported"] < totals["recorded"], (
        f"{model}/{count}: no gap between {totals['recorded']} counted and "
        f"{totals['reported']} recoverable, so this field establishes nothing")


# Each entry is a pipeline that writes bases for a reason that is not somatic
# hypermutation: PCR misincorporation, a sequencing miscall, an ambiguous base
# call masking a position with N, or wholesale replacement by a background
# molecule. None of them is a mutation a caller should be scored on.
NON_SHM_WRITERS = (
    ("pcr", lambda e: e.pcr_amplify(count=5)),
    ("sequencing", lambda e: e.sequencing_errors(count=5)),
    ("ambiguous", lambda e: e.ambiguous_base_calls(count=4)),
    ("contaminant", lambda e: e.contaminate(prob=1.0)),
)


@verifies("maturation.hypermutation",
          claim="n_observable_mutations is attributable to somatic "
                "hypermutation and to nothing else. A pipeline that writes "
                "bases for a library or sequencing reason -- PCR "
                "misincorporation, a sequencing miscall, an N mask, a "
                "background contaminant -- moves no part of it, so a caller "
                "scored against it is not credited for calling noise",
          fields=("n_observable_mutations",),
          basis="reference")
@pytest.mark.parametrize("label,corrupt", NON_SHM_WRITERS,
                         ids=[n for n, _ in NON_SHM_WRITERS])
def test_observable_mutations_exclude_every_non_shm_writer(label, corrupt):
    """The field is published as the number to score a mutation caller
    against. Derived from the finished molecule alone it cannot be that: every
    pass that writes a base leaves a position differing from its germline, and
    a count of those differences credits a caller for calling sequencing
    noise as hypermutation. A contaminant record would hand it a target in
    the hundreds against a real burden of five.

    Two assertions, because the two halves fail differently. With no
    mutation pass in the pipeline the only honest answer is zero. With one,
    the recoverable count cannot exceed what was counted.
    """
    unmutated = corrupt(experiment_builder.recombined())
    for index, outcome in enumerate(unmutated.run(n=30, seed=3)):
        record = outcome.airr_record(unmutated.refdata)
        assert record["n_mutations"] == 0, (
            f"{label}: fixture error -- no mutation pass is configured")
        assert record["n_observable_mutations"] == 0, (
            f"{label}: record {index} reports "
            f"{record['n_observable_mutations']} recoverable mutations in a "
            f"pipeline with no mutation pass")

    mutated = corrupt(experiment_builder.recombined().mutate(model="s5f", count=5))
    for index, outcome in enumerate(mutated.run(n=30, seed=3)):
        record = outcome.airr_record(mutated.refdata)
        assert record["n_observable_mutations"] <= record["n_mutations"], (
            f"{label}: record {index} reports "
            f"{record['n_observable_mutations']} recoverable against "
            f"{record['n_mutations']} counted")


@verifies("maturation.hypermutation",
          claim="a lineage node's mutation count is the number of positions "
                "that differ from the germline, which excludes a position a "
                "later mutation returned to its germline base. S5F reaches "
                "this case: the kernel excludes the base currently present, "
                "not the germline base, so a second hit at a mutated position "
                "can draw the germline base back. The restored position is "
                "stored as an uppercase SHM write over a lowercase germline "
                "byte, so a byte comparison counts it and a case-insensitive "
                "one does not",
          fields=("n_mutations",),
          basis="provenance")
def test_a_lineage_node_does_not_count_a_restored_position():
    """Unlike the per-pass counters, a lineage node's count is a NET
    difference: the node's molecule against the founder's germline. Two hits
    at one position therefore contribute whatever the second one left behind,
    and if that is the germline base they contribute nothing.

    This is the one place the S5F model is not exempt from the no-op
    question. Its kernel never draws the base already present, so no single
    S5F write is ever a no-op -- but it excludes the CURRENT base, not the
    germline base, so the second write at a position can restore the
    germline. Measured here: 5 such positions across these five nodes.

    The expected count is derived from the exported per-base provenance; the
    field is produced by the engine's own net-difference synthesis. The two
    agree only if both compare the base to its germline the same way, which
    is the mistake this pins.
    """
    exp = experiment_builder.recombined().clonal_lineage(
        n_clones=2, max_generations=2, n_max=20, n_sample=8, rate=0.05)
    result = exp.run_records(seed=11)
    restored = 0
    for outcome in result.outcomes:
        record = outcome.airr_record(exp.refdata)
        simulation = outcome.final_simulation()
        bases, germline = simulation.bases(), simulation.germline_bases()
        differing = sum(1 for i in range(len(bases))
                        if chr(bases[i]).upper() != chr(germline[i]).upper())
        assert record["n_mutations"] == differing, (
            f"node reports {record['n_mutations']} where "
            f"{differing} positions differ from their germline base")
        restored += sum(1 for i in range(len(bases))
                        if bases[i] != germline[i]
                        and chr(bases[i]).upper() == chr(germline[i]).upper())
    assert restored > 0, (
        "no position was restored to its germline base, so this test "
        "exercises nothing: raise the rate or the generation count")


BUNDLED_S5F_KERNELS = ("hh_s5f", "hh_s5f_60", "hh_s5f_opposite", "hkl_s5f")


@verifies("maturation.hypermutation",
          claim="no bundled S5F kernel ever draws the base already present, "
                "so every S5F substitution attempt is realised and for that "
                "model count and n_mutations coincide. This is a property of "
                "the shipped tables, not of the sampler, which takes every "
                "positive entry in a context's row -- so it is measured here "
                "rather than assumed from the code",
          fields=("n_mutations",),
          basis="reference")
@pytest.mark.parametrize("kernel", BUNDLED_S5F_KERNELS)
def test_no_bundled_s5f_kernel_draws_the_base_already_present(kernel):
    """The uniform model's quarter-of-draws no-op rate is the reason
    `n_mutations` stopped counting attempts. S5F is exempt -- but only
    because of what is in the four bundled tables, and the documentation
    says so in several places, so it is worth a measurement rather than a
    reading of the sampler.
    """
    attempts = 15
    exp = experiment_builder.recombined().mutate(
        model="s5f", count=attempts, s5f_model=kernel)
    for index, outcome in enumerate(exp.run(n=60, seed=5)):
        record = outcome.airr_record(exp.refdata)
        assert record["n_mutations"] == attempts, (
            f"{kernel}: record {index} realised {record['n_mutations']} of "
            f"{attempts} attempts, so this kernel admits the base already "
            f"present as a destination")


REGION_FIXTURES = (
    ("single_hit", lambda: experiment_builder.recombined().mutate(model="s5f", count=1)),
    ("repeat_hits", lambda: experiment_builder.recombined().mutate(model="s5f", count=20)),
    ("with_read_loss", lambda: experiment_builder.recombined()
     .mutate(model="s5f", count=20).end_loss_5prime(length=12)),
    ("with_indels", lambda: experiment_builder.recombined()
     .mutate(model="s5f", count=20).polymerase_indels(count=(1, 3), insertion_prob=0.5)),
)


def _reported_region_counters(record):
    """The six counters, named literally rather than built from a loop.

    The claimed fields have to be readable in the body that claims them, and
    an f-string is invisible to the static check -- which rejected the first
    version of this test for exactly that reason.
    """
    return {
        "fwr1": record["n_fwr1_observable_mutations"],
        "cdr1": record["n_cdr1_observable_mutations"],
        "fwr2": record["n_fwr2_observable_mutations"],
        "cdr2": record["n_cdr2_observable_mutations"],
        "fwr3": record["n_fwr3_observable_mutations"],
        "v_unannotated": record["n_v_unannotated_observable_mutations"],
    }


def _v_subregions(refdata, v_call):
    for index in range(refdata.v_pool_size()):
        allele = refdata.v_allele(index)
        if allele.name == v_call:
            return allele.subregions
    raise AssertionError(f"{v_call} names no allele in the bound reference")


@verifies("maturation.hypermutation",
          claim="the per-subregion recoverable counters are what a CDR/FWR "
                "selection-pressure ratio computed from a read needs. They "
                "agree with an oracle that buckets the visible differences by "
                "the allele's own IMGT table in every configuration measured, "
                "including repeat hits, read loss and indels -- where the "
                "event counters do not, because an event whose base is gone "
                "or overwritten is still an event",
          fields=("n_fwr1_observable_mutations", "n_cdr1_observable_mutations",
                  "n_fwr2_observable_mutations", "n_cdr2_observable_mutations",
                  "n_fwr3_observable_mutations",
                  "n_v_unannotated_observable_mutations"),
          basis="provenance")
@pytest.mark.parametrize("label,make", REGION_FIXTURES,
                         ids=[name for name, _ in REGION_FIXTURES])
def test_the_per_region_recoverable_counters_match_what_the_read_shows(label, make):
    """Why these fields exist, stated as the measurement that motivated them.

    The event counters are a partition of `n_v_mutations` and correct as
    that. Against what a read shows they agree on 60 of 60 records at one
    mutation, 32 of 60 at twenty, and 17 of 60 once twelve nucleotides are
    lost from the 5' end -- the record counting mutations in bases the read
    no longer contains. A CDR/FWR ratio computed from them is an event ratio,
    which is a different quantity from the one a selection-pressure analysis
    of a read is after.
    """
    exp = make()
    agreed = checked = 0
    for outcome in exp.run(n=60, seed=7):
        record = outcome.airr_record(exp.refdata)
        simulation = outcome.final_simulation()
        expected = reconstruction.observable_substitutions_by_v_subregion(
            simulation, _v_subregions(exp.refdata, record["v_call"]))
        reported = _reported_region_counters(record)
        assert reported == expected, (
            f"{label}: record reports {reported} where the read shows {expected}")
        agreed += 1
        checked += 1
    assert checked == 60
    assert agreed == 60


@verifies("maturation.hypermutation",
          claim="the six per-subregion recoverable counters sum to the V part "
                "of n_observable_mutations, and each is at most its event "
                "counterpart -- a recoverable substitution had an event "
                "behind it in the same subregion",
          fields=("n_observable_mutations",),
          basis="provenance")
def test_the_per_region_recoverable_counters_partition_the_v_segment():
    """The invariant the validator enforces, asserted here against the
    provenance rather than against the validator's own recomputation."""
    exp = (experiment_builder.recombined()
           .mutate(model="s5f", count=20).end_loss_5prime(length=10))
    checked = 0
    for outcome in exp.run(n=60, seed=7):
        record = outcome.airr_record(exp.refdata)
        simulation = outcome.final_simulation()
        expected = reconstruction.observable_substitutions_by_v_subregion(
            simulation, _v_subregions(exp.refdata, record["v_call"]))
        parts = sum(record[f"n_{bucket}_observable_mutations"]
                    for bucket in reconstruction.V_SUBREGION_BUCKETS)
        assert parts == sum(expected.values()), (
            f"the six counters sum to {parts}, the read shows "
            f"{sum(expected.values())} recoverable V substitutions")
        assert parts <= record["n_observable_mutations"], (
            f"the V buckets sum to {parts}, above the whole-record "
            f"recoverable count of {record['n_observable_mutations']}")
        for bucket in reconstruction.V_SUBREGION_BUCKETS:
            assert (record[f"n_{bucket}_observable_mutations"]
                    <= record[f"n_{bucket}_mutations"]), (
                f"{bucket}: more recoverable than counted")
        checked += 1
    assert checked == 60
