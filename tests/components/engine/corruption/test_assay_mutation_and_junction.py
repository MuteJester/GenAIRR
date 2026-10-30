"""Independent regressions for read-length mutation bias and missing junction AA."""

import pytest
import GenAIRR as ga
from GenAIRR import _engine

from tests.support.cases.productive import translate
from tests.support.fixtures.minimal_refdata import _vj_basic

pytestmark = pytest.mark.behavior


def experiment(np_length=3):
    return (
        ga.Experiment.on(_vj_basic())
        .recombine(np1_lengths=[(np_length, 1.0)]).trim(enabled=False)
    )


@pytest.mark.parametrize("model", ["uniform", "s5f"])
@pytest.mark.parametrize("loss", [0, 6, 100])
def test_mutation_rate_uses_stage_molecule_length_even_after_complete_loss(model, loss):
    exp = experiment().mutate(model=model, count=3).end_loss_3prime(length=loss)
    result = exp.run_records(n=3, seed=5)
    for record in result:
        # The subject is the DENOMINATOR: the rate is over the stage molecule
        # (18), not the shorter read left after loss. Asserting it against the
        # record's own count states that and nothing else -- the count itself
        # is bounded below, because `count=3` is three attempts and a uniform
        # draw that lands on the base already there is no longer counted.
        assert 0 < record["n_mutations"] <= 3
        assert record["sequence_length"] == 18 - min(loss, 18)
        assert record["mutation_rate"] == pytest.approx(record["n_mutations"] / 18)
    assert result.validate_records(exp.refdata).ok


def test_repeated_mutation_stages_preserve_cumulative_burden_and_zero_pass_invariance():
    exp = experiment().mutate(model="uniform", count=2).mutate(count=3).mutate(count=0)
    exp.end_loss_5prime(length=6)
    result = exp.run_records(n=4, seed=9)
    for record in result:
        # Three stages of 2 + 3 + 0 attempts; the zero-count pass must not
        # perturb the cumulative burden, which is the invariance under test.
        assert 0 < record["n_mutations"] <= 5
        assert record["mutation_rate"] == pytest.approx(record["n_mutations"] / 18)
    assert result.validate_records(exp.refdata).ok


@pytest.mark.parametrize("post_mutations", [0, 2])
def test_clones_keep_parent_mutation_counters_and_rate_through_observation_loss(post_mutations):
    from GenAIRR._compiled.clonal import CompiledClonalExperiment
    # The public DSL intentionally places SHM after the fork. Exercise the
    # lower-level continuation API with an already mutated founder.
    parent = experiment().mutate(model="uniform", count=3).compile()
    template = (
        experiment().expand_clones(n_clones=3, per_clone=2)
        .mutate(model="uniform", count=post_mutations).end_loss_3prime(length=6).compile()
    )
    exp = CompiledClonalExperiment(parent.simulator, template._post, template._fork, parent.refdata)
    result = exp.run_records(seed=4)
    counts = set()
    for record in result:
        # `count=N` is N attempts; a uniform draw landing on the base already
        # there is not counted, so the realised total is at most the attempts.
        assert record["n_mutations"] <= 3 + post_mutations
        # The subject: the rate's denominator is the stage molecule (18), not
        # the shorter read left after the 3' loss.
        assert record["mutation_rate"] == pytest.approx(record["n_mutations"] / 18)
        assert sum(record[f"n_{segment}_mutations"] for segment in ("v", "d", "j", "np")) == record["n_mutations"]
        counts.add(record["n_mutations"])
    if post_mutations == 0:
        # With no post-fork pass every record of one clone carries its
        # founder's counters unchanged. Grouped by clone, because the three
        # founders mutate independently and their realised counts differ.
        by_clone = {}
        for record in result:
            by_clone.setdefault(record["parent_id"], set()).add(record["n_mutations"])
        for parent, group in sorted(by_clone.items()):
            assert len(group) == 1, (
                f"clone {parent} diverged from its founder: {sorted(group)}")
    assert result.validate_records(exp.refdata).ok


def test_repertoire_preserves_founder_burden_with_and_without_postfork_pass():
    """A zero-count post-fork pass must not perturb the founders' burden.

    Both variants are built here rather than parametrized, because that is
    the claim the name makes and a parametrized run can only ever see one of
    them. The previous version asserted every record reported exactly 3
    mutations, which held only while the engine counted substitutions that
    changed nothing: `count=3` is three attempts, and the three founders now
    realise 1, 2 and 2. A repertoire collapses each clone to one record with
    a `duplicate_count`, so there is nothing to compare within a clone -- the
    comparison that matters is between the two variants.
    """
    from GenAIRR._compiled.repertoire import CompiledRepertoireExperiment

    def burdens(post_pass):
        parent = experiment().mutate(model="uniform", count=3).compile()
        template = experiment().clonal_repertoire(n_clones=3, max_size=2)
        if post_pass:
            template.mutate(count=0).end_loss_3prime(length=6)
        template = template.compile()
        exp = CompiledRepertoireExperiment(
            parent.simulator, template._post, template._step, parent.refdata)
        result = exp.run_records(seed=9)
        assert result.validate_records(exp.refdata).ok
        for record in result:
            assert record["n_mutations"] <= 3
            # The denominator is the stage molecule (18), not the shorter
            # read the 3' loss leaves in the post-pass variant.
            assert record["mutation_rate"] == pytest.approx(record["n_mutations"] / 18)
        return {record["clone_id"]: record["n_mutations"] for record in result}

    without, with_pass = burdens(False), burdens(True)
    assert without == with_pass, (
        f"the zero-count post-fork pass changed the founders' burden: "
        f"{without} -> {with_pass}")
    assert len(without) == 3 and any(without.values())


def test_lineage_net_difference_rate_is_not_inflated_by_read_loss():
    exp = (
        ga.Experiment.on("human_igh").recombine()
        .clonal_lineage(n_clones=1, max_generations=3, n_max=40, n_sample=8, rate=0.03)
        .end_loss_3prime(length=30)
    )
    result = exp.run_records(seed=9)
    assert result.records
    for record in result:
        assert record["mutation_rate"] == pytest.approx(
            record["n_mutations"] / (record["sequence_length"] + record["end_loss_3_length"])
        )
    assert result.validate_records(exp.refdata).ok


def test_lineage_corruption_trace_spans_follow_founder_choices_without_duplicate_shm():
    def lineage():
        return (
            ga.Experiment.on("human_igh").recombine()
            .clonal_lineage(n_clones=1, max_generations=3, n_max=40, n_sample=8, rate=0.03)
        )

    pristine_exp = lineage()
    corrupted_exp = (
        lineage().polymerase_indels(count=2, insertion_prob=0.5)
        .end_loss_3prime(length=30)
    )
    pristine = pristine_exp.run_records(seed=9)
    corrupted = corrupted_exp.run_records(seed=9)
    assert len(pristine) == len(corrupted) > 0
    for before, after, base, merged in zip(
        pristine, corrupted, pristine.outcomes, corrupted.outcomes
    ):
        assert before["n_mutations"] == after["n_mutations"]
        assert before["mutation_rate"] == after["mutation_rate"]
        assert after["mutation_rate"] == pytest.approx(
            before["n_mutations"] / before["sequence_length"]
        )
        choices = merged.trace().choices()
        prefix = base.trace().choices()
        assert [(c.address, c.value) for c in choices[:len(prefix)]] == [
            (c.address, c.value) for c in prefix
        ]
        events = merged.events()
        assert sum(event.pass_name.startswith("mutate.") for event in events) == 1
        assert events[0].trace_start == events[0].trace_end == 0
        cursor = len(prefix)
        for event in events[1:]:
            assert event.trace_start == cursor
            assert event.trace_end > cursor
            assert all(
                choice.address.startswith(event.pass_name)
                for choice in choices[event.trace_start:event.trace_end]
            )
            cursor = event.trace_end
        assert cursor == len(choices)
    assert corrupted.validate_records(corrupted_exp.refdata).ok


def test_merged_clonal_history_replays_against_equivalent_complete_plan():
    split = experiment().expand_clones(n_clones=1, per_clone=1).mutate(model="uniform", count=3)
    original = split.run_records(seed=7).outcomes[0]
    complete = experiment().mutate(model="uniform", count=3).compile()
    trace = complete.simulator.trace_file_from(original, seed=7)
    replayed = complete.simulator.replay_from_trace_file(trace)
    assert original.airr_record(complete.refdata) == replayed.airr_record(complete.refdata)


@pytest.mark.parametrize("np_length", [1, 2])
def test_nonproductive_physical_junction_keeps_complete_codon_translation(np_length):
    exp = experiment(np_length=np_length)
    result = exp.run_records(n=20, seed=9)
    for record in result:
        assert record["vj_in_frame"] is False
        assert record["productive"] is False
        assert record["junction_aa"] == translate(record["junction"])
        assert len(record["junction_aa"]) == len(record["junction"]) // 3
    assert result.validate_records(exp.refdata).ok


def test_nonproductive_junction_preserves_stop_translation():
    ref = _engine.RefDataConfig.vj()
    ref.add_v_allele("v1*01", "v1", b"AAACCCTGTTAA", anchor=6)
    ref.add_j_allele("j1*01", "j1", b"TGGAAA", anchor=0)
    exp = ga.Experiment.on(ref).recombine(np1_lengths=[(0, 1)]).trim(enabled=False)
    result = exp.run_records(n=1, seed=0)
    assert result[0]["junction_aa"] == "C*W"
    assert result[0]["productive"] is False
    assert result[0]["stop_codon"] is True
    assert result.validate_records(exp.refdata).ok


def test_zero_mutations_and_zero_requested_records_are_well_defined():
    exp = experiment().mutate(count=0).end_loss_3prime(length=100)
    result = exp.run_records(n=1, seed=0)
    assert result[0]["n_mutations"] == 0
    assert result[0]["mutation_rate"] == 0.0
    with pytest.raises(ValueError, match="n must be at least 1"):
        exp.run_records(n=0, seed=0)
