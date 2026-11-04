"""Two run seeds draw two repertoires, not one repertoire twice.

Clone seeds used to be ``run_seed + clone_index * 1_000_000``, which
made clone 1 of a run at seed 0 the same simulation as clone 0 of a run
at seed 1,000,000 -- same founder, and with the read offsets agreeing
too, the same corrupted read. Every clonal DSL derived its seeds that
way, so this is asserted through each of them rather than through the
schedule object alone: the schedule is only correct if the DSLs use it.
"""

from collections import defaultdict

import pytest

import GenAIRR as ga
from GenAIRR._seed_schedule import CURRENT, _Additive

pytestmark = pytest.mark.behavior

#: Clone 1 under the old schedule took exactly the seed clone 0 took at
#: this run seed. Everything below is that one coincidence, observed.
TRANSLATED_RUN_SEED = 1_000_000


def _by_clone(result, *fields):
    grouped = defaultdict(list)
    for record in result.records:
        grouped[record["clone_id"]].append(tuple(record[field] for field in fields))
    return grouped


def test_the_old_schedule_translated_and_the_current_one_does_not():
    """The arithmetic the observations below are observations of."""
    old = _Additive()
    assert old.clone(0, 1) == old.clone(TRANSLATED_RUN_SEED, 0)
    assert CURRENT.clone(0, 1) != CURRENT.clone(TRANSLATED_RUN_SEED, 0)


def test_a_repertoire_clone_is_not_another_runs_clone():
    def run(seed, n_clones):
        return (ga.Experiment.on("human_igh").recombine()
                .clonal_repertoire(n_clones=n_clones, max_size=5)
                .run_records(seed=seed))

    a = _by_clone(run(0, 3), "sequence", "v_call", "d_call", "j_call", "junction")
    b = _by_clone(run(TRANSLATED_RUN_SEED, 1), "sequence", "v_call", "d_call", "j_call", "junction")
    assert set(a) == {0, 1, 2} and set(b) == {0}, "premise: both runs produced their clones"
    assert a[1] != b[0]
    assert a[0] != b[0]


def test_a_repertoire_read_is_not_another_runs_read():
    """With one read per clone the founder *and* the read seed coincided,
    so whole corrupted reads repeated across runs."""
    def run(seed, n_clones):
        return (ga.Experiment.on("human_igh").recombine()
                .clonal_repertoire(n_clones=n_clones, max_size=1)
                .sequencing_errors(rate=0.02)
                .run_records(seed=seed))

    a = _by_clone(run(0, 3), "sequence")
    b = _by_clone(run(TRANSLATED_RUN_SEED, 1), "sequence")
    assert [len(reads) for reads in (a[0], a[1], b[0])] == [1, 1, 1], "premise: one read each"
    assert a[1] != b[0]


def test_an_expanded_clones_descendant_is_not_another_runs_descendant():
    def run(seed, n_clones):
        return (ga.Experiment.on("human_igh").recombine()
                .expand_clones(n_clones=n_clones, per_clone=2)
                .sequencing_errors(rate=0.02)
                .run_records(seed=seed))

    a = _by_clone(run(0, 3), "sequence")
    b = _by_clone(run(TRANSLATED_RUN_SEED, 1), "sequence")
    assert [len(d) for d in (a[0], a[1], b[0])] == [2, 2, 2], "premise: two descendants each"
    assert a[1] != b[0]


def test_a_lineage_family_is_not_another_runs_family():
    """The founder, every survival retry and every observed node all hung
    off the same clone seed, so a whole family repeated."""
    def run(seed, n_clones):
        return (ga.Experiment.on("human_igh").recombine()
                .clonal_lineage(n_clones=n_clones, max_generations=4, n_max=100,
                                n_sample=10, rate=0.03, lambda_base=1.5)
                .run_records(seed=seed))

    a = _by_clone(run(0, 3), "sequence")
    b = _by_clone(run(TRANSLATED_RUN_SEED, 1), "sequence")
    assert len(a[1]) > 1 and len(b[0]) > 1, "premise: both families grew past their founder"
    assert a[1] != b[0]
