"""The named schemes a clonal run derives its sub-simulation seeds under.

The seeds in here are a wire value: a family-replay bundle records which
schedule drew it, so the first schedule's arithmetic has to stay exactly
what it was while the current one has to stop sharing streams between
runs.
"""

import pytest

from GenAIRR._seed_schedule import BY_LABEL, CURRENT, _Additive, _Keyed

pytestmark = pytest.mark.behavior


# The first schedule's arithmetic, as `origin/master` writes it at
# `_compiled/repertoire.py:106`, `_compiled/lineage.py:164,178,298` and
# `_compiled/repertoire.py:134`. Quoted here rather than read back off
# the schedule object: a replay has to reproduce the seeds the capture
# used, so this is the contract, not an implementation detail.
def _master_clone(run_seed, clone_index):
    return int(run_seed) + clone_index * 1_000_000


def _master_attempt(clone_seed, attempt_index):
    return clone_seed + attempt_index * 7_919


def _master_read_or_node(clone_seed, index):
    return clone_seed + 1 + index


class TestTheScheduleThatBundlesRecord:
    def test_its_label_is_the_one_written_into_existing_bundles(self):
        assert _Additive.label == "clone_1000000_attempt_7919_node_plus_one_v1"

    def test_it_derives_exactly_what_it_always_derived(self):
        v1 = _Additive()
        for run_seed in (0, 1, 7, 2**40):
            for clone_index in (0, 1, 2, 37):
                clone_seed = v1.clone(run_seed, clone_index)
                assert clone_seed == _master_clone(run_seed, clone_index)
                for index in (0, 1, 49, 7918):
                    assert v1.attempt(clone_seed, index) == _master_attempt(clone_seed, index)
                    assert v1.read(clone_seed, index) == _master_read_or_node(clone_seed, index)
                    assert v1.node(clone_seed, index) == _master_read_or_node(clone_seed, index)

    def test_its_defects_are_still_present_because_a_replay_needs_them(self):
        """Not a schedule anyone should draw under -- a schedule a bundle may
        have been drawn under. The three collisions are asserted so that
        the contrast below is a measurement rather than a claim."""
        v1 = _Additive()
        assert v1.clone(0, 1) == v1.clone(1_000_000, 0)
        assert v1.read(4242, 7918) == v1.attempt(4242, 1)
        assert v1.node(4242, 7918) == v1.read(4242, 7918)


class TestTheScheduleRunsDrawUnderNow:
    def test_it_is_the_current_one_and_says_so(self):
        assert CURRENT is BY_LABEL[_Keyed.label]
        assert _Keyed.label.endswith("_v2")

    def test_a_different_run_seed_is_a_different_lattice(self):
        seen = {}
        for run_seed in (0, 1, 2, 1_000_000, 2_000_000, 2**63 - 1):
            for clone_index in range(48):
                word = CURRENT.clone(run_seed, clone_index)
                assert word not in seen, f"{(run_seed, clone_index)} reuses {seen[word]}"
                seen[word] = (run_seed, clone_index)

    def test_the_axes_under_one_clone_stay_apart(self):
        clone_seed = CURRENT.clone(0, 0)
        seen = {}
        for axis in (CURRENT.attempt, CURRENT.read, CURRENT.node):
            for index in list(range(128)) + [7918, 7919, 15837, 999_999]:
                word = axis(clone_seed, index)
                assert word not in seen, f"{(axis.__name__, index)} reuses {seen[word]}"
                seen[word] = (axis.__name__, index)

    def test_every_seed_it_derives_is_one_the_engine_can_be_given(self):
        """A seed is handed straight to the engine, which takes a u64."""
        words = [CURRENT.clone(2**64 - 1, 2**32 - 1)]
        words.append(CURRENT.read(words[0], 2**32 - 1))
        words.append(CURRENT.node(words[0], 0))
        words.append(CURRENT.attempt(words[0], 49))
        assert all(0 <= word <= 2**64 - 1 for word in words)

    def test_an_index_it_cannot_name_is_refused_not_folded(self):
        with pytest.raises((OverflowError, ValueError)):
            CURRENT.clone(0, 2**32)


def test_a_recorded_label_names_exactly_one_schedule():
    assert set(BY_LABEL) == {_Additive.label, _Keyed.label}
    assert all(label == type(schedule).label for label, schedule in BY_LABEL.items())
