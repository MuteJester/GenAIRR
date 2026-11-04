"""Where a subject's seed and a genotype attempt's seed come from.

Neither of these is the collision the clonal schedule had: Python's
Mersenne Twister mixes, so `random.Random(seed).getrandbits(63)` did not
make subject 1 of one cohort into subject 0 of another. What it did was
put a second generator in the derivation chain -- one the engine's
`rng::discipline` guard cannot see, one whose stream depends on CPython's
`random` module rather than on the engine's determinism version, and one
that means "how GenAIRR turns a seed into a sub-simulation's seed" has
two answers.

So these two axes derive through the same keyed primitive as a clone, a
retry and a read, and the seeds are pinned here against that primitive
rather than against a recomputation of the old chain.
"""

import pytest

import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR._seed_schedule import sub_seed
from GenAIRR.genotype import Genotype

pytestmark = pytest.mark.behavior


def _cfg():
    return gdata.HUMAN_IGH_OGRDB


class TestACohortSubjectSeed:
    def test_it_is_a_named_substream_of_the_run_seed(self):
        cfg = _cfg()
        gts = [Genotype.sample(cfg, seed=s, subject_id=f"S{s}") for s in (0, 1, 2)]
        cohort = ga.Experiment.on(cfg).recombine().run_cohort(
            gts, n_per_subject=1, seed=404)
        assert len(cohort.subjects) == 3, "premise: the cohort ran every subject"
        assert [s.seed for s in cohort.subjects] == [
            sub_seed(404, "subject", 0),
            sub_seed(404, "subject", 1),
            sub_seed(404, "subject", 2),
        ]

    def test_the_run_seed_is_mixed_in_rather_than_handed_to_a_second_generator(self):
        """The property that makes this one mechanism and not two: the
        seed a subject gets is the one the engine's primitive derives,
        not the one CPython's `random` module happens to produce."""
        import random

        superseded = random.Random(404)
        would_have_been = [superseded.getrandbits(63) for _ in range(3)]
        now = [sub_seed(404, "subject", i) for i in range(3)]
        assert would_have_been != now, (
            "premise: the two chains disagree, or this test cannot tell them apart")
        assert len(set(now)) == 3

    def test_two_cohorts_at_different_run_seeds_share_no_subject(self):
        assert not (
            {sub_seed(0, "subject", i) for i in range(64)}
            & {sub_seed(1, "subject", i) for i in range(64)}
        )


class TestAGenotypeAttemptSeed:
    @staticmethod
    def _attempt_seeds(monkeypatch, **kwargs):
        """Every seed `sample()` hands to a draw, in order."""
        seen = []
        original = Genotype._draw_one.__func__

        def recording(cls, cfg, seed, *rest):
            seen.append(seed)
            return original(cls, cfg, seed, *rest)

        monkeypatch.setattr(Genotype, "_draw_one", classmethod(recording))
        Genotype.sample(_cfg(), **kwargs)
        return seen

    def test_it_is_a_named_substream_of_the_requested_seed(self, monkeypatch):
        seen = self._attempt_seeds(monkeypatch, seed=7, subject_id="A")
        assert seen, "premise: the recording wrapper saw a draw at all"
        assert seen[0] == sub_seed(7, "genotype_attempt", 0)

    def test_a_resample_cannot_land_on_a_neighbouring_draw(self, monkeypatch):
        """What the base RNG was there for, kept.

        A failed attempt at `seed` must not re-draw what a direct draw at
        `seed + 1` would, or raising `max_resamples` would quietly turn
        one caller's genotype into another's.
        """
        assert sub_seed(7, "genotype_attempt", 1) != sub_seed(8, "genotype_attempt", 0)
        assert sub_seed(7, "genotype_attempt", 1) != sub_seed(7, "genotype_attempt", 0)
        seen = self._attempt_seeds(monkeypatch, seed=7, subject_id="A")
        assert seen[0] != sub_seed(8, "genotype_attempt", 0)

    def test_a_subject_axis_and_an_attempt_axis_do_not_meet(self):
        """Both axes now start from a plain integer the caller chose, so
        they would collide if they shared a role."""
        assert not (
            {sub_seed(5, "subject", i) for i in range(64)}
            & {sub_seed(5, "genotype_attempt", i) for i in range(64)}
        )
