"""What the live population does as a lineage fills its germinal center.

None of this was measured anywhere. The whole clonal suite passed both before
and after the regulation law was replaced, because every test asked what a
lineage's records looked like and none asked how many cells were alive. The law
that decides that was wrong in three ways at once, and the one test that
mentioned carrying capacity asserted a bound on the peak -- which an
oscillating population satisfies, and so does a population that reached
capacity and died twenty-nine generations before the run ended.
"""
from collections import Counter

import pytest

import GenAIRR as ga
from GenAIRR import _engine
from GenAIRR._s5f_loader import load_builtin_s5f_kernel

pytestmark = pytest.mark.behavior


def _founder():
    compiled = ga.Experiment.on("human_igh").recombine().compile()
    return compiled.run(n=1, seed=0)[0].final_simulation()


@pytest.fixture(scope="module")
def grow():
    """Grow a neutral, mutation-free lineage and return cells per generation."""
    founder = _founder()
    mut, sub = load_builtin_s5f_kernel("hh_s5f")

    def _grow(*, lambda_base, n_max, generations=40, seed=7):
        tree = _engine.simulate_lineage(
            founder, mut, sub, 0.0, lambda_base, 0.0, generations, n_max, 1, seed)
        per_generation = Counter(node.generation for node in tree.nodes())
        return [per_generation.get(g, 0) for g in range(generations + 1)]

    return _grow


def test_a_lineage_that_reaches_capacity_is_still_alive_at_the_end(grow):
    """Reaching the carrying capacity used to be fatal with certainty.

    The hard cap truncated a generation to exactly ``n_max``, and the logistic
    law's effective rate at exactly ``n_max`` was exactly zero, so every cell
    in the next generation drew zero offspring. Across 180 engine runs, 67
    reached capacity and all 67 were extinct one generation later.
    """
    reached = died = 0
    for seed in range(1, 25):
        population = grow(lambda_base=4.0, n_max=50, seed=seed)
        if max(population) >= 50:
            reached += 1
            if population[-1] == 0:
                died += 1
    # Guard the premise: seeds that never fill the center say nothing about
    # what filling it costs. This is the assertion the original test lacked.
    assert reached >= 20, f"only {reached} of 24 seeds reached capacity"
    assert died == 0, f"{died} of {reached} lineages died after reaching capacity"


@pytest.mark.parametrize("lambda_base", [1.5, 2.0, 4.0])
def test_the_settled_population_is_the_capacity_itself(grow, lambda_base):
    """``n_max`` is the population a lineage settles at, for every rate.

    Under the discrete logistic law the settling point was
    ``n_max * (1 - 1/lambda_base)``, so the parameter named a population the
    lineage never reached: 125 of a requested 500 at ``lambda_base=1.5``, 250
    at 2.0, 375 at 4.0. Beverton-Holt's fixed point is ``n_max``, and because
    the law is monotone the fixed point does not move with the rate.
    """
    # Seed 22 is one of the seeds whose founder establishes at all three
    # rates. A founder dies out with probability q solving
    # q = exp(lambda_base(q - 1)) -- about 0.42 at 1.5 -- so a seed picked
    # without checking tests founder luck rather than regulation.
    n_max = 500
    population = grow(lambda_base=lambda_base, n_max=n_max, seed=22)
    settled = population[-1]
    # Guard the premise: a dead lineage settles at zero and would pass a
    # one-sided bound.
    assert settled > 0, "the lineage died, so there is no settled population"
    assert 0.8 * n_max < settled < 1.2 * n_max, (
        f"lambda_base={lambda_base} settled at {settled}, not at the capacity "
        f"{n_max}; the law it replaced settled at "
        f"{round(n_max * (1 - 1 / lambda_base))}")


def test_a_fast_growing_population_does_not_swing(grow):
    """Above ``lambda_base=3`` the old law had no settled population at all.

    The discrete logistic map is past its period-doubling threshold there: at
    3.5 it cycled between 383 and 875, and at 4.0 it moved chaotically between
    218 and 942 of a capacity of 1000, while the documentation promised a
    plateau. Beverton-Holt is monotone in the live population, so it cannot
    oscillate at any rate.
    """
    n_max = 500
    population = grow(lambda_base=4.0, n_max=n_max, generations=40, seed=22)
    tail = population[-10:]
    assert min(tail) > 0, "the lineage died before the window being measured"
    # The old law's tail spanned better than 4x. Poisson noise on a population
    # of 500 at replacement is a few percent per generation.
    assert max(tail) / min(tail) < 1.5, (
        f"the population is still swinging over the last ten generations: {tail}")


def test_the_extinction_error_does_not_advise_the_thing_that_caused_it():
    """The old message blamed the founder and advised raising ``lambda_base``.

    Under the old law that advice made the failure more likely: a faster
    growing population reaches the cap sooner, and reaching it was fatal. The
    message also asserted that each founder "drew 0 offspring", which was
    false -- they had grown to capacity first.
    """
    experiment = (ga.Experiment.on("human_igh").recombine()
                  .clonal_lineage(n_clones=1, max_generations=6, n_max=1000,
                                  n_sample=5, lambda_base=0.05)
                  .compile())
    with pytest.raises(ValueError) as caught:
        experiment.run_records(seed=3)
    message = str(caught.value)
    assert "went extinct on every one of" in message
    assert "allow_extinction=True" in message
    # It must not state a cause it has not established.
    assert "drew 0 offspring" not in message
