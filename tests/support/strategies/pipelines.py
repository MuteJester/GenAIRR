"""Draw a pipeline the compiler will accept.

The generator's job is to spend its examples on programs that run. A sweep
that draws uniformly over method names and keyword values spends most of its
budget on ValueError -- hypermutation on a TCR locus, a D inversion where
there is no D, a paired-end insert longer than the molecule -- and a campaign
measured in examples rather than in records is a campaign that mostly tested
the DSL's own argument checking, which is already tested.

So validity is built in rather than filtered out: the locus is drawn first,
only passes legal for it are offered, and the drawn subset is ordered by the
stage it declares. What remains invalid after that is a finding about the
registry, not noise -- `test_every_drawn_pipeline_builds_and_runs` is what
says so.
"""
from __future__ import annotations

from hypothesis import strategies as st

from tests.support.pipelines import Pipeline, for_locus
from tests.support.pipelines.registry import LOCI, STAGES


def loci() -> st.SearchStrategy[str]:
    return st.sampled_from(sorted(LOCI))


def _kwargs(one) -> st.SearchStrategy[dict]:
    """Draw a value for every parameter the pass declares.

    All of them, every time, rather than a drawn subset: a pass's defaults are
    exercised by the rest of the suite, and the sweep is here for the
    combinations nobody configured on purpose.
    """
    if not one.parameters:
        return st.just({})
    keys = sorted(one.parameters)
    return st.tuples(*(one.parameters[key] for key in keys)).map(
        lambda values: dict(zip(keys, values)))


@st.composite
def pipelines(draw, *, min_passes: int = 1, max_passes: int = 5,
              locus: str | None = None) -> Pipeline:
    """A locus and an ordered, legal stack of passes with drawn parameters.

    ``max_passes`` is five because that is roughly what a real library-prep
    pipeline looks like and because the battery's cost is per record, not per
    pass; deeper stacks are reachable by raising it in a campaign that has the
    time.
    """
    chosen_locus = locus or draw(loci())
    available = for_locus(chosen_locus)
    names = draw(st.lists(st.sampled_from([one.name for one in available]),
                          min_size=min_passes, max_size=max_passes,
                          unique=True))
    by_stage = {one.name: STAGES.index(one.stage) for one in available}
    ordered = sorted(names, key=lambda name: (by_stage[name], name))
    # Mutually exclusive passes are dropped rather than rejected. Rejecting
    # the example would spend the budget on the DSL's argument checking, which
    # is tested elsewhere; dropping keeps the rest of the stack, which is what
    # the campaign is for. The first kept pass wins, so the choice is a
    # function of the draw and not of dict order.
    offered_all = {one.name: one for one in available}
    kept: list[str] = []
    for name in ordered:
        if any(name in offered_all[other].excludes for other in kept):
            continue
        kept.append(name)
    ordered = kept
    steps = []
    offered = {one.name: one for one in available}
    for name in ordered:
        steps.append((name, draw(_kwargs(offered[name]))))
    return Pipeline(locus=chosen_locus, steps=tuple(steps))


def pipelines_with(name: str, **kwargs) -> st.SearchStrategy[Pipeline]:
    """Pipelines that definitely contain ``name``.

    For a campaign aimed at one pass's interactions rather than at the space
    as a whole -- the shape every finding so far would have been caught by.
    """
    return pipelines(**kwargs).filter(lambda pipeline: pipeline.has(name))
