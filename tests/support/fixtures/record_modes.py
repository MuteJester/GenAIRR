"""One record from every mode GenAIRR can run in.

The catalog's field-existence check needs to know what fields exist, and a
plain record does not carry all of them: clonal runs add lineage truth, and a
hand-maintained exemption list for the difference is exactly the kind of thing
that rots -- the first version of it named a field (`generation`) that exists
nowhere in GenAIRR, in the one place the check stops looking.
"""
from __future__ import annotations

import functools

import GenAIRR as ga

FULL_STACK = dict(model="uniform", count=2)


def _plain():
    exp = (ga.Experiment.on("human_igh").recombine().mutate(**FULL_STACK)
           .polymerase_indels(count=1).sequencing_errors(count=1)
           .paired_end(r1_length=150, insert_size=300))
    return set(exp.run(n=1, seed=1)[0].airr_record(exp.refdata))


def _from_records(exp):
    rows = list(exp.compile().run_records(seed=1))
    return set(rows[0])


@functools.lru_cache(maxsize=1)
def fields_by_mode() -> dict[str, frozenset[str]]:
    """Mode name -> the record fields that mode produces."""
    return {
        "plain": frozenset(_plain()),
        "clonal_lineage": frozenset(_from_records(
            ga.Experiment.on("human_igh").recombine()
            .clonal_lineage(n_clones=1, max_generations=3, n_sample=5))),
        "clonal_repertoire": frozenset(_from_records(
            ga.Experiment.on("human_igh").recombine()
            .clonal_repertoire(n_clones=2, max_size=4))),
    }


@functools.lru_cache(maxsize=1)
def every_known_field() -> frozenset[str]:
    """Every field any mode produces, measured rather than listed."""
    return frozenset().union(*fields_by_mode().values())
