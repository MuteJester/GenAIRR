"""Deprecation contract for ``Experiment.expand_clones``.

Verifies that:

1. Calling ``expand_clones()`` emits a ``DeprecationWarning`` that
   mentions ``clonal_lineage``.
2. Despite the warning, ``expand_clones()`` continues to produce the
   correct star-topology clonal output (n_clones * per_clone records
   with the expected ``clone_id`` values).
"""

import warnings

import pytest

import GenAIRR as ga


def test_expand_clones_emits_deprecation_warning_but_still_works():
    exp = ga.Experiment.on("human_igh").recombine()
    with pytest.warns(DeprecationWarning, match="clonal_lineage"):
        exp = exp.expand_clones(n_clones=2, per_clone=3)
    # Pipeline continues to work after the deprecation warning.
    exp = exp.mutate(count=5)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        result = exp.run_records(seed=0)
    # 2 clones × 3 descendants = 6 records.
    assert len(result) == 6
    # Exactly the expected clone IDs are present.
    clone_ids = {r["clone_id"] for r in result}
    assert clone_ids == {0, 1}
