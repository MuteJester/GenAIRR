"""Can a dataset say which sampling model produced it?

`test_reproducibility` asserts that one seed gives one answer. This module
covers what happens when that answer changes, which it has now done once: the
same seed, on the same platform, through the same API, produces a different
repertoire before and after per-pass random substreams.

Without a version stamp the two are indistinguishable. Someone holding a TSV
has no way to tell which engine drew it, a benchmark comparing a prediction to
a dataset cannot tell a model change from a regression, and "reproduce it with
seed 42" is not a reproducible instruction. So the number travels with the
data.

**On the record rather than only beside it.** A GenAIRR dataset leaves this
library as a file of rows and is then passed on by itself -- attached to a
paper, dropped in a bucket, handed to a collaborator. A sidecar is the part of
that journey that gets lost. One constant column per row is the price of a file
that still answers the question after it has been separated from everything
that was shipped with it.

A dataset with no such column predates the stamp, which is why the numbering
starts at 2: numbering the current model 1 would make 1 mean both "the first
versioned model" and "unknown".
"""
from __future__ import annotations

import csv

import pytest

import GenAIRR as ga
from tests.support import golden

pytestmark = pytest.mark.composition

COLUMN = "determinism_version"


def _records(n: int = 3):
    from tests.support.builders import experiment as experiment_builder

    return (experiment_builder.recombined(locus="human_igh")
            .compile().run_records(n=n, seed=7).records)


@pytest.mark.parametrize("name,locus,n,variant", golden.every_fixture())
def test_the_exported_file_carries_the_version(name, locus, n, variant, tmp_path):
    """Read back through the export path, because that is the artifact.

    Parametrized over every reproducibility fixture rather than one
    configuration: the column has to survive each locus and each feature
    combination the goldens cover, and a stamp applied in one projection path
    and not another would be worse than no stamp at all.
    """
    path = golden.write_reproducibility_tsv(
        tmp_path / f"{locus}.tsv", locus, n=n, variant=variant)
    with open(path, encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows, f"{name} exported no rows"
    assert all(row[COLUMN] == str(ga.DETERMINISM_VERSION) for row in rows), (
        f"{name}: a row reports a version other than the engine's "
        f"{ga.DETERMINISM_VERSION}")


def test_the_version_is_readable_without_running_anything():
    """The question "which model is installed?" must be answerable before a
    run, not only from a run's output -- that is how a script decides whether
    the dataset it is about to compare against is comparable at all."""
    assert isinstance(ga.DETERMINISM_VERSION, int)
    assert ga.DETERMINISM_VERSION >= 2, (
        "1 is reserved for the unnumbered model that predates the stamp")


def test_every_record_reports_the_engine_that_drew_it():
    """On the record dict, not only in the file, so the field is present for
    the in-memory and streaming consumers too."""
    records = _records()
    assert all(record[COLUMN] == ga.DETERMINISM_VERSION for record in records)


def test_the_version_is_a_column_the_exporter_knows_about():
    """Ordered deliberately rather than appended as an unknown extra.

    A field that only appears because the exporter tacks unrecognised keys on
    the end would move whenever another extra is added, and `with_metadata`
    would not know the name is taken.
    """
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    assert COLUMN in _DEFAULT_COLUMN_ORDER


def test_a_user_cannot_overwrite_the_version_with_metadata():
    """The stamp is the engine's claim about itself. A dataset whose version
    column says whatever its author wanted answers nothing."""
    from tests.support.builders import experiment as experiment_builder

    with pytest.raises(ValueError):
        (experiment_builder.recombined(locus="human_igh")
         .with_metadata(**{COLUMN: 1})
         .compile().run_records(n=1, seed=7))
