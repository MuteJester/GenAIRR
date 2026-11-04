"""Can a dataset -- or a trace file -- say which sampling model produced it?

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

**Two artifacts carry a version, for two different reasons.** A record
carries it so a reader can tell what drew the rows in front of them. A trace
file carries it so the engine can refuse to rerun something it would not
reproduce: a trace records a seed and a plan signature, which between them say
what to draw and from where, and neither says which model does the drawing.
Both halves are here because they are one claim -- the model that drew an
artifact travels with the artifact -- and splitting them across modules is how
one of them quietly stops being true.
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


def test_a_clonal_dataset_says_which_model_drew_it_too():
    """Version 2 moved the clonal seed schedule as well as the passes.

    A clone, a survival retry and a read used to take seeds offset from
    the run seed, which made clone 1 of one run the same simulation as
    clone 0 of another; they are keyed substreams now. That is the same
    kind of change as the pass substreams and carries the same number,
    so a clonal dataset has to carry it too -- and clonal records leave
    through their own export route, with `clone_id` and
    `duplicate_count`, rather than the plain one asserted above.
    """
    result = (ga.Experiment.on("human_igh").recombine()
              .clonal_repertoire(n_clones=4, max_size=5).run_records(seed=3))
    assert result.records, "premise: the repertoire produced records"
    assert {record["clone_id"] for record in result.records} == {0, 1, 2, 3}
    assert all(record[COLUMN] == ga.DETERMINISM_VERSION for record in result.records)


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


# ──────────────────────────────────────────────────────────────────
# The other artifact a version has to travel on
# ──────────────────────────────────────────────────────────────────


@pytest.fixture
def vj_compiled():
    return ga.Experiment.on("human_igk").recombine().compile()


def _emit(compiled, path, seed=4242):
    """Write a trace file and return `(path, its json)`."""
    import json

    outcome = compiled.simulator.run(seed=seed)
    compiled.simulator.trace_file_from(outcome, seed=seed).write_to(str(path))
    return path, json.loads(path.read_text(encoding="utf-8"))


def test_a_trace_file_says_which_model_drew_it(vj_compiled, tmp_path):
    """On the artifact, not inferred from the engine version beside it.

    `engine_version` is already in the file and does not answer this: two
    releases can ship the same sampling model, and a determinism change is not
    tied to a version bump.
    """
    from GenAIRR import _engine as genairr_engine

    path, body = _emit(vj_compiled, tmp_path / "t.trace.json")
    assert body[COLUMN] == ga.DETERMINISM_VERSION
    assert genairr_engine.TraceFile.read_from(str(path)).determinism_version == (
        ga.DETERMINISM_VERSION)


def test_a_trace_drawn_by_another_model_is_refused_rather_than_rerun(
        vj_compiled, tmp_path):
    """The whole point of the field: the rerun fails instead of succeeding
    with something else.

    The failure has to name the mismatch. A rerun that diverges silently is
    the defect; a rerun that fails with an unrelated message is the defect
    with a worse diagnosis.
    """
    import json

    from GenAIRR import _engine as genairr_engine

    path, body = _emit(vj_compiled, tmp_path / "t.trace.json")
    body[COLUMN] = ga.DETERMINISM_VERSION + 1
    path.write_text(json.dumps(body), encoding="utf-8")

    loaded = genairr_engine.TraceFile.read_from(str(path))
    assert loaded.determinism_version == ga.DETERMINISM_VERSION + 1, (
        "the doctored file must load, or this tests the loader instead")
    with pytest.raises(ValueError) as caught:
        vj_compiled.simulator.rerun_from_trace_file(loaded)
    assert "determinism" in str(caught.value).lower()


def test_a_trace_written_before_the_stamp_still_reruns(vj_compiled, tmp_path):
    """The one case where the gate yields, and why.

    A file with no such field cannot say which model drew it, so the engine
    cannot compare. Refusing those would break every trace written before this
    version, including the committed v1 fixtures whose whole purpose is to
    prove the on-disk format still loads. They keep the treatment they have
    always had: the signature gates and nothing more.
    """
    import json

    from GenAIRR import _engine as genairr_engine

    path, body = _emit(vj_compiled, tmp_path / "t.trace.json")
    del body[COLUMN]
    path.write_text(json.dumps(body), encoding="utf-8")

    loaded = genairr_engine.TraceFile.read_from(str(path))
    assert loaded.determinism_version == 0, (
        "an absent field reads as 0, which means unknown, not version zero")
    vj_compiled.simulator.rerun_from_trace_file(loaded)
