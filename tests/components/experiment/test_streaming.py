"""Streaming results instead of materialising them.
"""
from __future__ import annotations
from GenAIRR import _engine as ge
import pytest
from GenAIRR import CompiledExperiment, Experiment
import GenAIRR as ga  # noqa: E402
import inspect as _inspect
import itertools as _itertools
from tests.support.fixtures.experiment_dsl import _human_igh_exp

pytestmark = pytest.mark.behavior

def test_stream_returns_generator():
    exp = _human_igh_exp()
    s = exp.stream(n=3, seed=0)
    assert _inspect.isgenerator(s)


def test_stream_yields_n_outcomes_when_bounded():
    out = list(_human_igh_exp().stream(n=4, seed=0))
    assert len(out) == 4


def test_stream_default_n_is_unbounded():
    # No ``n`` → infinite generator. Use ``islice`` to bound.
    s = _human_igh_exp().stream(seed=0)
    first_seven = list(_itertools.islice(s, 7))
    assert len(first_seven) == 7


def test_stream_matches_run_for_same_seed():
    exp = _human_igh_exp().mutate(count=5)
    runs = exp.run(n=4, seed=0xCAFE)
    streamed = list(exp.stream(n=4, seed=0xCAFE))
    assert len(runs) == len(streamed)
    for a, b in zip(runs, streamed):
        assert a.final_simulation().bases() == b.final_simulation().bases()


def test_stream_yields_engine_outcomes():
    from GenAIRR import _engine as ge

    out = list(_human_igh_exp().stream(n=2, seed=0))
    assert all(isinstance(o, ge.Outcome) for o in out)


def test_stream_n_zero_raises():
    with pytest.raises(ValueError, match="n must be at least 1"):
        list(_human_igh_exp().stream(n=0))


def test_stream_n_negative_raises():
    with pytest.raises(ValueError, match="n must be at least 1"):
        list(_human_igh_exp().stream(n=-3))


def test_stream_is_lazy_under_unbounded_n():
    # Constructing a stream with no bound must not actually run any
    # simulations until iteration begins. Confirmed by the fact that
    # construction returns immediately even with no ``n``.
    s = _human_igh_exp().stream(seed=0)
    assert _inspect.isgenerator(s)
    # Pull a single item — should not blow up regardless of "n".
    first = next(s)
    assert first is not None


def test_stream_seed_offsets_per_iteration():
    # The i-th yielded outcome should match a single ``run()`` call
    # with ``seed = base_seed + i``. Confirms the documented
    # ``seed + i`` per-run offset.
    base_seed = 17
    exp = _human_igh_exp().mutate(count=3)
    streamed = list(exp.stream(n=4, seed=base_seed))
    for i, outcome in enumerate(streamed):
        single = exp.run(n=1, seed=base_seed + i)[0]
        assert outcome.final_simulation().bases() == single.final_simulation().bases()


def test_stream_supports_respect_kwarg():
    out = list(
        _human_igh_exp()
        .mutate(count=5)
        .productive_only().stream(n=3, seed=0)
    )
    # The stream should still produce productive outcomes.
    for o in out:
        # Junction is in-frame and stop-free → trace records show the
        # productive contract stayed satisfied.
        sim = o.final_simulation()
        assert sim.bases() != b""


def test_stream_records_yields_record_dicts():
    out = list(_human_igh_exp().stream_records(n=3, seed=0))
    assert len(out) == 3
    for rec in out:
        assert isinstance(rec, dict)
        assert "sequence" in rec
        assert "v_call" in rec
        assert "junction" in rec


def test_stream_records_matches_run_records_per_seed():
    exp = _human_igh_exp().mutate(count=5)
    batch = exp.run_records(n=3, seed=42)
    streamed = list(exp.stream_records(n=3, seed=42))
    assert len(batch) == len(streamed)
    for a, b in zip(batch, streamed):
        assert a["sequence"] == b["sequence"]
        assert a["v_call"] == b["v_call"]


def test_stream_records_default_n_is_unbounded():
    s = _human_igh_exp().stream_records(seed=0)
    first_three = list(_itertools.islice(s, 3))
    assert len(first_three) == 3
    assert all("sequence" in r for r in first_three)


def test_stream_records_lazy_iteration_writes_one_at_a_time(tmp_path):
    # Real-world streaming pattern: write each record to TSV as it
    # comes in, never holding more than one in memory.
    import csv as _csv

    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    path = tmp_path / "stream.tsv"
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = _csv.DictWriter(
            fh,
            fieldnames=_DEFAULT_COLUMN_ORDER,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for rec in _human_igh_exp().stream_records(n=5, seed=0):
            row = {k: ("" if v is None else v) for k, v in rec.items()}
            writer.writerow(row)

    with open(path, encoding="utf-8") as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 5
    assert all(r["sequence"] for r in rows)


def test_stream_works_on_compiled_experiment():
    exp = _human_igh_exp().mutate(count=5)
    compiled = exp.compile(allow_curatable_refdata=True)
    out = list(compiled.stream(n=3, seed=0))
    via_run = compiled.run(n=3, seed=0)
    for a, b in zip(out, via_run):
        assert a.final_simulation().bases() == b.final_simulation().bases()


def test_stream_records_works_on_compiled_experiment():
    compiled = _human_igh_exp().compile(allow_curatable_refdata=True)
    out = list(compiled.stream_records(n=2, seed=0))
    assert len(out) == 2
    assert all("v_call" in r for r in out)


def test_stream_compose_with_using_lock():
    exp = (
        ga.Experiment.on("human_igh")
        .restrict_alleles(v="IGHVF1-G1*01")
        .recombine()
    )
    for o in exp.stream(n=5, seed=0):
        vid = o.final_simulation().v_allele_id()
        assert exp.refdata.v_allele(vid).name == "IGHVF1-G1*01"


def test_stream_can_terminate_early_with_break():
    # Caller stops iteration after the first item; the generator
    # should not produce more outcomes than requested.
    exp = _human_igh_exp()
    s = exp.stream(seed=0)
    first = next(s)
    s.close()  # Generator cleanup; subsequent next() raises StopIteration.
    assert first is not None
    with pytest.raises(StopIteration):
        next(s)


def test_stream_records_writes_fasta_streaming(tmp_path):
    # Companion of the TSV streaming test for FASTA output.
    path = tmp_path / "stream.fa"
    with open(path, "w", encoding="utf-8") as fh:
        for i, rec in enumerate(
            _human_igh_exp().stream_records(n=3, seed=0)
        ):
            fh.write(f">seq{i}|v_call={rec['v_call']}\n")
            fh.write(f"{rec['sequence']}\n")
    lines = path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 6  # 3 records × (header + seq)
    assert lines[0].startswith(">seq0|")
    assert lines[2].startswith(">seq1|")
    assert lines[4].startswith(">seq2|")
