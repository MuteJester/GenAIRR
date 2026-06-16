"""Tests for the non-tree clonal repertoire DSL (``clonal_repertoire``).

``clonal_repertoire`` generalizes the deprecated ``expand_clones``: per
clone a size is drawn from a heavy-tailed distribution (with an
unexpanded-singleton fraction), that many reads are emitted through the
post-fork library-prep / sequencing passes, and identical reads are
genotype-collapsed into AIRR records carrying ``duplicate_count``.
"""

from collections import Counter

import pytest

import GenAIRR as ga


def test_bcr_repertoire_with_corruption_has_clone_ids_and_duplicate_count():
    r = (ga.Experiment.on("human_igh").recombine()
         .clonal_repertoire(n_clones=20, max_size=50, unexpanded_fraction=0.3)
         .sequencing_errors(rate=0.005)
         .run_records(seed=0))
    assert len(r.records) > 0
    assert all("duplicate_count" in rec and rec["duplicate_count"] >= 1 for rec in r.records)
    assert all("clone_id" in rec for rec in r.records)
    per_clone = Counter(rec["clone_id"] for rec in r.records)
    assert min(per_clone.values()) >= 1


def test_pure_copy_no_postfork_one_record_per_clone():
    r = (ga.Experiment.on("human_igh").recombine()
         .clonal_repertoire(n_clones=10, max_size=100).run_records(seed=1))
    # no post-fork passes -> each clone collapses to ONE record with
    # duplicate_count = size
    per_clone = Counter(rec["clone_id"] for rec in r.records)
    assert all(c == 1 for c in per_clone.values())
    assert all(rec["duplicate_count"] >= 1 for rec in r.records)


def test_tcr_repertoire_runs_no_shm():
    r = (ga.Experiment.on("human_tcrb").allow_curatable_refdata().recombine()
         .clonal_repertoire(n_clones=10, max_size=50).run_records(seed=0))
    assert len(r.records) > 0
    assert all(rec.get("n_mutations", 0) == 0 for rec in r.records)
    assert all("duplicate_count" in rec for rec in r.records)


def test_tcr_repertoire_rejects_mutate():
    with pytest.raises(Exception):
        (ga.Experiment.on("human_tcrb").allow_curatable_refdata().recombine()
         .clonal_repertoire(n_clones=5, max_size=10).mutate(rate=0.05).compile())


def test_determinism():
    mk = lambda: (ga.Experiment.on("human_igh").recombine()
                  .clonal_repertoire(n_clones=8, max_size=40).run_records(seed=3))
    a, b = mk(), mk()
    assert [r["sequence"] for r in a.records] == [r["sequence"] for r in b.records]
    assert [r["duplicate_count"] for r in a.records] == [r["duplicate_count"] for r in b.records]


def test_validation_works():
    r = (ga.Experiment.on("human_igh").recombine()
         .clonal_repertoire(n_clones=5, max_size=30).sequencing_errors(rate=0.003)
         .run_records(seed=0, validate_records=True))
    assert len(r.records) > 0


def test_clonal_repertoire_rejects_premature_mutate():
    # .mutate() before clonal_repertoire is a descendant-phase step and
    # must be rejected (SHM is descendant-specific).
    with pytest.raises(ValueError):
        (ga.Experiment.on("human_igh").recombine()
         .mutate(rate=0.05).clonal_repertoire(n_clones=5, max_size=10))


def test_clonal_repertoire_rejects_double_fork():
    with pytest.raises(ValueError):
        (ga.Experiment.on("human_igh").recombine()
         .clonal_repertoire(n_clones=5, max_size=10)
         .clonal_repertoire(n_clones=3, max_size=5))


@pytest.mark.parametrize("method", ["invert_d", "receptor_revision"])
def test_recombination_time_edits_reject_after_clonal_repertoire(method):
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .clonal_repertoire(n_clones=5, max_size=10)
    )
    with pytest.raises(ValueError, match="before the clonal fork"):
        getattr(exp, method)(prob=0.05)


@pytest.mark.parametrize("method", ["invert_d", "receptor_revision"])
def test_recombination_time_edits_reject_after_clonal_lineage(method):
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .clonal_lineage(n_clones=1, max_generations=2, n_sample=2)
    )
    with pytest.raises(ValueError, match="before the clonal fork"):
        getattr(exp, method)(prob=0.05)


def test_clonal_repertoire_validates_args():
    base = ga.Experiment.on("human_igh").recombine()
    with pytest.raises(ValueError):
        base.clonal_repertoire(n_clones=0)
    with pytest.raises(ValueError):
        ga.Experiment.on("human_igh").recombine().clonal_repertoire(
            n_clones=5, size_distribution="bogus")
    with pytest.raises(ValueError):
        ga.Experiment.on("human_igh").recombine().clonal_repertoire(
            n_clones=5, exponent=0.0)
    with pytest.raises(ValueError):
        ga.Experiment.on("human_igh").recombine().clonal_repertoire(
            n_clones=5, sigma=-1.0)
    with pytest.raises(ValueError):
        ga.Experiment.on("human_igh").recombine().clonal_repertoire(
            n_clones=5, max_size=0)
    with pytest.raises(ValueError):
        ga.Experiment.on("human_igh").recombine().clonal_repertoire(
            n_clones=5, unexpanded_fraction=1.5)
