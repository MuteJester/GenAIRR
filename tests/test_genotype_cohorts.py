"""Cohorts: Experiment.run_cohort + CohortResult."""
import pytest

import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR.genotype import Genotype
from GenAIRR.cohort import CohortResult, CohortSubjectResult


def _cfg():
    return gdata.HUMAN_IGH_OGRDB


def test_cohort_result_core_accessors():
    cfg = _cfg()
    g0 = Genotype.sample(cfg, seed=0, subject_id="A")
    g1 = Genotype.sample(cfg, seed=1, subject_id="B")
    r0 = ga.Experiment.on(cfg).with_genotype(g0).recombine().run_records(n=3, seed=1)
    r1 = ga.Experiment.on(cfg).with_genotype(g1).recombine().run_records(n=2, seed=2)
    subjects = [
        CohortSubjectResult(subject_id="A", genotype=g0, result=r0, refdata=object(), seed=11, count=3),
        CohortSubjectResult(subject_id="B", genotype=g1, result=r1, refdata=object(), seed=22, count=2),
    ]
    c = CohortResult(subjects)
    assert c.subject_ids == ["A", "B"]
    assert c.genotypes == [g0, g1]
    assert c.results == [r0, r1]
    assert c.result_for("B") is r1
    assert c.refdata_for("A") is subjects[0].refdata
    assert len(c) == 5                      # total records
    assert len(c.records) == 5              # fresh concatenated list
    c.records.append({"x": 1})              # mutating the returned list...
    assert len(c.records) == 5              # ...does not affect the cohort
    with pytest.raises(KeyError):
        c.result_for("NOPE")
