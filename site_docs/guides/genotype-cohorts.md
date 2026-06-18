# Genotype cohorts

A single genotype models one subject. To simulate a **cohort** - many subjects,
each with their own genotype - use `run_cohort`. It runs the single-subject path
once per subject and collects the results:

```python
import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR.genotype import Genotype

cfg = gdata.HUMAN_IGH_OGRDB
# one sampled genotype per donor
donors = [Genotype.sample(cfg, seed=s, subject_id=f"donor_{s}") for s in range(5)]

cohort = ga.Experiment.on(cfg).recombine().run_cohort(
    donors, n_per_subject=200, seed=0, expose_provenance=True)

cohort.subject_ids            # ['donor_0', ..., 'donor_4']
len(cohort)                   # 1000 total records
cohort.result_for("donor_2")  # that donor's SimulationResult
cohort.to_csv("cohort.csv")   # combined, subject-tagged, unique sequence_id
```

`run_cohort` returns a `CohortResult`:

- `.subject_ids` / `.genotypes` / `.results` - per-subject, in input order.
- `.result_for(sid)` / `.refdata_for(sid)` - one subject's `SimulationResult` and
  the reference it ran against. Each subject can carry different novel/private
  alleles, so its *effective* germline (base catalogue + that subject's novels)
  differs; keeping it per subject is what lets `validate_records` and
  novel-allele truth calls stay correct across a heterogeneous cohort.
- `.records` - a fresh, subject-tagged, `sequence_id`-namespaced concatenation
  (each `sequence_id` is `"{subject_id}_{...}"`, so combined AIRR/FASTA export
  never collides).
- `.to_dataframe()` / `.to_csv()` / `.to_fasta()` - combined export over a stable
  union of columns.

**Record counts.** `n_per_subject` applies to all subjects; pass `counts` (a
parallel list, same length as the genotypes) to vary per-subject repertoire
sizes. A count of `0` is allowed - that subject appears in the cohort with zero
records.

```python
cohort = ga.Experiment.on(cfg).recombine().run_cohort(
    donors, counts=[500, 200, 0, 1000, 300], seed=0)
```

**Subject IDs.** Taken from each genotype's `subject_id`; if none are set they are
auto-assigned `subject_0..N-1`. Mixed (some set, some not) or duplicate IDs raise.

**Determinism.** Each subject gets an independent sub-seed derived from `seed`, so
a cohort is fully reproducible and subjects are independent.

`run_cohort` is mutually exclusive with `with_genotype`, `restrict_alleles`, and
`recombine(*_allele_weights=...)` (the genotype owns allele expression). It
**supports** `receptor_revision` (each subject's replacement V is restricted to
its own carried alleles); clonal forks are not combined with a cohort in this
release.

