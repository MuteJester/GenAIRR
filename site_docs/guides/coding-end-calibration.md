# Retained coding end models

`CodingEndOutcomeModel` stores finite canonical `(five, three, weight)` outcomes
for V, J and optionally D. Positive outcomes delete reference bases; negative
outcomes retain original-tip P bases. D rows preserve correlation between ends.
These are retained outcomes of processing, rather than hairpin opening rates.

```python
from GenAIRR import CodingEndOutcomeModel, Experiment

model = CodingEndOutcomeModel(
    v=((0, -1, 1), (0, 0, 3)),
    d=((-1, 1, 2), (1, -1, 1)),
    j=((0, 0, 1),),
)
experiment = (
    Experiment.on("human_igh").recombine()
    .process_coding_ends(model=model, reference_ends_are_coding_ends=True)
)
```

The D distribution above contains two possible pairs. It never expands into
four combinations of independently sampled ends. The direct API also accepts
`d_joint=[(five, three, weight), ...]`, excluding separate `d_5`/`d_3` inputs.

The model is immutable and supports `to_dict()`/`from_dict()` for JSON storage.
Optional `reference_hashes` restrict compatibility to native reference content
identities. Compatibility is checked during configuration and again against
the final compiled reference after curation or genotype changes.

## Supervised estimation

`fit_generating_records()` requires completed `original_coding_end_v1` generating
provenance. It reads historical choices, so subsequent mutation and read loss
do not become recombination processing. It does not infer P from palindromes,
NP lengths or aligner annotations.

```python
fitted = CodingEndOutcomeModel.fit_generating_records(
    records,
    chain_type="vdj",
    reference_hash=compiled.refdata.content_hash(),
    source="independent native simulation founders",
    training_sampling_law="unconditioned_generating_choices",
    independent_generating_events=True,
)
report = fitted.provenance
held_out = fitted.score_generating_records(held_out_records)
```

The sampling-law and independent-event declarations are input contracts. Known
clone/lineage rows and collapsed duplicates are rejected; independently select
founders before fitting. Productive or otherwise selected samples need a
selection/ascertainment likelihood and are rejected by this supervised method.
Pooled frequencies also reflect allele-specific feasibility filtering and need
not recover the original pre-filter weights across heterogeneous alleles.

The report records counts, source, training assumptions and marginal 95% Wilson
intervals for independent observed events. These are not simultaneous intervals
or estimates of enzyme kinetics. Explicit prior support and nonnegative
`pseudocount` values are supported. Unobserved support is never invented.

## Compatible history likelihoods

`model.fit_compatible_histories()` provides an in-house finite mixture fitter.
Each observation supplies complete compatible histories with their forward
likelihoods. It sums explanations during expectation maximization rather than
assigning a single trim/P history to each observation. D outcomes remain joint.

```python
observations = [[
    {"outcomes": {"V": (0, 0), "D": (-1, 1), "J": (0, 0)},
     "likelihood": 0.8},
    {"outcomes": {"V": (0, -1), "D": (1, -1), "J": (0, 0)},
     "likelihood": 0.2},
]]
fitted = model.fit_compatible_histories(
    observations,
    source="specified synthetic forward observation model",
    training_sampling_law="unconditioned_generating_choices",
    independent_observations=True,
    complete_scenarios=True,
)
```

The caller's forward model must account for the remaining sequence-generating
processes, allele feasibility and data ascertainment. Generic AIRR fields do
not provide these likelihoods. The fitter reports convergence, likelihood
trajectory and residual history uncertainty. It does not establish parameter
identifiability or produce biochemical confidence intervals.

Validation includes supervised joint-outcome recovery, held-out scoring, an
identifiable ambiguous observation mixture and equivalent-history fixtures.
Measured biochemical rate calibration requires appropriate data and an
observation model.

The reproducible native recovery experiment uses 2,000 independent training
events and 500 held-out events, one feasible allele per segment, correlated D
outcomes and both D orientations:

```bash
python -m tools.metadata_biology_audit.coding_end_synthetic_recovery --output recovery.json
```

It reports every target and recovered probability, raw counts, marginal
intervals and held-out scores. Fixed tolerances are checked without choosing
seeds based on the fitted result. Empirical calibration and donor/clone
validation await measured data.
