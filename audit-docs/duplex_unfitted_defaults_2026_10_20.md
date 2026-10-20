# Unfitted duplex defaults

The finite duplex repair model can run without measured biochemical data.
Every `DuplexHistory` has a default `weight=1.0`. Leaving these weights out
therefore defines an equal-weight control over the explicitly supplied history
IDs. It does not require an external simulator, fitting service or additional
dependency.

## What is uniform

After choosing the generating V, D and J alleles and D orientation, GenAIRR
selects among the positive, eligible authored rows for that context. If the
context has `k` eligible rows and all weights are omitted, each row has
probability `1/k`. With `productive_only()`, eligibility additionally requires
a positive productivity assessment; the surviving rows retain equal weights.
Physically invalid or unsupported positive histories fail compilation before
this productivity filter. No invalid history becomes valid through a default
weight, and failed histories are not hidden by retries.

Uniformity applies to **history IDs**, including histories that produce
identical DNA. A context containing three equal-weight rows that make the same
sequence and one row that makes another sequence assigns those sequences
probabilities `3/4` and `1/4`. Duplicating a row under another ID intentionally
adds another history atom and changes that context's probability distribution.
There is no automatic deduplication by final sequence.

Upstream choices retain their own laws. Existing reference usage weights are
used when supplied; otherwise ordinary allele sampling is uniform over its
selectable support. Explicit allele, genotype, chromosome and orientation
controls remain effective. Giving one allele context more history rows does
not increase that allele's upstream probability. The current law conditions
locally within a selected context and must not be described as globally
conditioning the full joint distribution on productivity.

## What these defaults establish

Equal weights are transparent, unfitted simulation controls. They do not
estimate the frequency of particular hairpin openings, strand resection,
polymerase actions or joining orders. The current interface requires complete
authored scripts, with bounded physical operations and reference bindings.
It does not generate arbitrary stochastic chemistry histories automatically.
Uniform defaults over a declared finite support make that support usable;
they cannot establish that the support is biologically exhaustive.

Family integration inherits the selected repair once from its founder. Each
descendant carries that founder history while acquiring its own observed
changes and coverage. Redrawing a fresh repair history for a descendant would
change the founding rearrangement and is outside this inheritance contract.
The existing lineage controls, including neutral `selection_strength=0.0`,
remain simulation settings. Their existing defaults are not newly inferred
biochemical repair parameters. Existing bundled mutation kernels keep their
declared meanings; equal repair weights do not imply that every other model
component is uniform.

## Provenance and later data

The serialized model stores the ordered row IDs, scripts, weights and exact
reference identity. Its digest and sampling law are bound into durable source
traces. Each record retains the selected history ID/index, natural weight,
model digest and local conditioning law. Omitting a weight and explicitly
supplying `1.0` produce the same serialized model and reproducible run.

Measured or fitted nonnegative finite weights can later replace the row
weights through the existing `DuplexHistory(..., weight=...)` API. The finite
support, reference context, ascertainment and sampling law must match what is
being estimated. Generic final AIRR sequences alone do not identify unique
strand-by-strand histories. The older `CodingEndOutcomeModel` deliberately
rejects duplex histories; a future estimator must consume their actual
history-aware likelihood or generating provenance. Until such data and a
suitable estimator exist, the default rows remain authored, unfitted controls.

## Regression evidence

`tests/test_duplex_uniform_defaults.py` checks omitted and explicit equal
weights through canonical serialization, complete source choices, records,
durable choice replay and seed reruns. It also checks equal mass for distinct
history IDs with identical final DNA, and confirms that changing the number
of histories in a context preserves earlier allele choices under both
uniform and explicitly weighted allele sampling.

See the [finite experiment implementation](duplex_experiment_implementation_2026_10_20.md)
for the physical, productivity, resource and replay contracts shared by these
defaults.
