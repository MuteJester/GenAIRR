# Duplex repair through clonal families

Slice D extends the finite authored repair model from checkpoint `a75db8c`
through `expand_clones`, `clonal_repertoire` and `clonal_lineage`. The
[design](duplex_family_design_2026_10_20.md) defines founder inheritance,
observation, abundance and replay. All simulation remains inside GenAIRR;
there is no external simulator or new dependency.

## Defaults without measured data

Omitting `DuplexHistory.weight` retains the existing value `1.0`. Selection is
uniform over positive eligible authored history IDs within the chosen source
context when their weights are equal. Local productivity conditioning and
physical support checks still apply. Different histories that produce the same
DNA retain separate probability mass. Giving one context more history rows
does not increase its upstream allele probability.

These are explicit unfitted controls over a supplied finite vocabulary. They
do not estimate biochemical rates or automatically enumerate all possible
hairpin-processing histories. Existing reference usage weights, uniform allele
fallbacks, mutation models and neutral lineage controls keep their established
meanings. Later fitted repair weights can use the same immutable model API.
The [defaults note](duplex_unfitted_defaults_2026_10_20.md) details this law;
the guide demonstrates an optional user metadata tag marking unfitted weights.

## Inheritance and actual branch history

Each clone samples and repairs its founder once. Duplex lineage survival
retries reuse that founder and retry growth, preserving the founding
rearrangement. Descendants share immutable repair evidence and inherit the
generating alleles, selected history, physical ledger and certified founder
placements. Each branch develops independent observation maps and edit
histories. Founder repair identities describe shared ancestry after cell
division.

Duplex lineage outcomes compose the founder and actual incoming-edge outcomes.
They retain cumulative division mutations, including reversions and repeated
changes at the same site. They are not reconstructed from a net difference
between final DNA and germline. Mechanism counters and mutation denominators
therefore remain tied to the actual counted events and their molecule.
`lineage_mutations_from_parent` reports the actual incoming-edge mutation count.
Downstream structural edits and read corruption retain their distinct counters.
Continuation composition validates its initial and final state boundaries.

Affinity scoring translates the contiguous phase-certified coding span,
including codons that cross junction boundaries. Neutral growth without a
target can proceed when coding phase is unresolved. Requested affinity scoring
in that situation raises an explicit error. Productive constraints use the
existing duplex observed-VDJ assessment at the founder, division and observed
read stages; no unknown productivity is silently promoted to productive.

## Records and abundance

Duplex lineage sampling collapses repeated draws of the same cell ID. Distinct
cells with identical DNA remain separate observations. Their `duplicate_count`
and lineage abundance count sample draws, conserving the requested sample
size. Duplex repertoire collapse compares complete generated record truth
apart from `sequence_id`, so equal DNA cannot hide distinct repair or mutation
histories. Fixed expansion continues to emit its descendants individually.
Legacy collapse and record shapes remain unchanged.

All three duplex workflows expose dense founder `parents` addressed by
`parent_id`. Lineage results additionally expose aligned `families` and
`lineage_trees`; accepted extinct trees retain their founder positions.
Nullable AIRR fields and canonical JSON ancestry cells survive ordinary
DataFrame and delimited exports; record dictionaries remain JSON serializable.
User metadata is applied after truth-based collapse through the ordinary
experiment API.

Family validation compares the complete immutable founder envelope, including
bindings, physical product, evidence, selected history and identity. Current
branch maps and edits can differ. Parent validation catches missing sidecars
even if all sibling dictionaries were edited. Duplex lineage validation also
checks actual record addressing, abundance, affinity and edge counts against
cached original native nodes. Changing both a clone and parent tag to another
equal-DNA tree does not evade this check. Exact value types are checked.

## Replay and engine changes

`result.families[i].replay_node(node_id, refdata)` reconstructs any retained
cell, including roots and unobserved nodes, using recorded division choices
without drawing new randomness. It checks reference identity, contracts,
mutator configuration, sequence, event counters and ancestry. The returned
state precedes downstream read observation. Incoming-edge outcomes are kept
for duplex nodes; full ancestor compositions are materialized for observed
cells or on demand. Immutable founder sharing is checked directly in native
tests.

This is an in-memory node replay interface. Generic durable plan traces still
reject merged family outcomes. Ordinary duplex source traces retain their
existing durable replay support; a persisted whole-family package is separate
work.

S5F lineage dispatch now returns structured replay errors instead of panicking.
Replay recognizes a missing site/base pair only when effective support is
actually empty. In duplex mode a selected positive site with an empty
destination row records a site-only choice; replay verifies that no-op and
rejects a truncated positive destination pair or out-of-range site. Fresh
sampling and legacy trace shapes remain unchanged. Finite duplex plan
identities now hash the exact numeric S5F kernel, preventing custom kernels
with the same display name from being substituted during source replay.

## Verification and practical scope

The [verification artifact](duplex_family_verification_2026_10_20.json) records
the final changed-source and extension hashes, complete suite results and
commands. Supporting controls use the final release extension:

- All 2,995 nonperformance Python tests pass, along with all eight existing
  performance checks under unchanged budgets. Full native debug and release
  suites pass 1,630 and 1,629 tests respectively, with the same three inherited
  ignored cases. Separately, 50 focused family/default integration checks and
  105 contract/family checks pass. The strict documentation build also passes.
- The independent [family sweep](duplex_family_sweep_2026_10_20.json) covers
  96 plans across all three workflows, both join orders, both D orientations,
  strict/permissive sampling, productivity constraints and read edits. Its
  192 founders produce 896 records; 928 node replays include unobserved cells
  and roots. Literal DNA, ancestry, cumulative edits, counters, mutation
  denominators and abundance checks report zero errors.
- The [projection regression](duplex_family_projection_regression_2026_10_20.json)
  covers 32,768 cases, including 15,360 intended incompatible projections;
  the [observation regression](duplex_family_observation_regression_2026_10_20.json)
  covers 2,048 founders and 18,432 continuations. Both report zero errors.
- All six actual Python guide examples execute successfully, including the
  unfitted provenance tag, family exports and exact node replay.

The first full Python run passed 2,994 tests and failed one obsolete
source-text architecture pin that forbade any parent-validator call in shared
clonal code. The replacement checks actual legacy calls and retains the
duplex-only stronger gate. After that test and its architecture note were
updated, the focused contracts and complete Python repeat passed. Production
sources and the release extension did not change between those runs.

The [isolated benchmark](duplex_family_benchmark_2026_10_20.json) uses six
processes, two clones per run and depths 2, 4 and 6. It retains 664 nodes and
exports 117 observed records, all replayed exactly. Median generation,
validation and metadata export takes approximately 77, 101 and 130 ms at the
three depths. Native node replay takes 9, 12 and 16 ms; replay plus validation,
export and equality checks takes 41, 55 and 72 ms. Compilation takes about
0.43 ms.

Maximum current RSS increases are 648 KiB during compilation, 10,380 KiB
during generation/export, and 172 KiB during replay. The largest increase from
before compilation through replay is 11,196 KiB. Peak high-water increments
are also recorded separately; startup high-water marks can mask later
allocations. These short 15-base controls and two repeats per depth are
descriptive, not full-length performance guarantees. Process RSS does not
establish individual allocation sharing, and this benchmark retains family
state rather than streaming it.

## Remaining work

Durable whole-family and external initial-snapshot formats, phased genotype
families, VJ repair, receptor revision, contaminants and paired-end adapters
remain explicit unsupported combinations. General flap/mismatch processing,
reopening completed joins and unrestricted overlapping layouts still need
their own physical and projection contracts. Automatic stochastic history
generation and measured biochemical calibration remain future work; equal
authored weights do not imply either is complete. Synthetic recovery is the
current calibration check while measured data are unavailable.
