# Finite authored duplex repair in Experiment

Slice C follows the native seed checkpoint `6413a6b` and the
[approved integration design](duplex_airr_projection_and_integration_design_2026_10_20.md).
The opt-in entry point is `Experiment.recombine(duplex_repair=model)`, with
`DuplexHistory` and `DuplexRepairModel` in `GenAIRR.experimental`.

## Configuration and lowering

An immutable model pins exact reference content and allele evidence, authored
row order, stable unique history IDs, finite nonnegative weights, canonical
source sequences and generating names, D orientation, ordered operations and
work budgets. `bind` captures the reference identity; compilation validates
actual bindings and histories against the final reference catalogue, including
an explicitly authored effective genotype catalogue. Serialization contains
plain data and no Python execution callback.

One atomic native source pass follows existing allele/genotype and optional
D-orientation sampling. It owns the completed V-D and D-J repair product.
Automatic cartridge trim/P/NP assembly defaults are bypassed. Explicit legacy
NP settings, even an explicit zero distribution, and explicit trim calls, even
`trim(enabled=False)`, are rejected. Original-end processing is a distinct
model and cannot be combined with this mode.

All naturally reachable positive contexts must have positive authored support.
Explicit allele restrictions narrow that domain. Phased contexts are derived
from expressible complete chromosomes rather than the Cartesian product of
marginal allele sets. Zero chromosome, gene-usage and copy-expression weights
cannot establish natural reachability. Every active positive row is checked
through the shared physical kernel, exact bindings and seed capability before
productivity filtering. Valid overlapping/nonrepresentable products remain
available through standalone projection; they reject the requested ordinary
experiment with row-specific diagnostics. No failed row is dropped or retried.

The complete plan is checked before sampling. The supported suffix is uniform
or S5F mutation, PCR/quality substitution, N masking, structural indels, terminal
loss and read orientation. VJ, family simulation, receptor revision,
contamination, paired-end windows, arbitrary custom contracts and the legacy
preservation-only `AnchorPreserved` require further adapters and fail explicitly.

## Probability and bounded work

The declared law is `finite_authored_history_local_conditioning_v1`.
Within a selected allele/orientation context, unconstrained history probability
is its authored positive weight divided by that context's weight sum. Under
`productive_only`, only rows assessed as definitely productive are eligible,
and their original weights are normalized locally. False and unresolved
assessments are ineligible.

Earlier allele, gene, haplotype and orientation choices retain natural relative
weights among choices with a feasible completion. They receive an existential
completion mask, never a multiplier derived from the number or success mass of
repair histories. This is sequential local support conditioning, not globally
conditioning the joint generation law. Different rows with identical products
remain separate choices with their own mass and IDs. Zero-weight rows are
never sampled and cannot satisfy missing-context coverage.

For a phased genotype, an eligible gene retains its natural mass: gene usage
times the total positive carried-copy dosage on the selected chromosome.
Only the subsequent within-gene copy draw is normalized over eligible copies.
For example, a gene with copy masses 1 and 9 retains gene mass 10 when only
its mass-1 copy has a productive completion. Against a second eligible gene
with natural mass 1, the first gene retains 10:1 odds. Completion filtering
cannot turn those odds into 1:1.

`max_preflight_operations` charges reachable-context checks, authored
operations and source weight entries inspected for natural and distinct
productivity-filtered flat supports. Genotype checks charge chromosome, explicit
usage, copy and derived active gene-mass entries; their realized filtered CDFs
are also checked before drawing or replaying. `max_preflight_bases` reserves each active
script's declared base budget before preparing it. These limits control
compilation work and allocation reservations; they carry no biochemical rates
or wall-clock guarantee. Exhaustion is a resource error, not proof of impossible
productive support. Max-weight scaling
prevents total-weight overflow; numerical loss of positive support during
scaling or cumulative addition is rejected explicitly. Flat source masks are
checked during compilation; genotype's realized chromosome/gene/copy CDFs
also fail explicitly before the affected draw or replay. The numerical domain
is binary64 authored weights and the native RNG, not arbitrary real numbers.

## Truth, output and replay

The persistent seed and observation model from slice B preserves the founder
ledger, physical origins, template ancestry, certified placements, structural
versus observation loss, native event counters and current coding assessment.
The selected history ID/index, natural weight, model digest and conditioning
law travel with founder provenance. Current allele re-calls remain separate
from the standard generating calls and alignment reference.

Normal `run_records`, compiled execution and streaming use nullable duplex
record dispatch. Each record carries self-contained ancestry and current
alignment paths. Strict export converts only AIRR coordinates and uses the
oriented D offset/zero-width query boundary for CIGARs. Nested provenance stays
in its declared native coordinate systems and serializes as canonical JSON in
delimited exports. Nullable DataFrames preserve missing integers and booleans.
Reserved generated fields cannot be replaced by custom metadata.

Validation checks the actual mutable record, including nested JSON value types,
against validated native repair/observation evidence and independent productivity
checks. Legacy scalar estimators continue to reject duplex provenance rather
than treating unknown P values as zero or NP geometry as pure N addition.

Source-backed trace files start from empty IR and bind the full model and
upstream sampling law, reference identity, contracts and duplex policy. Replay
consumes the recorded positive row in its exact context and re-executes the
shared kernel without a replacement draw. Wrong rows, unsupported or zero-mass
choices, edited bindings, model operations/weights and policy downgrades fail.
A private cached execution identity prevents a trace writer from relabeling an
outcome generated under different source weights as its own. External-seed
`run_from` traces still need a durable initial-snapshot format and remain gated.

Old-mode seeded behavior and legacy trace signatures remain unchanged. The
stronger complete sampling-law signature is selected for finite duplex source
plans. Family integration, general flap/mismatch repair, reopening completed
joins and fitting measured biochemical rates remain separate milestones.

## Verification

Verification of the final source and release extension:

- 2,953 Python tests, plus all 8 existing performance checks with unchanged
  budgets; 68 focused integration checks also passed separately.
- Full native debug and release suites: 1,618 and 1,617 passes respectively,
  with the same three inherited ignored cases.
- 32,768 projection controls, including 15,360 intended incompatibility
  rejections; 2,048 founders and 18,432 observation continuations; zero errors.
- 6,144 exact finite-law/metadata controls across 12 plans, plus 512 exact
  genotype gene-choice controls under strict and permissive policies.
- Strict documentation build and all five actual Python guide examples,
  including nullable TSV export and durable replay.

An initial W6 performance run exceeded its one-second budget at 1.187 seconds.
The isolated profiled case and the complete eight-test repeat passed without
code or budget changes; the initial slowdown's cause was not established.

The compilation benchmark used 18 isolated processes and streamed 18,000
records. Median compilation took approximately 1–14 ms for 8–128 authored
histories across 1 or 8 contexts. Streaming 1,000 records took 22–26 ms.
Maximum observed RSS increases were 4,544 KiB during compilation and 1,024 KiB
during streaming. These are short 15-base synthetic controls; RSS does not prove
individual immutable allocations are shared, and each stream measured only
one full metadata export.

The [verification artifact](duplex_experiment_verification_2026_10_20.json)
pins the final binary and changed-source hashes, commands, results and
supporting artifacts. Existing A/B verification artifacts retain their
historical results. The [guide](../site_docs/guides/experimental-duplex-repair.md)
contains the executable public API examples.
