# GenAIRR productivity and indel correctness plan

Proposed on 2026-10-20 for `fix/biological-metadata-correctness`.
This is the approved implementation plan. Current results and remaining model
boundaries are in the [implementation report](productivity_indel_repair_implementation_2026_10_20.md).

The next release should first make productivity decisions, constrained indel
sampling and replay trustworthy. Reference annotations need a reproducible repair
process. Detailed duplex repair and empirical calibration can then build on those
foundations as separately validated models.

All implementation, simulation, annotation resolution and estimation will remain
inside GenAIRR. Reference releases, scientific papers and experimental data are
inputs, not external software dependencies. This plan adds no dependency on IGoR,
OLGA, external annotation programs or hosted inference services.

## Evidence and priorities

The [original limitation report](original_end_existing_limits_2026_10_20.md)
contains the mouse reference failures and the fixed-precursor demonstration that
strict indel failures occur despite nonempty support. Additional read-only probes
and source inspection during this planning review found related problems:

| Finding | Evidence | Consequence |
| --- | --- | --- |
| Invalid annotated anchors can pass as productive | Legacy productive generation accepts synthetic V `AAACCCGGG`, anchor 6, and J `TTTGGG`, anchor 0; its junction is `GF`. Record validation accepts the result. | Reference admission, anchor preservation and biological acceptance are currently conflated. |
| Ambiguous coding sequence can be declared stop-free and productive | V `AAANNNTGT`, anchor 6, with J `TTTGGG`, anchor 0, yields `KXCFG`, `productive=True`, `stop_codon=False`. | Translation to `X` is being treated as evidence of no stop. |
| Reference rules and productive solvers disagree | A custom J rule allowing Y accepts J `TATGGG`, anchor 0. Legacy productive generation succeeds; original-end productive compilation rejects it. | A shared, explicit productivity policy is needed. |
| Indel reachability can exclude valid histories | Starting from V `TGT`, J `TGG`, no N, insert A at position 3 and delete that inserted base at position 3. Replay accepts this productive pair, but the fresh sampler's stationary continuation calculation excludes that first insertion. | The existing calculation is not a sound proof of empty support. A new event can create a later legal site. |
| Replay permits impossible or malformed events | An insertion coordinate `2**32 + 6` replays at 6; an insertion replays with insertion probability zero; a fixed-count-one plan accepts two events; mixed no-op and real events can replay on a nonempty sequence with insertion probability one. | Replay must validate typed coordinates, positive proposal support and explicit no-op provenance. |
| Reference construction loses evidence | The Python bridge forwards the anchor but discards `anchor_meta`; native coding phase is derived from the anchor. A small FASTA build passes checksum integrity with missing anchors. | Checksums do not certify biological annotations. |
| Reference authoring has an incomplete native migration | The FASTA/allele path refers to the absent `GenAIRR._native._anchor` module and can fall back to missing anchors. The inspected fixture also loses header functionality, and omitted identity inference defaults to Human with a warning. | Repair the authoring path before rebuilding affected cartridges. |

These are separate from the implemented original-tip P materializer. Its
ancestry checks remain useful, but cannot certify the biological interpretation
of invalid reference landmarks.

## 1. Define productivity once and use it throughout the engine

Introduce a versioned `ProductivityPolicy` and a pure `ProductivityAssessment`.
The assessment should contain individual predicates, their evidence status,
reason codes and the evaluated sequence scope. Production consumers should use
this result rather than reimplementing the final Boolean.

Keep three concepts distinct:

- **Reference admission:** whether a catalogue may be loaded under the chosen
  validation mode.
- **Reference readiness:** whether the particular allele has enough trustworthy
  landmark, phase and coverage information for the requested calculation.
- **Sequence acceptance:** whether the realized sequence satisfies the declared
  productivity policy.

`allow_curatable_refdata()` changes admission only. It cannot certify a bad
anchor, invent a coding phase or relax the productive policy.

The default policy should be explicitly named as a GenAIRR sequence-level V(D)J
policy: usable V/J landmarks, V Cys and J Phe/Trp, compatible coding frames,
intact required anchor provenance, and no coding stop in the evaluated span.
This does not establish receptor expression, correct splicing, regulatory
function or antigen binding. AIRR permits nullable productivity and does not
mandate one strict productivity definition. [AIRR Rearrangement specification](https://docs.airr-community.org/en/stable/datarep/rearrangements.html)

The canonical anchor convention follows the conserved V-domain positions
identified by IMGT, rather than an arbitrary terminal cysteine or phenylalanine.
[IMGT V-domain numbering](https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html)

Catalogue validation rules can retain locus-specific expectations. A synthetic
or nonstandard productivity policy must be explicit and included in the model
identity; changing catalogue validation rules alone must not silently change
what `productive=True` means. Preserve `AnchorPreserved` as a standalone
optional preservation contract. The default productive bundle should instead
use reference-readiness, anchor-provenance integrity and live residue-class
acceptance. In particular, a J W-to-F or F-to-W change passes the canonical
residue criterion when all other predicates pass. An explicitly requested
`AnchorPreserved` contract additionally requires the original amino acid.
This removes a remaining ambiguity between canonical sequence acceptance and
exact reference-amino-acid conservation; it is part of the versioned behavior
change.

### Truth values and coverage

Use three outcomes. A known defect takes precedence over uncertainty only when
the evidence supporting that defect is itself reliable.

| Available evidence | Standard productivity result |
| --- | --- |
| All required predicates pass with sufficient observed coverage | `True` |
| A trusted frame, stop or required-anchor defect proves policy failure | `False` |
| Required reference landmark, coding phase, base identity or coverage is unresolved | `None` |

An invalid annotation is not evidence that the biological molecule is
nonproductive. Suppress standard junction coordinates/content when the reference
landmark cannot define a defensible junction. Conversely, retain a junction
defined by trustworthy original landmarks when mutation changes its amino acid;
report its actual translation and assess it under the selected policy.
If source evidence establishes a genuine noncanonical residue at the correctly
identified landmark, its failure of the canonical policy is determinate.
Uncertainty about where the landmark is located remains a different condition.

For ambiguous codons, evaluate possible translations under the supported
alphabet. `NNN` permits both stop and sense codons and cannot establish
`stop_codon=False`. An ambiguous codon whose possible resolutions are all sense
codons can still establish no stop at that codon. Unknown reading phase prevents
a stop conclusion unless independent evidence establishes it.

Standard flags should describe the final observed evidence. Terminal loss of
required coding evidence produces uncertainty, not an automatic biological
failure. `stop_codon=False` may describe the observed complete codons without
certifying absent portions. Loss outside the evaluated coding span need not
change the result. If molecular productivity before observation is exposed, it
must be a separate, stage-labeled assessment taken before read corruption.
Terminal clipping that retains the required landmarks and sufficient phase
evidence can retain a determinate observed-span assessment; this does not
certify the absent portion of the source molecule.

This intentionally changes some current `True` values to `None`. Preserve the
existing public promise that observed-scope `productive_only()` emits records
with `productive is True`: unknown-producing observations must be excluded or
rejected. A future explicit molecule-scope option can allow productive source
molecules to yield partial or damaged reads without mislabeling those reads.
Do not silently change the scope of the existing method.

### Integration

Wire the assessment and policy through:

- `engine_rs/src/coding.rs`, a small new productivity module, and reference
  readiness helpers.
- Both legacy and P-aware feasibility paths, allele/genotype selection,
  revision and the productive contract bundle.
- Mutation support, complete structural-event checks, native final-output
  checks and Python repertoire/clonal execution.
- `airr_record/builder.rs`, `airr_record/junction.rs` and validation.
- Nullable Python/JSON/DataFrame/TSV conversion and reserved diagnostic fields.

An intermediate assembly may be incomplete without constituting a contract
violation. Distinguish “deferred until assembly completes” from “passed.”
At the final emission boundary, the productive bundle must verify that the
assessment is definitively accepted, even if a permissive pass was skipped.
A failed final assertion is an engine error, never a reason to relabel or
silently discard the record.

Generation and projection should share policy semantics. The record validator
should reconstruct evidence from sequence/provenance independently of the
builder, including checking required absences when a junction is undefined.
Annotation limitations and record inconsistencies need separate diagnostics:
an honestly unknown record can still be internally valid. Tests should also
independently reconstruct landmarks, frames and translations from literal
fixtures. Comparing two callers of the same flawed helper is insufficient.

## 2. Repair and version the reference authoring pipeline

First restore a single GenAIRR-native anchor resolver used by both allele
construction and the cartridge builder. Missing resolver functionality must be
reported as an authoring failure, rather than silently producing an apparently
finished functional cartridge. Diagnostic partial catalogues remain possible
with explicit annotation-incomplete status.

Add an evidence record per allele, reusing the existing manifest/report system:

- Source accession/release and sequence checksum.
- Species, locus, functionality, aliases and source orientation.
- Anchor coordinate, resolution method, evidence and confidence.
- Independently supported coding phase and coding interval, where available.
- Completeness of each junction-facing coding boundary.
- Any annotation correction: previous value, replacement and supporting source.

Carry this evidence across the Python/Rust boundary. The existing Python
`JAllele.frame` is derived from the anchor and is not independent phase evidence.
Validated legacy anchors may provide a clearly labeled compatibility phase;
questionable anchors must not become a phase source. A known phase does not
establish that a junction landmark is correct.

Rebuild affected cartridges from pinned source snapshots and an annotation
patch ledger. Automatic motif searches may suggest candidates, but cannot
authorize a correction. Keep unresolved annotations explicit and exclude them
from certified productive support. Preserve real nonfunctional alleles for
unconstrained simulations.

Audit all bundled cartridges, starting with mouse TCRA/TCRG and the separate
mouse IGK/IGL duplicate-name/alphabet failures. Resolve duplicates using source
identity and allele semantics; do not silently drop rows or substitute ambiguous
bases to make validation pass. Preserve functional-status headers and require
explicit or supported species identity; unknown species must remain unknown.

Separate structural integrity, annotation completeness and suitability for a
named simulation capability in build reports. Keep empirical parameter overlays
distinct from structural reference releases and declare their compatibility.
A rebuilt reference receives a new content identity. Extend versioned hashing
and loading for new fields while retaining historical hash verification for old
trace formats.

Release gate: every changed annotation is traceable; unchanged allele sequences
remain byte-identical; any removals or renames are listed; unresolved data remain
visible; independent checks validate the certified subset.

## 3. Replace constrained indel sampling with a defined probability model

The present sampler has two problems: its continuation calculation assumes a
future site distribution based on the current state, and its bounded final
rejection loop equates failure to find a tuple with proof of impossibility.
Increasing the retry budget fixes neither definition.

Adopt a versioned natural sequential law for the indel pass. For a fixed input
state and an already drawn event count `k`, each event has its original
probability:

```text
insertion at site i of base b: p_insert / (current_length + 1) * p_base(b)
deletion at site i:          (1 - p_insert) / current_length
```

Preserve the explicit empty-pool deletion semantics. Let `A` denote all declared
path requirements and the final tuple contract. The target distribution is:

```text
P(events | input, k, A)
  proportional to
  product_i q(event_i | actual state after preceding events) * indicator(A)
```

This conditions this indel batch only. It does not change recombination's
existing local support conditioning, condition the entire experiment globally,
or reweight the user's count distribution. An impossible already-selected count
must not be silently redrawn.

Before implementation, pin which constraints are path requirements and which
are final-state requirements. Under the recommended final-state productive
policy, temporary frame shifts and temporary stops are allowed inside an
uncommitted tuple if the final state satisfies the policy. Inserting between
anchor bases and subsequently removing that insertion can also restore the
anchor. Therefore, the current per-event `Forbidden` classification is not
automatically a sound pruning rule. Deleting an original anchor base is
irreversible under the existing provenance-preservation contract and may be
pruned when that contract is active.

### Exact support for short tuples

Implement an exhaustive reference sampler first for small fixtures, then a
production exact solver for counts 0, 1 and 2. Compute continuation mass from
each candidate's actual successor state:

```text
Z(state, 0) = 1 if the complete tuple is accepted, otherwise 0
Z(state, r) = sum_event q(event | state) * Z(apply(state, event), r - 1)
sample weight(event) = q(event | state) * Z(successor, r - 1)
```

The initial state, policy and any relevant history are fixed inputs to this
recurrence. They must participate in memoization identity when they affect
acceptance. Equal sequence strings alone do not imply equal provenance or
contract state. Multiple histories reaching the same sequence retain their
combined probability mass; histories are not deduplicated as sequences.

Separate Boolean reachability from floating-point mass. Numerical underflow
must not produce an empty-support verdict. Use stable scaling/log sums and
diagnose invalid weights separately.

For the practical count-two path, use a temporary edit view over the original
sequence with coordinate/provenance mapping. Precompute codon summaries for all
three phases and compose unchanged intervals plus edited boundaries. Check
anchors and the J tail against their actual transformed positions. This avoids
cloning and scanning a full `Simulation` for every candidate pair. Validate this
optimized evaluator against both the exhaustive oracle and fully materialized
native contract checks. Commit only the chosen tuple through
`MutationTransaction`, preserving the existing observer/provenance machinery.

The count-two release gate includes the documented 370-nt precursor and the
minimal insertion-then-deletion counterexample. Any claimed guaranteed operating
range must be backed by complete enumeration and measured cost; do not advertise
an arbitrary “up to six indels” guarantee.

### Higher counts and work limits

Use the same natural law for larger tuples. A practical first implementation
can use independent complete-tuple rejection sampling followed by exact
continuation search within an explicit work budget. A successful rejection
draw has the intended conditional distribution. A completed fallback must
sample by the same law; “return the first valid tuple” is biased.

Budget complete independent attempts. Interrupting expensive sampled tuples
partway through validation can preferentially retain cheap histories. Exact
search may have a deterministic node budget, but an unfinished search must
discard partial masses and report exhaustion, rather than sampling its explored
subset.

Do not combine the historical approximate frame-weighted proposal with a
natural-law fallback. Those are different distributions. Rejecting an early
prefix is permitted only if it provably cannot complete, or if the path policy
itself rejects that prefix.

Expose distinct outcomes:

| Outcome | Required evidence or behavior |
| --- | --- |
| `EmptyAdmissibleSupport` | Exhaustive enumeration or a sound impossibility proof |
| `SamplingBudgetExceeded` | Search ended without a sample or proof; support remains unresolved |
| `SupportUnavailable` | Required proposal support cannot be enumerated |
| Invalid support/numerical error | Malformed probabilities or a computation failure |
| Replay violation | The recorded history is invalid for the declared model |

Extend the native fallible execution path and Python exception mapping so budget
exhaustion never causes the permissive `.expect(...)` path to panic. Recommended
default: report resource exhaustion in either execution mode. If a user
explicitly chooses to skip such a batch, expose its reason and requested/applied
counts; it is not a successful realization of the requested indel count.

Record search strategy, version, work performed and completion status in
diagnostics. Use deterministic work limits rather than wall-clock termination
to preserve reproducibility. Batch APIs must not hide state-dependent failures
by discarding failed records, which would change the emitted repertoire.

General exact acceleration beyond two events is a benchmarked follow-up using
memoized edit states and sound summaries. Arbitrarily large counts with custom
contracts cannot honestly be promised cheap, guaranteed exact sampling.

## 4. Make replay enforce the same model

Replay should consume recorded choices without new randomness, validate them
against the declared proposal law and contract semantics, then commit atomically.

Use checked integer conversions before narrowing sites or counts. Validate the
recorded count against the configured count source, kind probability, positive
finite insertion-base support, evolving bounds and final tuple predicates. Use
the new model's path rules rather than retaining old prefix restrictions that
exclude recoverable intermediate states.

Distinguish an empty-pool deletion from an explicitly skipped batch. Preserve
older permissive trace sentinels through a documented legacy decoding path,
but do not allow a mixed fabricated no-op/event history to masquerade as a
natural new-model draw. Record skip reasons outside biological event provenance.

Include sampler law, productivity policy, scope and behavior-affecting settings
in signatures. Existing refdata hashes already include rules, functionality and
subregions; extend that system rather than building a competing identity scheme.
Version any changed hash encoding and trace compatibility checks.

New constrained runs will intentionally change some seeded results. Valid old
event traces may still reproduce their sequence through explicit compatibility
handling, but current metadata must not revive the old false productivity
claims. Do not promise bit-for-bit old annotations after correcting their
meaning.

Correct stale documentation at the same time, including the public
`productive_only()` docstring's claim that TCR compilation is unsupported.

## 5. Acceptance tests and release evidence

The correctness release must include the cases below.

1. **Reference and productivity:** literal valid/invalid/unknown anchors;
   custom admission rules versus explicit policies; wrong-anchor phase;
   ambiguous codons; trustworthy mutated anchors; J-tail shifts; internal and
   terminal coverage loss; reverse observations; full read loss.
2. **Cross-path agreement:** legacy, P-aware and original-end recombination;
   genotype restrictions; V replacement; mutation; native and Python execution;
   clonal branches; fresh and replayed output. Every observed-scope productive
   emission is definitively productive under its recorded policy.
3. **Exact indel oracle:** enumerate small sequences, counts and skewed weights;
   compare legal histories and normalized masses, not just successful outputs.
   Include insertion cancellation, temporary anchor interruption, same-site
   edits, terminal edits, noncoding frame-neutral edits and empty-pool events.
4. **Probability validation:** compare optimized and exhaustive continuation
   masses numerically, then use seeded distribution tests with predetermined
   tolerances. Frequency checks supplement deterministic checks.
5. **Failure taxonomy:** force tiny budgets with known support; genuine empty
   support; extreme weights; unavailable distributions; zero-event cases.
   Verify error class, no partial commit, and explicit skip accounting.
6. **Replay rejection:** oversized and negative sites, impossible event kinds,
   unsupported bases, mixed sentinels, missing choices, changed policies and
   changed reference evidence. Verify accepted replay preserves P ancestry and
   all observed coordinate/counter metadata.
   Separately tamper with record `True`/`False`/`None` values, reuse an invalid
   reference phase, and populate junction fields whose boundaries are undefined.
   Validation must catch these cases even when normal projection would return
   early.
7. **Cross-locus sweep:** restore strict productive indels to the earlier stress
   matrix; retain the original seeds and fixed precursors; add broad seed
   coverage. Report invalid references and budget limits separately from
   biological outcomes. An uncorrected cartridge is not a successful test.
8. **Compatibility and performance:** run the full native debug/release and
   Python suites against the rebuilt extension. Benchmark no-contract baseline,
   counts 1/2/4/6, realistic read lengths and skewed probabilities; record
   throughput, tail latency, memory, fallback frequency and exhausted budgets.
   Publish the exact build/reference identities with results.

Promote temporary reproducers into permanent, minimal regression fixtures when
implementation starts. The planning probes are not a replacement for these
release tests, and the previous green suites predate the proposed fixes.

## 6. Extend biochemical repair as a separate model

Keep `original_coding_end_v1` available with its stated retained-outcome
semantics. Develop an opt-in duplex model after the reference and correctness
work above.

Primary biochemical studies support sequence-dependent hairpin opening and
joint effects of end processing, polymerases and microhomology. They do not
supply a universal repertoire-wide frequency table. [Lu, Schwarz and Lieber, 2007](https://pubmed.ncbi.nlm.nih.gov/17932067/),
[Gu and colleagues, 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2564285/)

The proposed implementation is a bounded local state machine for the two coding
ends at each join. Represent both oriented strands, base pairing, nicks/gaps,
terminus chemistry and original-coordinate ancestry. Specify permitted
transitions for opening, strand-specific resection, TdT extension,
annealing, template-dependent synthesis and ligation. Each transition must
validate its substrate and preserve a replayable history. A polymerase extension
requires an appropriate primer/template configuration; unresolved gaps or
mismatches cannot silently become a completed sequence.

Keep nick offset, overhang length, physical resection, final net deletion,
retained P and other templated synthesis as separate quantities. Treat the two
D ends and successive joins consistently under inversion. Ensure one model owns
junction additions so duplex-generated N/P cannot be added again by legacy
passes.

Extend ancestry to represent shared templates and both contributors to an
annealed overlap. A single germline coordinate or a P/N flag cannot express
every such history. Project an AIRR sequence only after the modeled product is
resolved; preserve richer ancestry in versioned GenAIRR metadata.

Start with deterministic opening and paired-end primer/template completion
fixtures, then add erosion and joining, then bounded microhomology/polymerase
alternatives. An isolated 3′ overhang must not acquire a fictitious fill-in from
its recessed 5′ terminus. Validate short substrates
exhaustively, orientation symmetry, strand connectivity, exact replay and
independently specified experimental examples. The same pure transition kernel
must serve feasibility, runtime and replay. Synthetic transition weights are
declared model assumptions, not measured enzyme rates.

## 7. Calibrate identifiable outcomes before microscopic rates

Introduce a native `CodingEndOutcomeModel` carrying joint retained outcomes,
conditioning variables, support, model version, compatible reference hashes and
training provenance. Expose joint D rows in Python as well as Rust; do not force
independence where data support correlation.

First support supervised estimation from records with genuine generating
history, using simulated truth for recovery tests. Biological rearrangements
require a different likelihood that sums compatible allele, trim, P and N
explanations and accounts for read errors, clipping and mutation. A palindromic
sequence match alone is not observed biochemical P provenance.

Current NP estimators accepting missing P fields does not establish that the
input contains only N. Fit consistent joint outcomes rather than separately
counting P as N and treating observed net trim as physical strand erosion.
Select and document training cohorts, genotype information, clone weighting,
SHM and sequencing uncertainty. Nonproductive data can reduce some selection
effects but do not eliminate all ascertainment bias.

Declare the target probability law before estimation. Current productive
recombination applies local feasibility conditioning; its output is not the
unconstrained generative model globally conditioned on productivity. The
likelihood, simulation mode and training selection model must agree.

Release criteria are synthetic parameter recovery where identifiable, held-out
donor/clone predictive checks, calibration uncertainty, sensitivity to annotation
and selection assumptions, and explicit reporting of indistinguishable
parameters. Keep unsupported parameters fixed or expose ranges instead of
reporting spuriously precise estimates.

Fitting detailed opening/resection/polymerase rates is a later milestone.
Final junction sequences alone generally cannot distinguish all microscopic
histories. Intermediate-product or controlled perturbation data are needed
before claiming biochemical calibration of those mechanisms.

## Implementation sequence

| Slice | Deliverable | Dependency and completion criterion |
| --- | --- | --- |
| 1 | Permanent reproducers, precise errors and replay input checks | Reproduce the reported defects; distinguish exhaustion from proof; reject malformed traces |
| 2 | Shared productivity policy, nullable assessment and final guarantees | Independent truth fixtures agree across all consumers; document changed unknown/coverage behavior |
| 3 | Native annotation resolver and reference evidence transport | No silent anchor-resolution fallback; preserve identity, functionality and independently known phase |
| 4 | Versioned natural-law indel solver and edit-view evaluator | Exact counts 0–2, sound higher-count strategy, correct weighted oracle and strict precursor regressions |
| 5 | Audited cartridge releases and integrated engine validation | Source-backed corrections or explicit unresolved status; full suites, expanded sweep and benchmarks |
| 6 | Native retained-outcome estimation | Probability law and input provenance fixed; recovery and held-out evaluation pass |
| 7 | Experimental duplex transitions and richer ancestry | Independent substrate fixtures, topology invariants and replay pass |
| 8 | Empirical duplex calibration | Appropriate data identify the fitted quantities; uncertainty and limitations are reported |

Slices 2 and 3 can progress in parallel after the behavioral fixtures are pinned.
Slice 4 can build its sampler infrastructure concurrently, but its final
predicate must use slice 2. Retained-outcome estimation does not need to wait
for duplex simulation. Duplex rates remain a research deliverable with separate
evidence requirements.

The immediate engine changes are narrow: a common biological assessment,
reference capabilities, a sound constrained indel solver, typed failure
diagnostics and strict replay validation. Retain the persistent simulation,
transaction and observer architecture. Optimize the new edit evaluator against
measured costs before considering broader engine refactoring.
