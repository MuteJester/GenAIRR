# Durable family replay design

Proposed on 2026-10-21 against `18ecdf6e9bb180712e0f61f12cdbe55706775ee5`
on `fix/biological-metadata-correctness`. This document specifies future work;
the APIs and artifact described here are not implemented.

## 1. Outcome and first-release boundary

A newly captured duplex lineage run can be saved, moved to a fresh process,
and reconstructed using only the bundle and a compatible GenAIRR installation.
The original Python program, compiled object, cartridge installation and
mutation-model pickle are unnecessary. Replay reconstructs the founder, every
attempted tree, accepted cellular histories, sampled abundances and final
observed records. It validates their generating choices and full ancestry.

All execution and validation remain native parts of GenAIRR. This work adds no
external simulator, service, executable Python payload or runtime dependency.
It needs no measured biochemical data and does not fit or change repair rates.

| First release includes | Remains separate |
| --- | --- |
| Source-backed duplex VDJ `clonal_lineage`, including multiple clones | Legacy lineage, fixed expansion and repertoire bundle adapters |
| Neutral growth, explicit targets and automatically generated targets | New selection, repair or growth models |
| Strict/permissive runs under currently audited duplex contracts | Arbitrary custom passes or contracts |
| All attempted trees, including failed survival attempts and accepted extinction | Export of partially failed runs |
| Existing supported post-lineage substitutions, masking, indels, clipping and reversal | VJ, phased genotype families, revision, contaminants and paired-end adapters |
| Exact choice-based reconstruction, including growth and sampling | Proof of the originating seed stream or a seed-rerun API |
| Founder reconstructed from an empty-state source plan | Arbitrary external initial snapshots or `run_from` founders |

The artifact is named `genairr_family_replay`, schema version 1. Its initial
workflow is `duplex_lineage_v1`. Other family types do not masquerade as trees
or receive partial support under this workflow. Their future adapters must
declare their own capabilities and compatible schema evolution.

## 2. Guarantees and scientific meaning

Three checks have distinct meanings:

1. **Structural/integrity checks:** the file is bounded, well formed and its
   typed content commitments agree. This alone does not verify a simulation.
2. **Semantic replay:** GenAIRR reconstructs all declared choices through the
   versioned source, target, growth, mutation, sampling and observation rules;
   all state, event, topology and record comparisons pass.
3. **Seed rerun:** reproducing the original PRNG stream is a different claim.
   Version 1 records seed labels for provenance but provides semantic choice
   replay, without substituting new draws or claiming seed-stream proof.

An offspring count, mutation or sampled index must have positive support in
its reconstructed context. Admissibility does not certify that the choices
were randomly sampled with the declared frequencies. Content hashes do not
authenticate authorship: someone can construct another internally consistent
realization and recompute its hashes. No signature or empirical validity is
implied.

The replayed metadata preserves the existing distinctions between generating
alleles and re-calls, repair and downstream events, physical and information
ancestry, certified placement and sequence coincidence, missing values and
zeros, cellular state and read observations. It does not independently certify
reference annotations or receptor expression/function.

## 3. Why saving the current family object is insufficient

| Existing foundation | Additional durable evidence required |
| --- | --- |
| `DuplexFamilyReplay` keeps founder/edge outcomes, S5F kernel/rate, contracts, policy and reference hash | Source execution binding, reconstructible plans/references, growth/target settings and attempt history |
| `replay_mutation_edge` validates actual mutation choices, events and ancestry without RNG | Versioned durable edge receipts and semantic state commitments |
| `grow_core` creates the tree | Ordered proposed offspring counts, visited parents, cap decisions, child seed labels and termination |
| `sample_cells` accumulates abundance by cell identity | Ordered sample draws against the derived final-live set |
| The wrapper retries survival and merges read corruption | All rejected attempts, accepted-attempt identity and unmerged observation receipts |
| Result validation retains original node addressing | A sealed output manifest and reconstructed immutable origin registry |
| `TraceFile` v4 binds choices to signatures | Portable typed plan, reference and numeric-model encodings |

Current results cannot recover discarded growth/sample/read evidence. Export
requires a new run with capture enabled. The writer must reject uncaptured
results; it must not silently rerun their seed to manufacture a history.

The existing reference content hash uses delimiter-concatenated strings.
Different free-form identity fields containing delimiters/newlines can encode
the same byte message. This is an encoding alias, not a SHA-256 collision.
The new bundle therefore binds a typed, length-framed reference encoding in
addition to the legacy hash required by an embedded v4 source trace. Historical
trace hash meanings remain unchanged. The analogous rule applies to plan
signatures: human-readable concatenated signatures are comparators, not a
portable executable specification or the new authoritative identity.

## 4. Artifact structure and executable specifications

Use one UTF-8 JSON file initially, with no archive extraction, external file
references, compression dependency or implicit network lookup. A content table
deduplicates shared references, kernels, repair models and plans. References
between entries use their typed content digests; each reference must resolve
inside the file.

The top-level payload contains:

- Artifact kind/schema and the complete required capability/version set.
- Engine/package/build provenance, a tested numeric compatibility profile,
  and recorded seed labels/seed-scheduling version.
- Ordered typed content tables and their digests.
- One execution specification and ordered per-clone receipts.
- The immutable output manifest and expected record dictionaries.
- A root commitment covering every semantic field and every table reference.

Required component versions include the plan, reference, canonical encoding,
state/event commitment, family choice vocabulary, growth, target-generation,
sampling, record projection and productivity policies. Unknown required
versions or capabilities reject. A package version string alone is insufficient
to establish execution compatibility.

### Reference and model specifications

`ReferenceSpecV1` stores the complete effective, ordered V/D/J/C catalogues,
identity, rules, allele IDs/names/genes/sequences, functionality, anchors,
subregions and sequence-bound annotation/phase/boundary evidence. Preserve pool
order: IDs, source support and re-call ties depend on it. Storing only the
generating alleles would change both generation and inference.

Keep reference admission options and annotation/productivity policy explicit.
Loading the embedded reference does not upgrade uncertain evidence. An
admission override cannot establish biological readiness.

`KernelSpecV1` stores all 1,024 mutability and 4,096 substitution values in
their exact native order and binary64 representation. A model name is provenance,
not a lookup instruction. Include complete numeric kernels used by division or
any supported source/observation stage, even if their display names coincide.

`DuplexModelSpecV1` retains all ordered authored histories, operations, bindings,
weights, conditioning law and resource/preflight settings. It must not retain
only the selected history: the alternatives define its sampling context.
Existing unit defaults and all actual effective weights remain unchanged.

### Native plan codec

Add a closed `ReplayPlanSpecV1` codec for the concrete native passes needed by
the supported source and observation paths. Serialize normalized effective
parameters at compilation, after cartridge overlays/defaults and locks have
been resolved. Do not serialize Python class names/private dataclasses and
later reinterpret their defaults.

Each audited native pass supplies a typed replay descriptor; the default for
an unsupported/custom pass is no descriptor. A native registry reconstructs
only those closed variants through their ordinary constructors and compiler.
All distributions, ordered supports/weights, rates, count sources, selection
settings and meaningful execution/preflight limits must round-trip. No `eval`,
pickle, callbacks or arbitrary module imports are allowed.

The compile-time descriptor and the running plan receive the same immutable
execution identity. A writer cannot pair an outcome with a lookalike plan or
different kernel. Codec tests require exact descriptor round trips and
agreement with every covered pass's complete sampling-law signature. Signatures
remain useful diagnostics; the typed descriptor digest is authoritative.

Contracts have typed variants and parameters, not just names. Persist actual
strict/permissive execution policy and admission settings from the run, not a
later caller-supplied value. A continuation's policy is bound explicitly to its
duplex initial state; source-less post plans must not fall back to ordinary
policy/signature selection simply because they lack a source pass.

## 5. Canonical commitments and numeric compatibility

Define one native `typed_canonical_v1` encoder for hashing, shared by the Python
API through the engine. Use domain/type tags and length framing, ordered arrays
and maps sorted by UTF-8 key bytes. Preserve strings without Unicode
normalization. Reject duplicate map keys. Distinguish null, Boolean, integer,
float, string, array and object; `false`, `0` and `0.0` must remain distinct.

Encode scientific/configuration binary64 values as fixed-width hexadecimal bit
strings in JSON. Decode with exact bit reconstruction and validate finiteness
and the parameter's ordinary domain. This retains signed zero and tiny positive
weights without relying on a JSON client's number conversions. Fixed-width
integer IDs/counts use checked decoders; u64 seed labels use canonical decimal
strings. Boolean-to-integer coercion and narrowing casts are errors.

Dynamic values in expected records and user annotations use a lossless tagged
wire codec, not untyped JSON number guessing. Every value has an explicit tag:
null; Boolean; string; integer with a canonical decimal string in
`[-2**63, 2**64-1]`; binary64 with exactly 16 hex digits; ordered array of tagged
values; or object with ordered, unique string-key/tagged-value entries. An object
is wrapped as an object value even if its own keys resemble these codec tags.
For example, a literal hex-looking string is tagged as a string, never inferred
to be a float. Preserve `False`, `0`, `0.0`, `-0.0` and tiny finite floats through
both hashing and readback. API results decode to ordinary Python dictionaries.

Thus supported user metadata is a finite JSON subset with string keys and the
integer range above. Unsupported Python objects or out-of-range integers reject
capture explicitly; normal uncaptured runs retain their existing behavior.

`SemanticStateCommitmentV1` covers ordered nucleotide fields/flags/provenance,
regions, generating assignments/orientations/trims, inference support, actual
running counters and burden, plus the complete duplex founder and current edit/
coverage history. Shared founder evidence can be referenced by digest. Define
an explicit inventory from native fields; adding a semantic field requires a
new supported state version rather than silently omitting it.

Exclude pointer identity, allocator layout, cache version counters and transient
dirty-window logs. Rebuild live-call caches under the versioned refresh policy;
recomputed re-calls and all generated record fields still must match. Event
commitments cover ordered committed events and their stage-local choice spans,
not just final DNA or a net mutation count. Source and edge receipts store
expected state/event digests and counts; native replay regenerates full outcomes.
No unrestricted `Deserialize<Simulation>` or trusted final-pool import is needed.

Affinity uses floating-point arithmetic and `exp`; universal bit-identical
replay across every platform is not promised. Version 1 declares a tested
numeric profile and checks actual numerical commitments exactly. Unsupported
profiles or changed results raise a numerical compatibility error. Do not use
tolerances to relabel a different history as exact replay. Add a deterministic
portable math implementation later if wider compatibility requires it.

## 6. Execution receipts

IDs are scoped by `(clone_id, attempt_id, node_id)`; observation IDs additionally
identify their read stage. Dense `parent_id` indices are separate from clone
IDs. No receiver infers scope from equal DNA, a shared repair ID or matching
integer node IDs in different trees.

### One source receipt per clone

Capture the source-bound v4 trace and its immutable execution identity while
the result of `_pre.run` still carries that binding, before native family
construction copies its raw `Outcome`. Bind this receipt to the typed source
plan/reference/model, actual policy and the empty initial-state commitment.
Replay the source exactly once per clone. All attempts refer to that same
founder state/event history; survival retries never draw a replacement repair.

Direct native families without this source evidence cannot be exported by the
first version. Merged outcomes must continue to fail the ordinary v4 writer's
source-binding checks. The new format does not weaken those gates.

### Every growth attempt

An `AttemptReceiptV1` contains its ordinal, recorded seed label, target mode,
target choices/resolved target, ordered growth choices, child edge receipts,
termination summary and ordered sampling choices. Keep rejected candidates
before the wrapper discards them. The accepted-attempt index is a derived
check, never a substitute for earlier attempts.

Targets follow the current rules:

- Neutral/no-target mode has no target draws or invented affinity evidence.
- An explicit target is bound as configuration, with no target-generation draws.
- An automatic target records ordered positions and replacement residues,
  including repeated positions and reversions. Validate each step against the
  evolving target and its actual destination support. Preserve the empty-input
  and zero-substitution cases without fabricated choices.

Currently an automatic target is generated anew from each attempt's family
seed. Preserve that attempt-local target law. A fixed target across retries
would be a separate simulation-policy change, not a replay implementation fix.

For each generation, replay visits live parents in their actual order and
consumes a `proposed_offspring` count, including zero. Derive live-set size,
logistic damping, fitness and effective offspring rate from the reconstructed
state. Require finite, nonnegative rates and admissible counts: zero rate
permits only zero offspring; a positive rate permits the declared Poisson
support within integer/resource limits. This checks semantic support, not
originating PRNG words or exact finite-PRNG probability mass.

For every realized child, consume its child-seed label and addressed division
trace, derive its ID/parent/generation, and replay mutation before computing
that child's affinity. Preserve actual cumulative events and burdens, including
repeated mutations and reversions. Child seed labels are retained provenance;
division replay consumes choices, not that seed's RNG stream.

Capacity follows the current loop exactly. Offspring are proposed before the
inner capacity check. If a parent exactly fills capacity, later parents can
still have zero proposals before a positive proposal triggers the stop.
Record only actually visited parents and instantiated children; cap-suppressed
offspring are not births. Unvisited parents receive no invented draws.
Retain current generation replacement semantics, including zero next-generation
growth at full logistic capacity. Derive final-live IDs and the stopping reason
instead of trusting a saved tree or treating every leaf as alive.

The first capture profile remains in the current Knuth sampler's small-rate
domain: conservatively require finite
`lambda_base * (1 + selection_strength) <= 64` before capture. This bounds all
effective rates for the current nonnegative-distance affinity model and avoids
the large-rate threshold-underflow regime. It is a numerical/computational
guard, not a biochemical prior. Check realized rates too. Reject unsupported
capture configurations explicitly; do not clamp rates, replace the sampler or
change uncaptured runs. Wider-rate support needs its own versioned sampler work.

### Sampling, survival and observed reads

Each sampling receipt records an ordered index into the derived final-live
array, with the selected node ID as a redundant check. Draw count and eligibility
follow the existing with-replacement rule. Extinct attempts have no sample
draws. Derive cell abundance and observation flags; distinct equal-DNA cells
remain distinct. Repeated draws of one cell aggregate that cell's abundance.

Record the exact survival policy, retry bound and schedule. With the existing
default policy, all earlier attempts must fail its acceptance predicate and
the accepted attempt must be the first that passes. With `allow_extinction`,
accept the first attempt, including extinction. Preserve an accepted extinct
founder/tree/family position while emitting zero records. Exhausted or aborted
runs do not produce an apparently complete result bundle. DSL validation still
requires positive generation, capacity and sample controls.

Sampling draw order and emitted record order are different. Current records
are emitted in observed-node arena order. The post-read plan executes once
per observed cell after repeated sample draws have been aggregated, not once
per sample draw. Preserve that behavior and its original sequence-ID rules.

Capture each unmerged post-read continuation as a `StageReceiptV1`, with the
actual input node commitment, typed plan/context, local addressed choices,
diagnostics and expected state/event commitments. Replay from that cell's
reconstructed state and merge through the ordinary checked continuation path.
No trace is inferred by subtracting common prefixes from merged outcomes.
Pristine reads explicitly use an identity observation stage, with no fabricated
executed passes. Stage-local addresses are namespaced by the enclosing receipt,
so repeated addresses in distinct cells cannot consume each other's choices.

## 7. Capture, sealing and public API

Proposed opt-in API:

```python
from GenAIRR.experimental import FamilyReplayBundle

# experiment is a supported source-backed duplex clonal_lineage experiment.
result = experiment.run_records(
    seed=42, strict=True, validate_records=True,
    capture_family_replay=True,
)
result.save_family_replay("family.gafamily.json")

# A later process needs no original experiment, cartridge or kernel loader.
bundle = FamilyReplayBundle.read_from("family.gafamily.json")
restored = bundle.replay_records()
cell = bundle.replay_node(clone_id=0, node_id=3)
```

`FamilyReplayBundle.read_from(path, limits=...)` and replay accept caller-owned
`FamilyReplayLimits`; captured runs accept the corresponding
`family_replay_limits` option. `save_family_replay(path, replace=False)` exposes
the explicit file-replacement choice. These optional arguments must not change
the recorded stochastic configuration or relax semantic verification.

All names above are proposals. `replay_records` always performs full validation
before returning `SimulationResultWithLineages`, with native parents, cellular
families, trees, merged observed outcomes and the original record order. Rebuild
the private node-origin registry from verified scopes and native trees rather
than trusting mutable saved record tags. The loaded result supports the existing
record/family/parent validators and ordinary exports.

`replay_node` returns an accepted attempt's cellular state before observation,
including for root/unsampled nodes. It requires the full bundle to have passed
validation first. The existing in-memory `families[i].replay_node(node, refdata)`
continues to work with its current contract. Accepted extinct families retain
their root even though they have no observed records.

Replay accepts resource ceilings, not seed, source-weight, policy, metadata or
sampling overrides. A caller may provide expected reference identity for an
additional exact typed-match gate; it is not needed for self-contained recovery.
Inspecting a parsed manifest must not return a trusted native family or claim
that semantic replay passed. Cache only verified immutable evidence; mutable
returned records cannot become the cache's new ground truth.

Capture defaults to off. Before drawing a founder, capture preflight checks
workflow support, all plan/reference/model codecs, metadata representability,
policy/numeric profile, limits and conservative growth bounds. Unsupported
fixed/repertoire/native-external-seed/legacy/genotype requests reject before
sampling. Capture observes existing draws without adding any, reordering them
or changing conditional support; failures return structured errors without
retrying a different biological realization.

The compiled and `Experiment` entry points share receipt capture. Preserve their
current metadata semantics: the direct compiled path does not suddenly acquire
the top-level `Experiment` metadata stamps. A private entry-point-aware capture
context carries the actually applicable annotations and a deferred-finalization
flag. An `Experiment` call validates its intended annotations before sampling,
delegates capture with sealing deferred, applies its metadata, and seals exactly
once. A direct compiled call seals its actual output with no top-level stamps.
Metadata declared on a compiled wrapper but unused by that entry point is not
part of its executable bundle and does not cause capture rejection merely by
being non-JSON. Capture records which annotations were actually applied.

The finalization step deep-freezes/canonically encodes generated records,
wrapper fields and actual JSON annotations; it must not retain shallow mutable
aliases. It never reopens an already sealed envelope to accommodate later
metadata. Normal result dictionaries remain mutable, and subsequent nested
changes cannot alter the frozen expected evidence.

Saving validates the actual mutable result against this immutable envelope
and its native evidence before taking the final file snapshot. Removed rows,
changed generated fields, altered addressing or changed annotations reject;
they are not blessed by recomputing a digest. Users who need annotations should
declare them before the captured run. Manually constructed records-only results,
old uncaptured results and mutated private evidence do not acquire a generating
history through export. Capture may perform extra validation even when ordinary
`validate_records` is false; it must never alter sampled choices.

## 8. Bounded loading and replay

Implement the file reader separately from the existing whole-string v4 reader.
Use a byte-limited reader and bounded typed visitors before allocating tables
or arenas. Reject duplicate keys/IDs, unknown executable variants, invalid
types/nonfinite values, count/size overflows, missing or surplus receipts,
unresolved component references and cycles. Announced counts are cross-checks,
not permission to allocate. No partial native result is published.

`FamilyReplayLimits` is caller owned. Initial software defaults are a 256 MiB
file limit, depth 64, one million total nodes, five million total choices,
fifty million cumulative replayed nucleotide instances and ten million repair/
replay work units. Also start with five million cumulative allocated/copied
history entries and a 256 MiB ceiling on accounted retained provenance storage.
Enforce separate reference, model-operation, string/metadata and attempt counts
within those global ceilings. The total node/choice/base
limits include rejected attempts, source stages and observed-read stages, not
only returned records. Checked arithmetic applies before multiplication or
capacity reservation. These are overridable computational ceilings, not
biological limits or wall-clock guarantees; benchmark them before release and
document changes independently of scientific defaults.

Specify accounting units in the codec implementation before release: base work
charges source construction and each stage's input/output length, plus any
intermediate materialization; work units charge model operations, candidate/
search steps, visited parents, target/sample decisions and fresh capture RNG
words. In particular, a small Poisson rate is not a worst-case word-count bound.
Capture must fence that sampler's actual loop work without extra RNG draws,
and fail the whole captured run on exhaustion. Replay consumes semantic counts
and never reruns that sampling loop. Internal sampler budgets remain effective
and cannot exceed the caller's overall ceiling.

History budgets charge cumulative edits, removed links, current links, template/
certificate data and event history before copying/extending them, including
candidate/temporary states and rejected attempts. Charge hashing/validation of
history prefixes to work too. Accounted retained storage includes container
capacity under the qualified native layout; these deterministic software
ceilings do not assert a hard process-RSS cap or eliminate allocator overhead.
Checking file size or nucleotide counts alone is insufficient.

Do not let bundle-provided preflight/work settings increase the reader's local
ceilings. If the original specification needs more resources, reject with a
resource error; never truncate the history, reinterpret budget exhaustion as
empty support/extinction, or return a partial weighted sample.

Full replay proceeds as follows:

1. Check format, required capabilities, typed commitments and local limits.
2. Reconstruct full references/models/plans through audited constructors,
   perform ordinary compile/admission/support checks and verify their identities.
3. Replay each clone's source from empty state through its v4 choices, with the
   additional new typed bindings. Verify founder state/event commitments.
4. For every attempt in order, replay the target and shared growth control flow.
   Consume actual visited-parent/child decisions and replay each division once
   before later affinity-dependent decisions. Derive topology and final-live
   membership. Verify all intermediate and final commitments.
5. Consume ordered sample indices, derive observation flags/abundances, and
   enforce the retry/acceptance policy across all attempts.
6. Replay saved post-read stages from each accepted observed cellular state;
   merge source/division/observation outcomes using checked boundaries and
   stage-local event-span rebasing.
7. Recompute native records and wrapper metadata, apply exactly the captured
   declared annotations once, and compare every expected field/value/type and
   row position. Rebuild native family/result/origin objects and run all record,
   family and parent validators.
8. Require every receipt/choice cursor to be consumed exactly. Publish the
   verified result only after every clone and output passes.

Use shared record/consume interfaces for normal target/growth/sampling control
flow; do not maintain an independent approximate replay growth algorithm.
Existing mutation/pass replay remains the state-transition engine. Reject all
fresh RNG use during semantic replay, including hidden fallback draws. Keep a
zero-word fence for ordinary pass replay and an equivalent no-generator fence
around the family decision provider. Stored seed labels never authorize a
fallback draw. Sampling diagnostics are nonbiological evidence and retain their
declared interpretation instead of becoming events.

Version 1 replays complete attempted trees once, discarding rejected attempt
arenas after validation and retaining accepted cellular state as needed. The
serialized format stores one edge receipt per division and shared immutable
founder/model data; it does not duplicate full ancestral paths or serialize
every IR revision. Its serialized size should scale with nodes, recorded
edits/choices and the complete exported record payload. This is an artifact-size
target, not a claim of linear replay memory: current duplex states copy
cumulative history vectors on branch mutation, and edge outcomes retain
parent/child states. Enforce the history-copy/retention budgets above. Persistent
prefix sharing or lazy materialization can reduce those costs later without
changing semantic identity. Profile real-length and deep-tree references as
well as short controls before claiming practical scaling.

### Failure and file-publication rules

Expose structured errors with component/clone/attempt/node/stage paths, with
distinct categories for unsupported schema/capability/codec, missing capture,
reference/plan/model mismatch, invalid choice/order/scope, state/event/record
mismatch, numerical incompatibility, resource exhaustion and I/O failure.
Unsupported evidence is distinct from proven impossible biological support.
An internally valid file does not pass semantic validation merely because its
root hash matches.

Write through a same-directory temporary file, flush/sync and close it, then
publish atomically. Default publication must use a race-safe no-clobber operation,
not an existence check followed by an overwriting rename. Explicit replacement
uses an atomic replacement path. Sync the containing directory on qualified
platforms where required for crash persistence. Failure cleans the temporary
file and leaves the previous artifact intact. A failed or cancelled run cannot
leave an apparently complete family file. No multi-file transaction or archive
path handling is needed.

## 9. Compatibility and extension boundaries

Existing `TraceFile` schemas 1–4, signatures, address vocabulary and ordinary
source/continuation writer gates retain their meanings. The family artifact
has its own required family-choice vocabulary; it nests existing addressed
source/pass traces without renaming their local addresses. New component
semantics or required wire fields require the corresponding version change.
Unknown semantics fail explicitly rather than degrading to a legacy comparator.

Old uncaptured results cannot be upgraded to this bundle. A records-only JSON
file, a saved tree, and a merged ordinary trace are not alternate valid inputs.
No arbitrary IR object graph is accepted as an external initial snapshot.

Later external-snapshot support needs a separately versioned snapshot schema,
an explicit evidence/validation boundary and a declaration of which pre-snapshot
history is asserted versus reproducible. It must not claim source generation
that was not recorded. Source-backed families do not require that work first.

Later fixed expansion/repertoire adapters need their own control receipts.
Repertoire collapse must retain every raw read's observation receipt and ordered
collapse membership: identical full records may still arise from distinct
no-op choice histories. Its no-post bulk-copy shortcut records multiplicity
without inventing executed passes. Legacy lineage needs a decision about its
current synthesized net-mutation history before receiving an actual-event
replay claim. None is silently enabled by this format's first release.

## 10. Verification and release gates

Use independent literal oracles as well as round trips. Comparing a writer and
reader that share the same faulty helper is insufficient.

| Area | Required checks |
| --- | --- |
| Self-contained execution | Save, end the process, load/replay in a fresh subprocess with named cartridge and bundled kernel loaders disabled; no original experiment/compiled object; exact records and native validation |
| Source truth | Both D orientations/join orders, same-DNA distinct histories, unequal weights, productivity conditioning, restored landmarks, source-reference mismatch and missing evidence |
| Actual mutations | Repeated hits/reversions; literal edge history and cumulative counts/burden, independent sibling maps, roots and unsampled cells |
| Growth | Hand-authored neutral and selected examples, zero offspring, generation termination, cap-suppressed proposals, capacity filled exactly followed by zero/positive proposals, true final-live eligibility |
| Targets/retries | Explicit/automatic/inactive targets; repeated target positions; attempt-local target changes; first accepted attempt; all earlier failed attempts; accepted extinction and aligned empty family positions |
| Observations | Ordered with-replacement samples, equal-DNA distinct cells, repeated-cell abundance, arena-order output, one post stage per observed cell, pristine identity stages |
| Read truth | PCR/quality substitutions, N masking, indels, terminal loss and reversal; literal current/removed ancestry and correct pre-loss mutation denominator; null fields survive ordinary exports |
| Output/sealing | Exact generated types; lossless tagged values including Boolean/int/float distinctions, signed zero and hex-looking strings; deferred Experiment versus direct-compiled sealing; unused non-JSON compiled metadata; nested-object mutation cannot change the envelope and causes rejection if it changes actual records; removed/reordered/edited rows and cross-clone/node namespace swaps |
| Strong identity | Same-name different S5F arrays; changed source weights/model/locks; same-size different references; changed annotation/phase/rules; delimiter/newline reference-identity aliases distinguished by the new typed digest |
| Hostile/malformed input | Duplicate keys/IDs, future versions, extra/missing/reordered stages, cyclic/dangling references, Boolean IDs, overflow, nonfinite or malformed bit encodings, truncated choices and digests recomputed over semantically invalid payloads |
| Resource behavior | Byte/depth/node/choice/base/history/work limits fail before excess charged allocation/execution; deep-tree copied ancestor prefixes, temporary states and rejected attempts charged; interrupted/failed atomic writes preserve the previous target |
| Numerical/compatibility | Qualified numeric profile fixtures and exact mismatch errors; current v1–v4 trace goldens unchanged; capture-on/off seeded choices and records equal |

Tamper tests must exercise semantic checks after recomputing all content hashes,
not only corruption detected by the root digest. Valid alternative realizations
are not universally distinguishable from originals without external authenticity
evidence; tests must not promise that.

Run full native debug/release and Python correctness suites, existing unchanged
performance budgets, strict docs and executable examples at the final binary.
Benchmark capture overhead, file size, parsing peak memory and full replay
separately across short and realistic references, tree sizes, read corruption
and survival retries. Record source/binary hashes, qualified numeric profile,
commands and independent sweep counts in a new dated verification artifact.

## 11. Implementation sequence

| Slice | Deliverable | Gate |
| --- | --- | --- |
| F1 | Typed reference/kernel/model/plan codecs, canonical encoding, state/event commitments and capability checks | Exact normalized round trips; full semantic field inventory; delimiter alias and same-name kernel tests; supported source/read plans reconstruct identically |
| F2 | Native target/growth/sampling record/consume providers and attempt receipts | Independent tiny-family/cap/target/sample oracles; unchanged seeded generation; zero fresh RNG on replay |
| F3 | Source-bound wrapper capture, all survival attempts, unmerged read receipts and sealed output manifest | One founder per clone; full retry evidence; exact post-read records and metadata/origin validation; unsupported capture fails before sampling |
| F4 | Bounded versioned bundle reader/writer and complete public replay APIs | Fresh-process self-contained recovery, exact native result reconstruction, semantic tamper and resource/error tests, atomic writes |
| F5 | End-to-end audit, compatibility, real-length benchmarks and guide | All release gates and a committed source/binary-pinned verification report |

Intermediate slices are internal foundations, not a released whole-family
replay feature. Complete every source, attempt, cell and read stage before
advertising the first workflow as durable. Keep later workflow adapters,
arbitrary snapshots, seed reruns and additional biochemistry explicitly outside
this implementation's completion criterion.

## 12. Code anchors

- [Ordinary durable trace format and hashes](../engine_rs/src/trace_file.rs)
- [Source binding and trace writer gates](../engine_rs/src/python/compiled.rs)
- [Current in-memory family evidence and replay](../engine_rs/src/python/lineage.rs)
- [Generation loop and mutation-edge replay](../engine_rs/src/lineage/branching.rs)
- [Family growth/sampling orchestration](../engine_rs/src/lineage/family.rs)
- [Observed-cell sampling](../engine_rs/src/lineage/sampling.rs)
- [Wrapper retry/read/record integration](../src/GenAIRR/_compiled/lineage.py)
- [Checked outcome composition](../engine_rs/src/pass/outcome.rs)
- [Current result/origin registry](../src/GenAIRR/result.py)
- [Completed family milestone](duplex_family_implementation_2026_10_20.md)
