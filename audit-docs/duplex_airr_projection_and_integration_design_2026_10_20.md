# Duplex AIRR projection and experiment integration design

Proposed on 2026-10-20 for `fix/biological-metadata-correctness`, following
the coupled repair checkpoint `c5ad816`. This is a design for the next
implementation, not a report of implemented functionality.

The next milestone is a reference-bound projection of repaired duplexes into
trustworthy sequence metadata, followed by opt-in `Experiment` integration.
The repair product must remain the source of truth. Fitting it into today's
single-segment sequence representation must never erase conflicting ancestry,
manufacture a germline deletion or silently change the sampled history.

All repair, projection, sampling, estimation and validation remain inside
GenAIRR. This design adds no external simulator, annotation service or runtime
dependency. Authored weights are model controls; measured biochemical rates
will require later data. Synthetic recovery remains the calibration target for
now.

## Existing foundations and the missing contract

The [coupled implementation](coupled_duplex_implementation_2026_10_20.md)
already processes V-D and D-J on one shared D duplex, in either authored order,
including inverted D. It preserves original physical nucleotide identities,
template links, removed material, origin anchors and ligation seams. It does
not yet produce standard AIRR records or participate in normal experiments.

Three facts must stay distinct:

| Fact | Example | Representation |
| --- | --- | --- |
| Physical ancestry | This top nucleotide was synthesized; its predecessor was resected. | Immutable repair ledger and operation events |
| Information ancestry | Its template ultimately descends from D reference position 4. | Template graph, including removed templates |
| Reference placement | It restores D position 4 in the D body, or lies elsewhere after a hairpin transfer. | Certified alignment placement, separate from physical ancestry |

A template coordinate alone cannot establish placement. A sequence match
cannot establish any of these facts.

The AIRR specification allows missing field values and custom extensions;
its exported interval coordinates are one based and closed. Its NP fields
describe the combined N/P intervals between gene alignments. Those rules
permit an honest partial record, but do not define how a simulator should
assign causal duplex ancestry. The policies below are GenAIRR conventions.
[AIRR Rearrangement specification](https://docs.airr-community.org/en/stable/datarep/rearrangements.html)

## Decisions for the first implementation

1. Use the completed product's coding-oriented top strand as the query.
   Preserve both physical strand contributions. Reverse-complementing a read
   is an observation transform and does not change this founder convention.
2. Derive alignment placement from replayable body certificates, considering
   both paired strands. Chemical origin and alignment role remain separate.
3. Retain every supported body claim. A column shared by V and D must not be
   assigned to one of them by a tie breaker or a best sequence match.
4. Return partial projections with explicit reasons when a scalar field is
   not determined. An explicitly empty tract has length zero; an unresolved
   tract has a missing length.
5. Convert into the current execution IR only when its single-segment layout
   can represent the product faithfully. Check this capability for the whole
   authored support before running an experiment. Never resample to hide a
   conversion failure.
6. Add a named duplex productivity policy for certified restored landmarks.
   Keep `canonical_observed_vdj_v1` unchanged.
7. Integrate finite, authored, reference-bound histories first. Do not invent
   enzyme-rate defaults or make arbitrary Python callbacks part of replay.

## Reference binding and coordinate systems

A proposed `DuplexReferenceBindings` value binds each source substrate to:

- Reference content identity, generating allele ID and canonical sequence
  checksum, species/locus and source instance ID.
- The exact supplied coding interval within that allele, with canonical
  forward coordinates and boundary evidence or an explicit boundary assertion.
- D orientation, landmark and coding-phase evidence, including their existing
  readiness status.

Matching literal sequence strings to an allele catalogue is insufficient.
Identical alleles remain different generating choices. Standalone products
made from raw strings can retain substrate labels but cannot acquire allele
calls without explicit bindings. A binding mismatch is an error, not an
uncertain biological annotation.

Known partial coding boundaries cannot be overridden as complete. Preserve
the existing distinction between an asserted supplied coding end and a source
verified boundary. Productivity additionally requires its own trustworthy
landmark and phase evidence; boundary readiness does not establish either.

Maintain three separate coordinate spaces:

| Space | Convention | Purpose |
| --- | --- | --- |
| Repair axis | Signed product coordinates; original oriented origins remain stable | Physical operations and seams |
| Canonical allele | Zero based, forward allele coordinates | Source identity, body placement and D provenance |
| Current coding query | Zero based, half open intervals | Python records, observation links and ordinary projection |

Convert to AIRR one based closed intervals only at strict export. For inverted
D, canonical coordinates remain forward-reference coordinates while alignment
traversal and reference bases are reversed and complemented. Preserve the
existing explicit D orientation metadata; do not reinterpret D inversion as
whole-read `rev_comp`. Strict CIGAR prefixes and reference offsets must be
derived from the same oriented path, including for asymmetric D intervals.

## Placement certificates and operation evidence

Extend the shared native executor with a deterministic annotation observer.
Each operation receives a stable `(join ordinal, operation ordinal)` ID.
Record nucleotide creation, removal, transfer, template use and physical tip
alongside existing IDs and coordinates. The observer must see both joins in
one execution; it must not reconstruct the second join from fresh germline
material.

Each certificate relates a reference interval, a physical strand placement
and the operations that justify it. It records source instance, segment,
canonical position, current repair-axis position and supporting nucleotide or
template links. The rules are:

| Transition | Certificate behavior |
| --- | --- |
| Initialize a source duplex | Each physical strand receives body certificates for its canonical reference positions. |
| Invert D or translate a component | Transform physical placement and strand; retain canonical identity. |
| Resect | Remove that nucleotide's live placement; preserve its historical certificate and removal event. |
| Transfer around a hairpin | Retain physical/template ancestry, revoke its previous body placement, and record a hairpin-transfer placement at the new tip. |
| Add N | Record a new untemplated nucleotide and tip; create no body certificate. |
| Fill from a body template in its native column | Certify restoration at that same body position; retain synthesized physical kind and the template edge. |
| Fill from transferred or N material | Preserve that nonbody placement ancestry; a distant reference ancestor does not create a body certificate. |
| Pair and ligate | Preserve every supported body claim in the column, including claims from different segments. |

Use the source origin plus oriented reference offset as a geometric invariant,
not as sufficient evidence: a coincidental return to that coordinate does not
erase a transfer history. Already established certificates may support later
in-place synthesis through multiple template generations. Removal of a
template after copying does not invalidate that historical evidence.

The query is always top-oriented, but a retained complementary body strand can
justify the column's reference placement. For example, top V-derived hairpin
material paired with bottom D-body can carry a D body alignment certificate
while retaining V-derived physical ancestry. This is a declared reference
information convention, not a claim that the top nucleotide was an original
D nucleotide. By contrast, V-body opposite D-body supplies two body claims;
both must survive projection.

Version the evidence envelope separately from the physical repair model and
the projection policy. Existing serialized coupled scripts remain accepted.
Old product snapshots without operation evidence require verified replay with
their exact substrates and histories before projection. Do not silently
reinterpret incomplete snapshots under a new version.

## Projection results and representability

Introduce a pure native projection boundary. Names below are proposed API
shapes, not existing public methods:

```rust
project_coupled(product, bindings, policy) -> Result<DuplexProjection, BindingOrIntegrityError>
try_into_simulation_seed(projection) -> Result<DuplexSeed, UnsupportedLayout>
```

`DuplexProjection` contains the query, per-segment certified alignment paths,
landmark evidence, per-field availability and reasons, a versioned rich repair
extension, and a capability result for installation into the current IR.
Binding corruption is an error. A valid molecular product with incomplete or
overlapping conventional annotations is a successful, explicitly partial
projection.

Projection first collects body certificates without choosing an owner. It
deduplicates agreeing support by exact source instance, allele and canonical
position: top and bottom evidence for the same position is one claim, not two
body contributions. Different canonical positions of the same allele remain
distinct claims, just as positions from different segments do.

It then constructs a path for each generating allele. Each individual path
must have in-bounds, unique, orientation-consistent reference placement, with
at most one canonical position per query column in that path. Every
aligned query column and every reference gap must be accounted for. Repeated
or nonmonotone placements remain explicit graph evidence rather than being
collapsed into a misleading CIGAR.

For a conventional combined alignment and a normal IR seed, additionally
require an ordered V, NP1, D, NP2, J layout, no overlapping body ownership and
no segment re-entry that the current region structure cannot express. Every
nonbody column must have an unambiguous geometric interval or a supported
internal insertion path. Unsupported layouts receive specific reason codes,
not a biological-impossibility label.

Independent per-gene alignments may remain available when their combined
layout cannot be flattened. In that case preserve valid per-gene coordinates
and CIGARs, but leave the global alignment pair and its alignment-relative
coordinates missing. Do not return unequal or mutually inconsistent global
alignment strings. This restriction reflects this GenAIRR representation;
it is not a claim that AIRR forbids overlapping gene annotations.

### Field semantics

| Field group | Duplex rule |
| --- | --- |
| `v_call`, `d_call`, `j_call` | The explicitly bound generating alleles. Keep `recalled_*_call` separate. |
| Segment query/reference coordinates and CIGAR | Derived from that allele's certified path, never from a re-call or a motif search. Missing when the path is not representable. |
| `sequence_alignment`, `germline_alignment` | A consistent combined projection when the layout is supported; otherwise both missing. |
| Segment identity | The existing declared mismatch/gap calculation over an available certified alignment. Synthesis does not itself count as a mismatch. |
| Net terminal trims | Missing terminal body positions at repair completion under the placement policy, in canonical allele orientation. Not physical resection totals. |
| `np1`, `np2` and lengths | Exact intervals between supported body paths, measured in query columns. A known adjacent boundary gives length zero. |
| Existing P length scalars | Historical per-end P emission/assembly counts. Missing unless a lossless mapping to that exact existing contract is established. Final retained tract length alone is insufficient. |
| `junction`, `junction_aa` | Actual observed sequence and translation between trusted canonical landmark placements, including nonproductive products. Unresolved landmarks leave the junction missing. |
| Mutation and indel counters | Existing downstream event definitions. Repair synthesis is not SHM; physical repair resection is not a later polymerase indel. |

Keep the D generating call if D was sampled but no D body survives repair
completion. Report its contribution as absent and unsupported alignment fields
as missing. This differs from a VJ model in which no D was sampled. If no
certified boundary separates two merged junction tracts, expose the combined
nonbody tract in the extension and leave the NP1/NP2 split missing. Do not
invent a D boundary from the ligation seam.

Later loss of a previously certified D body is a different case. Preserve its
founder reference path and any still-determined junction boundaries. Complete
structural deletion can retain a certified deletion-only path; export each
field according to its own available evidence, without fabricating a nonempty
query interval for zero query consumption. Clipping narrows observation
coverage rather than becoming that deletion. Neither event automatically
makes a known NP split unknown or changes founder trims.

Overlapping segment paths have an overlap interval, not a negative NP length.
The standalone projection may report the known intervening gap as empty when
both boundaries are known, together with explicit overlap evidence. It must
not claim that the product had no N/P processing. More generally NP intervals
are geometry; chemical N/P counts are separate history measurements.

For trims, calculate each end from an available terminal body envelope. If no
body survives, do not distribute complete loss arbitrarily between two ends.
Internal missing body positions require their own alignment gap evidence and
must not be absorbed into terminal trims. Later clipping changes observed
reference coverage, not these founder net trims. A canonical interval omitted
from the supplied substrate is source-coverage loss, not evidence of repair
trimming; bindings must distinguish it from a complete end actually processed
by the history.

### Nullable schema changes

The current native record has nonnullable trim, NP-length and P-length
scalars, and string fields commonly use an empty string. Before exposing
partial duplex records, add explicit optional values for affected fields and
audit their Rust, bridge, validator and Python consumers. `None` means
unresolved or unavailable; `""` with length zero means a known empty tract
where applicable. Strict TSV represents missing values as empty cells.
Use nullable integer columns where DataFrame conversion needs them.

Legacy and original-end records must retain their existing determined values,
record shape and seeded output. Only the new duplex mode may introduce these
new unknown cases. Provide structured availability reasons so that an absent
D, unresolved boundary and unsupported combined layout are distinguishable.

Carry availability into `DuplexSeed` and persistent processing state as well
as the output record. The current assignment P/trim integers cannot express
unknown values. Introduce a typed processing-metrics view that distinguishes
legacy determined totals from duplex evidence, and migrate every consumer of
those totals to it. Compatibility integer caches must not become evidence:
zero cannot stand in for an unknown P count or trim. Seed installation and
record validation must check this invariant.

Audit retained-outcome, P-length, NP-length/base and trim estimators explicitly.
They must reject unsupported duplex provenance or consume a documented
eligible subset with missingness diagnostics. Missing P cannot be treated as
zero, and a combined NP interval cannot be relabeled as purely untemplated N.

## Repair history and observation provenance

Add a reserved structured `duplex_repair` extension; do not relabel the existing
`coding_end_processing.model = original_coding_end_v1` payload. The proposed
extension includes:

- Schema, kernel, placement and productivity policy versions, bindings,
  authored history ID, executed join order and both ligation seams.
- Immutable nucleotide ledger, operation events and template graph, retaining
  removed nodes and both physical strand contributors per founder column.
- Body placement certificates, competing claims and projection availability
  reasons.
- Separate physical operation totals, founder query contributions after both
  joins and current observation coverage, with explicit units and strand/end
  attribution.

Count physical nucleotides per operation and physical strand; count query
contributions once per selected query column. A nucleotide may have multiple
historical transfer events without becoming multiple retained nucleotides.
Keep untemplated addition, hairpin transfer and templated synthesis as separate
axes; a synthesized mate of an N addition is not a second untemplated event.
Classify copies of nonoriginal hairpin material explicitly, rather than calling
every transferred base a canonical germline-derived P nucleotide. Do not force
ambiguous or mixed tracts into a legacy per-end P scalar.

Name the count timepoint explicitly: final founder query contributions are
measured after both joins, not immediately after the first local ligation.
Record physical tip, canonical source end and source orientation separately,
including their transformation under D inversion. A first-join emission can
be removed, moved or copied during the second join; its historical count need
not equal either the final founder or observed query-column count.

Store the detailed ancestry behind
`Simulation::duplex_provenance: Option<Arc<DuplexProvenance>>`. Separate its
immutable founder repair snapshot from the persistent current-position map:

```text
founder column -> top and bottom repair IDs, body claims, creation evidence
current query position -> founder column or downstream insertion origin
removed current link -> structural deletion, observation loss or replacement
```

After replication, repair IDs identify founder ancestry, not the same atoms
in descendant cells. Clone branches share the immutable snapshot and maintain
independent current maps. A current substitution changes observed bytes, not
the historical complementary strand. A current indel does not mutate the old
repaired duplex into a new biochemical intermediate.

| Downstream operation | Required consequence |
| --- | --- |
| SHM, PCR or quality substitution | Preserve placement and founder links; update observed bytes and the appropriate existing event counters. |
| Structural insertion | Add a new observation origin without a founder repair ID; include its actual alignment insertion. |
| Structural deletion | Preserve founder history and mark the removed link as a structural deletion, retaining the appropriate CIGAR gap. |
| Terminal loss or cropping | Mark observation loss and narrow coverage; never manufacture an alignment deletion or biochemical resection. |
| Reverse-complement output | Transform read presentation, retaining coding-oriented maps and AIRR's existing `rev_comp` convention. |
| Clone fork, checkpoint or replay | Preserve immutable ancestry and independently update descendant coverage. |

All installations and edits must go through `SimulationBuilder` and typed
events. Audit every persistent `Simulation` constructor, clone and replacement
path. A new optional sidecar that disappears in one constructor is a metadata
bug even when the emitted sequence still looks correct.

The existing `Nucleotide.segment` and `germline_pos` can cache a certified
single-path alignment role for representable seeds. They cannot be treated as
physical-original identity. Duplex-aware consumers must use the placement
evidence; physical P/N ancestry stays in the sidecar. Audit productivity,
mutation subregions/rates, indel attribution, live re-calls and all validators
for that distinction. No existing consumer may silently read a fake original
base from a newly synthesized cache entry.

## Duplex productivity policy

Add `duplex_observed_vdj_v1`, with the existing
`observed_coding_span` scope, selected explicitly by duplex mode. Its canonical
residue, phase, definite-stop and coverage criteria reuse the current policy
semantics; its landmark evidence accepts certified in-place restoration.
Do not silently broaden `canonical_observed_vdj_v1`.

A restored landmark is locatable only when all three canonical reference
positions have a unique, contiguous, orientation-consistent placement and
trusted positional landmark evidence. Coding phase is a separate requirement
for the relevant frame/stop predicates and an overall `True` assessment; its
absence does not erase a verified landmark or the junction defined by two
such landmarks. A verified noncanonical live residue can establish a policy
failure even when unrelated phase evidence is missing.

A copied Cys or Phe/Trp elsewhere, including a P tract, cannot become the
landmark. Evaluate the current query codon, not the historical template
codon. A reliable failed predicate gives `False`; unresolved required evidence
gives `None`; all predicates passing gives `True`.

Unavailable combined D alignment does not automatically make a well-supported
V/J junction or coding-span assessment unknown. Conversely, a full-length
looking query does not compensate for unresolved landmark evidence.
Reference admission overrides never certify the new placements or landmarks.

Preserve `productive_only()` as a final observed-scope promise that every
emitted record has `productive is True`. Use the same assessment for compile
feasibility, strict/permissive execution, downstream contracts and final
emission. Unknown and false completions must have distinct diagnostics even
when both are ineligible. A resource limit is not evidence of empty support.

Carry the active policy in the compiled plan, simulation and checkpoint.
Replace hard-coded policy assumptions in AIRR building, final execution
assertions and trace validation with explicit dispatch. The indel sampler's
optimized complete-bundle assessment must use a validated duplex placement
view or decline that optimization; it cannot silently apply the old
single-source view. Stage contracts still need the realized downstream state;
precomputing productive repair rows alone does not establish productive
output after mutation or clipping.

The public preservation-only `AnchorPreserved` contract has a distinct
retained-provenance and exact-reference-amino-acid promise. A synthetic
`germline_pos` cache must not make it pass by accident. Initially reject this
custom contract with duplex mode until a ledger-aware preservation adapter is
validated; restoration-aware canonical acceptance belongs to the new named
productive bundle. Other custom or subset contracts require an explicit
capability check rather than inheriting the old optimized view.

## Opt in integration and probability law

The proposed public entry point is
`Experiment.recombine(duplex_repair=bound_model)`, backed by an immutable
native-serializable model. A model contains finite history rows with stable
unique IDs, finite nonnegative weights, exact source bindings and D orientation.
It has no sequence-matching fallback. Missing coverage for any selectable
allele/orientation context is a compile error unless the user explicitly
restricts that context out of the requested model. Every selectable context
must have positive total natural weight; an all-zero context is invalid.

For a selected context `c`, normalize its authored positive weights `w_i`.
Duplicate products and different histories remain separate rows; do not
deduplicate them and change their mass. Define
`finite_authored_history_local_conditioning_v1` as:

```text
unconstrained: P(history i | c) = w_i / sum_j w_j
productive:   P(history i | c, accepted) =
                  w_i * I[assessment_i is True]
                  / sum_j(w_j * I[assessment_j is True])
```

For productive generation, retain the engine's local conditioning convention:
upstream feasible allele/orientation choices retain their natural relative
weights, and history is conditioned within the selected context. This is not
global conditioning of the whole joint generative law; that would require
weighting contexts by their completion probabilities and must be a separately
named future option. Genotype/haplotype eligibility remains an upstream rule.
Document the actual selection hierarchy and test each level's normalization.

Compile every positive row through the shared kernel and projection. Invalid
authored operations are configuration errors. Exhausted compilation work has
its own resource error. A physically valid but nonrepresentable row is a
capability error for normal `Experiment` in this first release: reject the
whole requested plan, naming the row and unsupported layout. Keep that row
available to the standalone partial projection API. Never silently drop it,
change its weight or retry a different row. This gate must run before any
productivity conditioning so that conditioning cannot conceal a layout limit.

Lower duplex recombination as one atomic assembly pass after allele/context
selection. It owns opening, resection, synthesis and both junctions. Reject
conflicting explicit legacy trim, P and NP settings; bypass default legacy
assembly/P/NP passes. Do not silently apply two junction-generation models.
Track whether these settings were explicitly supplied when constructing the
immutable DSL; comparing their values with defaults cannot establish intent.
Neither original-end retained-outcome fitting nor its parameters automatically
become a probability distribution over these microscopic histories.

Bind plan signatures and replay to reference/allele identities, sequence and
boundary evidence, landmark status and coding-phase evidence, orientation,
ordered operations, weights/history IDs, budgets, policy and capability
versions, genotype controls and conditioning law. Replay validates positive
selected-row support and re-executes the same native kernel without a new draw.
It must reject edited bindings, unsupported histories, missing
template links and malformed coordinates. Record any RNG consumption/version
changes for the new mode; preserve old-mode seeded behavior.

Serialize `duplex_repair` as a structured value in JSON and one canonical JSON
cell in TSV. Strict AIRR coordinate conversion applies only to declared AIRR
fields; nested signed repair coordinates and zero-based provenance positions
retain their explicitly declared coordinate systems.

### Capability boundaries

The first experiment slice supports reference-bound VDJ products with
representable layouts, both D orientations and join orders, followed by the
already supported substitution, indel, read-loss and read-orientation passes.
Clonal integration is released only after every continuation/export path
preserves founder ancestry and observed-scope productivity. Until then,
reject the combination at compilation rather than dropping provenance.

Receptor revision needs explicit source-instance replacement semantics for
columns/templates that cross the replaced boundary. Keep it unavailable with
duplex mode until that design and its tests are implemented. VJ needs a
two-source adapter and matching conditioning tests; do not pretend a dummy D
is equivalent. Resolved A/C/G/T substrates remain the kernel requirement;
IUPAC references cannot be silently resolved into a sampled concrete allele.

General overlapping/re-entrant execution layouts, flap or mismatch repair,
reopening protected ligation seams, RSS/RAG chronology and biochemical kinetics
remain separately named model extensions. None can be claimed from the AIRR
adapter or from synthetic recovery of authored weights.

## Concrete acceptance cases

These are geometry and annotation controls, not measured biochemical examples.
Positions in this table are zero based in the coding query/canonical allele.

| Case | Required result |
| --- | --- |
| Blunt V `ACG`, D `AAT`, J `TCC` | Query `ACGAATTCC`, exact generating calls and contiguous body paths, known empty NP intervals. |
| Same inputs, inverted D | Query `ACGATTTCC`; canonical D identities unchanged; reversed/complemented alignment traversal and strict offsets verified independently. |
| Resect D-top position 0 and fill from surviving D-bottom | Synthesized physical top nucleotide; restored D body alignment at position 0; resection retained in history, no invented final deletion. |
| Later opening transfers an earlier synthesis | Same ledger ID and template ancestry; new transferred placement, no stale body coordinate. |
| N additions from opposing ends pair in one column | One query column, distinct physical additions and template histories; no count of two query bases. |
| V `ACG`, D `GTA`, J `CCC`; VD blunt, resect V-bottom 1 and D-top 1, synapse at 2; DJ blunt at D offset 3 | Query `ACGTACCC`; column 2 supports V position 2 and D position 0. Preserve both paths and overlap; normal single-owner IR installation reports unsupported. |
| V `AAATGT`, D `CCC`, J `TTTAAA`; VD blunt, resect V-top 1, synapse at 6, fill top; DJ blunt at 3 | Query `AAATGTCCCTTTAAA`; V-anchor last base is synthetic, certified in place. With trusted anchor/phase evidence, the duplex policy can accept it. |
| A P tract happens to contain a canonical anchor codon | No landmark certification by motif; no invented productive result. |
| Verified noncanonical landmark without coding-phase evidence | Preserve its locatable junction evidence and determinate residue failure; missing phase cannot turn that failure into acceptance. |
| Fully absent D body or unresolved tract split | Generating D call retained; no fabricated coordinates or NP1/NP2 split. |
| Structural deletion versus read clipping of the same founder positions | Same founder ledger; different current loss reasons and CIGAR/coverage consequences. |
| Complete downstream D deletion and terminal reference-gap loss | Preserve the certified deletion path/boundaries where still known, with no invented nonempty D query interval; distinguish uncovered terminal gaps. |
| SHM, indels, clipping and reverse output after repair | Recompute observed fields while preserving founder repair history and mutation-stage denominators. |
| Clone branches with different changes and exact replay | Shared founder ancestry, independent observed maps; identical reconstructed metadata on replay. |

Validate physical transitions and annotation with independent literal oracles.
A projector and validator calling the same placement helper do not establish
correctness. Extend the existing short-substrate exhaustive sweep with
placement conservation, inverted/asymmetric intervals, removed-template
chains and both join orders. Do not assume the orders commute for all
histories.

Export tests must reconstruct sequence/reference consumption from CIGAR and
compare all coordinate spaces, test null versus empty through JSON, DataFrame
and strict TSV, and verify reserved extensions cannot be overwritten by custom
metadata. Weighted small-support tests must enumerate expected row/context
probabilities, including duplicate products, missing contexts and partial
layouts. Replay corruption tests must fail before installation.

## Engine integration points

| Area | Required change |
| --- | --- |
| `engine_rs/src/duplex_repair.rs` and `duplex_repair/coupled.rs` | Shared deterministic operation evidence and certificate observation; no parallel repair executor. |
| `src/GenAIRR/experimental.py` | Reference-bound standalone projection and exact versioned serialization. |
| `engine_rs/src/ir/simulation.rs`, `builder.rs` and `sim_event.rs` | Persistent ancestry, atomic seed installation and typed loss/update events. |
| `engine_rs/src/airr_record/` | Certified path consumption, nullable values, reserved extension and independent consistency checks. |
| `engine_rs/src/productivity.rs`, `compiled/execute.rs` and `trace_file.rs` | Active policy dispatch through feasibility, final output and replay validation. |
| Native mutation/indel paths, clonal continuation and record validators | Preserve current maps; distinguish placement role, chemical ancestry and observed loss. |
| `src/GenAIRR/_pipeline_ir.py`, `_lowering.py` and `_experiment/recombination.py` | Immutable finite model, configuration conflicts, one assembly owner and full-support capability checks. |
| `src/GenAIRR/_result_export.py` and native bridge | Preserve missing versus empty, strict coordinate conversion and extension serialization. |

Hash-cons or memoize immutable reference/context evidence where useful, with
keys covering all identities and policies. Keep detailed ledger serialization
out of per-event hot paths. Lightweight final-record references are acceptable
only if the corresponding immutable repair snapshot is retained in an
explicit export bundle; an exported truth record must not contain dangling
ancestry IDs after its in-memory outcome is released.

## Implementation order and release gates

| Slice | Deliverable | Gate before proceeding |
| --- | --- | --- |
| A | Native event evidence, reference bindings, placement certificates and standalone partial AIRR projection | Hand-derived restoration, transfer, overlap and inversion cases pass; rich histories remain conserved; null fields survive every export. |
| B | Atomic representable seed installation and persistent edit/observation maps | Mutation, indels, clipping, reversal and persistence preserve both ordinary metadata and ancestry; existing modes keep their seeded results. |
| C | Opt-in finite-history experiment lowering, duplex productivity, exact feasibility and replay | Declared local probability law verified, full-support capability preflight, final productive assertions, reference/admission gates and tamper checks pass. |
| D | Clonal continuation and export | Branch-specific coverage, event counters, summaries and replay pass independent end-to-end ancestry checks. |

Keep slices independently reviewable. Do not announce general experiment
support after only slice A or an isolated native seed test. Before completing
the milestone, run the affected Python suites, full native debug/release
suites, existing performance budgets, legacy-output checks and strict docs.
Benchmark new-mode ledger/certificate memory and compilation work separately;
immutable sharing and bounded memoization must not change probabilities or
hide exhausted work.

The completion criterion is an ordinary opt-in GenAIRR experiment that emits,
exports and replays defensible AIRR metadata after downstream edits, while
preserving the full repair ancestry and declaring every unresolved field.
Reference curation with pinned source evidence can proceed in parallel.
Learning a distribution over histories from measured data follows this
annotation contract; it does not substitute for it.
