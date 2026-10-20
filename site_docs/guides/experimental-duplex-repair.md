# Experimental duplex repair

`GenAIRR.experimental.DuplexRepair` is an explicitly authored local experiment
for two complete coding-end hairpins. The Rust kernel represents both strands
and preserves their separate origins. `CoupledDuplexRepair` extends this to two
successive joins sharing one D molecule. A finite reference-bound authored
model can feed the normal `Experiment` pipeline as described below.

```python
from GenAIRR.experimental import DuplexRepair

script = (
    DuplexRepair("AAG", "CAA")
    .open("left", offset=1)
    .open("right", offset=1)
    .synapse(right_start=3)
    .ligate()
)
product = script.run()
assert product["sequence"] == "AAGCAA"
replayed = DuplexRepair.from_dict(script.to_dict()).run()
assert replayed == product
```

Top strands follow the coding coordinate axis from 5′ to 3′. Bottom strands
follow the same axis from 3′ to 5′. `right_start` positions the original right
coding origin on the left coding axis. Opening at offset `h` transfers original
bases across the hairpin turn and produces a `2*h` 3′ overhang; opening itself
does not synthesize nucleotides. The paired annealed columns can retain origins
from both input ends.

## Supported transitions

- `open(side, offset=...)` opens each original hairpin once while retaining a
  paired stem.
- `resect(side, strand, count=...)` removes bases from the exposed end of the
  selected strand and preserves the removed origins.
- `add_n(side, sequence)` extends the exposed 3′ terminus. Its sequence follows
  polymerization order, 5′ to 3′, for either side.
- `synapse(right_start=...)` aligns the ends. Opposing bases must complement;
  overlapping same-direction strands requiring flap processing are rejected.
- `fill(strand)` extends a paired 3′ primer across a gap using a continuous
  opposing template. Each new base records its template identifier.
- `ligate()` requires two resolved 3′-OH/5′-phosphate nicks and a fully paired
  product.

Opening, resection and N extension precede synapsis in this bounded model.
Ligation completes the history. Unsupported transitions, ambiguous input,
unresolved mismatches, gaps and work-budget violations raise explicit errors.
An isolated 3′ overhang cannot be filled from a recessed 5′ end. Discontinuous
or unpaired-primer polymerase mechanisms require additional explicit operations.

Physical resection and final net loss are distinct. For example, resecting the
last top base of a blunt left end and filling it from the retained bottom strand
can restore the same final coding sequence while recording the resection and
new synthesis separately.

The result contains both original references, both final strand sequences,
paired columns, nucleotide identities, original strand coordinates, template
links, removed bases and the operation history. Only completed duplexes produce
a sequence. `max_operations` and `max_bases` bound the experiment.

## Model scope

The operations are specified mechanisms with authored choices. They carry no
fitted enzyme rates, universal opening probabilities or inferred unique repair
history. `DuplexRepair` handles a single local join, compatible annealing
and template completion; the coupled API below shares material across two joins.
Flap processing, mismatch resolution, general repair
kinetics and experimental calibration remain explicit extensions.

Verification covers nucleotide conservation during opening, paired ancestry,
opposing N additions sharing a product column, resection/restoration, strand
symmetry, exact script replay and rejection of invalid primer/template geometry.

## Coupled V-D and D-J repair

`CoupledDuplexRepair` initializes the two strands of V, D and J once. Each
`DuplexJoin` acts on the actual connected components left by preceding joins.
It preserves the original nucleotide identities and any earlier synthesis.
The following inverted-D example performs D-J first:

```python
from GenAIRR.experimental import CoupledDuplexRepair, DuplexJoin

def blunt(name, right_start):
    return (
        DuplexJoin(name)
        .open("left", offset=0)
        .open("right", offset=0)
        .synapse(right_start=right_start)
        .ligate()
    )

script = (
    CoupledDuplexRepair("ACG", "AAT", "TCC", d_inverted=True)
    .join(blunt("dj", 3))
    .join(blunt("vd", 3))
)
product = script.run()
assert product["sequence"] == "ACGATTTCC"
assert CoupledDuplexRepair.from_dict(script.to_dict()).run() == product
```

Exactly one `vd` and one `dj` join must complete. Either order is supported.
The order is an authored part of the history, without a fitted probability.
Each join supports the same opening, resection, N extension, synapsis, fill and
ligation operations as a local script. Joining consumes existing components;
it never constructs a new germline substrate from an intermediate sequence.

### Ends, coordinates and inversion

For `vd`, `left` means the exposed V tip and `right` means the exposed left tip
of oriented D. For `dj`, `left` means the exposed right tip of oriented D and
`right` means the exposed J tip. D's canonical ends exchange assembled roles
when inverted. `top`/`bottom` always describe the current physical coding axis;
an original base separately retains its canonical source strand.

For `vd`, `right_start` locates D's original oriented origin relative to V's
origin. For `dj`, it locates J's original origin relative to oriented D's origin.
These origin anchors survive resection and move with the whole component at
synapsis. They are not recomputed from the first surviving base. The result's
`origin_positions` gives their final common-axis positions, and
`ligation_seams` records both strand junction boundaries.

Canonical V/D/J references are supplied as input. Inverting D exchanges its
physical strand maps and reverses their placement. It preserves each original
base's identity, canonical strand, canonical reference coordinate and character.
It does not create an independent reverse-complement copy of D.

### Shared substrate and ancestry

Both D ends share finite physical material. The first join can alter the
substrate available to the second, including supplying paired material or
template-dependent synthesis. Opening transfers the bases actually present;
their identity and original/template ancestry survive movement between strands.

A sum of the two opening offsets is insufficient to decide support. For
example, V `AAG`, D `CA` and J `TAA`, with offset one at all four facing tips,
can produce `AAGCATAA` successively. The first partner supplies pairing retained
through the second opening; after both openings no original-original D pair
remains.
Transition checks use the current component geometry and the protected prior
ligation bonds.

Each `columns` entry contains both physical strand contributors. Its bases use
canonical `original` records with `segment`, `strand` and `reference_position`.
The global `ledger` contains every original or synthesized nucleotide once.
Every allocated ID is present either in final columns or in `removed`, and
every template ID resolves in the ledger even if that template was later
removed. New synthesis is never relabeled as germline by a subsequent join.

`max_operations` is shared across both scripts. `max_bases` counts the two
original strands of all three references and every newly synthesized base;
removal refunds none. Partial, incompatible, over-budget or unresolved products
raise errors and export no completed sequence.

### Standalone AIRR projection

Completed coupled scripts can produce a reference-bound AIRR record and
detailed ancestry with `script.project(...)`. Supply the actual generating
allele names and a `RefDataConfig`; all three substrates must match those
whole canonical alleles. Names are never inferred from sequence similarity.
This standalone adapter remains separate from `Experiment`.

The following uses deliberately short synthetic references as a geometry
control. Its explicit coding-end assertion is appropriate to these authored
substrates; it does not certify biological reference annotations.

```python
from GenAIRR import RefDataConfig
from GenAIRR.experimental import export_duplex_projection

reference = RefDataConfig.vdj()
reference.add_v_allele("V*01", "V", b"ACG")
reference.add_d_allele("D*01", "D", b"AAT")
reference.add_j_allele("J*01", "J", b"TCC")

projection = script.project(
    reference, v_call="V*01", d_call="D*01", j_call="J*01",
    reference_ends_are_coding_ends=True,
)
assert projection["record"]["sequence"] == "ACGATTTCC"
assert projection["record"]["d_germline_start"] == 0
export_duplex_projection(projection, "duplex.json")
export_duplex_projection(projection, "duplex.tsv", format="tsv")
```

Unknown junction-facing boundary evidence requires the explicit assertion.
A known partial boundary cannot be overridden. Landmark and phase readiness
remain separate from boundary admission.

`record` contains the ordinary fields and a reserved `duplex_repair` extension
with both strands, operation events, placement certificates and the original
repair ledger. `paths` contains independent per-gene alignments. The query is
the coding-oriented top strand; a body position supported on either paired
strand can justify its alignment. Synthesized restoration stays chemically
synthetic, while a transferred hairpin base loses its former body placement.

A valid V-body/D-body overlap retains both gene paths. The combined alignment
and its relative coordinates remain missing, with reasons in `availability`;
the adapter never chooses one owner to hide the overlap. Empty NP intervals
have length zero. Unknown values remain `None`. Existing P scalars remain
missing because this kernel has no equivalent legacy P-emission contract;
physical events and query contributions after both joins have separate counts.
The legacy P, NP, trim and retained-outcome fitters reject duplex provenance;
they cannot reinterpret these new histories as their old scalar models.

`capabilities.simulation_seed` reports whether the bound product can enter
the native single-segment execution IR, including its coordinate limits.
`experiment_integration` is false for standalone products. Native projection verifies the snapshot and its evidence
by exact deterministic replay before returning a record.

Productivity uses the explicit `duplex_observed_vdj_v1` policy. It accepts
certified in-place landmarks restored by synthesis and evaluates their actual
codons. A matching Cys or Phe/Trp elsewhere cannot substitute for a lost
landmark. A verified rejected residue can establish `False` even when coding
phase is unknown. A defined junction keeps its translation on nonproductive
products. This predicts sequence-level acceptance in the observed coding span,
not expression or receptor function.

JSON exports retain the complete envelope in native coordinates. Strict TSV
exports one AIRR record, converts its interval coordinates, and serializes the
ancestry extension as one JSON cell. Nested repair/provenance coordinates
retain their declared systems. Inverted D retains canonical forward reference
bounds; its CIGAR uses the oriented reference offset.

### Persistent native seeds and observations

A representable completed script can now create a persistent native simulation.
The bound references, immutable repair snapshot and certified body placements
travel with that simulation through audited native continuations:

```python
from GenAIRR import _engine
from GenAIRR.experimental import project_duplex_simulation

seed = script.seed(
    reference, v_call="V*01", d_call="D*01", j_call="J*01",
    reference_ends_are_coding_ends=True,
)
plan = _engine.PassPlan()
plan.push_mutate_uniform([(1, 1.0)])
plan.push_corrupt_3prime_loss([(1, 1.0)])
compiled = plan.compile(
    refdata=reference, strict=True, allow_curatable_refdata=True,
)
outcome = compiled.run_from(seed, 123)
observation = project_duplex_simulation(outcome.final_simulation(), reference)
export_duplex_projection(observation, "duplex-observation.json")
```

This is an experimental native continuation API. The finite-history
`Experiment` integration below provides ordinary generation, family
continuation and replay. It adds no external tools or dependencies.

Supported continuations are the existing uniform and S5F mutation, PCR and
quality substitution, N masking, structural indel, terminal loss and read
reverse-complement passes. Every pass in a continuation is checked before any
pass runs, including unsupported suffixes. Source resampling, trimming,
new P/NP generation, contamination, echo, paired-end windowing and custom
passes are rejected explicitly. An empty continuation is valid. References
must match the founder binding; generating alleles and net repair trims cannot
be replaced. Legacy preservation-only `AnchorPreserved` is rejected; productive
contracts use the named duplex policy and current certified landmarks.

`record.duplex_repair.product`, its physical ledger, operation evidence and
repair counts always describe the founder after the two joins. They are never
rewritten by downstream mutation or read loss. The `observation` extension
records current query links, cumulative edits, removal causes and sampled
PCR/quality stages. A founder-column link describes ancestry after copying;
it does not assert that descendants contain the founder's physical atoms.
Inserted observation identities are local to an observation branch. Independent
branches require their observation identity when comparing insertion links.

Current alignments use generating alleles and certified body coordinates.
Re-calls, when computed by the native continuation, appear only in separate
`recalled_v_call`, `recalled_d_call` and `recalled_j_call` fields. A structural
deletion remains a reference-consuming `D`; clipping reduces observed reference
coverage. Complete structural deletion can retain a deletion-only gene path
without inventing query interval endpoints. The path's `query_boundary` locates
that zero-width placement for strict CIGAR export. NP amino acids include only
complete codons inside the tract in the current molecule's coding frame.
Unknown phase gives a missing amino-acid value.

Mutation burden retains each stage's input molecule denominator even after
read loss. Mechanism counts accumulate over the observation's history. In
particular, `n_pcr_errors` and `n_quality_errors` count sampled attempts;
`duplex_repair.applied_error_substitutions` separately counts applied
substitutions. Attempt summaries include zero-effect and empty-query stages.
These cumulative counts remain available across successive continuations,
whereas a single legacy trace is not a complete descendant history. Low-level
unattributed edits make mechanism counters unavailable instead of assigning
a biological cause.

Reverse-complement output changes `record.sequence` and `rev_comp`; alignments,
landmarks, query links and their coordinates remain coding oriented. The
projector verifies the immutable repair snapshot, full edit and query replay,
source assignments, coverage and structural partition before emitting metadata.
Unknown legacy P lengths remain missing. Use this nullable projector and its
JSON/TSV exporter, or normal nullable `Outcome`/`SimulationResult` record
exports. Durable traces from externally supplied seeds still require an
initial-snapshot format. Source-backed experiment traces below bind
the generating model and support ordinary replay. The stored repair script
can also be replayed independently.

### Finite authored histories in Experiment

`DuplexRepairModel` pins a finite set of named histories to exact reference
content and evidence. `Experiment.recombine(duplex_repair=model)` checks all
naturally selectable positive histories before any draw, including source
coverage, physical validity and compatibility with the native execution layout.
It raises an error for an unsupported active row and never retries another
history to hide that limitation. Explicit allele restrictions narrow the
requested source contexts; genotype eligibility remains chromosome linked.

```python
from GenAIRR import Experiment, RefDataConfig
from GenAIRR.experimental import DuplexHistory, DuplexRepairModel

refs = RefDataConfig.vdj()
refs.add_v_allele("V*01", "V", b"AAATGT", anchor=3)
refs.add_d_allele("D*01", "D", b"CCC")
refs.add_j_allele("J*01", "J", b"TTTAAA", anchor=0)

blunt_vd = blunt("vd", 6)
restored_vd = (
    DuplexJoin("vd").open("left", offset=0).open("right", offset=0)
    .resect("left", "top", count=1).synapse(right_start=6)
    .fill("top").ligate()
)
first = CoupledDuplexRepair("AAATGT", "CCC", "TTTAAA").join(blunt_vd).join(blunt("dj", 3))
second = CoupledDuplexRepair("AAATGT", "CCC", "TTTAAA").join(restored_vd).join(blunt("dj", 3))
model = DuplexRepairModel.bind(
    refs,
    [DuplexHistory("blunt", first, "V*01", "D*01", "J*01", weight=1),
     DuplexHistory("restored", second, "V*01", "D*01", "J*01", weight=3)],
    reference_ends_are_coding_ends=True,
)
compiled = (
    Experiment.on(refs).recombine(duplex_repair=model)
    .mutate(model="uniform", count=2).end_loss_3prime(length=1)
    .productive_only().compile()
)
result = compiled.run_records(n=5, seed=42, strict=True, validate_records=True)
result.to_tsv("duplex-experiment.tsv")
outcome = compiled.run(seed=42, strict=True)[0]
trace = compiled.simulator.trace_file_from(outcome, 42)
replayed = compiled.simulator.replay_from_trace_file(trace, strict=True)
assert outcome.airr_record(refs, sequence_id="same") == replayed.airr_record(refs, sequence_id="same")
```

These short references are geometry/annotation controls, so their supplied
coding boundaries are explicitly asserted. The weights describe a 1:3 choice
between the two authored histories; they are not measured repair frequencies.
Their final DNA is identical, but the histories and physical ancestry remain
distinct. Models serialize with `to_dict`/`to_json` and reconstruct with
`DuplexRepairModel.from_dict`; row order and unique IDs are preserved.

The named law is `finite_authored_history_local_conditioning_v1`. Histories are
weighted within the selected V/D/J/orientation context. Under
`productive_only`, only definitely productive rows remain eligible. Earlier
allele/gene/chromosome/orientation choices retain their natural relative weights
among feasible continuations; repair success fractions never multiply those
weights. This is sequential local conditioning, rather than global conditioning
of the entire joint law. Downstream edits use the same duplex productivity
assessment, and every successful productive emission remains `True`.

In a phased genotype, each eligible gene keeps its natural usage times total
positive carried-copy dosage. The later copy draw is conditioned within that
gene. A failed copy therefore cannot reduce an otherwise eligible gene's
upstream weight.

Use `.invert_d(prob=...)` to select D orientation; provide positive history
coverage for every orientation with nonzero natural probability. Zero-weight
rows are disabled and cannot supply required context coverage. An all-zero
active context is invalid. Missing context rows raise a compilation error.

The authored model owns both junctions. Explicit legacy NP-length arguments,
`.trim(...)` (including `enabled=False`), and original-end processing conflict
with this mode. Default cartridge trim/P/NP assembly is bypassed; allele-usage
weights remain upstream. VJ, receptor revision, contamination,
paired-end read windows and arbitrary custom contracts require separate
adapters. They are rejected before generation, including zero-probability
configurations of unsupported passes.

Compilation has explicit `max_preflight_operations` and
`max_preflight_bases` limits. Operations include reachable-context checks,
authored transitions and inspected source weights, including distinct
productivity-filtered flat supports and natural genotype masses; bases reserve
the active scripts' declared allocation
budgets. Exhaustion is a resource error, not empty biological support. Binary64
weight normalization rejects numerical loss of positive support rather than
silently dropping a row. These controls carry no enzyme kinetics.

Ordinary result and stream exports preserve nullable values and full ancestry.
`duplex_repair.selection` describes the selected authored row;
`duplex_repair.current_paths` contains current native alignment paths. Strict
TSV converts only AIRR fields and keeps this extension as one JSON cell.
DataFrames retain nullable integer and boolean columns. `validate_records`
checks actual records and typed nested provenance. Generated fields, including
`duplex_repair`, cannot be overwritten by custom metadata.

Source-backed durable traces bind references, model content, the complete
upstream sampling law and active policy/contracts. Replay checks positive
support in the recorded context and executes the same kernel without a new
history draw. External-seed `run_from` traces still require an initial-snapshot
format; a durable full-family package and measured calibration remain later steps.

### Families with unfitted history weights

Omitting a history's weight uses `1.0`, an equal-weight control over the valid
authored history IDs in that source context. This supplies a useful starting
distribution without measured repair data. It does not assign enzyme rates or
make distinct histories with the same DNA one choice. Later fitted weights can
replace these values through the same model API.
The example labels their origin with sample metadata; equal numeric weights
alone do not establish how a model was obtained.

Repair is sampled once per clone. Flat expansion and lineage branches inherit
the immutable founder ledger while retaining independent observed histories.
Lineage survival retries reuse that founder and retry only growth. Existing
neutral growth and mutation parameters remain explicit model controls.

```python
equal_model = DuplexRepairModel.bind(
    refs,
    [DuplexHistory("blunt", first, "V*01", "D*01", "J*01"),
     DuplexHistory("restored", second, "V*01", "D*01", "J*01")],
    reference_ends_are_coding_ends=True,
)
family_result = (
    Experiment.on(refs).with_metadata(repair_weight_source="unfitted_uniform_authored")
    .recombine(duplex_repair=equal_model)
    .clonal_lineage(n_clones=1, max_generations=2, n_max=8, n_sample=6,
                    rate=0.05, lambda_base=3)
    .productive_only().run_records(seed=42, strict=True, validate_records=True)
)
assert len(family_result) > 0
assert all(r["repair_weight_source"] == "unfitted_uniform_authored" for r in family_result)
assert len({r["duplex_repair"]["selection"]["history_id"] for r in family_result}) == 1
family = family_result.families[0]
observed = next(node for node in family_result.lineage_trees[0].nodes() if node.observed)
replayed_node = family.replay_node(observed.id, refs)
original_node = family.node_outcomes()[observed.id]
assert replayed_node.airr_record(refs, sequence_id="same") == original_node.airr_record(refs, sequence_id="same")
family_result.to_tsv("duplex-family.tsv")
```

Duplex lineage records retain distinct sampled cells even when their DNA is
identical. Repeated draws of the same cell contribute to its
`lineage_abundance` and `duplicate_count`. Duplex repertoire collapse requires
all generated truth to agree, so equal DNA cannot hide different mutation,
indel or ancestry histories. Existing legacy collapse behavior is unchanged.

`parents` and `parent_id` connect records to their actual founder. Lineage
results also retain native `families` aligned with `lineage_trees`.
`lineage_mutations_from_parent` counts actual events on the incoming edge;
cumulative mutation counts include reversions. Family validation compares
immutable repair evidence while allowing current maps and paths to diverge.
Lineage record validation checks node metadata against its retained origin.
Local insertion identities require their enclosing clone/node/read context.

`replay_node` accepts any retained node, including roots and unsampled cells.
It consumes the recorded division choices without replacement draws
and returns the cellular state before any per-read observation passes. It is
an in-memory API requiring the matching bound reference; it does not claim a
durable whole-family format. Ordinary compiled trace writers continue to reject
merged family outcomes. Affinity scoring, when requested, translates the
contiguous certified coding span across junction boundaries; unresolved coding
phase cannot be replaced with an invented affinity.

### Bounded substrate contract

The input supplies three coding duplexes and four junction-facing hairpins.
The model does not simulate RSS recognition, RAG cleavage activation or a
detailed cleavage chronology. Both join orders are mechanistic controls under
this supplied-substrate contract, rather than a repertoire-wide claim about
their prevalence.

Completed ligation bonds and their two adjacent physical nucleotides are
protected. Processing that would consume or cross those bonds is outside this
bounded model and raises an explicit error. It is not classified as biologically
impossible in every repair mechanism. General flap/mismatch repair, reopening
completed joins, stochastic enzyme kinetics and biochemical calibration remain
extensions. Fully paired completion is required for the final VDJ product.

Hairpin opening and end processing are supported by primary biochemical work,
including sequence-dependent opening in
[Lu, Schwarz and Lieber (2007)](https://pubmed.ncbi.nlm.nih.gov/17932067/)
and coding-joint reconstitution in
[Lu and colleagues (2008)](https://pubmed.ncbi.nlm.nih.gov/18722175/).
The coupled short-substrate controls above are GenAIRR geometry tests, not
measurements from those studies. No universal enzyme rates or empirical
frequencies are assigned to the authored scripts.
