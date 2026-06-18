# The AIRR record

<p class="lead">An AIRR record is one row of GenAIRR's output -
one simulated molecule (or one paired-end read pair) with its
ground truth, its observed sequence, and the counters that
describe what happened on the way through the pipeline. This
page is the field-by-field tour, organised by category, with the
biology-vs-artefact distinctions made explicit so you can pick
the right field for the question you're asking.</p>

## Derived, not accumulated

A GenAIRR record is **not** an accumulator of fields filled in
across the pipeline. It's computed once, at the end, by reading
the final intermediate-representation pool. Coordinates, CIGARs,
junction, productivity, calls - everything is derived from one
source of truth, so nothing can drift out of sync.

Most simulators carry an AIRR record alongside the sequence and
update fields as the pipeline progresses. An indel happens, and
they scramble to decrement coordinates downstream; a trim
arrives, and they bookkeep the gene end. Every update is a
chance to introduce a desynchronisation bug.

GenAIRR keeps no such running record. Coordinates exist only on
the final pool. CIGARs are walked once at the end. `productive`
is decided by reading the junction codons from the final state.
There is no parallel state to fall out of sync - there is
nothing to fall out of sync *with*.

That property is also why the [`validate_records`](../validation/validate-records.md)
gate works: the validator re-runs the same projection and asserts
agreement with the record. If projection is pure (and it is),
the validator can independently re-derive every field and catch
any drift introduced by a code-path change.

## What an AIRR record is

A record is a `dict` (or one row of `result.to_dataframe()`).
Field names follow the AIRR-C v1 schema where one exists, with
GenAIRR-specific additions for ground truth, mutation
provenance, corruption counters, clonal structure, and
paired-end layout. The same column order ships in
`SimulationResult.to_tsv()` / `.to_csv()` and is reachable
programmatically as `result.to_dataframe().columns`.

The canonical groups, in the order they appear on the record:

1. AIRR metadata (sequence + length + locus)
2. V / D / J calls + alignment coordinates + trim
3. Junction
4. NP regions
5. Productivity flags
6. SHM and corruption counters
7. Per-segment indel + mutation partitions
8. V-subregion mutation partition
9. End-loss artefacts
10. Contamination flag
11. D-inversion provenance
12. Receptor-revision provenance
13. Read layout (paired-end)
14. Clonal stamping (when a clonal workflow runs)
15. Truth columns (when `expose_provenance=True`)

The remaining sections walk through each.

## Core sequence fields

| Field | Type | Meaning |
|---|---|---|
| `sequence_id` | str | `"{prefix}{i}"` - defaults to `"seq0"`, `"seq1"`, … |
| `sequence` | str | The observed nucleotide sequence; lower-case bases mark corrupted positions (low quality / PCR / sequencing error) |
| `sequence_aa` | str | Amino-acid translation in the V reading frame |
| `sequence_alignment` | str | Gapped sequence aligned to the germline (when present) |
| `germline_alignment` | str | The germline alignment of `sequence_alignment` |
| `germline_alignment_d_mask` | str | Germline alignment with the D region masked |
| `sequence_length` | int | Length of `sequence` in nucleotides |
| `rev_comp` | bool | `True` when the projection reverse-complemented the molecule (set by `random_strand_orientation`) |
| `locus` | str | Cartridge identity (e.g. `"IGH"`, `"TRB"`) |

The lower-case base convention on `sequence` is load-bearing -
the FASTQ exporter (`to_fastq` / `to_paired_fastq`) routes those
positions to the low-quality bucket. See the
[Corruption + sequencing artefacts guide](../guides/corruption-sequencing.md)
for what each corruption pass writes into `sequence`.

## Allele calls and truth

```text
v_call, d_call, j_call    ← the call surface
truth_v_call, truth_d_call, truth_j_call  ← provenance (opt-in)
```

`v_call` (and `d_call`, `j_call`) carry the **resolved allele
name** committed during recombination - or the post-revision
identity after `receptor_revision()`. The call field is a
comma-separated **tie set** when the engine sampled an
indistinguishable group:

```text
v_call: "IGHV3-23*01"
v_call: "IGHV3-23*01,IGHV3-23*04"   ← tie set: two-way ambiguity
```

Downstream code that wants the single committed allele should
split on `","` and pick the first entry; downstream code that
wants to credit ambiguity fractionally should split and divide.
Either policy is valid; just pick one.

**`d_call`** is empty (`""`) on VJ-chain records (kappa, lambda,
TCR alpha). The empty string - not `None` - is the canonical
absent value, mirroring how missing fields work in AIRR TSV.

### Truth columns (opt-in)

When the experiment runs with `expose_provenance=True`:

```python
result = exp.run_records(n=10, seed=42, expose_provenance=True)
```

the records also carry `truth_v_call`, `truth_d_call`, and
`truth_j_call`. These are the **same** call values as
`v_call` / `d_call` / `j_call` in the ordinary case (GenAIRR's
records *are* truth by construction). The truth columns become
load-bearing in two cases:

- **Family validation** uses them to confirm every descendant of
  a clonal family shares the recombination-time identity (see
  [Clonal simulation overview](../guides/clonal-families.md)).
- **Aligner benchmarking** treats `v_call` as the aligner's
  prediction and `truth_*_call` as the ground truth, even when
  both came from the same simulator - the column split makes the
  benchmark script's join symmetric.

Without `expose_provenance=True`, the truth columns are absent
entirely.

## Junction and productivity

| Field | Type | Meaning |
|---|---|---|
| `junction` | str | Junction sequence (V end + NP1 + D + NP2 + J start, conserved Cys to Trp/Phe) |
| `junction_aa` | str | Amino-acid translation of `junction` |
| `junction_start`, `junction_end` | int | Junction coordinates within `sequence` |
| `junction_length` | int | Length of `junction` in nucleotides (always a multiple of 3 when `productive: True`) |
| `productive` | bool | `True` iff in-frame junction + no junction stop codon + V Cys preserved + J Trp / Phe preserved |
| `vj_in_frame` | bool | V and J in the same reading frame |
| `stop_codon` | bool | `True` iff a stop codon exists in `sequence_aa` |

When the pipeline includes `.productive_only()`, every record
has `productive: True` by construction - the constraint masks
the sampling support before the draw, so the engine never
produces an unproductive record in the first place. See
[Recombination and junction biology](../guides/recombination-junction.md#productivity)
for the four-clause definition.

## Recombination structure

### Per-segment alignment coordinates

For each of V, D, J:

| Field | Meaning |
|---|---|
| `{v,d,j}_sequence_start` | Position in `sequence` where this segment begins (0-based by default; 1-based when `airr_strict=True`) |
| `{v,d,j}_sequence_end` | Exclusive end position in `sequence` |
| `{v,d,j}_alignment_start` | Position in `sequence_alignment` |
| `{v,d,j}_alignment_end` | Exclusive alignment end |
| `{v,d,j}_germline_start` | Position in the germline allele |
| `{v,d,j}_germline_end` | Exclusive germline end |

### Trim fields (recombination-stage)

| Field | Meaning |
|---|---|
| `v_trim_3` | Bases removed from the V allele's 3′ end during recombination |
| `d_trim_5`, `d_trim_3` | Bases removed from the D allele's 5′ and 3′ ends |
| `j_trim_5` | Bases removed from the J allele's 5′ end |
| `v_trim_5`, `j_trim_3` | Always `0` - these positions aren't trimmed during recombination (the canonical biology) |

These four trim fields are the recombination-stage diet. They
are **not** the observation-stage length loss - that's
`end_loss_5_length` and `end_loss_3_length` (next section).
Mixing them up is the most common confusion on the record
surface; the trim fields describe biology (the recombinase
chewed bases off allele ends before junction joining), while
end-loss describes the sequencer (the read ran short).

### NP regions

| Field | Meaning |
|---|---|
| `np1` | Non-templated bases between V and D ends - **P-clean** (V–D junction in VDJ; V–J junction in VJ) |
| `np2` | Non-templated bases between D and J ends - **P-clean** (VDJ only; empty on VJ) |
| `np1_aa`, `np2_aa` | Amino-acid translations |
| `np1_length`, `np2_length` | Lengths in nucleotides |

`np1` and `np2` are the **non-templated** strings only. When the
engine claims a P-nucleotide span back as a templated extension
of V, D, or J, those positions drop out of `np1` / `np2` - the
NP strings are P-clean by construction.

### P-nucleotide lengths

Four per-record fields carry the palindromic-insertion lengths
sampled during recombination:

| Field | Meaning |
|---|---|
| `p_v_3_length` | Number of P bases off the V allele's 3′ end (V → NP1 side) |
| `p_d_5_length` | Number of P bases off the D allele's 5′ end (NP1 → D side) - VDJ only, `0` on VJ |
| `p_d_3_length` | Number of P bases off the D allele's 3′ end (D → NP2 side) - VDJ only, `0` on VJ |
| `p_j_5_length` | Number of P bases off the J allele's 5′ end (NP2 → J side on VDJ; NP1 → J side on VJ) |

P bases **contribute to `sequence` and `junction`** (they are
real palindromic nucleotides in the assembled molecule) but
`np1` / `np2` remain N-only. **GenAIRR exposes P lengths, not P
strings** - there is no per-base P field. If you need the actual
P-nucleotide bases for a record, slice them from `sequence` using
the per-segment coordinates plus the four length fields.

The engine samples P-insertion lengths from
`cfg.reference_models.p_nucleotide_lengths`. See
[Junction N/P additions](../guides/junction-additions.md) for the
layout diagrams.

### CIGAR fields

| Field | Meaning |
|---|---|
| `{v,d,j}_cigar` | Per-segment CIGAR string against the assigned allele |
| `{v,d,j}_score` | Per-segment alignment score (when reported by the projector) |
| `{v,d,j}_identity` | Per-segment percent identity |
| `{v,d,j}_support` | Per-segment support value (when reported) |
| `c_call` | Constant-region call (empty when no C-region biology is modeled) |

## Mutation and artefact counters

### Biological SHM

```text
n_mutations            ← total biological SHM events (canonical)
n_v_mutations          ← V-segment SHM
n_d_mutations          ← D-segment SHM
n_j_mutations          ← J-segment SHM
n_np_mutations         ← NP1 + NP2 combined SHM
mutation_rate          ← realised per-base rate for this record
```

These five sum cleanly:
`n_v + n_d + n_j + n_np == n_mutations` on every record. The
[`validate_records`](../validation/validate-records.md) gate
checks this equality and fires `MutationCountSumMismatch` if
it ever breaks.

**`n_mutations` is biology only.** PCR errors, sequencing
errors, indel-pass errors, and end-loss never increment these
counters - they have their own.

### V-subregion mutation partition

Six fields that partition `n_v_mutations` by the assigned V
allele's IMGT subregion intervals:

```text
n_fwr1_mutations
n_cdr1_mutations
n_fwr2_mutations
n_cdr2_mutations
n_fwr3_mutations
n_v_unannotated_mutations
```

Sum equals `n_v_mutations` by construction. The unannotated
bucket catches three legitimate non-zero cases (the
V-side CDR3 stretch is the most common); see
[SHM and mutation targeting](../guides/shm-targeting.md#reading-mutation-counters)
for the full rules.

### Library + sequencing artefacts (non-biological)

```text
n_pcr_errors           ← bases mutated by PcrAmplifyPass
n_quality_errors       ← bases corrupted by SequencingErrorsPass / AmbiguousBaseCallsPass
n_indels               ← total indels across all PolymeraseIndelsPass passes
n_v_indels             ← indels landing in the V segment
n_d_indels             ← indels landing in the D segment
n_j_indels             ← indels landing in the J segment
end_loss_5_length      ← bases lost from the 5′ end (EndLossPass / primer_trim_5prime)
end_loss_3_length      ← bases lost from the 3′ end (EndLossPass / primer_trim_3prime)
is_contaminant         ← True when this record is a contaminant (set by `contaminate`)
```

`n_v_indels + n_d_indels + n_j_indels ≤ n_indels` - indels that
land in NP1 or NP2 are counted in `n_indels` but not in any
per-segment bucket (NP indels don't belong to a germline
segment).

`primer_trim_*prime` is a backwards-compatibility alias for
`end_loss_*prime` - both write the same `end_loss_*_length`
field.

## Advanced mechanism provenance

Three fields surface the engine's recombination-stage editing
decisions:

| Field | Type | Meaning |
|---|---|---|
| `d_inverted` | bool | `True` when `invert_d()` committed the D allele in reverse-complement orientation; `False` otherwise (VJ chains, VDJ without `invert_d`, inversion that landed on the forward branch) |
| `receptor_revision_applied` | bool | `True` when `receptor_revision()` fired and replaced the committed V; `False` otherwise |
| `original_v_call` | str | When `receptor_revision_applied: True`, the V allele name the recombine pass originally committed (before revision). Empty string `""` otherwise - never `None` |

When `receptor_revision_applied: True`, `v_call` reports the
**post-revision identity** and `original_v_call` carries the
pre-revision name. See
[Recombination editing (D inversion + receptor revision)](../guides/recombination-editing.md)
for the biology and the rules.

## Read layout

Eight fields populated only when `.paired_end(...)` is in the
pipeline:

| Field | Type | Default when absent |
|---|---|---|
| `read_layout` | str | `""` |
| `r1_sequence` | str | `""` |
| `r2_sequence` | str | `""` |
| `r1_start`, `r1_end` | int | `0` |
| `r2_start`, `r2_end` | int | `0` |
| `insert_size` | int | `0` |

On single-molecule pipelines (no `.paired_end(...)` call) all
eight default to their sentinel values, but the columns are
still present in the record. Treat a non-empty `r1_sequence` as
the canonical "this is a paired-end record" check; `read_layout`
also carries the layout label (`"paired_end"` when set).

See [Paired-end reads and FASTQ](../guides/paired-end-fastq.md)
for the quality model and the layout coordinates.

## Clonal fields

All modern clonal workflows stamp `clone_id`, but the rest of the
surface depends on which clonal model produced the record:

| Field | Type | Produced by | Meaning |
|---|---|---|---|
| `clone_id` | int | `clonal_lineage`, `clonal_repertoire`, `expand_clones` | Planted clone / family label (0-based) |
| `duplicate_count` | int | `clonal_lineage`, `clonal_repertoire` | AIRR-standard abundance after genotype collapse |
| `parent_id` | int | legacy `expand_clones` | Identifier of the ancestor `Outcome`; equals `clone_id` for records expanded from that ancestor |
| `lineage_node_id` | int | `clonal_lineage` | Node id of the observed cell in the ground-truth lineage tree |
| `lineage_parent_id` | int | `clonal_lineage` | Parent node id in the lineage tree (`-1` for founder) |
| `lineage_generation` | int | `clonal_lineage` | Generation depth of the observed cell |
| `lineage_abundance` | int | `clonal_lineage` | Observation count after final-cell sampling and genotype collapse |
| `lineage_affinity` | float | `clonal_lineage` | Sequence-distance proxy to the target; 0 only when no affinity model is active |

On non-clonal runs, these fields are **absent** from the record dict. Don't write
code that assumes they're always present; check with `"clone_id" in rec`.
`validate_families()` is a safe no-op on non-clonal batches because it handles
the absent case explicitly. See
[Clonal simulation overview](../guides/clonal-families.md).

## Validation

The records page composes cleanly with the validation surface:

- **[`validate_records(refdata)`](../validation/validate-records.md)**
  re-derives every counter, every coordinate, and every truth
  field from the underlying engine `Outcome` and checks that
  the record agrees. This is the load-bearing AIRR-output gate.
- **`validate_families(refdata=None)`** groups records by
  `clone_id` and asserts each family agrees on `truth_v_call`,
  `truth_d_call`, `truth_j_call` (when present).
- **`validate_families_with_parents(refdata)`** is for legacy
  `expand_clones` results with `result.parents`; it compares each
  descendant against its actual parent `Outcome`.

The full validation picture lives at the
[Validation hub](../validation/index.md).

## Common mistakes

A handful of issues that show up repeatedly with the record
surface.

**Treating all `n_*` counters as biological mutations.** Only
the SHM partition (`n_mutations`, `n_v_mutations`,
`n_d_mutations`, `n_j_mutations`, `n_np_mutations`, the six
V-subregion fields) describes biology. `n_pcr_errors`,
`n_quality_errors`, `n_indels`, and the end-loss lengths
describe library / sequencer artefacts. They live on the same
record by design - they don't share a counter.

**Assuming call fields are single alleles.** `v_call`,
`d_call`, and `j_call` can be comma-separated tie sets. Split on
`","` before parsing; even a single-allele call can become
ambiguous in a future cartridge revision.

**Confusing trim with end-loss.** `v_trim_3` is recombination
biology (recombinase chewed bases off the V allele's 3′ end
before junction joining). `end_loss_3_length` is sequencing
artefact (the read ran short / was 3′ end-loss-clipped). The
field names are similar; the biology is different.

**Inferring P-nucleotides from `np1` / `np2`.** Don't. The NP
strings are P-clean - P bases that have been claimed back as
templated extensions of V / D / J are deliberately excluded
from `np1` / `np2`. Use the four `p_*_length` fields
(`p_v_3_length`, `p_d_5_length`, `p_d_3_length`, `p_j_5_length`)
to read palindromic-insertion lengths. GenAIRR exposes lengths,
not P strings; slice `sequence` with the per-segment coordinates
if you need the actual P bases.

**Expecting paired-end fields without `.paired_end(...)`.** The
eight read-layout fields are present in every record but default
to empty / zero on single-molecule pipelines. `r1_sequence ==
""` is the canonical "no paired-end projection ran" check. Don't
try to write FASTQ from a single-molecule pipeline - `to_fastq`
emits the assembled sequence; `to_paired_fastq` raises if
`read_layout != "paired_end"`.

**Expecting `truth_*_call` columns to always be present.** They
appear only when `expose_provenance=True` is passed to
`run_records(...)`. Without the flag, the columns are absent
entirely - not `None`-valued, absent.

**Expecting `clone_id` on non-clonal records.** Without a clonal
workflow (`clonal_lineage`, `clonal_repertoire`, or legacy
`expand_clones`), clonal fields are not stamped on the record dict at
all. Check for presence with `"clone_id" in rec`, not
`rec.get("clone_id") is not None`.

## Where to go next

- **[Your first AIRR record](../getting-started/first-airr-record.md)**
  a worked walk-through of one record, field by field.
- **[Export the results](../getting-started/export-results.md)**
  how records become TSV / CSV / FASTA / FASTQ.
- **[SHM and mutation targeting](../guides/shm-targeting.md)**
  the SHM counters in depth.
- **[Corruption and sequencing artefacts](../guides/corruption-sequencing.md)**
  the artefact counters in depth.
- **[Clonal simulation overview](../guides/clonal-families.md)** -
  `clone_id`, `duplicate_count`, lineage metadata, and family validation.
- **[Validation hub](../validation/index.md)** - re-deriving every
  field from the underlying `Outcome`.
