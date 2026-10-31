# The AIRR record

<p class="lead">An AIRR record is one row of GenAIRR's output -
one simulated molecule (or one paired-end read pair) with its
ground truth, its observed sequence, and the counters that
describe what happened on the way through the pipeline. This
page is the field-by-field tour, organised by category, with the
biology-vs-artefact distinctions made explicit so you can pick
the right field for the question you're asking.</p>

## Derived, not accumulated

GenAIRR builds the record at the end from the final nucleotide
pool, allele assignments, live calls, trace, and event history.
Coordinates and CIGARs describe the final alignment; productivity
checks the observed coding span, reading frame, and conserved
anchors. Persistent assignments retain recombination provenance
across clonal forks.

[`validate_records`](../validation/validate-records.md) checks the
record against engine state and reference data. Passing confirms
agreement with the implemented checks; it does not prove that every
biological assumption or reference annotation is correct. Projection
and validation share some helpers, so independent checks remain useful.

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

Biological annotations use coding orientation. When `rev_comp=True`,
`sequence` contains the observed reverse-strand molecule, while
`sequence_aa`, alignments, junction/NP strings, and biological query
coordinates refer to its reverse complement. Paired-read windows
instead index the observed `sequence`.

The lower-case base convention on `sequence` is load-bearing -
the FASTQ exporter (`to_fastq` / `to_paired_fastq`) routes those
positions to the low-quality bucket. See the
[Corruption + sequencing artefacts guide](../guides/corruption-sequencing.md)
for what each corruption pass writes into `sequence`.

## Allele calls and truth

`v_call`, `d_call`, and `j_call` name the final generating alleles. They are
single-valued and define the references used by the germline alignments,
CIGARs, identities, and germline coordinates. Receptor revision updates
`v_call` to the replacement allele; `original_v_call` preserves its predecessor.
`d_call` is the empty string on VJ chains.

Post-mutation inference is exported separately as `recalled_v_call`,
`recalled_d_call`, and `recalled_j_call`. These may contain comma-separated
ties or differ from the generating allele. When the generating allele is
in a tie, it appears first; otherwise calls follow allele-ID order.
Donor-genotype experiments restrict re-calls to alleles carried by that donor.
Re-calls do not change the ground-truth reference or absorb junction bases
into a segment's alignment.
When no inferred evidence is available, for example after complete clipping,
the corresponding `recalled_*_call` is empty; the standard generating call remains.

```python
rec["v_call"]           # generating allele
rec["recalled_v_call"]  # evidence-based call; split on commas for ties
```

### Truth columns (opt-in)

`run_records(..., expose_provenance=True)` also emits `truth_v_call`,
`truth_d_call`, and `truth_j_call`. These compatibility aliases equal the
standard calls. They remain available for existing family-validation and
benchmarking workflows.

## Junction and productivity

| Field | Type | Meaning |
|---|---|---|
| `junction` | str | Junction sequence (V end + NP1 + D + NP2 + J start), from the first surviving base of the conserved Cys codon to the last surviving base of the Trp/Phe codon; empty when a whole conserved codon is gone or not covered by the read |
| `junction_aa` | str | Amino-acid translation of `junction` |
| `junction_start`, `junction_end` | int | Junction coordinates within the coding-oriented sequence |
| `junction_length` | int | Length of `junction` in nucleotides (always a multiple of 3 when `productive: True`) |
| `productive` | bool or null | Acceptance under the named sequence-level policy for the observed coding span |
| `vj_in_frame` | bool or null | Junction length and junction start agree with the retained V coding frame |
| `stop_codon` | bool or null | Whether a complete observed coding codon is definitely a stop; unknown if phase is unavailable or a possible stop remains unresolved |
| `productivity_policy` | str | `canonical_observed_vdj_v1` |
| `productivity_scope` | str | `observed_coding_span` |
| `productivity_reasons` | list[str] | Stable diagnostic codes, or `["accepted"]`; JSON in TSV/CSV cells |

The policy requires original V/J landmark provenance, V Cys and J Phe/Trp
residue classes, V/J and J-tail frame agreement, and definitely stop-free
complete coding codons. It does not establish expression, splicing, folding,
or receptor function. J Phe and Trp are both accepted; exact preservation of
the reference residue is an optional separate native contract.

A known defect yields `False` even when another check is unresolved. Otherwise,
missing required evidence yields `None`, and only all accepted checks yield
`True`. For example, `NNN` can hide a stop, `GCN` is definitely sense, and `TAR`
is definitely stop. An invalid reference landmark is unresolved annotation,
not proof that the molecule is biologically nonproductive. Reference-admission
exceptions never repair landmarks or invent a translation frame.

Terminal read loss may still permit `True` when the required landmarks and
coding evidence remain observed. Loss of a required landmark yields unknown;
a structural deletion of an original landmark is a known defect. This scope
makes no claim about sequence outside the observation. `junction_aa` translates
the physical junction in complete local codons even if `vj_in_frame` is false.

`.productive_only()` conditions sampling on accepted support and asserts the
final observed assessment in both execution modes. It raises if that guarantee
cannot be met; no failed or unresolved final record is emitted. The record
validator reconstructs productivity independently of the projection helper.

## Recombination structure

### Per-segment alignment coordinates

Python record dicts use 0-based, half-open intervals. TSV defaults
to `airr_strict=True`, which converts biological coordinates to
1-based closed intervals. CSV and DataFrame exports keep Python
conventions unless `airr_strict=True` is requested.

For each of V, D, J, the Python fields are:

| Field | Meaning |
|---|---|
| `{v,d,j}_sequence_start` | Position in the coding-oriented sequence where this segment begins |
| `{v,d,j}_sequence_end` | Exclusive end position in the coding-oriented sequence |
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
| `np1` | Bases between the final V and D alignment bounds (V–J on VJ chains), including intervening N and P bases |
| `np2` | Bases between the final D and J alignment bounds, including intervening N and P bases (VDJ only; empty on VJ) |
| `np1_aa`, `np2_aa` | Complete codons within those intervals, translated in the retained V coding frame |
| `np1_length`, `np2_length` | Lengths in nucleotides |

The NP strings are slices of the coding-oriented sequence between
the final segment alignments. Their lengths equal the lengths of
those strings. Bases assigned to an adjacent germline alignment
are outside these intervals, regardless of their original N/P
provenance.

### P-nucleotide lengths

Four per-record fields carry the palindromic-insertion lengths
sampled during recombination:

| Field | Meaning |
|---|---|
| `p_v_3_length` | Number of P bases off the V allele's 3′ end (V → NP1 side) |
| `p_d_5_length` | Number of P bases off the D allele's 5′ end (NP1 → D side) - VDJ only, `0` on VJ |
| `p_d_3_length` | Number of P bases off the D allele's 3′ end (D → NP2 side) - VDJ only, `0` on VJ |
| `p_j_5_length` | Number of P bases off the J allele's 5′ end (NP2 → J side on VDJ; NP1 → J side on VJ) |

P bases contribute to the assembled molecule and to `np1` / `np2`
when they lie between the final alignment bounds. The four P
lengths record emitted recombination history and survive clonal
forks or later read loss. They do not necessarily equal the number
of P bases still observed. There is no separate P-sequence field;
alignment bounds and P lengths alone cannot recover P provenance
after arbitrary mutation, indels, or truncation.

The engine samples P-insertion lengths from
`cfg.reference_models.p_nucleotide_lengths`. See
[Junction N/P additions](../guides/junction-additions.md) for the
layout diagrams.

### CIGAR fields

| Field | Meaning |
|---|---|
| `{v,d,j}_cigar` | Local per-segment CIGAR against the generating allele named by the standard call |
| `{v,d,j}_score` | Per-segment alignment score (when reported by the projector) |
| `{v,d,j}_identity` | Fraction from 0 to 1: matching bases divided by alignment columns, including gaps |
| `{v,d,j}_support` | Per-segment support value (when reported) |
| `c_call` | Constant-region call (empty when no C-region biology is modeled) |

Python dicts carry local `M` / `I` / `D` operation runs. Strict
exports add leading query-offset `S` and reference-offset `N` runs
when nonzero. End loss shortens observed reference coverage; it does not
create alignment deletions. A structural deletion inside the observed span
still produces a `D` operation. Strict TSV/CSV also write booleans as `T` / `F`;
DataFrames retain boolean values.

## Mutation and artefact counters

### Biological SHM

```text
n_mutations             ← biological SHM substitutions that changed a base (canonical)
n_observable_mutations  ← the subset a mutation caller can recover from `sequence`
n_v_mutations           ← V-segment SHM
n_d_mutations           ← D-segment SHM
n_j_mutations           ← J-segment SHM
n_np_mutations          ← NP1 + NP2 combined SHM
mutation_rate           ← cumulative realized SHM event burden per input-molecule base
```

`n_mutations` counts substitutions that changed a base, not substitution
attempts. A uniform replacement is drawn from A/C/G/T including the base
already present, and such a draw is not a mutation: `mutate(model="uniform",
count=7)` typically reports five or six.

### Scoring a mutation caller

Use `n_observable_mutations`, not `n_mutations`. It is the number of
hypermutation substitutions still visible in `sequence`, and it is below
`n_mutations` by three intended effects:

- a mutation on an untemplated junction base (an N or P addition) changed a
  base, so it is counted, but it has no germline base to differ from, so
  nothing can recover it — about 5% of records that carry a mutation;
- two mutations can hit one position, the later one restoring the germline
  base, which leaves nothing to see;
- a position lost to end loss or a deletion is no longer in the molecule.

It is attributable to hypermutation and to nothing else. A PCR
misincorporation, a sequencing miscall or an `N` mask leaves a position
differing from its germline too, and none of them moves this count — they
are reported on `n_pcr_errors` and `n_quality_errors`, which is what a
caller's false positives should be measured against. Two scoping notes:
the count is over `sequence`, which stays the whole molecule under
`paired_end()` while `r1_sequence` / `r2_sequence` cover less; and on a
contaminant record the molecule is a background sequence, so no mutation
field on it describes hypermutation.

The four per-segment fields sum cleanly:
`n_v + n_d + n_j + n_np == n_mutations` on every record. The
[`validate_records`](../validation/validate-records.md) gate
checks this equality and fires `MutationCountSumMismatch` if
it ever breaks.

`mutation_rate` sums each mutation stage's applied events divided by that
stage's input molecule length. Later clipping leaves those contributions
unchanged. It counts repeated hits and reversions, so it is not a mismatch
fraction measured on the observed V segment. Lineage summary records instead
use net differences divided by the node's pre-observation molecule length.

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


### Per-subregion recoverable counters

Six more fields, the regional counterpart of `n_observable_mutations`:

```text
n_fwr1_observable_mutations
n_cdr1_observable_mutations
n_fwr2_observable_mutations
n_cdr2_observable_mutations
n_fwr3_observable_mutations
n_v_unannotated_observable_mutations
```

**Use these for a CDR/FWR selection-pressure ratio computed from a read.**
The event counters above are a correct partition of `n_v_mutations`, and a
ratio built from them is an *event* ratio — which is a different quantity
from what a read shows, because an event whose base was later overwritten or
lost is still an event.

Measured against an independent oracle that buckets the visible differences
by the same IMGT table, the event counters agree with the read on 60 of 60
records at one mutation per record, 32 of 60 at twenty, and 17 of 60 once
twelve nucleotides are lost from the 5' end. These six agree on 60 of 60 in
every one of those configurations.

They are derived like `n_observable_mutations` — the hypermutation passes'
own sites resolved against the finished molecule — so a position restored by
a repeat hit or removed by read loss is absent from them for the same reason
it is absent from the whole-record count. They sum to the V part of
`n_observable_mutations`, and each is at most its event counterpart;
`validate_records` checks both.

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

`r1_start/end` and `r2_start/end` are always 0-based, half-open
windows on the observed `sequence`, even with `rev_comp=True`
or strict export. `r2_sequence` is the reverse complement of
its window.

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
  checks supported counters, coordinates, calls, and productivity
  against the underlying engine `Outcome` and reference data.
  It tests consistency with the implemented model, not universal
  biological correctness or arbitrary custom metadata schemas.
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

**Scoring a mutation caller against `n_mutations`.** It counts
every substitution that changed a base, including ones on
untemplated bases that leave no germline base to differ from
and ones a later mutation reversed. `n_observable_mutations`
is the recoverable number; see
[Scoring a mutation caller](#scoring-a-mutation-caller).

**Reading `count=N` as a promise of `n_mutations == N`.** For
`model="uniform"`, `count` is the number of attempts and
roughly a quarter of them draw the base already present.

**Treating all `n_*` counters as biological mutations.** Only
the SHM partition (`n_mutations`, `n_observable_mutations`,
`n_v_mutations`, `n_d_mutations`, `n_j_mutations`,
`n_np_mutations`, the six V-subregion fields) describes biology. `n_pcr_errors`,
`n_quality_errors`, `n_indels`, and the end-loss lengths
describe library / sequencer artefacts. They live on the same
record by design - they don't share a counter.

**Confusing inferred calls with generating alleles.** Standard calls are
single generating alleles. The `recalled_*_call` fields describe inference
and can contain comma-separated ties.

**Confusing trim with end-loss.** `v_trim_3` is recombination
biology (recombinase chewed bases off the V allele's 3′ end
before junction joining). `end_loss_3_length` is sequencing
artefact (the read ran short / was 3′ end-loss-clipped). The
field names are similar; the biology is different.

**Treating `np1` / `np2` as N-only strings.** They contain the
combined gap between adjacent alignments, including intervening
P bases. Use the four `p_*_length` fields
(`p_v_3_length`, `p_d_5_length`, `p_d_3_length`, `p_j_5_length`)
to read emitted palindromic-insertion lengths. These provenance
counters do not identify which P bases remain after read corruption.

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
