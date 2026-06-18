# Clonal simulation overview

<p class="lead">GenAIRR has three clonal surfaces. Use
<code>clonal_lineage</code> when you need B-cell affinity-maturation trees,
<code>clonal_repertoire</code> when you need TCR or flat-BCR clone-size /
abundance repertoires, and legacy <code>expand_clones</code> only when you
need the older fixed-size star model. All three stamp planted clone labels so
AIRR clone-calling and ML benchmarks can compare inferred groups against the
truth the simulator created.</p>

## Choose the right clonal model

| Use case | DSL | Biology | Output truth |
|---|---|---|---|
| **BCR lineage reconstruction / affinity maturation** | [`clonal_lineage(...)`](clonal-lineage.md) | Generation-synchronous B-cell tree, per-division S5F SHM, optional sequence-distance selection, final live-cell sampling | AIRR records with `clone_id`, `lineage_*`, `duplicate_count`; one `LineageTree` per clone |
| **TCR clone-size / abundance benchmark** | [`clonal_repertoire(...)`](clonal-repertoire.md) | One rearranged T cell copied to a heavy-tailed clone size; no SHM; optional per-read technical noise | AIRR records with `clone_id` and AIRR `duplicate_count` |
| **Flat BCR abundance without genealogy** | [`clonal_repertoire(...)`](clonal-repertoire.md) | One BCR rearrangement copied to a heavy-tailed size; optional flat post-fork SHM per copy | AIRR records with `clone_id` and `duplicate_count` |
| **Legacy fixed-size star** | `expand_clones(...)` | One parent rearrangement and a fixed `per_clone` descendant count | AIRR records with `clone_id`, `parent_id`; parent `Outcome`s on `result.parents` |

`expand_clones` is still supported for old scripts, but new clone-related
benchmarks should usually start with `clonal_lineage` or `clonal_repertoire`.
Those two encode the distinction AIRR users usually care about: BCR lineages are
mutation trees, while TCR clones are abundance groups with technical read noise.

## Shared DSL shape

Clonal workflows all start by creating one founder rearrangement per clone:

```python
ga.Experiment.on("human_igh").recombine()
```

Everything before the clonal fork runs once per clone. Everything after a flat
fork (`clonal_repertoire` or `expand_clones`) runs once per emitted read/copy.
`clonal_lineage` is different: it grows the SHM tree internally, then optional
library-prep / sequencing artefact passes run once per observed cell.

Only one clonal fork is allowed in a pipeline:

```python
.clonal_lineage(...)
.clonal_repertoire(...)
.expand_clones(...)
```

are mutually exclusive.

`run_records(n=...)` is not the way to set clone output size for modern clonal
models. Use the clonal parameters instead:

| Model | Output-size knobs |
|---|---|
| `clonal_lineage` | `n_clones`, `max_generations`, `n_max`, `n_sample`, extinction/survival, genotype collapse |
| `clonal_repertoire` | `n_clones`, `size_distribution`, `max_size`, `unexpanded_fraction`, genotype collapse |
| `expand_clones` | `n_clones × per_clone` fixed descendants |

## BCR lineage trees

Use `clonal_lineage` when you need a real B-cell genealogy:

```python
import GenAIRR as ga

result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .clonal_lineage(
          n_clones=20,
          max_generations=6,
          n_max=300,
          n_sample=30,
          rate=0.01,
          lambda_base=1.6,
          selection_strength=10.0,
      )
      .sequencing_errors(rate=0.005)
      .run_records(seed=0, validate_records=True)
)

rec = result.records[0]
print(rec["clone_id"], rec["lineage_node_id"], rec["lineage_generation"])
print(rec["lineage_abundance"], rec["duplicate_count"])

tree = result.lineage_trees[rec["clone_id"]]
newick = tree.to_newick()
node_table = tree.to_node_table_tsv()
```

What happens:

1. `recombine()` creates one naive BCR founder per clone.
2. The Rust lineage engine grows a tree for `max_generations`, with live-cell
   carrying capacity `n_max`.
3. Each child receives per-division S5F SHM at `rate`.
4. If selection is enabled, offspring rates are modulated by a BLOSUM62
   sequence-distance proxy, not a physical antigen-binding model.
5. `n_sample` cells are sampled from the living final-generation population and
   identical genotypes are collapsed into `lineage_abundance` /
   `duplicate_count`.

`clonal_lineage` is BCR-only. Calling it on TCR refdata raises `ValueError`
because T cells do not somatically hypermutate.

Deep dive: [Clonal lineage trees](clonal-lineage.md).

## TCR and flat-BCR abundance repertoires

Use `clonal_repertoire` when the clone truth is membership and abundance, not a
lineage tree:

```python
import GenAIRR as ga

result = (
    ga.Experiment.on("human_tcrb")
      .allow_curatable_refdata()
      .recombine()
      .clonal_repertoire(
          n_clones=200,
          size_distribution="power_law",
          exponent=2.0,
          max_size=500,
          unexpanded_fraction=0.5,
      )
      .sequencing_errors(rate=0.005)
      .run_records(seed=0, validate_records=True)
)

for rec in result.records[:5]:
    print(rec["clone_id"], rec["duplicate_count"], rec["v_call"], rec["j_call"])
```

What happens:

1. `recombine()` creates one rearrangement per clone.
2. A clone size is drawn from a heavy-tailed distribution: rounded continuous
   power-law by default, or log-normal.
3. That many copies pass through post-fork per-read passes such as
   `sequencing_errors`, `pcr_amplify`, `polymerase_indels`, `end_loss_*`,
   `ambiguous_base_calls`, `random_strand_orientation`, or `paired_end`.
4. Identical output sequences collapse into one AIRR record whose
   `duplicate_count` is the represented abundance.

For TCR, do not add `.mutate(...)`; the API rejects it. TCR within-clone sequence
diversity should come from technical artefact passes only. For flat BCR
abundance, you may add post-fork `.mutate(...)`, but that is independent SHM off
the founder, not a tree.

Deep dive: [Clonal repertoires](clonal-repertoire.md).

## Legacy fixed-size stars

`expand_clones` remains available for old fixed-size star benchmarks:

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .expand_clones(n_clones=5, per_clone=10)
      .mutate(model="s5f", rate=0.02)
      .run_records(seed=1)
)

print(len(result.records))                         # 50
print(result.records[0]["clone_id"])               # 0
print(result.records[0]["parent_id"])              # 0
parent = result.parents[result.records[0]["parent_id"]]
```

`expand_clones` records carry `parent_id` and `result.parents` because the old
star model keeps an explicit parent `Outcome`. Modern `clonal_repertoire` does
not expose `parents`; it collapses abundance into records. `clonal_lineage`
exposes lineage truth through `lineage_*` fields and `result.lineage_trees`
instead.

## Output fields by model

| Field / object | `clonal_lineage` | `clonal_repertoire` | `expand_clones` |
|---|---|---|---|
| `clone_id` | Yes: planted family label | Yes: planted clone label | Yes: planted family label |
| `duplicate_count` | Yes: alias of `lineage_abundance` after final-cell sampling | Yes: collapsed abundance | No standard abundance field in the legacy star model |
| `lineage_node_id` / `lineage_parent_id` / `lineage_generation` | Yes | No | No |
| `lineage_abundance` / `lineage_affinity` | Yes | No | No |
| `parent_id` | No; use `lineage_parent_id` for tree parent | No | Yes |
| `result.parents` | No | No | Yes |
| `result.lineage_trees` | Yes | No | No |
| `result.outcomes` | One per observed record | One per emitted/collapsed record | One per descendant record |

For external clone-calling tools, keep the planted label under a non-AIRR name so
it does not collide with the tool's inferred `clone_id`:

```python
import pandas as pd

df = pd.DataFrame(result.records).rename(columns={"clone_id": "true_clone_id"})
df.to_csv("repertoire.tsv", sep="\t", index=False)
```

`duplicate_count` is the AIRR-standard abundance column consumed by
abundance-aware workflows. `clonal_repertoire` and `clonal_lineage` emit it
directly.

## Validation

Use record validation on every clonal workflow:

```python
result = exp.run_records(seed=42, validate_records=True)
```

or explicitly:

```python
record_report = result.validate_records(refdata)
assert record_report, record_report.summary()
```

Family validation is records-only and works across all clonal models that carry
`clone_id`:

```python
family_report = result.validate_families()
assert family_report, family_report.summary()
```

Currently `validate_families()` checks that every record in a clonal batch has a
`clone_id` and, when `truth_v_call` / `truth_d_call` / `truth_j_call` are present
from `expose_provenance=True`, that those truth calls are invariant within each
clone. It does not validate lineage topology, clone-size priors, or biological
realism.

`validate_families_with_parents(refdata)` is specific to legacy
`expand_clones`, because it compares descendant records against
`result.parents`. For `clonal_lineage`, validate the tree objects directly with
`tree.validate()` and use the exported Newick/FASTA/node table for lineage-tool
scoring.

## Ordering rules

Ancestor-phase steps go before the clonal fork:

| Step | Why |
|---|---|
| `.recombine()` | Defines the clone's V/D/J, trims, NP sequence, and junction |
| `.invert_d(...)` | Recombination-time D orientation; inherited by the clone |
| `.receptor_revision(...)` | Recombination/development-time V replacement; inherited by the clone |
| `.productive_only()` / `.restrict_alleles(...)` | Constraints on the founder draw |

Descendant/read-phase steps go after a flat fork (`clonal_repertoire` or
`expand_clones`):

| Step | Notes |
|---|---|
| `.mutate(...)` | BCR-only flat SHM; not allowed on TCR |
| `.pcr_amplify(...)`, `.sequencing_errors(...)`, `.polymerase_indels(...)` | Per-read technical artefacts |
| `.ambiguous_base_calls(...)`, `.end_loss_*prime(...)`, `.random_strand_orientation(...)` | Per-read observation artefacts |
| `.paired_end(...)` | Supported after legacy `expand_clones`; accepted after `clonal_repertoire` with abundance-collapse caveats; not yet supported after `clonal_lineage` |

For `clonal_lineage`, do not add `.mutate(...)` afterward: SHM is internal to the
tree and controlled by `clonal_lineage(rate=...)`. Library-prep and sequencing
artefact passes may follow; `paired_end` is still a future addition for lineage
output.

When `paired_end` follows `clonal_repertoire`, records are still collapsed by
assembled `sequence` and carry `duplicate_count`. TSV/DataFrame output preserves
that abundance. FASTQ exporters do not expand `duplicate_count` back into multiple
read pairs, so use this path for paired fields on collapsed records, not for exact
per-copy paired FASTQ depth.

## Common mistakes

**Using `clonal_lineage` for TCR.** TCR clones do not SHM. Use
`clonal_repertoire` for TCR clone-size and abundance benchmarks.

**Expecting exact discrete Zipf from `clonal_repertoire`.** The default
`power_law` sampler is a rounded continuous inverse-CDF draw. It is heavy-tailed
and Zipf-like, but not an exact discrete Zipf PMF.

**Expecting `parent_id` on every clonal model.** `parent_id` belongs to legacy
`expand_clones`. Use `lineage_parent_id` for BCR lineage-tree parentage, and use
`clone_id` + `duplicate_count` for `clonal_repertoire`.

**Using `n=` with modern clonal models.** `clonal_lineage` and
`clonal_repertoire` compute record counts from their own parameters and genotype
collapse. Passing `n` to `run_records` raises.

## Where to go next

- **[Clonal lineage trees](clonal-lineage.md)** - full BCR lineage model,
  affinity-selection proxy, tree exporters, and Change-O validation example.
- **[Clonal repertoires](clonal-repertoire.md)** - TCR and flat-BCR abundance
  model, clone-size parameters, `duplicate_count`, and tool export notes.
- **[Validation & reproducibility](../validation/index.md)** - record and family
  validation layers.
- **[Corruption + sequencing artefacts](corruption-sequencing.md)** - technical
  noise passes that compose with clonal workflows.
- **[Paired-end reads and FASTQ](paired-end-fastq.md)** - paired-end output;
  currently available for non-lineage clonal workflows.
