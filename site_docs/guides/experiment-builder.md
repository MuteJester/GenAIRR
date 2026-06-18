# The Experiment builder

<p class="lead">Every simulation you write in GenAIRR is built with
the <code>Experiment</code> DSL: a fluent pipeline of biological,
library-prep, and sequencing mechanisms that compiles to a
deterministic plan. This page is the control panel - what every
method does, the order they go in, and which guide to read for
each one.</p>

!!! tip "Your learning path"
    You're at the hub of the **"I want to simulate sequences"**
    path. The mechanism guides linked below this page (SHM,
    recombination, clonal families, corruption, paired-end) all
    plug into this DSL. New here? Start with the
    [Quick start](../getting-started/quick-start.md) instead;
    come back when you need a specific mechanism.
    [See all paths →](../learn.md)

## What the Experiment builder does

`Experiment` is a *builder*. Each chained method returns the same
`Experiment` object, extended by one more pass:

```python
exp = (
    ga.Experiment.on("human_igh")     # bind to a cartridge
      .recombine()                    # append V(D)J recombination
      .mutate(rate=0.05)              # append biological SHM
      .pcr_amplify(count=(0, 3))      # append PCR-error corruption
)
```

The pipeline only runs when you call `run_records(...)` (records the
output to AIRR format) or `run(...)` (returns full `Outcome` objects
for deep introspection). Until then, you're just declaring the
plan.

| Call | Returns | Use it for |
|---|---|---|
| `exp.run_records(n=..., seed=...)` | `SimulationResult` (list of AIRR dicts) | Routine simulation work |
| `exp.run(n=..., seed=...)` | `list[Outcome]` (full IR + trace + revision history) | Advanced introspection, debugging, replay |
| `exp.compile()` | `CompiledExperiment` (compiled plan, reusable across batches) | Hot loops where you want compile-once / run-many |

## A complete full-stack example

A realistic pipeline exercising biology, observation-stage
corruption, paired-end read layout, and runtime validation:

```python
import GenAIRR as ga

result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .productive_only()
      .mutate(model="s5f", rate=0.03)
      .polymerase_indels(count=1)
      .end_loss_5prime(length=(0, 3))
      .end_loss_3prime(length=(0, 3))
      .random_strand_orientation(prob=0.5)
      .paired_end(r1_length=150, insert_size=(250, 450))
      .run_records(n=100, seed=1, validate_records=True)
)

print(len(result))                                     # 100
print(result[0]["productive"], result[0]["v_call"])    # True 'IGHVF10-G38*04'
print(result[0]["read_layout"])                        # 'paired_end'
print(result[0]["r1_sequence"][:30])                   # 'gaggtgcagctg...'
```

That's a finished simulation: every record is productive, has
realistic S5F mutations, may carry an indel and/or end-loss, is
flipped strand at random, and ships with R1/R2 paired-end windows
sliced to a 150-bp read length over a 250-450-bp insert. The
`validate_records=True` flag opts the postcondition validator into
the run.

## Pipeline stage map

The Experiment surface partitions into seven groups by what the
pass actually does:

| Stage | Methods | Effect |
|---|---|---|
| **Recombination** | `.recombine()` | Sample alleles, trim, fill NP1/NP2, assemble |
| **Recombination-stage mechanisms** | `.invert_d(prob=...)`, `.receptor_revision(prob=...)` | D in reverse-complement orientation; post-recombine V replacement |
| **Constraints** | `.productive_only()`, `.restrict_alleles(v=[...], d=[...], j=[...])` | Constrain the sample space (productive triad; allele subsetting) |
| **Biological mutation** | `.mutate(model="s5f"\|"uniform", rate=..., segment_rates={...}, v_subregion_rates={...})` | Somatic hypermutation per descendant |
| **Clonal structure** | `.clonal_lineage(...)`, `.clonal_repertoire(...)`, legacy `.expand_clones(...)` | BCR trees, TCR / flat-BCR abundance repertoires, or fixed-size star families |
| **Library / sequencing corruption** | `.pcr_amplify(count=...)`, `.polymerase_indels(count=...)`, `.ambiguous_base_calls(count=...)`, `.sequencing_errors(count=...)`, `.end_loss_5prime(length=...)`, `.end_loss_3prime(length=...)` | Library-prep + sequencer artefacts |
| **Read layout** | `.paired_end(r1_length=..., r2_length=..., insert_size=...)`, `.random_strand_orientation(prob=...)` | R1/R2 windows; strand flips |
| **Bookkeeping** | `.with_metadata(experiment_id=..., tissue=...)`, `.contaminate(prob=...)` | Stamp user fields onto every record; inject background contaminants |

**Junction mechanisms (P and N additions) are NOT separate Experiment
methods.** They're driven by the cartridge's empirical models plane
`ReferenceEmpiricalModels.np_lengths` / `np_bases` /
`p_nucleotide_lengths`. Authoring those values on your cartridge
controls how much N is added and what bases get drawn; the
`.recombine()` pass consumes them automatically.

## Order matters

GenAIRR's API rejects pipelines whose steps biologically cannot
compose. The clonal methods are the main ordering boundary:

| Method | Use it for | Post-fork behavior |
|---|---|---|
| `clonal_lineage(...)` | BCR affinity-maturation trees | SHM is internal to the tree; optional library-prep / sequencing artefacts run once per observed cell |
| `clonal_repertoire(...)` | TCR or flat-BCR abundance repertoires | Copies are emitted through post-fork per-read passes and collapsed into `duplicate_count` |
| `expand_clones(...)` | Legacy fixed-size star families | Fixed `n_clones × per_clone` descendants |

For flat clonal models (`clonal_repertoire` and `expand_clones`), steps before
the fork run once per clone; steps after run once per emitted read/copy.

```python
exp = (
    ga.Experiment.on("human_igh")
      # ──── ancestor phase (runs ONCE per clonal parent) ────
      .recombine()
      .invert_d(prob=0.05)
      .receptor_revision(prob=0.05)
      # ──── fork ─────────────────────────────────────────────
      .clonal_repertoire(n_clones=50, max_size=100, unexpanded_fraction=0.3)
      # ──── descendant phase (runs ONCE per descendant) ──────
      .mutate(model="s5f", rate=0.05)
      .pcr_amplify(count=(0, 3))
      .end_loss_5prime(length=(0, 8))
      .paired_end(r1_length=150, insert_size=300)
)

result = exp.run_records(seed=42)
# Each clone shares the same V(D)J recombination + D orientation +
# receptor revision; each emitted read has independent SHM, PCR errors,
# end-loss, and R1/R2 windows. Identical reads collapse into duplicate_count.
```

**Ancestor-phase passes** (anything before a flat clonal fork) run once per
clonal parent and propagate to every emitted copy. Use them for the V(D)J
recombination event itself, recombination-stage mechanisms (D inversion,
receptor revision), and constraints on the founder draw (`productive_only`,
`restrict_alleles`).

**Descendant/read-phase passes** (anything after a flat clonal fork) run
independently per emitted copy. Use them for BCR flat SHM, all library /
sequencer corruption, and read layout. TCR rejects `.mutate(...)`; T cells do
not SHM. If `paired_end` follows `clonal_repertoire`, remember the result is
still abundance-collapsed by assembled sequence: `duplicate_count` carries copy
number, and FASTQ export does not expand it into multiple read pairs.

`clonal_lineage` handles its own tree-internal SHM through
`clonal_lineage(rate=...)`; do not add `.mutate(...)` after it. Library-prep and
sequencing artefacts may follow lineage output, but `paired_end` is not wired
through `clonal_lineage` yet.

If no clonal method is present, every pass runs per record.

## Choosing counts vs rates

Many mechanisms accept either a fixed count, a (min, max) tuple
for uniform sampling, or an empirical list of (value, weight)
pairs for arbitrary distributions. The pattern is consistent
across passes:

```python
# Fixed count - every record gets exactly this many.
exp.polymerase_indels(count=2)

# Uniform tuple - every record draws a count in [0, 3].
exp.pcr_amplify(count=(0, 3))

# Empirical list - every record draws from this empirical distribution.
exp.ambiguous_base_calls(count=[(0, 0.6), (1, 0.3), (2, 0.1)])
```

For `mutate`, `count=` is the absolute number of mutation events,
and `rate=` is a per-base SHM rate that scales with sequence
length:

```python
exp.mutate(model="s5f", rate=0.03)      # ~3% per-base S5F SHM
exp.mutate(model="uniform", count=12)   # exactly 12 uniform mutations per record
```

`end_loss_*prime(length=...)` follows the same `int` / `(low, high)` /
`[(value, weight), ...]` shape. Read [`Tune corruption
rates`](tune-corruption.md) for the per-mechanism shape reference;
the API reference page has the formal signatures.

## How `productive_only()` works

`.productive_only()` is **constraint-aware sampling, not
reject-and-retry**. Most simulators build a sequence first and
check it afterwards; if it's broken they discard and resample.
That's slow, statistically biased toward weakly-constrained
edges of the distribution, and offers no way to reason about
*which* draws are admissible. GenAIRR inverts the order: at
every draw point the engine consults the active contract bundle
and only candidates that keep it intact stay in the support.

`productive_only()` registers four predicates as the contract
bundle:

| Predicate | What it enforces | Where it fires |
|---|---|---|
| **AnchorPreserved (V)** | V's conserved Cysteine codon survives the 3′ V trim | At V-trim sampling |
| **AnchorPreserved (J)** | J's conserved Trp (heavy / κ / TR) or Phe (λ) codon survives the 5′ J trim | At J-trim sampling |
| **ProductiveJunctionFrame** | Junction length (Cys-codon start through anchor-codon end) is divisible by 3 | At NP-region length sampling |
| **NoStopCodonInJunction** | No `TAA` / `TAG` / `TGA` codon in the assembled junction | At NP base-draw sampling |

Bases the engine touches that the contracts don't gate (mutation
events, all corruption passes) can still degrade the result
downstream - that's the intentional split. Heavy SHM or 5′ / 3′
loss can still produce a non-productive *observed* sequence;
the recombination-time identity stays productive.

### Strict mode on empty support

Some configurations are over-constrained - a particular allele
pair may have no admissible NP length that keeps the junction
in-frame. The `strict=` flag on `run_records(...)` picks the
behaviour:

- **`strict=False` (default)** - fall back to the unconstrained
  distribution at that single draw and continue. Records can
  contain rare non-productive sequences when the constraint
  couldn't be satisfied.
- **`strict=True`** - raise `StrictSamplingError(pass_name,
  address, reason)` and stop. Use this for training datasets,
  validation runs, and formal benchmarks where 100 % contract
  compliance is required.

`StrictSamplingError` is **not** a `ValueError` subclass; catch
it explicitly. See the
[Validation hub](../validation/index.md#strict-vs-permissive)
for the recovery pattern.

## Validation and reproducibility

```python
result = exp.run_records(n=1000, seed=42, validate_records=True)
```

Two flags carry most of the reproducibility story:

- **`seed=`** - same seed → byte-identical output across runs and
  platforms. `n` consecutive seeded records use seeds
  `[seed, seed+1, …, seed+n-1]`, so batches stitch together if you
  offset the starting seed.
- **`validate_records=True`** - runs the postcondition validator
  inline on every record after the run. Off by default;
  release-tier loops opt in. See
  [`validate_records`](../validation/validate-records.md).

The full trace of every random draw is on each `Outcome`'s `trace()`
method (`result.outcomes[i].trace()` after `run_records`). The
trace is reproducible: replaying a recorded outcome at the same
seed produces the same draws, which is how golden tests work.

## Where to go next

| Topic | Guide |
|---|---|
| Biology → API surface lookup | [Biology map](biology-map.md) |
| Per-segment + per-V-subregion SHM targeting | [Targeted SHM rates](shm-targeting.md) |
| R1/R2 windows + insert sizes + FASTQ output | [Paired-end reads and FASTQ](paired-end-fastq.md) |
| Choosing a clonal model | [Clonal simulation overview](clonal-families.md) |
| BCR lineage trees | [Clonal lineage trees](clonal-lineage.md) |
| TCR / flat-BCR abundance repertoires | [Clonal repertoires](clonal-repertoire.md) |
| Authoring or tuning a custom reference cartridge | [Reference cartridge](../concepts/reference-cartridge.md) |
| Confirming output integrity post-run | [`validate_records`](../validation/validate-records.md) |
| Per-record AIRR field catalogue | [Your first AIRR record](../getting-started/first-airr-record.md) |
| The mental model behind the DSL | [The simulation pipeline](../concepts/pipeline.md) |
| Every public symbol's signature | [API reference](../reference/index.md) |
