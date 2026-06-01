# Biology map

<p class="lead">One page that maps every biological mechanism
GenAIRR models to the API surface that controls it, the stage
where it fires, the AIRR output fields it affects, and the guide
that explains it. Useful when you know what biology you want to
model and need to find the GenAIRR knob — or when you're reading
a record and want to know what produced a specific field.</p>

## How to use this map

You know the biology term, you want the GenAIRR surface. Find
the row in the table below, read the surface column, click the
guide for context. Going the other direction — you know the
field name and want to know which mechanism produced it — the
[AIRR record concept page](../concepts/airr-record.md) is the
better starting point.

## Mechanism map

| Biology | GenAIRR surface | Stage | Main output fields | Guide |
|---|---|---|---|---|
| **V(D)J recombination** | `.recombine()` | Recombination | `v_call`, `d_call`, `j_call`, `junction`, `np1`, `np2`, all `*_sequence_*` coords, `productive` | [Recombination + junction biology](recombination-junction.md) |
| **Exonuclease trimming** | `cfg.reference_models.trims` (cartridge) + `.trim(v_3=..., d_5=..., d_3=..., j_5=..., enabled=...)` (experiment override) | Recombination | `v_trim_3`, `d_trim_5`, `d_trim_3`, `j_trim_5` | [Recombination + junction biology](recombination-junction.md) |
| **N-addition (NP)** | `cfg.reference_models.np_lengths` + `np_bases` (cartridge) | Recombination | `np1`, `np2`, `np1_length`, `np2_length`, `np1_aa`, `np2_aa` | [Junction N/P additions](junction-additions.md) |
| **P-nucleotide insertions** | `cfg.reference_models.p_nucleotide_lengths` (cartridge) | Recombination | `p_v_3_length`, `p_d_5_length`, `p_d_3_length`, `p_j_5_length` (P bases contribute to `sequence` + `junction`, not to `np1`/`np2`) | [Junction N/P additions](junction-additions.md) |
| **D inversion** | `.invert_d(prob=...)` | Recombination editing | `d_inverted` | [D inversion + receptor revision](recombination-editing.md) |
| **Receptor revision** | `.receptor_revision(prob=...)` | Recombination editing | `receptor_revision_applied`, `original_v_call`, post-revision `v_call` | [D inversion + receptor revision](recombination-editing.md) |
| **Productivity constraint** | `.productive_only()` | Constraint-aware sampling | `productive=True` (by construction) | [Experiment builder](experiment-builder.md#how-productive_only-works) |
| **Allele restriction** | `.restrict_alleles(v=..., d=..., j=...)` | Constraint-aware sampling | Sampled `v_call`/`d_call`/`j_call` restricted to the listed alleles | [Recombination + junction biology](recombination-junction.md) |
| **Somatic hypermutation (SHM)** | `.mutate(model="s5f", rate=...)` or `.mutate(model="uniform", count=...)` | Biology — descendant phase | `n_mutations`, `mutation_rate`, `n_v_mutations`, `n_d_mutations`, `n_j_mutations`, `n_np_mutations` | [SHM and mutation targeting](shm-targeting.md) |
| **Targeted SHM** | `.mutate(..., segment_rates=..., v_subregion_rates=...)` | Biology — descendant phase | Same plus the six V-subregion counters (`n_fwr1_mutations` … `n_v_unannotated_mutations`) | [SHM and mutation targeting](shm-targeting.md) |
| **PCR substitution errors** | `.pcr_amplify(count=...)` or `.pcr_amplify(rate=...)` | Library / sequencing artefact — descendant | `n_pcr_errors`, lowercase corruption markers in `sequence` | [Corruption + sequencing artefacts](corruption-sequencing.md) |
| **Sequencing errors** | `.sequencing_errors(count=...)` or `.sequencing_errors(rate=...)` | Library / sequencing artefact — descendant | `n_quality_errors`, lowercase corruption markers in `sequence` | [Corruption + sequencing artefacts](corruption-sequencing.md) |
| **Ambiguous base calls (N)** | `.ambiguous_base_calls(count=...)` | Library / sequencing artefact — descendant | `n_quality_errors`, `N` characters in `sequence` | [Corruption + sequencing artefacts](corruption-sequencing.md) |
| **Polymerase indels** | `.polymerase_indels(count=..., insertion_prob=0.5)` | Library / sequencing artefact — descendant | `n_indels`, `n_v_indels`, `n_d_indels`, `n_j_indels` (NP indels count toward total only) | [Corruption + sequencing artefacts](corruption-sequencing.md) |
| **End-loss (5′ / 3′)** | `.end_loss_5prime(length=...)`, `.end_loss_3prime(length=...)` (or `primer_trim_*prime` aliases) | Library / sequencing artefact — descendant | `end_loss_5_length`, `end_loss_3_length` | [Corruption + sequencing artefacts](corruption-sequencing.md) |
| **Random strand orientation** | `.random_strand_orientation(prob=0.5)` | Read layout — descendant | `rev_comp` | [Corruption + sequencing artefacts](corruption-sequencing.md) |
| **Paired-end layout** | `.paired_end(r1_length=..., insert_size=...)` | Read layout — descendant | `read_layout`, `r1_sequence`, `r2_sequence`, `r1_start`, `r1_end`, `r2_start`, `r2_end`, `insert_size` | [Paired-end reads and FASTQ](paired-end-fastq.md) |
| **Clonal expansion** | `.expand_clones(n_clones=..., per_clone=...)` | Ancestor / descendant fork | `clone_id`, `parent_id` (stamped Python-side) | [Clonal families](clonal-families.md) |
| **Contamination** | `.contaminate(prob=...)` | Library / sequencing artefact — descendant | `is_contaminant` | [Experiment builder](experiment-builder.md) |
| **Sample metadata** | `.with_metadata(**fields)` | Bookkeeping — post-run | Arbitrary user-stamped columns | [Experiment builder](experiment-builder.md) |

## Stage ordering

GenAIRR's pipeline splits into five conceptual stages. The
engine enforces the ordering at compile time — any out-of-order
call raises `ValueError`.

**1. Recombination.** `.recombine()` runs the V(D)J join,
trim-and-fill, NP-region generation, and P-nucleotide
insertion. Productivity constraints (`productive_only`,
`restrict_alleles`) mask the sampling support inside this stage.

**2. Recombination editing.** `.invert_d()` and
`.receptor_revision()` edit the just-recombined molecule. Each
can fire at most once per record.

**3. Ancestor / descendant fork (clonal pipelines only).**
`.expand_clones()` partitions the pipeline. Everything before
the fork runs once per ancestor; everything after fires per
descendant.

**4. Biology — descendant phase.** `.mutate(...)` accumulates
biological SHM on top of recombination. On clonal pipelines
this fires *after* `expand_clones`; SHM is per-descendant, not
shared across the family.

**5. Library / sequencing artefacts + read layout —
descendant phase.** All corruption passes (`pcr_amplify`,
`sequencing_errors`, `ambiguous_base_calls`, `polymerase_indels`,
`end_loss_5prime`, `end_loss_3prime`, `random_strand_orientation`,
`contaminate`) plus the read-layout projection
(`paired_end`). On clonal pipelines these must come after
`.expand_clones()`; calling any of them *before* the fork raises
`ValueError`.

**Per-batch bookkeeping** (`.with_metadata(...)`) stamps the
result after every other stage has run.

The two main ordering invariants:

- **Library artefacts never precede biology.** SHM is a *biological*
  mutation; the corruption passes model the wet lab. Reversing
  the order would model SHM mutating an already-corrupted
  sequence, which doesn't match reality.
- **All descendant-phase passes follow `expand_clones`.** That's
  what makes them per-descendant. Putting them earlier would
  share their effects across the whole family.

## Cartridge-controlled vs Experiment-controlled

A clean partition between what's biology (cartridge) and what's
experimental design (Experiment).

### Cartridge-controlled (`DataConfig`)

The reference cartridge carries the immutable biological priors:

- **Allele universe** — `cfg.alleles.v`, `cfg.alleles.d`,
  `cfg.alleles.j`, `cfg.alleles.c`. The sequences and metadata
  the recombinase has access to.
- **Empirical recombination distributions** —
  `cfg.reference_models.allele_usage`,
  `cfg.reference_models.trims`,
  `cfg.reference_models.np_lengths`,
  `cfg.reference_models.np_bases`,
  `cfg.reference_models.p_nucleotide_lengths`. The per-segment
  draws the recombination pass samples from.
- **SHM kernel** — the cartridge's S5F mutability table (used
  when `.mutate(model="s5f", ...)` runs).
- **V-subregion annotations** — `cfg.alleles.v[i].subregions`.
  Required for `v_subregion_rates`.
- **Rules / anchors** — `cfg.reference_rules` (V Cys + J anchor
  expectations, allowed bases, severity).

The cartridge says **what** biology is available. The four
empirical-model planes participate in the plan signature; the
allele catalogue + rules participate in the refdata content
hash. See [Inspect manifest + build report](cartridge-manifest-report.md).

### Experiment-controlled (`Experiment`)

The `Experiment` DSL carries the experimental design:

- **Which mechanisms to enable** — the chained method calls.
  Omit a method, that mechanism doesn't fire.
- **Rates and counts** — kwargs on the methods (e.g.
  `mutate(rate=0.05)`, `pcr_amplify(count=(0, 3))`,
  `invert_d(prob=0.05)`).
- **Constraints** — `productive_only()`,
  `restrict_alleles(v=..., d=..., j=...)`. Both
  constraint-aware: they prune the sampling support at relevant
  draw points rather than rejecting after the fact.
- **Targeting overrides** — `segment_rates`, `v_subregion_rates`
  on `mutate(...)`; per-experiment `trim(v_3=..., d_5=..., ...)`
  distributions that override the cartridge defaults.
- **Clonal structure** — `expand_clones(n_clones, per_clone)`.
- **Read layout** — `paired_end(...)`,
  `random_strand_orientation(...)`.
- **Run-time flags** — `strict`, `expose_provenance`,
  `validate_records` on `run_records(...)`.

The Experiment says **how** the available biology is exercised.
Experiment knobs participate in the plan signature so any change
flips replay safety.

## What's validated

GenAIRR's two-layer validation surface — `validate_records` and
`validate_families` — covers the engine's *internal* consistency.
It does NOT validate the biological realism of your chosen
priors.

**Validated:**

- **Projection consistency** — every AIRR field is independently
  re-derived from the underlying `Outcome` events and compared
  with the projected record. Bugs in projection / live-call
  cache / counter aggregation fire here.
- **Counter partitions** — `n_v_mutations + n_d_mutations +
  n_j_mutations + n_np_mutations == n_mutations`. The six
  V-subregion counters sum to `n_v_mutations`. Indel counters
  partition correctly (with NP indels counted toward the total
  but not the per-segment partition). Mismatches fire
  `MutationCountSumMismatch` and friends.
- **Junction + productivity** — junction coordinates re-derived
  from anchor codons; `productive` flag re-derived from the
  four-clause definition; `vj_in_frame` and `stop_codon` checked.
- **Paired-end geometry** — when `read_layout == "paired_end"`,
  R1/R2 coordinates checked against `insert_size`; reads are
  consistent with their parent assembled sequence.
- **Family invariants** — `validate_families` and
  `validate_families_with_parents` assert each `clone_id` group
  agrees on `truth_v_call` / `truth_d_call` / `truth_j_call`.
  The parent-aware form additionally compares descendants
  against their actual parent `Outcome`.

**NOT validated:**

- **Biological realism of the chosen priors.** A cartridge that
  ships an unrealistic SHM rate or a wrong NP-length distribution
  will produce internally-consistent records that don't match
  real biology. Pick your priors with the same care you'd use
  for any simulator.
- **The cartridge's identity claims.** The validator trusts the
  cartridge's species / locus / reference-set declaration. If
  the cartridge is mis-labelled, the records will faithfully
  reflect the mis-labelled biology.
- **Compatibility with downstream tooling.** AIRR-strict
  coordinate conventions (`airr_strict=True` on the exporters)
  are an export-time setting, not a validator-enforced invariant.

For the canonical reproducibility surface — plan signatures,
refdata content hashes, trace replay — see
[Trace, replay, reproducibility](trace-replay.md) and
[Validation hub](../validation/index.md).

## Where to go next

- **Designing a new simulation** → start with
  [The Experiment builder](experiment-builder.md). It's the
  control panel every mechanism plugs into.
- **Matching a real dataset's empirical distributions** →
  [Estimate cartridge models from real data](estimate-cartridge-models.md).
  Allele usage, trims, NP lengths, NP bases, P lengths — all
  estimable from AIRR-like records.
- **Debugging an unexpected output** → [Validation hub](../validation/index.md).
  Start with `validate_records(refdata)`; escalate to trace
  replay if a specific record is suspect.
- **Adding a new biological mechanism** →
  [Architecture (Contributor)](../architecture/index.md). The
  audit-first workflow, the engine invariants, and the
  before-you-add checklist live there.
