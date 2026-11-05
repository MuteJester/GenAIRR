<p align="center">
  <img src="https://raw.githubusercontent.com/MuteJester/GenAIRR/master/site_docs/assets/genairr_logo.png" alt="GenAIRR" width="200">
</p>

<h1 align="center">GenAIRR</h1>

<p align="center">
  <b>Synthetic Adaptive Immune Receptor Repertoire Generator</b>
</p>

<p align="center">
  <a href="https://pypi.org/project/GenAIRR/"><img src="https://img.shields.io/pypi/v/GenAIRR.svg?logo=pypi&logoColor=white" alt="PyPI"></a>
  <a href="https://github.com/MuteJester/GenAIRR/actions/workflows/test.yml"><img src="https://github.com/MuteJester/GenAIRR/actions/workflows/test.yml/badge.svg" alt="Tests"></a>
  <a href="https://pypi.org/project/GenAIRR/"><img src="https://img.shields.io/pypi/pyversions/GenAIRR.svg?logo=python&logoColor=white" alt="Python"></a>
  <a href="https://github.com/MuteJester/GenAIRR/blob/master/LICENSE"><img src="https://img.shields.io/github/license/MuteJester/GenAIRR" alt="License"></a>
</p>

<p align="center">
  High-performance BCR and TCR sequence simulation with full ground-truth annotations.<br/>
  Rust kernel &middot; 23 species &middot; constraint-aware sampling &middot; cross-platform wheels
</p>

<p align="center">
  <a href="https://mutejester.github.io/GenAIRR/"><b>📖 Documentation</b></a>
</p>

---

## Installation

```bash
pip install GenAIRR
```

GenAIRR ships as a single wheel that bundles both the Python API and the Rust simulation kernel — no extra packages, no compiler needed. Pre-built wheels are published for **Linux** (x86_64, aarch64), **macOS** (Intel + Apple Silicon), and **Windows** (x64), supporting Python **3.9+**.

Building from source needs a stable Rust toolchain (`rustup install stable`) — see [CONTRIBUTING.md](CONTRIBUTING.md).

---

## Quick Start

```python
import GenAIRR as ga

# Generate 1,000 productive human heavy-chain sequences. Every sequence
# comes back with the full AIRR-format annotation block — gene calls,
# junction, productive flag, identity, mutation counts.
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .productive_only().run_records(n=1000, seed=42)
)

# `result` is a SimulationResult — list-like over AIRR record dicts.
# Each dict has the 50+ standard AIRR fields per row.
len(result)                 # 1000
rec = result[0]

rec["sequence"]             # 'caggtgcagctggtgcagtctggggctgaggtg...' (nucleotide)
rec["sequence_aa"]          # 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYY...' (translated)
rec["locus"]                # 'IGH'
rec["v_call"]               # 'IGHVF6-G23*01'    (generating allele)
rec["recalled_v_call"]      # post-mutation inference; comma-separated when tied
rec["d_call"]               # 'IGHD1-20*01'
rec["j_call"]               # 'IGHJ1*01'
rec["junction_aa"]          # 'CARCRTRLTRAYHITGTAHFQHW'
rec["productive"]           # True
rec["v_identity"]           # 1.0  (matches/total over the V segment)
rec["n_mutations"]          # 0
rec["determinism_version"]  # 4    (which sampling model drew this record)

# Export in any of the standard formats. TSV/FASTA/FASTQ are dependency-free;
# to_dataframe() needs pandas (pip install pandas).
result.to_tsv("repertoire.tsv")        # AIRR-spec TSV (50+ columns)
result.to_fasta("sequences.fasta")     # FASTA with v_call/j_call in the headers
result.to_fastq("sequences.fastq")     # FASTQ with illumina-shaped quality scores
df = result.to_dataframe()             # one row per record, AIRR columns
```

Python record coordinates are zero-based and half-open. `to_tsv()` converts
them to AIRR's one-based closed intervals, adds CIGAR placement offsets and
writes booleans as `T`/`F`. Use `to_tsv(..., airr_strict=False)` for the legacy
Python-style TSV. With `rev_comp=True`, biological annotations index the
coding-oriented reverse complement of `sequence`; paired-read windows index
the emitted `sequence` directly.

For paired-end sequencing pipelines, add `.paired_end(...)` to the experiment and export to two synchronized FASTQ files:

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .paired_end(r1_length=150, insert_size=300)
      .run_records(n=100, seed=1)
)
result.to_paired_fastq("reads_R1.fastq", "reads_R2.fastq")
```

Headers use the AIRR record's `sequence_id` with the universally-portable `/1` and `/2` suffix (`@seq0/1` and `@seq0/2`) — no `|`-pipe metadata that some aligners (STAR < 2.7) split on. R1 / R2 sequences are written verbatim from the AIRR `r1_sequence` / `r2_sequence` fields (R2 is already reverse-complemented at projection time; the writer does NOT apply a second flip). Default `overwrite=False` refuses to clobber existing files; pass `overwrite=True` to replace. Quality strings use the same pluggable models as `to_fastq` — `quality="illumina"` (default trapezoid shape, applied per read) or `quality="constant"` with `q=...`. The requested `insert_size` has to fit inside the molecule: on human IGH a record is roughly 350–440 nt, so an insert of 450 leaves the R2 window off the end of the sequence and `validate_records` reports `ReadWindowOutOfBounds`. See [`audit-docs/fastq_export_design.md`](audit-docs/fastq_export_design.md).

`Experiment.on(...)` accepts **a config-name string** (e.g. `"human_igh"`, `"mouse_tcrb"`), **a `DataConfig`** loaded from the bundled species pickles, or **a `RefDataConfig`** for [custom reference data](#custom-reference-data). `.productive_only()` is the constraint-aware bundle — covered in the next section. Drop it to sample without the productivity constraint; records then report true, false, or unresolved productivity according to their observed evidence.

> See the full walkthrough in the docs: [Quick Start](https://mutejester.github.io/GenAIRR/getting-started/quick-start.html) · [Your first AIRR record](https://mutejester.github.io/GenAIRR/getting-started/first-airr-record.html) · [Interpreting Results](https://mutejester.github.io/GenAIRR/concepts/airr-record.html)

---

## A realistic pipeline — everything in one place

The Experiment DSL is a fluent builder. Each step appends to the pipeline; the same `Experiment` is returned so calls chain. The example below uses every major feature GenAIRR offers — recombination, BCR clonal lineage trees with affinity maturation, per-division somatic hypermutation, primer-trimming, structural indels, PCR errors, N-base injection, custom metadata, and the productive constraint:

```python
import GenAIRR as ga

result = (
    ga.Experiment.on("human_igh")
      # 1. V(D)J recombination — sample alleles, trim, fill NP1/NP2, assemble.
      .recombine()
      # 2. Clonal lineage trees — grow 50 real B-cell lineages from the
      #    recombined founder: generation-synchronous birth–death over
      #    5 generations, per-division S5F somatic hypermutation at 0.75%
      #    per base, affinity selection, then sample 20 cells per clone.
      #    Passes BEFORE this point shape the founder rearrangement;
      #    passes AFTER apply per sampled cell (library prep). Every clone
      #    shares the founder's V(D)J but carries its own SHM + errors,
      #    and you get a ground-truth tree per clone (see lineage_trees).
      .clonal_lineage(n_clones=50, max_generations=5, n_sample=20,
                      rate=0.0075, selection_strength=10.0)
      # 3. Sequencing artefacts per sampled cell: primer trimming, structural
      #    indels, PCR substitution errors, quality-driven N injection.
      .primer_trim_5prime(length=(0, 8))
      .primer_trim_3prime(length=(0, 4))
      .polymerase_indels(count=(0, 2), insertion_prob=0.5)
      .pcr_amplify(count=(0, 3))
      .ambiguous_base_calls(count=(0, 2))
      # 4. Stamp arbitrary metadata onto every record.
      .with_metadata(experiment_id="exp001", tissue="peripheral_blood")
      # The constraint is carried, not just applied to the founder: each
      # lineage division and each corruption pass draws only substitutions,
      # indels and masks that keep the cell productive, so every emitted
      # record is productive by construction.
      .productive_only().run_records(seed=42)
)

len(result)                                  # 440   (sampled live cells across 50 trees)
sum(1 for r in result if r["productive"])    # 440   (all of them — see the note above)
len(result.lineage_trees)                    # 50    (one ground-truth tree per clone)

# Same clone, different cells — shared founder V(D)J,
# independent SHM + errors:
result[0]["clone_id"], result[1]["clone_id"]              # (0, 0)
result[0]["v_call"],   result[1]["v_call"]                # both 'IGHVF6-G21*11'
result[0]["n_mutations"], result[1]["n_mutations"]        # (11, 17) — independent SHM
result[0]["n_pcr_errors"], result[1]["n_pcr_errors"]      # (1, 0)   — independent errors

# Custom metadata propagated:
result[0]["experiment_id"], result[0]["tissue"]           # ('exp001', 'peripheral_blood')

result.to_tsv("repertoire.tsv")
newick = result.lineage_trees[0].to_newick()  # ground-truth lineage for clone 0
```

## Clonal lineages & repertoires

GenAIRR models clonal structure three ways. Pick by what you need:

| Method | Model | Use for |
|---|---|---|
| **`clonal_lineage`** | BCR affinity-maturation **trees** — generation-synchronous birth–death, per-division S5F SHM, optional affinity selection, sample + genotype-collapse | B-cell lineage trees with **ground truth** (Newick/FASTA per clone), benchmarking lineage/clone inference, ML training data |
| **`clonal_repertoire`** | Non-tree **abundance** repertoires — heavy-tailed clone sizes (power-law/lognormal) + unexpanded-singleton fraction, reads through library-prep, collapsed to `duplicate_count` | **TCR** repertoires (no SHM) and flat-BCR abundance; clone-calling benchmarks |
| `expand_clones` *(deprecated)* | Fixed-size **star**: one founder + `per_clone` independent descendants | legacy; still works and still the only way to get an exact per-clone count, but it raises a `DeprecationWarning` — prefer the two above |

```python
import GenAIRR as ga

# BCR clonal lineage trees — affinity maturation, with ground truth
bcr = (ga.Experiment.on("human_igh").recombine()
       .clonal_lineage(n_clones=50, max_generations=6, n_sample=30,
                       rate=0.01, selection_strength=10.0)   # 0 = neutral
       .sequencing_errors(rate=0.001)
       .run_records(seed=0))
bcr.records              # per-cell AIRR records: clone_id, lineage_node_id, lineage_affinity, ...
bcr.lineage_trees[0].to_newick()   # ground-truth tree (branch length = per-edge mutations)

# TCR clone-size repertoire — proliferated rearrangements, no SHM, abundance
tcr = (ga.Experiment.on("human_tcrb").allow_curatable_refdata().recombine()
       .clonal_repertoire(n_clones=200, size_distribution="power_law",
                          exponent=2.0, max_size=500, unexpanded_fraction=0.5)
       .sequencing_errors(rate=0.005)
       .run_records(seed=0))
tcr.records              # AIRR records: clone_id (membership) + duplicate_count (abundance)
```

`clonal_lineage` is BCR-only (it applies S5F SHM and rejects TCR loci); `clonal_repertoire` covers TCR and flat BCR. See the guides:
[Clonal lineage trees](https://mutejester.github.io/GenAIRR/guides/clonal-lineage.html) ·
[Clonal repertoires (TCR & abundance)](https://mutejester.github.io/GenAIRR/guides/clonal-repertoire.html).

Other feature flags worth knowing:

| Step | What it does |
|-----|-----|
| `.contaminate(prob=0.02)` | Replace ~2% of records with unrelated background sequences. |
| `.sequencing_errors(count=(0, 5))` | Lowercase 0–5 bases per sequence to mark sequencer-low-quality positions. |
| `.random_strand_orientation(prob=0.5)` | Flip ~50% of records to the reverse strand (with the `rev_comp` flag set). |
| `.restrict_alleles(v=[...], d=[...], j=[...])` | Restrict allele sampling to a specific subset — useful for benchmarking against a known repertoire. |
| `.mutate(model="uniform", rate=0.03)` | Use a uniform-rate mutation model instead of S5F. |
| `.mutate(model="s5f", rate=0.03, segment_rates={"V": 1.0, "D": 0.2, "J": 0.5, "NP": 0.0})` | Restrict SHM targeting by biological region class. Buckets are `"V"`, `"D"`, `"J"`, `"NP"` (the NP entry covers both NP1 and NP2). Omitted buckets default to `1.0`; `0.0` disables a region entirely so sites in it drop out of proposal support before contract admissibility. Default (no kwarg) is byte-identical to the pre-slice engine. Composes with `productive_only()`; replay reproduces exactly when the same rate vector is used. The rate vector itself is **not** part of cartridge identity / Rust `content_hash` — it's a per-experiment parameter (audit's documented v1 boundary). |
| `.mutate(model="s5f", rate=0.03, segment_rates={"V": 1.0, "NP": 0.0}, v_subregion_rates={"CDR": 2.0, "FWR": 0.5})` | Refine V-segment SHM targeting by IMGT subregion on top of `segment_rates`. Accepted keys are the five canonical labels `"FWR1"` / `"CDR1"` / `"FWR2"` / `"CDR2"` / `"FWR3"` plus the two aliases `"FWR"` (expands to FWR1 / FWR2 / FWR3) and `"CDR"` (expands to CDR1 / CDR2). **Alias expansion runs first, then explicit labels override** — `{"FWR": 0.5, "FWR2": 2.0}` resolves to `FWR1=0.5, FWR2=2.0, FWR3=0.5, CDR1=1.0, CDR2=1.0`. The V-site weight becomes `base × segment_rate(V) × v_subregion_rate(label)`. Non-V sites and V alleles without subregion annotations receive identity factor `1.0`, so the kwarg composes cleanly with mixed cartridges and with `productive_only()`. Default (no kwarg, `{}`, or explicit all-ones) is byte-identical to the pre-slice engine. Requires the cartridge to carry V-subregion annotations — the bundled human OGRDB cartridges do; calling `v_subregion_rates` against an unannotated cartridge raises `ValueError` at builder time. Like `segment_rates`, the rate vector is part of the plan signature (replay against a different vector fails before consuming choices) but **not** part of `refdata_content_hash`. See [`audit-docs/v_subregion_shm_rate_design.md`](audit-docs/v_subregion_shm_rate_design.md). |
| **`n_mutations` / `n_v_mutations` / `n_d_mutations` / `n_j_mutations` / `n_np_mutations`** (AIRR fields) | `n_mutations` is the **global biological SHM** counter (uniform + S5F substitutions only, sourced from `Simulation.mutation_count`). The four per-segment fields **partition** it by carried event segment — `n_v + n_d + n_j + n_np == n_mutations` by construction, enforced by the validator's `MutationCountSumMismatch` cross-check. NP1 and NP2 events roll into `n_np_mutations` (matches the `segment_rates` DSL grouping). **PCR / quality / N-corruption / receptor revision / D inversion / contaminant are intentionally excluded** — observation-stage and recombination-stage changes do not count as SHM. PCR + quality artefacts surface separately as `n_pcr_errors` / `n_quality_errors`. See [`audit-docs/mutation_provenance_audit.md`](audit-docs/mutation_provenance_audit.md). |
| **`n_fwr1_mutations` / `n_cdr1_mutations` / `n_fwr2_mutations` / `n_cdr2_mutations` / `n_fwr3_mutations` / `n_v_unannotated_mutations`** (AIRR fields) | **Per-V-subregion SHM counters** — a partition of **`n_v_mutations`** (not of the global `n_mutations`). The five canonical IMGT labels bucket V SHM events by the assigned V allele's `FWR1` / `CDR1` / `FWR2` / `CDR2` / `FWR3` interval; `n_v_unannotated_mutations` catches everything else under V. Three "unannotated" cases: (1) the user's cartridge has at least one V allele without subregion annotations (legacy / hand-authored), (2) the SHM event lands in the **V-side CDR3 contribution stretch** between `FWR3.end` and `len(V_allele.seq)` — the five canonical labels deliberately stop at FWR3, so junctional V-tail SHM goes here, (3) an SHM pass runs after an indel pass and mutates an inserted V base (non-canonical pass order). On bundled human OGRDB cartridges the unannotated bucket is a small but non-zero baseline (case 2 — V-tail / CDR3 contribution). Partition invariant: `n_fwr1 + n_cdr1 + n_fwr2 + n_cdr2 + n_fwr3 + n_v_unannotated == n_v_mutations` on every record, enforced by the validator's six `N<Region>MutationsMismatch` checks plus the `VSubregionMutationCountSumMismatch` cross-field invariant. Two-bucket aggregates (`n_cdr_mutations` / `n_fwr_mutations`) are deliberately NOT exposed — derive them downstream as `df["n_cdr_mutations"] = df["n_cdr1_mutations"] + df["n_cdr2_mutations"]`. See [`audit-docs/v_subregion_mutation_counters_audit.md`](audit-docs/v_subregion_mutation_counters_audit.md). |
| `.invert_d(prob=0.05)` | Heavy-chain only. With probability `prob`, sample the D allele in reverse-complement orientation (V(D)J inversion event). Records the decision under `sample_allele.d.inverted`; the assembled D bytes are the WC complement of the original allele; AIRR records carry `d_inverted: bool` for provenance. |
| `.receptor_revision(prob=0.05)` | Heavy-chain only. With probability `prob`, replace V after initial recombination with a different germline allele (B-cell receptor revision). Records `receptor_revision.applied` (+ replacement allele id / 3' trim on `True`); the V slice in the assembled pool is rewritten; AIRR records carry `receptor_revision_applied: bool` and `original_v_call: str` (the pre-revision V; empty when no revision happened). `v_call` continues to report the post-revision identity. |
| `.paired_end(r1_length=150, insert_size=300)` | Both VDJ and VJ. Add Illumina-style paired-end read-layout projection to every record: AIRR fields `r1_sequence` (forward window) + `r2_sequence` (reverse-complement of the 3' window) + `r1_start/end` / `r2_start/end` / `insert_size` / `read_layout="paired_end"`. Optional `r2_length` (defaults to `r1_length`). Each length accepts `int` / `(low, high)` / empirical `[(value, weight)]`. Trace records `paired_end.r1_length` / `.r2_length` / `.insert_size`; AIRR builder reads them back to populate the windows. Projection-only — no IR mutation, no live-call invalidation, runs after end-loss + rev-comp. |
| `compile()` then `compiled.run_records(...)` | Compile the plan once, reuse it across many batches — see [Compile once](#compile-once-run-many-times). |

---

## Per-individual genotypes

A `Genotype` is one subject's diploid germline complement — which V/D/J alleles
they carry, on which chromosome, with which genes deleted or duplicated. Attach
one and recombination becomes **haplotype-phased**: the V, D and J of a
rearrangement are drawn from a single chromosome, so the linkage that
haplotype-inference methods exploit is actually there in the data.

```python
import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR.genotype import Genotype

cfg = gdata.HUMAN_IGH_OGRDB
g = (
    Genotype.from_dataconfig(cfg)
      .complete_from_reference("homozygous_first_reference")
      .heterozygous("IGHVF1-G1", "IGHVF1-G1*01", "IGHVF1-G1*02")
      .delete_gene("IGHVF3-G7", haplotype="both")
      .with_subject("DONOR01")
)

result = (
    ga.Experiment.on(cfg)
      .with_genotype(g)              # recombination is now haplotype-phased
      .recombine()
      .run_records(n=1000, seed=7)
)

result[0]["subject_id"]              # 'DONOR01' — on every record
result[0]["haplotype"]               # 0 or 1 — which chromosome this read used

# The genotype leaves the library with the records it explains:
result.to_genotype_tsv("genotype_truth.tsv")
# subject_id  locus  segment  gene  allele  chromosome  novel  gene_state
```

`to_genotype_tsv()` writes one row per carried allele, so a dataset you hand to
someone else arrives with the truth it is meant to be scored against. A deleted
gene keeps its row with the literal allele `(none)`; a run with no genotype, or a
genotype that leaves reference genes unspecified, raises rather than writing a
file that looks complete and is not. Whether a given read drew a **novel** allele
is a join of `truth_v_call` against that table on `(subject_id, allele)` —
`run_records(..., expose_provenance=True)` adds the `truth_*_call` columns.

`run_cohort(genotypes, n_per_subject=...)` runs many subjects in one call and
namespaces their records; `CohortResult.to_genotype_tsv(path)` writes the whole
cohort's truth table.

> Docs: [Genotypes](https://mutejester.github.io/GenAIRR/guides/genotype.html) · [Genotype sampling & population priors](https://mutejester.github.io/GenAIRR/guides/genotype-priors.html) · [Genotype cohorts](https://mutejester.github.io/GenAIRR/guides/genotype-cohorts.html) · [Benchmarking genotype inference](https://mutejester.github.io/GenAIRR/guides/genotype-benchmarking.html)

---

## Constraint-aware sampling

GenAIRR's signature feature is **constraint-aware sampling**: contracts that prune the candidate distribution at sample time, not retries after the fact. The canonical bundle is `productive()` (preserved coding frame, no stops across the observed coding sequence, and intact V/J anchors):

```python
import GenAIRR as ga

# Every sequence is productive by construction. No retry loops, no
# post-hoc filtering — the engine only ever picks NP lengths, NP bases,
# and mutation substitutions that satisfy the bundle.
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .productive_only().run_records(n=1000, seed=42)
)
assert all(rec["productive"] for rec in result)
```

The named policy is `canonical_observed_vdj_v1`, scoped to the observed coding
span. It requires V Cys and J Phe/Trp original landmarks, consistent frame, and
definitely stop-free complete codons. Missing evidence or possible hidden stops
remain nullable in unconstrained records. Reference admission exceptions never
change this policy. Both execution modes require an accepted final assessment;
this is a sequence-level prediction, not a claim of receptor function.

`.productive_only()` supports fixed and mixed nonzero P-nucleotide lengths.
The engine keeps a candidate only when the remaining allele, trim, N/P, and
D-orientation choices can complete a productive sequence. It renormalizes the
original weights among those candidates at each draw. Fixed P lengths are
preserved; impossible support raises an error even in permissive mode instead
of dropping requested additions. This conditions each local draw, rather than
sampling from a globally conditioned joint probability law.

> Docs: [Sample only productive sequences](https://mutejester.github.io/GenAIRR/guides/productive.html)

### Strict vs permissive mode

Initial allele/genotype assignment and P-aware productive completion fail explicitly when no admissible support exists, including in permissive mode. Other passes can use their documented no-op or sentinel behavior; `strict=True` surfaces those sampling failures as exceptions:

```python
import GenAIRR as ga

try:
    ga.Experiment.on("human_igh").recombine().productive_only().run_records(n=10, seed=42, strict=True)
except ga.StrictSamplingError as e:
    pass_name, address, reason = e.args
    # pass_name e.g. "generate_np.np1", address e.g. "np.np1.length",
    # reason in {"empty_admissible_support", "support_unavailable", ...}
    print(f"{pass_name} could not satisfy the contract at {address}: {reason}")
```

---

## Reproducibility

The promise is narrow and worth stating exactly: **one seed and one
determinism version give one answer.** Two runs at the same
`GenAIRR.DETERMINISM_VERSION` agree byte for byte; two runs at different
versions promise each other nothing — not "mostly the same", not
"statistically equivalent". Every exported record carries the number as a
`determinism_version` column, and so does every trace file, so a dataset says
which sampling model drew it.

```python
import GenAIRR as ga

ga.DETERMINISM_VERSION          # 4 — the installed engine's sampling model

# Same seed + same determinism version → byte-identical records.
a = ga.Experiment.on("human_igh").recombine().run_records(n=100, seed=42)
b = ga.Experiment.on("human_igh").recombine().run_records(n=100, seed=42)
assert a[0]["sequence"] == b[0]["sequence"]
assert a[0]["determinism_version"] == ga.DETERMINISM_VERSION

# `n` runs use seeds [seed, seed+1, ..., seed+n-1], so a seed names a
# *record*, not a run: consecutive batches stitch together by offsetting
# the starting seed by at least `n`.
batch_a = ga.Experiment.on("human_igh").recombine().run_records(n=100, seed=0)
batch_b = ga.Experiment.on("human_igh").recombine().run_records(n=100, seed=100)
# batch_a[50] is byte-equal to a one-off run at seed=50.
# Offsetting by less than `n` overlaps: seeds 0, 1, 2 are not three
# independent replicates, they are one dataset read from three
# starting points.
```

The reproducibility goldens run in CI on **Linux, macOS and Windows across
Python 3.10–3.13** and assert that the same seed reproduces the same bytes on
each. Linux aarch64 ships a wheel but is not covered by those goldens.

> Docs: [Trace, replay, reproducibility](https://mutejester.github.io/GenAIRR/guides/trace-replay.html)

---

## Validation posture

GenAIRR's release readiness rests on seven guardrails, each backed by
an audit doc and a golden-test file. Every PR runs the full suite;
every release additionally runs the golden trace compatibility
fixtures so on-disk formats stay byte-stable.

- **Productive validity** — junction-frame, no-stop, anchor-preserved bundle, with strict / permissive failure modes pinned per pass.
- **Provenance correctness** — indels, end-loss, allele-call ambiguity, junction fields, and per-pass event ledger each have isolated-scenario tests.
- **Reproducibility** — every audit slice includes trace-replay round-trips that reproduce sequence, AIRR coordinates, and per-pass event counts.
- **Distribution invariants** — Monte-Carlo tests (±5σ) prove the constrained samplers draw from `natural_weight × admissibility`, with explicit negative controls against renormalization bugs.
- **Performance budgets** — wall-time regression guards on nine representative workloads catch ~10× slowdowns before they ship.
- **Postcondition validator** — `result.validate_records(refdata)` runs the engine's own truth oracle over every projected AIRR record (re-derives counters, junction, allele tie-set, structural coords from outcome state). Returns a `ValidationReport` you can `assert` on as a one-line CI guard. **This is the public output-correctness contract.**
- **Live-call cache parity** — `outcome.check_live_call_cache_parity(refdata)` compares the cached `SegmentLiveCall` on the final simulation against a from-scratch recompute, per V/D/J. **This is the internal cache-correctness check on the state that feeds projection.**

Two independent layers, each one a single call:

```python
result = exp.run_records(n=1000, seed=0)

# Layer 1 — public AIRR output correctness.
report = result.validate_records(refdata)
assert report, report.summary()

# Layer 2 — internal live-call cache correctness.
for outcome in result.outcomes:
    for p in outcome.check_live_call_cache_parity(refdata):
        assert p["tie_set_matches"], p
```

**Troubleshooting rule.** If a CI run has both layers failing on the same batch, fix the **parity** divergence first: a stale cache leaks into projection and produces spurious validator failures downstream. Once parity is green, rerun the validator — remaining failures then point at a genuine projection-layer bug rather than cache fallout.

The navigable index — guarantees → audit docs → test files — lives in
[`audit-docs/validation_matrix.md`](audit-docs/validation_matrix.md). Use it to
locate the right test for a change you're making, or to scope a
follow-up that touches an open drift item in any audit's §6.

```bash
make validate-fast      # correctness only (~60s)
make validate-full      # + performance budgets (regression guard)
make validate-release   # + golden trace compat + wheel build
```

---

## Compile once, run many times

For a hot loop, `compile()` once and reuse the plan. Contracts (`respect=`) are baked into the compiled plan, so they only need to be passed once:

```python
import GenAIRR as ga

compiled = (
    ga.Experiment.on("human_igk")
      .recombine()
      .productive_only().compile()
)

# Run 10 batches of 100, seeded so they don't overlap.
for batch in range(10):
    result = compiled.run_records(n=100, seed=batch * 100)
    result.to_tsv(f"batch_{batch:02d}.tsv")
```

---

## What you get back

`.run_records(...)` returns a `SimulationResult` — a list-like wrapper around a batch of AIRR record dicts:

| Method / attribute | Returns | Description |
|-----|-----|-----|
| `len(result)` | `int` | Number of records in the batch. |
| `result[i]` | `dict` | The i-th AIRR record. Standard 0-based indexing + slicing. |
| `for rec in result:` | iterates `dict`s | Records in `[seed, seed+1, …, seed+n-1]` order. |
| `result.records` | `list[dict]` | The underlying list. Mutate-through is fine. |
| `result.to_tsv(path, *, airr_strict=True)` | — | TSV with 1-based closed biological coordinates, CIGAR placement offsets, and `T` / `F` booleans. |
| `result.to_csv(path, *, airr_strict=False)` | — | Comma-separated; keeps Python record conventions unless `airr_strict=True`. |
| `result.to_fasta(path, *, prefix="seq")` | — | FASTA. Headers include `v_call` and `j_call`. |
| `result.to_fastq(path, *, quality="illumina", **kw)` | — | FASTQ. Quality models: `"illumina"` (smoothed trapezoid) or `"constant"`. |
| `result.to_dataframe(*, airr_strict=False)` | `pandas.DataFrame` | One row per record. Requires pandas (`pip install pandas`). |
| `result.to_genotype_tsv(path)` | — | The genotype truth sidecar: one row per carried allele. Raises when the run carried no genotype. |
| `result.outcomes` | `list[Outcome] \| None` | The underlying `Outcome` objects, for advanced introspection (see below). |

Each record dict has 50+ AIRR fields. The most commonly used:

| Field | Example value | Description |
|-----|-----|-----|
| `sequence` | `'gaggtgcagctggtg…'` | Assembled nucleotide sequence (uppercase + lowercase corruption markers). |
| `sequence_aa` | `'EVQLVESGGG…'` | Codon-rail translation. Stops emit `*`, ambiguous codons emit `X`. |
| `locus` | `'IGH'` | Locus code derived from `v_call` / `j_call`. |
| `v_call` / `d_call` / `j_call` | `'IGHV3-23*01'` | Generating alleles used by the truth alignments. |
| `recalled_v_call` / `recalled_d_call` / `recalled_j_call` | `'IGHV3-23*01,IGHV3-23*02'` | Post-mutation inferred calls, including ties. |
| `junction` / `junction_aa` | `'TGC…GAC'` / `'CAR…D'` | Junction nucleotide + AA. AA includes the V Cys (anchor) through J W/F+3. |
| `productive` | `True` / `False` / `None` | In-frame junction AND no stop codons AND anchors preserved. `None` when undefined (e.g. junction not present). |
| `v_identity` / `d_identity` / `j_identity` | `0.987` | Fraction from 0 to 1: matching bases divided by alignment columns, including gaps. |
| `v_cigar` / `d_cigar` / `j_cigar` | `'17D279M'` | Local per-segment M/I/D runs in Python dicts; strict exports add leading query `S` and reference `N` offsets when nonzero. |
| `n_mutations` / `n_pcr_errors` / `n_quality_errors` / `n_ambiguous_bases` / `n_indels` | `4` / `0` / `2` / `3` / `1` | Per-record event counts from the trace. `n_quality_errors` comes from `sequencing_errors()`; `ambiguous_base_calls()` reports on `n_ambiguous_bases`. |
| `n_observable_mutations` | `3` | The SHM substitutions a mutation caller can **recover from `sequence`**. Below `n_mutations`: a mutation on an untemplated junction base has no germline base to differ from, a second hit can restore the germline base, and a position lost to end loss is gone. **Score a caller against this field, not `n_mutations`.** |
| `determinism_version` | `4` | Which sampling model drew the row. Two datasets with different numbers are not comparable. |

The full schema (plus the `*_sequence_start/end`, `*_alignment_start/end`, `*_germline_start/end` coordinate fields, `vj_in_frame`, `stop_codon`, `rev_comp`, and others) is documented at [The AIRR record](https://mutejester.github.io/GenAIRR/concepts/airr-record.html).

### Advanced: full pipeline state via `Outcome`

When you need step-by-step IR history or the raw trace of every random draw — debugging an engine bug, building a custom alignment tool, replaying a specific seed — use `.run()` instead of `.run_records()`. It returns a list of `Outcome` objects that carry the full pipeline state:

| Accessor | Returns | Description |
|----------|---------|-------------|
| `outcome.final_simulation()` | `Simulation` | End-of-pipeline IR snapshot. |
| `outcome.revision(i)` | `Simulation` | IR after the i-th pass — full step-by-step history. |
| `outcome.revision_after(name)` | `Simulation \| None` | First revision produced by the named pass. |
| `outcome.pass_names()` | `list[str]` | Names of every pass that ran, in order. |
| `outcome.trace()` | `Trace` | Addressed log of every random draw. |

Each `Simulation` exposes `len(sim)` (pool length), `sim.bases() → bytes`, `sim.regions() → list[Region]`, `sim.germline_position(i)`, `sim.v_allele_id() / .d_allele_id() / .j_allele_id()`. Each `Region` carries `segment` (`"V"`/`"D"`/`"J"`/`"NP1"`/`"NP2"`), `start`/`end`/`len()`, `frame_phase`, and `amino_acids() → bytes` (codon-rail translation, including stop markers and ambiguous codons).

`outcome.trace()` supports `find(address)`, `prefix_query(prefix)`, and `prefix_count(prefix)` — every random draw is keyed by a hierarchical address (`"sample_allele.v"`, `"np.np1.length"`, `"np.np1.bases[3]"`, …). This is the same trace the engine uses internally for replay determinism.

`.run_records(...)` also exposes these via `result.outcomes[i]` — so you can have both the AIRR records *and* the deep introspection from a single call.

> Docs: [The Experiment builder](https://mutejester.github.io/GenAIRR/guides/experiment-builder.html) · [Trace, replay, reproducibility](https://mutejester.github.io/GenAIRR/guides/trace-replay.html) · [The AIRR record](https://mutejester.github.io/GenAIRR/concepts/airr-record.html)

---

## Supported Species & Chains

GenAIRR ships with **106 built-in configurations** covering 23 species (sourced from IMGT and OGRDB).

```python
import GenAIRR as ga
print(ga.list_configs())  # all available configs
```

| Species | BCR | TCR |
|---------|-----|-----|
| Human | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| Mouse | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| Rat | IGH, IGK, IGL | &mdash; |
| Rabbit | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| Dog | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| Cat | IGK, IGL | TCRA, TCRB, TCRD, TCRG |
| Rhesus | IGH, IGK, IGL | TCRA, TCRB, TCRD, TCRG |

<details>
<summary>All 23 species</summary>

Alpaca, Cat, Chicken, Cow, Cynomolgus, Dog, Dromedary, Ferret, Goat, Gorilla,
Horse, Human, Mouse (generic + C57BL/6J), Pig, Platypus, Rabbit, Rat, Rhesus,
Salmon, Sheep, Trout, Zebrafish.

</details>

```python
import GenAIRR as ga

ga.Experiment.on("rabbit_tcrb").recombine().run_records(n=500)

# 80 of the 106 bundled catalogues carry unresolved anchor or
# pseudogene/ORF findings, so strict compilation quarantines them. Opt in
# explicitly to simulate from one — the records are then drawn from a
# catalogue whose landmarks the project has not been able to confirm.
ga.Experiment.on("mouse_igh").allow_curatable_refdata().recombine().run_records(n=500)
ga.Experiment.on("rhesus_igk").allow_curatable_refdata().recombine().run_records(n=500)
```

The 26 catalogues that compile strictly today are the three human OGRDB
cartridges (`human_igh`, `human_igk`, `human_igl`), `HUMAN_IGH_EXTENDED`, and
22 IMGT catalogues whose anchors all resolved. Everything else — including
every mouse catalogue and `human_tcrb` — needs `allow_curatable_refdata()`.
The compile error lists the findings allele by allele, so you can also filter
the catalogue or curate it with `Experiment.curate_refdata(...)` instead.

> Docs: [Choose your path](https://mutejester.github.io/GenAIRR/learn.html) · [Reference evidence and quarantine](https://mutejester.github.io/GenAIRR/cartridges/reference-evidence.html) · [API Reference](https://mutejester.github.io/GenAIRR/reference/index.html)

---

## Custom reference data

For non-builtin alleles (custom IMGT pulls, in-house references, etc.) you can build a `RefDataConfig` directly and pass it to `Experiment.on(...)`:

```python
import GenAIRR as ga

# Real bytes, so the example runs as written. Substitute your own.
igk = ga.HUMAN_IGK_OGRDB
v1, v2 = [a for g in igk.v_alleles.values() for a in g][:2]
j1 = next(a for g in igk.j_alleles.values() for a in g)

cfg = ga.RefDataConfig.vj()
cfg.add_v_allele("v_custom*01", "v_custom", v1.ungapped_seq.upper().encode(), anchor=v1.anchor)
cfg.add_v_allele("v_custom*02", "v_custom", v2.ungapped_seq.upper().encode(), anchor=v2.anchor)
cfg.add_j_allele("j_custom*01", "j_custom", j1.ungapped_seq.upper().encode(), anchor=j1.anchor)

result = (
    ga.Experiment.on(cfg)
      .recombine(np1_lengths=[(n, 1.0) for n in range(7)])  # no empirical NP model
      .trim(enabled=False)                                  # no trim distributions
      .run_records(n=100, seed=42)
)
result[0]["v_call"], result[0]["junction_aa"]   # ('v_custom*01', 'CLQHNSYPPWRFTF')
```

`RefDataConfig.vdj()` builds a heavy-chain-shaped refdata (with a D pool); `add_d_allele(...)` populates it. Anchors are 0-based offsets of the V Cys / J W or F codon's first base, used to keep the junction frame-aligned during recombination — the engine checks that the codon at the anchor really translates to Cys (V) or Trp/Phe (J) and refuses the reference otherwise. Sequence bytes must be A/C/G/T/N; a raw `RefDataConfig` carries no empirical NP-length or trim distributions, so pass them explicitly (as above) or accept the two `UserWarning`s that say the engine fell back to a uniform NP length and a no-op trim.

For non-human or otherwise non-standard references — custom J anchor amino acids, extended sequence alphabets, hand-authored NP-length or trim distributions, pseudogene-aware allele filtering — configure a **reference cartridge** instead of relying on the built-in locus defaults. The cartridge has four typed planes (identity, catalogue, rules, empirical models) plus an explicit curation policy. See [`audit-docs/reference_cartridge.md`](audit-docs/reference_cartridge.md) for the model, the authoring surface (`ReferenceRulesSpec`, `ReferenceEmpiricalModels`, `Experiment.curate_refdata`), and end-to-end examples.

Cartridges can also carry **V-region substructure annotations** — per-V-allele IMGT `FWR1` / `CDR1` / `FWR2` / `CDR2` / `FWR3` intervals derived from the IMGT-gapped reference sequence. The bundled OGRDB cartridges populate them at load time for all human IGH / IGK / IGL V alleles; user cartridges can override the dict explicitly. Subregion intervals are inspectable on the bridged `RefDataConfig`, surfaced in `DataConfig.cartridge_manifest()["models"]["shm"]["v_subregion_support"]`, and folded into `refdata_content_hash` so two cartridges that differ only in subregion boundaries hash differently. **They also drive SHM targeting** via `Experiment.mutate(v_subregion_rates={…})` — see the targeted-SHM example in the pipeline-steps table above — and they bucket the per-region AIRR mutation counters (`n_cdr1_mutations` …) and their recoverable counterparts (`n_cdr1_observable_mutations` …), both described in the pipeline-steps table. See [`audit-docs/v_region_substructure_audit.md`](audit-docs/v_region_substructure_audit.md) and [`audit-docs/v_subregion_shm_rate_design.md`](audit-docs/v_subregion_shm_rate_design.md) for the audits and per-slice scope.

### Typed NP base models (uniform / empirical / Markov)

Cartridges can also author **per-NP-region base distributions** via the typed `ReferenceEmpiricalModels.np_bases` plane. Three `kind` values are supported end-to-end:

- `"uniform"` — equivalent to the engine default (4-way A/C/G/T, byte-identical to the pre-typed-model baseline).
- `"empirical_first_base"` — position-independent weighted categorical. Every NP slot samples from the same `first_base` distribution.
- `"markov"` — 1-step previous-base-conditional sampling. Position 0 uses `first_base`; positions 1+ select a transition row keyed by the previous emitted base.

```python
from GenAIRR.reference_models import (
    NpBaseModelSpec,
    ReferenceEmpiricalModels,
)

cfg.reference_models = ReferenceEmpiricalModels(
    np_bases={
        "NP1": NpBaseModelSpec(
            kind="markov",
            first_base={"G": 1.0},
            transitions={
                "A": {"T": 1.0},
                "C": {"G": 1.0},
                "G": {"A": 1.0},
                "T": {"C": 1.0},
            },
        )
    }
)
```

Both rows must cover every canonical from-base — the Python validator rejects partial matrices at construction time. Plan signatures fold the full Markov payload (first-base row + 4 transition rows), so replay against a different matrix fails the signature gate before any choice is consumed; same-cartridge + same-seed replay is byte-identical. Productive-only sampling composes through the existing admit-mask intersection.

**Legacy `NP_transitions` / `NP_first_bases` are NOT auto-lifted yet.** The bundled cartridges still carry those dicts as documented orphan fields, and the manifest reports `legacy_fallback=False`. Authors who want Markov today populate `ReferenceEmpiricalModels.np_bases` explicitly; an opt-in auto-lift slice is a separate cartridge-migration decision. See [`audit-docs/np_markov_base_generator_design.md`](audit-docs/np_markov_base_generator_design.md) and [`audit-docs/junction_n_addition_audit.md`](audit-docs/junction_n_addition_audit.md) for the full design.

### Templated P-nucleotide (palindromic) additions

For an opt-in model tied to the **original coding tips**, configure
`.process_coding_ends(v_3=[(-1, 1.0)], reference_ends_are_coding_ends=True)`
after `.recombine()`. Positive signed outcomes delete reference bases,
negative outcomes retain original-tip P extensions, and omitted ends are
fixed at zero. This replaces automatic trim/P settings and cannot be combined
with an explicit `.trim()` call. The boundary confirmation requires complete
junction-facing coding ends in the reference. These are caller-supplied net
outcome probabilities, not calibrated hairpin-opening probabilities.
See [original coding-end processing](site_docs/guides/recombination-junction.md#processing-original-coding-ends-opt-in)
for VDJ orientation, provenance, and productive-sampling semantics.

The existing cartridge P model below remains available when
`process_coding_ends()` is not configured.

Cartridges can also author **per-end P-nucleotide length distributions** via `ReferenceEmpiricalModels.p_nucleotide_lengths`, keyed by junction side label (`"V_3"`, `"D_5"`, `"D_3"`, `"J_5"`). P-bases are **templated palindromic complements** of the source allele's post-trim coding flank — only the per-end length is sampled; the bytes themselves derive deterministically from `(allele, trim, orientation, length)` via `complement_base`. The four ends interleave with NP1 / NP2 in pool order between V/D and D/J coding bases:

```python
from GenAIRR.reference_models import (
    EmpiricalDistributionSpec,
    ReferenceEmpiricalModels,
)

cfg.reference_models = ReferenceEmpiricalModels(
    p_nucleotide_lengths={
        "V_3": EmpiricalDistributionSpec([(0, 0.8), (1, 0.2)]),
        "J_5": EmpiricalDistributionSpec([(0, 0.9), (1, 0.1)]),
    }
)
```

Empty dict (the bundled-cartridge default) means the pipeline omits every P-addition pass — byte-identical to the pre-slice baseline. VJ cartridges reject D-end keys at the validation layer (no D segment to extend); VDJ cartridges accept all four ends. Each end's length distribution folds into the plan signature, so replay against a different P-length distribution fails the signature gate before any choice is consumed.

P-nucleotide provenance is surfaced on every AIRR record via four int fields:

```text
p_v_3_length, p_d_5_length, p_d_3_length, p_j_5_length
```

…each counting the number of templated P-bytes emitted at that V(D)J coding-end junction side. Zero on records from cartridges without a typed P-plane. The validator's `PLengthMismatch` issue kind catches downstream record tampering; cache parity is unaffected (P-bytes don't move structural region boundaries).

**Legacy `p_nucleotide_length_probs` is NOT auto-lifted yet.** The bundled cartridges still carry the orphan dict; the manifest reports `legacy_fallback=False`. Per-base P strings (`p_v_3`, …) and an aggregate `n_p_nucleotides` are still not exported, but each NP region now reports how its bases arose: `np1_n_tdt` / `np1_n_p` / `np1_n_inserted` (and the `np2_*` triple) partition `np1_length` into TdT N additions, templated P bases and bases an indel inserted. See [`audit-docs/p_nucleotide_design.md`](audit-docs/p_nucleotide_design.md) for the full design.

### Build a cartridge from FASTA

Beyond the 106 bundled cartridges, the 4-path manual construction surfaces, and the typed `ReferenceEmpiricalModels` / `ReferenceRulesSpec` authoring layer, GenAIRR ships an audit-trail builder for new cartridges from raw FASTA. **`ReferenceCartridgeBuilder`** stages FASTA parsing, identity attribution, V-subregion derivation, and rules/models attachment into a single fluent chain; every stage writes a structured entry to the build report so downstream consumers can audit how the cartridge was constructed.

```python
import GenAIRR as ga

builder = (
    ga.ReferenceCartridgeBuilder
    .from_fasta(
        v_fasta="v.fa",
        d_fasta="d.fa",
        j_fasta="j.fa",
        chain_type="BCR_HEAVY",
    )
    .infer_identity(
        species="human",
        locus="IGH",
        reference_set="custom",
        name="my_igh",
    )
    .infer_v_subregions()
)

cfg = builder.build()
print(cfg.cartridge_manifest())
print(builder.report().to_dict())
```

`.from_fasta(...)` accepts file paths, open text files, or raw FASTA strings (string with `\n` and leading `>`). Duplicate allele names, empty sequences, and constructor failures land in `report.rejected` with structured `{stage, segment, allele_name, reason}` entries rather than corrupting the cartridge silently. `.infer_v_subregions()` reuses the same `compute_v_region_boundaries` helper the bundled-cartridge bridge uses, so a builder-produced cartridge with IMGT-gapped V FASTA gets the five canonical subregion intervals for free.

`.build()` finalises the cartridge: it stamps the canonical checksum onto `schema_sha256`, attaches the `build_report` as a typed `CartridgeBuildReport` dataclass, runs `verify_integrity()`, and returns a plain `DataConfig` you can drop straight into `ga.Experiment.on(cfg)`. The build report's `.to_dict()` round-trips through `json.dumps` so CI artifacts can capture the full provenance.

**Estimate allele usage from observed rearrangements.** The first data-derived estimator landed: `ReferenceCartridgeBuilder.estimate_allele_usage(rearrangements, *, min_count=1.0, ambiguous="fractional", replace=True)` reads AIRR `v_call` / `d_call` / `j_call` columns from a `list[dict]`, AIRR-C TSV path, or open text handle, and writes per-segment allele weights into the typed `ReferenceEmpiricalModels.allele_usage` plane.

```python
builder = (
    ga.ReferenceCartridgeBuilder
    .from_fasta(v_fasta="v.fa", d_fasta="d.fa", j_fasta="j.fa", chain_type="BCR_HEAVY")
    .infer_identity(species="human", locus="IGH", reference_set="custom", name="my_igh")
    .estimate_allele_usage(rearrangements, ambiguous="fractional")
)
cfg = builder.build()
```

Recombination on the resulting cartridge **uses the estimated weights by default** — `ga.Experiment.on(cfg).recombine()` now samples alleles according to the typed plane unless an explicit `recombine(v_allele_weights=..., d_allele_weights=..., j_allele_weights=...)` kwarg is passed (kwarg > cartridge plane > uniform). Ambiguity policy is `"fractional"` (splits credit across tie-set entries), `"truth_first"` (uses the first call only), or `"reject"` (drops ambiguous rows into `report.rejected`). Unknown alleles, missing-D-on-VDJ rows, and below-`min_count` drops all surface as structured entries on the build report. Per-segment weights are normalised to sum to 1.0; the `manifest["models"]["allele_usage"]` block surfaces `available` / `nonempty_segments` / `legacy_gene_use_dict_present` / `in_plan_signature=False` (inherited soft gap from `v_allele_weights`). See [`audit-docs/allele_usage_estimation_design.md`](audit-docs/allele_usage_estimation_design.md).

**The other estimators have landed since.** `estimate_trim_distributions`, `estimate_np_length_distributions`, `estimate_np_base_model`, `estimate_p_nucleotide_lengths` and `estimate_genotype_priors` all ship alongside `estimate_allele_usage` and append their own stage entries to the build report; `set_genotype_priors` attaches a population model you already have. Only `estimate_shm_rates` is still deferred. See [`audit-docs/reference_cartridge_authoring_audit.md`](audit-docs/reference_cartridge_authoring_audit.md) for the full design and [Estimate models from data](https://mutejester.github.io/GenAIRR/guides/estimate-cartridge-models.html) for the per-estimator walkthrough.

---

## Key Features

- **Rust simulation kernel** &mdash; persistent IR with full revision history, addressed-trace introspection, `cargo test`-grade unit coverage.
- **Constraint-aware sampling** &mdash; contracts prune candidate distributions at sample time so productive sequences come out of the engine by construction; no retry loops.
- **Strict-mode opt-in** &mdash; surface unsatisfiable plans as `StrictSamplingError` instead of silently relaxing the bundle.
- **Deterministic seeds** &mdash; one seed and one `determinism_version` reproduce every byte of the pool and every entry of the trace; the goldens check that on Linux, macOS and Windows across Python 3.10&ndash;3.13.
- **Full revision history** &mdash; `outcome.revision(i)` exposes the IR after each pass for fine-grained debugging.
- **Addressed trace** &mdash; every random draw is keyed by a hierarchical string (`"np.np1.bases[3]"`) and survives end-to-end into the returned `Outcome`.
- **23 species, 106 configs** &mdash; built-in IMGT + OGRDB reference pickles ship with the wheel.
- **Zero mandatory Python dependencies** &mdash; one wheel, everything in the box.

---

## Optional Extras

```bash
pip install GenAIRR[all]          # numpy, scipy, graphviz, tqdm, fastmcp
pip install GenAIRR[dataconfig]   # numpy + scipy (custom DataConfig analysis)
pip install GenAIRR[viz]          # graphviz
pip install GenAIRR[mcp]          # fastmcp (for the MCP server, see next section)
```

---

## MCP server — drive GenAIRR from an LLM agent

GenAIRR ships an MCP server that exposes 14 tools an LLM agent (Claude, Cursor, etc.) can call to discover configs, simulate repertoires, validate AIRR records, and replay specific seeds — all without writing Python. Install the extra, then point your MCP client at `python -m GenAIRR.mcp_server`:

```bash
pip install GenAIRR[mcp]
```

### Config snippets

**Claude Code** — `.mcp.json` in the project root (or `~/.claude/mcp.json` globally):

```json
{
  "mcpServers": {
    "genairr": {
      "type": "stdio",
      "command": "python",
      "args": ["-m", "GenAIRR.mcp_server"]
    }
  }
}
```

**Claude Desktop** — `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS) or `%APPDATA%\Claude\claude_desktop_config.json` (Windows):

```json
{
  "mcpServers": {
    "genairr": {
      "command": "/path/to/venv/bin/python",
      "args": ["-m", "GenAIRR.mcp_server"]
    }
  }
}
```

**Cursor** — `~/.cursor/mcp.json` or per-project:

```json
{
  "mcpServers": {
    "genairr": {
      "command": "python",
      "args": ["-m", "GenAIRR.mcp_server"]
    }
  }
}
```

Use the full path to the venv's `python` if GenAIRR isn't installed in the system interpreter — the MCP server inherits the launching process's Python environment.

### What you get

After reloading the client, the agent has 14 tools available under the `genairr` namespace:

| Category | Tools |
|----------|-------|
| Discovery | `list_configs`, `config_info`, `list_alleles`, `inspect_allele` |
| Simulation | `simulate_repertoire`, `simulate_preset`, `simulate_allele` |
| Analysis | `validate_records`, `align_to_germline`, `score_allele_calls`, `analyze_mutations`, `classify_regions`, `summarize_dataset` |
| Reproducibility | `replay_seed` |

Every tool returns a uniform `{ok, tool, elapsed_ms, result | error}` envelope; failures carry a stable error-code token (`config_not_found`, `allele_not_found`, `invalid_preset`, `invalid_parameter`, `malformed_record`, `seed_replay_mismatch`) the agent can branch on.

A quick smoke-test prompt to verify the install: *"List the available GenAIRR configs, then simulate 100 productive human heavy-chain sequences with moderate SHM and summarise the V-gene usage."* — the agent should chain `list_configs` → `simulate_repertoire(config="human_igh", n=100, productive_only=true, mutation_model="s5f", mutation_count_min=5, mutation_count_max=15)` and read `v_usage_top` from the result.

---

## Documentation

The full documentation site is at **[mutejester.github.io/GenAIRR](https://mutejester.github.io/GenAIRR/)**. Useful starting points:

- **Getting started** — [Installation](https://mutejester.github.io/GenAIRR/getting-started/installation.html) · [Quick start](https://mutejester.github.io/GenAIRR/getting-started/quick-start.html) · [Your first AIRR record](https://mutejester.github.io/GenAIRR/getting-started/first-airr-record.html) · [Export the results](https://mutejester.github.io/GenAIRR/getting-started/export-results.html)
- **Concepts** — [AIRR Record](https://mutejester.github.io/GenAIRR/concepts/airr-record.html) · [Reference Cartridge](https://mutejester.github.io/GenAIRR/concepts/reference-cartridge.html)
- **Guides** — [The Experiment builder](https://mutejester.github.io/GenAIRR/guides/experiment-builder.html) · [Productive sampling](https://mutejester.github.io/GenAIRR/guides/productive.html) · [Clonal lineage trees](https://mutejester.github.io/GenAIRR/guides/clonal-lineage.html) · [Clonal repertoires](https://mutejester.github.io/GenAIRR/guides/clonal-repertoire.html) · [Genotypes](https://mutejester.github.io/GenAIRR/guides/genotype.html) · [Paired-end reads and FASTQ](https://mutejester.github.io/GenAIRR/guides/paired-end-fastq.html) · [Corruption + sequencing artefacts](https://mutejester.github.io/GenAIRR/guides/corruption-sequencing.html) · [Targeted SHM rates](https://mutejester.github.io/GenAIRR/guides/shm-targeting.html) · [Trace, replay, reproducibility](https://mutejester.github.io/GenAIRR/guides/trace-replay.html)
- **Validation** — [validate_records](https://mutejester.github.io/GenAIRR/validation/validate-records.html) · [Validation overview](https://mutejester.github.io/GenAIRR/validation/index.html)
- **Reference** — [API + Configs + AIRR fields](https://mutejester.github.io/GenAIRR/reference/index.html)

---

## Citing GenAIRR

If GenAIRR is useful in your research, please cite:

> Konstantinovsky T, Peres A, Polak P, Yaari G. An unbiased comparison of immunoglobulin sequence aligners. *Briefings in Bioinformatics*. 2024;25(6):bbae556. [doi:10.1093/bib/bbae556](https://doi.org/10.1093/bib/bbae556)

---

## Contributing

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup and guidelines. If you're adding a pass, a contract, or any code that mutates the simulation IR:

- **Read first**: [audit-docs/engine_architecture.md](audit-docs/engine_architecture.md) — codifies the invariants (contracts constrain support, trace = choices, events = consequences, live-call refresh follows events) and lists the anti-patterns CI catches.
- **Then copy from**: [audit-docs/adding_a_pass.md](audit-docs/adding_a_pass.md) — minimal pass template, three required test patterns with `passes::test_support` helpers, and a crib sheet of which existing pass to model the new one on.

## License

GPL-3.0. See [LICENSE](LICENSE).
