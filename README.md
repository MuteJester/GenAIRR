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

rec["sequence"]             # 'gaggtgcagctggtggagtctgggggaggc...' (nucleotide)
rec["sequence_aa"]          # 'EVQLVESGGGLVQPGGSLRLSCSAS...'      (translated)
rec["locus"]                # 'IGH'
rec["v_call"]               # 'IGHVF10-G38*04'   (comma-separated if the call ties)
rec["d_call"]               # 'IGHD2-15*01'
rec["j_call"]               # 'IGHJ2*01'
rec["junction_aa"]          # 'CVKDDGNRGYCSGGSCYGRCCALDYWYFDLW'
rec["productive"]           # True
rec["v_identity"]           # 1.0  (matches/total over the V segment)
rec["n_mutations"]          # 0

# Export in any of the standard formats. TSV/FASTA/FASTQ are dependency-free;
# to_dataframe() needs pandas (pip install GenAIRR[all]).
result.to_tsv("repertoire.tsv")        # AIRR-spec TSV (50+ columns)
result.to_fasta("sequences.fasta")     # FASTA with v_call/j_call in the headers
result.to_fastq("sequences.fastq")     # FASTQ with illumina-shaped quality scores
df = result.to_dataframe()             # one row per record, AIRR columns
```

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

Headers use the AIRR record's `sequence_id` with the universally-portable `/1` and `/2` suffix (`@seq0/1` and `@seq0/2`) — no `|`-pipe metadata that some aligners (STAR < 2.7) split on. R1 / R2 sequences are written verbatim from the AIRR `r1_sequence` / `r2_sequence` fields (R2 is already reverse-complemented at projection time; the writer does NOT apply a second flip). Default `overwrite=False` refuses to clobber existing files; pass `overwrite=True` to replace. Quality strings use the same pluggable models as `to_fastq` — `quality="illumina"` (default trapezoid shape, applied per read) or `quality="constant"` with `q=...`. See [`docs/fastq_export_design.md`](audit-docs/fastq_export_design.md).

`Experiment.on(...)` accepts **a config-name string** (e.g. `"human_igh"`, `"mouse_tcrb"`), **a `DataConfig`** loaded from the bundled species pickles, or **a `RefDataConfig`** for [custom reference data](#custom-reference-data). `.productive_only()` is the constraint-aware bundle — covered in the next section. Drop it to allow non-productive sequences (~30% of records will then have stop codons in the junction).

> See the full walkthrough in the docs: [Quick Start](https://mutejester.github.io/GenAIRR/lesson-1.html) · [Interpreting Results](https://mutejester.github.io/GenAIRR/concept-airr-record.html)

---

## A realistic pipeline — everything in one place

The Experiment DSL is a fluent builder. Each step appends to the pipeline; the same `Experiment` is returned so calls chain. The example below uses every major feature GenAIRR offers — recombination, clonal expansion, per-descendant somatic hypermutation, primer-trimming, structural indels, PCR errors, N-base injection, custom metadata, and the productive constraint:

```python
import GenAIRR as ga

result = (
    ga.Experiment.on("human_igh")
      # 1. V(D)J recombination — sample alleles, trim, fill NP1/NP2, assemble.
      .recombine()
      # 2. Clonal structure — 50 lineages × 20 sister sequences each.
      #    Passes BEFORE this point apply to the parent rearrangement;
      #    passes AFTER apply per-descendant. So each clone shares the
      #    same V(D)J recombination but accumulates its own SHM + errors.
      #    NOTE: expand_clones is the legacy fixed-size *star* model. For
      #    real clonal trees and repertoires, see "Clonal lineages &
      #    repertoires" below (clonal_lineage / clonal_repertoire).
      .expand_clones(n_clones=50, per_clone=20)
      # 3. Somatic hypermutation per descendant — S5F context-dependent
      #    model at 5% per-base rate (matches memory-B-cell SHM).
      .mutate(rate=0.05)
      # 4. Sequencing artefacts per descendant: primer trimming, structural
      #    indels, PCR substitution errors, quality-driven N injection.
      .primer_trim_5prime(length=(0, 8))
      .primer_trim_3prime(length=(0, 4))
      .polymerase_indels(count=(0, 2), insertion_prob=0.5)
      .pcr_amplify(count=(0, 3))
      .ambiguous_base_calls(count=(0, 2))
      # 5. Stamp arbitrary metadata onto every record.
      .with_metadata(experiment_id="exp001", tissue="peripheral_blood")
      # Constraint-aware sampling: the productive() bundle is enforced at
      # rearrangement + SHM time. Corruption passes can still introduce
      # stop codons / frameshifts post-hoc, so expect ~70% productive
      # when aggressive corruption is in the chain — that mirrors real
      # wet-lab data, where a productive B-cell can sequence as a
      # non-productive read because of an indel during library prep.
      .productive_only().run_records(seed=42)
)

len(result)                                  # 1000  (= n_clones × size)
sum(1 for r in result if r["productive"])    # 697   (~70% under this corruption load)

# Same clone, different descendants — same V(D)J recombination,
# independent SHM + errors:
result[0]["clone_id"], result[1]["clone_id"]              # (0, 0)
result[0]["v_call"],   result[1]["v_call"]                # both 'IGHVF10-G38*04'
result[0]["n_mutations"], result[1]["n_mutations"]        # (13, 15) — independent SHM
result[0]["n_pcr_errors"], result[1]["n_pcr_errors"]      # (1, 1)   — independent errors

# Custom metadata propagated:
result[0]["experiment_id"], result[0]["tissue"]           # ('exp001', 'peripheral_blood')

result.to_tsv("repertoire.tsv")
```

## Clonal lineages & repertoires

GenAIRR models clonal structure three ways. Pick by what you need:

| Method | Model | Use for |
|---|---|---|
| **`clonal_lineage`** | BCR affinity-maturation **trees** — generation-synchronous birth–death, per-division S5F SHM, optional affinity selection, sample + genotype-collapse | B-cell lineage trees with **ground truth** (Newick/FASTA per clone), benchmarking lineage/clone inference, ML training data |
| **`clonal_repertoire`** | Non-tree **abundance** repertoires — heavy-tailed clone sizes (power-law/lognormal) + unexpanded-singleton fraction, reads through library-prep, collapsed to `duplicate_count` | **TCR** repertoires (no SHM) and flat-BCR abundance; clone-calling benchmarks |
| `expand_clones` *(deprecated)* | Fixed-size **star**: one founder + `per_clone` independent descendants | legacy; kept working, but prefer the two above |

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
[Clonal lineage trees](site_docs/guides/clonal-lineage.md) ·
[Clonal repertoires (TCR & abundance)](site_docs/guides/clonal-repertoire.md).

Other feature flags worth knowing:

| Step | What it does |
|-----|-----|
| `.contaminate(prob=0.02)` | Replace ~2% of records with unrelated background sequences. |
| `.sequencing_errors(count=(0, 5))` | Lowercase 0–5 bases per sequence to mark sequencer-low-quality positions. |
| `.random_strand_orientation(prob=0.5)` | Flip ~50% of records to the reverse strand (with the `rev_comp` flag set). |
| `.restrict_alleles(v=[...], d=[...], j=[...])` | Restrict allele sampling to a specific subset — useful for benchmarking against a known repertoire. |
| `.mutate(model="uniform", rate=0.03)` | Use a uniform-rate mutation model instead of S5F. |
| `.mutate(model="s5f", rate=0.03, segment_rates={"V": 1.0, "D": 0.2, "J": 0.5, "NP": 0.0})` | Restrict SHM targeting by biological region class. Buckets are `"V"`, `"D"`, `"J"`, `"NP"` (the NP entry covers both NP1 and NP2). Omitted buckets default to `1.0`; `0.0` disables a region entirely so sites in it drop out of proposal support before contract admissibility. Default (no kwarg) is byte-identical to the pre-slice engine. Composes with `productive_only()`; replay reproduces exactly when the same rate vector is used. The rate vector itself is **not** part of cartridge identity / Rust `content_hash` — it's a per-experiment parameter (audit's documented v1 boundary). |
| `.mutate(model="s5f", rate=0.03, segment_rates={"V": 1.0, "NP": 0.0}, v_subregion_rates={"CDR": 2.0, "FWR": 0.5})` | Refine V-segment SHM targeting by IMGT subregion on top of `segment_rates`. Accepted keys are the five canonical labels `"FWR1"` / `"CDR1"` / `"FWR2"` / `"CDR2"` / `"FWR3"` plus the two aliases `"FWR"` (expands to FWR1 / FWR2 / FWR3) and `"CDR"` (expands to CDR1 / CDR2). **Alias expansion runs first, then explicit labels override** — `{"FWR": 0.5, "FWR2": 2.0}` resolves to `FWR1=0.5, FWR2=2.0, FWR3=0.5, CDR1=1.0, CDR2=1.0`. The V-site weight becomes `base × segment_rate(V) × v_subregion_rate(label)`. Non-V sites and V alleles without subregion annotations receive identity factor `1.0`, so the kwarg composes cleanly with mixed cartridges and with `productive_only()`. Default (no kwarg, `{}`, or explicit all-ones) is byte-identical to the pre-slice engine. Requires the cartridge to carry V-subregion annotations — the bundled human OGRDB cartridges do; calling `v_subregion_rates` against an unannotated cartridge raises `ValueError` at builder time. Like `segment_rates`, the rate vector is part of the plan signature (replay against a different vector fails before consuming choices) but **not** part of `refdata_content_hash`. See [`docs/v_subregion_shm_rate_design.md`](audit-docs/v_subregion_shm_rate_design.md). |
| **`n_mutations` / `n_v_mutations` / `n_d_mutations` / `n_j_mutations` / `n_np_mutations`** (AIRR fields) | `n_mutations` is the **global biological SHM** counter (uniform + S5F substitutions only, sourced from `Simulation.mutation_count`). The four per-segment fields **partition** it by carried event segment — `n_v + n_d + n_j + n_np == n_mutations` by construction, enforced by the validator's `MutationCountSumMismatch` cross-check. NP1 and NP2 events roll into `n_np_mutations` (matches the `segment_rates` DSL grouping). **PCR / quality / N-corruption / receptor revision / D inversion / contaminant are intentionally excluded** — observation-stage and recombination-stage changes do not count as SHM. PCR + quality artefacts surface separately as `n_pcr_errors` / `n_quality_errors`. See [`docs/mutation_provenance_audit.md`](audit-docs/mutation_provenance_audit.md). |
| **`n_fwr1_mutations` / `n_cdr1_mutations` / `n_fwr2_mutations` / `n_cdr2_mutations` / `n_fwr3_mutations` / `n_v_unannotated_mutations`** (AIRR fields) | **Per-V-subregion SHM counters** — a partition of **`n_v_mutations`** (not of the global `n_mutations`). The five canonical IMGT labels bucket V SHM events by the assigned V allele's `FWR1` / `CDR1` / `FWR2` / `CDR2` / `FWR3` interval; `n_v_unannotated_mutations` catches everything else under V. Three "unannotated" cases: (1) the user's cartridge has at least one V allele without subregion annotations (legacy / hand-authored), (2) the SHM event lands in the **V-side CDR3 contribution stretch** between `FWR3.end` and `len(V_allele.seq)` — the five canonical labels deliberately stop at FWR3, so junctional V-tail SHM goes here, (3) an SHM pass runs after an indel pass and mutates an inserted V base (non-canonical pass order). On bundled human OGRDB cartridges the unannotated bucket is a small but non-zero baseline (case 2 — V-tail / CDR3 contribution). Partition invariant: `n_fwr1 + n_cdr1 + n_fwr2 + n_cdr2 + n_fwr3 + n_v_unannotated == n_v_mutations` on every record, enforced by the validator's six `N<Region>MutationsMismatch` checks plus the `VSubregionMutationCountSumMismatch` cross-field invariant. Two-bucket aggregates (`n_cdr_mutations` / `n_fwr_mutations`) are deliberately NOT exposed — derive them downstream as `df["n_cdr_mutations"] = df["n_cdr1_mutations"] + df["n_cdr2_mutations"]`. See [`docs/v_subregion_mutation_counters_audit.md`](audit-docs/v_subregion_mutation_counters_audit.md). |
| `.invert_d(prob=0.05)` | Heavy-chain only. With probability `prob`, sample the D allele in reverse-complement orientation (V(D)J inversion event). Records the decision under `sample_allele.d.inverted`; the assembled D bytes are the WC complement of the original allele; AIRR records carry `d_inverted: bool` for provenance. |
| `.receptor_revision(prob=0.05)` | Heavy-chain only. With probability `prob`, replace V after initial recombination with a different germline allele (B-cell receptor revision). Records `receptor_revision.applied` (+ replacement allele id / 3' trim on `True`); the V slice in the assembled pool is rewritten; AIRR records carry `receptor_revision_applied: bool` and `original_v_call: str` (the pre-revision V; empty when no revision happened). `v_call` continues to report the post-revision identity. |
| `.paired_end(r1_length=150, insert_size=300)` | Both VDJ and VJ. Add Illumina-style paired-end read-layout projection to every record: AIRR fields `r1_sequence` (forward window) + `r2_sequence` (reverse-complement of the 3' window) + `r1_start/end` / `r2_start/end` / `insert_size` / `read_layout="paired_end"`. Optional `r2_length` (defaults to `r1_length`). Each length accepts `int` / `(low, high)` / empirical `[(value, weight)]`. Trace records `paired_end.r1_length` / `.r2_length` / `.insert_size`; AIRR builder reads them back to populate the windows. Projection-only — no IR mutation, no live-call invalidation, runs after end-loss + rev-comp. |
| `compile()` then `compiled.run_records(...)` | Compile the plan once, reuse it across many batches — see [Compile once](#compile-once-run-many-times). |

---

## Constraint-aware sampling

GenAIRR's signature feature is **constraint-aware sampling**: contracts that prune the candidate distribution at sample time, not retries after the fact. The canonical bundle is `productive()` (in-frame junction + no stop codons + V/J anchors preserved):

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

> Docs: [Productive sequences](https://mutejester.github.io/GenAIRR/guide-productive.html)

### Strict vs permissive mode

By default, if a contract can't admit any candidate at a sampling step the runtime falls back to unconstrained sampling and the run continues. Pass `strict=True` to surface the failure as an exception instead — useful for catching unsatisfiable plans early during development:

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

```python
import GenAIRR as ga

# Same seed → byte-identical records across runs and platforms.
a = ga.Experiment.on("human_igh").recombine().run_records(n=100, seed=42)
b = ga.Experiment.on("human_igh").recombine().run_records(n=100, seed=42)
assert a[0]["sequence"] == b[0]["sequence"]

# `n` runs use seeds [seed, seed+1, ..., seed+n-1] so consecutive
# batches stitch together by offsetting the starting seed.
batch_a = ga.Experiment.on("human_igh").recombine().run_records(n=100, seed=0)
batch_b = ga.Experiment.on("human_igh").recombine().run_records(n=100, seed=100)
# batch_a[50] is byte-equal to a one-off run at seed=50.
```

> Docs: [Reproducibility](https://mutejester.github.io/GenAIRR/guide-reproduce.html)

---

## Validation posture

GenAIRR's release readiness rests on five guardrails, each backed by
an audit doc and a golden-test file. Every PR runs the full suite;
every release additionally runs the golden trace compatibility
fixtures so on-disk formats stay byte-stable.

- **Productive validity** — junction-frame, no-stop, anchor-preserved bundle, with strict / permissive failure modes pinned per pass.
- **Provenance correctness** — indels, end-loss, allele-call ambiguity, junction fields, and per-pass event ledger each have isolated-scenario tests.
- **Reproducibility** — every audit slice includes trace-replay round-trips that reproduce sequence, AIRR coordinates, and per-pass event counts.
- **Distribution invariants** — Monte-Carlo tests (±5σ) prove the constrained samplers draw from `natural_weight × admissibility`, with explicit negative controls against renormalization bugs.
- **Performance budgets** — wall-time regression guards on seven representative workloads catch ~10× slowdowns before they ship.
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
[`docs/validation_matrix.md`](audit-docs/validation_matrix.md). Use it to
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
| `result.to_tsv(path, *, airr_strict=False)` | — | AIRR-format TSV. `airr_strict=True` converts coordinates to 1-based-inclusive per spec. |
| `result.to_csv(path, *, airr_strict=False)` | — | Comma-separated. Same options as `to_tsv`. |
| `result.to_fasta(path, *, prefix="seq")` | — | FASTA. Headers include `v_call` and `j_call`. |
| `result.to_fastq(path, *, quality="illumina", **kw)` | — | FASTQ. Quality models: `"illumina"` (smoothed trapezoid) or `"constant"`. |
| `result.to_dataframe(*, airr_strict=False)` | `pandas.DataFrame` | One row per record. Requires pandas (`pip install GenAIRR[all]`). |
| `result.outcomes` | `list[Outcome] \| None` | The underlying `Outcome` objects, for advanced introspection (see below). |

Each record dict has 50+ AIRR fields. The most commonly used:

| Field | Example value | Description |
|-----|-----|-----|
| `sequence` | `'gaggtgcagctggtg…'` | Assembled nucleotide sequence (uppercase + lowercase corruption markers). |
| `sequence_aa` | `'EVQLVESGGG…'` | Codon-rail translation. Stops emit `*`, ambiguous codons emit `X`. |
| `locus` | `'IGH'` | Locus code derived from `v_call` / `j_call`. |
| `v_call` / `d_call` / `j_call` | `'IGHV3-23*01'` | Gene calls. Comma-separated tie set when the evidence walker can't disambiguate. |
| `junction` / `junction_aa` | `'TGC…GAC'` / `'CAR…D'` | Junction nucleotide + AA. AA includes the V Cys (anchor) through J W/F+3. |
| `productive` | `True` / `False` / `None` | In-frame junction AND no stop codons AND anchors preserved. `None` when undefined (e.g. junction not present). |
| `v_identity` / `d_identity` / `j_identity` | `0.987` | Match rate over each segment's CIGAR M/D ops. |
| `v_cigar` / `d_cigar` / `j_cigar` | `'17D279M'` | CIGAR strings. Only M/I/D ops are emitted — no soft-clips. |
| `n_mutations` / `n_pcr_errors` / `n_quality_errors` / `n_indels` | `4` / `0` / `2` / `1` | Per-record error counts from the trace. |

The full schema (plus the `*_sequence_start/end`, `*_alignment_start/end`, `*_germline_start/end` coordinate fields, `vj_in_frame`, `stop_codon`, `rev_comp`, and others) is documented at [Interpreting Results](https://mutejester.github.io/GenAIRR/concept-airr-record.html).

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

> Docs: [Simulation Pipeline](https://mutejester.github.io/GenAIRR/concept-pipeline.html) · [Ground-truth Contracts](https://mutejester.github.io/GenAIRR/concept-contracts.html) · [Interpreting Results](https://mutejester.github.io/GenAIRR/concept-airr-record.html)

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

ga.Experiment.on("mouse_igh").recombine().run_records(n=500)
ga.Experiment.on("rabbit_tcrb").recombine().run_records(n=500)
ga.Experiment.on("rhesus_igk").recombine().run_records(n=500)
```

> Docs: [Learn GenAIRR](https://mutejester.github.io/GenAIRR/learn.html) · [Reference (Configs)](https://mutejester.github.io/GenAIRR/reference.html)

---

## Custom reference data

For non-builtin alleles (custom IMGT pulls, in-house references, etc.) you can build a `RefDataConfig` directly and pass it to `Experiment.on(...)`:

```python
import GenAIRR as ga

cfg = ga.RefDataConfig.vj()
cfg.add_v_allele("v_custom*01", "v_custom", b"GAAGTACAGCTGGTGCAG...", anchor=288)
cfg.add_v_allele("v_custom*02", "v_custom", b"GAAGTACAGCTAGTGCAG...", anchor=288)
cfg.add_j_allele("j_custom*01", "j_custom", b"TGGGGCCAAGGG...",       anchor=10)

result = ga.Experiment.on(cfg).recombine().run_records(n=100, seed=42)
```

`RefDataConfig.vdj()` builds a heavy-chain-shaped refdata (with a D pool); `add_d_allele(...)` populates it. Anchors are 0-based offsets of the V Cys / J W or F codon's first base, used to keep the junction frame-aligned during recombination.

For non-human or otherwise non-standard references — custom J anchor amino acids, extended sequence alphabets, hand-authored NP-length or trim distributions, pseudogene-aware allele filtering — configure a **reference cartridge** instead of relying on the built-in locus defaults. The cartridge has four typed planes (identity, catalogue, rules, empirical models) plus an explicit curation policy. See [`docs/reference_cartridge.md`](audit-docs/reference_cartridge.md) for the model, the authoring surface (`ReferenceRulesSpec`, `ReferenceEmpiricalModels`, `Experiment.curate_refdata`), and end-to-end examples.

Cartridges can also carry **V-region substructure annotations** — per-V-allele IMGT `FWR1` / `CDR1` / `FWR2` / `CDR2` / `FWR3` intervals derived from the IMGT-gapped reference sequence. The bundled OGRDB cartridges populate them at load time for all human IGH / IGK / IGL V alleles; user cartridges can override the dict explicitly. Subregion intervals are inspectable on the bridged `RefDataConfig`, surfaced in `DataConfig.cartridge_manifest()["models"]["shm"]["v_subregion_support"]`, and folded into `refdata_content_hash` so two cartridges that differ only in subregion boundaries hash differently. **They also drive SHM targeting** via `Experiment.mutate(v_subregion_rates={…})` — see the targeted-SHM example in the pipeline-steps table above. Per-region AIRR mutation counters (`n_cdr1_mutations` etc.) are still deferred — that's a separate future slice. See [`docs/v_region_substructure_audit.md`](audit-docs/v_region_substructure_audit.md) and [`docs/v_subregion_shm_rate_design.md`](audit-docs/v_subregion_shm_rate_design.md) for the audits and per-slice scope.

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

**Legacy `NP_transitions` / `NP_first_bases` are NOT auto-lifted yet.** The bundled cartridges still carry those dicts as documented orphan fields, and the manifest reports `legacy_fallback=False`. Authors who want Markov today populate `ReferenceEmpiricalModels.np_bases` explicitly; an opt-in auto-lift slice is a separate cartridge-migration decision. See [`docs/np_markov_base_generator_design.md`](audit-docs/np_markov_base_generator_design.md) and [`docs/junction_n_addition_audit.md`](audit-docs/junction_n_addition_audit.md) for the full design.

### Templated P-nucleotide (palindromic) additions

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

**Legacy `p_nucleotide_length_probs` is NOT auto-lifted yet.** The bundled cartridges still carry the orphan dict; the manifest reports `legacy_fallback=False`. v1 ships *lengths-only* — per-base P strings (`p_v_3`, ...) and the aggregate `n_p_nucleotides` are deferred. See [`docs/p_nucleotide_design.md`](audit-docs/p_nucleotide_design.md) for the full design.

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

Recombination on the resulting cartridge **uses the estimated weights by default** — `ga.Experiment.on(cfg).recombine()` now samples alleles according to the typed plane unless an explicit `recombine(v_allele_weights=..., d_allele_weights=..., j_allele_weights=...)` kwarg is passed (kwarg > cartridge plane > uniform). Ambiguity policy is `"fractional"` (splits credit across tie-set entries), `"truth_first"` (uses the first call only), or `"reject"` (drops ambiguous rows into `report.rejected`). Unknown alleles, missing-D-on-VDJ rows, and below-`min_count` drops all surface as structured entries on the build report. Per-segment weights are normalised to sum to 1.0; the `manifest["models"]["allele_usage"]` block surfaces `available` / `nonempty_segments` / `legacy_gene_use_dict_present` / `in_plan_signature=False` (inherited soft gap from `v_allele_weights`). See [`docs/allele_usage_estimation_design.md`](audit-docs/allele_usage_estimation_design.md).

**Remaining estimators stay deferred.** `estimate_trim_distributions` / `estimate_np_length_distributions` / `estimate_np_base_model` / `estimate_p_nucleotide_lengths` / `estimate_shm_rates` are scoped to follow-up slices — each will append new stage entries to the build report without changing the facade. See [`docs/reference_cartridge_authoring_audit.md`](audit-docs/reference_cartridge_authoring_audit.md) for the full design.

---

## Key Features

- **Rust simulation kernel** &mdash; persistent IR with full revision history, addressed-trace introspection, `cargo test`-grade unit coverage.
- **Constraint-aware sampling** &mdash; contracts prune candidate distributions at sample time so productive sequences come out of the engine by construction; no retry loops.
- **Strict-mode opt-in** &mdash; surface unsatisfiable plans as `StrictSamplingError` instead of silently relaxing the bundle.
- **Deterministic seeds** &mdash; same seed reproduces every byte of the pool and every entry of the trace, across runs and platforms.
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

- **Learn** — [Lesson 1: V(D)J recombination](https://mutejester.github.io/GenAIRR/lesson-1.html) · [Lesson 2: Pipeline scrubber](https://mutejester.github.io/GenAIRR/lesson-2.html) · [Lesson 3: S5F SHM](https://mutejester.github.io/GenAIRR/lesson-3.html) · [Lesson 4: Sequencing artifacts](https://mutejester.github.io/GenAIRR/lesson-4.html) · [Lesson 5: Ground-truth payoff](https://mutejester.github.io/GenAIRR/lesson-5.html)
- **Concepts** — [Simulation Pipeline](https://mutejester.github.io/GenAIRR/concept-pipeline.html) · [Persistent IR](https://mutejester.github.io/GenAIRR/concept-persistent-ir.html) · [Contracts](https://mutejester.github.io/GenAIRR/concept-contracts.html) · [AIRR Record](https://mutejester.github.io/GenAIRR/concept-airr-record.html) · [Live Calls](https://mutejester.github.io/GenAIRR/concept-live-call.html)
- **Guides** — [Build a config](https://mutejester.github.io/GenAIRR/guide-build-config.html) · [Productive sampling](https://mutejester.github.io/GenAIRR/guide-productive.html) · [Clonal families](https://mutejester.github.io/GenAIRR/guide-clonal-families.html) · [Export](https://mutejester.github.io/GenAIRR/guide-export.html) · [Replay](https://mutejester.github.io/GenAIRR/guide-replay.html) · [Reproduce a dataset](https://mutejester.github.io/GenAIRR/guide-reproduce.html) · [Compare SHM models](https://mutejester.github.io/GenAIRR/guide-compare-shm.html) · [Tune corruption](https://mutejester.github.io/GenAIRR/guide-tune-corruption.html)
- **Reference** — [API + Configs + AIRR fields](https://mutejester.github.io/GenAIRR/reference.html)

---

## Citing GenAIRR

If GenAIRR is useful in your research, please cite:

> Konstantinovsky T, Peres A, Polak P, Yaari G. An unbiased comparison of immunoglobulin sequence aligners. *Briefings in Bioinformatics*. 2024;25(6):bbae556. [doi:10.1093/bib/bbae556](https://doi.org/10.1093/bib/bbae556)

---

## Contributing

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup and guidelines. If you're adding a pass, a contract, or any code that mutates the simulation IR:

- **Read first**: [docs/engine_architecture.md](audit-docs/engine_architecture.md) — codifies the invariants (contracts constrain support, trace = choices, events = consequences, live-call refresh follows events) and lists the anti-patterns CI catches.
- **Then copy from**: [docs/adding_a_pass.md](audit-docs/adding_a_pass.md) — minimal pass template, three required test patterns with `passes::test_support` helpers, and a crib sheet of which existing pass to model the new one on.

## License

GPL-3.0. See [LICENSE](LICENSE).
