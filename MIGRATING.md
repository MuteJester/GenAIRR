# Migrating

Two migrations live in this file. Start with the one that matches the version
you are coming from.

- [2.3.1 → 2.4.0](#migrating-from-231-to-240) — stored data and clonal
  parameters.
- [1.x → 2.0](#migrating-from-genairr-1x-to-20) — the Experiment DSL rename.

---

# Migrating from 2.3.1 to 2.4.0

2.4.0 does not rename anything you call. What it changes is **what a seed
produces**, **what `n_max` means**, and **which bundled catalogues compile
without an opt-in**. Nothing below is a find/replace; each item is a decision
about data or parameters you already have.

---

## 1. Your stored data is from a different sampling model

`GenAIRR.DETERMINISM_VERSION` is now `4`. Everything GenAIRR released up to and
including 2.3.1 was the unnumbered model that preceded the stamp. A seed that
produced a dataset under 2.3.1 does **not** reproduce it under 2.4.0 — not
approximately, not statistically. Three separate increments moved the draws:
per-pass random substreams and keyed sub-simulation seeds (2),
Beverton-Holt lineage regulation (3), and affinity selection scoring the
protein the cell expresses (4).

**What to do**

- **Keep the data, not the seed.** A dataset you still need is a file of rows;
  keep the file. Do not plan to regenerate it from the seed.
- **Record the version alongside any dataset you publish.** Every record now
  carries a `determinism_version` column and every trace file carries the same
  number. Read the installed engine's number without running anything:

  ```python
  import GenAIRR as ga
  ga.DETERMINISM_VERSION          # 4
  ```

- **Check the column before scoring against someone else's dataset.** A
  dataset with no such column was drawn by the pre-2.4.0 model. Two datasets
  with different numbers are not comparable, and a trace file whose number
  disagrees with the installed engine now fails replay by name instead of
  quietly reproducing something else.
- **Re-pin golden files.** If your test suite asserts simulated bytes for a
  fixed seed, regenerate the fixtures once and commit the new values together
  with the determinism version you generated them at.

---

## 2. `n_max` delivers about 3x more cells

`clonal_lineage(n_max=...)` is now the live population a lineage **settles
at**. Under 2.3.1 it was the pivot of a discrete logistic map that settled at
`n_max * (1 - 1/lambda_base)` and a hard cap truncated any generation that
exceeded it — so a requested 1000 delivered about 333 cells at the default
`lambda_base=1.5`, and reaching capacity killed the lineage outright.

**What to do**

- If you tuned `n_max` to get a population you were happy with, divide it by
  `1 / (1 - 1/lambda_base)` to recover roughly the old plateau. At the default
  `lambda_base=1.5` that is `n_max / 3`.
- Budget memory and runtime for the larger population if you keep your old
  `n_max`.
- Stop reading `lambda_base` as "mean offspring per cell". It is the Poisson
  mean for a founder-baseline cell in an **empty** germinal center. The mean a
  cell actually draws against is that damped by the live population and, under
  selection, scaled by its fitness relative to the founder.

```python
# 2.3.1: a settled population of ~333
exp.clonal_lineage(n_clones=20, n_max=1000, lambda_base=1.5, ...)

# 2.4.0, same population
exp.clonal_lineage(n_clones=20, n_max=333, lambda_base=1.5, ...)

# 2.4.0, if you wanted 1000 all along: no change needed
exp.clonal_lineage(n_clones=20, n_max=1000, lambda_base=1.5, ...)
```

Affinity numbers from before this release are not comparable either, and
neither are the genealogies selection shaped: under 2.3.1 a cell grown without
mutation scored `3.57e-67` against its own protein instead of `1.0`, and a stop
codon scored as a perfect match. A `target_aa` whose length does not match the
receptor's protein is now refused rather than scored over the overlapping
prefix.

---

## 3. `expand_clones()` is deprecated

It still works, and it is still the only surface that delivers an exact
`n_clones x per_clone` count. It now raises a `DeprecationWarning`.

**What to do**

- **Keeping it:** nothing, unless your test suite turns warnings into errors.
  Then filter it:

  ```python
  import warnings
  with warnings.catch_warnings():
      warnings.filterwarnings("ignore", category=DeprecationWarning)
      result = exp.expand_clones(n_clones=10, per_clone=50).run_records(seed=0)
  ```

- **Moving off it:** pick by what the clone means in your benchmark.

  | You want | Use |
  |---|---|
  | A BCR genealogy with ancestral nodes, per-division SHM and a ground-truth tree per clone | `clonal_lineage(n_clones=..., max_generations=..., n_sample=..., rate=...)` |
  | TCR or flat-BCR clone **abundance** — membership plus a size | `clonal_repertoire(n_clones=..., size_distribution=..., max_size=...)` |

  Neither is a drop-in replacement for a fixed `per_clone`. A lineage run
  returns at most `n_sample` records per clone after identical genotypes are
  collapsed; a repertoire run collapses identical copies into one record
  carrying `duplicate_count`.

---

## 4. 80 of the 106 bundled catalogues need an explicit opt-in

Strict compilation — the default — now refuses a bundled catalogue that carries
unresolved anchor or pseudogene/ORF findings. Chains that ran under 2.3.1 can
raise `ValueError` at compile time with the findings listed allele by allele.

Affected: every mouse catalogue, `human_tcrb`, `rhesus_igk`, and 77 others.
Unaffected: the three human OGRDB cartridges (`human_igh`, `human_igk`,
`human_igl`), `HUMAN_IGH_EXTENDED`, and 22 IMGT catalogues whose anchors all
resolved.

**What to do** — one of:

```python
# Sample from the catalogue as it is, acknowledging the findings.
ga.Experiment.on("mouse_igh").allow_curatable_refdata().recombine()

# Or drop the offending alleles first.
ga.Experiment.on("mouse_igh").curate_refdata("functional_anchors_only").recombine()
```

`run_records(..., allow_curatable_refdata=True)` is the per-run equivalent of
the first form.

---

## 5. Counters that now report what you asked for

`mutate(model="uniform")`, `pcr_amplify` and `sequencing_errors` drew their
replacement base from all four bases, including the one already at the
position, so about a quarter of draws wrote back what was there and were
correctly not counted: `pcr_amplify(count=8)` delivered about six errors. The
replacement is now drawn from the three bases the position does not hold, so
`count=N` delivers `N`. **If you calibrated a corruption rate against the old
behaviour, it is now about 32% higher than you intended** — multiply by 0.757
to recover the realized rate you had.

Two further counter changes to know about:

- `n_ambiguous_bases`, not `n_quality_errors`, is what `ambiguous_base_calls()`
  reports on. `n_quality_errors` belongs to `sequencing_errors()`.
- `n_observable_mutations` is the SHM a mutation caller can recover from
  `sequence`. **Score a caller against that field, not `n_mutations`**, which
  counts mutations on untemplated junction bases and positions later lost to
  end loss.

---

## 6. Smaller things that can break a script

- `result.to_dataframe()` needs `pandas`, which is **not** in any extra. The
  docs used to say it shipped under `[all]`; install it yourself.
- `result.validate_records(refdata)` wants the engine-side `RefDataConfig`.
  `exp.compile().refdata` gives you one; passing a `DataConfig` raises
  `TypeError`. `report.summary()` returns a `dict[str, int]`, not a string.
- `mutate(segment_rates={"V": 0, "D": 0, "J": 0, "NP": 0})` is refused with
  `ValueError` instead of silently making the pass a no-op.
- `Experiment.on("x.json")` — any string ending in `.json` — raises `TypeError`
  where it raised `ValueError`. Read a portable reference with
  `GenAIRR.Reference.read(path)`; a path string is never read implicitly.
- Artifact files GenAIRR writes (portable references, family-replay bundles)
  follow your umask instead of being created owner-only `0600`.
- `validate_records` / `validate_families_with_parents` against a reference
  other than the one the run used raises `ValueError` instead of misreporting.
  A genotype run with novel alleles must be validated against
  `compiled.refdata`.

## 7. What you may now want to use

Everything above is about breakage. These are the 2.4.0 additions most
likely to change how you work, in case the migration is the only part of
the release you read:

- **Pin a published germline release, and ship the evidence.**
  `GenAIRR.reference_sources` locks one named OGRDB release and writes a
  dossier a reviewer can verify offline, with no network and no trust in
  you. `ReferenceBuild.bundled("human_igh_ogrdb_10_observed_v1")` needs
  no network at all. See
  [Published germline sources](https://mutejester.github.io/GenAIRR/guides/published-references.html).
- **Real clonal genealogies.** `clonal_lineage()` replaces the star
  topology with a birth-death process and hands back the ground-truth
  tree per clone; `clonal_repertoire()` covers TCR and flat-BCR
  abundance. Both are the intended replacements for `expand_clones()`
  (§3).
- **Per-individual diploid genotypes**, with population priors, cohorts
  and a genotype truth sidecar for scoring inference.
- The full list is the `### Added` section of
  [CHANGELOG.md](CHANGELOG.md).

---

# Migrating from GenAIRR 1.x to 2.0

GenAIRR 2.0 is a breaking release. The Python `Experiment` DSL was reworked
around one principle: a fluent chain should read like a wet-lab protocol,
not a configuration file. Scripts written against 1.x **will not run
unmodified**. This guide is the migration map.

The Rust engine (`engine_rs/`) is unchanged in shape — only the Python
DSL surface moved. Pre-built `_engine.RefDataConfig` / `ContractSet`
objects you may have constructed manually still work the same way.

---

## TL;DR — find/replace

| v1.x | v2.0 |
|---|---|
| `respect=ga.productive()` | `.productive_only()` |
| `.using(v=...)` | `.restrict_alleles(v=...)` |
| `.recombine(trim=True)` | `.recombine()` (default on) |
| `.recombine(trim=False)` | `.recombine().trim(enabled=False)` |
| `mutate(count=...)` | `mutate(rate=...)` *(preferred)* or `mutate(count=...)` |
| `corrupt_pcr(count=...)` | `pcr_amplify(count=...)` |
| `corrupt_quality(count=...)` | `sequencing_errors(count=...)` |
| `corrupt_contaminants(prob=...)` | `contaminate(prob=...)` |
| `corrupt_indels(count=..., insertion_prob=...)` | `polymerase_indels(count=..., insertion_prob=...)` |
| `corrupt_ns(count=...)` | `ambiguous_base_calls(count=...)` |
| `corrupt_reverse_complement(prob=...)` | `random_strand_orientation(prob=...)` |
| `corrupt_5prime_loss(length=...)` | `primer_trim_5prime(length=...)` |
| `corrupt_3prime_loss(length=...)` | `primer_trim_3prime(length=...)` |
| `with_clonal_structure(n_clones=N, size=K)` | `expand_clones(n_clones=N, per_clone=K)` |

The parameter rename inside `expand_clones` is `size` → `per_clone`;
`n_clones` keeps its name.

---

## Side-by-side example

A realistic memory-B-cell chain in **v1.x**:

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine(trim=True)
      .expand_clones(n_clones=50, per_clone=20)
      .mutate(model="s5f", count=(5, 15))
      .corrupt_5prime_loss(length=(0, 8))
      .corrupt_indels(count=(0, 2), insertion_prob=0.5)
      .corrupt_pcr(count=(0, 3))
      .corrupt_ns(count=(0, 2))
      .with_metadata(donor="P1", tissue="peripheral_blood")
      .run_records(seed=42, respect=ga.productive())
)
```

The same chain in **v2.0**:

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .expand_clones(n_clones=50, per_clone=20)
      .mutate(rate=0.05)
      .primer_trim_5prime(length=(0, 8))
      .polymerase_indels(count=(0, 2), insertion_prob=0.5)
      .pcr_amplify(count=(0, 3))
      .ambiguous_base_calls(count=(0, 2))
      .with_metadata(donor="P1", tissue="peripheral_blood")
      .productive_only()
      .run_records(seed=42)
)
```

Try `exp.describe()` (new in v2.0) to read the chain back as biology —
if it reads weird, the chain itself is wrong.

---

## Detailed notes

### `respect=ga.productive()` → `.productive_only()`

In v1 you passed a contract bundle as a runtime kwarg on `.run()` /
`.compile()`. In v2 it's a chainable method that declares the constraint
on the experiment itself. The method is **order-independent** —
`.productive_only()` anywhere in the chain produces the same compiled
simulator.

```python
# v1
exp.run_records(n=100, seed=0, respect=ga.productive())

# v2
exp.productive_only().run_records(n=100, seed=0)
```

The `respect=` kwarg is removed from `compile()`, `run()`, `run_records()`,
`stream()`, and `stream_records()`. Pass `strict=True` if you want
runtime-impossible contracts to raise instead of fall back to permissive
sampling (unchanged from v1).

### `.using(...)` → `.restrict_alleles(...)`

Same semantics (narrows allele sampling to a named subset, sampling
stays uniform within the subset), but the new name actually describes
what's happening. v1's "using" was ambiguous — using what for what?

### `.recombine(trim=…)` → `.trim(enabled=…)`

Trim is now its own chainable method, called *after* `.recombine()`.
Default-on stays the default-on; you only need to call `.trim()` if you
want to disable it or supply custom distributions.

```python
# v1
exp.recombine(trim=False).mutate(...)

# v2
exp.recombine().trim(enabled=False).mutate(...)
```

`.trim()` is **position-guarded**: calling it before `.recombine()` or
after a mutation / corruption step raises `ValueError`. This catches
the confusion where v1 users sometimes thought the kwarg controlled
post-recombination trimming.

For custom distributions:

```python
exp.recombine().trim(
    v_3=[(0, 1.0), (1, 1.0), (2, 1.0)],
    j_5=[(0, 1.0), (1, 1.0)],
)
```

Omitted segments keep their empirical defaults.

### `mutate(count=…)` → `mutate(rate=…)` (preferred)

`mutate(rate=0.03)` is the canonical v2 form. Per record, the count is
drawn from `Poisson(rate × pool_len)` — so a 5% rate on a 350-bp record
expects ~17 mutations, but the actual count varies stochastically and
*scales with each record's actual sequence length*. This matches how
immunologists report SHM in the literature.

`mutate(count=…)` is still accepted for benchmark scripts that need a
deterministic count independent of record length. You **must** pass
exactly one of `rate` or `count` — both, or neither, raises
`ValueError`.

### `corrupt_*` → wet-lab-stage names

The v1 `corrupt_*` umbrella lumped PCR amplification, sequencing
readout, library prep, and contamination under one engineering verb.
v2 splits them by wet-lab stage. See the table at the top of this guide.

A note on `primer_trim_5prime` / `primer_trim_3prime` vs `.trim()`:
**these are biologically distinct operations** at different stages.

- `.trim()` configures **exonuclease trim** during V(D)J recombination
  (RAG-mediated, before NP insertion). Only one per experiment.
- `primer_trim_5prime` / `primer_trim_3prime` model **read-end loss**
  during library prep / sequencing. They drop bases off the assembled
  sequence and can be called alongside other library-stage steps.

### `with_clonal_structure` → `expand_clones`

The new name is the verb of the biological process (B-cell clonal
expansion). Parameter rename inside: `size=20` becomes
`per_clone=20`. (`n_clones` stays as-is.) Both must be positive
ints; calling `expand_clones()` twice still raises.

---

## Things that didn't change

These continue to work identically:

- `Experiment.on(config_name | DataConfig | RefDataConfig)`.
- `.recombine(np1_lengths=..., np2_lengths=..., v_allele_weights=...)`.
- `.with_metadata(**fields)`.
- `.run(n=..., seed=..., strict=...)` returning `Outcome` objects.
- `.run_records(n=..., seed=..., strict=..., expose_provenance=...)`
  returning `SimulationResult`.
- Engine-direct `ga._engine.PassPlan` + `ga._engine.run(...)`.
- The S5F kernels and their string IDs (`hh_s5f`, `hkl_s5f`, etc.).
- `ga.productive()` still returns a `ContractSet` — it's just no
  longer the way you wire it into an `Experiment` (use
  `.productive_only()` instead).

---

## Pre-flight checklist

Before opening a PR with your migrated scripts:

1. `grep -rn 'corrupt_\|respect=\|\.using\|with_clonal_structure\|trim=' your_scripts/` returns nothing relevant.
2. Each migrated experiment prints a clean `exp.describe()` — read it back as a wet-lab protocol.
3. Your test suite passes against the new method shapes.

If something still feels wrong, file an issue with a v1 / v2 chain pair
and the expected behavior — the v2 surface is designed to read
naturally, and discrepancies with the wet-lab mental model are bugs we
want to know about.
