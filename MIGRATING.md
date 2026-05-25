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
| `corrupt_indels(count=..., insertion_prob=...)` | `library_indels(count=..., insertion_prob=...)` |
| `corrupt_ns(count=...)` | `mask_low_quality(count=...)` |
| `corrupt_reverse_complement(prob=...)` | `random_strand_orientation(prob=...)` |
| `corrupt_5prime_loss(length=...)` | `primer_trim_5prime(length=...)` |
| `corrupt_3prime_loss(length=...)` | `primer_trim_3prime(length=...)` |
| `with_clonal_structure(n_clones=N, size=K)` | `expand_clones(n=N, per_clone=K)` |

Parameter renames are inside `expand_clones`: `n_clones` → `n`, `size` → `per_clone`.

---

## Side-by-side example

A realistic memory-B-cell chain in **v1.x**:

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine(trim=True)
      .with_clonal_structure(n_clones=50, size=20)
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
      .expand_clones(n=50, per_clone=20)
      .mutate(rate=0.05)
      .primer_trim_5prime(length=(0, 8))
      .library_indels(count=(0, 2), insertion_prob=0.5)
      .pcr_amplify(count=(0, 3))
      .mask_low_quality(count=(0, 2))
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
expansion). Parameter renames inside: `n_clones=10, size=20` becomes
`n=10, per_clone=20`. Both must be positive ints; calling
`expand_clones()` twice still raises.

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
