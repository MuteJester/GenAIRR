---
title: API Reference
sidebar_label: Overview
---

# API Reference

This page is a comprehensive catalog of the GenAIRR Python API.

---

## The `Experiment` Class

The `Experiment` class is a fluent builder used to construct simulation pipelines.

### Initialization

#### `Experiment.on(config)`
Starts a simulation on a specific reference configuration.
*   **`config`** (`str`, `DataConfig`, or `RefDataConfig`): The name of a built-in config (e.g., `"human_igh"`), a loaded `DataConfig` object, or a native `RefDataConfig`.

### Biological Recombination

#### `.recombine(*, np1_lengths, np2_lengths, trim, v_allele_weights, d_allele_weights, j_allele_weights)`
Adds a V(D)J or VJ recombination pass.
*   **`np1_lengths` / `np2_lengths`** (`Iterable[Tuple[int, float]]`): Custom (length, weight) distribution for NP regions. Defaults to empirical data.
*   **`trim`** (`bool`): Whether to apply exonuclease trimming. Default `True`.
*   **`v/d/j_allele_weights`** (`Dict[str, float]`): Bias sampling by providing specific weights for named alleles.

#### `.using(*, v, d, j)`
Restricts the allele pool for the *next* recombination step.
*   **`v`, `d`, `j`** (`str` or `List[str]`): A single allele name or a list of names. Lists are sampled uniformly.

### Somatic Hypermutation

#### `.mutate(*, model, count, s5f_model)`
Adds a somatic hypermutation pass. (B-cells only).
*   **`model`** (`str`): `"s5f"` (default) or `"uniform"`.
*   **`count`**: The number of mutations to apply (see [Distributions](#distributions) below).
*   **`s5f_model`** (`str`): The specific S5F kernel to use. Options: `"hh_s5f"` (default), `"hh_s5f_60"`, `"hh_s5f_opposite"`, `"hkl_s5f"`.

### Technical Artifacts (Corruption)

All corruption methods accept a `count` or `length` [Distribution](#distributions).

| Method | Argument | Description |
|--------|----------|-------------|
| `.corrupt_pcr()` | `count` | Random base substitutions mimicking PCR errors. |
| `.corrupt_quality()` | `count` | Position-dependent sequencing quality substitutions. |
| `.corrupt_5prime_loss()` | `length` | Deletes bases from the 5' end. |
| `.corrupt_3prime_loss()` | `length` | Deletes bases from the 3' end. |
| `.corrupt_indels()` | `count` | Random insertions and deletions. |
| `.corrupt_ns()` | `count` | Replaces bases with `N`. |
| `.corrupt_rev_comp()` | `prob` | Probability (0.0 - 1.0) of reverse-complementing the read. |
| `.corrupt_contaminants()` | `prob` | Probability of replacing the entire sequence with noise. |

### Advanced Structure

#### `.with_clonal_structure(*, n_clones, size)`
Forks the pipeline into clonal families.
*   **`n_clones`** (`int`): Number of independent recombination events.
*   **`size`** (`int`): Number of descendant reads per clone.

#### `.with_metadata(**fields)`
Injects custom fields into every output record.
*   **`**fields`**: Keyword arguments (e.g., `sample_id="P1"`) to be included in the AIRR record.

### Execution

#### `.run(n, seed, respect, strict)`
Compiles and runs the simulation. Returns a list of native `Outcome` objects.
*   **`n`** (`int`): Number of sequences to generate.
*   **`seed`** (`int`): PRNG seed.
*   **`respect`** (`ContractSet`): Biological constraints (e.g., `ga.productive()`).
*   **`strict`** (`bool`): If `True`, raises an error if a contract cannot be satisfied.

#### `.run_records(n, seed, respect, strict, expose_provenance)`
Same as `.run()`, but returns a `SimulationResult` wrapper.
*   **`expose_provenance`** (`bool`): If `True`, includes `truth_v_call`, etc.

#### `.stream_records(n, seed, respect, strict)`
Lazily yields record dictionaries (ideal for millions of sequences).

---

## The `SimulationResult` Class

A list-like container for simulated records.

| Method | Description |
|--------|-------------|
| `.to_csv(path, airr_strict)` | Exports to TSV (default) or CSV. Use `airr_strict=True` for 1-based coordinates. |
| `.to_fasta(path)` | Exports sequences to FASTA format. |
| `.to_dataframe(airr_strict)`| Returns a pandas DataFrame (requires `pandas`). |
| `.records` | Access the underlying list of dictionaries. |
| `.outcomes` | Access the underlying list of native Rust `Outcome` objects. |

---

## Distributions

Many GenAIRR methods accept a flexible distribution shape for counts and lengths:

1.  **Fixed:** `10` — Every sequence gets exactly 10.
2.  **Uniform Range:** `(5, 15)` — A random integer between 5 and 15 inclusive.
3.  **Empirical:** `[(0, 0.5), (5, 0.4), (10, 0.1)]` — A list of `(value, weight)` tuples.

---

## Global Functions

*   **`ga.list_configs()`**: Returns a list of all 106 built-in configuration names.
*   **`ga.productive()`**: Returns the standard contract bundle for productivity.
*   **`ga.set_seed(seed)`**: Sets the global PRNG seed for the current session.
