---
title: Biological Events
sidebar_label: Biological Events
---

# Biological Events

V(D)J recombination is not a simple cut-and-paste operation. Several biological phenomena can modify the rearranged sequence during or shortly after recombination, before the cell ever encounters antigen. GenAIRR models three of these through the `.recombine()` phase: D-gene inversion, receptor revision, and allele locking (for controlled experiments).

---

## D-gene inversion

### What it models

During V(D)J recombination, the D gene segment is flanked by recombination signal sequences (RSSs) on both sides. Because the D segment has RSSs on both ends, it can be incorporated in either orientation — sense or antisense. When the D is inserted in the reverse-complement orientation, this is called **D-gene inversion**.

D-gene inversion is a well-documented biological phenomenon. It's estimated to occur in ~5-15% of heavy chain rearrangements. The inverted D produces a different amino acid sequence in the CDR3 region, contributing to antibody diversity. This is one of the mechanisms that expands the CDR3 repertoire beyond what forward-orientation D genes alone could provide.

D-gene inversion only applies to VDJ chains (IGH, TCRB, TCRD) — VJ chains have no D segment.

### How it works

After the V, D, and J segments are assembled and trimmed, the C engine checks whether this sequence should have its D inverted. It draws against the configured probability. If triggered:

1. The engine walks the D-segment nodes in the ASeq linked list and collects all bases and germline positions
2. It writes them back in reverse order, complementing each base (A↔T, C↔G)
3. Both the `current` and `germline` fields are updated — the germline reference for this sequence is now the reverse complement of the D allele

This happens **before** somatic hypermutation, so any subsequent S5F mutations are applied to the inverted D sequence and will use the inverted 5-mer contexts.

The `d_call` field in the output still reports the original D allele name (e.g., `IGHD2-21*02`), since the gene identity hasn't changed — only its orientation.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `prob` | 0.15 | Probability of D-gene inversion per sequence (0-1) |

### Usage

```python
from GenAIRR import Experiment
from GenAIRR.ops import with_d_inversion

result = (
    Experiment.on("human_igh")
    .recombine(
        with_d_inversion(0.15),
    )
    .run(n=1000, seed=42)
)
```

---

## Receptor revision

### What it models

Receptor revision (also called **receptor editing** in the B cell context) is a secondary rearrangement event. After an initial V(D)J rearrangement, if the resulting receptor is autoreactive (binds self-antigens), the cell can attempt a rescue by replacing the V gene with a different one. The new V gene replaces most of the original V segment, but a short stretch at the 3' end of the original V — the **footprint** — is preserved at the V-D junction.

This creates a chimeric V region: the bulk of the sequence comes from the new V allele, but the last 5-20 bases before the NP1 region still match the original V. In real data, this signature is detectable by aligning the V region and finding a breakpoint where the best-matching allele changes.

Receptor revision is primarily a B cell phenomenon but has been reported in T cells as well.

### How it works

When triggered, the C engine:

1. **Selects a new V allele** from the config's allele pool. It picks a V allele from a different gene family than the original — the replacement must come from a distinct V gene (e.g., IGHV1-2 → IGHV3-11, never IGHV1-2 → IGHV1-69). It makes up to 50 attempts to find a suitable replacement.

2. **Determines the footprint length** by drawing uniformly from `[footprint_min, footprint_max]`. The footprint is the number of bases at the 3' end of V that are preserved from the original allele.

3. **Replaces the upstream V bases** — all V-segment nodes before the footprint are overwritten with the new allele's sequence. Both `current` and `germline` are updated to the new allele's bases, so the germline alignment reflects the replacement allele.

4. **Updates the V call** — `v_call` in the output reports the new (replacement) V allele, not the original.

The replacement is **length-preserving**: it overwrites V bases in-place without adding or removing nodes, so all downstream positions remain stable. If the new V allele is shorter than the region being replaced, some positions may not be overwritten.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `prob` | 0.05 | Probability of receptor revision per sequence (0-1) |
| `footprint` | (5, 20) | (min, max) length of the preserved original-V footprint in nucleotides |

### Usage

```python
from GenAIRR import Experiment
from GenAIRR.ops import with_receptor_revision

result = (
    Experiment.on("human_igh")
    .recombine(
        with_receptor_revision(prob=0.05, footprint=(5, 20)),
    )
    .run(n=1000, seed=42)
)
```

### The footprint signature

In a receptor-revised sequence, the `v_call` reports the replacement allele. If you align the full V region to this allele, you'll see a good match everywhere except the last few bases (the footprint), which come from the original V. This mismatch pattern is how receptor revision is detected in real repertoire data.

The footprint length is typically short — 5-20 nucleotides — because it represents the physical overlap between the old and new V gene segments at the recombination site.

---

## Allele locking

### What it models

In normal simulation, GenAIRR samples V, D, and J alleles from the config's allele frequency distribution. Allele locking overrides this — you specify exactly which allele(s) to use. This is useful for:

- **Benchmarking**: test a tool's accuracy on a known allele
- **Controlled experiments**: isolate the effect of a specific V gene on downstream analysis
- **Generating training data**: create datasets balanced across specific alleles

### How it works

The `using()` op accepts keyword arguments for each segment type: `v`, `d`, `j`, and `c`. Each can be a single allele name (string) or a list of allele names. When a list is provided, the engine samples uniformly from that list (ignoring the config's frequency distribution).

Allele names must exactly match the names in the config. These follow the format used by the config's source database (typically IMGT nomenclature).

### Parameters

| Keyword | Description |
|---------|-------------|
| `v` | V allele name(s) to use |
| `d` | D allele name(s) to use (VDJ chains only) |
| `j` | J allele name(s) to use |
| `c` | C (constant region) allele name(s) to use |

### Usage

```python
from GenAIRR import Experiment
from GenAIRR.ops import using

# Lock to a specific V allele
result = (
    Experiment.on("human_igh")
    .recombine(
        using(v="IGHVF10-G50*04"),
    )
    .run(n=100, seed=42)
)

# All sequences use the specified V allele
assert all(rec["v_call"] == "IGHVF10-G50*04" for rec in result)
```

You can lock any combination of segments:

```python
result = (
    Experiment.on("human_igh")
    .recombine(
        using(
            v="IGHVF10-G50*04",
            j="IGHJ4*02",
        ),
    )
    .run(n=100, seed=42)
)
```

Or provide a list to sample from a subset:

```python
result = (
    Experiment.on("human_igh")
    .recombine(
        using(v=["IGHVF10-G50*04", "IGHVF10-G35*03"]),
    )
    .run(n=100, seed=42)
)

# Each sequence uses one of the two specified V alleles
v_calls = set(rec["v_call"] for rec in result)
print(v_calls)  # e.g. {"IGHVF10-G50*04", "IGHVF10-G35*03"}
```

:::note
The allele names must match exactly what the config provides. Use a small test run to see available allele names for your config:
```python
result = Experiment.on("human_igh").run(n=50, seed=1)
v_alleles = sorted(set(r["v_call"] for r in result))
print(v_alleles)
```
:::

---

## Combining recombine ops

All `.recombine()` ops can be combined freely:

```python
from GenAIRR import Experiment
from GenAIRR.ops import with_d_inversion, with_receptor_revision, using

result = (
    Experiment.on("human_igh")
    .recombine(
        with_d_inversion(0.15),
        with_receptor_revision(0.05),
        using(v="IGHVF10-G50*04"),
    )
    .run(n=1000, seed=42)
)
```

Note that `using(v=...)` sets the **initial** V allele. If receptor revision is triggered, the V allele may be replaced with a different one — so `v_call` in the output may differ from what you locked.
