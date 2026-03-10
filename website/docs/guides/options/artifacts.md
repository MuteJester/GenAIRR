---
title: Sequencing Artifacts
sidebar_label: Artifacts
---

# Sequencing Artifacts

Real AIRR-seq data passes through a long experimental pipeline before you see it in a FASTQ file. Each step — library preparation, PCR amplification, sequencing, and base calling — can introduce artifacts that distort the original biological sequence. If you're training a model or benchmarking a tool on synthetic data, your simulation should include these artifacts, otherwise you're testing on unrealistically clean data.

GenAIRR simulates artifacts across three DSL phases, matching the order in which they occur in a real experiment:

| Phase | What it models | Ops |
|-------|---------------|-----|
| `.prepare()` | Library preparation: primer masking, UMI attachment, PCR amplification | `with_primer_mask`, `with_umi`, `with_pcr` |
| `.sequence()` | The sequencing run: 5'/3' signal loss, quality profile, reverse complement | `with_5prime_loss`, `with_3prime_loss`, `with_quality_profile`, `with_reverse_complement` |
| `.observe()` | Post-sequencing noise: indels, N-bases, contaminants | `with_indels`, `with_ns`, `with_contaminants` |

The C engine applies these in pipeline order, so a PCR error introduced in `.prepare()` can later be compounded by a sequencing error at the same position in `.sequence()`.

---

## 5' and 3' signal loss

### What it models

In Sanger and short-read sequencing, the ends of reads are often unreliable. The 5' end may lose bases due to primer mispriming, RNA degradation, or reverse transcription drop-off. The 3' end may be truncated when the sequencer runs out of reagent, or may read into the adapter sequence. These artifacts are a major source of V-gene and J-gene call errors in real data — the missing or corrupted bases at the ends make it harder to align back to the germline.

### How it works

For each sequence, the C engine randomly selects one of three event types:

| Event | Probability | What happens |
|-------|-------------|-------------|
| **Remove only** | 40% | Bases are deleted from the end of the sequence |
| **Add only** | 30% | Random bases (uniform ACGT) are prepended/appended |
| **Remove + Add** | 30% | Bases are removed, then random bases are added in their place |

The "remove + add" event is the most disruptive — it simulates adapter read-through where the true biological bases are gone and replaced by nonsense. The amount removed and added are each drawn uniformly from their respective `[min, max]` ranges.

The engine applies a safety floor: it will never remove enough bases to leave the sequence shorter than 5 nucleotides. Added bases are capped at 50. Added bases are tagged with the `SEG_ADAPTER` segment type internally, so they are distinguishable from biological sequence in the ASeq linked list, though this distinction is not directly exposed in the AIRR output.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_remove` | 1 | Minimum bases to remove |
| `max_remove` | 20 | Maximum bases to remove |
| `min_add` | 1 | Minimum random bases to add |
| `max_add` | 10 | Maximum random bases to add |

### Usage

```python
from GenAIRR import Experiment
from GenAIRR.ops import with_5prime_loss, with_3prime_loss

result = (
    Experiment.on("human_igh")
    .sequence(
        with_5prime_loss(min_remove=5, max_remove=30),
        with_3prime_loss(min_remove=5, max_remove=20),
    )
    .run(n=1000, seed=42)
)
```

5' and 3' corruption are independent — you can use one without the other. In practice, 5' loss is more common in RNA-based protocols (due to reverse transcription), while 3' loss is more common in Sanger sequencing.

### Effect on annotations

After corruption, all AIRR coordinate fields (`v_sequence_start`, `j_sequence_end`, etc.) are adjusted to reflect the new sequence boundaries. This means your annotations remain correct even after bases are removed or random bases are added at the ends. The `sequence` field contains the corrupted version — this is what a real sequencer would give you.

---

## PCR amplification errors

### What it models

Before sequencing, target DNA is amplified through polymerase chain reaction (PCR). Each cycle of PCR copies every molecule, but the polymerase occasionally introduces errors — typically at a rate of 10^-4 to 10^-5 per base per cycle for high-fidelity polymerases. Over 25-35 cycles, these errors accumulate. Unlike SHM mutations (which are context-dependent and biased toward specific motifs), PCR errors are position-independent with uniform substitution probabilities.

### How it works

The C engine computes an **effective error rate** from the per-cycle rate and cycle count:

```
effective_rate = 1 - (1 - error_rate) ^ n_cycles
```

For example, with `error_rate=1e-4` and `cycles=30`:
```
effective_rate = 1 - (1 - 0.0001)^30 = 1 - 0.997 = 0.003
```

So about 0.3% of positions will have a PCR error — roughly 1 error per 300-base sequence.

The engine then walks every position in the sequence (including NP regions, unlike SHM which only targets V/D/J positions). At each position, it draws against the effective rate. If a position is selected for error, it substitutes a uniformly random different base (equal probability among the three alternatives). Positions that are already `N` are skipped.

PCR errors are flagged with `NUC_FLAG_PCR_ERROR` internally, keeping them distinguishable from SHM mutations and sequencing errors.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `error_rate` | 1e-4 | Per-base per-cycle error probability |
| `cycles` | 30 | Number of PCR cycles |

### Usage

```python
from GenAIRR.ops import with_pcr

result = (
    Experiment.on("human_igh")
    .prepare(
        with_pcr(error_rate=1e-4, cycles=30),
    )
    .run(n=1000, seed=42)
)

rec = result[0]
print(rec["n_pcr_errors"])  # number of PCR errors introduced
print(rec["pcr_errors"])    # position:from>to record of each error
```

### Typical values

| Polymerase | Per-cycle rate | Typical cycles |
|-----------|----------------|----------------|
| Taq (standard) | ~2 x 10^-4 | 25-35 |
| High-fidelity (Phusion, Q5) | ~1 x 10^-6 | 25-30 |
| Pfu | ~1 x 10^-6 | 25-30 |

---

## Sequencing quality profile

### What it models

Illumina sequencers produce reads where base-call quality degrades along the read. The first bases are called with high confidence (Q30+, error rate ~0.001), but quality drops toward the 3' end as phasing and signal decay accumulate. This creates a characteristic position-dependent error profile that's visible in any FastQC report.

### How it works

The C engine uses a **linear interpolation** between a base error rate (at position 0) and a peak error rate (at the last position):

```
error_rate(pos) = base_rate + (pos / seq_length) * (peak_rate - base_rate)
```

At each position, the engine draws against this local error rate. If an error occurs, it applies a **transition/transversion bias**: 70% of errors are transitions (A↔G, C↔T) and 30% are transversions (A↔C, A↔T, G↔C, G↔T). This matches the known chemistry of Illumina sequencing, where the most common base-call errors are confusing purines with purines or pyrimidines with pyrimidines.

Positions that are already `N` are skipped. Errors are flagged with `NUC_FLAG_SEQ_ERROR`.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `base` | 0.001 | Error rate at the 5' end of the read (~Q30) |
| `peak` | 0.02 | Error rate at the 3' end of the read (~Q17) |

### Usage

```python
from GenAIRR.ops import with_quality_profile

result = (
    Experiment.on("human_igh")
    .sequence(
        with_quality_profile(base=0.001, peak=0.02),
    )
    .run(n=1000, seed=42)
)

rec = result[0]
print(rec["n_sequencing_errors"])  # total sequencing errors
print(rec["sequencing_errors"])    # position:from>to record
```

### Interpreting the error gradient

For a 350-nucleotide sequence with default parameters (`base=0.001`, `peak=0.02`):
- Position 0: error rate = 0.001 (Q30, very reliable)
- Position 175 (midpoint): error rate = 0.0105
- Position 350 (end): error rate = 0.02 (Q17, noticeable noise)

This means the 3' half of the read accumulates most of the errors, which is exactly what real Illumina data looks like.

---

## Indels

### What it models

Small insertions and deletions (indels) can arise from sequencing errors, PCR slippage (especially in homopolymer runs), or during library preparation. Unlike SHM point mutations, indels shift the reading frame and can make downstream annotation much harder. They're relatively rare in Illumina data (< 0.1% per position) but more common in long-read technologies.

### How it works

The engine walks every position in the sequence that belongs to a germline segment (V, D, J, or C). At each eligible position, it draws against the per-position `prob`. If triggered, it flips a coin (50/50) to decide insertion vs deletion:

- **Insertion**: a random base (uniform ACGT) is inserted after the current position. The new node is tagged with `NUC_FLAG_INDEL_INS` and inherits the segment type of its neighbor.
- **Deletion**: the current node is removed from the linked list.

Anchor positions (the conserved Cysteine and Tryptophan/Phenylalanine that define the junction boundaries) are protected — they are never targeted for indels, since losing them would make the sequence impossible to annotate.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `prob` | 0.01 | Per-position probability of an indel event |

### Usage

```python
from GenAIRR.ops import with_indels

result = (
    Experiment.on("human_igh")
    .observe(
        with_indels(prob=0.01),
    )
    .run(n=1000, seed=42)
)
```

### Effect on annotations

Indels shift all downstream positions. The AIRR coordinate fields (`v_sequence_start`, `d_sequence_start`, etc.) reflect the pre-indel positions, so they may not exactly match the post-indel sequence. This is intentional — it mirrors the situation in real data where annotations refer to the germline alignment, not the observed read.

---

## N-base insertions

### What it models

Ambiguous base calls appear as `N` in sequencing output. This happens when the sequencer cannot confidently distinguish the fluorescent signal at a given cycle — common in low-quality regions, mixed templates, or degraded samples. A sprinkling of Ns throughout a read is typical in lower-quality runs.

### How it works

The engine walks every position in the sequence. At each position (regardless of segment type), it draws against `prob`. If triggered, the current base is replaced with `N` and flagged with `NUC_FLAG_IS_N`. Positions that are already `N` are skipped.

Unlike indels, N-insertion doesn't change the sequence length — it's a substitution, not an insertion. The original base is lost.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `prob` | 0.01 | Per-position probability of replacing with N |

### Usage

```python
from GenAIRR.ops import with_ns

result = (
    Experiment.on("human_igh")
    .observe(
        with_ns(prob=0.005),
    )
    .run(n=1000, seed=42)
)
```

---

## Primer masking

### What it models

Many AIRR-seq protocols use V-gene-specific primers that anneal to the FR1 region. When the primer binds, it overwrites the first ~78 bases of the V gene with the primer sequence, which is identical to the germline. Any SHM mutations that were present in FR1 are effectively erased — the primer sequence replaces whatever was there.

This is important for SHM analysis: if your real data was generated with FR1 primers, any mutations you see in FR1 are primer artifacts (or PCR errors), not genuine SHM. Enabling primer masking in your simulation produces data that matches this protocol.

### How it works

The engine walks the V-segment nodes from the 5' end. For each node with a germline position less than `length` (default: 78, which is the IMGT FR1 boundary), if the current base differs from the germline base, it reverts the base to germline. Any error flags (`NUC_FLAG_SEQ_ERROR`, `NUC_FLAG_PCR_ERROR`) are also cleared.

The result is that FR1 looks pristine — as if no mutations ever occurred there.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `length` | 0 | Number of bases to mask (0 = full FR1 = 78 bases) |

### Usage

```python
from GenAIRR.ops import with_primer_mask, rate

result = (
    Experiment.on("human_igh")
    .mutate(rate(0.02, 0.08))
    .prepare(
        with_primer_mask(),
    )
    .run(n=1000, seed=42)
)
```

---

## UMI barcodes

### What it models

Unique Molecular Identifiers (UMIs) are short random DNA sequences (typically 8-16 nucleotides) that are ligated to each molecule before PCR amplification. They allow deduplication of PCR duplicates in downstream analysis — reads with the same UMI originated from the same pre-amplification molecule.

### How it works

The engine generates a random nucleotide barcode of the specified length (uniform ACGT at each position) and prepends it to the sequence. The UMI bases are tagged with `SEG_UMI` segment type internally. The UMI appears as the first `length` characters of the output `sequence` field.

UMI length is capped at 30 nucleotides.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `length` | 12 | UMI length in nucleotides |

### Usage

```python
from GenAIRR.ops import with_umi

result = (
    Experiment.on("human_igh")
    .prepare(
        with_umi(12),
    )
    .run(n=1000, seed=42)
)

rec = result[0]
# The first 12 bases of the sequence are the UMI barcode
print(rec["sequence"][:12])  # e.g. "TTCACACTAGAA"
```

Common UMI lengths in practice: 8 (10x Genomics), 10 (MIGEC), 12 (standard), 16 (high-diversity libraries).

---

## Reverse complement

### What it models

During library preparation, molecules can be captured in either orientation. Depending on the protocol, a fraction of reads may be reverse-complemented relative to the sense strand. This is especially common in protocols without strand-specific adapters.

### How it works

For each sequence, the engine draws against `prob`. If triggered, the entire ASeq linked list is reversed and each base is complemented (A↔T, C↔G). All flags, germline annotations, and segment tags are preserved through the reversal. The `is_reverse_complement` field is set to `True` in the output record.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `prob` | 0.5 | Probability of reverse-complementing each sequence |

### Usage

```python
from GenAIRR.ops import with_reverse_complement

result = (
    Experiment.on("human_igh")
    .sequence(
        with_reverse_complement(0.3),
    )
    .run(n=100, seed=42)
)

rc_count = sum(1 for r in result if r["is_reverse_complement"])
print(f"{rc_count}/{len(result)} reverse-complemented")
```

---

## Contaminants

### What it models

In high-throughput sequencing, a small fraction of reads may come from non-target sources: cross-sample contamination during multiplexed runs, phiX control spike-in that bleeds through demultiplexing, or environmental DNA. These reads have no relation to any immune receptor — they're pure noise.

### How it works

For each sequence, the engine draws against the contamination `rate`. If triggered, the **entire sequence is replaced** — all biological content is discarded and replaced with a contaminant. All AIRR annotations are zeroed out (no V/D/J calls, not productive, note set to `"contaminant"`).

Two contaminant sources are available:

| Source | Description |
|--------|-------------|
| `"random"` | Random ACGT nucleotides (300-500 bp) |
| `"phix"` | Random substring of the phiX174 genome (~600 bp fragment) |

The `is_contaminant` field is set to `True` in the output record.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `rate` | 0.01 | Per-sequence contamination probability |
| `source` | `"random"` | `"random"` or `"phix"` |

### Usage

```python
from GenAIRR.ops import with_contaminants

result = (
    Experiment.on("human_igh")
    .observe(
        with_contaminants(rate=0.01, source="phix"),
    )
    .run(n=1000, seed=42)
)

contam_count = sum(1 for r in result if r["is_contaminant"])
print(f"{contam_count}/{len(result)} contaminants")
```

---

## Pipeline order

The C engine applies artifact steps in a fixed order that mirrors a real experiment:

```
1. Primer mask         (.prepare)   — Revert FR1 mutations to germline
2. UMI                 (.prepare)   — Prepend random barcode
3. PCR amplification   (.prepare)   — Introduce polymerase errors
4. 5' corruption       (.sequence)  — Remove/add bases at 5' end
5. 3' corruption       (.sequence)  — Remove/add bases at 3' end
6. Quality profile     (.sequence)  — Position-dependent sequencing errors
7. Reverse complement  (.sequence)  — Flip read orientation
8. Indels              (.observe)   — Random insertions/deletions
9. N-bases             (.observe)   — Ambiguous base calls
10. Contaminants       (.observe)   — Replace with non-immune sequence
```

This order matters. For example:
- Primer masking reverts SHM mutations **before** PCR errors are added, so PCR errors can appear in FR1 even after masking
- PCR errors are introduced **before** sequencing errors, so both can accumulate at the same position
- 5'/3' corruption happens **before** quality profile, so sequencing errors are distributed across the already-corrupted sequence
- Contaminants replace the **entire** sequence, so any upstream artifacts are irrelevant for contaminated reads

## Full realistic pipeline

Combine all artifact types for a complete wet-lab simulation:

```python
from GenAIRR import Experiment
from GenAIRR.ops import (
    rate, model,
    with_primer_mask, with_umi, with_pcr,
    with_5prime_loss, with_3prime_loss, with_quality_profile,
    with_reverse_complement,
    with_indels, with_ns, with_contaminants,
)

result = (
    Experiment.on("human_igh")

    .mutate(
        model("s5f"),
        rate(0.02, 0.08),
    )

    .prepare(
        with_primer_mask(),
        with_umi(12),
        with_pcr(error_rate=1e-4, cycles=30),
    )

    .sequence(
        with_5prime_loss(min_remove=5, max_remove=30),
        with_3prime_loss(min_remove=5, max_remove=20),
        with_quality_profile(base=0.001, peak=0.02),
        with_reverse_complement(0.3),
    )

    .observe(
        with_indels(prob=0.005),
        with_ns(prob=0.005),
        with_contaminants(rate=0.01),
    )

    .run(n=1000, seed=42)
)
```

## Output fields

| Field | Type | Description |
|-------|------|-------------|
| `n_pcr_errors` | int | Number of PCR polymerase errors |
| `pcr_errors` | str | Position:from>to record of PCR errors |
| `n_sequencing_errors` | int | Number of sequencing quality errors |
| `sequencing_errors` | str | Position:from>to record of sequencing errors |
| `is_reverse_complement` | bool | Whether the read was reverse-complemented |
| `is_contaminant` | bool | Whether the sequence was replaced by a contaminant |
