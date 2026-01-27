# Step Deep-Dives

This page explains the internal logic of every built-in pipeline step — what problem it solves, how it works algorithmically, and what data it reads from and writes to the `SimulationContainer`.

For parameter details and constructor signatures, see the [API Reference](../api_reference.md) and [Parameter Reference](../parameter_reference.md).

---

## SimulateSequence

### Problem It Solves

Creates an immunoglobulin or TCR sequence from scratch by simulating V(D)J recombination with somatic hypermutation. This is the foundational step — every pipeline starts here.

### Algorithm

1. **Select alleles**: Choose V, D (if applicable), J, and C alleles from the DataConfig. Selection is weighted by empirical gene usage frequencies, unless specific alleles are provided via `specific_v`, `specific_d`, `specific_j`.

2. **Trim segments**: For each segment, sample trim amounts from the DataConfig's trimming distributions:
    - V: trimmed at the 3' end only
    - D: trimmed independently at both 5' and 3' ends (validated so combined trims don't exceed allele length)
    - J: trimmed at the 5' end only

3. **Generate NP junctions**: Create non-templated nucleotide regions between segments using a first-order Markov chain model:
    - Sample the NP region length from `NP_lengths`
    - Sample the first base from `NP_first_bases`
    - Generate subsequent bases using position-specific transition probabilities from `NP_transitions`
    - NP1: between V and D (or V and J for light chains)
    - NP2: between D and J (heavy chains and TCR-beta only)

4. **Assemble**: Concatenate trimmed-V + NP1 + trimmed-D + NP2 + trimmed-J into the complete sequence.

5. **Validate**: Run the [FunctionalityValidator](productivity.md) to check reading frame, stop codons, and anchor residues.

6. **Retry** (if `productive=True`): If the sequence is non-functional, generate a completely new random sequence. Retry up to 25 times. If all fail, issue a warning and proceed.

7. **Apply mutations**: Pass the assembled sequence to the mutation model (S5F or Uniform), which introduces substitutions at the sampled rate. If the mutation model has `productive=True`, mutations that would create stop codons or destroy anchors are rejected.

8. **Re-validate**: Run the validator again on the mutated sequence to update productivity flags.

### Container Writes

| Field | Value |
|-------|-------|
| `sequence` | The mutated nucleotide sequence |
| `v_call`, `d_call`, `j_call`, `c_call` | Lists containing the selected allele(s) |
| `v_sequence_start/end`, `d_sequence_start/end`, `j_sequence_start/end` | Segment positions in the final sequence |
| `v_germline_start/end`, `d_germline_start/end`, `j_germline_start/end` | Positions within the original germline allele |
| `v_trim_5/3`, `d_trim_5/3`, `j_trim_5/3` | Trim amounts applied |
| `junction_sequence_start/end` | CDR3 region boundaries |
| `mutations` | Dictionary of `{position: "X>Y"}` |
| `mutation_rate` | Fraction of mutated positions |
| `productive`, `stop_codon`, `vj_in_frame`, `note` | Validation results |

---

## FixVPositionAfterTrimmingIndexAmbiguity

### Problem It Solves

When the V segment is trimmed at its 3' end, random nucleotides are added in the NP1 junction. By chance, some of those NP nucleotides may be identical to the bases that were trimmed away. This makes the true V/NP boundary ambiguous — an alignment tool cannot tell where V ends and NP begins.

This step resolves the ambiguity by expanding the V boundary to absorb any matching bases.

### Algorithm

1. Extract the NP region immediately after the current V end (V_end to D_start, or V_end to J_start for light chains).
2. Retrieve the trimmed-off portion of the V germline (the last `v_trim_3` bases of the V reference).
3. Walk through the trimmed region character by character, comparing with the NP region:
    - If the first NP base matches the last trimmed base, expand V by 1 position.
    - Continue comparing until a mismatch is found or one sequence is exhausted.
4. Update `v_sequence_end`, `v_germline_end`, and `v_trim_3` to reflect the expansion.

### Example

```
Before:  V segment ends at position 280, v_trim_3 = 5
         Trimmed V tail:  ...ACGTC
         NP region start: ACGT...

         The first 4 NP bases (ACGT) match the trimmed tail.

After:   V segment ends at position 284, v_trim_3 = 1
         The 4 matching bases are now attributed to V, not NP.
```

!!! note "Why this matters for benchmarking"
    If ground-truth V boundaries don't account for this ambiguity, an alignment tool that correctly identifies the overlap would appear to disagree with the ground truth — producing a false-negative evaluation.

---

## FixDPositionAfterTrimmingIndexAmbiguity

### Problem It Solves

Same concept as the V fix, but applied to both ends of the D segment. D segments are short (often 10–40 bp) and heavily trimmed, so NP nucleotides frequently match the trimmed D sequence on either side.

### Algorithm

1. Extract NP1 (V_end to D_start) and NP2 (D_end to J_start).
2. Retrieve the trimmed portions:
    - 5' side: first `d_trim_5` bases of D reference
    - 3' side: last `d_trim_3` bases of D reference
3. For the **5' end**: Walk backwards through the trimmed-5 region and NP1, comparing from the end. Each match moves D_start one position earlier.
4. For the **3' end**: Walk forward through the trimmed-3 region and NP2. Each match moves D_end one position later.
5. Update `d_sequence_start/end`, `d_germline_start/end`, and `d_trim_5/3`.

!!! info "Only for heavy chains and TCR-beta"
    Light chains have no D segment. Skip this step for kappa and lambda chains.

---

## FixJPositionAfterTrimmingIndexAmbiguity

### Problem It Solves

Same concept applied to the J segment's 5' end. After J is trimmed at the 5' end, NP nucleotides before J may match the trimmed-away J sequence.

### Algorithm

1. Extract the NP region before J (D_end to J_start, or V_end to J_start for light chains).
2. Retrieve the trimmed 5' portion of J (first `j_trim_5` bases of J reference).
3. Walk backwards through both sequences comparing characters. Each match moves J_start one position earlier.
4. Update `j_sequence_start`, `j_germline_start`, and `j_trim_5`.

---

## CorrectForVEndCut

### Problem It Solves

After trimming, the remaining V sequence may be identical across multiple V alleles. For example, if two V alleles differ only in their last 10 bases and 12 bases are trimmed, both alleles produce the same observed sequence. The original V allele is no longer uniquely identifiable.

This step documents all V alleles that are indistinguishable given the observed trim amount.

### Algorithm

1. Load the **V_3_TRIM_SIMILARITY_MAP** from the DataConfig. This precomputed map stores, for each V allele and trim amount, the list of other V alleles that become identical.
2. Look up the current V allele (`v_call[0]`) and current `v_trim_3` amount.
3. Retrieve all equivalent alleles from the map.
4. Add them to `v_call`, keeping the original allele first.

### Container Writes

| Field | Change |
|-------|--------|
| `v_call` | Extended from `[original]` to `[original, equivalent1, equivalent2, ...]` |

!!! info "This doesn't change the sequence"
    This step only updates metadata. The sequence itself is unchanged.

---

## CorrectForDTrims

### Problem It Solves

Same concept as CorrectForVEndCut, but for D alleles. Since D segments are trimmed from both ends, the similarity lookup uses a 2D key: the `(d_trim_5, d_trim_3)` tuple.

### Algorithm

1. Load the **D_5_3_TRIM_SIMILARITY_MAP** from the DataConfig.
2. Look up the current D allele with both trim amounts as a tuple key.
3. Add all equivalent alleles to `d_call`.

!!! info "Only for heavy chains and TCR-beta"
    Skip for light chains.

---

## DistillMutationRate

### Problem It Solves

The raw mutation dictionary from `SimulateSequence` includes mutations at every position — including NP junctions. But NP regions are non-templated random nucleotides with no germline reference, so "mutations" there are meaningless. This step removes them and recalculates a clean mutation rate.

### Algorithm

1. Identify NP region positions:
    - For heavy chains: `[v_sequence_end, d_sequence_start)` and `[d_sequence_end, j_sequence_start)`
    - For light chains: `[v_sequence_end, j_sequence_start)`
2. Find all mutations whose positions fall within NP regions.
3. Remove those mutations from the `mutations` dictionary.
4. Recalculate: `mutation_rate = len(remaining_mutations) / len(sequence)`

### Container Writes

| Field | Change |
|-------|--------|
| `mutations` | NP-region mutations removed |
| `mutation_rate` | Recalculated (lower or equal) |

!!! warning "Must run before artifact steps"
    Corruption, N-insertion, and indels modify the sequence without being biological mutations. If DistillMutationRate runs after them, positions will be shifted and the rate will be wrong.

---

## CorruptSequenceBeginning

### Problem It Solves

Real NGS data often has imperfect 5' ends due to PCR artifacts, primer degradation, or incomplete reverse transcription. This step simulates these artifacts by modifying the beginning of the sequence.

### Algorithm

With probability `probability` (default 0.7), one of three corruption events is selected based on `event_weights`:

#### Event 1: Add Random Nucleotides

1. Sample the number of bases to add from a Beta distribution scaled by `nucleotide_add_coefficient`.
2. Generate a random nucleotide string of that length.
3. Prepend it to the sequence.
4. Shift all position annotations **right** by the added amount.
5. Shift all mutation positions right.

#### Event 2: Remove Nucleotides

1. Sample the number of bases to remove from a Beta distribution scaled by `nucleotide_remove_coefficient` (capped at V region length).
2. Remove that many bases from the sequence start.
3. Shift all position annotations **left** by the removed amount.
4. Discard mutations that were in the removed region; shift the rest left.
5. Look up V alleles that become equivalent given the now-larger effective V trim.
6. Recalculate the mutation rate.

#### Event 3: Remove Then Add

1. Perform the remove operation (as above).
2. Sample a smaller number of bases to add (using `nucleotide_add_after_remove_coefficient`).
3. Prepend random bases.
4. Check if the added bases partially recreate the removed V sequence — if so, absorb the overlap back into V.
5. Update equivalent V alleles based on the final effective trim.

After any corruption event, the step **recalculates productivity** (stop codons, reading frame, anchors) because the sequence may have changed in ways that affect these.

### Container Writes

| Field | Change |
|-------|--------|
| `sequence` | Modified (bases added/removed at 5' end) |
| `corruption_event` | `'add'`, `'remove'`, or `'remove_before_add'` |
| `corruption_add_amount` / `corruption_remove_amount` | How many bases |
| `corruption_added_section` / `corruption_removed_section` | The actual bases |
| All position fields | Shifted left or right |
| `mutations` | Shifted; removed if in deleted region |
| `v_call` | Extended with equivalent alleles |
| `productive`, `stop_codon`, `vj_in_frame` | Recalculated |

---

## EnforceSequenceLength

### Problem It Solves

Sequencing platforms have fixed read lengths (e.g., Illumina MiSeq: ~576 bp for paired 2×300). Sequences longer than the read length need to be truncated.

### Algorithm

If `len(sequence) > max_length`:

1. Calculate `trim_amount = len(sequence) - max_length`.
2. Remove `trim_amount` bases from the **5' end** (since sequencing reads from 5' to 3', 5' truncation simulates incomplete reads).
3. Shift all position annotations left by `trim_amount`.
4. Remove mutations, Ns, and indels that fall in the trimmed region; shift the rest left.
5. If trimming extends into the V region, adjust `v_germline_start` to reflect the portion of V that was lost.

### Container Writes

| Field | Change |
|-------|--------|
| `sequence` | Truncated to `max_length` |
| All position fields | Shifted left |
| `mutations`, `Ns`, `indels` | Entries removed or shifted |
| `v_germline_start` | Increased if V was trimmed |

---

## InsertNs

### Problem It Solves

Real sequencing data contains positions where the base caller could not determine the nucleotide with confidence. These positions are reported as 'N'. This step simulates that.

### Algorithm

With probability `probability` (default 0.5):

1. Calculate `num_replacements = int(len(sequence) * n_ratio)`.
2. Select that many random positions in the sequence.
3. Replace each selected base with 'N'.
4. Log each replacement in the `Ns` dictionary as `{position: "original_base > N"}`.
5. Remove any mutations at N-replaced positions (since the original base is now unknown).
6. **V allele disambiguation**: For N positions within the V region, map them back to germline positions and query the **V_N_AMBIGUITY_CORRECTION_GRAPH**. This graph identifies which V alleles become indistinguishable when specific discriminating positions are obscured by Ns. Add those alleles to `v_call`.
7. Recalculate the mutation rate (counting Ns within the V-J region as additional uncertainty).

### Container Writes

| Field | Change |
|-------|--------|
| `sequence` | Selected positions replaced with 'N' |
| `Ns` | Dictionary of N positions and original bases |
| `mutations` | Entries at N positions removed |
| `v_call` | Extended with N-ambiguous equivalent alleles |
| `mutation_rate` | Increased to account for N uncertainty |

---

## InsertIndels

### Problem It Solves

Sequencing and PCR errors can introduce small insertions and deletions. This step simulates those artifacts.

### Algorithm

With probability `probability` (default 0.5):

1. **Find valid positions**: All positions from 0 to `j_sequence_end`, excluding:
    - NP junction regions (non-germline; indels there are not meaningful)
    - Positions already occupied by Ns
    - Positions already mutated
2. **Sample indel count**: Uniformly from 1 to `max_indels`.
3. **For each indel** (processed sequentially to handle shifting):
    - Choose insertion or deletion based on `insertion_probability` / `deletion_probability` weights.
    - **Insertion**: Pick a random nucleotide, insert at the position. Shift all downstream positions right by 1. Log as `"I < {base}"`.
    - **Deletion**: Record the deleted base, remove it from the sequence. Shift all downstream positions left by 1. Log as `"D > {base}"`.
    - Update mutation positions, N positions, and remaining indel positions to account for the shift.
4. **Recalculate productivity**: Indels can introduce frameshifts that destroy the reading frame.

### Container Writes

| Field | Change |
|-------|--------|
| `sequence` | Modified (bases inserted/deleted) |
| `indels` | Dictionary of `{position: "I < X" or "D > X"}` |
| All position fields | Shifted per indel |
| `mutations`, `Ns` | Positions shifted |
| `productive`, `stop_codon`, `vj_in_frame` | Recalculated |

---

## ShortDValidation

### Problem It Solves

D segments can be trimmed so aggressively that only a few bases remain. When the D region is very short (e.g., 1–4 bp), the D allele assignment is essentially random — many D alleles share such short subsequences. This step flags those cases.

### Algorithm

1. Calculate D length: `d_sequence_end - d_sequence_start`.
2. If length < `short_d_length` (default: 5 bp): replace `d_call` with `["Short-D"]`.

### Container Writes

| Field | Change |
|-------|--------|
| `d_call` | Replaced with `["Short-D"]` if D is too short |

!!! info "Only for heavy chains and TCR-beta"
    Light chains have no D segment.

---

## FilterTCRDJAmbiguities

### Problem It Solves

In TCR-beta, there are biological constraints on which D and J families can appear together. After position corrections and allele equivalence expansion, the `d_call` list may contain D alleles from families that are incompatible with the selected J allele. This step removes those invalid combinations.

### Algorithm

1. Extract the family numbers from all J alleles in `j_call`.
2. For each D allele in `d_call` (except the original, which is always kept):
    - Extract its family number.
    - If the family number is not compatible with any J family, remove it.
3. Update `d_call` with the filtered list.

### Container Writes

| Field | Change |
|-------|--------|
| `d_call` | Filtered to remove incompatible D-J family combinations |

!!! info "TCR-beta only"
    This step is specific to T-cell receptor beta chains and should not be included in BCR pipelines.
