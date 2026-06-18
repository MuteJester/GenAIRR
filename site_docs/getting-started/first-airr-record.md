# Your first AIRR record

<p class="lead">Every record GenAIRR emits is a dict with 50+
AIRR-format fields. This page walks the most important ones on a
real record so you can read your own output confidently.</p>

## One record

```python
import GenAIRR as ga

result = ga.Experiment.on("human_igh").recombine().run_records(n=1, seed=0)
rec = result[0]
print(rec)
```

The output is a Python dict. A trimmed view (the full dict has 50+
keys; we'll cover them by category below):

```python
{
    "sequence_id":   "seq0",
    "sequence":      "gaggtgcagctggtg...",
    "sequence_aa":   "EVQLVESGGGLVQPG...",
    "locus":         "IGH",
    "v_call":        "IGHVF10-G38*04",
    "d_call":        "IGHD2-15*01",
    "j_call":        "IGHJ2*01",
    "junction":      "TGTGCGAGAGATGATGGAAATAGAGGCTACTGCAGT...",
    "junction_aa":   "CARDDGNRGYCSGGSCYGRCCALDYW",
    "junction_length": 78,
    "productive":    True,
    "vj_in_frame":   True,
    "stop_codon":    False,
    "v_identity":    1.0,
    "n_mutations":   0,
    ...
}
```

## Sequence

```python
rec["sequence"]       # 'gaggtgcagctggtg…' - assembled nucleotide
rec["sequence_aa"]    # 'EVQLVESGGGLVQPG…' - codon-rail translation
rec["sequence_id"]    # 'seq0' - auto-assigned; carries to FASTQ headers
rec["locus"]          # 'IGH' - derived from v_call / j_call
```

`sequence` is what a sequencer would have produced. Uppercase
bytes are clean; lowercase ones are corruption-pass markers
(quality errors, N-base injection). `sequence_aa` is a codon-rail
translation - stops emit `*`, ambiguous codons emit `X`.

## V / D / J calls

```python
rec["v_call"]    # 'IGHVF10-G38*04'
rec["d_call"]    # 'IGHD2-15*01'
rec["j_call"]    # 'IGHJ2*01'
```

Each call is a single allele name or a comma-separated *tie set*
when an independent walker can't disambiguate from the sequence
alone. The first entry of the tie set is the truth allele - the
allele the engine actually used - so consumers that want a single
canonical name can read `rec["v_call"].split(",", 1)[0]`.

### Truth vs evidence-derived calls

GenAIRR's calls are **derived from the same evidence an aligner
would see**, not handed out from the engine's internal state. The
engine knows which allele it sampled (the *truth*), but the AIRR
projection runs an independent walker over the assembled sequence
and reports what that walker would conclude. Under high SHM or
short trims the walker can land on a tie set that doesn't include
the truth - that's biologically correct: the alignment evidence is
genuinely ambiguous. The truth-allele-first ordering convention
preserves the ground truth when it's in the set; the tie set
expresses the residual ambiguity downstream tools must handle.

## Junction + productivity

```python
rec["junction"]         # 'TGTGCGAGAGATGAT…' - V Cys through J W/F + 3
rec["junction_aa"]      # 'CARDDGNRGYCSGGSCYGRCCALDYW'
rec["junction_length"]  # 78
rec["productive"]       # True
rec["vj_in_frame"]      # True
rec["stop_codon"]       # False
```

`junction` is the canonical V Cys → J W/F+3 nucleotide window;
`junction_aa` is its in-frame translation. The productivity
flags decompose as:

- `vj_in_frame` - `junction_length % 3 == 0`.
- `stop_codon` - any stop in `junction_aa` (meaningful only when in-frame).
- `productive` - full triad: in-frame ∧ no junction stop ∧ V Cys
  preserved ∧ J W-or-F preserved.

When `productive_only()` is in the pipeline, every record carries
`productive: True` by construction.

## Identity, mutation, and error counters

```python
rec["v_identity"]        # 1.0  (matches / total over V's CIGAR M/D)
rec["d_identity"]        # 1.0
rec["j_identity"]        # 1.0
rec["n_mutations"]       # 0  (global SHM count)
rec["n_v_mutations"]     # 0  (SHM partitioned by carried segment)
rec["n_d_mutations"]     # 0
rec["n_j_mutations"]     # 0
rec["n_np_mutations"]    # 0
rec["n_pcr_errors"]      # 0  (from .pcr_amplify pass)
rec["n_quality_errors"]  # 0  (from .ambiguous_base_calls pass)
rec["n_indels"]          # 0  (from .polymerase_indels pass)
```

`n_mutations` is the biological SHM total; the per-segment fields
partition it exactly (`n_v + n_d + n_j + n_np == n_mutations` by
construction). PCR / quality / indel artefacts are reported on
separate counters - they're library/sequencing-stage noise, not
biological mutation.

## Paired-end fields (when requested)

```python
rec["read_layout"]    # 'single' by default; 'paired_end' after .paired_end(...)
rec["r1_sequence"]    # '' unless paired_end is in the pipeline
rec["r2_sequence"]    # ''
rec["r1_start"]       # 0
rec["r1_end"]         # 0
rec["r2_start"]       # 0
rec["r2_end"]         # 0
rec["insert_size"]    # 0
```

These fields are always present on the record dict for schema
stability, but they're populated only when `Experiment.paired_end(
r1_length=..., r2_length=..., insert_size=...)` is in the pipeline.
Otherwise they carry their zero / empty defaults.

## The full schema

The fields above cover the most common analysis surfaces. A few
others worth knowing about as you go deeper:

- **CIGAR strings** - `v_cigar` / `d_cigar` / `j_cigar`. Only M / I /
  D ops are emitted; no soft-clips.
- **Coordinate fields** - `v_sequence_start/end`,
  `v_germline_start/end`, etc. The sequence-frame coordinates plus
  the matching germline-frame coordinates per segment.
- **AIRR provenance flags** - `d_inverted` (D in reverse-complement),
  `receptor_revision_applied` (post-recombine V replacement),
  `is_contaminant` (replaced by background).
- **Per-V-subregion counters** - `n_fwr1_mutations` / `n_cdr1_mutations`
  / `n_fwr2_mutations` / `n_cdr2_mutations` / `n_fwr3_mutations`
  / `n_v_unannotated_mutations`. They partition `n_v_mutations`.
- **Custom metadata** - `Experiment.with_metadata(experiment_id=..., tissue=...)`
  stamps arbitrary user fields onto every record.

See the [AIRR Record concept page](../concepts/airr-record.md) for
the full field catalogue and derivation rules.

---

## Next step

→ [Export the results](export-results.md) - write your records to
TSV, FASTA, FASTQ, paired-end FASTQ, or pandas DataFrame.
