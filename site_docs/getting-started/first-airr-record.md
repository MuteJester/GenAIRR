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

Each standard call names the generating allele used by the corresponding
alignment. On VJ chains, `d_call` is empty.

### Truth and evidence-based re-calls

`recalled_v_call`, `recalled_d_call`, and `recalled_j_call` contain GenAIRR's
post-mutation inference. They may differ from the generating allele or
contain comma-separated ties. They do not change the alignment's reference.
Optional `truth_*_call` fields remain compatibility aliases for the standard
calls when `expose_provenance=True`.

## Junction + productivity

```python
rec["junction"]         # 'TGTGCGAGAGATGAT…' - V Cys through J W/F + 3
rec["junction_aa"]      # 'CARDDGNRGYCSGGSCYGRCCALDYW'
rec["junction_length"]  # 78
rec["productive"]       # True
rec["vj_in_frame"]      # True
rec["stop_codon"]       # False
```

`junction` is the nucleotide window from the V Cys codon through the J W/F codon;
`junction_aa` is its in-frame translation. The productivity
flags decompose as:

- `vj_in_frame` - the junction is a whole number of codons **and** begins
  on a codon boundary of the reading frame:

  ```python
  (rec["junction_end"] - rec["junction_start"]) % 3 == 0 \
      and rec["junction_start"] % 3 == rec["coding_frame_start"] % 3
  ```

  The second clause is why this is not simply `junction_length % 3 == 0`.
  Without an indel the two agree on every record; with one they disagree
  on about a quarter, because an indel upstream of the junction shifts the
  frame while leaving the window a whole number of codons. Measured on
  human_igk: 200 of 200 agreement without indels, 54 of 200 disagreement
  with one, always with the flag reporting out-of-frame.
- `stop_codon` - any stop in the translated coding span, which is
  `sequence[coding_frame_start:coding_frame_end]`. **Not** limited to the
  junction: a record whose `junction_aa` translates cleanly can still
  carry `stop_codon: True`, with `coding_stop` in `productivity_reasons`.
  Measured on 46 of 299 hypermutated records. Three-valued over ambiguous
  bases: a codon containing an `N` is a stop if every base it could be
  makes one, not a stop if none does, and `None` otherwise - so
  `ambiguous_base_calls()` can leave this null with
  `coding_stop_unresolved` in the reasons.
- `productive` - the conjunction of **five** predicates, each of them
  exported, under three-valued logic: `False` if any is false, `True` if
  all are true, `None` otherwise.

  | predicate | field |
  |---|---|
  | the V anchor was accepted | `v_anchor_accepted` |
  | the J anchor was accepted | `j_anchor_accepted` |
  | the junction is in frame | `vj_in_frame` |
  | no stop in the coding span | `stop_codon` is `False` |
  | the J tail joins the constant region in frame | `j_tail_in_frame` |

  The last one is easy to miss and is why a record can satisfy everything
  the other four ask for and still be unproductive, with
  `j_tail_out_of_frame` as the only reason - 5 of 200 IGH records carrying
  one polymerase indel. A frame-shifted J tail cannot splice to a constant
  region, so the rearrangement does not express whatever the junction
  looks like.

  The two anchor predicates are tri-state for reasons a reader cannot see
  from the bases: an ambiguous base inside an anchor codon leaves it
  *unresolved* rather than rejected, and an insertion inside the codon
  disrupts its provenance without changing the base. That is why they are
  fields rather than something to recompute.

**Every one of these three flags is reproducible from the record.** If you
are scoring a tool against GenAIRR and it disagrees with one of them, you
can recompute ours from the fields above and see which of you is wrong -
which is the point of exporting the frame and the predicates rather than
only the verdict.

When `productive_only()` is in the pipeline, every record carries
`productive: True` by construction.

## Identity, mutation, and error counters

```python
rec["coding_frame_start"] # 0    (first base of the reading frame)
rec["coding_frame_end"]   # 316  (one past the last base assessed)
rec["j_tail_in_frame"]    # True (the fifth predicate of `productive`)
rec["v_anchor_accepted"]  # True (the first two)
rec["j_anchor_accepted"]  # True
rec["v_identity"]        # 1.0  (matches / total over V's CIGAR M/D)
rec["d_identity"]        # 1.0
rec["j_identity"]        # 1.0
rec["n_mutations"]       # 0  (global SHM count)
rec["n_observable_mutations"]  # 0  (the subset a caller can recover)
rec["n_v_mutations"]     # 0  (SHM partitioned by carried segment)
rec["n_d_mutations"]     # 0
rec["n_j_mutations"]     # 0
rec["n_np_mutations"]    # 0
rec["n_pcr_errors"]      # 0  (from .pcr_amplify pass)
rec["n_quality_errors"]  # 0  (from .sequencing_errors pass)
rec["n_indels"]          # 0  (from .polymerase_indels pass)
```

`v_identity` is the fraction of that segment's alignment matching the
germline — matches over every column the segment's CIGAR covers, so an
inserted or deleted column counts against it. A segment with no bases in
the read reports `None` rather than `1.0`.

`v_score` / `d_score` / `j_score` and `v_support` / `d_support` /
`j_support` are **always `None`**. They are an aligner's outputs: how well
the chosen allele explains the read, and how confident the choice is.
GenAIRR assembles from a known allele, so there is no hypothesis to score
and no call to be uncertain about. The columns are present because the AIRR
schema has them — if you are scoring an aligner, compute these from its
alignment, not from ours.

`n_mutations` is the biological SHM total; the per-segment fields
partition it exactly (`n_v + n_d + n_j + n_np == n_mutations` by
construction). PCR / quality / indel artefacts are reported on
separate counters - they're library/sequencing-stage noise, not
biological mutation.

Like `n_mutations`, the two error counters report the substitutions that
**changed a base**, not the number requested: `pcr_amplify(count=5)`
typically reports three or four, because a replacement is drawn from
A/C/G/T including the base already present. `ambiguous_base_calls()` masks
positions with `N` and is not reported on either counter; its count is
available from the trace (`corrupt.ns.count`).

It counts the substitutions that **changed a base**. A uniform
replacement is drawn from A/C/G/T including the base already present,
and those draws are not counted — so `mutate(model="uniform", count=7)`
typically reports five or six, because `count` is the number of attempts.
No single `model="s5f"` write is ever a no-op: that kernel does not draw
the base already present.

`n_observable_mutations` is the number a mutation caller can recover from
`sequence`, and it is the field to score a caller against. It is below
`n_mutations` by three intended effects: a mutation on an untemplated
junction base (an N or P addition) changes a base and is counted, but has
no germline base to differ from, so nothing can recover it — about 5% of
records that carry a mutation; two mutations can hit one position with the
later one restoring the germline base; and a position lost to end loss or
a deletion is no longer in the molecule.

It is attributable to hypermutation alone. A PCR misincorporation, a
sequencing miscall or an `N` mask also leaves a position differing from
its germline, and none of them moves this count — those are reported on
`n_pcr_errors` and `n_quality_errors`.

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
- **Per-subregion recoverable counters** - the same six regions,
  counting what a caller can recover rather than what happened:
  `n_fwr1_observable_mutations` and friends. Use these for a
  CDR/FWR selection-pressure ratio computed from a read; the event
  counters above agree with a read on 17 of 60 records once twelve
  nucleotides are lost from the 5' end.
- **Custom metadata** - `Experiment.with_metadata(experiment_id=..., tissue=...)`
  stamps arbitrary user fields onto every record.

See the [AIRR Record concept page](../concepts/airr-record.md) for
the full field catalogue and derivation rules.

---

## Next step

→ [Export the results](export-results.md) - write your records to
TSV, FASTA, FASTQ, paired-end FASTQ, or pandas DataFrame.
