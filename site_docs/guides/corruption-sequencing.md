# Corruption and sequencing artefacts

<p class="lead">Real sequencing data doesn't ship clean. PCR
amplifies errors, sequencers miscall bases, polymerases insert
and delete, adapter loss clips the ends, and library prep flips
the strand half the time. GenAIRR models each of these as a
distinct observation-stage pass — separate from biological SHM,
separately counted, and individually tunable to match the
artefact profile of the sequencer you're benchmarking against.</p>

## What counts as observation-stage corruption

Seven passes belong in the observation-stage block. Each models a
distinct artefact mechanism, each has its own AIRR provenance, and
each is independently tunable:

| Pass | What it does | Counter / field |
|---|---|---|
| `.pcr_amplify(...)` | PCR-cycle substitution errors | `n_pcr_errors` |
| `.sequencing_errors(...)` | Sequencer base-call errors (lowercases bytes) | `n_quality_errors` |
| `.ambiguous_base_calls(...)` | N-base injection (replaces bytes with `N`) | (no dedicated counter — see below) |
| `.polymerase_indels(...)` | Polymerase-stage insertions and deletions | `n_indels`, `n_v_indels`, `n_d_indels`, `n_j_indels` |
| `.end_loss_5prime(...)` | 5' adapter / primer loss | `end_loss_5_length` |
| `.end_loss_3prime(...)` | 3' adapter / primer loss | `end_loss_3_length` |
| `.random_strand_orientation(...)` | Strand flip during library prep | `rev_comp` |

Plus `.paired_end(...)`, which isn't corruption — it's a read-layout
projection — but it sits in the same descendant-phase block of the
pipeline and is usually configured alongside the corruption passes.
See [Paired-end reads and FASTQ](paired-end-fastq.md).

## Biological mutation vs artefact

The single most important distinction: **`.mutate(...)` is the
only pass that increments `n_mutations`**. Every corruption pass
above is observation-stage; their counts surface on dedicated
fields. A downstream pipeline that conflates them will misattribute
real biology to library prep and vice versa.

| Mechanism | Pipeline stage | Counter |
|---|---|---|
| Somatic hypermutation | Biology — `.mutate(...)` | `n_mutations` (plus per-segment `n_v_mutations` etc.) |
| PCR substitution | Library — `.pcr_amplify(...)` | `n_pcr_errors` |
| Sequencer base-call error | Sequencer — `.sequencing_errors(...)` | `n_quality_errors` |
| N-base injection | Sequencer / quality flagging — `.ambiguous_base_calls(...)` | (none on the record) |
| Polymerase indel | Library — `.polymerase_indels(...)` | `n_indels` + per-segment partition |
| Adapter / primer loss | Library — `.end_loss_*prime(...)` | `end_loss_*_length` |
| Strand flip | Library — `.random_strand_orientation(...)` | `rev_comp` |

`n_mutations` is biology, every other field is artefact. If you
want a "total perturbation budget" you sum them downstream — they
don't aggregate automatically.

See [SHM and mutation targeting](shm-targeting.md) for the biology
side of the boundary in detail.

## PCR and sequencing errors

Two substitution-style corruption passes — same shape, different
biology.

### PCR errors

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .pcr_amplify(count=(0, 3))
      .run_records(n=100, seed=1)
)

rec = result[0]
rec["n_pcr_errors"]    # int — number of PCR substitutions on this record
```

Each PCR error is a single-base substitution drawn uniformly across
the molecule. The `count` argument follows the standard
`int` / `(low, high)` / `[(value, weight), ...]` pattern.

### Sequencing errors

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .sequencing_errors(count=(0, 5))
      .run_records(n=100, seed=1)
)

rec["n_quality_errors"]    # int — number of sequencer base-call errors
```

`sequencing_errors` models sequencer-quality miscalls. The pass
**lowercases the affected bytes** in the assembled sequence as a
visible quality marker. So `sequence` may end up with mixed-case
bytes: uppercase = clean, lowercase = sequencer-flagged
low-quality position. Downstream pipelines that want to mask
quality-flagged bases can do
`sequence.upper()` to normalise.

## N corruption / ambiguous bases

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .ambiguous_base_calls(count=(0, 2))
      .run_records(n=100, seed=1)
)

rec["sequence"]    # may contain 'N' characters at corrupted positions
```

N-base injection replaces individual bytes with `N`. Useful when
the downstream tool you're benchmarking handles ambiguity codes
specifically (allele-call walkers, junction translators).

**There is no `n_ns` counter on the AIRR record.** The pass leaves
no provenance field; the only signal is the `N` characters in
`sequence`. If you need a count, derive it: `rec["sequence"].count("N")`.

A few facts about how `N` corruption composes with downstream
analysis:

- **Allele-call walkers handle `N` as a wildcard.** The GenAIRR
  walker classifies `N` as matching any canonical reference base
  (per the standard convention), so an N-injected record can still
  produce a meaningful `v_call`. Tie sets may widen.
- **Junction translation produces `X` codons.** An `N`-containing
  codon translates to `X`; `junction_aa` and `sequence_aa` may
  carry `X` characters. The `productive` predicate still evaluates
  correctly — `X` codons don't satisfy the stop-codon check but
  also don't satisfy the V-Cys or J-W/F anchor preservation checks
  if they land at the anchor positions.

## Polymerase indels

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .polymerase_indels(count=(0, 2), insertion_prob=0.5)
      .run_records(n=100, seed=1)
)

rec["n_indels"]      # total polymerase indel events
rec["n_v_indels"]    # partitioned by carried segment
rec["n_d_indels"]
rec["n_j_indels"]
```

Each event is either an insertion (random base inserted at a
random position) or a deletion (single base removed), drawn per
`insertion_prob`. The `n_indels` total partitions exactly into the
per-segment counters: `n_v + n_d + n_j == n_indels` for events
that land in V / D / J. Events in NP1 / NP2 roll into `n_indels`
only (they're not allocated to V/D/J).

### Indels vs end-loss

A common confusion. Both shorten / extend the sequence in some
way; the biology is completely different:

| Mechanism | What it does | AIRR fields |
|---|---|---|
| **Polymerase indel** | Random single-base insertions and deletions *inside* the molecule | `n_indels`, per-segment partition |
| **End-loss** | Block deletion from the molecule's 5' or 3' end | `end_loss_5_length`, `end_loss_3_length` |

If your downstream pipeline conflates them, the per-segment
identity and coordinate calculations will be wrong. Use
`n_indels` for "did the polymerase mis-step in the middle of the
molecule"; use `end_loss_*_length` for "how much got chewed off
the ends".

## End-loss

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .end_loss_5prime(length=(0, 8))
      .end_loss_3prime(length=(0, 4))
      .run_records(n=100, seed=1)
)

rec["end_loss_5_length"]    # bases removed from 5' end
rec["end_loss_3_length"]    # bases removed from 3' end
```

End-loss models adapter / primer loss during library prep —
the read starts a few bases in from the molecule's actual 5' end
and / or terminates a few bases short of the 3' end. The
`length` argument follows the standard
`int` / `(low, high)` / `[(value, weight), ...]` shape.

### `primer_trim_*prime` is a backwards-compatible alias

```python
# These two are exactly equivalent:
.end_loss_5prime(length=(0, 8))
.primer_trim_5prime(length=(0, 8))    # legacy alias — same pass, same fields
```

The legacy `primer_trim_*prime` names predate the engine's
distinction between recombination-stage exonuclease trim
(`TrimPass`, the `v_trim_3` family) and observation-stage end
loss (`EndLossPass`, the `end_loss_*_length` family). New code
should use `end_loss_*prime` for clarity.

### End-loss vs recombination trim

Another common confusion. Recombination trims and end-loss are
two separate passes with two separate sets of fields:

| Mechanism | Engine pass | AIRR fields | Biology stage |
|---|---|---|---|
| **Recombination trim** | `TrimPass` | `v_trim_3`, `d_trim_5`, `d_trim_3`, `j_trim_5` | Recombination (germline exonuclease) |
| **End loss** | `EndLossPass` | `end_loss_5_length`, `end_loss_3_length` | Observation (library prep / sequencing) |

Don't conflate them — they capture different biology and a
downstream pipeline that adds them together gets the wrong total.
See the [Recombination + junction biology](recombination-junction.md#trims-vs-end-loss)
guide for the detailed comparison.

## Random strand orientation

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .random_strand_orientation(prob=0.5)
      .run_records(n=100, seed=1)
)

rec["rev_comp"]    # bool — True if this record was flipped
```

At `prob=0.5`, half the records get reverse-complemented in place.
The `rev_comp` flag tells you which way the strand went;
`sequence` and `sequence_aa` reflect the post-flip orientation;
`junction`, `np1`, `np2` similarly reflect the post-flip bytes.

A few facts to know:

- **The flip is projection-level, not biology.** Real biology
  isn't being inverted; the read just happens to come off the
  sequencer in the other orientation. Allele truth and all the
  per-segment provenance fields stay correct; coordinates remap
  to the flipped sequence frame.
- **Paired-end windows are sliced AFTER the strand flip.** If
  `.paired_end(...)` follows `.random_strand_orientation(...)`,
  R1 and R2 are sliced from the post-flip molecule, so they
  reflect the strand decision automatically. Don't apply a
  second flip downstream.

## Ordering with clonal families

**All seven corruption passes are descendant-phase.** They model
per-read artefacts — every descendant of a clone gets independent
PCR errors, independent indel events, independent end-loss
draws, etc.

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .expand_clones(n_clones=10, per_clone=20)
      .mutate(model="s5f", rate=0.03)        # biology — descendant-phase
      .pcr_amplify(count=(0, 3))             # corruption — descendant-phase
      .polymerase_indels(count=(0, 2))
      .ambiguous_base_calls(count=(0, 2))
      .sequencing_errors(count=(0, 5))
      .end_loss_5prime(length=(0, 8))
      .end_loss_3prime(length=(0, 4))
      .random_strand_orientation(prob=0.5)
      .paired_end(r1_length=150, insert_size=300)
      .run_records(seed=1)
)
```

Calling any of them *before* `.expand_clones(...)` raises
`ValueError` at chain time. See [Clonal families](clonal-families.md)
for the ancestor / descendant phase discipline in full.

## Per-platform calibrated profiles

Different sequencing platforms break sequences in characteristic
ways: MiSeq adds quality noise at the 3′ end, NovaSeq drops
bases uniformly with a slight 5′ bias, PacBio HiFi inserts
homopolymer indels. Three starting profiles for the platforms
researchers ask about most often:

```python
# Counts are per-record ranges. For a ~350bp read, count=(0, 2) Ns
# is roughly 0.1% N-density.
PROFILES = {
    "miseq_300pe": dict(
        end_loss_5=(5, 25),
        end_loss_3=(8, 45),       # 3′ decay heavier than 5′
        ambig_ns=(0, 2),
        indels=(0, 1),
    ),
    "novaseq_150pe": dict(
        end_loss_5=(3, 15),
        end_loss_3=(3, 15),
        ambig_ns=(0, 1),
        indels=(0, 1),
    ),
    "pacbio_hifi": dict(
        end_loss_5=(0, 5),
        end_loss_3=(0, 5),
        ambig_ns=(0, 1),
        indels=(0, 3),            # HiFi: low subs, higher indels
    ),
}

def apply_profile(exp, name):
    p = PROFILES[name]
    return (
        exp.end_loss_5prime(length=p["end_loss_5"])
           .end_loss_3prime(length=p["end_loss_3"])
           .ambiguous_base_calls(count=p["ambig_ns"])
           .polymerase_indels(count=p["indels"])
    )

result = apply_profile(
    ga.Experiment.on("human_igh").recombine().mutate(count=(5, 20)),
    "miseq_300pe",
).run_records(n=5000, seed=42)
```

The set of passes you call **is** the profile — these are
starting points, not exhaustive: add `pcr_amplify` or
`sequencing_errors` rates when the platform calls for them.

### Calibrating from your own reads

Don't trust a preset blindly. Run a small panel under the
profile and match two distributions to a sample of real reads
from your sequencer:

**Match the N-density.**

```python
n_rate_real = (real_reads.str.count("N") / real_reads.str.len()).mean()
n_rate_sim  = (sim["sequence"].str.count("N") / sim["sequence_length"]).mean()
# Adjust the ambiguous_base_calls count range until these match
```

**Match the length distribution.**

```python
# End-loss is the dominant driver of read-length variance
sim_lens  = sim["sequence_length"]
real_lens = real_reads.str.len()
# Match means + variances by widening or tightening the (min, max)
# ranges on end_loss_5prime / end_loss_3prime
```

## Validation and replay

`validate_records` re-derives every corruption counter from the
event ledger and trace addresses, then compares against the
record's reported values. Issue kinds the catalogue covers:

| Counter | Issue kind on mismatch |
|---|---|
| `n_pcr_errors` | `NPcrErrorsMismatch` |
| `n_quality_errors` | `NQualityErrorsMismatch` |
| `n_indels` | `NIndelsMismatch` |
| `n_v_indels` / `n_d_indels` / `n_j_indels` | per-segment indel-count mismatches |
| `end_loss_5_length` / `end_loss_3_length` | `EndLossLengthMismatch` |

The `rev_comp` flag is derived from the IR; mismatch surfaces
similarly. Every corruption mechanism records its choices on
stable trace addresses (`corrupt.pcr.count`, `corrupt.quality.count`,
`corrupt.indel.count`, `corrupt.end_loss.5`, `corrupt.end_loss.3`,
etc.), so trace replay reproduces the same artefacts deterministically.

Paired-end FASTQ export writes the projected `r1_sequence` /
`r2_sequence` fields verbatim — by AIRR-projection time the
corruption + strand-flip + paired-end pipeline has already
finished, so the exported reads carry all the artefacts the
record's counters report. See [Validation hub](../validation/index.md)
for the validator's full scope.

## Common mistakes

A handful of issues that show up repeatedly with corruption.

**Expecting PCR or sequencing errors to increment `n_mutations`.**
They never do. `n_mutations` is biology only. PCR substitution
surfaces on `n_pcr_errors`; sequencer miscalls surface on
`n_quality_errors`. If your downstream pipeline reads
`n_mutations` as a "total perturbation budget", it's leaving the
artefact counters off.

**Confusing end-loss with recombination trim.** Recombination
trim is biology (germline exonuclease, fields `v_trim_3` etc.);
end-loss is library-prep artefact (`end_loss_*_length`). They're
separate passes, separate fields, separate biology stages. The
[Recombination + junction biology](recombination-junction.md#trims-vs-end-loss)
guide has the side-by-side comparison.

**Putting corruption before `.expand_clones()`.** All seven
corruption passes are descendant-phase — they're per-read
artefacts that need to vary within a clonal family. The DSL
rejects this at chain time with the uniform message: "<method>
must be called after expand_clones(); it is descendant-specific
and must be sampled independently for each clone member."

**Reverse-complementing R2 again after random strand orientation.**
When `.random_strand_orientation(...)` is in the pipeline,
`sequence` (and `r1_sequence` / `r2_sequence` if paired-end is
on) already reflect the post-flip orientation. R2 is also
already reverse-complemented at projection time (that's the
paired-end-output convention). Applying any further flips
downstream produces wrong reads. The validator's
`PairedEndWindowMismatch` check enforces the per-record R2 RC
state.

**Looking for an `n_ns` counter.** `ambiguous_base_calls` doesn't
ship one. The only signal is the `N` characters in `sequence` —
count them yourself with `rec["sequence"].count("N")`.

## Where to go next

- **[SHM and mutation targeting](shm-targeting.md)** — the biology
  side of the boundary; the per-segment / per-V-subregion
  partition.
- **[Paired-end reads and FASTQ](paired-end-fastq.md)** — the
  read-layout projection that usually sits alongside corruption.
- **[Clonal families](clonal-families.md)** — the
  ancestor / descendant phase discipline that gates every
  corruption pass.
- **[Recombination + junction biology](recombination-junction.md)**
  — the trim vs end-loss boundary in detail.
- **[Validation & reproducibility](../validation/index.md)** — the
  validator's issue catalogue covering every corruption counter.
- **[Your first AIRR record](../getting-started/first-airr-record.md)**
  — the field-by-field catalogue including every counter named on
  this page.
