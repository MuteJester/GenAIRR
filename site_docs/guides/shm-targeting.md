# SHM and mutation targeting

<p class="lead">Somatic hypermutation is GenAIRR's richest knob — you
can pick a per-base rate, choose a uniform or context-aware model,
target specific V/D/J/NP buckets, drill further into FWR vs CDR
sub-regions, and read every realised mutation back off the AIRR
record. This guide is the one place that ties it all together.</p>

## What SHM means in GenAIRR

GenAIRR treats **biological SHM** as a first-class mechanism distinct
from sequencing and library artefacts. Three things follow from
that:

- `.mutate(...)` is the only pass that produces biological SHM. Its
  events show up on every record's `n_mutations` counter and on the
  per-segment / per-V-subregion partitions.
- The library and sequencing corruption passes
  (`.pcr_amplify`, `.ambiguous_base_calls`, `.sequencing_errors`,
  `.polymerase_indels`, `.end_loss_*prime`) **do not** increment
  `n_mutations`. They surface separately on `n_pcr_errors`,
  `n_quality_errors`, `n_indels`, and the `end_loss_*_length`
  fields.
- Recombination-stage mechanisms (D inversion, receptor revision)
  rewrite the assembled sequence but are *also* excluded from
  `n_mutations` — they happen before any biology-stage mutation.

If a downstream tool wants the count of "biology that actually
changed the sequence under SHM", `n_mutations` is the answer. For
anything else, look at the corruption counters.

## Minimal mutation

A working SHM pass is one method call:

```python
import GenAIRR as ga

result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .mutate(model="s5f", rate=0.03)
      .run_records(n=100, seed=1)
)

rec = result[0]
print(rec["n_mutations"], rec["n_v_mutations"], rec["n_d_mutations"], rec["n_j_mutations"])
```

`rate=0.03` is a per-base SHM rate that scales with sequence
length; `model="s5f"` selects context-aware mutability.

## Uniform vs S5F

Two models ship.

### `model="uniform"`

Every base proposes a mutation with equal probability per site.
Simple, fast, useful when you want a baseline SHM with no
biological context bias.

```python
exp.mutate(model="uniform", rate=0.05)        # ~5% per-base mutation
exp.mutate(model="uniform", count=12)         # exactly 12 mutations per record
```

`rate=` is per-base; `count=` is the absolute number of mutation
events per record. Use whichever shape matches your intent.

### `model="s5f"`

The canonical Shapiro 5-mer model. Each site's mutability is
weighted by its 5-base context, so hot spots concentrate in the
CDRs (where SHM actually accumulates in vivo) rather than scattering
uniformly across the V region.

```python
exp.mutate(model="s5f", rate=0.03)
```

S5F draws its 5-mer mutability table and substitution preferences
from the cartridge's SHM model block. The bundled human cartridges
all ship S5F kernels; you can inspect the support on a custom
cartridge through the manifest:

```python
manifest = cfg.cartridge_manifest()
manifest["models"]["shm"]                     # full SHM model inventory
manifest["models"]["shm"]["v_subregion_support"]
# {"available": True, "annotated_v_count": 320, "total_v_count": 320}
```

The `v_subregion_support.available` flag tells you whether the
cartridge carries the FWR/CDR boundaries you'll need for
subregion-targeted SHM (next two sections).

## Controlling where SHM can land

`segment_rates={...}` multiplies the per-site mutability by a
per-bucket factor. Four buckets — V, D, J, NP — and every site
falls into exactly one:

```python
.mutate(
    model="s5f",
    rate=0.03,
    segment_rates={"V": 1.0, "D": 0.2, "J": 0.5, "NP": 0.0},
)
```

In this configuration:

- The V region keeps the cartridge's S5F mutability untouched
  (factor 1.0).
- D bases get a 0.2× scaling — SHM is biologically rare in D.
- J bases get a 0.5× scaling.
- **NP regions are excluded entirely.** A `0.0` factor doesn't
  just make NP sites unlikely to mutate — it removes them from the
  proposal support before constraint admissibility runs. Zero
  means *cannot land here*, not *unlikely to land here*.

Omitted buckets default to `1.0`. The default (no kwarg at all, or
`{}`, or explicit all-ones) is byte-identical to the pre-targeting
engine — opt-in only.

## V-subregion targeting

For V, you can drill further: the five canonical IMGT subregions
— `FWR1` / `CDR1` / `FWR2` / `CDR2` / `FWR3` — each get their own
multiplicative factor on top of `segment_rates["V"]`.

```python
.mutate(
    model="s5f",
    rate=0.03,
    v_subregion_rates={"CDR": 2.0, "FWR": 0.5, "FWR2": 1.0},
)
```

Two aliases:

- `"FWR"` expands to `FWR1`, `FWR2`, `FWR3`.
- `"CDR"` expands to `CDR1`, `CDR2`.

**Alias expansion runs first, then explicit labels override.** The
example above resolves to:

```python
{"FWR1": 0.5, "FWR2": 1.0, "FWR3": 0.5, "CDR1": 2.0, "CDR2": 2.0}
```

— FWR2 wins its own value because the explicit label beats the
alias.

A few hard rules:

- **Subregion annotations are required.** Calling
  `v_subregion_rates={...}` against a cartridge whose V alleles
  don't carry `subregions` raises `ValueError` at builder time.
  The bundled human OGRDB cartridges have them; check the
  manifest's `v_subregion_support` block before assuming.
- **Non-V sites are unaffected.** D, J, and NP sites multiply by
  identity `1.0` — `v_subregion_rates` only modulates V positions.
- **Zero excludes.** `{"CDR": 0.0}` removes CDR1 and CDR2 from the
  proposal support, exactly like a zero `segment_rates` entry.
- The composition is multiplicative: a V-site SHM weight is
  `base_mutability × segment_rates["V"] × v_subregion_rates[label]`.

## Reading mutation counters

Every record carries three layers of mutation counters that
partition exactly.

### Global

```python
rec["n_mutations"]    # total biological SHM events
```

This is the canonical "how many mutations did this record
accumulate" number.

### Per-segment partition

```python
rec["n_v_mutations"]   # V-region SHM
rec["n_d_mutations"]   # D-region SHM
rec["n_j_mutations"]   # J-region SHM
rec["n_np_mutations"]  # NP1 + NP2 combined
```

`n_v + n_d + n_j + n_np == n_mutations` by construction. The
validator's `MutationCountSumMismatch` issue fires if this ever
breaks — it's a load-bearing invariant of the SHM partition.

### Per-V-subregion partition

```python
rec["n_fwr1_mutations"]
rec["n_cdr1_mutations"]
rec["n_fwr2_mutations"]
rec["n_cdr2_mutations"]
rec["n_fwr3_mutations"]
rec["n_v_unannotated_mutations"]
```

These six partition `n_v_mutations` (not `n_mutations`):
`n_fwr1 + n_cdr1 + n_fwr2 + n_cdr2 + n_fwr3 + n_v_unannotated ==
n_v_mutations` on every record.

`n_v_unannotated_mutations` catches V-bucket SHM events that
don't land in any of the five canonical labels. Three legitimate
cases produce non-zero values:

- A V allele in your cartridge has no subregion annotations (legacy
  or hand-authored — the bundled OGRDB cartridges annotate
  everything).
- The SHM event lands in the **V-side CDR3 contribution stretch**
  between FWR3.end and the V allele's end. The five canonical
  labels deliberately stop at FWR3; junctional V-tail SHM goes
  into the unannotated bucket.
- An SHM pass runs after an indel pass and mutates an inserted V
  base.

The first case is rare; the second is normal background on
bundled cartridges and explains why `n_v_unannotated_mutations`
isn't always zero.

!!! info "Two-bucket aggregates are not exposed"
    `n_cdr_mutations` and `n_fwr_mutations` are deliberately NOT
    AIRR fields. Derive them downstream when you need them:

    ```python
    df["n_cdr_mutations"] = df["n_cdr1_mutations"] + df["n_cdr2_mutations"]
    df["n_fwr_mutations"] = (
        df["n_fwr1_mutations"] + df["n_fwr2_mutations"] + df["n_fwr3_mutations"]
    )
    ```

## Interaction with `productive_only`

GenAIRR follows a **constrain-before-propose** discipline: when
`productive_only()` is in the pipeline, the mutation pass never
proposes substitutions that would break the productive triad
(in-frame junction, no junction stop, V Cys preserved, J W-or-F
preserved). The constraint masks the proposal support; rejected
sites simply don't participate in the draw.

What you observe as a user:

- Every record still has `productive: True`.
- `n_mutations` is the count of **accepted biological mutations**
  — substitutions that landed and survived the constraint mask.
- Asking for `rate=0.1` on a cartridge whose junction would
  saturate at 5% productive-safe sites is **not an error**; the
  draw simply caps at what the constraint allows.

This composes cleanly with `segment_rates` and `v_subregion_rates`
— targeting an SHM rate at CDR2 while requiring productivity
simply means the CDR2-targeted proposals that would break the
junction frame are masked out. There's no retry loop; the
constraint prunes the candidate set at sample time.

## Replay and validation

Two reproducibility surfaces stay in your favour automatically:

- **The rate vectors fold into the plan signature.** Both
  `segment_rates` and `v_subregion_rates` are part of the plan
  signature, so replaying a trace against an `Experiment` with a
  different rate vector fails *before* any choices are consumed
  — you get a clear plan-signature mismatch, not silent divergence.
- **The counter partition is validator-enforced.** When you call
  `result.validate_records(refdata)`, the validator independently
  re-derives every counter from the event ledger and checks
  every partition equality (`n_v + n_d + n_j + n_np == n_mutations`,
  the six per-V-subregion counters summing to `n_v_mutations`,
  etc.). A discrepancy fires `MutationCountSumMismatch` or
  `VSubregionMutationCountSumMismatch` issues with the per-record
  breakdown.

Neither rate vector enters `refdata_content_hash` — they're
per-experiment knobs, not cartridge identity.

## Comparing two SHM models

Holding the recombination seed fixed and varying *only* the
mutation pass gives you a clean A/B isolation of what the SHM
choice does to the downstream distribution. Useful when picking
between S5F variants, calibrating a custom model, or
stress-testing an aligner against non-canonical mutation
patterns.

The pattern: two panels, identical seed, swapped model. Bases
the mutation pass doesn't touch are byte-identical across the
runs — only the SHM events differ.

```python
import GenAIRR as ga

def panel(model: str):
    return (
        ga.Experiment.on("human_igh")
          .recombine()
          .mutate(model=model, count=(10, 20))
          .run_records(n=5_000, seed=42, expose_provenance=True)
    )

a = panel("s5f").to_dataframe()
b = panel("uniform").to_dataframe()

# Sanity check: every row pairs by sequence_id with the same truth_v_call
assert (a["truth_v_call"] == b["truth_v_call"]).all()
```

Two diagnostics are most informative:

### V-identity densities

`v_identity` is the strongest single axis — every record yields
one value in [0, 1] and the shape is sensitive to whether the
model concentrates mutations in motifs or scatters them.

```python
import matplotlib.pyplot as plt
import seaborn as sns

fig, ax = plt.subplots()
sns.kdeplot(data=a["v_identity"], label="S5F",     ax=ax)
sns.kdeplot(data=b["v_identity"], label="uniform", ax=ax)
ax.set_xlabel("V identity")
ax.legend()
```

### 5-mer motif distribution at mutated positions

Density plots can hide structural differences. Bucket every
mutated position by its 5-mer context (centred on the mutated
base) and the gap becomes obvious — S5F over-represents AID
hotspots like `WRC` / `GYW`, uniform scatters flat across
contexts.

```python
from collections import Counter

# Mutated positions aren't a record field — derive them by diffing
# sequence_alignment against germline_alignment, then take the 5-mer
# context around each mutated position from germline_alignment.
def motif_counts(df, k=5):
    counts = Counter()
    half = k // 2
    for _, row in df.iterrows():
        sa, germ = row["sequence_alignment"], row["germline_alignment"]
        for p in range(half, len(germ) - half):
            if sa[p] != germ[p] and sa[p] != "-" and germ[p] not in ("-", "N"):
                counts[germ[p - half : p + half + 1]] += 1
    return counts

ma = motif_counts(a)
mb = motif_counts(b)
# Diff the top 20 motifs — S5F should over-represent WRC, GYW
```

## Common mistakes

A handful of issues that show up repeatedly with the SHM surface.

**Expecting PCR or quality errors to increase `n_mutations`.**
They never do. `n_mutations` is biology only. Use
`n_pcr_errors` / `n_quality_errors` / `n_indels` for
library-and-sequencing artefacts; use `n_mutations` only for SHM.

**Setting every rate to zero.** `segment_rates={"V": 0, "D": 0,
"J": 0, "NP": 0}` doesn't produce zero mutations with a "no
sites available" error — it makes the proposal support empty, so
the pass draws nothing. Every record will have `n_mutations: 0`.
If you wanted that, fine; if you didn't, check which buckets you
zeroed.

**Using `v_subregion_rates` on an unannotated cartridge.** Raises
`ValueError` at builder time. The bundled OGRDB cartridges
annotate every human V allele; if you authored a custom cartridge
via `ReferenceCartridgeBuilder.from_fasta`, run
`.infer_v_subregions()` (or supply the intervals manually) before
adding `v_subregion_rates` to your pipeline.

**Expecting CDR / FWR two-bucket counters.** They're not AIRR
fields. The catalogue stops at the five canonical IMGT labels
plus `n_v_unannotated_mutations`. Combine them downstream when
you need the two-bucket view (see the admonition above).

**Mistaking `n_v_unannotated_mutations` for a cartridge bug.**
On bundled cartridges this counter is usually non-zero but small
— it's the V-side CDR3 contribution stretch (everything past
FWR3.end is unannotated by design). Only a *large* unannotated
count on a cartridge whose `v_subregion_support.available` is
`True` is a real signal worth investigating.

## Where to go next

- **[Corruption and sequencing artefacts](corruption-sequencing.md)**
  — the observation-stage mechanisms (PCR / indel / end-loss /
  N-injection / strand) that surface on separate counters from
  `n_mutations`.
- **[The Experiment builder](experiment-builder.md)** — how SHM
  composes with the other passes in the pipeline.
- **[`validate_records`](../validation/validate-records.md)** — what
  the counter-partition invariants actually check.
- **[Reference cartridge](../concepts/reference-cartridge.md)** —
  the V-subregion annotation surface that drives subregion
  targeting.
- **[Your first AIRR record](../getting-started/first-airr-record.md)**
  — the field-by-field catalogue these counters live on.
