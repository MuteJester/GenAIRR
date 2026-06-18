# Recombination and junction biology

<p class="lead">Recombination is the single pass that creates the
V(D)J receptor - it picks alleles, trims their coding ends, fills
NP regions, optionally adds palindromic P-bases, and assembles
the molecule. This guide is the conceptual bridge between the
high-level Experiment builder and the detailed Junction N/P
guide: what the recombination pass actually produces, how to read
the resulting AIRR fields, and where the boundaries sit between
recombination biology and observation-stage corruption.</p>

## What recombination creates

A single `.recombine()` call produces everything a downstream
analysis would call "the receptor":

| Surface | What it carries |
|---|---|
| **V/D/J assignment** | The truth alleles the engine sampled, exposed on `v_call` / `d_call` / `j_call` (and `truth_*_call` when provenance exposure is on) |
| **Trims** | Per-end exonuclease byte counts on `v_trim_3`, `d_trim_5`, `d_trim_3`, `j_trim_5` |
| **N additions** | The NP1 (and NP2 on VDJ) sequences sampled from the cartridge's typed N-base / N-length models |
| **P additions** | The optional palindromic bases derived from each coding end, counted on `p_v_3_length`, `p_d_5_length`, `p_d_3_length`, `p_j_5_length` |
| **Assembled regions** | The final structural regions (`V`, `NP1`, `D`, `NP2`, `J`) carrying their coordinates on the record |
| **Junction** | The canonical V Cys → J W/F + 3 window exposed on `junction` / `junction_aa` / `junction_length` |

Everything before `.mutate()`, a clonal fork (`clonal_lineage`,
`clonal_repertoire`, legacy `expand_clones`), or any corruption pass is
recombination's responsibility. After the pass, the molecule is "finished" in
the recombination-biology sense - SHM, clonal expansion, and library prep happen
on top of it.

## A minimal recombination

```python
import GenAIRR as ga

result = ga.Experiment.on("human_igh").recombine().run_records(n=5, seed=1)
rec = result[0]

print(rec["v_call"], rec["d_call"], rec["j_call"])
print(rec["junction_aa"])
print(rec["v_trim_3"], rec["d_trim_5"], rec["d_trim_3"], rec["j_trim_5"])
print(rec["np1"], rec["np2"])
print(rec["productive"], rec["vj_in_frame"], rec["stop_codon"])
```

Five records each go through one recombination pass - sample
alleles, trim, fill NP regions, assemble. No SHM, no PCR errors,
no end-loss. Every field above is recombination's output.

## VJ vs VDJ chains

GenAIRR splits loci into two shapes:

| Chain type | Locus | Shape |
|---|---|---|
| `BCR_HEAVY` | IGH | VDJ |
| `TCR_BETA` | TRB | VDJ |
| `TCR_DELTA` | TRD | VDJ |
| `BCR_LIGHT_KAPPA` | IGK | VJ |
| `BCR_LIGHT_LAMBDA` | IGL | VJ |
| `TCR_ALPHA` | TRA | VJ |
| `TCR_GAMMA` | TRG | VJ |

On a VJ chain there is no D segment and no NP2 region. The
record dict still has `d_call` / `d_trim_5` / `d_trim_3` / `np2`
/ `np2_length` / `p_d_5_length` / `p_d_3_length` fields for
schema stability, but they carry their empty/zero defaults on
every record. Code that may receive either shape can safely read
the fields without branching first; the values are just zero or
empty string when D and NP2 don't exist.

Authoring a VJ cartridge with D-end keys on `trims` or
`p_nucleotide_lengths`, or with `NP2` on `np_lengths`, raises
`ValueError` at the cartridge attach time - the validator catches
the impossible shape before the engine gets the chance.

## Allele calls and ambiguity

Each call (`v_call` / `d_call` / `j_call`) is either a single
allele name or a comma-separated **tie set** when an independent
walker can't disambiguate between alleles from the assembled
sequence alone:

```python
rec["v_call"]   # 'IGHV1-2*02,IGHV1-2*04'  ← tie set: two alleles equally consistent
```

Two facts to know about call provenance:

- **The truth allele lives at position 0 of the tie set.** When
  GenAIRR projects an AIRR record, the truth allele (the one the
  engine actually sampled) is placed first if it's in the tie
  set. Consumers that want a single canonical call read
  `rec["v_call"].split(",", 1)[0]`.
- **The calls are evidence-derived, not handed out.** The engine
  knows what it sampled, but the AIRR projection runs an
  independent walker over the assembled sequence and reports
  what *that walker* would conclude. Under high SHM or short
  trims the walker can genuinely lose the truth allele from the
  tie set - that's not a bug, it's the alignment evidence being
  ambiguous.

When you want the original sampled allele explicitly, opt into
provenance exposure:

```python
result = exp.run_records(n=5, seed=1, expose_provenance=True)
rec = result[0]
rec["truth_v_call"], rec["truth_d_call"], rec["truth_j_call"]
```

`truth_*_call` is single-valued (no tie sets) and carries the
allele the engine sampled regardless of what the walker would
report.

### Controlling the allele universe

The cartridge defines which alleles exist - bundled cartridges
ship 100+ V alleles per locus, and the recombination pass
samples uniformly across them by default. Two surfaces narrow
the universe:

```python
# Per-experiment subset:
ga.Experiment.on("human_igh") \
   .recombine() \
   .restrict_alleles(v=["IGHV1-2*02", "IGHV3-23*01"])

# Per-cartridge weights:
cfg.reference_models = ReferenceEmpiricalModels(
    allele_usage=AlleleUsageSpec(v={"IGHV1-2*02": 100.0, "IGHV3-23*01": 60.0}),
)
```

`.restrict_alleles(...)` is a hard subset for the experiment;
allele-usage weights are an empirical bias from the cartridge.
See [Reference cartridge](../concepts/reference-cartridge.md)
for the full surface.

## Trims

The recombination pass models exonuclease trimming at each
biologically-active coding end:

```python
rec["v_trim_3"]   # bases trimmed from V's 3' coding end
rec["d_trim_5"]   # bases trimmed from D's 5' coding end (VDJ only)
rec["d_trim_3"]   # bases trimmed from D's 3' coding end (VDJ only)
rec["j_trim_5"]   # bases trimmed from J's 5' coding end
```

These are the four ends GenAIRR's recombination model trims.
Two fields exist on the AIRR record but are **always zero** on
every record GenAIRR produces:

```python
rec["v_trim_5"]   # always 0  - no recombination-stage 5' V trim in v1
rec["j_trim_3"]   # always 0  - no recombination-stage 3' J trim in v1
```

The four populated fields are driven by the cartridge's
`ReferenceEmpiricalModels.trims` plane:

```python
cfg.reference_models = ReferenceEmpiricalModels(
    trims={
        "V_3": EmpiricalDistributionSpec([(0, 5.0), (1, 3.0), (2, 1.0)]),
        "J_5": EmpiricalDistributionSpec([(0, 5.0), (1, 2.0)]),
    },
)
```

For per-experiment overrides, use `Experiment.trim(v_3=..., d_5=...,
d_3=..., j_5=...)` chained after `.recombine()`.

### Trims vs end-loss

A common confusion. Two completely separate mechanisms:

| Mechanism | Engine pass | AIRR fields | Biology stage |
|---|---|---|---|
| **Recombination trim** | `TrimPass` | `v_trim_3`, `d_trim_5`, `d_trim_3`, `j_trim_5` | Recombination (germline exonuclease) |
| **End-loss** | `EndLossPass` (corruption-stage) | `end_loss_5_length`, `end_loss_3_length` | Observation (5'/3' adapter / primer loss during library prep + sequencing) |

The trim fields capture biology; the end-loss fields capture
library-prep artefacts. If your downstream pipeline conflates
them, the per-allele identity calculations will be wrong.

## N and P additions

NP regions and P-nucleotides are also recombination biology, but
they're rich enough to deserve their own page -
[**Junction N/P additions**](junction-additions.md) covers the
typed cartridge planes that drive them, the engine pass order,
and the AIRR fields where realised state lands. Quick summary
here:

- **N additions** - random bases sampled from the cartridge's
  `np_lengths` + `np_bases` planes. Realised on `np1` /
  `np1_length`, `np2` / `np2_length` (VDJ only).
- **P additions** - deterministic palindromic bases derived from
  each coding end. Realised lengths on `p_v_3_length` /
  `p_d_5_length` / `p_d_3_length` / `p_j_5_length`.

Crucially: `np1` and `np2` strings are P-clean. The structural
NP region span excludes P bases by construction, so
`len(rec["np1"]) == rec["np1_length"]` on every record. Don't
re-derive lengths by string counting - they may include neither
P nor NP attribution as you'd expect.

## Productivity

`productive` is a triad with four predicates that must all hold:

| Predicate | Check |
|---|---|
| `vj_in_frame` | `junction_length % 3 == 0` |
| (no junction stop) | No `*` in `junction_aa` (only meaningful when in-frame) |
| V Cys preserved | The V Cys anchor is the right amino acid after recombination |
| J W/F preserved | The J W-or-F anchor is the right amino acid after recombination |

`productive == True` iff all four hold. When `productive == False`,
the issue payload from `validate_records` names which predicate
fired - `OutOfFrame`, `JunctionStopCodon`, `VAnchorAaChanged`,
or `JAnchorAaChanged`. The four flags on the record decompose the
verdict so you can debug a non-productive storm without re-running.

### Constrained sampling

The simplest way to ensure every record is productive is
`.productive_only()`:

```python
result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .productive_only()
      .run_records(n=1000, seed=42)
)
```

This isn't a retry loop. The engine *narrows the sample space*
at recombination time - trims, NP bases, NP lengths, P lengths
are all proposed against an admissible-only mask. The engine
never proposes a candidate that would break productivity, so
every record carries `productive: True` by construction.

When productivity can't be satisfied at a given site (rare -
usually means an unsatisfiable plan, like a cartridge whose
trim distribution leaves no in-frame option), default
permissive mode falls back to unconstrained sampling and the
record may end up non-productive. Pass `strict=True` to
`run_records(...)` to surface this as `StrictSamplingError`
instead.

## Reading recombination outputs

Every recombination-produced field on the record dict, organised
by concern:

```python
# Allele calls (tie sets when ambiguous)
rec["v_call"], rec["d_call"], rec["j_call"]
rec["truth_v_call"], rec["truth_d_call"], rec["truth_j_call"]  # only with expose_provenance=True
rec["locus"]                                                    # derived

# Trims (recombination-stage; v_trim_5 and j_trim_3 are always 0)
rec["v_trim_3"], rec["d_trim_5"], rec["d_trim_3"], rec["j_trim_5"]

# NP regions
rec["np1"], rec["np1_length"]
rec["np2"], rec["np2_length"]      # VDJ only

# P-nucleotide lengths
rec["p_v_3_length"], rec["p_d_5_length"], rec["p_d_3_length"], rec["p_j_5_length"]

# Junction
rec["junction"], rec["junction_aa"], rec["junction_length"]

# Productivity
rec["productive"], rec["vj_in_frame"], rec["stop_codon"]

# Coordinates + CIGAR
rec["v_sequence_start"], rec["v_sequence_end"]
rec["v_germline_start"], rec["v_germline_end"]
rec["v_cigar"]                      # M/I/D ops only
# (same shape for d_*, j_*)

# Assembled molecule
rec["sequence"], rec["sequence_aa"]
```

CIGAR strings only emit canonical M / I / D / S / N / P / X / = ops
no soft-clips. Coordinates are 0-based exclusive end by default;
pass `airr_strict=True` to TSV/CSV/DataFrame exports for the
AIRR-C 1-based-inclusive convention.

## Common mistakes

A handful of issues that show up repeatedly when reading
recombination output.

**Confusing trim with end-loss.** Recombination trims
(`*_trim_*`) are biology - the germline exonuclease activity
during recombination. End-loss (`end_loss_*_length`) is a
corruption-stage artefact - 5'/3' adapter or primer loss
during library prep / sequencing. They're separate fields,
separate engine passes, separate biology stages. A downstream
pipeline that conflates them will misattribute the source of
sequence shortening.

**Expecting D on VJ chains.** A VJ cartridge has no D segment.
The record dict still has `d_call` / `d_trim_5` / `d_trim_3` /
`np2` / `np2_length` / `p_d_5_length` / `p_d_3_length` fields
for schema stability, but they all carry their empty/zero
defaults. Code that reads them shouldn't branch on `chain_type`;
it should just see the empty/zero shape.

**Assuming one allele call when tie sets are possible.** Every
call (`v_call`, `d_call`, `j_call`) can be a comma-separated
list. The first entry is the truth allele when present in the
tie set; consumers that want one canonical name read
`call.split(",", 1)[0]`. If you split downstream code on a
single-call assumption, high-SHM datasets will surface tie sets
and your join will be silently wrong.

**Treating `productive` as one boolean without checking
components.** When `productive == False`, the four decomposition
flags (`vj_in_frame`, `stop_codon`, V Cys preserved, J W/F
preserved) name exactly which predicate fired. The validator's
issue payload also exposes them. Reading only `productive`
hides why a record failed; checking the components lets you
debug a non-productive storm in one pass.

## Where to go next

- **[Junction N/P additions](junction-additions.md)** - the
  detailed cartridge planes that drive N and P biology, plus the
  estimator-specific quirks.
- **[D inversion + receptor revision](recombination-editing.md)**
  advanced recombination-stage mechanisms that edit V/D truth.
- **[Reference cartridge](../concepts/reference-cartridge.md)** -
  the four-plane cartridge model that controls allele universe,
  trims, NP lengths/bases, and P lengths.
- **[SHM and mutation targeting](shm-targeting.md)** - the
  biology stage that runs *after* recombination is complete.
- **[Your first AIRR record](../getting-started/first-airr-record.md)**
  the field-by-field catalogue.
- **[Validation & reproducibility](../validation/index.md)** -
  the postcondition validator that checks every recombination
  field against the simulated outcome.
- **[The Experiment builder](experiment-builder.md)** - where
  `.recombine()` and `.productive_only()` sit in the full
  pipeline.
