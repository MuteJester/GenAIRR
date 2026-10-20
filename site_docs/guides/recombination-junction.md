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

## Processing original coding ends (opt in)

`process_coding_ends()` models the net outcome of processing the original
junction-facing coding ends. Supply a finite distribution of signed outcomes
for each end:

| Outcome | Meaning |
|---|---|
| Positive `d` | Remove `d` original reference bases; retain no P extension at that end |
| Zero | Preserve the original end without a P extension |
| Negative `-p` | Preserve the original reference end and retain `p` palindromic bases templated from that original tip |

The following small synthetic example declares complete coding ends and fixed
outcomes. Its weights are demonstration inputs, not fitted biological defaults.

```python
import GenAIRR as ga
from GenAIRR import _engine

reference = _engine.RefDataConfig.vj()
reference.add_v_allele("V-demo*01", "V-demo", b"AAACCCTGT", anchor=6)
reference.add_j_allele("J-demo*01", "J-demo", b"TTTGGG", anchor=0)

experiment = (
    ga.Experiment.on(reference)
    .recombine(np1_lengths=[(0, 1.0)])
    .process_coding_ends(
        v_3=[(-1, 1.0)],
        j_5=[(-2, 1.0)],
        reference_ends_are_coding_ends=True,
    )
)
record = experiment.run_records(n=1, seed=1)[0]
# V + P(V3) + P(J5) + J = AAACCCTGT + A + AA + TTTGGG
assert record["junction_aa"] == "CKF"
```

Passing `reference_ends_are_coding_ends=True` explicitly confirms that the
provided junction-facing allele boundaries are complete biological coding
ends. GenAIRR cannot establish this from a truncated reference sequence.
Omitted end distributions are fixed at zero. They do not inherit empirical
trim or P distributions. N-addition length and base models remain configured
by `recombine()` and the cartridge.

Configure processing after `recombine()` and before mutation, corruption, or
a clonal fork. It replaces automatic cartridge trim/P settings and rejects
an explicit `trim()` call in either order, including `trim(enabled=False)`.
Repeated processing configuration is rejected. VJ chains reject `d_5` and
`d_3`, even when their outcomes are zero.

For VDJ chains, the independent `d_5` and `d_3` inputs are combined into a
finite joint outcome table using products of their weights. Those ends refer
to the canonical reference orientation. D inversion reverses the processed
segment and its retained extensions; the physical D5/D3 P counters can
therefore swap relative to the canonical outcomes. Joint deletions must fit
the original allele and extensions must have enough original template.

With `productive_only()`, candidate outcomes retain their supplied relative
weights among choices with a feasible productive completion. Subsequent N
and orientation choices use the same completion-support rule. This is
sequential restriction to feasible support; it is not exact sampling from the
unconstrained model conditioned on productivity, which would require weighting
by the probability mass of every future completion. Impossible support raises
an error rather than substituting zero P.

Records from this mode include a `coding_end_processing` dictionary. Its
`choices` describe canonical signed outcomes and generating alleles;
`retained_p_at_join` and `observed_p` distinguish emitted extensions from
surviving observed P bases. `p_origins` records original template positions,
assembled ends and current query positions; removed bases retain their origin
entry with a null query position. These query coordinates use the coding
orientation. Standard AIRR NP fields continue to combine intervening N and P.
`recombination_complete` is true after V, J and any assigned D have been
assembled; it remains true after later observation loss. On partial native
outcomes this flag is false, and `retained_p_at_join` contains provisional
emitted-to-date counts, not a claim that joining has completed.

CSV and TSV exports serialize `coding_end_processing` as one JSON cell;
use `json.loads(row["coding_end_processing"])` after reading the table.
The nested dictionary declares `coordinate_system="zero_based_coding_oriented"`.
Its template and query positions remain zero-based and coding-oriented,
including in strict AIRR TSV exports and reverse-strand observations. AIRR
conversion of the standard top-level coordinate columns does not alter these
nested extension coordinates. Python records and DataFrames retain the
dictionary value.

The existing empirical NP plane is still interpreted as N-only input to the
sampler. Activating this mode does not prove that a cartridge's historical NP
estimates excluded P. Training combined AIRR NP lengths or strings as N and
then adding retained P would count those bases twice. Current NP estimators
reject explicit nonzero P provenance; missing P fields do not establish that
the training data were N-only. This mode does not infer a new N model or
calibrate end-processing weights from repertoire records.

This is an effective model of retained terminal outcomes, compatible with the
negative-deletion representation used by probabilistic repertoire models.
It does not infer or claim to simulate exact Artemis nick positions,
strand-specific resection, fill-in timing, or microhomology-mediated joining.
The signed weights are supplied by the caller; no universal biochemical
calibration is assumed. Pipelines that omit this method retain the existing
trim/P model and seeded behaviour.

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

`v_call`, `d_call`, and `j_call` identify the final generating alleles and
provide the reference frame for all alignment metadata. Post-mutation
inference is retained separately:

```python
rec["v_call"]           # 'IGHVF1-G1*01' — generating allele
rec["recalled_v_call"]  # may be a different allele or comma-separated tie
```

Re-calls put the generating allele first only when it belongs to the tied
set; otherwise names follow allele-ID order. A re-call never moves a true
junction boundary or changes the reference used to represent an indel.
`expose_provenance=True` retains the `truth_*_call` compatibility aliases,
which equal the standard generating calls.

### Controlling the allele universe

The cartridge defines which alleles exist - bundled cartridges
ship 100+ V alleles per locus, and the recombination pass
samples uniformly across them by default. Two surfaces narrow
the universe:

```python
# Per-experiment subset:
ga.Experiment.on("human_igh") \
   .recombine() \
   .restrict_alleles(v=["IGHVF1-G1*01", "IGHVF1-G2*01"])

# Per-cartridge weights:
cfg.reference_models = ReferenceEmpiricalModels(
    allele_usage=AlleleUsageSpec(v={"IGHVF1-G1*01": 100.0, "IGHVF1-G2*01": 60.0}),
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
  `np_lengths` + `np_bases` planes. The final `np1` / `np2`
  intervals combine intervening N and P bases between alignments.
- **P additions** - deterministic palindromic bases derived from
  each coding end. Realised lengths on `p_v_3_length` /
  `p_d_5_length` / `p_d_3_length` / `p_j_5_length`.

`np1` and `np2` follow the final segment alignment bounds, so
`len(rec["np1"]) == rec["np1_length"]` and likewise for NP2.
They include P bases remaining in those gaps. The P-length fields
instead preserve emitted recombination history, including bases
later removed by read loss.

## Productivity

When the required coding and junction evidence is present,
`productive` requires all of these checks:

| Predicate | Check |
|---|---|
| `vj_in_frame` | Junction length is a multiple of three and its V anchor starts in the retained V coding frame |
| J-tail frame | The observed J tail and endpoint retain the coding phase, accounting for terminal read loss |
| No coding stop | No stop codon in the complete observed V(D)J coding span, including outside the junction |
| Canonical V landmark | The original three landmark bases remain adjacent and encode Cys |
| Canonical J landmark | The original three landmark bases remain adjacent and encode Phe or Trp |

The versioned policy is `canonical_observed_vdj_v1`, with scope
`observed_coding_span`. All checks must be accepted for `True`. A known defect
produces `False`; otherwise unresolved required evidence produces `None`.
Ambiguous codons are evaluated over their possible canonical resolutions:
`NNN` can hide a stop, `GCN` cannot, and `TAR` is definitely stop.
`productivity_reasons` gives stable diagnostic codes.

Invalid reference annotations never supply a guessed coding frame. Independent
phase evidence can support a stop assessment even when a landmark is unknown.
Terminal clipping retains an observed-span prediction if all required landmarks
and phase remain; absent required landmarks yield unknown. The policy does not
predict expression, folding, or functionality outside that observed span.

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

The engine narrows the sample space
at recombination time - allele choices, trims, N bases, and N/P lengths
are checked against the active productivity constraints. Admissible
choices preserve the coding frame, stop-free coding span, and anchors.
Fixed and mixed nonzero P lengths are supported. Their deterministic
bases and future D orientation enter the completion check alongside
N-addition choices. Original local weights are renormalized over
candidates that leave a productive completion.

If no valid initial productive assignment or P-aware completion
exists, execution raises an error even in permissive mode. Requested
P additions are not silently replaced with zero. Other later steps
can use their declared no-op or sentinel behavior when no admissible edit exists.
Pass `strict=True` to `run_records(...)` to surface sampling
failures rather than those permissive fallbacks.

## Reading recombination outputs

Every recombination-produced field on the record dict, organised
by concern:

```python
# Generating alleles and separate inferred calls
rec["v_call"], rec["d_call"], rec["j_call"]
rec["recalled_v_call"], rec["recalled_d_call"], rec["recalled_j_call"]
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

Python record CIGARs contain local `M` / `I` / `D` runs, and
coordinates are 0-based with exclusive ends. TSV defaults to
`airr_strict=True`: biological coordinates become 1-based inclusive,
CIGARs gain leading query `S` and reference `N` offsets when nonzero,
and booleans are written as `T` / `F`. CSV and DataFrame exports
default to `airr_strict=False`; strict DataFrames retain booleans.

Biological query coordinates refer to the coding-oriented sequence.
For `rev_comp=True`, reverse-complement the observed `sequence` before
using these coordinates to slice it. Paired-read windows always
index the observed `sequence` and remain 0-based, half-open.

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

**Treating re-calls as generating truth.** Standard calls name the generating
alleles. Split `recalled_*_call` on commas when analyzing inferred ambiguity;
its first entry need not be the generating allele.

**Treating `productive` as one boolean without checking
components.** When `productive == False`, inspect `vj_in_frame`,
`stop_codon`, the J tail, and the conserved anchor amino acids.
A validation mismatch can also carry a decision reason. Reading only `productive`
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
