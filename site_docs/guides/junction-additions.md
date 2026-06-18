# Junction N and P additions

<p class="lead">The junction between V, D, and J segments isn't just
the gap between coding ends - it's where N-addition and P-nucleotide
biology gets written. GenAIRR drives both from typed cartridge
planes, and every base ends up labelled on the AIRR record so you
can see exactly where it came from.</p>

## What gets added at the junction

Two distinct mechanisms write bases into the junction. They look
similar on a sequencer trace, but they come from different places.

| Addition | Source | Determinism | Cartridge plane |
|---|---|---|---|
| **N addition** | TdT-like template-free polymerase | Random - bases are sampled from a model | `np_lengths`, `np_bases` |
| **P addition** | Hairpin opening of the coding end | **Deterministic** - bases are the palindromic complement of the trimmed coding-end flank | `p_nucleotide_lengths` |

Both contribute bases that affect the junction sequence, junction
length, and therefore the productive-triad evaluation. The
distinction matters because:

- **N additions** are sampled from a distribution. You author the
  distribution; the engine draws values.
- **P additions** are derived bytes - once the trim and length are
  fixed, the bases are forced by the coding-end sequence. You only
  control the length, not the bases.

## Where the regions sit

GenAIRR's pass order interleaves the additions with the V/D/J
assembles so the geometry is fixed:

**VJ chains** (no D):

```text
V coding | P_V_3 | NP1 | P_J_5 | J coding
```

**VDJ chains**:

```text
V coding | P_V_3 | NP1 | P_D_5 | D coding | P_D_3 | NP2 | P_J_5 | J coding
```

A few rules to keep in mind when you reason about the geometry:

- **D inversion fires before P_D_5 is derived.** If `.invert_d(...)`
  is in the pipeline and a record draws an inverted D, the P_D_5
  P-bases are the palindromic complement of the *reverse-
  complemented* D's 5' coding end. The pipeline order is
  `recombine → invert_d → P_D_5 derivation → D assemble`, so
  P-bases always reflect the final D orientation.
- **P bases live BETWEEN structural regions.** They aren't part of
  the NP1 / NP2 region span. `np1` and `np2` strings on the AIRR
  record don't include P bases; the dedicated `p_*_length`
  fields count them.

## N length models

`ReferenceEmpiricalModels.np_lengths` controls how many bases each
NP region holds. Keys are `NP1` (always) and `NP2` (VDJ only - VJ
cartridges reject `NP2` keys at attach time):

```python
from GenAIRR import EmpiricalDistributionSpec, ReferenceEmpiricalModels

cfg.reference_models = ReferenceEmpiricalModels(
    np_lengths={
        "NP1": EmpiricalDistributionSpec([(0, 1.0), (3, 4.0), (6, 2.0)]),
        "NP2": EmpiricalDistributionSpec([(0, 1.0), (2, 3.0), (4, 2.0)]),
    },
)
```

Each `EmpiricalDistributionSpec` is a flat `[(value, weight), ...]`
list - non-negative integer values, finite positive weights. The
engine normalises internally; weights don't need to sum to 1.

Per-experiment overrides go through `Experiment.recombine`:

```python
ga.Experiment.on(cfg).recombine(
    np1_lengths=[(0, 1.0), (1, 4.0), (2, 2.0)],
    np2_lengths=[(0, 1.0), (3, 2.0)],
)
```

The kwarg wins over the cartridge plane wins over the legacy
`cfg.NP_lengths` dict wins over a uniform fallback. **Passing
`np2_lengths` on a VJ cartridge raises `ValueError`** - there's
no NP2 region to extend.

## N base models

`ReferenceEmpiricalModels.np_bases` controls **which bases** get
drawn at each NP position. Three kinds, each shipped end-to-end:

```python
from GenAIRR import NpBaseModelSpec

# Uniform - every position samples A/C/G/T uniformly.
NpBaseModelSpec(kind="uniform")

# Empirical first-base - full base composition over every position.
NpBaseModelSpec(
    kind="empirical_first_base",
    first_base={"A": 0.4, "C": 0.1, "G": 0.4, "T": 0.1},
)

# Markov - 1-step previous-base-conditional.
NpBaseModelSpec(
    kind="markov",
    first_base={"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
    transitions={
        "A": {"A": 0.1, "C": 0.4, "G": 0.4, "T": 0.1},
        "C": {"A": 0.1, "C": 0.1, "G": 0.4, "T": 0.4},
        "G": {"A": 0.4, "C": 0.1, "G": 0.1, "T": 0.4},
        "T": {"A": 0.4, "C": 0.4, "G": 0.1, "T": 0.1},
    },
)
```

Notes on each kind:

- **`uniform`** is the engine baseline - byte-identical to the
  pre-typed-model behaviour when no spec is authored.
- **`empirical_first_base`** uses the supplied categorical at
  *every* position (the historical name is preserved; the
  biologically correct estimate is the full base composition).
- **`markov`** uses `first_base` at position 0 and the transition
  matrix for position 1+. **All four A/C/G/T from-base rows must
  be populated** - the spec validator rejects partial matrices at
  attach time.

Attach the spec to either NP region independently:

```python
cfg.reference_models = ReferenceEmpiricalModels(
    np_bases={
        "NP1": NpBaseModelSpec(kind="markov", first_base=..., transitions=...),
        # NP2 omitted - falls back to uniform
    },
)
```

!!! warning "Legacy `NP_transitions` / `NP_first_bases` are not auto-lifted"
    Bundled cartridges still carry the legacy dicts as documented
    orphan fields. The engine no longer reads them; they're
    inspectable via the manifest's `legacy_np_first_bases_present`
    / `legacy_np_transitions_present` flags but ignored at run
    time. If you want Markov behaviour today, populate
    `ReferenceEmpiricalModels.np_bases` explicitly.

## P-nucleotide length models

`ReferenceEmpiricalModels.p_nucleotide_lengths` controls **how many
P bases** get derived at each junction end. Four keys on VDJ
chains (`V_3`, `D_5`, `D_3`, `J_5`), two on VJ chains (`V_3`,
`J_5`). The validator rejects D-end keys on VJ at attach time.

```python
cfg.reference_models = ReferenceEmpiricalModels(
    p_nucleotide_lengths={
        "V_3": EmpiricalDistributionSpec([(0, 0.7), (1, 0.2), (2, 0.1)]),
        "D_5": EmpiricalDistributionSpec([(0, 0.8), (1, 0.2)]),
        "D_3": EmpiricalDistributionSpec([(0, 0.8), (1, 0.2)]),
        "J_5": EmpiricalDistributionSpec([(0, 0.85), (1, 0.15)]),
    },
)
```

**Bases are deterministic.** Once a record draws a P length at a
given end, the bases are the reverse-complement palindrome of the
adjacent coding-end's trimmed flank. You don't get to author the
bases - only the length distribution and (indirectly via the trim
distribution) which bases the palindrome will reflect.

Empty / omitted P-end keys mean no P-pass runs at that end -
byte-identical to the pre-P-plane baseline. Authoring
`p_nucleotide_lengths` is opt-in.

## Reading outputs

Every record carries the realised N and P state on dedicated AIRR
fields:

```python
rec["np1"]              # 'ACGT'  - NP1 nucleotide string
rec["np1_length"]       # 4
rec["np2"]              # 'TAGC'  - NP2 string (VDJ only; empty on VJ)
rec["np2_length"]       # 4

rec["p_v_3_length"]     # 1      - number of P bases at V_3 end
rec["p_d_5_length"]     # 0
rec["p_d_3_length"]     # 1
rec["p_j_5_length"]     # 0
```

Two things to know:

- **`np1` / `np2` are P-clean.** The structural NP region span
  excludes P bases by construction, so `len(rec["np1"]) ==
  rec["np1_length"]` on every record GenAIRR produces - even
  under a max-P plane.
- **P bases don't have their own string fields in v1.** Per-base P
  strings (`p_v_3`, `p_d_5`, etc.) and an aggregate
  `n_p_nucleotides` counter are deliberately deferred - the
  `p_*_length` fields are the only P provenance v1 ships.

## Interaction with `productive_only`

Both N and P additions affect junction length, frame, and content,
so they participate in the productive-triad constraint. When
`productive_only()` is in the pipeline:

- **N base draws are constraint-masked.** A position whose
  candidate bases would all push the junction out-of-frame, into
  a stop codon, or off-anchor drops out of the proposal support
  before the draw. The constraint never sees a bad candidate.
- **N length draws are constraint-masked.** Lengths whose addition
  would force a non-productive junction simply aren't proposed.
- **P length draws follow the same discipline.** A P length that
  would break productivity drops out of the support.

If the constraint admits *no* candidates at a given site,
GenAIRR's default permissive mode falls back to unconstrained
sampling and the run continues (the record may end up
non-productive at that site). Pass `strict=True` to
`run_records(...)` to surface the empty-support condition as
`StrictSamplingError` instead - useful during cartridge
development when you want to catch unsatisfiable plans loudly:

```python
result = exp.run_records(n=100, seed=0, strict=True)
```

## Estimating models from your data

If you have an AIRR-format dataset, the cartridge builder can
estimate all three planes for you:

```python
import csv

rearrangements = list(csv.DictReader(open("rearrangements.tsv"), delimiter="\t"))

cfg = (
    ga.ReferenceCartridgeBuilder
    .from_fasta(v_fasta="v.fa", d_fasta="d.fa", j_fasta="j.fa", chain_type="BCR_HEAVY")
    .infer_identity(species="human", locus="IGH", reference_set="custom", name="my_igh")
    .estimate_np_length_distributions(rearrangements)
    .estimate_np_base_model(rearrangements, kind="markov", pseudocount=0.5)
    .estimate_p_nucleotide_lengths(rearrangements)
    .build()
)
```

Each estimator reads only the AIRR fields it needs:

- **`estimate_np_length_distributions`** consumes `np1_length` /
  `np2_length` directly. GenAIRR's own records always carry
  these; external AIRR tools usually do too.
- **`estimate_np_base_model`** consumes the `np1` / `np2` *string*
  fields (not derived from length arithmetic). With `kind="markov"`,
  pass a non-zero `pseudocount` if the input is small - otherwise
  unobserved transitions leave rows the spec validator will
  reject.
- **`estimate_p_nucleotide_lengths`** consumes `p_v_3_length` /
  `p_d_5_length` / `p_d_3_length` / `p_j_5_length`. **External
  AIRR tools (IgBLAST, MiXCR, ...) do NOT populate these fields.**
  Run this estimator only against GenAIRR's own records, against a
  dataset where you populated the P-length columns yourself, or
  expect the estimator to warn (it emits a stage-level warning
  when ≥95% of contributing rows reported zero for a given key).

The build report on the resulting cartridge captures every
estimator's inputs and inferred distributions for provenance.

## Common mistakes

A handful of issues that show up repeatedly with the junction
addition surfaces.

**Expecting legacy `NP_transitions` / `NP_first_bases` to drive
Markov base sampling.** They don't. The engine no longer auto-
lifts the legacy dicts; if you want Markov, populate
`np_bases["NP1"]` (and / or `NP2`) with an explicit
`NpBaseModelSpec(kind="markov", ...)`. The manifest's
`legacy_np_*_present` flags tell you what's sitting unused on the
bundled cartridges.

**Confusing P and N bases.** Both end up in the junction, but they
come from different planes (`p_nucleotide_lengths` vs
`np_lengths` + `np_bases`) and surface on different AIRR fields
(`p_*_length` vs `np1` / `np2` / `np*_length`). P bases are
deterministic; N bases are sampled. If your downstream pipeline is
mis-counting one as the other, check which AIRR fields it's
reading.

**Expecting `estimate_p_nucleotide_lengths` to infer from generic
AIRR.** External AIRR tools don't model P-nucleotides, so their
output either omits the `p_*_length` columns or populates them as
zero. The estimator will produce a degenerate `[(0, 1.0)]`
distribution and emit the ≥95%-zero warning. Use this estimator
only against data that genuinely carries P-length annotations.

**Expecting `NP2` on a VJ cartridge.** VJ chains have only one NP
region. Attempting to author `np_lengths["NP2"]` or
`p_nucleotide_lengths["D_5"]` / `["D_3"]` on a VJ cartridge
raises `ValueError` at validate time; the engine never gets the
chance to silently lose your spec.

## Where to go next

- **[The Experiment builder](experiment-builder.md)** - how the
  recombination pipeline composes with corruption and clonal
  expansion.
- **[Reference cartridge](../concepts/reference-cartridge.md)** -
  the four typed planes (identity, catalogue, rules, empirical
  models) the engine reads.
- **[SHM and mutation targeting](shm-targeting.md)** - biology-stage
  mutation after junction biology is complete.
- **[Validate AIRR records](../validation/validate-records.md)** -
  what the validator checks about NP / P fields specifically.
