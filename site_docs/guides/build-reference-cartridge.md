# Build a reference cartridge

<p class="lead">When the 106 bundled cartridges don't match what you
need - a non-human locus, a curated allele subset, lab-specific
priors, or a simulation calibrated to a real dataset -
<code>ReferenceCartridgeBuilder</code> is the practical path. This
guide walks the builder pipeline end to end: FASTA inputs in,
validated cartridge + build report out.</p>

## When to build your own cartridge

Reach for the builder when any of these apply:

- **Non-human or non-standard reference.** Your species isn't in
  the bundled catalogue, or you have a custom IMGT pull, or your
  J anchors use a non-canonical amino acid.
- **Curated allele set.** You want pseudogenes excluded, or only
  a subset of alleles available for benchmarking, or you've
  hand-curated a list of known-functional alleles.
- **Lab-specific priors.** Your sequencer / library prep produces
  characteristic trim distributions or NP biases you want
  simulated faithfully.
- **Simulation benchmark matching a dataset.** You have a real
  AIRR dataset and want a simulator whose biology mirrors it -
  estimated allele usage, trims, NP lengths, and base composition
  drawn from the data itself.

For everything else, `ga.Experiment.on("human_igh")` (or one of the
106 other bundled shortcuts) is faster and lower-risk.

## Minimal FASTA build

The shortest path that produces a working cartridge:

```python
import GenAIRR as ga

builder = (
    ga.ReferenceCartridgeBuilder
    .from_fasta(
        v_fasta="V.fa",
        d_fasta="D.fa",
        j_fasta="J.fa",
        chain_type="BCR_HEAVY",
    )
    .infer_identity(
        species="human",
        locus="IGH",
        reference_set="custom",
        name="my_igh",
    )
    .infer_v_subregions()
)
cfg = builder.build()

# Drop straight into Experiment.
result = ga.Experiment.on(cfg).recombine().run_records(n=100, seed=0)
```

What each call does:

| Call | Effect |
|---|---|
| `from_fasta(v_fasta=..., d_fasta=..., j_fasta=..., chain_type=...)` | Classmethod constructor - parses V/J (and D on VDJ; C on optional) FASTA into per-segment allele buckets. Accepts file paths, open text handles, or raw FASTA strings. |
| `.infer_identity(species=..., locus=..., reference_set=..., name=...)` | Populates the cartridge's identity plane and tags every parsed allele with provenance. |
| `.infer_v_subregions()` | Derives IMGT `FWR1` / `CDR1` / `FWR2` / `CDR2` / `FWR3` intervals from gapped V sequences. Idempotent - calling twice re-derives. Required if you want `v_subregion_rates`-driven SHM targeting later. |
| `.build()` | Assembles a validated `DataConfig`, stamps the checksum, attaches the build report + manifest snapshot, runs `verify_integrity()` as the final gate. |

`from_fasta` is the only constructor. The builder is fluent -
every step returns the same builder so calls chain.

## Inspect the build report

Every call appends a structured entry to the build report so you
can audit exactly how the cartridge was constructed:

```python
report = builder.report()        # CartridgeBuildReport dataclass
d = report.to_dict()             # JSON-clean dict - pickleable, CI-artifact-safe
print(d.keys())
# dict_keys(['stages', 'warnings', 'rejected', 'manifest_snapshot', 'checksum_at_build_time'])
```

Five top-level fields:

| Field | Shape | What's in it |
|---|---|---|
| `stages` | `list[dict]` | One entry per builder call, each with `{stage, inputs, inferred, warnings}` |
| `warnings` | `list[str]` | Build-finalisation warnings that aren't tied to a single stage |
| `rejected` | `list[dict]` | Per-allele FASTA-parse drops AND per-row estimator rejections, mixed; filter by `stage` |
| `manifest_snapshot` | `dict` | The result of `cfg.cartridge_manifest()` captured at `build()` time |
| `checksum_at_build_time` | `str` | The `schema_sha256` value stamped onto the built cartridge |

A typical inspection workflow:

```python
# How many alleles got parsed vs rejected at FASTA time?
parse_stage = next(s for s in report.stages if s["stage"] == "from_fasta")
print(parse_stage["inputs"])
# {'chain_type': 'BCR_HEAVY', 'v_alleles_parsed': 285, 'v_alleles_rejected': 2,
#  'j_alleles_parsed': 13, 'j_alleles_rejected': 0,
#  'd_alleles_parsed': 27, 'd_alleles_rejected': 0}

# Which V alleles got dropped?
parse_drops = [r for r in report.rejected
               if r["stage"] == "from_fasta" and r.get("segment") == "V"]
for r in parse_drops:
    print(r["allele_name"], "-", r["reason"])

# What did each estimator see?
for s in report.stages:
    if s["stage"].startswith("estimate_"):
        print(s["stage"], s["inputs"])
```

The `manifest_snapshot` is a frozen capture of the cartridge's
state at build time - useful for cartridge diffs and CI provenance.

## Validate and compile

`build()` runs the validation gate as its final step, so a
successfully-returned cartridge is already valid. The next gate
runs at compile time, when an `Experiment` instance picks it up:

```python
exp = ga.Experiment.on(cfg).recombine()
compiled = exp.compile()         # raises on validation failure
```

GenAIRR's validation classifies issues into two severities:

- **Fatal** - duplicate allele names, invalid bytes, locus / chain
  mismatch, empty required pools. These always reject; neither
  curation nor opt-in helps. Fix them at the cartridge level.
- **Curatable** - pseudogene-shape anchor anomalies (missing anchor,
  out-of-bounds, AA outside the expected set). These reject under
  strict mode but pass under `allow_curatable_refdata()`.

When you're authoring against a catalogue that contains
pseudogenes, `allow_curatable_refdata()` is the right opt-in for
the compile path:

```python
compiled = (
    ga.Experiment.on(cfg)
    .allow_curatable_refdata()    # let pseudogenes through compile validation
    .recombine()
    .compile()
)
```

For the alternative (drop pseudogenes from the catalogue entirely
rather than opt-in), see [Curation](#curation) below.

## Add reference rules

If your locus uses non-canonical anchor amino acids or accepts an
extended sequence alphabet, attach a `ReferenceRulesSpec` before
`build()`:

```python
from GenAIRR import ReferenceRulesSpec, AnchorRuleSpec

rules = ReferenceRulesSpec(
    allowed_bases=["A", "C", "G", "T", "N", "R"],     # extended alphabet
    v_anchor=AnchorRuleSpec(expected_aa=["C"], required=True),
    j_anchor=AnchorRuleSpec(expected_aa=["Y"], required=True),  # non-standard J
)

cfg = (
    builder
    .with_rules(rules)
    .build()
)
```

`with_rules` validates the spec immediately - if the anchor /
severity / alphabet shape is wrong, you get the error here rather
than at compile time. Severity strings are `"fatal"` or
`"curatable"` (the same vocabulary the cartridge validator uses).

Calling `with_rules` again replaces the previously-attached rules
wholesale; the build report records both calls.

## Add empirical models

If you've authored empirical distributions externally (e.g. fitted
distributions, hand-tuned NP-length distributions, an
`AlleleUsageSpec` derived from a different tool), attach them
through `with_models`:

```python
from GenAIRR import (
    AlleleUsageSpec,
    EmpiricalDistributionSpec,
    NpBaseModelSpec,
    ReferenceEmpiricalModels,
)

models = ReferenceEmpiricalModels(
    allele_usage=AlleleUsageSpec(v={"MyV1*01": 100.0, "MyV2*01": 30.0}),
    trims={
        "V_3": EmpiricalDistributionSpec([(0, 5.0), (1, 3.0), (2, 1.0)]),
        "J_5": EmpiricalDistributionSpec([(0, 5.0), (1, 2.0)]),
    },
    np_lengths={
        "NP1": EmpiricalDistributionSpec([(0, 1.0), (3, 4.0), (6, 2.0)]),
        "NP2": EmpiricalDistributionSpec([(0, 1.0), (2, 3.0), (4, 2.0)]),
    },
    np_bases={
        "NP1": NpBaseModelSpec(
            kind="markov",
            first_base={"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
            transitions={...},
        ),
    },
    p_nucleotide_lengths={
        "V_3": EmpiricalDistributionSpec([(0, 0.7), (1, 0.3)]),
    },
)

cfg = builder.with_models(models).build()
```

The five typed sub-planes - `allele_usage` / `trims` /
`np_lengths` / `np_bases` / `p_nucleotide_lengths` - are validated
at attach time against the cartridge's chain type, so a D-end
key on a VJ cartridge fails immediately.

**SHM targeting is NOT a cartridge concern.** Per-segment SHM
rates (`segment_rates`) and per-V-subregion SHM rates
(`v_subregion_rates`) ride on `Experiment.mutate(...)` instead.
They fold into the plan signature but don't enter cartridge
identity - same cartridge, different SHM regimes per experiment.
See [SHM and mutation targeting](shm-targeting.md).

## Estimate models from AIRR-like data

If you have a representative AIRR-format dataset, the builder can
estimate every typed model plane for you. Five estimators have
shipped - and the focused **[Estimate cartridge models from real
data](estimate-cartridge-models.md)** guide is the deep dive on
the workflow, ambiguity handling, rejection reasons, and common
mistakes. The condensed reference below covers the same surface
at a glance:

```python
import csv

rearrangements = list(csv.DictReader(open("rearrangements.tsv"), delimiter="\t"))

cfg = (
    ga.ReferenceCartridgeBuilder
    .from_fasta(v_fasta="v.fa", d_fasta="d.fa", j_fasta="j.fa", chain_type="BCR_HEAVY")
    .infer_identity(species="human", locus="IGH", reference_set="custom", name="my_igh")
    .infer_v_subregions()
    .estimate_allele_usage(rearrangements, ambiguous="fractional")
    .estimate_trim_distributions(rearrangements)
    .estimate_np_length_distributions(rearrangements)
    .estimate_np_base_model(rearrangements, kind="markov", pseudocount=0.5)
    .estimate_p_nucleotide_lengths(rearrangements)
    .build()
)
```

Each estimator reads the AIRR fields it needs and writes into the
corresponding `reference_models.*` sub-plane:

| Method | AIRR columns consumed | Writes into |
|---|---|---|
| `estimate_allele_usage(rearrangements, *, min_count=1.0, ambiguous="fractional", replace=True)` | `v_call`, `d_call`, `j_call` | `reference_models.allele_usage` |
| `estimate_trim_distributions(rearrangements, *, min_count=1, pseudocount=0.0, replace=True)` | `v_trim_3`, `d_trim_5`, `d_trim_3`, `j_trim_5` | `reference_models.trims` |
| `estimate_np_length_distributions(rearrangements, *, min_count=1, pseudocount=0.0, replace=True)` | `np1_length`, `np2_length` (VDJ) | `reference_models.np_lengths` |
| `estimate_np_base_model(rearrangements, *, kind="markov", min_count=1, pseudocount=0.0, replace=True)` | `np1`, `np2` (VDJ) | `reference_models.np_bases` |
| `estimate_p_nucleotide_lengths(rearrangements, *, min_count=1, pseudocount=0.0, replace=True)` | `p_v_3_length`, `p_d_5_length`, `p_d_3_length`, `p_j_5_length` | `reference_models.p_nucleotide_lengths` |

Each estimator is idempotent under the `replace` flag: `replace=True`
(default) overwrites the previous spec and stamps `replaced=True`
on the new stage entry; `replace=False` raises `ValueError` if a
prior estimator (or `with_models`) already populated the same
sub-plane.

A few notes:

- **`estimate_allele_usage`'s `ambiguous` kwarg** picks the
  tie-set policy: `"fractional"` (split the row's credit across
  tie-set entries), `"truth_first"` (use the first call only), or
  `"reject"` (drop ambiguous rows into `report.rejected`).
- **`estimate_np_base_model`'s `kind` kwarg** picks the output
  model - `"markov"` (default; first-base row + 4×4 transitions)
  or `"empirical_first_base"` (position-independent categorical).
  Markov on sparse input needs `pseudocount > 0` to avoid leaving
  unfillable transition rows.
- **`estimate_p_nucleotide_lengths` requires P-length fields in
  the input.** External AIRR tools (IgBLAST, MiXCR, ...) do NOT
  populate `p_*_length`. Run this estimator against GenAIRR's own
  output or a dataset where you populated the P-length columns
  yourself - the estimator emits a stage-level warning when ≥95%
  of contributing rows reported zero.

**`estimate_shm_rates` is intentionally deferred.** SHM rate
estimation (per-segment + per-V-subregion) is the only estimator
from the original audit list still pending; per-experiment SHM
targeting via `Experiment.mutate(segment_rates=...,
v_subregion_rates=...)` is the workaround until it ships.

## Curation

Curation is *which subset of the catalogue participates in
sampling* - orthogonal to the four cartridge planes. Two policies
in v1:

```python
exp = (
    ga.Experiment.on(cfg)
    .curate_refdata("functional_anchors_only")   # drop V/J alleles failing anchor rules
    .recombine()
)
```

- **`"raw"`** - identity. No filtering.
- **`"functional_anchors_only"`** - drop V/J alleles that fail the
  active `AnchorRule` (missing anchor, anchor out of bounds, AA
  outside `expected_aa`). D and C pools pass through unchanged.

**Curation vs `allow_curatable_refdata`:**

- `curate_refdata(...)` *removes* the non-canonical alleles from
  the cartridge before sampling. Strict validation passes by
  construction. This is the professional model - the cartridge
  itself identifies which alleles participate.
- `allow_curatable_refdata()` *keeps* the catalogue intact and
  relaxes the validator at compile time. Curatable issues pass;
  Fatal issues still reject.

Curation never silences Fatal structural corruption (duplicate
names, invalid bytes, locus/chain mismatch) - those continue to
surface from the validator on the curated cartridge.

The manifest is the canonical introspection surface for any
authored cartridge:

```python
manifest = cfg.cartridge_manifest()

manifest["identity"]                          # species, locus, reference_set, name, source
manifest["catalogue"]["v_allele_count"]       # post-curation count
manifest["rules"]["allowed_bases"]            # what the alphabet is
manifest["models"]["allele_usage"]            # available + nonempty segments + soft-gap flag
manifest["models"]["trim_models"]             # keys + source + plan-signature flag
manifest["models"]["np_length_models"]        # keys + source + plan-signature flag
manifest["models"]["np_base_models"]          # per-key kind + supported kinds
manifest["models"]["p_nucleotide_models"]     # length_keys + supported_ends
manifest["models"]["shm"]["v_subregion_support"]  # annotation coverage
manifest["hashes"]["data_config_checksum"]    # canonical schema sha256
manifest.get("errors", [])                    # bridge-failure messages (key is conditional)
```

Two cartridges whose manifests differ on any content-hash-
participating field will produce different trace files for the
same seed - that's how reproducibility crosses machines and code
versions.

## Common mistakes

A handful of issues that show up repeatedly when authoring custom
cartridges.

**No gapped V sequence → no V subregions.** `infer_v_subregions`
derives the IMGT region intervals from each V allele's gapped
(IMGT-aligned) sequence. If your FASTA carries ungapped sequences
only, `v_subregion_support.annotated_v_count` will be 0 and
`v_subregion_rates`-driven SHM targeting will refuse to compile.
Either include the IMGT-gapped sequence in the FASTA headers or
supply the `(start, end)` intervals manually before `build()`.

**Expecting legacy fields to auto-lift.** Bundled cartridges carry
several legacy nested dicts (`NP_lengths`, `trim_dicts`,
`NP_first_bases`, `NP_transitions`, `gene_use_dict`,
`p_nucleotide_length_probs`) that the modern engine no longer
reads. Authoring on your own cartridge requires writing into the
typed planes - the legacy dicts are not auto-lifted. The manifest's
`legacy_*_present` flags tell you which legacy fields are sitting
unused.

**Wrong chain type.** Declaring `chain_type="BCR_HEAVY"` while
populating the catalogue with kappa alleles produces a
`LocusChainTypeMismatch` Fatal issue at validate time. The
validator catches it; the engine never gets the chance to
misbehave. Always make `chain_type` match the alleles you're
loading.

**Missing D FASTA for VDJ.** `from_fasta` requires `d_fasta` when
`chain_type` is one of the VDJ shapes (BCR_HEAVY, TCR_BETA,
TCR_DELTA). Omitting it raises immediately. For VJ shapes
(BCR_LIGHT_KAPPA, BCR_LIGHT_LAMBDA, TCR_ALPHA, TCR_GAMMA), `d_fasta`
must NOT be supplied - passing it on VJ raises the same way.

**Using `estimate_p_nucleotide_lengths` on generic AIRR.** External
AIRR tools don't model P-nucleotides; their output either omits
the `p_*_length` columns or populates them as zero. The estimator
will produce a degenerate `[(0, 1.0)]` distribution and emit the
≥95%-zero warning. Run this estimator only against data that
genuinely carries P-length annotations - typically GenAIRR's own
output, or a dataset where you populated the P-length columns
yourself.

## Where to go next

- **[Inspect manifest + build report](cartridge-manifest-report.md)**
  the focused guide for auditing what's in a cartridge and
  pinning it in CI.
- **[Reference cartridge concept](../concepts/reference-cartridge.md)**
  the four-plane model behind the builder API.
- **[SHM and mutation targeting](shm-targeting.md)** - per-segment
  + per-V-subregion SHM rates that ride on `Experiment.mutate`,
  not the cartridge.
- **[Junction N/P additions](junction-additions.md)** - the
  cartridge planes that drive NP and P biology, plus the
  estimator-specific quirks.
- **[Validate AIRR records](../validation/validate-records.md)** -
  the postcondition validator that gates every cartridge's
  output at run time.
- **[The Experiment builder](experiment-builder.md)** - what to do
  with the cartridge once it's built.
