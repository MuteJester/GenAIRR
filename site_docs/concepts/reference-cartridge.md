# Reference cartridge

<p class="lead">A reference cartridge is the sealed, typed description
of the biological universe a simulation runs against. Every record
GenAIRR produces is attributable to one cartridge - its identity,
its rules, and its catalogue all participate in the simulation's
content hash and trace metadata.</p>

!!! tip "Your learning path"
    This is the conceptual entry point for the
    **"I want to build a cartridge"** path. Once the four-plane
    model below is clear, the
    [Build a reference cartridge](../guides/build-reference-cartridge.md)
    guide walks you through the practical builder workflow, and
    [Estimate models from data](../guides/estimate-cartridge-models.md)
    covers fitting empirical biology to a dataset.
    [See all paths →](../learn.md)

## What is a reference cartridge?

If you've used the simulator with `Experiment.on("human_igh")`, you
already used a cartridge - the bundled `HUMAN_IGH_OGRDB` cartridge,
loaded for you. The cartridge concept makes that "reference data"
explicit and typed: instead of an opaque pickle that the engine
inspects with hardcoded assumptions, every piece of biology is
declared on one of four planes.

```python
import GenAIRR as ga

cfg = ga.HUMAN_IGH_OGRDB
manifest = cfg.cartridge_manifest()
# JSON-clean inventory of what's on the cartridge - identity, plane keys,
# legacy state, supported model kinds, content-hash participation, ...
```

## The four planes

A cartridge has four typed planes plus one orthogonal concept
(curation, covered below). Each plane answers exactly one question.

| Plane | Question | Lives on |
|---|---|---|
| **Identity** | What cartridge is this? | `cfg.metadata` + `cfg.name` |
| **Catalogue** | Which alleles exist? | `cfg.v_alleles` / `d_alleles` / `j_alleles` |
| **Rules** | How should the engine interpret alleles? | `cfg.reference_rules` |
| **Empirical models** | What are the default sampling distributions? | `cfg.reference_models` |

### Identity

Identity declares *what cartridge this is* - species, locus,
reference set, name, and curation source. Two cartridges with
identical catalogues but different declared identity hash
differently; the trace files for their runs are not confusable.

```python
cfg.metadata.species         # Species.HUMAN
cfg.metadata.chain_type      # ChainType.BCR_HEAVY
cfg.metadata.reference_set   # "OGRDB V8"
cfg.name                     # "human_igh"
```

The cartridge validator cross-checks `chain_type` against the locus
in `identity`: IGH / TRB / TRD must be VDJ; IGK / IGL / TRA / TRG
must be VJ. A mismatch is a *fatal* issue - neither curation nor
opt-in can mask it.

### Catalogue

The set of V / D / J / C alleles available to the simulator. Each
allele carries a name, gene, sequence, segment, and an optional
anchor position (V Cys / J W-or-F codon offset). The catalogue is
read through `cfg.v_alleles` / `cfg.d_alleles` / `cfg.j_alleles`
(dicts keyed by gene) - the cartridge concept is unchanged from
earlier GenAIRR versions.

V alleles can additionally carry **subregion annotations** - the
five canonical IMGT region intervals (`FWR1` / `CDR1` / `FWR2` /
`CDR2` / `FWR3`) that drive per-V-region SHM targeting and per-V-
region mutation counters. When a V allele has a populated IMGT-
gapped sequence, the bridge derives these intervals for you;
otherwise you can supply them explicitly. See
[Build a reference cartridge](../guides/build-reference-cartridge.md).

### Rules

The **interpretation layer**. Tells the engine how anchor codons,
sequence alphabets, and severity classifications should be read.

```python
from GenAIRR import ReferenceRulesSpec, AnchorRuleSpec

rules = ReferenceRulesSpec(
    allowed_bases=["A", "C", "G", "T", "N"],
    v_anchor=AnchorRuleSpec(expected_aa=["C"], required=True),
    j_anchor=AnchorRuleSpec(expected_aa=["W"], required=True),
)

cfg.reference_rules = rules
```

Severity strings are `"fatal"` or `"curatable"` - the same
vocabulary the cartridge validator uses. When a `ReferenceRulesSpec`
is attached to a cartridge, it overrides the engine's bundled-locus
defaults (IGH → W anchor, IGK/IGL/TR* → F anchor). Non-canonical
species with custom J anchor amino acids or extended sequence
alphabets are first-class authoring scenarios - not "hack the
engine" workarounds.

### Empirical models

The **defaults plane** - distributions the engine samples from when
the user doesn't override at recombine time. v1 carries five typed
sub-planes plus the SHM-targeting kwargs that ride on
`Experiment.mutate`:

```python
from GenAIRR import (
    AlleleUsageSpec,
    EmpiricalDistributionSpec,
    NpBaseModelSpec,
    ReferenceEmpiricalModels,
)

cfg.reference_models = ReferenceEmpiricalModels(
    allele_usage=AlleleUsageSpec(
        v={"IGHVF1-G1*01": 100.0, "IGHVF1-G2*01": 60.0},
    ),
    trims={
        "V_3": EmpiricalDistributionSpec([(0, 5.0), (1, 3.0), (2, 1.0)]),
        "J_5": EmpiricalDistributionSpec([(0, 5.0), (1, 2.0)]),
    },
    np_lengths={
        "NP1": EmpiricalDistributionSpec([(0, 1.0), (3, 4.0), (6, 2.0)]),
    },
    np_bases={
        "NP1": NpBaseModelSpec(
            kind="markov",
            first_base={"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
            transitions={...},
        ),
    },
    p_nucleotide_lengths={
        "V_3": EmpiricalDistributionSpec([(0, 0.8), (1, 0.2)]),
    },
)
```

| Sub-plane | What it drives | Keys |
|---|---|---|
| `allele_usage` | Per-segment allele sampling weights | `v`, `d`, `j` |
| `trims` | Recombination-stage exonuclease trim distributions | `V_3`, `D_5`, `D_3`, `J_5` |
| `np_lengths` | NP1 / NP2 region length distributions | `NP1`, `NP2` |
| `np_bases` | NP-region base sampling - uniform / empirical / Markov | `NP1`, `NP2` |
| `p_nucleotide_lengths` | Templated P-nucleotide length distributions | `V_3`, `D_5`, `D_3`, `J_5` |

**SHM targeting** is a per-experiment surface, not a cartridge
plane: `Experiment.mutate(segment_rates=..., v_subregion_rates=...)`
takes per-V/D/J/NP bucket rates and per-IMGT-subregion rates,
folds them into the plan signature, and is not part of the
cartridge identity. See [Targeted SHM rates](../guides/shm-targeting.md).

Empirical-model precedence at run time is: explicit
`recombine(...)` kwarg first, then the typed cartridge plane, then
the legacy nested-dict fallback (`cfg.NP_lengths` /
`cfg.trim_dicts`), then a uniform placeholder. So a cartridge can
ship typed defaults that users override per-experiment.

## Build a cartridge from FASTA

For new cartridges with an audit trail, use `ReferenceCartridgeBuilder`.
The page that follows is a concept overview; for the practical
end-to-end builder workflow (FASTA in → validated cartridge out, with
rules / models / estimators / curation / common-mistakes) see
[**Build a reference cartridge**](../guides/build-reference-cartridge.md).

```python
import GenAIRR as ga

builder = (
    ga.ReferenceCartridgeBuilder
    .from_fasta(
        v_fasta="v.fa",
        d_fasta="d.fa",
        j_fasta="j.fa",
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
report = builder.report()       # CartridgeBuildReport - pickleable, JSON-clean
print(report.to_dict())          # Capture for CI artefacts
```

`build()` finalises the cartridge, stamps the canonical
`schema_sha256` checksum, attaches the typed `CartridgeBuildReport`
audit trail, runs `verify_integrity()`, and returns a plain
`DataConfig` you can drop straight into `ga.Experiment.on(cfg)`.

### Estimate biology from your own AIRR data

If you have a representative AIRR-format dataset, the builder can
estimate empirical models for you instead of authoring them by
hand:

```python
import csv

rearrangements = list(csv.DictReader(open("rearrangements.tsv"), delimiter="\t"))

cfg = (
    ga.ReferenceCartridgeBuilder
    .from_fasta(v_fasta="v.fa", d_fasta="d.fa", j_fasta="j.fa", chain_type="BCR_HEAVY")
    .infer_identity(species="human", locus="IGH", reference_set="custom", name="my_igh")
    .estimate_allele_usage(rearrangements, ambiguous="fractional")
    .estimate_trim_distributions(rearrangements)
    .estimate_np_length_distributions(rearrangements)
    .estimate_np_base_model(rearrangements, kind="markov", pseudocount=0.5)
    .estimate_p_nucleotide_lengths(rearrangements)
    .build()
)
```

Each estimator reads the AIRR fields it needs (`v_call` for allele
usage, `*_trim_*` for trims, `np1_length` / `np2_length` for NP
lengths, `np1` / `np2` for base composition, `p_*_length` for
P-nucleotide lengths) and writes the result into the corresponding
typed empirical-model sub-plane. The build report's
`stages` list captures every estimator's inputs and inferred
distributions for provenance. See
[Estimate models from data](../guides/estimate-cartridge-models.md) for
per-estimator detail.

## The manifest

`cfg.cartridge_manifest()` returns a JSON-clean inventory of
everything on the cartridge - identity, catalogue counts, rules,
typed-plane keys, legacy state, model-kind support, and
content-hash participation. It's the canonical surface for
introspection:

```python
manifest = cfg.cartridge_manifest()

manifest["identity"]                       # name, species, locus, reference_set, source
manifest["catalogue"]["v_allele_count"]    # int
manifest["rules"]                          # anchor expectations, allowed bases
manifest["models"]["allele_usage"]         # available / nonempty_segments / soft-gap flag
manifest["models"]["trim_models"]          # keys / source / plan-signature flag
manifest["models"]["np_length_models"]     # keys / source / plan-signature flag
manifest["models"]["np_base_models"]       # per-key kind + supported kinds
manifest["models"]["p_nucleotide_models"]  # length_keys / supported_ends
manifest["models"]["shm"]                  # v_subregion_support coverage
```

CI dashboards, cartridge-diff tools, and replay attribution all
read the manifest. Two cartridges whose manifests differ on any
field that participates in the content hash will produce different
trace files for the same seed. The focused
[Inspect manifest + build report](../guides/cartridge-manifest-report.md)
guide is the practical deep dive on what each key means and how
to pin invariants in CI.

## Curation policy

Curation is *not validation* and *not catalogue authoring*. It's
the answer to a third question: *given a catalogue and a set of
rules, which subset of alleles do we let the engine sample from?*

Two policies in v1:

- `"raw"` - identity (no filtering).
- `"functional_anchors_only"` - drop V/J alleles that fail the
  active `AnchorRule`: missing anchor, anchor out of bounds, or
  anchor codon AA outside the expected set. D and C pools pass
  through unchanged.

```python
ga.Experiment.on("mouse_igh") \
    .curate_refdata("functional_anchors_only") \
    .recombine() \
    .compile()
```

Curation **never** silences structural corruption (duplicate
names, invalid bytes, locus/chain mismatch) - those continue to
surface from the validator on the curated cartridge.

### Curation vs `allow_curatable_refdata`

Two ways to handle pseudogene-bearing catalogues:

- `curate_refdata("functional_anchors_only")` **removes** the
  non-canonical alleles. The cartridge is clean; strict validation
  passes by construction. This is the professional model - the
  cartridge identifies which alleles participate.
- `allow_curatable_refdata()` **keeps** the catalogue as-is and
  relaxes the validator at compile time. Strict validation passes
  Curatable issues (pseudogene-shape anchor anomalies) but still
  rejects Fatal ones (empty pools, duplicates, invalid bytes,
  locus/chain mismatch).

Recommended progression: start strict; if you want to *exclude*
pseudogenes use `curate_refdata`; if you want to *sample from* them
explicitly use `allow_curatable_refdata`.

## How validation, curation, and compile interact

Three layers, in the order they fire:

1. **Validation** describes the catalogue against the cartridge's
   rules. Returns a list of issue dicts each tagged with severity
   (`"fatal"` or `"curatable"`).
2. **Curation** selects which alleles participate. Re-runs validation
   on the curated cartridge. Fatal structural issues are NOT fixed
   by curation - they still surface.
3. **Compile** runs validation under either strict mode (default)
   or `AllowCuratable` mode (after `.allow_curatable_refdata()`).
   Fatal issues always reject; Curatable issues reject only under
   strict mode.

```python
# Direct authoring
issues = cfg.validate()
cfg.validate_strict()                          # raises on any issue
cfg.validate_with_mode(mode="allow_curatable") # raises only on Fatal

# Through Experiment
ga.Experiment.on(cfg).curate_refdata("functional_anchors_only").compile()
ga.Experiment.on(cfg).allow_curatable_refdata().compile()
```

When authoring a new cartridge, the practical workflow is: build the
catalogue + identity + rules + models, run `cfg.validate()` to see
what the validator says, then either curate or opt-in depending on
what you want for pseudogenes - and fix anything Fatal at the
cartridge level (neither curation nor opt-in helps with those).

## Common mistakes

A handful of issues that show up repeatedly when authoring custom
cartridges. None are subtle once you know to look for them.

**Missing anchors.** A V allele without an `anchor` (Cys codon
offset) or a J allele without an `anchor` (W/F codon offset) will
surface as a Fatal validation issue when you try to compile.
`ReferenceCartridgeBuilder.from_fasta` populates anchors when the
native resolver is available; without it, you'll see warnings in
the build report and need to supply anchors via `anchor_override`
at allele construction time.

**Wrong locus / chain combination.** Declaring `metadata.chain_type =
ChainType.BCR_HEAVY` while populating the catalogue with kappa
alleles produces a `LocusChainTypeMismatch` Fatal issue. The
validator catches it; the engine never gets the chance to misbehave.

**Unannotated V subregions.** Passing `v_subregion_rates={...}` to
`Experiment.mutate` against a cartridge whose V alleles lack
`subregions` annotations raises at builder time. Either populate
the `gapped_seq` field on each V allele so the bridge derives
subregions automatically, or supply the `(start, end)` intervals
explicitly.

**Expecting legacy fields to auto-lift.** The legacy
`DataConfig.NP_lengths` / `trim_dicts` / `gene_use_dict` /
`NP_transitions` / `NP_first_bases` / `p_nucleotide_length_probs`
dicts continue to drive the engine when no typed plane is
authored - but they are NOT lifted into the typed
`ReferenceEmpiricalModels` planes automatically. If you populate
the typed plane, the typed plane wins and the legacy dict is
ignored. The manifest's `legacy_*_present` flags tell you which
legacy fields are sitting on the cartridge unused.

---

## Deep architecture notes

For the original design audit and the engine-side
`RefDataConfig` bridge, see
[`docs/reference_cartridge.md`](https://github.com/MuteJester/GenAIRR/blob/master/docs/reference_cartridge.md)
and
[`docs/reference_cartridge_authoring_audit.md`](https://github.com/MuteJester/GenAIRR/blob/master/docs/reference_cartridge_authoring_audit.md).
For per-estimator audits (allele usage, trims, NP length, NP base,
P-nucleotide length), see the matching design docs in `docs/`. For
the full engine-wide validation matrix, see
[`docs/validation_matrix.md`](https://github.com/MuteJester/GenAIRR/blob/master/docs/validation_matrix.md).
