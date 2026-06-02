# Reference Cartridge

A **reference cartridge** is the sealed, typed description of the
biological universe a simulation runs against. Every record GenAIRR
produces is attributable to one cartridge — its identity, its rules,
and its catalogue all participate in the simulation's content hash
and trace metadata.

This document describes the cartridge model, the planes that make it
up, and how the user-facing authoring surface
(`DataConfig`, `Experiment`) plugs into the engine-native
`RefDataConfig`.

---

## What is a reference cartridge?

In the previous generation of GenAIRR, "reference data" meant a
loosely-typed `DataConfig` (a Python pickle) plus some implicit
biology — anchor expectations were inferred from allele names,
"is this a heavy chain?" was deduced from `metadata.chain_type`,
and the engine carried hardcoded knowledge about which catalogue
shapes were biologically sensible.

The cartridge model replaces that with **four explicit planes** plus
one orthogonal concept:

| Plane | Question | Lives on |
| ----- | -------- | -------- |
| **identity** | What cartridge is this? | `RefDataConfig.identity` |
| **catalogue** | Which alleles exist? | `RefDataConfig.v_pool` / `d_pool` / `j_pool` / `c_pool` |
| **rules** | How should the engine interpret alleles? | `RefDataConfig.rules` / `DataConfig.reference_rules` |
| **empirical models** | What are the default sampling distributions? | `DataConfig.reference_models` (Python only in v1) |

**Curation** is orthogonal: given a cartridge, *which subset of the
catalogue participates in simulation?* See
[Curation policy](#curation-policy) below.

---

## The four planes

### Identity

```python
refdata.identity()
# {
#   "species": "Human",
#   "locus": "IGH",
#   "reference_set": "OGRDB V8",
#   "name": "Heavy Chain",
#   "source": "DataConfig",
# }
```

Every field is `Optional[str]`. Synthetic test cartridges may carry
none of it; bundled `Experiment.on("human_igh")` data has all five
fields populated by the loader.

Identity participates in the cartridge content hash. Two cartridges
with identical catalogues + rules but different declared identity
hash differently — trace files attribute outputs to the cartridge,
not just to the catalogue.

The validator also cross-checks `identity.locus` against
`chain_type`: IGH / TRB / TRD must be VDJ; IGK / IGL / TRA / TRG
must be VJ. A mismatch raises `LocusChainTypeMismatch` (Fatal —
never opt-outable).

### Catalogue

The set of V / D / J / C alleles available to the simulator. Each
allele carries `name`, `gene`, `seq`, `segment`, and an optional
`anchor` position (V Cys / J W-or-F codon offset). Population is
exactly what `cfg.add_v_allele(...)` etc. have always done — the
catalogue concept is unchanged from earlier GenAIRR versions.

**V-region substructure annotations.** V alleles can additionally
carry an `Allele.subregions` mapping from IMGT region label to
`(start, end)` ungapped nucleotide intervals — the five canonical
labels `FWR1` / `CDR1` / `FWR2` / `CDR2` / `FWR3`. When a V allele
has a populated `gapped_seq` (the IMGT-aligned sequence with `.`
characters at gap positions), the bridge derives these intervals
automatically via the `compute_v_region_boundaries` helper. A user
can override the dict explicitly; the bridge validates the
overrides at load time (canonical labels, no duplicates,
`start < end`, in-bounds, no overlap) and rejects malformed input
with a clear error. Subregions cross the bridge as
`Vec<VSubregion>` on the Rust `Allele`, fold into the cartridge
content hash, and surface coverage statistics in the manifest's
`models.shm.v_subregion_support` block. **Today they are
inspectable + hashed only — they do NOT yet drive SHM targeting
(no `v_subregion_rates`) or appear in per-region AIRR counters
(no `n_cdr1_mutations` etc.).** Those are separate future slices.
See [`v_region_substructure_audit.md`](v_region_substructure_audit.md).

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
```

Severity strings are `"fatal"` or `"curatable"` — the same vocabulary
the validator's issue dicts use.

When a `ReferenceRulesSpec` is attached to a `DataConfig`, the
loader ships it verbatim into `RefDataConfig.rules` and the
engine's bundled-locus inference (IGH → W, IGK/IGL/TR* → F) is
overridden.

### Empirical models

The **defaults plane** — distributions the engine samples from when
the user doesn't override at recombine time.

```python
from GenAIRR import ReferenceEmpiricalModels, EmpiricalDistributionSpec

models = ReferenceEmpiricalModels(
    np_lengths={
        "NP1": EmpiricalDistributionSpec([(0, 1.0), (3, 4.0), (6, 2.0)]),
    },
    trims={
        "V_3": EmpiricalDistributionSpec([(0, 5.0), (1, 3.0), (2, 1.0)]),
        "J_5": EmpiricalDistributionSpec([(0, 5.0), (1, 2.0)]),
    },
)
```

In v1 this plane lives only on the Python side. The Rust engine
consumes already-lowered pass distributions; the Python loader
preferred order is:

1. `cfg.reference_models[...]` (typed, validated).
2. `cfg.NP_lengths` / `cfg.trim_dicts` (legacy nested-dict
   extraction, marginalised across genes).
3. Uniform `[(0, 1.0), ..., (6, 1.0)]` placeholder
   (raw-`RefDataConfig` fallback).

Distributions are validated shape-only Python-side: non-negative
int values, finite positive weights, known keys, no D-trim on a VJ
cartridge.

---

## Rust `RefDataConfig` vs Python `DataConfig`

| | `DataConfig` (Python) | `RefDataConfig` (Rust) |
| - | - | - |
| Identity | `metadata.species`, `metadata.chain_type`, `metadata.reference_set`, `name` | `identity.species`, `identity.locus`, `identity.reference_set`, `identity.name`, `identity.source` |
| Catalogue | `v_alleles`, `d_alleles`, `j_alleles` (nested `{gene: [Allele]}`) | `v_pool`, `d_pool`, `j_pool`, `c_pool` (flat `AllelePool`) |
| Rules | `reference_rules: Optional[ReferenceRulesSpec]` | `rules: ReferenceRules` |
| Empirical models | `reference_models: Optional[ReferenceEmpiricalModels]` | (not yet — engine consumes lowered pass distributions) |
| Validation | shape-only on the specs | full `validate()` / `validate_with_mode()` |
| Curation | (none — see Rust side) | `curated("functional_anchors_only")` |

The bridge is `dataconfig_to_refdata(cfg)`: alleles + identity +
rules cross the PyO3 boundary; empirical models stay Python-side
and feed `Experiment.recombine`'s default extraction.

---

## Curation policy

Curation is **not validation** and **not catalogue authoring**. It's
the answer to: *given a catalogue and a set of rules, which subset
of alleles do we let the engine sample from?*

Two policies in v1:

- `"raw"` — identity (no filtering).
- `"functional_anchors_only"` — drop V/J alleles that fail the
  active `AnchorRule`: missing anchor, anchor out of bounds, or
  anchor codon AA outside `expected_amino_acids`. D and C pools
  pass through unchanged.

```python
# Strict compile passes after curation; raw catalogue would fail.
ga.Experiment.on("mouse_igh") \
    .curate_refdata("functional_anchors_only") \
    .recombine() \
    .compile()
```

Curation **never silences structural corruption** (duplicate names,
invalid bytes, locus/chain mismatch) — those continue to surface
from the validator on the curated cartridge. If curation empties a
required pool, compile fails with `EmptyRequiredPool`.

The curated cartridge's `identity.source` is tagged
`|curated:functional_anchors_only` so trace files distinguish raw
from curated artefacts and content hashes diverge.

### Curation vs `allow_curatable_refdata`

Two ways to handle pseudogene-bearing catalogues:

- `curate_refdata("functional_anchors_only")` **removes** the
  non-canonical alleles. Cartridge is clean; strict validation
  passes by construction. This is the professional model — the
  cartridge identifies which alleles participate.
- `allow_curatable_refdata()` **keeps** the catalogue as-is and
  relaxes the validator at compile time. Strict validation passes
  Curatable issues (pseudogene-shape anchor anomalies) but still
  rejects Fatal ones (empty pools, duplicates, invalid bytes,
  locus/chain mismatch).

Recommended progression: start strict; if you want to *exclude*
pseudogenes, use `curate_refdata`; if you want to *sample from*
them explicitly, use `allow_curatable_refdata`.

---

## Three cartridge creation paths

A `DataConfig` can reach the simulator from three sources, each
with a distinct provenance signal. Pick the path that matches
how you arrived at the cartridge:

| Path | Use case | `build_report` | `schema_sha256` |
|---|---|---|---|
| **Bundled cartridge** (`ga.HUMAN_IGH_OGRDB`, …) — 106 ship with the wheel via lazy `__getattr__` | Production users with a supported species/locus combo | `None` (predate the builder; legacy provenance) | populated at wheel-build time |
| **Manual `DataConfig(...)`** — direct dataclass construction + field-by-field population (see [Minimal custom cartridge](#minimal-custom-cartridge) below) | Tests, lightweight programmatic creation, edge-case authoring | `None` (no audit-trail layer) | NOT populated by default; call `cfg.schema_sha256 = cfg.compute_checksum()` before `verify_integrity()` |
| **`ReferenceCartridgeBuilder`** — staged FASTA-driven authoring with a build report | New cartridges from raw FASTA where you want a reproducibility/audit record | populated with per-stage entries + manifest snapshot + checksum | populated automatically by `build()` |

### Pick the builder when you want an audit trail

`ReferenceCartridgeBuilder` is the recommended path when you
need provenance about *how* the cartridge was produced —
which FASTA inputs went in, which stages ran, which alleles
were rejected at parse time and why, what the manifest
looked like at the moment `build()` finalised. The build
report rides through pickle, so a teammate loading your
cartridge can replay the audit trail without re-running the
builder.

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
report = builder.report()
```

`cfg.build_report` is a `CartridgeBuildReport` dataclass; its
`.to_dict()` method returns a JSON-clean dict so CI
artifacts / audit dashboards can capture the full provenance.
v1 statistical estimators (`estimate_allele_usage`,
`estimate_trim_distributions`, `estimate_np_length_distributions`,
`estimate_np_base_model`, `estimate_p_nucleotide_lengths`,
`estimate_shm_rates`) are intentionally deferred to a
follow-up slice; the v1 builder covers FASTA → identity →
V-subregions → rules/models attachment → `build()`. See
[`docs/reference_cartridge_authoring_audit.md`](reference_cartridge_authoring_audit.md)
for the full design.

### When NOT to use the builder

- **Bundled cartridge already covers your case** — just
  `ga.HUMAN_IGH_OGRDB`; the build report stays `None` because
  bundled cartridges predate the builder.
- **One-off test fixture** — manual `DataConfig(name="…")`
  with a hand-stitched catalogue is faster to author and
  doesn't need an audit trail.
- **You need parameters the builder doesn't ship yet** —
  statistical estimators are deferred; attach
  `ReferenceEmpiricalModels` directly via `.with_models(...)`
  if you've estimated rates / NP / P distributions outside
  GenAIRR.

---

## Minimal custom cartridge

Author a `DataConfig` with explicit rules + models and run it
through `Experiment`:

```python
from datetime import date

import GenAIRR as ga
from GenAIRR import (
    AnchorRuleSpec,
    ConfigInfo,
    DataConfig,
    EmpiricalDistributionSpec,
    ReferenceEmpiricalModels,
    ReferenceRulesSpec,
    Species,
    ChainType,
)
from GenAIRR.alleles.allele import VAllele, DAllele, JAllele

# (1) Catalogue — V/D/J alleles
v1 = VAllele("MyV1*01", "TGT...", length=300, anchor_override=288)
d1 = DAllele("MyD1*01", "GGGCCC...", length=15)
j1 = JAllele("MyJ1*01", "TGG...", length=60, anchor_override=10)

# (2) Identity
metadata = ConfigInfo(
    species=Species.HUMAN,
    chain_type=ChainType.BCR_HEAVY,
    reference_set="custom",
    last_updated=date(2026, 1, 1),
    has_d=True,
)

# (3) Rules — defaults are fine for canonical biology
rules = ReferenceRulesSpec(
    v_anchor=AnchorRuleSpec(expected_aa=["C"]),
    j_anchor=AnchorRuleSpec(expected_aa=["W"]),
)

# (4) Empirical models — defaults the engine samples from
models = ReferenceEmpiricalModels(
    np_lengths={
        "NP1": EmpiricalDistributionSpec([(0, 1.0), (3, 4.0), (6, 2.0)]),
        "NP2": EmpiricalDistributionSpec([(0, 1.0), (2, 3.0), (4, 2.0)]),
    },
)

cfg = DataConfig(
    name="my_custom_igh",
    metadata=metadata,
    v_alleles={"MyV1": [v1]},
    d_alleles={"MyD1": [d1]},
    j_alleles={"MyJ1": [j1]},
    gene_use_dict={"V": {}, "D": {}, "J": {}},
    trim_dicts={"V_3": {}, "D_5": {}, "D_3": {}, "J_5": {}},
    reference_rules=rules,
    reference_models=models,
)

result = ga.Experiment.on(cfg).recombine().run_records(n=100, seed=42)
```

---

## Non-standard species example

A custom J anchor expectation (`Y` instead of `W`), an extended
alphabet that admits `R` (purine ambiguity code), and a tailored
NP-length distribution:

```python
from GenAIRR import (
    AnchorRuleSpec,
    EmpiricalDistributionSpec,
    ReferenceEmpiricalModels,
    ReferenceRulesSpec,
)

rules = ReferenceRulesSpec(
    allowed_bases=["A", "C", "G", "T", "N", "R"],
    v_anchor=AnchorRuleSpec(expected_aa=["C"]),
    j_anchor=AnchorRuleSpec(expected_aa=["Y"]),  # non-standard locus
)

models = ReferenceEmpiricalModels(
    np_lengths={
        "NP1": EmpiricalDistributionSpec(
            [(0, 0.1), (5, 0.4), (10, 0.3), (15, 0.2)]
        ),
    },
)

cfg.reference_rules = rules
cfg.reference_models = models
```

The cartridge now declares: *my J alleles' anchor codons translate
to tyrosine, my catalogue may contain purine-ambiguity bases, my NP1
lengths are long-tailed.* No engine-side code needs to change — the
rules and models are data, not behaviour.

---

## How validation, curation, and compile interact

Three layers, in order of when they fire:

1. **Validation** (`RefDataConfig.validate()` /
   `validate_with_mode(mode)`) describes the catalogue against the
   cartridge's rules. Returns a list of issue dicts (each tagged
   with `severity`: `"fatal"` or `"curatable"`).

2. **Curation** (`RefDataConfig.curated(policy)`) selects which
   alleles participate. Re-runs validation on the curated cartridge.
   Fatal structural issues are NOT fixed by curation — they still
   surface.

3. **Compile** (`Experiment.compile()`) runs validation under either
   strict mode (default) or `AllowCuratable` mode (after
   `.allow_curatable_refdata()`). Fatal issues always reject;
   Curatable issues reject only under strict mode.

```python
# Direct authoring:
issues = cfg.validate()
cfg.validate_strict()                          # raises on any issue
cfg.validate_with_mode(mode="allow_curatable") # raises only on Fatal
curated = cfg.curated("functional_anchors_only")

# Through Experiment:
ga.Experiment.on(cfg).curate_refdata("functional_anchors_only").compile()
ga.Experiment.on(cfg).allow_curatable_refdata().compile()
```

If you're authoring a new cartridge, the workflow is usually:

1. Build the catalogue + identity + rules + models.
2. Call `cfg.validate()` (or `Experiment.on(cfg).compile()`) to see
   what the validator says.
3. If it complains about pseudogene-shape anchors and you want them
   excluded: `.curate_refdata("functional_anchors_only")`.
4. If you want them included anyway: `.allow_curatable_refdata()`.
5. If it complains about something Fatal (duplicate names, invalid
   bytes, locus/chain mismatch): fix the cartridge itself; neither
   curation nor opt-in helps.

---

## Why this matters

Before the cartridge model, biology lived implicitly in:

- engine code (locus-prefix inference from allele names),
- DataConfig keys (loose nested dicts), and
- "well, the user just knows" tribal knowledge.

After:

- biology lives **on the cartridge as data** — identity, rules,
  empirical models;
- the engine code that consumes the cartridge is generic and
  data-driven;
- trace files attribute outputs to a precise cartridge identity, so
  reproducibility crosses machines and code versions;
- non-standard species (custom alphabet, non-W/F J anchor,
  empirical distributions you measured) are first-class authoring
  scenarios, not "hack the engine" workarounds.
