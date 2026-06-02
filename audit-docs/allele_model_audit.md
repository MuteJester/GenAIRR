# Allele Model Audit

A read-only architecture review of the Python `Allele` layer that
sits beneath the now-complete reference cartridge (identity / rules
/ catalogue / curation / empirical models).

The cartridge model is sealed at the planes-and-rules level. The
allele *objects* underneath those planes are still a mix of T2-8
provenance metadata, legacy build-time helpers, and orphaned trim
logic. Before any biology / DSL slice runs on top of the cartridge,
this audit documents:

1. which legacy paths still fire,
2. what crosses the Python → Rust bridge,
3. what is silently dropped,
4. what should eventually become part of the cartridge catalogue,
5. how reliable the string-split gene / family parsing is,
6. whether the `anchor_override` provenance gap matters,
7. whether trimming still belongs on allele objects.

The companion test file
[`tests/test_allele_model_audit.py`](../tests/test_allele_model_audit.py)
pins every behavioural claim made here. A future refactor that
changes these behaviours should update both files in lockstep — the
final test in that module asserts every section header below
exists.

---

## Production paths that still use legacy methods

The following Allele methods exist on every subclass:

- `VAllele._find_anchor`, `JAllele._find_anchor` — anchor-codon
  resolution.
- `_get_trim_length`, `get_trimmed` — sample a trim amount and
  return the trimmed sequence.

### `_find_anchor`

**Build-time only.** `_find_anchor` runs inside `Allele.__init__`
when one of the bundled-data builders in `.private/scripts/`
constructs `VAllele("name", gapped_seq, length)` — the constructor
calls `self._find_anchor()` unless `anchor_override=...` was passed.

At simulation time, alleles arrive via *unpickling* of
`data/builtin_dataconfigs/*.pkl`. Python's pickle protocol bypasses
`__init__`, so `_find_anchor` is never invoked on a loaded
bundled `DataConfig`. Anchors live on the pickled instances as
plain integers.

The current `_find_anchor` implementation delegates to the C
resolver in `GenAIRR._native._anchor` (T2-8). That module is built
by a separate process; in a vanilla dev install it may not be
present, and the audit tests deliberately avoid importing it.

### `_get_trim_length` / `get_trimmed`

**Orphaned.** A repo-wide grep for `.get_trimmed(` and
`._get_trim_length(` finds:

- self-dispatch inside `src/GenAIRR/alleles/allele.py`
  (`get_trimmed` calls its own `_get_trim_length`);
- this audit's test file (referencing the names in code, not as a
  caller);
- no other production source file or test invokes them.

The actual simulation trim path is:

```
Experiment.recombine()
  └→ _dataconfig_extract.extract_recombine_defaults(cfg)
       ├→ cfg.reference_models.trims[...] (preferred, typed)
       └→ cfg.trim_dicts[...] (legacy nested dict, marginalised)
  └→ engine PassPlan trim distribution
```

The Allele-resident trim methods predate that cartridge path. They
are dead code as far as the current engine is concerned.

---

## Fields that cross into Rust

Confirmed in
[`src/GenAIRR/_refdata_resolver.py::_push_alleles`](../src/GenAIRR/_refdata_resolver.py):

| Python `Allele` field | Rust `Allele` field | Notes |
| --- | --- | --- |
| `name` | `name` | the user-visible identifier |
| `gene` | `gene` | string-split derivative of `name` |
| `ungapped_seq` (`str`) | `seq` (`bytes`) | ASCII-encoded into the Rust pool |
| `anchor` (`int` or `None`) | `anchor` (`u16` or `None`) | the resolved or override-supplied position |
| `functional_status` (enum or `str` or `None`) | `functional_status` (`Option<FunctionalStatus>`) | normalised to `"functional"` / `"orf"` / `"pseudogene"` / `"unknown"`; unknown labels collapse to `None` |

That is the complete contract today. The audit test
`test_bridge_transfers_expected_fields_per_allele` greps the source
of `_push_alleles` to enforce this invariant in both directions:
each expected read must be present, and none of the dropped-field
reads (below) may sneak in.

**`functional_status` history.** This field was added to the
cross-bridge set by the *functional-status curation v1* slice. It
was previously documented under "Fields that are dropped" — the
slice that promoted it updated this section, the dropped-fields
section, and the audit tests in lockstep, exactly as the
"Candidates for the cartridge catalogue" section originally
intended.

---

## Fields that are dropped

The following Python `Allele` fields exist but never reach Rust:

| Field | Class default | Populated when | Used downstream |
| --- | --- | --- | --- |
| `gapped_seq` | n/a (always set) | always | `utilities/imgt_regions.py` only (codon-rail helpers) |
| `length` | n/a (always set) | always | unused in production after construction |
| `family` | n/a (always set) | always | unused in production after construction |
| `anchor_meta` | `None` | when `_find_anchor` ran (T2-8) | inspectable in Python; no downstream consumer |
| `locus` | `None` | never | nowhere |
| `aliases` | `()` | never | nowhere |
| `species` | `None` | never (bundled data leaves it at default) | nowhere |
| `source` | `None` | never | nowhere |
| `frame` (JAllele only) | not set | when J `_find_anchor` ran | nowhere |

The first three (`gapped_seq` / `length` / `family`) are
*derivable* — `length` from `ungapped_seq.length`, `family` from
`name.split(...)`, `gapped_seq` only matters if a future analysis
needs IMGT-numbered positions.

The remaining six are *aspirational*: the class declares them so
loaders can populate them, but the bundled builders never do.

---

## Candidates for the cartridge catalogue

The cartridge architecture already has identity (species, locus,
source) at the *cartridge level*. Per-allele extensions worth
considering, in priority order:

1. ~~**`functional_status`**~~ — **promoted** in *functional-status
   curation v1*. The Rust `Allele` now carries an optional
   `FunctionalStatus` enum (`Functional` / `Orf` / `Pseudogene` /
   `Unknown`), wired through `cfg.curated("functional_status",
   allowed=…, keep_unannotated=…)` and `Experiment.curate_refdata`.
   Bundled `.pkl` builders still don't populate the field — the
   default `keep_unannotated=True` preserves their behaviour. A
   later data-migration slice can rebuild the bundled catalogues to
   carry IMGT's labels and tighten the default.

2. **`anchor_meta.confidence`** (`HIGH_CANONICAL` / `MEDIUM` /
   `LOW` / `REJECTED`). Would let the validator distinguish
   "missing anchor" from "rejected by resolver", and let
   curation drop the latter while keeping the former under a
   permissive rule.

3. **`aliases`**. Would let the OGRDB paralog graph round-trip
   through Rust for cross-reference lookups (currently any
   `aliases` query is Python-only).

Recommend *against* moving:

- `gapped_seq` — IMGT-numbering is a Python loader concern; the
  Rust catalogue stores positional anchors which already encode
  the only semantically-load-bearing fact.
- `length` — redundant.
- `family` — derivable; if needed by Rust, derive once at the
  bridge instead of carrying both copies.
- `species` / `locus` / `source` — already on `ReferenceIdentity`
  at the cartridge level, where they actually belong (not
  per-allele).

---

## Is gene / family parsing reliable?

The parsing logic in `Allele.__init__`:

```python
self.family = self.name.split("-")[0] if "-" in self.name else self.name.split("*")[0]
self.gene = self.name.split("*")[0]
```

`test_gene_family_split_pins_documented_cases` pins the behaviour
across canonical IMGT, OGRDB, and edge-case names. Summary:

| Name | family | gene |
| --- | --- | --- |
| `IGHV1-2*01` | `IGHV1` | `IGHV1-2` |
| `IGHJ4*02` | `IGHJ4` | `IGHJ4` |
| `IGKV3-20*01` | `IGKV3` | `IGKV3-20` |
| `TRBV20-1*01` | `TRBV20` | `TRBV20-1` |
| `IGHVF1-G3*01` (OGRDB) | `IGHVF1` | `IGHVF1-G3` |
| `IGHV1-2-3*01` (multi-hyphen) | `IGHV1` | `IGHV1-2-3` |
| `IGHV*01` (no hyphen) | `IGHV` | `IGHV` |
| `IGHV1-2` (no allele suffix) | `IGHV1` | `IGHV1-2` |

The split is **reliable for canonical IMGT and OGRDB names** that
populate the bundled catalogues — every bundled allele has at
least one `-` between family and gene, and at least one `*` between
gene and allele.

It is **fragile for synthetic / custom catalogues** in two ways:

- Names without `-` collapse `family == gene`. Acceptable for
  toy fixtures; surprising for real custom species where the
  user expected family extraction to work.
- Multi-hyphen names lose the back half from `family`
  (`IGHV1-2-3*01` → `family="IGHV1"`). This is the documented
  behaviour for OGRDB-style family-numbering and IMGT
  duplication suffixes; the `gene` field preserves the full
  pre-`*` string regardless.

Recommendation: **acceptable as-is for v1**. If a future custom-
species cartridge demands explicit family/gene metadata, add
optional `family=` / `gene=` constructor kwargs that override the
split; don't change the default parsing (it's stable for bundled
data and the existing builders depend on the current values).

---

## Is `anchor_override` provenance acceptable?

Current behaviour pinned by
`test_anchor_override_sets_anchor_and_leaves_meta_none`:

```python
if anchor_override is not None:
    self.anchor = anchor_override   # provenance: None
else:
    self._find_anchor()              # populates anchor_meta
```

After override: `anchor` is set, `anchor_meta` is `None`. Two
alleles with the same `anchor` value can have wildly different
trust levels — one validated by the C resolver (codon confirmed
Cys), one supplied as an opaque integer the user vouches for.

**Acceptable for v1, with a known limit.** The Rust-side rules
validator re-checks the codon at compile time (`VAnchorNotCys` /
`JAnchorUnexpectedAa` issues), so an override pointing at a non-Cys
codon is caught downstream. What's lost is *provenance*: a trace
file or audit can't tell whether the anchor was validated at load
or supplied verbatim.

**Recommended remediation (future slice)**: when `anchor_override`
is supplied, populate `anchor_meta` with a synthetic
`AnchorResult { confidence: EXPLICIT, source: "override" }` so the
override path leaves the same `anchor_meta` shape as the resolver
path. Cheap fix, no API change, closes the provenance gap.

---

## Should trimming still belong on allele objects?

**No.** The trim methods on `Allele` are orphaned, and the
cartridge architecture already provides the right home for trim
distributions:

- `ReferenceEmpiricalModels.trims` for typed, validated,
  cartridge-authored defaults;
- `DataConfig.trim_dicts` for the legacy nested-dict shape (still
  consumed by the extractor as a fallback);
- the engine's `PassPlan` carries the lowered distribution as a
  pass parameter, sampled at simulation time.

Allele-resident trim methods serve no path in any of these layers.
They can be removed in a future cleanup slice with no behavioural
change. The cleanup is non-trivial because:

- `_get_trim_length` and `get_trimmed` are declared `@abstractmethod`
  on the `Allele` ABC. Removing the abstract declarations is
  fine; removing the concrete subclass implementations is fine.
- The build-time anchor resolution path (`_find_anchor`) still
  belongs on Allele subclasses, so the class hierarchy itself
  isn't deletable.

Recommendation: **deprecate now, remove in a follow-up slice.**
Marking them as deprecated docs the intent without breaking any
hypothetical downstream user who imported them from
`GenAIRR.alleles.allele`.

---

## Summary

| Item | Status | Action before biology / DSL work |
| --- | --- | --- |
| `_find_anchor` build-time-only | Working as designed | none |
| `_get_trim_length` / `get_trimmed` | Orphaned | deprecate, schedule removal |
| Rust bridge transfers 5 fields | Stable, pinned by tests (was 4 pre-status slice) | none |
| `functional_status` not used | **Promoted** (functional-status curation v1) | bundled rebuild = future data slice |
| `anchor_meta.confidence` not used | Documented gap | candidate for validator severity refinement |
| `aliases` / `locus` / `species` / `source` not used | Documented gap | low priority |
| `gene` / `family` split parsing | Reliable for IMGT/OGRDB | optional kwargs in future |
| `anchor_override` leaves `anchor_meta=None` | Documented provenance gap | fixable in one short slice |
| Trim on alleles | Legacy / orphaned | deprecate in same slice as the removal |

The audit surfaces no concrete correctness bug. The legacy allele
layer has documented gaps but doesn't compromise the cartridge
architecture or the engine's invariants. Greenlight biology / DSL
work; revisit the deprecation and provenance-fix slices in
parallel.
