# Reference Cartridge — Completeness Audit

**Status: audit only.** Pins today's reference-cartridge surface
and the gaps a future "complete the cartridge" slice will need to
close. No implementation is proposed; the deliverable is the
shared vocabulary plus contract pins so any later slice lands
against an audited baseline.

This is the bridge audit between engine architecture and user-
facing reproducibility. The engine now depends on the cartridge
as a programmable reference object — four planes (identity,
catalogue, rules, models) plus curation. The remaining risk is
whether the Python `DataConfig` / `Allele` layer and the Rust
`RefDataConfig` still drift in shape and meaning, and whether
users have enough metadata to reproduce or compare simulations
across runs.

Companion to
[`tests/test_reference_cartridge_completeness_contract.py`](../tests/test_reference_cartridge_completeness_contract.py)
which freezes every claim here as either a `pin_scaffold_*`
(today's behaviour) or a `pin_absence_*` (deferred surface that a
later slice flips).

This audit deliberately starts narrow: it documents what exists
and where the boundaries are. The proposed next slice
(`cartridge_manifest()` exporter) is sketched in §12 but explicitly
not committed to — the audit's job is to fix the contract first.

---

## Existing scaffolding the contract pins

| Surface | Where | What it pins |
|---|---|---|
| `DataConfig` dataclass | [`src/GenAIRR/dataconfig/data_config.py:32`](../src/GenAIRR/dataconfig/data_config.py#L32) | Top-level pickle container; all cartridge fields live here. |
| `DataConfig.compute_checksum()` | [`src/GenAIRR/dataconfig/data_config.py:127`](../src/GenAIRR/dataconfig/data_config.py#L127) | Python-side SHA-256 over `pickle.dumps(self)` with documented surgery. |
| `DataConfig.verify_integrity()` | [`src/GenAIRR/dataconfig/data_config.py:181`](../src/GenAIRR/dataconfig/data_config.py#L181) | Load-time pickle integrity check; raises `DataConfigError` on mismatch. |
| Cartridge plane views | [`src/GenAIRR/dataconfig/cartridge_views.py`](../src/GenAIRR/dataconfig/cartridge_views.py) | Four frozen dataclasses (`CartridgeIdentityView` / `CartridgeCatalogueView` / `CartridgeRulesView` / `CartridgeModelsView`). |
| `Allele` ABC | [`src/GenAIRR/alleles/allele.py:29`](../src/GenAIRR/alleles/allele.py#L29) | Per-allele Python record; many fields beyond what crosses to Rust. |
| `dataconfig_to_refdata` bridge | [`src/GenAIRR/_refdata_resolver.py:210`](../src/GenAIRR/_refdata_resolver.py#L210) | Python→Rust translation; reads only a subset of Allele + DataConfig fields. |
| `RefDataConfig` Rust struct | [`engine_rs/src/refdata.rs:188`](../engine_rs/src/refdata.rs#L188) | Engine-side cartridge with `Allele` records carrying `functional_status`. |
| `RefDataConfig.content_hash()` | [`engine_rs/src/python/refdata.rs:279`](../engine_rs/src/python/refdata.rs#L279) | Rust-side SHA-256 over `chain + identity + rules + pools(name/gene/seg/seq/anchor/functional_status)`. |
| `refdata.curated(policy, …)` | [`engine_rs/src/python/refdata.rs:348`](../engine_rs/src/python/refdata.rs#L348) | Curation entry point; tags `identity.source` with `|curated:<policy>` so the content hash differs from raw. |

---

## 1. Q1 — Which `DataConfig` fields still affect simulation but
##        are not in the cartridge views?

Today's four cartridge views cover:

| View | Surfaces |
|---|---|
| `CartridgeIdentityView` | `name`, `metadata` |
| `CartridgeCatalogueView` | `v_alleles`, `d_alleles`, `j_alleles`, `c_alleles` |
| `CartridgeRulesView` | `reference_rules` |
| `CartridgeModelsView` | `reference_models`, `legacy_np_lengths` (= `NP_lengths`), `legacy_trim_dicts` (= `trim_dicts`) |

### DataConfig fields NOT in any cartridge view

Looking at the full `DataConfig` field list against the views:

| Field | In view? | Simulation effect | Audit note |
|---|---|---|---|
| `name` | identity | — | |
| `metadata` | identity | identity for locus cascade | |
| `v_alleles` / `d_alleles` / `j_alleles` / `c_alleles` | catalogue | — (catalogue is THE simulation input) | C dropped at bridge (engine has no C-segment passes). |
| `gene_use_dict` | **NOT in any view** | Used by recombine sampling weights when no explicit allele weights override it | **Audit gap** — affects sampling but not surfaced as a plane. |
| `trim_dicts` | models (legacy) | Trim distribution defaults when `reference_models` is None | Covered. |
| `NP_transitions` | **NOT in any view** | Per-base NP-content transitions | **Audit gap** — affects pool-byte sampling but not in a view. |
| `NP_first_bases` | **NOT in any view** | First-base distribution for NP regions | **Audit gap** — same shape as NP_transitions. |
| `NP_lengths` | models (legacy) | NP-length distribution defaults | Covered. |
| `correction_maps` | **NOT in any view** | Validator helper for allele-call corrections | **Audit gap** (likely Python-only consumer; verify before pinning as required). |
| `asc_tables` | **NOT in any view** | Same as correction_maps | **Audit gap**. |
| `p_nucleotide_length_probs` | **NOT in any view** | P-nucleotide length distribution | **Audit gap**. |
| `dj_pairing_map` | **NOT in any view** | D-J pairing constraint when set | **Audit gap** — affects recombination admissibility. |
| `schema_version` / `schema_sha256` | none | Pickle integrity machinery | Not a "plane" — meta. |
| `build_report` | none | Diagnostic provenance, removed from checksum | Not a plane — diagnostic. |
| `reference_rules` | rules | Engine interpretation layer | Covered. |
| `reference_models` | models | Typed distribution overrides | Covered. |

### Verdict

**Six "still affects simulation" DataConfig fields have no
cartridge view**: `gene_use_dict`, `NP_transitions`,
`NP_first_bases`, `correction_maps`, `asc_tables`,
`p_nucleotide_length_probs`, `dj_pairing_map`. These are mostly
legacy / Python-internal surfaces — they don't cross to the Rust
`RefDataConfig` (the bridge only reads alleles + metadata +
reference_rules) — but they DO influence the Python-side
recombine sampling path when `reference_models` is None.

This is the audit's primary "completeness" gap: until either
they're absorbed into `CartridgeModelsView` (as a typed
`SamplingDefaults` sub-plane) or formally pinned as deprecated,
their inclusion in `compute_checksum` but absence from the views
means a future contributor could easily change a cartridge's
behaviour without updating any view-level documentation or test.

Pinned by `pin_scaffold_dataconfig_fields_not_in_any_view` —
the contract test enumerates these so an absence becomes a
deliberate documented gap, not an oversight.

---

## 2. Q2 — Which `Allele` fields are dropped at the Python→Rust
##        bridge?

The bridge function `_push_alleles`
([`_refdata_resolver.py:153`](../src/GenAIRR/_refdata_resolver.py#L153))
extracts exactly:

```python
allele.name           # → engine Allele.name
allele.gene           # → engine Allele.gene
allele.ungapped_seq   # → engine Allele.seq (bytes)
allele.anchor         # → engine Allele.anchor (Option<u16>)
allele.functional_status  # → engine Allele.functional_status (Option<FunctionalStatus>)
```

### Promoted (the recent change)

`functional_status` was promoted from "Python-only" to "crosses
to Rust" in the cartridge-curation work — its presence on the
Rust-side `Allele` enables `refdata.curated("functional_status", …)`
to filter pools server-side and contributes to the content hash.

### Still dropped (intentional)

| Allele field | Dropped because | Status |
|---|---|---|
| `aliases` | OGRDB paralog graph; Python-only lookup surface | Intentional — engine doesn't need allele aliases for simulation. |
| `anchor_meta` | T2-8 provenance dataclass (codon/residue/confidence/method/rejection reason); Python-only | Intentional — engine consumes the resolved `anchor: Option<u16>` directly; the resolution metadata is diagnostic. |
| `gapped_seq` | Only `ungapped_seq` crosses (Rust engine works on ungapped bytes) | Intentional — gap dots have no engine-side meaning. |
| `family` | Re-derivable from `name` (e.g. `"IGHV1-2*01"` → family `"IGHV1"`); Rust derives it on construction | Intentional — redundant. |
| `locus` | At the cartridge level via `identity.locus`, not per-allele | Intentional — locus is a cartridge-wide property, not per-allele. |
| `species` | At the cartridge level via `identity.species` | Intentional — same. |
| `source` | At the cartridge level via `identity.source` | Intentional — same. |
| `length` / `ungapped_len` | Re-derivable from `seq.len()` | Intentional — redundant. |

### Verdict

**Seven Allele fields are intentionally dropped** (`aliases`,
`anchor_meta`, `gapped_seq`, `family`, `locus`, `species`,
`source`, plus the two redundant lengths). Three of these
(`locus`, `species`, `source`) are explicitly cartridge-level —
duplicated per-allele in Python but consolidated to identity in
Rust. Two (`aliases`, `anchor_meta`) are Python-only diagnostic
surfaces. The rest are pure redundancy.

Pinned by `pin_scaffold_allele_fields_dropped_at_bridge_documented` —
the contract test enumerates exactly which fields are dropped so a
future contributor adding a new Allele field has to consciously
decide whether it crosses the bridge or stays Python-side.

---

## 3. Q3 — Does `compute_checksum()` align with
##        `refdata.content_hash()`?

The two hashes serve **different roles** and consequently cover
different surfaces:

### `DataConfig.compute_checksum()` — Python-side cartridge integrity

SHA-256 over `pickle.dumps(self, protocol=4)` with this surgery:
- `schema_sha256` zeroed (can't include itself).
- `build_report` removed (diagnostic; legacy compatibility).
- `reference_rules` removed only when `None` (soft-transition).
- `reference_models` removed only when `None` (soft-transition).

This covers **everything** in `__dict__`: gene_use_dict,
NP_transitions, NP_first_bases, NP_lengths, trim_dicts,
correction_maps, asc_tables, p_nucleotide_length_probs,
dj_pairing_map, reference_models, reference_rules, identity
fields, the catalogue, etc. Any change to any of these — even
fields that don't cross to Rust — produces a different checksum.

### `RefDataConfig.content_hash()` — Rust-side cartridge identity

SHA-256 over a canonical text format covering exactly:
- `chain` (vdj / vj).
- Identity: `species`, `locus`, `reference_set`, `name`, `source`.
- Rules: `alphabet` + `v_anchor` + `j_anchor` (each with
  `required`, `expected`, `missing_severity`, `mismatch_severity`).
- Pools (v/d/j/c): for every allele, `name`, `gene`, `segment`,
  `seq`, `anchor`, `functional_status`.

This covers **only what crosses the bridge**. Fields that live
on `DataConfig` but never reach `RefDataConfig` don't contribute.

### Confirmed divergences

| Field swapped | `compute_checksum` changes? | `content_hash` changes? |
|---|---|---|
| `reference_rules` (truly different from defaults) | Yes | Yes |
| `reference_models` (non-None) | **Yes** | **No** ⚠️ |
| `functional_status` on alleles | Yes (via catalogue pickle bytes) | Yes |
| Curation applied via `refdata.curated(...)` | n/a (Python `DataConfig` unchanged) | Yes (identity.source tagged) |
| `gene_use_dict` mutation | Yes | No |
| `NP_transitions` mutation | Yes | No |
| Identity `species` / `locus` / `name` mutation | Yes | Yes (after bridge re-runs) |

### Are the differences intentional?

Mostly yes. `compute_checksum` is the **pickle-integrity** check —
it catches any byte-level mutation of the saved file, including
legacy fields the engine no longer consumes. `content_hash` is
the **trace-attribution identity** — it answers "is this the same
cartridge the trace was recorded against?" Two cartridges with
identical catalogues + identity + rules SHOULD share a content
hash even if their legacy field tails differ.

The one row marked ⚠️ — `reference_models` changes
`compute_checksum` but not `content_hash` — is the user-spec's
"current v1 boundary." Because `reference_models` is consumed by
the Python-side recombine sampler (not the Rust bridge), swapping
it changes simulation output but not the cartridge's Rust-side
identity. The audit pins this **as the v1 contract** without
calling it a bug; a future slice may want to fold
`reference_models` digest into `content_hash` so trace replay
catches model swaps.

Pinned by `pin_scaffold_reference_models_changes_compute_checksum_not_content_hash` and `pin_scaffold_reference_rules_changes_both_hashes`.

---

## 4. Q4 — Are `reference_rules` and `reference_models`
##        consistently part of cartridge identity?

### `reference_rules`

**Yes, fully part of identity** — when set to a value that
materially differs from the locus-derived defaults, both
`compute_checksum` and `content_hash` change. The Rust bridge
transfers it verbatim into `RefDataConfig.rules` and
`content_hash` includes the alphabet + per-anchor rule
specification.

Edge: a `reference_rules` spec that matches the locus-derived
defaults exactly produces the **same** content_hash (because the
Rust-side rules end up identical regardless of how they got
there). This is correct behaviour — identity is over the
*effective* rules, not their provenance.

### `reference_models`

**Partially part of identity** — `compute_checksum` covers them;
`content_hash` does NOT.

**Pin this as the v1 contract.** A `reference_models` swap
changes the Python-side sampling distributions and therefore
changes simulation output, but the Rust engine's `RefDataConfig`
doesn't see them (they're consumed Python-side before the bridge).
Trace replay against the Rust `RefDataConfig` would succeed across
the swap because the recorded choices are the source of truth in
replay — but a fresh run with a different `reference_models` and
the same seed produces different output that the content_hash
wouldn't flag.

**This is the audit's marked-acceptable boundary.** Pinned by
`pin_scaffold_reference_models_changes_compute_checksum_not_content_hash`.
A future slice that wants symmetric identity would either:
(a) digest `reference_models` and fold the digest into the
identity passed to the Rust bridge, or (b) move
`reference_models` consumption to the Rust side (out of scope for
the current architecture).

---

## 5. Q5 — Are curation policies fully reproducible?

Curation runs on the Rust side via `refdata.curated(policy, …)`.
Two policies ship:

| Policy | What it does |
|---|---|
| `"raw"` | Identity — unchanged clone. |
| `"functional_anchors_only"` | Drops V/J alleles whose anchor doesn't satisfy the active rule. |
| `"functional_status"` | Filters V/D/J by IMGT functional status; `allowed` + `keep_unannotated` parameters. |

### Does identity source tag enough?

**Yes.** The curated cartridge's `identity.source` is appended
with `|curated:<policy_tag>`. The tag is canonical:
`functional_status:functional|keep_unannotated=true`,
`functional_anchors_only`, etc. Two cartridges curated identically
produce identical content hashes; two cartridges curated
differently produce different hashes.

Verified: `HUMAN_IGH_OGRDB` curated as
`functional_status(allowed=["functional"], keep_unannotated=True)`
ends up with
`identity.source == "DataConfig|curated:functional_status:functional|keep_unannotated=true"`,
and `content_hash` differs from the raw cartridge.

### Is the allowed functional-status set represented anywhere
inspectable?

`FunctionalStatus` (a Rust enum: `Functional`, `Orf`, `Pseudogene`,
`Unknown`) is the canonical set. It crosses to Python as a
`Optional[str]` accessor (`allele.functional_status()`) on the
PyO3 `Allele` binding. The Python `_KNOWN_STATUSES` constant
([`_refdata_resolver.py:189`](../src/GenAIRR/_refdata_resolver.py#L189))
mirrors the lowercase set: `"functional"`, `"orf"`, `"pseudogene"`,
`"unknown"`.

The two sets agree by construction (Rust normalises to lowercase
in `functional_status_to_str`). No drift today.

Pinned by `pin_scaffold_curation_tags_identity_source` and
`pin_scaffold_functional_status_set_documented`.

---

## 6. Q6 — Do bundled cartridges expose enough metadata for
##        users to reproduce/compare simulations?

For `HUMAN_IGH_OGRDB` after `dataconfig_to_refdata(cfg)`:

```python
refdata.identity() == {
    "species": "Human",
    "locus": "IGH",
    "reference_set": "OGRDB V8",
    "name": "Heavy Chain",
    "source": "DataConfig",
}
```

Plus the rules plane (locus-derived J anchor `["W"]` for IGH;
default V anchor + alphabet); the catalogue (V/D/J pools);
and `content_hash()` for stable identity.

### What's present

- **species / locus / reference_set / name / source** — yes, on
  every bundled cartridge.
- **rules** — yes, derived from locus when no explicit spec.
- **catalogue** — yes, V/D/J pools.
- **curation policy** — yes, via `identity.source`'s
  `|curated:…` suffix when applicable.
- **empirical model provenance** — partial: `reference_models`
  is the typed surface but bundled cartridges currently set it
  to `None` (legacy nested-dict fallback). A user inspecting a
  bundled cartridge sees `cartridge_models.reference_models is
  None` and `legacy_np_lengths` / `legacy_trim_dicts` populated.

### What's missing (the audit's "complete" question)

- **A single inspectable manifest.** Today a user has to call
  `cfg.cartridge_identity`, `cfg.cartridge_rules`,
  `cfg.cartridge_catalogue`, `cfg.cartridge_models`, plus
  `dataconfig_to_refdata(cfg).content_hash()`, plus
  `dataconfig_to_refdata(cfg).identity()` — five separate
  surfaces — to assemble the full picture. No
  `cartridge_manifest()` method exists.
- **Provenance for the empirical-model digest.** Even
  cartridges that DO use `reference_models` produce a manifest
  with no digest of those models (because no
  manifest exists).
- **Counts by functional-status bucket.** A user comparing two
  cartridges can count alleles total but not
  "how many are functional / ORF / pseudogene." The Rust side
  knows; no Python accessor surfaces the histogram.

Pinned by `pin_absence_no_cartridge_manifest_method` and
`pin_absence_no_functional_status_histogram_accessor`.

---

## 7. Q7 — What must be fixed before calling the cartridge
##        model "complete"?

The audit's verdict — three concrete gaps to close, **in this
order**:

### Gap A — `cartridge_manifest()` exporter (most user-facing)

A single method that returns a stable JSON-serialisable dict:

```python
{
    "schema_version": 1,
    "identity": {species, locus, reference_set, name, source},
    "rules": {alphabet, v_anchor, j_anchor},
    "catalogue_counts": {v: int, d: int, j: int, c: int},
    "functional_status_histogram": {
        "v": {functional: int, orf: int, pseudogene: int, unknown: int, unannotated: int},
        "d": {...},
        "j": {...},
    },
    "models": {
        "reference_models_present": bool,
        "legacy_np_lengths_keys": [...],
        "legacy_trim_dicts_keys": [...],
    },
    "curation": {
        "applied": str | None,  # e.g. "functional_status:functional|keep_unannotated=true"
    },
    "hashes": {
        "compute_checksum": str,
        "content_hash": str,
    },
}
```

This is what the user spec sketches at the end ("A
`cartridge_manifest()` method that returns a stable
JSON-serializable summary…"). The audit recommends but doesn't
implement.

### Gap B — `reference_models` digest in identity (the v1 boundary)

A future slice that wants `content_hash` to catch
`reference_models` swaps should produce a stable digest of the
spec (e.g. `sha256` over `(np_lengths_keys + trims_keys +
per-distribution-canonical-values)`) and fold it into either:
- `identity.source` (cheap; matches the curation tag pattern), or
- A new identity field `models_digest` (cleaner).

Out of scope for the audit; pin the v1 boundary so the slice
lands deliberately.

### Gap C — Cartridge-view coverage of orphan simulation fields

The six "still affects simulation" `DataConfig` fields without a
view (Q1) should each either:
- Be folded into `CartridgeModelsView` as a typed sub-plane
  (e.g. `legacy_np_transitions`), OR
- Be formally classified as deprecated and pinned for removal.

The decision matters because a contributor today can change
`gene_use_dict` and produce a different `compute_checksum` (so
pickle integrity is preserved) but no view's docstring or test
points at the change. The audit recommends the fold-in; absence
pin documents the current state.

### Deferred (NOT in scope for "complete")

- Cross-platform reproducibility beyond what `content_hash`
  already provides (locale, endianness, etc.).
- Per-allele provenance (`anchor_meta`, `aliases`) reaching the
  Rust side. These are diagnostic surfaces; if a downstream tool
  needs them, surface them through Python-only accessors rather
  than the bridge.

---

## 8. Edge cases the contract covers

1. **Legacy pickle without `reference_rules` / `reference_models`.**
   The fields fall through to defaults via `__getattr__`; the
   soft-transition checksum policy pops them only when `None` so
   legacy pickles continue to verify against existing checksums.
   Pinned by `pin_scaffold_legacy_pickle_compute_checksum_stable`.

2. **`reference_models` spec that matches legacy defaults.** Two
   cartridges with the same effective NP/trim distributions
   produce different `compute_checksum` if one uses
   `reference_models` and the other uses the legacy dicts.
   `content_hash` is unaffected by either (neither crosses to
   Rust). Pin this as a known v1 quirk.

3. **Identity-only swap.** Renaming `cfg.metadata.species` from
   "Human" to "Homo sapiens" changes both hashes — and changes
   downstream cartridge identity attribution in traces. Pin the
   expected behaviour.

4. **Functional-status normalisation.** "F" / "f" / "Functional"
   / "functional" all collapse to `"functional"` via
   `_normalise_functional_status`. Unknown strings drop to
   `None`. This means a cartridge with an exotic functional
   status string normalises silently and the bridge produces an
   `Allele` with `functional_status = None`. Pinned so a future
   contributor adding a new IMGT status string updates both the
   Python set and the Rust enum in lockstep.

5. **Multiple curation steps in sequence.** Today's
   `refdata.curated(p1).curated(p2)` produces an
   `identity.source` with two `|curated:…` suffixes, and a
   distinct content hash. Pinned.

6. **C-segment dropped at bridge.** The engine has no C-segment
   passes, so `cfg.c_alleles` doesn't cross. `content_hash`'s
   `c_pool` section is always `c_pool:0` for cartridges built via
   `dataconfig_to_refdata`. Pinned as the current contract.

---

## 9. Backwards compatibility

This audit does **not** propose any user-visible change. The
gaps it identifies (manifest exporter, models digest, view
coverage of orphan fields) are slices a future contributor can
land independently. The contract pins make today's behaviour the
documented architecture; flipping any of them requires a slice +
a doc update in lockstep.

The one place a careful contributor must watch: the
`_KNOWN_STATUSES` constant ↔ Rust `FunctionalStatus` enum
correspondence. They're parallel sets in two languages; pinning
their content explicitly means a refactor that adds a new IMGT
status (e.g. `"P_partial"`) updates both sides at the same time.

---

## 10. Validator integration

(No validator changes proposed. The cartridge integrity
validator is `DataConfig.verify_integrity()`; the engine-side
validator is `refdata.validate(mode="strict"|"curatable")` —
both unchanged by this audit.)

---

## 11. Implementation order (recommended)

Following the audit's "narrow next slice" pattern:

1. **Slice 1 — `cartridge_manifest()` exporter** (Gap A).
   Python-only. Builds the JSON-serialisable summary from
   existing surfaces; no new Rust accessors needed if the
   functional-status histogram is computed by walking the Python
   allele list. Smallest user-facing improvement; closes the
   "five separate calls to inspect a cartridge" friction.
2. **Slice 2 — Functional-status histogram accessor**
   (refinement of Slice 1). Add `cfg.functional_status_counts()`
   or similar so the manifest's histogram has a load-bearing
   computation path that's independently testable.
3. **Slice 3 — `reference_models` digest in identity** (Gap B).
   Heavier — touches both Python (digest production) and Rust
   (identity.source format or new field). Probably needs its own
   audit pass.
4. **Slice 4 — View coverage of orphan simulation fields**
   (Gap C). Decision-heavy: each orphan needs a "fold-in vs
   deprecate" call.

This audit only proposes the architecture and pins absences;
each slice flips the relevant `pin_absence_*` to
`pin_present_*` in lockstep.

---

## 12. Test surface — what the audit pins

Mirrored in
[`tests/test_reference_cartridge_completeness_contract.py`](../tests/test_reference_cartridge_completeness_contract.py).

### `pin_scaffold_*` — today's contract

1. The four cartridge views expose their documented planes.
2. `compute_checksum` is stable for the bundled cartridge across
   no-op operations (pickle round-trip).
3. `reference_rules` (truly different from defaults) changes
   both `compute_checksum` AND `content_hash`.
4. `reference_models` (non-None) changes `compute_checksum` but
   does NOT change `content_hash` — pinned as the v1 boundary.
5. `functional_status` on alleles changes `content_hash`.
6. Curation policy tags `identity.source` with `|curated:…` and
   changes `content_hash`.
7. The 7 dropped Allele fields are documented (`aliases`,
   `anchor_meta`, `gapped_seq`, `family`, `locus`, `species`,
   `source`) plus `length` / `ungapped_len` redundancy.
8. The bundled `HUMAN_IGH_OGRDB` exposes the expected identity
   fields after `dataconfig_to_refdata`.
9. Functional-status normalisation: "F" → "functional", "ORF"
   → "orf", unknown → None.
10. The Python `_KNOWN_STATUSES` set matches the Rust
    `FunctionalStatus` enum's lowercase names.
11. Six DataConfig fields are not covered by any view today
    (the orphan list).

### `pin_absence_*` — gaps Slice 1+ closes

12. No `DataConfig.cartridge_manifest()` method.
13. No `DataConfig.functional_status_counts()` (or any
    `*histogram*` accessor on `DataConfig` / cartridge views).
14. No `models_digest` field in `RefDataConfig.identity()`.
15. No new `CartridgeSamplingDefaultsView` / equivalent covering
    the orphan simulation fields.

### Doc anchor

16. The audit doc continues to exist and references the contract
    file; the 14-section structure is intact.

---

## 13. Out of scope

Documented here so a future contributor doesn't accidentally
expand the audit.

- **Refactoring `DataConfig` storage** — the audit only documents
  the current dataclass layout; restructuring is a separate
  slice with its own migration story.
- **Eliminating the bridge's "dropped" fields** — the seven
  dropped Allele fields are intentionally Python-only; moving
  any of them to the Rust side requires its own design.
- **Cross-process cartridge identity beyond `content_hash`** —
  e.g. signing manifests, distributing cartridge bundles. Out
  of scope.
- **Per-allele provenance migration** — `anchor_meta` is
  Python-only by design; a slice that wants it engine-side is
  separate.
- **Replacing pickle with a stable on-disk format** — would
  change `compute_checksum` definition. Separate audit.

---

## 14. Summary table

| Concern | Today's contract |
|---|---|
| Cartridge views | Four: identity, catalogue, rules, models. |
| `DataConfig` fields not in any view | `gene_use_dict`, `NP_transitions`, `NP_first_bases`, `correction_maps`, `asc_tables`, `p_nucleotide_length_probs`, `dj_pairing_map`. |
| `Allele` fields dropped at bridge | `aliases`, `anchor_meta`, `gapped_seq`, `family`, `locus`, `species`, `source` (`locus`/`species`/`source` re-exposed at cartridge level via `identity`). |
| `compute_checksum` coverage | Everything in `DataConfig.__dict__` minus zeroed `schema_sha256`, removed `build_report`, and `None`-popped `reference_rules` / `reference_models`. |
| `content_hash` coverage | `chain` + identity (5 fields) + rules (alphabet + 2 anchor rules) + pools (name/gene/segment/seq/anchor/functional_status). |
| `reference_rules` in identity | Yes — both hashes. |
| `reference_models` in identity | **Python only** — `compute_checksum` covers, `content_hash` does NOT. **v1 boundary.** |
| Functional status in identity | Yes — `content_hash`. |
| Curation policy in identity | Yes — `identity.source` tagged. |
| Bundled cartridge metadata | Identity, rules, catalogue, content hash — all present. Manifest exporter and functional-status histogram — absent. |
| Required for "complete" | (1) `cartridge_manifest()` exporter, (2) `reference_models` digest in identity, (3) view coverage of orphan simulation fields. |

The reference cartridge today is **operationally complete** —
every cartridge that ships exercises all four planes correctly,
the Rust bridge transfers what the engine needs, and identity is
stable for trace attribution. What's **not** complete is the
*inspectability* surface: users have to assemble cartridge
provenance from five separate calls, the orphan simulation
fields lack a documented home, and `reference_models` swaps
don't show up in the Rust-side identity.

The recommended next slice (`cartridge_manifest()` exporter)
closes the most user-visible gap with the smallest Python-only
diff. The deeper gaps (models digest, orphan-field absorption)
each warrant their own slice + audit cycle.
