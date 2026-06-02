# V-Region Substructure / CDR-FR Targeting — Audit + Slice 1 Shipped

**Status: Slice 1 (V-Subregion Cartridge Annotation Surface)
shipped; Slices 2 and 3 deferred.** The annotation plumbing
described below — Python `Allele.subregions`, Rust
`VSubregion` / `VSubregionLabel`, bridge-time derivation,
manifest `v_subregion_support` block, content-hash
participation, bridge-time validation — is now landed and
pinned by [`tests/test_v_subregion_annotations.py`](../tests/test_v_subregion_annotations.py).
Per-region SHM **sampling** (`v_subregion_rates` kwarg) and
per-region **mutation counters** (`n_cdr1_mutations`, …) remain
explicitly deferred — they are the Slice 2 / Slice 3 work and
are still pinned as `pin_absence_*` in
[`tests/test_v_region_substructure_contract.py`](../tests/test_v_region_substructure_contract.py).

This audit was originally pre-implementation. The body below
preserves the original framing for traceability; sections marked
**[Shipped — Slice 1]** describe how the recommendation actually
landed.

The audit pins today's V-subregion metadata surface — what
bundled cartridges already carry, what's derivable, what's
dropped at the bridge, what's absent — and specifies the next
slice. The Slice 1 deliverable is the shared vocabulary plus
contract pins so any future slice (CDR/FR rates, per-region
counters) lands against an audited baseline.

This audit is the natural follow-up to the targeted-SHM rates +
provenance counter slices: per-segment rates work at the V vs D
vs J vs NP granularity, but biological SHM is usually discussed
across V-region framework (FR) and complementarity-determining
regions (CDR). Before adding rates / counters at that finer
resolution, we need to know whether the cartridge and allele
model can support V subregion annotations cleanly.

Companion to
[`tests/test_v_region_substructure_contract.py`](../tests/test_v_region_substructure_contract.py)
which freezes today's surfaces (`pin_scaffold_*`) and the gaps
that an implementation slice would close (`pin_absence_*`).

**Pre-flight finding:** bundled cartridges already carry the
metadata needed (`gapped_seq` on every V allele, IMGT-gap
convention), AND a derivation helper exists
([`src/GenAIRR/utilities/imgt_regions.py`](../src/GenAIRR/utilities/imgt_regions.py))
that produces clean ungapped subregion boundaries on all 198
bundled V alleles. The next slice is a **cartridge annotation
surface**, not new data collection.

---

## Existing scaffolding the contract pins

| Surface | Where | What it pins |
|---|---|---|
| `Allele.gapped_seq` | [`src/GenAIRR/alleles/allele.py`](../src/GenAIRR/alleles/allele.py) | IMGT-gapped nucleotide sequence with `.` characters at gap positions. Carried on every bundled V allele. **Dropped at the Python→Rust bridge** ([`_refdata_resolver.py`](../src/GenAIRR/_refdata_resolver.py)). |
| `Allele.anchor_meta` | [`src/GenAIRR/alleles/allele.py:44-66`](../src/GenAIRR/alleles/allele.py#L44-L66) | T2-8 anchor provenance (codon, residue, confidence, method, rejection reason). Field exists on the Python `Allele` class, but is `None` on every bundled allele today — the resolver builds it on construction but bundled pickles predate the field. **Dropped at the bridge.** |
| `Allele.anchor: Optional[int]` | [`src/GenAIRR/alleles/allele.py:92`](../src/GenAIRR/alleles/allele.py#L92) | Ungapped position of the conserved Cys (V) / Trp/Phe (J) codon. **Crosses the bridge** as `RefDataConfig.Allele.anchor: Option<u16>`. |
| `IMGT_GAPPED_BOUNDARIES` + `compute_v_region_boundaries` | [`src/GenAIRR/utilities/imgt_regions.py:23-69`](../src/GenAIRR/utilities/imgt_regions.py#L23-L69) | Maps fixed IMGT gapped nucleotide positions (FWR1=0–78, CDR1=78–114, FWR2=114–165, CDR2=165–195, FWR3=195–312) to ungapped per-allele positions by counting non-gap characters. Pure Python. |
| `classify_position` | [`src/GenAIRR/utilities/imgt_regions.py:72-119`](../src/GenAIRR/utilities/imgt_regions.py#L72-L119) | Assembled-sequence-position → region label (`FWR1` / `CDR1` / `FWR2` / `CDR2` / `FWR3` / `CDR3` / `FWR4` / `NP`). Consumed by the MCP analysis helpers; **not consumed by the simulation pipeline.** |
| `Rust RefDataConfig::Allele` | [`engine_rs/src/refdata.rs`](../engine_rs/src/refdata.rs) | Fields: `name`, `gene`, `seq`, `segment`, `anchor`, `functional_status`, **`subregions: Vec<VSubregion>` (Slice 1)**. |
| `cartridge_manifest()["models"]["shm"]` | [`src/GenAIRR/dataconfig/data_config.py`](../src/GenAIRR/dataconfig/data_config.py) | Carries `segment_rate_support` (the V/D/J/NP bucket capability) and **`v_subregion_support` (Slice 1)**. |
| `Experiment.mutate(...)` signature | [`src/GenAIRR/experiment.py`](../src/GenAIRR/experiment.py) | Accepts `model` / `count` / `rate` / `s5f_model` / `segment_rates`. **No `v_subregion_rates` (Slice 2).** |
| AIRR record fields | [`engine_rs/src/airr_record/record.rs`](../engine_rs/src/airr_record/record.rs) | Carry `n_v_mutations` / `n_d_mutations` / `n_j_mutations` / `n_np_mutations` (per-segment SHM counters). **No CDR/FR counters (Slice 3).** |

---

## 1. Q1 — What V-subregion metadata exists today?

### Per-allele Python surfaces

- **`Allele.gapped_seq`** — present on every bundled V allele
  (198/198 in `HUMAN_IGH_OGRDB`, 33/33 D alleles, 7/7 J
  alleles). IMGT-gap convention (`.` characters at gap
  positions). Sample for `IGHVF1-G1*01`: 322 gapped chars / 301
  ungapped / 21 gaps. The gapped sequence aligns canonical IMGT
  codon positions so the same gapped position maps to the same
  biological structural location across alleles.
- **`Allele.anchor: Optional[int]`** — ungapped position of the
  conserved CDR3-defining codon (Cys for V at IMGT codon 104).
  Present on every bundled allele.
- **`Allele.anchor_meta: Optional[AnchorResult]`** — anchor
  provenance metadata (codon bytes, residue, confidence
  level, derivation method, rejection reason on failures).
  Field exists on the Allele class but is **`None` on every
  bundled allele** in `HUMAN_IGH_OGRDB`. The
  [`VAllele._find_anchor`](../src/GenAIRR/alleles/allele.py#L168-L203)
  resolver populates it at runtime construction, but bundled
  pickles were built before the field was added so it doesn't
  round-trip.

### Derivation helper (already shipping)

[`src/GenAIRR/utilities/imgt_regions.py`](../src/GenAIRR/utilities/imgt_regions.py)
ships a complete derivation surface:

- `IMGT_GAPPED_BOUNDARIES` — fixed `(start, end)` gapped
  nucleotide positions for FWR1 / CDR1 / FWR2 / CDR2 / FWR3.
  IMGT-canonical, locus-independent.
- `compute_v_region_boundaries(v_allele)` — counts non-gap
  characters to produce `{region: (ungapped_start,
  ungapped_end)}` per allele.
- `classify_position(pos, v_seq_start, v_seq_end, v_boundaries,
  junction_start, junction_end, j_seq_end)` — maps an
  assembled-sequence position to one of `FWR1` / `CDR1` /
  `FWR2` / `CDR2` / `FWR3` / `CDR3` / `FWR4` / `NP`.

The helper runs cleanly across all 198 V alleles in
`HUMAN_IGH_OGRDB`. Sample for `IGHVF1-G1*01` (ungapped length
301):

| Region | Ungapped `(start, end)` | Length |
|---|---|---|
| FWR1 | (0, 75) | 75 |
| CDR1 | (75, 105) | 30 |
| FWR2 | (105, 156) | 51 |
| CDR2 | (156, 177) | 21 |
| FWR3 | (177, 291) | 114 |

Biologically sensible — matches IMGT codon-counted lengths.

### Builder scripts dropping subregion data

The cartridge build path produces `Allele` instances populated
from the loader's gapped sequence; `gapped_seq` survives onto
the Python `DataConfig.v_alleles` dict. **Nothing is dropped at
the Python cartridge layer.** What IS dropped — at the
Python→Rust bridge — is documented in §2.

### Currently surfaced where?

The derivation helper is consumed only by
[`mcp_helpers.py`](../src/GenAIRR/utilities/mcp_helpers.py)
(MCP analysis tooling). It is **not consumed by the simulation
pipeline**, not by the cartridge manifest, not by any AIRR
projection, not by the SHM pass.

Pinned by:
- `pin_scaffold_bundled_v_alleles_carry_gapped_seq` — 198/198.
- `pin_scaffold_anchor_meta_field_exists_on_python_allele` —
  class attribute present.
- `pin_absence_anchor_meta_populated_on_bundled_alleles` —
  none of the 198 V alleles in `HUMAN_IGH_OGRDB` have
  `anchor_meta` populated.
- `pin_scaffold_imgt_regions_helper_exists` — the derivation
  helper ships.
- `pin_scaffold_compute_v_region_boundaries_works_on_bundled_alleles` —
  all 198 alleles yield five clean ungapped intervals.

---

## 2. Q2 — Does Rust `RefDataConfig` receive any V subregion info?

**No.** Pinned source-level on the Rust struct.

The Python→Rust bridge in
[`_refdata_resolver.py::_push_alleles`](../src/GenAIRR/_refdata_resolver.py#L153-L178)
forwards exactly:

```python
allele.name           # → engine Allele.name
allele.gene           # → engine Allele.gene
allele.ungapped_seq   # → engine Allele.seq (bytes)
allele.anchor         # → engine Allele.anchor
allele.functional_status  # → engine Allele.functional_status
```

Dropped at the bridge:
- `gapped_seq` (Python-only; ungapped goes across).
- `anchor_meta` (Python-only, always `None` on bundled data
  anyway).
- `aliases`, `family`, `locus`, `species`, `source` (cartridge-
  level only).

The Rust `RefDataConfig::Allele` struct has no V-subregion
fields:

```rust
pub struct Allele {
    pub name: String,
    pub gene: String,
    pub seq: Vec<u8>,
    pub segment: Segment,
    pub anchor: Option<u16>,
    pub functional_status: Option<FunctionalStatus>,
}
```

If a future slice wants Rust-side subregion awareness, it must
either:
- Compute subregions in Python at cartridge build time and
  forward as new fields on the Rust `Allele` (a `Vec<(label,
  start, end)>` or four `Option<u16>` start positions), OR
- Send `gapped_seq` to Rust and re-derive at compile time (more
  bandwidth, less clear ownership).

The audit's §4 recommendation is the first option.

Pinned by:
- `pin_scaffold_v_subregion_absent_from_rust_allele_struct` —
  source-level grep on the field list.
- `pin_absence_v_subregion_not_forwarded_at_bridge` —
  source-level pin that `_push_alleles` doesn't forward any
  subregion info.

---

## 3. Q3 — Can we derive subregions from gapped IMGT sequences?

**Yes, and we already do.** The
[`imgt_regions.py`](../src/GenAIRR/utilities/imgt_regions.py)
helper:

- Uses IMGT's invariant gapped-nucleotide boundaries (FWR1
  starts at gapped position 0, CDR1 at 78, FWR2 at 114, CDR2 at
  165, FWR3 at 195, FWR3 ends at 312 — the position of the
  conserved Cys codon).
- Counts non-gap characters in `gapped_seq[:gapped_pos]` to
  convert each boundary to the per-allele ungapped position.
- Runs in O(allele_len) per allele.

The derivation produces clean intervals on all 198 bundled
`HUMAN_IGH_OGRDB` V alleles (audit probe).

### Where should the derivation live?

**Recommended: at cartridge build time, not at simulation
runtime.** Three reasons:

1. The boundaries are deterministic and per-allele. Computing
   them once at build time costs nothing; computing them
   per-record (e.g. inside a hypothetical CDR-aware SHM pass)
   would re-derive constants on every record.
2. The Rust runtime doesn't have access to `gapped_seq` (it's
   Python-only). Pushing the derivation to Rust requires either
   bridging `gapped_seq` or pre-computing in Python.
3. Other consumers (cartridge inspection, AIRR projection
   per-position labelling, MCP analysis) want the same
   pre-computed intervals.

The natural shape:

```python
# At cartridge build time (or on first cartridge access):
v_allele.subregions = {
    "FWR1": (start, end),
    "CDR1": (start, end),
    "FWR2": (start, end),
    "CDR2": (start, end),
    "FWR3": (start, end),
}
```

Then bridge `subregions` across as a new field on the Rust
`Allele` — see §4 below.

### When derivation fails

The helper's `compute_v_region_boundaries` clamps to the gapped
sequence length, so a short / truncated gapped_seq still
produces an answer (just with empty trailing intervals). No
allele in `HUMAN_IGH_OGRDB` triggers the clamp.

Pinned by:
- `pin_scaffold_imgt_gapped_boundaries_match_imgt_aa_positions` —
  the constants match the documented IMGT AA-position
  ranges (1-26 / 27-38 / 39-55 / 56-65 / 66-104).

---

## 4. Q4 — How should the cartridge represent subregions?

### Recommended shape: per-allele dict on `Allele`

```python
@dataclass
class VAllele(Allele):
    # … existing fields …
    subregions: Optional[Dict[str, Tuple[int, int]]] = None
```

Coordinates **ungapped, allele-relative, half-open**. Same
convention as the existing `anchor` field. The dict carries the
five canonical IMGT regions:
`{"FWR1", "CDR1", "FWR2", "CDR2", "FWR3"}`. CDR3 / FWR4 are
**not** per-V-allele — they're defined by the recombined
sequence (CDR3 spans the V anchor → J anchor; FWR4 is the
post-junction J region). The cartridge V allele's `subregions`
dict is V-only.

### Why not separate fields per region

Option A (per-region fields, `fwr1_start: int`, `fwr1_end: int`,
etc.) bloats the Allele class with ten new fields. Option B
(`subregions: Dict[str, Tuple[int, int]]`) is one field and
preserves the IMGT region-name strings users will read in docs.

Recommendation: Option B. Single field, JSON-clean, matches the
helper's existing return shape.

### Rust bridge shape

A new `Vec<(SubregionLabel, u16, u16)>` field on the Rust
`Allele` struct, where `SubregionLabel` is a tightly-packed
enum (`Fwr1 = 0`, `Cdr1 = 1`, … `Fwr3 = 4`). Positions stay
`u16` matching `anchor`. Empty `Vec` when no subregions were
attached (D / J / VJ-chain V alleles where the derivation
doesn't apply).

Alternative: four `Option<u16>` start positions (FWR1 always
starts at 0). Compacter but less self-describing. The audit
recommends the Vec form.

### Content hash

Subregion boundaries should fold into
`refdata.content_hash()` so two cartridges with the same
allele bytes but different region annotations (e.g. an updated
IMGT spec) produce different hashes. Same shape as the
`anchor` field's contribution today.

Pinned by:
- `pin_absence_no_subregions_field_on_python_v_allele` — the
  proposed field doesn't exist yet.
- `pin_absence_no_subregion_label_enum_in_rust` — the proposed
  Rust enum doesn't exist.

---

## 5. Q5 — How should SHM use them?

**Out of scope for this audit.** The audit recommends a two-
slice path:

### Slice 1 — Cartridge annotation surface (recommended next)

Add `subregions` field to Python `VAllele` populated at
cartridge build time via `compute_v_region_boundaries`. Bridge
across to Rust as a `Vec<(SubregionLabel, u16, u16)>` on
`Allele`. Extend `cartridge_manifest()` with coverage statistics
(see §6).

**No SHM API change in this slice.** Just data plumbing.

### Slice 2 — `v_subregion_rates` kwarg on `mutate`

Future shape:

```python
.mutate(
    model="s5f",
    rate=0.03,
    segment_rates={"V": 1.0, "D": 0.5, "J": 0.5, "NP": 0.0},
    v_subregion_rates={"FWR1": 1.0, "CDR1": 4.0, "FWR2": 1.0,
                       "CDR2": 4.0, "FWR3": 1.0},
)
```

Composes with `segment_rates`: the V site weight becomes
`segment_rates["V"] * v_subregion_rates[subregion_at(site)]`.
Sites in V-region-derived `FR` regions get the FR rate; CDR
sites get the CDR rate; sites outside V are unaffected.

The biological convention: CDR positions mutate at ~3–5× the
FR rate (Lefranc et al., Yaari et al.). The kwarg lets users
quantify the bias explicitly without changing the kernel.

### Slice 3 — CDR/FR mutation counters

Add four AIRR fields:

```text
n_cdr1_mutations
n_cdr2_mutations
n_fwr1_mutations
n_fwr2_mutations
n_fwr3_mutations
```

Aggregation: walk `outcome.events()` filtered to the SHM passes,
classify each `BaseChanged.handle` against the V-region
boundaries from the IR (via the cartridge-pushed subregion
data), bucket. Same pattern as the existing per-segment counter
slice.

### Why three slices, not one

The audit's reasoning: each slice has its own validator
surface, replay semantics, and manifest extension. Bundling
them risks shipping a half-typed CDR concept (rates without
counters, or vice versa). The three-slice progression matches
the SHM model evolution (targeting → counters → release
consolidation) the targeted-SHM family already used.

---

## 6. Q6 — Validator / manifest

### Validator extensions

Any V-subregion slice that adds AIRR fields gets the same
re-derivation + mismatch issue pattern:

- `SubregionsBoundariesMismatch { segment, reported, expected }`
  — fires when an allele's subregions don't match the
  cartridge's stored values.
- `NCdr1MutationsMismatch` / `NCdr2MutationsMismatch` /
  `NFwr1MutationsMismatch` etc. — per-region counter mismatches
  (Slice 3).
- `SubregionSumMismatch` — sum-invariant cross-check
  (`n_cdr1 + n_cdr2 + n_fwr1 + n_fwr2 + n_fwr3 == n_v_mutations`)
  for Slice 3.

The slice should also surface a validator-side rejection for
**malformed subregion annotations** — overlapping intervals,
out-of-bounds end positions, decreasing region order, etc.
Same shape as the existing `RefDataConfig::validate(...)`
allele-level checks.

### Manifest extension

Slice 1's manifest extension would surface coverage statistics:

```python
manifest["models"]["shm"]["v_subregion_support"] = {
    "available": bool,
    "regions": ["FWR1", "CDR1", "FWR2", "CDR2", "FWR3"],
    "coverage": {
        "v_alleles_total": int,
        "v_alleles_with_subregions": int,
        "per_region_coverage": {
            "FWR1": int,  # allele count with this region's interval
            "CDR1": int,
            ...
        },
    },
    "derivation": "imgt_gapped_boundaries",
    "in_content_hash": True | False,
}
```

`in_content_hash` defaults to `True` (subregion boundaries are
allele identity per §4 recommendation). The slice can decide
whether to expose it as a v1 boundary like
`segment_rate_support` does.

Pinned by:
- `pin_absence_no_v_subregion_support_in_manifest`.
- `pin_absence_no_subregion_mismatch_validator_kinds`.

---

## 7. Edge cases the implementation slice must handle

1. **VJ chains (no D).** Subregion derivation works on V
   independently of D presence — the helper consumes only
   `v_allele.gapped_seq`. No change needed.

2. **TCR loci.** Should subregions be applied to TCR V alleles?
   Yes, IMGT defines CDR/FR for TCR V genes too. The derivation
   helper today uses the same IMGT_GAPPED_BOUNDARIES for any V
   gene — locus-independent. Pin the locus-agnostic behaviour.

3. **C-region alleles.** No CDR/FR on C — the derivation is
   V-only. The `subregions` field would be `None` for D / J /
   C alleles. Pin this clearly.

4. **Alleles missing `gapped_seq`.** Currently 0 / 198 in the
   bundled cartridge but possible for user-built cartridges
   from FASTA without IMGT alignment. The derivation should
   return `None` rather than raise; the manifest's coverage
   statistics surface the gap.

5. **Truncated / non-canonical V alleles.** Some V alleles in
   the wild are shorter than the canonical IMGT length (e.g.
   partial sequences). The helper today clamps to gapped length
   — the trailing regions come back as empty intervals. Pin
   the clamp behaviour.

6. **Receptor revision.** Replaces V mid-pipeline; the post-
   revision V's subregions apply, not the pre-revision V's.
   `Simulation.assignments.v.allele_id` is read at SHM time, so
   the SHM pass would consult the current V's subregions. No
   special handling needed.

7. **D-inverted records.** D-inversion doesn't affect V
   subregions. No interaction.

8. **AIRR record at TCR.** TCR records get the same subregion
   fields (when Slice 3 lands), or stay at zero if the
   cartridge has no TCR V subregion data attached.

---

## 8. Performance

### Per-allele derivation cost

`compute_v_region_boundaries` is O(allele_len × n_regions) ≈
O(300 × 5) = O(1500) character scans. Computed once at
cartridge build time (or first cartridge access); cached
thereafter. **Zero per-record cost.**

### Per-record classification cost (Slice 3)

If a future CDR/FR counter walks `BaseChanged` events and
classifies each position, the cost is O(n_mutations × log
n_regions) per record (binary search over the five intervals).
Negligible vs the SHM pass itself.

### Bridge bandwidth

Pushing five `(start, end)` pairs per V allele to Rust is ~40
bytes per allele × 198 alleles = ~8 KB per cartridge. One-time
cost at `dataconfig_to_refdata` time.

---

## 9. Trace / replay

**No new trace addresses.** Subregion boundaries are allele
identity (derived once at cartridge build / load time); they
don't enter per-record sampling. Replay reproduces because the
same cartridge produces the same boundaries.

Future Slice 2 (`v_subregion_rates` kwarg) is a compile-time
pass parameter — same shape as `segment_rates`. The plan
signature catches rate-vector drift at trace-file replay; no
new addresses.

Future Slice 3 (CDR/FR counters) is event-derived — same
shape as the per-segment SHM counters. No new addresses.

---

## 10. Manifest integration

Covered in §6. The Slice 1 extension adds a
`v_subregion_support` block under `models.shm` with coverage
counts + derivation provenance.

---

## 11. Implementation order

After this audit, three slices in order of dependency. Slice 1
has shipped; Slices 2 and 3 are still deferred.

### Slice 1 — Cartridge annotation surface **[Shipped]**

Scope as landed:

1. **Python `Allele.subregions`** —
   [`src/GenAIRR/alleles/allele.py`](../src/GenAIRR/alleles/allele.py)
   gets a class-level `subregions = None` default. V-only by
   intent; D / J alleles ignore it. Users can override with a
   `{label: (start, end)}` dict; the bridge validates it.
2. **Bridge-time derivation** —
   [`src/GenAIRR/_refdata_resolver.py`](../src/GenAIRR/_refdata_resolver.py)
   gains `_resolve_v_subregions(allele)`. For each V allele the
   bridge honors a user-set `allele.subregions` dict first, then
   falls back to deriving from `gapped_seq` via the existing
   `compute_v_region_boundaries` helper. Failures are swallowed
   into an empty list — a single allele's missing metadata is
   coverage, not a fatal error. D / J skip the dispatch
   entirely. Bundled cartridges are NOT re-pickled.
3. **Rust types** —
   [`engine_rs/src/refdata.rs`](../engine_rs/src/refdata.rs)
   carries `VSubregionLabel { Fwr1, Cdr1, Fwr2, Cdr2, Fwr3 }`,
   a `VSubregion { label, start, end }` struct, and a
   `pub subregions: Vec<VSubregion>` field on `Allele`. Other
   segments always carry an empty `Vec`.
4. **PyO3 wire format** —
   [`engine_rs/src/python/refdata.rs`](../engine_rs/src/python/refdata.rs)
   extends `add_v_allele` with a `subregions=None` kwarg taking
   a `List[(str, u16, u16)]`. `parse_subregions` validates:
   canonical label set (case-sensitive), no duplicate labels,
   `start < end`, `end <= seq_len`, no overlapping intervals.
   Each violation raises at bridge time with a clear message.
5. **Content hash** —
   [`engine_rs/src/trace_file.rs`](../engine_rs/src/trace_file.rs)
   folds `label:start-end` into the per-allele hash bytes.
   Cartridges with different subregion intervals hash
   differently. Empty list (D/J + legacy V without metadata)
   hashes as `subregions=` with no payload — backwards
   compatible for pools where no V has annotations.
6. **Manifest extension** —
   [`src/GenAIRR/dataconfig/data_config.py`](../src/GenAIRR/dataconfig/data_config.py)
   gains `models.shm.v_subregion_support` with
   `{available, labels, annotated_v_count, total_v_count,
   derivation, in_content_hash}`. The bridge walks the V pool
   to compute coverage; `available=True` iff any V allele has
   non-empty subregions. `derivation="bridge_imgt_gapped_seq"`.
   `in_content_hash=True`.
7. **PyAllele accessor** — `PyAllele.subregions` returns
   `List[(label, start, end)]`. Empty for D / J and for V
   alleles without metadata.
8. **Bundled coverage** — IGH 198/198, IGK 168/168, IGL 181/181
   on the OGRDB cartridges. Full derivation from IMGT-gapped
   sequences; legacy alleles without `gapped_seq` load cleanly
   with empty subregions.

Validation discipline: the validator surface lives at **bridge
time** (`parse_subregions` raises during `dataconfig_to_refdata`),
**not** at AIRR record validation time. Slice 3 would add
record-time `*SubregionMismatch` validator issue kinds; those
are still absent and pinned.

What Slice 1 did NOT do (still deferred):

- `Experiment.mutate(v_subregion_rates=...)` kwarg → Slice 2.
- `n_cdr1_mutations` / `n_cdr2_mutations` / `n_fwr1_mutations`
  / `n_fwr2_mutations` / `n_fwr3_mutations` AIRR fields →
  Slice 3.
- `*SubregionMismatch` AIRR-record validator issue kinds →
  Slice 3.

### Slice 2 — `v_subregion_rates` SHM kwarg

Builds on Slice 1. Adds compile-time pass parameter, weights
SHM site selection by `segment_rate × subregion_rate` for V
positions. Same shape as `segment_rates`.

### Slice 3 — CDR/FR mutation counters

Builds on Slice 1. Five AIRR fields, event-derived
aggregation, validator extension, manifest extension. Same
shape as the per-segment SHM counter slice.

### Why this order

Slice 1 unblocked Slices 2 and 3 (both need cartridge
subregion data). Slice 2 + Slice 3 can ship in either order
now — they don't depend on each other.

The next biology slice is the V-subregion SHM rate audit →
Slice 2 implementation.

---

## 12. Test surface — what this audit pins

Mirrored in
[`tests/test_v_region_substructure_contract.py`](../tests/test_v_region_substructure_contract.py).

### `pin_scaffold_*` — pre-existing contract that Slice 1 builds on

1. Bundled V alleles carry `gapped_seq` (198/198 in
   `HUMAN_IGH_OGRDB`).
2. `Allele.anchor_meta` field exists on the Python class.
3. `Allele.anchor` (ungapped position) crosses the Rust
   bridge.
4. `imgt_regions.IMGT_GAPPED_BOUNDARIES` exists with the five
   documented intervals.
5. `compute_v_region_boundaries` works on all 198 bundled V
   alleles.
6. IMGT gapped boundaries match the documented IMGT AA
   positions (1-26 / 27-38 / etc.).

### `pin_present_*` — Slice 1's landed surface

7. **Bridge forwards the five direct fields + delegates to
   `_resolve_v_subregions`.** `anchor_meta` still does NOT
   cross the bridge module (a separate cartridge-refresh
   surface).
8. **Rust `RefDataConfig::Allele` carries a `subregions`
   field** (`Vec<VSubregion>`).
9. **Python `Allele.subregions` attribute** exists on bundled
   instances.
10. **Rust declares `VSubregionLabel` enum + `VSubregion`
    struct** in [`engine_rs/src/refdata.rs`](../engine_rs/src/refdata.rs).
11. **`cartridge_manifest()["models"]["shm"]` carries the
    `v_subregion_support` block** with the six documented
    keys, five canonical labels, full coverage on bundled
    cartridges, and `in_content_hash=True`.
12. **`imgt_regions` is consumed by both the MCP helpers AND
    `_refdata_resolver`** — the bridge derivation surface.

### `pin_absence_*` — gaps Slice 2 / Slice 3 still close

13. No `v_subregion_rates` kwarg on `Experiment.mutate`
    (Slice 2).
14. No CDR/FR mutation counter fields on AIRR records
    (Slice 3): `n_cdr1_mutations`, `n_cdr2_mutations`,
    `n_fwr1_mutations`, `n_fwr2_mutations`,
    `n_fwr3_mutations`, `n_cdr_mutations`, `n_fwr_mutations`.
15. No `*SubregionMismatch` validator issue kinds on the
    Rust AIRR-record validator (Slice 3 — Slice 1's
    validation is bridge-time, not record-time).
16. `anchor_meta` is `None` on every bundled V allele in
    `HUMAN_IGH_OGRDB` (cartridge-refresh task, orthogonal).

### Doc anchor

17. The audit doc references the contract file; 14-section
    structure intact.

---

## 13. Out of scope

Documented here so a future slice author doesn't accidentally
expand the work.

- **Per-region (CDR1/CDR2 vs CDR3) SHM kernels.** S5F kernels
  treat the whole V region uniformly; per-region kernels are a
  separate biology slice well past Slice 3.
- **AID hotspot motif targeting.** S5F's 5-mer context kernel
  approximates AID hotspots empirically; explicit motif-aware
  targeting is out.
- **CDR3 sub-substructure** (D-derived junctional N-additions
  vs P-additions). The audit's NP bucket already covers NP1 /
  NP2; finer CDR3 internal structure is separate.
- **VJ chain CDR3 vs VDJ CDR3.** The cartridge's V subregions
  are V-only; CDR3 length variation between chain types is
  out of scope for the V annotation surface.
- **C-region subregions.** C alleles have CH1 / CH2 / CH3 / CH4
  / hinge but the engine doesn't model C-region passes today.
  Out of scope.
- **Re-pickling bundled cartridges to populate `anchor_meta`.**
  A separate "cartridge refresh" maintenance task; orthogonal
  to this audit.
- **TCR-specific subregion variants.** IMGT TCR V boundaries
  match the heavy-chain Ig V positions; if locus-specific
  divergence emerges (e.g. invariant TCR-α regions), that's a
  follow-up.

---

## 14. Summary table

| Concern | Post-Slice-1 state | Status |
|---|---|---|
| V-subregion metadata on bundled cartridges | `gapped_seq` present on 198/198 IGH V, 168/168 IGK V, 181/181 IGL V; `anchor_meta` still `None` | Bridge-time derivation covers all three; `anchor_meta` refresh remains a separate task |
| Subregion derivation helper | `imgt_regions.compute_v_region_boundaries` consumed by MCP **and** `_refdata_resolver._resolve_v_subregions` | **Shipped** |
| Rust `Allele` carries subregion data | `subregions: Vec<VSubregion>` field; `VSubregionLabel { Fwr1, Cdr1, Fwr2, Cdr2, Fwr3 }` enum | **Shipped** |
| Bridge forwards subregion data | Derived `(label, start, end)` triples cross via `add_v_allele(subregions=…)`; D / J always empty | **Shipped** |
| Cartridge manifest reports coverage | `models.shm.v_subregion_support` block (6 keys, 5 canonical labels, derivation tag, hash flag) | **Shipped** |
| Content hash includes subregion intervals | `refdata_content_hash` folds `label:start-end` into per-allele bytes | **Shipped** |
| Validator catches malformed subregions | Bridge-time `parse_subregions` rejects unknown label, duplicate label, `start>=end`, out-of-bounds, overlap | **Shipped (bridge-time)** |
| `Experiment.mutate(v_subregion_rates=...)` | No kwarg | Deferred — Slice 2 |
| AIRR CDR/FR counter fields | None | Deferred — Slice 3 |
| AIRR-record `*SubregionMismatch` validator kinds | None | Deferred — Slice 3 |
| Pre-flight bugs found | **None.** The data + helper exist; the gap was plumbing. | — |

The V-region substructure annotation surface is **shipped**:
bundled cartridges expose per-allele subregion intervals through
the bridge, manifest, and content hash. Inspectability and
cartridge identity are live; SHM **sampling** (rates) and
mutation **counters** are the next two biology slices.

The natural next step is the **V-subregion SHM rate audit**
preceding Slice 2 implementation. Slice 3 (CDR/FR counters) can
follow in either order.
