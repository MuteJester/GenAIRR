# Junction-call provenance audit

**Status:** audit + golden tests. Drift items in §6.

Counterpart to the
[productive-failure-mode](productive_failure_mode_audit.md) and
[allele-call](allele_call_audit.md) audits. This audit catalogues
how the AIRR junction-related fields (`junction`, `junction_aa`,
`junction_length`, `junction_start`, `junction_end`,
`vj_in_frame`, `stop_codon`, `productive`) are derived from the
post-pipeline IR, and how they behave under each evidence-changing
mechanism (recombination trims, NPs, D segments, indels, end-loss,
reverse-complement projection, SHM).

Audit-first pattern matching the prior slices — document current
behaviour, pin with golden tests, identify drift only after tests
force it out.

---

## 1. Field inventory + non-exposure

The AIRR record exposes the junction window directly:

| Field             | Type            | Origin                                                              |
|-------------------|-----------------|---------------------------------------------------------------------|
| `junction`        | `str`           | `rec.sequence[junction_start:junction_end]` (post-rev-comp aware).  |
| `junction_aa`     | `str`           | `translate_seq(junction)` — empty string `""` when out-of-frame.    |
| `junction_length` | `int`           | `len(junction)` (byte count, not coordinate delta).                 |
| `junction_start`  | `Optional[int]` | V_anchor pool position; coordinates flip under rev-comp.            |
| `junction_end`    | `Optional[int]` | J_anchor pool position + 3; coordinates flip under rev-comp.        |
| `vj_in_frame`     | `Optional[bool]`| `junction_length % 3 == 0`.                                         |
| `stop_codon`      | `Optional[bool]`| Stop codon (TAA/TAG/TGA) in any junction codon. *Only checked when in-frame.* |
| `productive`      | `Optional[bool]`| Composite: in-frame ∧ no junction stop ∧ V/J anchor amino acids preserved. |

**Not exposed today** (pinned by §7 tests):

| Expected field    | Status                                                                                   |
|-------------------|------------------------------------------------------------------------------------------|
| `cdr3`            | Not exposed. The junction window IS the CDR3 window per IMGT convention. No separate field. |
| `cdr3_aa`         | Same — not exposed.                                                                       |
| `v_anchor_codon_pos` / `j_anchor_codon_pos` | Not exposed. Anchor pool positions are internal to the junction calculation. |

Whether to add `cdr3` aliases is a product question (§6.2).

---

## 2. Junction window boundaries — the canonical rule

The junction window is computed by
[`compute_junction`](../engine_rs/src/junction.rs#L83) using
**germline-position lookup** (not offset arithmetic), so it
remains correct under structural indels:

```
v_anchor_pool = anchor_pool_position(sim, V, v_anchor_offset_in_allele)
j_anchor_pool = anchor_pool_position(sim, J, j_anchor_offset_in_allele)

start  = v_anchor_pool                 // first base of V Cys codon
end    = j_anchor_pool + 3             // one past J W/F codon
length = end - start                   // = byte count when end >= start
```

The window is closed-open `[start, end)`; `junction_length` after
slicing is `end - start`.

### 2.1 Returns `None` (junction undefined) when

- V or J assignment missing (pre-recombination).
- V or J allele not in refdata (signature mismatch).
- V or J anchorless allele.
- `v_trim_5 > v_anchor` or `j_trim_5 > j_anchor` (trim ate past
  the anchor codon).
- V or J region not yet assembled.
- `end <= start` (malformed simulation).

The AIRR builder gracefully handles `None`: the junction fields
either stay `None` (when explicit `Optional[...]`) or get defaults
(`junction` empty string, `junction_length` 0, etc.). The
`productive` triad is reported as `None` in that case.

---

## 3. Translation — `junction_aa`

`junction_aa` is computed via
[`translate_seq(junction)`](../engine_rs/src/codon.rs#L83):

- Codons read in triplets from `junction[0]` — frame is **always
  zero-aligned to V_anchor**.
- Each codon → amino acid via `GENETIC_CODE[64]` lookup.
- Non-canonical bases (anything other than A/C/G/T/U,
  case-insensitive) → `X` (ambiguous).
- Lowercase letters are upper-cased implicitly via the lookup table
  (mirroring the [allele-call audit](allele_call_audit.md) §1.3
  rule: lowercase is presentation only, not wildcard).
- Stop codons → `*`.

### 3.1 Out-of-frame: empty string

When `junction_length % 3 != 0`, `junction_aa` is the empty
string `""` (NOT `None`). The translator skips the call rather
than returning partial codons. Pinned by
`test_out_of_frame_junction_emits_empty_aa_and_false_stop`.

Consequence: users filtering "records with translation" must check
`vj_in_frame is True` rather than `junction_aa != ""` — the empty
string survives `bool` truthiness checks but isn't a translation.

### 3.2 Length invariant when in-frame

When in-frame, `len(junction) / 3 == len(junction_aa)`. Pinned by
`test_junction_translation_invariant_codon_to_aa`.

---

## 4. The productive triad

`productive = in_frame ∧ no_junction_stop ∧ anchors_preserved`,
evaluated at AIRR projection time (post-pipeline). The three
predicates are independent:

| `vj_in_frame` | `stop_codon` | anchors preserved | `productive` |
|---------------|--------------|-------------------|--------------|
| `False`       | `False` (vacuous) | n/a            | `False`      |
| `True`        | `True`       | n/a               | `False`      |
| `True`        | `False`      | `False`           | `False`      |
| `True`        | `False`      | `True`            | `True`       |

### 4.1 `stop_codon` is `False` when out-of-frame

The semantically honest answer for "is there a stop codon in this
junction" when the junction isn't in-frame is **undefined** — but
the field is typed as `Optional[bool]` and reports `False`. The
audit's call is that this is a documentation gap, not a bug: the
field name suggests "True iff there's a stop codon in the
junction frame," and out-of-frame junctions have no frame so the
answer is vacuously False.

Pinned by `test_stop_codon_is_false_when_out_of_frame`.

### 4.2 Anchor preservation under SHM

Mutation inside the V or J anchor codon that changes the
translated amino acid drops `productive` to `False` even when
the junction is in-frame and contains no stop. The
[`anchor_amino_acid_preserved`](../engine_rs/src/airr_record/junction.rs)
check translates the current anchor codon and compares to the
reference allele's anchor codon's amino acid.

Pinned by `test_v_anchor_mutation_to_non_cys_breaks_productive_only`.

---

## 5. Behaviour under evidence-changing mechanisms

### 5.1 NP segments

- **NP1 = 0 (VJ)**: junction = V_anchor_codon + J_anchor_codon
  (the two anchors touch). Length = 6, two codons (Cys + W/F).
  Productive.
- **NP1 = 3 (VJ in-frame)**: junction = V_anchor + NP1(3) + J_anchor.
  Length = 9, three codons.
- **VDJ**: junction spans V_anchor → NP1 → D → NP2 → J_anchor.
  Length = `V_anchor_to_end + NP1 + D_len + NP2 + (J_anchor + 3)`.

Pinned by `test_vj_zero_np_junction_collapses_to_anchor_codons` and
`test_vdj_junction_spans_v_np1_d_np2_j`.

### 5.2 Recombination trims (`v_3`, `d_5/3`, `j_5`)

V/J anchor-side trims that *preserve* the anchor codon shift the
non-anchor end of the V/J region but don't change `junction_start`
/ `junction_end` arithmetic — the anchor's germline position scan
still finds it. The junction window is unchanged in content;
the V/J region's coordinates outside the junction shift.

D and NP trims affect junction length directly (they live inside
the junction window).

Pinned by `test_v_trim_3_inside_anchor_zone_leaves_junction_unchanged`.

### 5.3 Indels

Indels are classified per-site by
[`ProductiveJunctionFrame::admissible_indel_class_at`](../engine_rs/src/contract/productive_junction_frame.rs#L231):

| Indel site relative to V/J regions       | Junction effect                                        |
|------------------------------------------|--------------------------------------------------------|
| `s < V_region.start`                     | Both V and J region starts shift → junction window shifts as a whole; **length unchanged**. |
| `V_region.start ≤ s < J_region.start`    | V stays, J shifts → junction window grows / shrinks → **length ± 1**. |
| `s ≥ J_region.start`                     | Neither shifts → junction window unchanged; **length unchanged**. |

A single indel in the "between V and J" zone flips
`vj_in_frame` to `False` (junction length becomes 9±1). Two
indels (insertion + deletion in the middle zone) preserve
length. Three indels in the middle zone shift by ±3 = 0 mod 3 →
in-frame again.

This is the rule the productive contract uses to filter
candidates; it's also the rule the AIRR projection observes after
the fact.

Pinned by `test_single_indel_in_v_anchor_zone_breaks_in_frame` and
`test_indel_before_v_anchor_shifts_window_but_preserves_length`.

### 5.4 End-loss (`end_loss_5prime` / `end_loss_3prime`)

End-loss removes bytes from the pool endpoints **after** all
sampling passes have run. The junction window itself is unchanged
content-wise (the V and J anchor codons survive end-loss because
they're not at the pool boundary in any normal pipeline), but
`junction_start` and `junction_end` shift to reflect the new pool
coordinates.

(`primer_trim_5prime` / `primer_trim_3prime` are backwards-compat
aliases for the canonical `end_loss_*prime` methods — same trace
address, same AIRR field, byte-identical records.)

Example: baseline junction at `[6, 15)` in an 18-base sequence.
After `end_loss_5prime(length=2)` removes the leading 2 bases:
the post-trim sequence is 16 bases; junction at `[4, 13)`. The
junction *content* is identical (`TGTACATGG`), only the
coordinates shifted.

Pinned by `test_end_loss_5_prime_shifts_junction_coordinates_not_content`.

### 5.5 Reverse-complement projection (`random_strand_orientation`)

When `rev_comp=True`, the AIRR projector calls
[`apply_rev_comp_projection`](../engine_rs/src/airr_record/sequence.rs#L14):

- Sequences (`sequence`, `np1`, `np2`, `junction`) are
  reverse-complemented in place.
- All `*_start / *_end` coordinate pairs are flipped:
  `new_start = seq_len - old_end`, `new_end = seq_len - old_start`.
- AA fields (`junction_aa`, `sequence_aa`, `np1_aa`, `np2_aa`)
  are **re-translated** from the reverse-complemented sequences,
  with the codon frame recomputed as `(junction_start % 3)`.

Consequence: `junction_aa` under rev-comp is NOT the reverse of
the original `junction_aa` — it's the translation of the
reverse-complement of the original `junction`. Empirically,
`"TGTACATGG" → "CTW"` becomes `"CCATGTACA" → "PCT"` after rev-comp,
not `"WTC"`.

Pinned by `test_rev_comp_re_translates_junction_aa_from_rc_string`.

---

## 6. Drift identified

### 6.1 `junction_aa = ""` (not `None`) when out-of-frame

The empty-string convention for out-of-frame junctions is
defensible (translation didn't fail, it just had no codons), but
diverges from the AIRR-standard recommendation of leaving the
field unset / null when the junction isn't in-frame. Users
filtering `df[df["junction_aa"] != ""]` get a different set than
`df[df["vj_in_frame"] == True]` when junction is truly undefined
(both `""` AND `vj_in_frame=None`).

**Possible fix:** return `None` for `junction_aa` when
`vj_in_frame` is `False` or `None`. Breaking change for callers
that test `len(junction_aa)`.

**Tradeoff:** the Rust kernel populates a `String` field, so
`None` would require an `Option<String>` schema bump. Defer
unless a downstream consumer reports the divergence.

### 6.2 `cdr3` / `cdr3_aa` not exposed as aliases

The AIRR spec distinguishes `junction` and `cdr3` (CDR3 = junction
minus the two anchor codons). GenAIRR exposes only `junction` —
users who want IMGT-style CDR3 need to slice `junction[3:-3]` and
re-translate.

**Possible fix:** add `cdr3 = junction[3:-3]` and
`cdr3_aa = junction_aa[1:-1]` as denormalised fields, populated
when `vj_in_frame` is `True`.

**Tradeoff:** two more AIRR fields; the math is trivial in
Python. Could be deferred indefinitely.

### 6.3 `stop_codon` is `False` (not `None`) when out-of-frame

Same flavour as §6.1 — the field reads `False` when the question
is ill-defined. A consumer counting "non-productive records by
reason" must use the triad `(vj_in_frame, stop_codon, productive)`
together rather than `stop_codon` alone.

**Possible fix:** return `None` when `vj_in_frame is False`.

### 6.4 No per-record reason for `productive=False`

When `productive=False`, the user must reconstruct the cause by
combining `vj_in_frame`, `stop_codon`, and (implicitly) "anchor
mutated" (which has no explicit boolean). A
`productive_fail_reason` field — one of `"out_of_frame"`,
`"junction_stop_codon"`, `"v_anchor_broken"`, `"j_anchor_broken"`,
`"junction_undefined"` — would denormalise the diagnosis.

**Possible fix:** add an `Option<String>` field populated when
`productive=False` (or always populated with `"productive"` /
the failure reason).

**Tradeoff:** denormalisation; the audit's failure matrix
(§4) already gives users the recipe to compute this in Python.

---

## 7. Test coverage in this slice

[`tests/test_junction_call_provenance.py`](../tests/test_junction_call_provenance.py)
pins the §2–§5 rules:

- **Field exposure (§1)**: assert `junction`, `junction_aa`,
  `junction_length`, `junction_start`, `junction_end`,
  `vj_in_frame`, `stop_codon`, `productive` are on every record;
  pin absence of `cdr3` / `cdr3_aa` / `v_anchor_codon_pos`.
- **Junction window boundaries (§2)**: VJ baseline at `[V_anchor,
  J_anchor + 3)`, junction content matches the V/NP/J slice.
- **Translation (§3)**: codon-by-codon mapping, length invariant
  when in-frame, empty string when out-of-frame.
- **Productive triad (§4)**: each cause of `productive=False`
  isolated: out-of-frame, stop codon, V-anchor mutation,
  J-anchor mutation.
- **NP variations (§5.1)**: NP1=0 baseline, NP1=3, VDJ with D and
  NP2.
- **Indels (§5.3)**: indel in V_anchor zone breaks frame; indel
  before V_anchor preserves length.
- **End-loss (§5.4)**: junction coordinates shift, content
  unchanged.
- **Rev-comp (§5.5)**: junction reverse-complemented, coordinates
  flipped, AA re-translated.
- **Replay round-trip**: every junction field reproduces under
  trace replay.
- **Pin-the-drift (§6)**: assert current state of empty-string
  AA, missing `cdr3` field, missing fail-reason.
