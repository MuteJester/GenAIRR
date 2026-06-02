# Primer-trim / end-loss provenance audit

**Status:** audit + golden tests. No code changes proposed in this
slice — drift items are listed in §6 for a follow-up.

This audit catalogues how the engine currently handles
recombination trim, primer trim, and end-loss, what AIRR fields
each one affects, and where biological provenance is collapsed.
The companion test file
[`tests/test_primer_trim_end_loss_provenance.py`](../tests/test_primer_trim_end_loss_provenance.py)
pins current behaviour so a future fix doesn't silently change
observable semantics.

---

## 1. What's there today

Three nominally distinct mechanisms exist in the user-facing
vocabulary:

| User-facing surface                                       | Biological intent                              |
|-----------------------------------------------------------|------------------------------------------------|
| `TrimPass` via `Experiment.trim(...)`                     | Recombination exonuclease trim (V/D/J pre-NP). |
| `EndLossPass` via `Experiment.end_loss_5prime(...)`       | Observation-stage 5' read-end / primer loss.   |
| `EndLossPass` via `Experiment.end_loss_3prime(...)`       | Observation-stage 3' read-end / primer loss.   |

`Experiment.primer_trim_5prime(...)` / `primer_trim_3prime(...)`
are backwards-compatible aliases for the `end_loss_*prime`
methods — same trace addresses, same AIRR fields, byte-identical
records. New code should prefer the canonical `end_loss_*` names.

The DSL also exposes `EndLossPass` implicitly through corruption
flows that the species datasets configure.

### 1.1 Engine implementation

At the engine layer there are **two** distinct mechanisms, not three:

- **Recombination trim** — [`TrimPass`](../engine_rs/src/passes/trim.rs).
  Updates the allele-instance `trim_5` / `trim_3` fields on the
  simulation. *Does not* delete pool bytes (the bytes never get
  pushed in the first place because assembly reads the post-trim
  slice).
  - Trace addresses: `trim.{v,d,j}_{5,3}`.

- **End-loss / primer-trim** — [`EndLossPass`](../engine_rs/src/passes/corrupt/end_loss.rs).
  Physically deletes pool bytes from the 5' or 3' end via
  [`MutationTransaction::delete_base_admitting`]. *Does not*
  update `trim_5` / `trim_3` allele-instance metadata.
  - Trace addresses: `corrupt.end_loss.{5,3}`.

**Finding 1.** The DSL's `end_loss_5prime()` /
`end_loss_3prime()` are thin wrappers over the same
`EndLossPass`. There is no separate "primer-trim" pass in the
engine — the legacy `primer_trim_*prime` names that predate the
engine's vocabulary survive as backwards-compatible aliases for
the canonical `end_loss_*prime` methods. Both spellings map to
the same code path, the same trace addresses, and the same
emitted events.

### 1.2 Engine-level semantics already correct

The end-loss module docstring spells out the design intent and
the engine implements it:

> End-loss is the *observation-stage* primer / read-degradation
> artifact. It physically deletes pool bytes via
> `MutationTransaction::delete_base_admitting`, which routes
> through the bundle's `admits_post_event`. It does NOT update
> allele `trim_5` / `trim_3` metadata — that's the
> recombination-stage `TrimPass`'s job.

So at the **engine** level, the two biological provenances are
not conflated: recombination trim writes to allele metadata, end-
loss deletes pool bytes. Anchor preservation under
`productive_only` enforces that end-loss can't bite into V/J
anchor codons. ✅

---

## 2. AIRR fields affected

What a downstream analyst sees in the projected AIRR record:

| AIRR field                              | Recombination trim drives it          | End-loss / primer-trim drives it       |
|-----------------------------------------|----------------------------------------|-----------------------------------------|
| `sequence`, `sequence_aa`               | Yes (post-trim allele slice).          | Yes (further shortened post-pool-delete).|
| `v_trim_5`, `v_trim_3`, `j_trim_5`, `j_trim_3` | **Yes** (sole writer).            | **No** (intentionally — see Finding 1). |
| `d_trim_5`, `d_trim_3`                  | **Yes** (sole writer).                 | No (end-loss is 5'/3' of the whole pool, not per-segment trim). |
| `v_sequence_start`, `v_sequence_end`    | Indirect (V region position).         | Indirect (V region shifts in pool).    |
| `j_sequence_start`, `j_sequence_end`    | Indirect.                              | Indirect.                                |
| `v_germline_start`, `v_germline_end`    | Indirect via `v_trim_5/3`.            | **Indirect via live-call walker only.** |
| `j_germline_start`, `j_germline_end`    | Indirect via `j_trim_5/3`.            | **Indirect via live-call walker only.** |
| `n_indels`                              | No.                                    | **No** — `n_indels` comes from `corrupt.indel.count` (the `IndelPass`), not from end-loss. |
| `n_pcr_errors`, `n_quality_errors`      | No.                                    | No.                                      |
| `productive`, `vj_in_frame`, `stop_codon` | Indirect (via junction shape).      | Indirect (anchor preservation enforced under `productive_only`). |

### 2.1 Walker vs trim-only fallback for germline coords

[`engine_rs/src/airr_record/builder.rs:335-355`](../engine_rs/src/airr_record/builder.rs):

```rust
if let Some((g_start, g_end)) = walk.ref_ranges[0] {
    rec.v_germline_start = Some(g_start);          // walker-derived
    rec.v_germline_end = Some(g_end);
} else if let Some(allele) = lookup_allele(refdata, Segment::V, v_id) {
    rec.v_germline_start = Some(rec.v_trim_5);     // trim-only fallback
    rec.v_germline_end = Some(allele.seq.len() as i64 - rec.v_trim_3);
}
```

- **Live-call walker path (production)**: `v_germline_start/end`
  reflect actual post-end-loss reference coverage. End-loss IS
  surfaced.
- **Trim-only fallback (no allele assignment, raw RefDataConfig
  runs without live calls)**: `v_germline_start/end` are derived
  from `v_trim_5/3` only and miss any end-loss contribution.

In typical user-facing pipelines the walker fires, so germline
coords are end-loss-aware. The fallback is exercised by
synthetic test fixtures and isolated unit tests.

---

## 3. AIRR fields NOT affected — the provenance gap

**Finding 2.** No AIRR field directly carries the end-loss /
primer-trim *amount*. The trace records it (at
`corrupt.end_loss.{5,3}`), but the projection drops it.

A downstream analyst comparing two records cannot, from the AIRR
output alone, tell:

- Record A: `v_trim_3 = 5`, no end-loss.
- Record B: `v_trim_3 = 0`, end-loss-3' removed 5 bases of V.

The sequence and `v_sequence_end` look identical between A and
B. The recombination-trim fields differentiate them in this
example *only because* end-loss doesn't touch the trim fields —
but there's no positive signal "primer-trim shortened me by 5".

The user can recover the end-loss amount in three ways, none
ergonomic:

1. Read the `outcome.trace()` and look up
   `corrupt.end_loss.{5,3}` records directly.
2. Compare `v_germline_start` (walker-derived) to `v_trim_5`:
   the difference equals the end-loss into V (when end-loss bit
   into V).
3. Compare the assembled-before-corruption sequence length to
   the final sequence length (requires intermediate revisions).

None of these are AIRR-format-portable.

---

## 4. Expected biological semantics

The architecture distinguishes three nominally separable
mechanisms. The cleanest target state is:

| Mechanism             | Biological stage          | Engine pass         | AIRR signal                                 |
|-----------------------|---------------------------|---------------------|---------------------------------------------|
| Recombination trim    | V(D)J recombination       | `TrimPass`          | `v_trim_5/3`, `d_trim_5/3`, `j_trim_5/3`    |
| Primer trim           | Library preparation       | (separate pass)     | `primer_trim_5_length`, `primer_trim_3_length` |
| End-loss              | Sequencing read degradation | (separate pass)   | `n_end_loss_5`, `n_end_loss_3` (or similar) |

In practice the user-facing biology rarely distinguishes primer
trim from end-loss — both are read-stage length losses with
similar statistical models. Whether to actually split them
(versus keeping them merged as today and just renaming the
exposure) is a product question, not an architectural one.

The **minimum** target state from the audit's perspective:

- One observation-stage shortening mechanism (merged primer-trim /
  end-loss is fine — that's what the engine does today).
- A dedicated AIRR field for the amount removed at each end, so
  the provenance reaches downstream analysts.
- `v_germline_start/end` and `j_germline_start/end` consistent
  between the walker-derived and trim-only fallback paths.

---

## 5. Productive-contract interaction (already correct)

Under `productive_only`, [`AnchorPreserved`](../engine_rs/src/contract/anchor_preserved.rs)
narrows the end-loss length distribution so the loss can't bite
into the V-anchor codon (5' loss) or the J-anchor codon (3'
loss). The end-loss pass clamps to the admissible support in
permissive mode and surfaces a structured error in strict mode.

This is the right semantics: observation-stage length losses
respect biological constraints carried over from recombination-
stage choices. The productive-contract stress matrix already
exercises primer-trim under `productive_only` and confirms every
record stays productive. ✅

---

## 6. Drift identified

### 6.1 Provenance gap on AIRR record — **✅ Resolved**

End-loss / primer-trim amounts were recorded in the trace but
never surfaced to the projected AIRR record.

**Fix landed:** added `end_loss_5_length` and `end_loss_3_length`
fields to [`AirrRecord`](../engine_rs/src/airr_record/record.rs),
populated from `corrupt.end_loss.{5,3}` trace records in
[`builder.rs`](../engine_rs/src/airr_record/builder.rs), mirrored
through the Rust→Python dict projection in
[`python/outcome.rs`](../engine_rs/src/python/outcome.rs), and
added to the canonical column list in
[`result.py`](../src/GenAIRR/result.py). Counters default to `0`
when no end-loss pass ran. Test coverage:
`test_end_loss_5_length_field_equals_trace_value`,
`test_end_loss_3_length_field_equals_trace_value`,
`test_end_loss_fields_default_to_zero_when_no_end_loss_pass_ran`,
and `test_combined_end_loss_fields_match_per_end_trace_values` in
[`test_primer_trim_end_loss_provenance.py`](../tests/test_primer_trim_end_loss_provenance.py).

### 6.2 DSL surface naming — **✅ Resolved**

The canonical DSL methods are now :meth:`Experiment.end_loss_5prime`
and :meth:`Experiment.end_loss_3prime`, matching the engine pass
(`EndLossPass`), the trace addresses (`corrupt.end_loss.{5,3}`),
and the AIRR fields (`end_loss_{5,3}_length`).

`primer_trim_5prime()` / `primer_trim_3prime()` survive as
backwards-compatible aliases for the canonical methods. Same
trace address, same AIRR field, byte-identical records for the
same seed — pinned by `tests/test_end_loss_dsl_alias.py`.

### 6.3 Trim-only fallback for germline coords — **✅ No drift (parity invariant pinned)**

**Original audit hypothesis:** the fallback at
[`builder.rs:335-355`](../engine_rs/src/airr_record/builder.rs)
ignores end-loss when computing `{v,d,j}_germline_start/end`, so
a fallback-firing record might disagree with what the walker
would have produced.

**Investigation outcome:** the walker encodes a 5'-end-loss as
leading `D` (deletion) ops in V's CIGAR, keeping
`v_germline_start = v_trim_5` regardless of how many pool bytes
end-loss removed. The 3' J side mirrors. This is the canonical
AIRR convention: `germline_start/end` describe the alignment's
allele-coordinate span, not the read-coordinate coverage.

The fallback's existing trim-only math (`v_germline_start =
v_trim_5`, `v_germline_end = allele_len - v_trim_3`) **already
matches that semantic**. A shift-by-`end_loss_length`
formulation would diverge from the walker and break AIRR-tool
parity. No code change required.

**Test coverage pinning the parity:**
- [`fallback_v_germline_coords_unaffected_by_end_loss_match_walker`](../engine_rs/src/airr_record/tests/projection.rs)
- [`fallback_j_germline_coords_unaffected_by_end_loss_match_walker`](../engine_rs/src/airr_record/tests/projection.rs)
- [`walker_path_keeps_v_germline_start_at_zero_under_5prime_end_loss`](../engine_rs/src/airr_record/tests/projection.rs)
- [`fallback_and_walker_paths_agree_on_germline_coords_under_end_loss`](../engine_rs/src/airr_record/tests/projection.rs)

Together these pin: walker emits `germline_start=0` under
5'-end-loss; fallback emits the same; both surface the end-loss
amount through the dedicated `end_loss_5_length` provenance
field (audit §6.1 — the *actual* drift fix).

### 6.4 Unified counter for observation-stage length losses — **✅ Resolved (merged with §6.1)**

Resolved together with §6.1. The `end_loss_5_length` and
`end_loss_3_length` fields are the canonical counters for the
observation-stage `EndLossPass` mechanism — explicit names that
match the engine pass and the trace address (no
`primer_trim_*` aliasing in the AIRR surface, even though the
DSL still exposes `primer_trim_*` methods per §6.2 deferral).

---

## 7. Test coverage in this slice

[`tests/test_primer_trim_end_loss_provenance.py`](../tests/test_primer_trim_end_loss_provenance.py)
pins the **current** behaviour with golden tests, including the
known drift items above. When a future slice fixes the drift,
the relevant assertions will need to be flipped from "documents
the gap" to "asserts the fix" — that flip is itself a useful
signal that the fix touched the right surface.

The golden suite covers:

- 5'-loss shortens the sequence by exactly the requested length.
- 3'-loss shortens the sequence by exactly the requested length.
- Recombination-trim fields (`v_trim_5/3`, `j_trim_5/3`) are
  untouched by end-loss.
- The trace carries the `corrupt.end_loss.{5,3}` records with
  the realized count.
- `n_indels` is unchanged by end-loss (separate biology).
- Replay through `replay_from_trace_file` reproduces sequence
  and AIRR coordinates record-for-record.
- (Pin-the-drift) No `end_loss_5_length` / `n_end_loss_*` field
  on the AIRR record — when one is added, this test fails and
  the fix-slice author knows to update it.

The audit doc and the test file are designed to be read
together: doc explains why, tests show what.
