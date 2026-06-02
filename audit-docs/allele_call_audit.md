# Allele-call disambiguation audit

**Status:** audit + golden tests. Drift items in §6.

This audit catalogues how the engine derives per-segment allele
*live calls* (`v_call`, `d_call`, `j_call`) from the assembled
pool, how those calls change under evidence-altering mechanisms
(SHM, PCR errors, sequencing errors, N substitution, indels), and
how truth (`truth_v_call`, etc., under `expose_provenance=True`)
stays decoupled from the evidence-driven call.

Pattern: same audit-first shape as
[primer_trim_end_loss_audit.md](primer_trim_end_loss_audit.md) and
[indel_provenance_audit.md](indel_provenance_audit.md) — document
current behaviour, pin with golden tests, catalogue drift, defer
fixes to follow-up slices.

---

## 1. What's there today

### 1.1 Live-call model

The per-segment "set of alleles consistent with current evidence"
lives on the `Simulation` IR as a
[`SegmentLiveCall`](../engine_rs/src/live_call/model.rs)
per assignable segment (V / D / J). Each `SegmentLiveCall` holds:

- `allele_call: AlleleBitSet` — dense bitset of the alleles
  currently in the *tie-set at maximum score* (the live call).
- `hypotheses: Vec<PlacementHypothesis>` — per-hypothesis
  breakdown when multiple placements coexist.
- `confidence: LiveCallConfidence` — one of `Unresolved`,
  `Unsupported`, `ExactSingle`, `ExactAmbiguous`.
- `evidence_version: u64` — monotonic version for change
  detection across passes.

The bitset is *dense* and indexed by allele ID — every allele in
the refdata's segment pool has a bit; "in the live call" means
"its bit is set."

### 1.2 Walker logic (per-position match counting)

[`walker_observer`](../engine_rs/src/live_call/walker/mod.rs)
implements a per-position match counter per allele:

- For each observed pool byte at reference position `ref_pos`,
  the walker consults [`ReferenceIndex::compatible_alleles_at`](../engine_rs/src/live_call/reference_index.rs)
  to find which alleles have a matching germline byte at that
  reference position.
- Each matching allele's `scores[i]` increments by 1.
- After the segment's region is walked, the alleles at the
  **maximum** score form the `allele_call` tie-set. Sub-maximal
  alleles are dropped.

The walker also tracks the *informative* / *uninformative* status
of each column for diagnostics — a column where every allele
matches (N wildcard) is uninformative and doesn't break ties even
though it increments every allele's score.

### 1.3 Match semantics — the rules

[`canonical_base_index`](../engine_rs/src/live_call/reference_index.rs#L357)
and [`observed_base_kind`](../engine_rs/src/live_call/reference_index.rs#L351)
specify the per-byte match rule:

| Observed byte             | Treated as       | Match behaviour                                                  |
|---------------------------|------------------|------------------------------------------------------------------|
| `A` / `C` / `G` / `T`     | Canonical        | Matches **only** alleles whose germline byte equals it (case-folded). |
| `a` / `c` / `g` / `t`     | Canonical        | Lowercase is **upper-cased** before lookup — same as the uppercase letter. |
| `N` / `n`                 | Wildcard         | Matches **every** allele at that position (informative=false).   |
| Anything else             | Not a base       | `None` — column contributes no evidence.                         |

**Key consequence:** the lowercase convention used by
[`sequencing_errors`](../src/GenAIRR/experiment.py#L216) is a
*presentation* marker (so a human reader can tell which bases
were corrupted), **not** an evidence semantic. A lowercase `t`
substituted onto a germline `C` is treated by the walker exactly
the same as an uppercase `T`: it's a substitution to a different
canonical base, and it can flip a call.

### 1.4 AIRR projection — `v_call` / `d_call` / `j_call`

[`projected_call_name`](../engine_rs/src/airr_record/projection.rs#L73)
and
[`live_call_name`](../engine_rs/src/airr_record/projection.rs#L207)
turn the tie-set into a comma-separated string with a
**truth-preference order**:

1. If the truth allele ID (the one originally sampled in the
   recombination pass) is in the tie-set, list it **first**.
2. Append the remaining tied allele IDs in ascending order.
3. If the tie-set is empty (live call `Unsupported` / unresolved),
   fall back to the truth allele name verbatim.

Examples (V tie-set):
- `{v1*01}` truth=v1*01 → `"v1*01"`
- `{v1*01, v1*02}` truth=v1*01 → `"v1*01,v1*02"`
- `{v1*01, v1*02}` truth=v1*02 → `"v1*02,v1*01"`
- `{v_B*01, v_C*01}` truth=v_A*01 (dropped) → `"v_B*01,v_C*01"`

### 1.5 Identity fields

`v_identity` / `d_identity` / `j_identity` are computed by
[`walk_alignment_columns`](../engine_rs/src/airr_record/walk/mod.rs)
as `matches / total` aligning the assembled pool against the
**projected (first-in-CSV)** allele's germline — not against
truth. So when a mutation flips the call from truth to a different
allele, identity is computed against the *new* call, which is why
the projected call can show `v_identity = 1.0` even when truth
has changed (the engine's interpretation of the evidence is
self-consistent).

### 1.6 Dirty windows + refresh

[`DirtySignalObserver`](../engine_rs/src/live_call/dirty_signal_observer.rs)
records per-position dirty windows from runtime events:

- `BaseChanged` → `DirtyWindow { site, reason: BaseEdited }`
- `IndelInserted` / `IndelDeleted` →
  `DirtyWindow { site, reason: StructuralIndel { delta: ±1 } }`

[`LiveCallRefreshPlan`](../engine_rs/src/live_call/refresh_plan.rs)
maps event records to refresh steps:

- `BaseChanged` only → `EditedSegments` (dirty-window-narrowed
  re-evaluation).
- Any `IndelInserted` / `IndelDeleted` → `AllStructural` (full
  V/D/J refresh; subsumes the edit-narrowed path).

The refresh runs after each pass via
[`LiveCallRefreshHook`](../engine_rs/src/live_call/refresh_hook.rs).

---

## 2. AIRR fields affected

| AIRR field                                         | Population path                                                     |
|----------------------------------------------------|---------------------------------------------------------------------|
| `v_call` / `d_call` / `j_call`                     | Comma-separated tie-set (truth-first then asc-id).                  |
| `v_identity` / `d_identity` / `j_identity`         | Matches against the **projected** allele (first in CSV).            |
| `v_cigar` / `d_cigar` / `j_cigar`                  | Walker per-column M/X/I/D ops vs the projected allele.              |
| `v_germline_start` / `..end` (etc.)                | Coordinates in the projected allele's reference space.              |
| `v_sequence_start` / `..end` (etc.)                | Region coordinates in the assembled pool.                           |
| `truth_v_call` / `truth_d_call` / `truth_j_call`   | When `expose_provenance=True`: the originally-sampled allele name.  |
| `n_mutations`                                      | Distinct from PCR / sequencing error counters; counts SHM only.     |
| `n_pcr_errors` / `n_quality_errors`                | Per-mechanism error counts; do NOT shift `n_mutations`.             |

---

## 3. Behaviour under evidence-changing mechanisms

This is the load-bearing section. Each row pins what should
happen when the listed mechanism touches a distinguishing /
non-distinguishing / N position.

### 3.1 SHM (`mutate`)

| Position kind                  | Tie-set effect                       | Truth field                |
|--------------------------------|--------------------------------------|----------------------------|
| Non-distinguishing             | Unchanged (still same tie-set).      | Unchanged (truth retained).|
| Distinguishing → matches another allele | Call **switches** to the new allele; truth drops out. | Unchanged (still original). |
| Distinguishing → matches nothing | Score decreases for truth; depending on other alleles, may switch or stay. | Unchanged.        |

The walker has no truth bias — once a mutation creates evidence
against the truth allele, the truth is dropped from the tie-set
exactly the same as any other allele. `truth_v_call` (under
`expose_provenance=True`) remains the only fact that survives.

### 3.2 PCR errors (`pcr_amplify`)

Substitutes a uniform random canonical base (uppercase) at a
random position. **Semantically identical to SHM** from the
walker's perspective — both write a new canonical base into the
pool. The distinction is counter-only: PCR increments
`n_pcr_errors`, SHM increments `n_mutations`.

### 3.3 Sequencing errors (`sequencing_errors`)

Substitutes a uniform random canonical base in *lowercase*
(presentation marker). **Walker treats lowercase as uppercase**,
so semantically identical to PCR / SHM. Counter is `n_quality_errors`.

### 3.4 N substitution (`ambiguous_base_calls`)

Substitutes the literal byte `N` at a random position. The walker
treats N as a wildcard:

- At a **non-distinguishing** position: every allele gets a +1
  from the wildcard match; tie-set composition unchanged
  (everyone in the tie keeps tying).
- At a **distinguishing** position: the position no longer
  discriminates; alleles that were tied minus the truth now also
  match → tie-set **widens**, truth retained.

### 3.5 Indels (`polymerase_indels`)

Trigger `AllStructural` refresh — every segment's live call is
recomputed against the post-indel pool. Insertions inside a region
emit `I` CIGAR ops (no germline credit); deletions trigger `D`
ops. The tie-set effect depends on whether the indel hits a
distinguishing position (which it usually doesn't, since indels
are independent of position bias).

See [`indel_provenance_audit.md`](indel_provenance_audit.md) for
how indel events flow through `n_indels` / `n_v_indels` /
`n_d_indels` / `n_j_indels`.

---

## 4. Replay determinism

The trace records every RNG-bearing choice (allele id at
recombine, mutation sites and bases, error sites and bases, indel
kind/site/base, etc.). The walker is deterministic given the
final pool state. Therefore:

- `v_call` / `d_call` / `j_call` reproduce exactly under replay.
- `v_identity` etc. reproduce exactly.
- `truth_v_call` reproduces (it's just the recombine allele
  trace).
- Cross-mechanism interaction (SHM + PCR errors on the same
  record) reproduces.

Pinned by the parametrised
`test_allele_call_replay_round_trip` in
[`tests/test_allele_call_provenance.py`](../tests/test_allele_call_provenance.py).

---

## 5. Expected biological semantics

| Question                                                                | Answer                                                              |
|-------------------------------------------------------------------------|---------------------------------------------------------------------|
| Does the walker pick a single "best" allele under ambiguity?            | **No.** It returns the full tie-set at the max evidence score.       |
| Can a mutation move the call **away** from truth?                       | **Yes** — at a distinguishing position. Truth survives in `truth_v_call` (when exposed). |
| Does a sequencing-error lowercase letter count as wildcard / N?         | **No.** Lowercase is upper-cased before lookup; only literal `N` / `n` is wildcard. |
| Does an N substitution always widen the tie-set?                        | **Only if** the N lands at a distinguishing position. At a non-distinguishing position, tie-set composition is unchanged (all current ties stay tied). |
| Is `v_identity` computed against truth or against the called allele?    | **Against the called (projected) allele** — the first entry in `v_call`. |
| Does the CSV order in `v_call` carry semantic meaning?                  | **Yes.** First entry = the projected allele (truth if in tie-set, otherwise lowest-id tied). |

---

## 6. Drift identified

Items below are **not fixed in this slice**. Catalogued for
follow-up scoping.

### 6.1 No `*_call_count` ambiguity counter on the AIRR record

The CSV in `v_call` carries the ambiguity, but a user filtering
"records with unambiguous V call" must `len(v_call.split(","))`
the string. A denormalised `v_call_count` / `d_call_count` /
`j_call_count` would be cheaper for SQL / dataframe filtering.

**Possible fix:** add three i64 fields populated from
`live_call.allele_call.population_count()`.

**Cost:** three new AIRR fields. Minor schema expansion. Could be
deferred indefinitely — `len(call.split(","))` works.

### 6.2 No explicit "called allele ≠ truth allele" indicator

When a mutation flips the call away from truth, the divergence
is recoverable only by users who request `expose_provenance=True`
and compare `v_call.split(",")[0] != truth_v_call`. A boolean
`v_call_truthful` (or similar) would denormalise the comparison
for SQL filters and reduce the chance of users silently treating
`v_call` as truth.

**Possible fix:** add three booleans (`v_call_is_truthful`,
`d_call_is_truthful`, `j_call_is_truthful`) populated when
`expose_provenance=True`. Or always populate, with `None` when
truth isn't recoverable.

**Cost:** three new fields, depends on default vs opt-in.
Product question — surface noise vs. helpful safety rail.

### 6.3 `identity` against the called allele can be misleading
under switched calls

When the call switches away from truth (e.g. mutation makes the
sequence "look like" v2*01 even though truth=v1*01),
`v_identity` is computed against v2*01 — and can read very high
(even 1.0 if all mutations conspire to match v2*01 exactly).
This is mathematically correct (identity against the *called*
allele) but can mask the fact that the original sample was
v1*01.

**Possible fix:** add `v_identity_to_truth` populated when
`expose_provenance=True`. Or document the convention in the
projection docstring more prominently.

**Cost:** symmetric to §6.2 — three opt-in fields. Could be
denormalised via post-hoc Python.

### 6.4 No per-pass live-call evolution exposed

The live call's `evidence_version` bumps each time a pass mutates
the IR — a benchmarking user might want "the tie-set just after
recombine" vs. "after SHM" to study how each mechanism collapses
ambiguity. Today the trajectory isn't surfaced; only the final
state is projected to AIRR.

**Possible fix:** add a Python helper on `Outcome` that returns
the live-call sequence over the revisions stack. No AIRR-record
change.

**Cost:** small Python surface addition. Genuinely defer until a
benchmarking user asks.

---

## 7. Test coverage in this slice

[`tests/test_allele_call_provenance.py`](../tests/test_allele_call_provenance.py)
pins the current behaviour with golden tests. Coverage:

- **Identical alleles tie**: `v_call` is the CSV of all identical
  alleles with truth first; identity = 1.0.
- **Distinguishing alleles, no mutation**: `v_call` equals
  `truth_v_call` (singleton).
- **Mutation at distinguishing position can switch call**: under
  controlled-truth + targeted mutation, the call can flip to a
  different allele; `truth_v_call` preserves the original sample.
- **Mutation at non-distinguishing position preserves the call**:
  sub-distinguishing mutations don't drop truth from the
  tie-set.
- **N at distinguishing position widens tie-set**: pinned via
  `restrict_alleles(v="v1*01") + ambiguous_base_calls(count=1)`
  + a seed where N lands at the distinguishing position.
- **Lowercase ≡ uppercase from the walker's POV**:
  sequencing-error fixture where a lowercase substitution at a
  distinguishing position flips the call (proves lowercase is
  *not* treated as wildcard).
- **Per-mechanism counter isolation**: `mutate` only bumps
  `n_mutations`; `pcr_amplify` only bumps `n_pcr_errors`;
  `sequencing_errors` only bumps `n_quality_errors`. The walker
  doesn't care about the source.
- **D allele disambiguation under VDJ**: parallel of the V tests
  for the D segment.
- **Replay round-trip**: trace replay reproduces `v_call`,
  `d_call`, `j_call`, identity fields, and counters under SHM /
  PCR / sequencing / N / mixed pipelines.

The audit doc and the test file are designed to be read
together: doc explains why, tests show what.
