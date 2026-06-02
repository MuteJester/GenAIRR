# SHM Model / Context-Mutation Audit

**Status: audit only.** Pins today's somatic-hypermutation (SHM)
architecture and the gaps a future "SHM model becomes
cartridge-owned" slice will need to close. No implementation is
proposed; the deliverable is the shared vocabulary plus contract
pins so any later slice lands against an audited baseline.

This audit re-enters biology after a long architecture stretch
(clonal split, parent observability, cartridge inspectability).
The engine ships two SHM models — uniform and S5F — plus the
corruption family (PCR / quality / N / indel / end-loss /
strand / paired-end). The split between "biological SHM" and
"observation-stage corruption" is the load-bearing semantic the
audit pins.

Companion to
[`tests/test_shm_model_contract.py`](../tests/test_shm_model_contract.py)
which freezes every claim here as either `pin_scaffold_*`
(today's contract) or `pin_absence_*` (deferred surface).

The audit is deliberately narrow: it documents what exists and
the boundaries between biological mutation and engineering
artefact. The proposed next slice (S5F kernel metadata in the
cartridge manifest + per-segment biological counters) is
sketched in §11 but explicitly not committed to.

---

## Existing scaffolding the contract pins

| Surface | Where | What it pins |
|---|---|---|
| `Simulation.mutation_count: u32` | [`engine_rs/src/ir/simulation.rs:43`](../engine_rs/src/ir/simulation.rs#L43) | IR-side biological-mutation counter; survives parent→descendant boundary via the assignments clone. |
| `Simulation::with_mutation_count` | [`engine_rs/src/ir/simulation.rs:89`](../engine_rs/src/ir/simulation.rs#L89) | Persistent setter; preserves all other Simulation fields. |
| `MutationTransaction::add_to_mutation_count` | [`engine_rs/src/passes/mutation_transaction/mod.rs:169`](../engine_rs/src/passes/mutation_transaction/mod.rs#L169) | THE single point where biological SHM counts are incremented. Called only by uniform + S5F. |
| `UniformMutationPass` | [`engine_rs/src/passes/mutate/uniform.rs:46`](../engine_rs/src/passes/mutate/uniform.rs#L46) | Position-independent SHM model; bumps `mutation_count` by realised count. |
| `S5FMutationPass` | [`engine_rs/src/passes/mutate/s5f/execution.rs`](../engine_rs/src/passes/mutate/s5f/execution.rs) | Context-aware SHM model; also bumps `mutation_count` by realised count. |
| Corruption passes (`corrupt/pcr/quality/ncorrupt/indel/...`) | [`engine_rs/src/passes/corrupt/`](../engine_rs/src/passes/corrupt/) | Fire `BaseChanged` / `IndelInserted` / `IndelDeleted` events; do NOT call `add_to_mutation_count`. |
| `n_mutations` AIRR field | [`engine_rs/src/airr_record/builder.rs:211`](../engine_rs/src/airr_record/builder.rs#L211) | IR-sourced: `rec.n_mutations = sim.mutation_count as i64`. |
| `mutation_rate` AIRR field | [`engine_rs/src/airr_record/builder.rs:212-216`](../engine_rs/src/airr_record/builder.rs#L212-L216) | Derived: `n_mutations / sequence_length` (or 0 when sequence_length == 0). |
| S5F kernel loader | [`src/GenAIRR/_s5f_loader.py:82-100`](../src/GenAIRR/_s5f_loader.py#L82-L100) | Four bundled kernels: `hh_s5f`, `hh_s5f_60`, `hh_s5f_opposite`, `hkl_s5f`. Loaded at `Experiment.mutate(...)` time, not from `DataConfig`. |
| Trace addresses | [`engine_rs/src/address.rs`](../engine_rs/src/address.rs) | Per-pass count / per-site / per-base addresses: `mutate.{uniform,s5f}.{count,site[i],base[i]}`. |

---

## 1. Q1 — What mutation models exist today?

Two **biological SHM** models plus six **observation-stage
corruption** passes. The split is enforced by which passes call
`tx.add_to_mutation_count(applied)` and consequently which
contribute to `Simulation.mutation_count` → AIRR `n_mutations`.

### Biological SHM

| Model | DSL surface | Trace addresses | Behaviour |
|---|---|---|---|
| Uniform | `.mutate(model="uniform", count=N)` or `.mutate(model="uniform", rate=r)` | `mutate.uniform.count` + `mutate.uniform.site[i]` + `mutate.uniform.base[i]` | Position-independent: pick `N` positions uniformly across the assembled pool; draw a replacement base from `base_dist`. |
| S5F | `.mutate(model="s5f", count=N, s5f_model="hh_s5f"\|"hh_s5f_60"\|"hh_s5f_opposite"\|"hkl_s5f")` | `mutate.s5f.count` + `mutate.s5f.site[i]` + `mutate.s5f.base[i]` | Context-aware: 1024 5-mer contexts → mutability + per-source-base substitution rows. Weighted site selection by mutability mass. |

Both share:
- The same `MutationTransaction` scoped session.
- The same realised-count semantics (only successful substitutions bump `mutation_count`; contract-rejected sites don't double-count).
- The same constraint-aware filtering via `substitute_position_constrained` / `force_substitute_base`.

### Observation-stage corruption (NOT biological SHM)

| Pass | DSL surface | Trace addresses | What it does |
|---|---|---|---|
| PCR | `.pcr_amplify(count=N)` | `corrupt.pcr.count` + `corrupt.pcr.site[i]` + `corrupt.pcr.base[i]` | Substitute `N` bases; models polymerase errors during amplification. |
| Quality | `.sequencing_errors(count=N)` | `corrupt.quality.count` + per-site | Lowercase-mark `N` bases; models per-read quality drops. |
| N-corruption | `.ambiguous_base_calls(count=N)` | `corrupt.ns.count` + per-site | Replace `N` positions with literal `N` byte. |
| Indel | `.polymerase_indels(count=N)` | `corrupt.indel.count` + per-event | Insert / delete bases. |
| End-loss | `.end_loss_5prime/3prime(length=L)` | `corrupt.end_loss.5` / `corrupt.end_loss.3` | Truncate the pool. |
| Strand orientation | `.random_strand_orientation(prob=p)` | `corrupt.rev_comp.applied` | Trace-only flag; no IR mutation (projection-time reverse-complement). |
| Contaminant | `.contaminate(prob=p)` | `corrupt.contaminant.*` | Replace assembled pool with uniform random bases. |
| Paired-end | `.paired_end(...)` | `paired_end.r1_length` / `r2_length` / `insert_size` | Read-layout sampling. |

**None of the corruption passes call
`tx.add_to_mutation_count(...)`.** They fire `BaseChanged` /
`IndelInserted` / `IndelDeleted` events but those events go to
the per-pass `EventRecord` ledger; they do not update
`Simulation.mutation_count`. This is the architectural separation
the audit pins.

Pinned by `pin_scaffold_only_mutate_passes_call_add_to_mutation_count`
(source-level grep) and `pin_scaffold_corruption_does_not_bump_n_mutations`
(behavioural).

---

## 2. Q2 — Which mutations count as biological SHM?

The user's question — which mutations count toward `n_mutations`?

### Counts as SHM (bumps `Simulation.mutation_count`)

- Realised substitutions from `UniformMutationPass` and
  `S5FMutationPass`. Both passes call
  `tx.add_to_mutation_count(applied)` where `applied` is the
  number of substitutions that actually wrote (constraint-
  rejected sites are skipped in permissive mode and don't
  double-count).

### Does NOT count as SHM

- PCR errors (`pcr_amplify`) — modelled as polymerase artefacts.
- Quality errors (`sequencing_errors`) — read-stage degradation.
- N-corruption (`ambiguous_base_calls`) — observation-stage.
- Polymerase indels (`polymerase_indels`) — read-stage.
- End-loss (`end_loss_*prime`) — observation-stage truncation.
- Strand orientation (`random_strand_orientation`) — projection-
  time flag.
- Contaminant (`contaminate`) — non-receptor material.
- Paired-end (`paired_end`) — read layout.
- Receptor revision (`receptor_revision`) — recombination-time V
  replacement. The V slice is replaced via `replace_segment`,
  which is structurally distinct from SHM; the pass does NOT
  bump `mutation_count`.
- Initial assembly + trim — not mutations at all.

### Distinction is structural, not heuristic

The Audit pin is at the call-site level: source-grep for
`tx.add_to_mutation_count` in `engine_rs/src/passes/` returns
exactly two production hits (uniform.rs, s5f/execution.rs). A
future contributor adding a new biological mutation pass must
call `add_to_mutation_count`; a future contributor adding a new
artefact pass must NOT. The contract test source-greps for
exactly these two hits and a third hit would surface the
addition for review.

### Events vs counters

Every per-site substitution fires `BaseChanged` regardless of
which pass emitted it. The events stream on
`Outcome.events[pass_index].simulation_events` lets downstream
analysis distinguish per-pass contributions (the `pass_name`
tells you whether the event came from `mutate.uniform` /
`mutate.s5f` / `corrupt.pcr` / etc.). The IR counter only counts
the biological subset.

`MutationCountChanged` is fired by the
mutation-transaction's commit path
(`SimulationBuilder::set_mutation_count` → event observer)
when `add_to_mutation_count` was called. Corruption passes don't
fire `MutationCountChanged` either.

Pinned by `pin_scaffold_n_mutations_only_reflects_shm_passes`.

---

## 3. Q3 — Is the S5F model cartridge-owned?

**No, not today.** The audit's primary completeness gap.

### Where S5F tables live

[`src/GenAIRR/_s5f_loader.py`](../src/GenAIRR/_s5f_loader.py)
loads kernels from bundled `.pkl` files (`HH_S5F_META.pkl`,
`HH_S5F_60_META.pkl`, `HH_S5F_Opposite_META.pkl`,
`HKL_S5F_META.pkl`). The four short names map to four bundled
files; the DSL `Experiment.mutate(model="s5f", s5f_model="hh_s5f")`
resolves the short name via `_BUILTIN_S5F_MODELS` and loads
the kernel at compile time.

The loaded kernel becomes a `Box<dyn S5fKernel>` passed to
`S5FMutationPass::new`. It is **not** stored on `DataConfig`,
**not** transferred via `dataconfig_to_refdata` to
`RefDataConfig`, and **not** part of
`refdata.content_hash()`.

### What `DataConfig.asc_tables` is

`asc_tables` (allele-specific corrections) is a different
data structure — `DataConfig.asc_tables: Dict[str, Any]`. It's
already pinned as an **orphan field** in the Reference Cartridge
Completeness Audit (§1 Q1) — present on DataConfig, NOT in any
cartridge view, NOT crossed to Rust. It's not the S5F kernel
either; it's a separate allele-call-correction table used by
some legacy validator paths.

### Is S5F in `cartridge_manifest()`?

**No.** The bundled cartridge's manifest has no `s5f` key
anywhere in its JSON. `cartridge_manifest()["models"]` describes
the empirical NP-length / trim distributions (Python-side); S5F
kernel choice happens at `Experiment.mutate(...)` time, so it's
a per-experiment configuration, not a per-cartridge property.

### Is S5F hashed anywhere meaningful?

**No** — same v1 boundary as `reference_models`. Swapping
`s5f_model="hh_s5f"` for `s5f_model="hkl_s5f"` changes simulation
output materially (different mutability ratios per 5-mer
context) but does not change `compute_checksum`,
`content_hash`, or any visible manifest field.

### Verdict

This is a **documented v1 boundary**, not a bug. S5F kernels
were designed as a simulation-time parameter, not a cartridge-
attached property. A user comparing two simulation runs that
used different S5F kernels cannot tell from the cartridge
manifest alone.

The audit's §11 Slice 1 recommends a `models.s5f_kernel_name` /
`models.s5f_kernel_digest` field on the manifest (and
optionally folded into `identity.source` so `content_hash`
catches it) — but defers implementation.

Pinned by:
- `pin_scaffold_s5f_kernels_bundled_separately` — the four
  builtin names + their `.pkl` paths.
- `pin_absence_no_s5f_in_cartridge_manifest` — the gap.
- `pin_absence_no_s5f_in_content_hash` — same gap, Rust side.

---

## 4. Q4 — Is targeting restricted to V only, VDJ, junction, or
##        full assembled sequence?

**Full assembled sequence — both uniform and S5F mutate the
entire pool.**

### Mechanism

Both passes start with:

```rust
let pool_len = sim.pool.len() as u32;
let count = sample_validated_count(&self.count_source, ctx, pool_len, ...);
```

Then the site selection loop ranges over `[0, pool_len)` —
covering V, NP1, D, NP2, J indiscriminately. The S5F path
weights sites by mutability ratio but doesn't restrict the
candidate range.

### Behavioural confirmation

A 30-mutation SHM run on a ~379-bp IGH heavy-chain pool produces
~29 differences distributed across the full sequence (audit
probe; pinned). If targeting were V-only, the differences would
concentrate in the V region; if junction-only, in CDR3.

### Biological note

This is a known modeling simplification. Real B-cell SHM
concentrates on V (with rate spikes at CDR1/CDR2) and some on J,
with much less activity at NP regions and CDR3. The current
model treats the pool as a flat substrate.

A future "per-segment SHM rate" slice (audit §11 Slice 2) would
add per-region rate scalars; out of scope for this audit.

Pinned by `pin_scaffold_shm_targets_full_assembled_pool` (probe-
level: differences span the full sequence at high SHM count).

---

## 5. Q5 — How does productive-only interact with SHM?

`productive_only()` attaches a `ContractSet` bundle to the
compiled artifact. The bundle's contracts gate per-site
admissibility for every substitution proposed by the mutation
transaction:

- `NoStopCodonInJunction` — rejects bases that would introduce
  a stop codon inside the junction codon rail.
- `AnchorPreserved` — rejects bases that would mutate the V/J
  anchor codon.
- `ProductiveJunctionFrame` — pre-fork compile-time gate on
  the NP-length × trim × allele combinations.

### Constrain-before-propose flow

`substitute_position_constrained` filters the candidate base
distribution against the active `ContractSet` BEFORE drawing.
When the admissible support shrinks to empty:

- **Permissive mode**: the slot is silently skipped (the
  realised count doesn't include it; `n_mutations` reflects only
  realised substitutions).
- **Strict mode**: a `ConstraintSampling` error is raised.

### Replay safety

Rejected sites do not appear in the trace at the rejected slot;
the trace records the site + base that was actually written. So
replay re-fires the same write at the same slot and the same
base draws. Trace replay against `strict=True` continues to work
even when the original recording was made with `strict=False`
and silent-skipped some slots — replay consumes the recorded
values verbatim; it doesn't re-evaluate admissibility (this is
the documented `replay against permissive sentinels` policy,
[`engine_rs/src/python/compiled.rs:230-249`](../engine_rs/src/python/compiled.rs#L230-L249)).

### Behavioural verification

`recombine().mutate(count=20).productive_only()` on
`human_igh` produces records with `productive=True` for every
descendant even at `count=20` — the contract filtering
preserved the productive triad across all 20 substitution
attempts (some may have skipped silently in permissive mode).

Pinned by `pin_scaffold_productive_only_preserves_triad_under_shm`.

---

## 6. Q6 — Allele-call ambiguity under SHM

The audit's question: does SHM affect evidence-derived calls,
and do truth calls stay stable?

### Evidence-derived calls drift; truth calls don't

- `truth_v_call` / `truth_d_call` / `truth_j_call` are derived
  from `Simulation.assignments.{v,d,j}.allele_id` — the
  recombination-time identities. SHM doesn't write to
  assignments; only `receptor_revision` does (and it writes
  via `replace_segment` with provenance in
  `receptor_revision_original_id`, per Bug D fix). So truth
  calls are IR-stable under any SHM count.

- `v_call` / `d_call` / `j_call` (the live-call columns) are
  evidence-derived: the live-call walker scores allele
  candidates per pool byte and surfaces the tied-best set. Under
  heavy SHM, the original allele's distinguishing positions
  get rewritten, so the tied-best set widens to a comma-separated
  list. The truth allele continues to appear in the tied set
  (it scores from the un-mutated positions), but the live-call
  string can grow.

### Behavioural confirmation

`recombine().mutate(count=15)` with `expose_provenance=True` on
`human_igh`:

```
seed 0: truth_v=IGHVF6-G21*15, v_call=IGHVF6-G21*15
seed 1: truth_v=IGHVF8-G30*01, v_call=IGHVF8-G30*01,IGHVF8-G30*02
seed 2: truth_v=IGHVF9-G32*01, v_call=IGHVF9-G32*01
```

Sometimes the live call widens; sometimes it doesn't. The truth
call is always present and always a single allele name. This is
the design — already pinned by
`test_g8_truth_retained_in_live_call_under_heavy_shm` in
`tests/test_experiment.py`.

### Identity-to-projected vs identity-to-truth gaps

The audit names two distinct gap shapes a future analysis tool
might want to surface:

- **Identity-to-projected gap** — `truth_v_call != v_call`
  literal. Under heavy SHM with tie-set widening, this is
  expected and not an error. A future "SHM ambiguity report"
  could quantify it per cartridge / per rate.
- **Identity-to-truth gap** — `truth_v_call not in
  v_call.split(",")`. This would be a real bug — the live caller
  dropped the truth allele. The score-and-tie caller's design
  prevents this from happening; the family-aware validator
  (`validate_families_with_parents`) catches it as a
  `ParentTruthVCallMismatch` if it ever did.

Both gaps are documented but no per-segment analysis surface
exists today.

Pinned by `pin_scaffold_truth_calls_stable_under_shm` and
`pin_scaffold_live_call_can_widen_under_heavy_shm`.

---

## 7. Q7 — Validator coverage

### What the validator catches

`Outcome.validate_record(refdata, sequence_id=...)` re-derives
an AIRR record from the outcome's IR + trace and compares it
field-by-field against the supplied record. It catches:

- `n_mutations` divergence between the expected (re-derived from
  `sim.mutation_count`) and the supplied record's `n_mutations`.
- `mutation_rate` divergence (since it's derived from
  `n_mutations`).
- Field-level tampering of any AIRR projection field.

### What the validator does NOT catch

- **Dict-level mutation count tampering through
  `validate_records`**. The validator re-projects the record
  internally; the supplied record dict isn't compared to the
  outcome state directly. So writing `record["n_mutations"] =
  999` and calling `result.validate_records(refdata)` returns
  ok — the validator re-derives, gets the right count, and
  ignores the dict.

  This is the same documented behaviour as the receptor-
  revision and `d_inverted` cases. Not a bug — Python's
  `validate_records` is a "re-derive and compare against
  outcome state" surface, not an "audit record dict against
  outcome" surface. To audit dict-level tampers you'd call
  `validate_airr_record` on the Rust side directly with a
  custom record (the Rust unit tests cover this path).

- **Per-site mutation position trace consistency**. The
  validator doesn't verify "every `BaseChanged` event in the
  trace corresponds to a difference between the pre- and post-
  mutation pool bytes." A pass that fired the right
  `add_to_mutation_count` value but didn't actually write the
  bases would pass `n_mutations` validation. The
  `MutationTransaction` boundary lockdown test in
  `engine_rs/src/passes/mutation_transaction/mod.rs:boundary_lockdown`
  protects this at the structural level (no pass mutates via
  low-level builder calls outside the transaction), but no
  runtime validator re-derives the bases.

### Mutation-ledger validator — future surface

A future `Outcome.validate_mutation_ledger(refdata)` could:

1. For each `BaseChanged` event from a `mutate.*` pass, verify
   the post-state pool byte at `event.position` differs from
   the pre-state byte.
2. Count realised substitutions and assert
   `count == sim.mutation_count - pre_mutation_count`.
3. Verify per-pass realised counts match
   `trace.find(mutate.{model}.count)`.

Out of scope for this audit; pinned absent.

Pinned by `pin_scaffold_validate_records_redrives_n_mutations`
and `pin_absence_no_mutation_ledger_validator`.

---

## 8. Q8 — Performance

Not deeply probed in this audit — performance audits are a
separate slice. But the surface-level picture:

### S5F per-site cost

S5F enumerates 1024 5-mer contexts and computes a per-position
mutability profile across the pool. For a ~380 bp IGH heavy-
chain pool:

- Profile construction: O(pool_len) — one 5-mer lookup per
  position.
- Site sampling: O(pool_len) — one cumulative-weight pass per
  mutation.
- Total per mutation: O(pool_len).
- Total per record: O(count × pool_len).

At `count=20`, this is ~7600 operations per record; in the
release-tier benchmark suite this lands well inside the
budget. The current performance budget tests in
[`tests/test_performance_budgets.py`](../tests/test_performance_budgets.py)
cover SHM under `count=10` and `rate=0.03` cases.

### Uniform per-site cost

Uniform doesn't build a mutability profile — site selection is
O(1) via uniform draw. Per-record cost is O(count). Cheaper
than S5F by a constant factor.

### Constraint filtering overhead

`substitute_position_constrained` adds per-site contract
evaluation. Under `productive_only()`, every base draw incurs
the `NoStopCodonInJunction` codon-rail check. This is the
dominant cost at high mutation rates with productive
constraints active; the existing benchmark covers it.

### No new performance pins added by this audit

The audit doesn't add performance pins beyond what
`tests/test_performance_baseline.py` already covers. A future
"SHM performance audit" slice would carve out per-pass cost
attributions and per-region cost (since SHM hits the full pool).

Pinned by `pin_scaffold_performance_baseline_covers_shm` — the
release-tier baseline test for SHM exists and is invoked.

---

## 9. Replay determinism

SHM is fully deterministic under same-seed replay:

- Same `seed` → same `count` draw.
- Same RNG state → same site sequence.
- Same site + same context → same kernel row (S5F) or same
  uniform draw.
- Constraint filtering produces the same admissible mask given
  the same pre-state pool, so same draws even under
  `productive_only()`.

Behavioural pin: two runs of the same SHM-bearing experiment
with the same seed produce byte-identical sequences across
all records. Already covered by Slice-2 byte-reproducibility
tests; re-pinned here for SHM-specific cross-doc traceability.

### Replay against a different S5F kernel

If a trace recorded with `hh_s5f` is replayed via
`replay_from_trace_file` against a plan compiled with
`hkl_s5f`, the replay STILL produces the same sequences —
because replay consumes the recorded site + base values
verbatim. The kernel only affects the sampling distribution,
not the actual writes. This is documented as a known
"replay-consumes-recorded-values" property of the engine
([`compiled.rs:230-249`](../engine_rs/src/python/compiled.rs#L230-L249));
the kernel mismatch is silent because the trace already carries
the chosen sites + bases.

This is the same shape as the `reference_models` and `S5F
kernel in content_hash` boundaries — a v1 acceptance, not a
bug. Pinned absent: no kernel digest in `content_hash`.

Pinned by `pin_scaffold_shm_replay_byte_deterministic`.

---

## 10. Validator integration

Slice 1's `validate_families` and Slice 3's
`validate_families_with_parents` both stay applicable under
SHM. The parent-aware validator already pins
`truth_v_call` (IR-sourced) invariance across descendants of a
clone; descendants diverge through SHM but their truth alleles
stay locked to the parent's.

No SHM-specific validator changes proposed by this audit.

---

## 11. Implementation order (recommended, NOT committed)

After the audit, three concrete slices in order of leverage:

### Slice 1 — SHM model metadata in cartridge manifest

Python-only:
- Plumb the S5F kernel name (or "uniform" / "none") through
  the `Experiment` builder so a manifest call can surface it.
- Add `manifest["models"]["shm_model"] = "uniform" | "s5f"`
  and `manifest["models"]["s5f_kernel_name"] = str | None`.
- Optionally add a `manifest["models"]["s5f_kernel_digest"]`
  by hashing the loaded kernel's mutability + substitution
  arrays.

Does NOT touch `content_hash` — keeps the v1 boundary
explicit until a separate slice closes it. Closes the
"two simulation runs with different S5F kernels look identical
in the manifest" gap.

### Slice 2 — Per-segment SHM rate scalars (biology slice)

Add `mutate(per_segment_rates={"V": 1.0, "D": 0.5, "J": 0.5, "NP1": 0.1, "NP2": 0.1})`
or similar. The mutation pass would consult per-position
region membership and weight the site distribution by the
segment's rate. Closes the "uniform/S5F treat the pool as
flat" modeling simplification.

### Slice 3 — Mutation-ledger validator

Add `Outcome.validate_mutation_ledger(refdata)` per §7.
Closes the "n_mutations is correct but per-site bases aren't
verified" runtime-validation gap.

This audit only proposes the architecture and pins absences;
each slice flips the relevant `pin_absence_*` to
`pin_present_*` in lockstep.

---

## 12. Test surface — what this audit pins

Mirrored in
[`tests/test_shm_model_contract.py`](../tests/test_shm_model_contract.py).

### `pin_scaffold_*` — today's contract

1. **Two SHM models exist**: `UniformMutationPass` /
   `S5FMutationPass`. Identified by source-level grep.
2. **`add_to_mutation_count` called only from the two SHM
   passes** (source-level pin).
3. **Corruption passes don't bump `n_mutations`** (behavioural):
   PCR / quality / N-corruption / indel runs with zero SHM
   produce `n_mutations == 0`.
4. **SHM does bump `n_mutations`** (behavioural): a `count=N`
   mutate run produces `n_mutations == N` (realised count).
5. **Combined SHM + corruption** — `n_mutations` reflects only
   the SHM count.
6. **`mutation_rate = n_mutations / sequence_length`** —
   derived, not a separate sample.
7. **SHM targets the full assembled pool** — substituted
   positions span V, NP1, D, NP2, J at high mutation count.
8. **Productive-only preserves triad under SHM** — every
   descendant of `productive_only() + mutate(count=20)` is
   productive.
9. **Truth calls stable under SHM** —
   `truth_v_call != ""` and matches the recombination-time
   allele.
10. **Live call can widen under heavy SHM but always includes
    the truth allele**.
11. **`validate_records` re-derives `n_mutations`** — clean
    batch validates; dict-level tampers are silently ignored
    (documented behavior, not a bug).
12. **Replay is byte-deterministic for SHM** — same seed →
    same sequences.
13. **S5F kernels bundled separately** — four canonical names
    (`hh_s5f`, `hh_s5f_60`, `hh_s5f_opposite`, `hkl_s5f`),
    loaded from `.pkl` files in `_s5f_loader.py`.
14. **Performance baseline covers SHM** — the release-tier
    perf-baseline file exists and exercises SHM.

### `pin_absence_*` — gaps Slices 1-3 close

15. **No S5F kernel name in cartridge manifest** — `manifest`
    JSON dump contains no "s5f" string.
16. **No S5F digest in `content_hash`** — refdata's
    content_hash is identical whether the experiment uses
    `hh_s5f` or `hkl_s5f` (sampling choices live in the trace,
    not in identity).
17. **No `shm_model` AIRR field** — record dict has
    `n_mutations` and `mutation_rate` but no string naming the
    model.
18. **No per-segment mutation rate / count AIRR fields** —
    record has `n_mutations` (total) but no
    `v_n_mutations` / `j_n_mutations` / etc.
19. **No `validate_mutation_ledger` method** on `Outcome` or
    `SimulationResult`.
20. **No mutation-position provenance AIRR field** — record
    doesn't carry the list of mutated positions; consumers
    must read `outcome.events()` and filter for
    `mutate.uniform.site[i]` / `mutate.s5f.site[i]` records.

### Doc anchor

21. The audit doc continues to exist and references the
    contract file; the 14-section structure stays intact.

---

## 13. Out of scope

- **Targeting biology (per-segment rates)**. Slice 2 territory;
  documented as a known modeling simplification but not changed.
- **AID/UNG biology** (somatic hypermutation by activation-
  induced cytidine deaminase mechanism). The S5F kernel
  approximates the empirical outcome; the audit doesn't model
  AID directly.
- **Affinity maturation / selection coefficients.** Out of
  scope — would compose against a future `FamilyOutcome`
  aggregate, not against SHM itself.
- **Class switch recombination.** Separate biology; out of
  scope.
- **Indel SHM (vs polymerase indels).** The current model only
  produces SHM substitutions; biological insertion / deletion
  events during SHM are not modelled. Out of scope.
- **Per-base context window > 5-mer.** S5F's 5-mer context is
  the documented kernel surface; longer contexts (7-mer, etc.)
  would need a separate kernel format.

---

## 14. Summary table

| Concern | Today's contract |
|---|---|
| SHM models | Two: `UniformMutationPass`, `S5FMutationPass`. |
| Counts as SHM | Only the two SHM passes bump `Simulation.mutation_count`. |
| AIRR `n_mutations` source | IR-sourced (`sim.mutation_count`); survives parent→descendant boundary. |
| AIRR `mutation_rate` | Derived: `n_mutations / sequence_length`. |
| Corruption passes | Fire `BaseChanged` events but do NOT bump `n_mutations`. |
| SHM target region | Full assembled pool (V + NP1 + D + NP2 + J). |
| S5F kernel location | Bundled `.pkl` files in `_s5f_loader.py`; NOT on DataConfig. |
| S5F in `cartridge_manifest` | Absent (v1 boundary). |
| S5F in `content_hash` | Absent (v1 boundary). |
| Productive-only interaction | Constrain-before-propose; realised-count semantics for SHM. |
| Truth calls under SHM | Stable (IR-sourced). |
| Live calls under SHM | Can widen; truth allele always in the tie set. |
| Validator coverage | `validate_records` re-derives `n_mutations`; dict-level tampers ignored. |
| Replay determinism | Byte-deterministic; kernel swap on replay silent (trace consumes recorded sites/bases). |
| Performance | O(count × pool_len) for both models; existing baseline covers. |
| Required for "complete" SHM | (1) Kernel metadata in manifest, (2) per-segment rates, (3) mutation-ledger validator. |

The SHM model today is **architecturally clean** — biology
counts separately from artefacts via the
`add_to_mutation_count` single-source-of-truth, the IR-sourced
counter survives the clonal boundary, and replay is exactly
deterministic. What's **not yet complete** is the cartridge
ownership of the kernel choice, the per-segment biological
targeting, and the position-level mutation-ledger validation.

The recommended next slice (S5F kernel metadata in the manifest)
closes the most user-visible gap with the smallest Python-only
diff. The deeper biology slices warrant their own audits.
