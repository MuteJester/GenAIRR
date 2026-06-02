# Productive-contract failure-mode audit

**Status:** audit + golden tests. Drift items in §6.

Counterpart to
[`tests/test_productive_stress_matrix.py`](../tests/test_productive_stress_matrix.py),
which proves "valid stacks stay productive." This audit proves the
*other half*: invalid choices fail predictably, surface diagnoses
that name the offending pass/address/reason, and never silently
corrupt biology. Audit-first like the prior slices —
[primer_trim_end_loss](primer_trim_end_loss_audit.md),
[indel_provenance](indel_provenance_audit.md), and
[allele_call](allele_call_audit.md) audits.

---

## 1. Contracts in the productive bundle

The `productive()` factory in
[`contract/mod.rs`](../engine_rs/src/contract/mod.rs) bundles
four contracts:

| Contract                          | `name()`                            | Checks                                                                 |
|-----------------------------------|-------------------------------------|------------------------------------------------------------------------|
| [`ProductiveJunctionFrame`](../engine_rs/src/contract/productive_junction_frame.rs) | `"productive_junction_frame"`       | Junction length divisible by 3 (V_anchor → J_anchor + 3).               |
| [`NoStopCodonInJunction`](../engine_rs/src/contract/no_stop_codon_in_junction.rs) | `"no_stop_codon_in_junction"`       | No TAA / TAG / TGA codons inside the junction window.                  |
| [`AnchorPreserved::V`](../engine_rs/src/contract/anchor_preserved.rs)            | `"anchor_preserved.v"`              | V anchor codon (Cys) amino-acid identity preserved.                    |
| [`AnchorPreserved::J`](../engine_rs/src/contract/anchor_preserved.rs)            | `"anchor_preserved.j"`              | J anchor codon (Trp/Phe) amino-acid identity preserved.                |

A fifth contract,
[`JunctionStopState`](../engine_rs/src/contract/junction_stop_state.rs),
isn't user-facing — it precomputes per-NP-position admissible-base
masks so the NP-base sampler can do O(1) per-base filtering
instead of per-candidate `admits_typed` dispatch.

---

## 2. Execution policy + sentinel infrastructure

[`ExecutionPolicy`](../engine_rs/src/compiled/mod.rs) determines
runtime behaviour when a sampler's admissible support is empty:

| Policy        | Behaviour on empty support                                                      |
|---------------|---------------------------------------------------------------------------------|
| `Permissive`  | Apply the pass's `EmptySupport` policy: write a sentinel value or skip the slot. |
| `Strict`      | Raise `PassError::ConstraintSampling` → Python `StrictSamplingError(pass, address, reason)`. |

[`EmptySupport`](../engine_rs/src/dist/filtered.rs#L115) is the
per-pass declared fallback:

| Pass                  | Address                              | Policy                       | Sentinel value     |
|-----------------------|--------------------------------------|------------------------------|--------------------|
| GenerateNP (length)   | `np1.length` / `np2.length`          | `Sentinel(0)`                | `0` (zero-length)  |
| GenerateNP (base)     | `np1.base[i]` / `np2.base[i]`        | `Sentinel(b'N')`             | `N` (IUPAC ambiguous) |
| Trim                  | `recombine.v_trim_3` (etc.)          | `Sentinel(0)`                | `0` (no trim)      |
| IndelPass             | `corrupt.indel.kind/site/base[i]`    | NoOp (site = -1)             | `-1` sentinel      |
| MutationTransaction   | `mutate.<addr>.base`                 | `Skip`                       | (no trace, no edit) |

`PassError::ConstraintSampling` carries three fields:
`pass_name`, `address`, and a [`FilteredSampleError`](../engine_rs/src/dist/filtered.rs#L38)
reason. The Python exception
[`StrictSamplingError`](../engine_rs/src/python/contract.rs#L27)
exposes them as a 3-tuple. The reason variants serialise to:

- `"support_unavailable"` — distribution has no enumerable support.
- `"empty_admissible_support"` — every candidate rejected by the
  active contract bundle.
- `"invalid_filtered_support"` — filtered weights non-finite / ≤ 0.

---

## 3. Failure detection: compile-time vs runtime

A key structural finding: **failures bifurcate into compile-time
preconditions and runtime sample-time narrowing**, and they use
**different exception types**.

### 3.1 Compile-time precondition checks

When the *static* admissible support is empty — i.e., every
candidate in a sampling distribution violates an active contract
regardless of any runtime sampling outcome — the plan compiler
detects this in
[`Experiment.compile()`](../src/GenAIRR/experiment.py) and raises
**`ValueError`** with a message of the form:

```
plan compilation failed; contract <name> precondition failed: <reason>
```

This is **fail-fast**: the error is reported before any record
runs. It is *not* a `StrictSamplingError` — it is a Python
`ValueError` raised during plan compilation, independent of the
`strict` flag at runtime. Both `strict=False` and `strict=True`
calls raise the same error from `compile()`.

Pinned compile-time failures:

| Trigger                                                  | Reason string contains                                                  |
|----------------------------------------------------------|-------------------------------------------------------------------------|
| All NP length candidates violate frame.                  | `"NP1 length support has no in-frame mass for the active V/J anchor and trim supports"` |
| All V trim_3 candidates exceed V_anchor_to_end.          | `"V trim support removes every valid anchor"`                          |
| All J trim_5 candidates exceed J_anchor_offset.          | `"J trim support removes every valid anchor"`                          |

### 3.2 Runtime sample-time narrowing (mixed supports)

When a distribution has a *mix* of valid and invalid candidates,
the static check passes (some in-contract candidate exists), and
the runtime sampler narrows per-call to the in-contract subset.
This works identically in strict and permissive modes — neither
raises, both produce records where the chosen value is always
admissible.

Examples (pinned by `tests/test_productive_failure_modes.py`):

- `recombine(np1_lengths=[(1, 1.0), (3, 1.0)])` under
  `productive_only` — only length 3 ever lands in the trace.
- `trim(v_3=[(0, 1.0), (3, 1.0), (5, 1.0)])` — only length 0
  lands (length 3 removes the full anchor codon, length 5 exceeds
  the allele).

### 3.3 Runtime sample-time failures (strict raises, permissive sentinel)

When dynamic state — produced by an earlier pass within the same
run — makes a sampler's admissible support empty, the
compile-time check can't catch it. This is where the policy split
matters:

- **Strict**: raises `StrictSamplingError(pass_name, address,
  reason)`.
- **Permissive**: writes the declared sentinel (or skips the
  slot) and continues.

The canonical example is the indel pass under `productive_only`
on a "short-J" fixture (J = anchor codon only, no post-anchor
zone, delete-only): no admissible deletion site exists.
Pinned in [`indel_provenance_audit.md`](indel_provenance_audit.md)
§5.3 and re-pinned here in §4.

---

## 4. Failure matrix

The following matrix catalogues each detectable failure mode by
mechanism, detection time, and the exact diagnostic surface.
"Pinned" entries cite the test in
[`tests/test_productive_failure_modes.py`](../tests/test_productive_failure_modes.py).

| #  | Contract                    | Mechanism            | When             | Strict mode                                                                              | Permissive mode                  | Pinned by                                                              |
|----|-----------------------------|----------------------|------------------|------------------------------------------------------------------------------------------|----------------------------------|------------------------------------------------------------------------|
| F1 | ProductiveJunctionFrame     | NP length (all bad)  | Compile-time     | `ValueError("...NP1 length support has no in-frame mass...")`                            | Same `ValueError`                | `test_all_invalid_np_lengths_fails_at_compile_time`                    |
| F2 | AnchorPreserved::V          | V trim_3 (all bad)   | Compile-time     | `ValueError("...V trim support removes every valid anchor")`                             | Same `ValueError`                | `test_all_invalid_v_trim_3_fails_at_compile_time`                      |
| F3 | AnchorPreserved::J          | J trim_5 (all bad)   | Compile-time     | `ValueError("...J trim support removes every valid anchor")`                             | Same `ValueError`                | `test_all_invalid_j_trim_5_fails_at_compile_time`                      |
| F4 | ProductiveJunctionFrame     | NP length (mixed)    | Runtime narrow   | Narrows to in-frame lengths only. No error.                                              | Same.                            | `test_mixed_np_length_support_narrows_at_sample_time`                  |
| F5 | AnchorPreserved::V          | V trim_3 (mixed)     | Runtime narrow   | Narrows to anchor-preserving lengths. No error.                                          | Same.                            | `test_mixed_v_trim_support_narrows_at_sample_time`                     |
| F6 | NoStopCodonInJunction       | NP base mask         | Runtime narrow   | Per-position mask narrows. Always ≥ 1 admissible base (3 stop codons can ban at most 3 bases). | Same. | `test_np_base_mask_excludes_stop_completing_bases`                     |
| F7 | Productive bundle           | Indel tuple          | Runtime fail     | `StrictSamplingError("corrupt.indel", "corrupt.indel.site[0]", "empty_admissible_support")` | NoOp sentinel: `site = -1`, no event fires. | `test_indel_no_admissible_tuple_raises_strict_emits_no_op_permissive`  |
| F8 | NoStopCodonInJunction       | SHM site/base        | Runtime narrow   | Per-base narrowing; SHM at rate 0.99 still completes without raising (3 stop codons can never ban all 4 bases). | Same. | `test_shm_high_rate_under_productive_never_fails_strict`               |

### 4.1 Diagnostic quality (`StrictSamplingError` fields)

When a runtime failure raises, the Python exception's `.args`
tuple is `(pass_name, address, reason_string)` — a structured
diagnosis that names the offending pass, the choice address, and
the canonical reason. Example:

```python
("corrupt.indel", "corrupt.indel.site[0]", "empty_admissible_support")
```

Pinned by `test_strict_error_carries_pass_address_reason`.

---

## 5. Replay behaviour under failure modes

Replay is **strictly positional**: the i-th choice in the trace
file is consumed by the i-th sampling slot the pass would have
made. Replay validates:

- address consistency (slot expects `np1.length`, trace records
  `np1.length`),
- value-kind consistency (`Int` vs `Base` vs `Bool`),
- pass plan + refdata signatures match.

Replay does **not** re-validate the recorded value against the
active contract bundle. This has a notable consequence: a trace
recorded in permissive mode with sentinel values (e.g., indel
slot recorded as `site = -1` NoOp) replays in strict mode
**without** raising. The strict-mode error happens during
*sampling*, not during *replay* — once the choice is recorded,
replay consumes it verbatim.

This is by design (see
[`replay.rs`](../engine_rs/src/replay.rs#L20)) but has user-
visible implications, listed in §6.

Pinned by `test_replay_of_permissive_sentinel_trace_passes_under_strict`.

---

## 6. Drift identified

### 6.1 Two exception types for "contract violation" failures

Compile-time empty-support failures raise Python `ValueError`
from `Experiment.compile()`; runtime empty-support failures raise
`StrictSamplingError` from `run_records()`. A user catching one
won't catch the other. There's no shared base class, no shared
field set (`ValueError` has only the message string;
`StrictSamplingError` has structured `(pass, address, reason)`).

**Tradeoff:** the two error types reflect genuinely different
detection times and structural assumptions — compile-time errors
can't carry an address because no sampling slot has been
allocated yet. Unifying them would force one side to lose
information. But the inconsistency surprises users.

**Possible fix:** introduce a shared `ContractError` base, with
compile-time errors carrying `pass=None, address=None`. Or
document the two-error pattern prominently in the strict-mode
section of `productive_only()`'s docstring.

**Pinned absence by:** `test_drift_pinned_compile_and_runtime_failures_use_different_exception_types`.

### 6.2 Replay of permissive-sentinel trace under strict ignores strictness

A trace recorded under permissive mode (with sentinel values for
contract no-ops) replays under `strict=True` without raising —
the sentinel value is consumed verbatim and the same outcome
reproduces. A user who runs `compiled.replay_from_trace_file(tf,
strict=True)` expecting "rerun this with strict semantics" gets
the original outcome, not a strict-mode-fresh-run.

This is faithful reproduction of the original outcome (replay's
documented role) but elides the strict-vs-permissive distinction
*at replay time*. Users wanting strict re-execution from scratch
should use `compiled.simulator.run(seed=..., strict=True)`, not
`replay_from_trace_file(..., strict=True)`.

**Tradeoff:** re-validating contracts during replay would either
break the "positional verbatim" invariant or require a separate
validation pass. Both are expensive. The current design treats
replay as "rebuild this exact outcome" rather than "rerun the
sampler with this seed and these inputs."

**Possible fix:** rename `strict=True` on replay to something
like `validate_contracts=True` if it's ever meant to do contract
re-validation, or document that the flag is essentially a no-op
on replay (only affects the unlikely case of trace having
genuinely invalid records).

**Pinned by:** `test_replay_of_permissive_sentinel_trace_passes_under_strict`.

### 6.3 Some failure modes are structurally unreachable

`NoStopCodonInJunction` can never make all 4 NP bases inadmissible
at a position (only 3 stop codons exist, so the per-base ban list
is at most 3 elements). Same for SHM substitution. The
`EmptySupport::Sentinel(b'N')` and `EmptySupport::Skip` infra
exists for defensive completeness but is unreachable through
standard pipelines.

This isn't a bug — defensive infrastructure for an unreachable
failure is fine. But the `EmptySupport::Sentinel(b'N')` codepath
isn't exercised by end-to-end tests; if it ever becomes
reachable through a future contract addition, no existing test
would catch a regression.

**Possible fix:** add a unit-level test that synthesises an
artificial contract with empty admissible base support and
verifies the `N` sentinel emission. Out of scope here.

**Pinned absence by:** `test_drift_pinned_np_base_sentinel_path_is_structurally_unreachable`.

---

## 7. Test coverage in this slice

[`tests/test_productive_failure_modes.py`](../tests/test_productive_failure_modes.py)
pins the failure matrix from §4 with golden tests. Coverage:

- **Compile-time precondition fail-fast** (F1, F2, F3): three
  fixtures each trigger a distinct precondition `ValueError`;
  assertions check the reason substring.
- **Runtime narrowing on mixed supports** (F4, F5): empirical
  distribution checks confirm only in-contract values appear in
  the trace.
- **NP-base stop-codon masking** (F6): high-NP-length fixture
  with stop-completing V context — confirm no stop codons ever
  appear in the junction across many seeds.
- **Indel no-admissible-tuple** (F7): short-J fixture + delete-
  only count=1 + productive_only. Strict raises; permissive
  records site=-1 sentinel.
- **SHM under productive doesn't fail strict** (F8): high-rate
  mutation never empties the per-base support.
- **Diagnostic quality**: `StrictSamplingError.args` tuple
  matches `(pass_name, address, reason_string)`.
- **Replay drift** (§6.2): permissive-mode sentinel trace
  replays under strict without raising — pinned with explicit
  comparison to a fresh strict run.
- **Pin-the-drift** (§6.1, §6.3): assert absence of unified
  contract-error base and that the NP-base sentinel codepath is
  not exercised by the failure matrix.
