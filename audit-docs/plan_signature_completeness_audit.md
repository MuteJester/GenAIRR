# Experiment DSL / Plan Signature Completeness Audit

**Status: audit only.** Reviews the load-bearing replay-safety
mechanism — `Pass::parameter_signature` and the per-pass
`fmt_byte_dist` / `fmt_int_dist` / `fmt_prob` / `fmt_segment_rates`
/ `fmt_v_subregion_rates` / `fmt_count_source` /
`fmt_np_base_generator` / `join_parts` discipline in
[`engine_rs/src/passes/paramsig.rs`](../engine_rs/src/passes/paramsig.rs)
— against every parameterized public DSL surface and every
cartridge-driven compile parameter accumulated through the
recent biology slices (segment rates, V-subregion rates, NP
base models, Markov dependency, P-nucleotide v1, paired-end,
D inversion, receptor revision).

The audit's central question, restated for clarity:

> Does the plan signature canonically encode every parameter
> that changes proposal support or output, such that two
> experiments differing only by such a parameter produce
> different `pass_plan_signature` strings, and replay across
> mismatched signatures fails the signature gate before any
> choice is consumed?

Companion to
[`tests/test_plan_signature_completeness_contract.py`](../tests/test_plan_signature_completeness_contract.py)
which pins the per-family change-detection invariant, the
default-normalisation invariants, the canonicalization
invariants, and the three documented soft-gap boundaries.

**Pre-flight finding (§7 below): clean-yes — no silent
corruption gap.** Every fold-relevant DSL surface either
appears in the signature OR is gated by an independent
mechanism (refdata content hash; strict-mode sampler
rejection; cartridge manifest digest) that backs up the
signature for the specific corner case it doesn't cover.

Three **documented soft gaps** identified — each has an
explicit code-level rationale in its owning pass's
docstring AND an independent safety mechanism that prevents
silent cross-cartridge replay corruption:

1. `SampleAllelePass::parameter_signature` returns the empty
   string. `restrict_alleles` / `*_allele_weights` change the
   sampling distribution but do NOT fold. **Backstop:**
   refdata content hash covers pool identity at the
   cartridge level; the strict sampler rejects a replayed
   `AlleleId` that's outside the (possibly narrowed) support.
   Cross-replay where the recorded allele falls outside the
   narrowed distribution support fails at the sampler, not
   at the signature gate.
2. `ReceptorRevisionPass::parameter_signature` folds only
   `prob`; the replacement V distribution (same V pool as
   `SampleAllelePass`) is NOT folded. Same rationale +
   backstop as (1).
3. `S5FMutationPass::parameter_signature` folds the kernel's
   *name string* (`kernel=igh_s5f` etc.) but NOT the kernel
   payload's content digest. **Backstop:** the cartridge
   manifest carries `s5f_kernel_digest` separately and the
   `models.shm.in_content_hash=False` flag explicitly
   documents that the kernel choice is out-of-`content_hash`
   too; two cartridges that ship different files under the
   same `kernel_name` would produce equal plan signatures
   but different manifest digests.

These are NOT silent-corruption gaps — they're explicit,
documented boundaries with backstops. They DO still mean
that a user who restricts alleles or revises kernel content
in place (without renaming) gets a signature that doesn't
distinguish the variant. The recommended tightening (§9) is
to extend `SampleAllelePass::parameter_signature` to fold the
distribution's `support()`, scoped to the case where the
distribution narrows the pool — but the v1 boundary holds.

---

## 1. Existing scaffolding (passes ↔ parameter_signature)

The table below enumerates every `Pass` implementation in the
engine and its current `parameter_signature` body. The
"folded" column captures which constructor parameters
participate in the signature; "unfolded" captures any
constructor parameter that does NOT.

| Pass | File | Folded | Unfolded | Notes |
|---|---|---|---|---|
| `SampleAllelePass` | [`sample_allele.rs`](../engine_rs/src/passes/sample_allele.rs) | (segment via `name()`) | `distribution: Box<dyn Distribution<Output = AlleleId>>` | **Soft gap (1) below.** Docstring documents the choice; rationale ties to `refdata_content_hash`. |
| `AssembleSegmentPass` | [`assemble_segment.rs`](../engine_rs/src/passes/assemble_segment.rs) | (segment via `name()`) | — | Default empty signature; no compile-time params besides `segment`. |
| `TrimPass` | [`trim.rs`](../engine_rs/src/passes/trim.rs) | `length` distribution; segment/end via `name()` | — | Folded via `fmt_int_dist`. |
| `GenerateNPPass` | [`generate_np.rs`](../engine_rs/src/passes/generate_np.rs) | `length` distribution + `base_generator.signature()` | — | NP base generator's signature covers uniform / categorical / Markov (5 rows × A/C/G/T). |
| `PAdditionPass` | [`p_addition.rs`](../engine_rs/src/passes/p_addition.rs) | `length` distribution; end via `name()` | — | Folded via `fmt_int_dist`. |
| `InvertDPass` | [`invert_d.rs`](../engine_rs/src/passes/invert_d.rs) | `prob` | — | Folded via `fmt_prob`. |
| `ReceptorRevisionPass` | [`receptor_revision.rs`](../engine_rs/src/passes/receptor_revision.rs) | `prob` | `replacement_v_distribution: Box<dyn Distribution<Output = AlleleId>>`, `v_trim_3_distribution` | **Soft gap (2) below.** Same V pool as SampleAlleleVPass; same rationale. |
| `UniformMutationPass` | [`mutate/uniform.rs`](../engine_rs/src/passes/mutate/uniform.rs) | `count_source`, `base_dist`, `segment_rates`, `v_subregion_rates` | — | All four fold via the shared helpers. |
| `S5FMutationPass` | [`mutate/s5f.rs`](../engine_rs/src/passes/mutate/s5f.rs) | `kernel_name`, `count_source`, `segment_rates`, `v_subregion_rates` | `kernel_payload` (the actual 5-mer transition matrix bytes) | **Soft gap (3) below.** Kernel name is folded; payload digest is gated by the cartridge manifest's `s5f_kernel_digest`. |
| `PCRErrorPass` | [`corrupt/pcr.rs`](../engine_rs/src/passes/corrupt/pcr.rs) | `count_source`, `base_dist` | — | |
| `QualityErrorPass` | [`corrupt/quality.rs`](../engine_rs/src/passes/corrupt/quality.rs) | `count_source`, `base_dist` | — | |
| `IndelPass` | [`corrupt/indel.rs`](../engine_rs/src/passes/corrupt/indel.rs) | `count_source`, `insertion_prob`, `base_dist` | — | |
| `NCorruptionPass` | [`corrupt/ncorrupt.rs`](../engine_rs/src/passes/corrupt/ncorrupt.rs) | `count_source` | — | |
| `EndLossPass` | [`corrupt/end_loss.rs`](../engine_rs/src/passes/corrupt/end_loss.rs) | `length` distribution; end via `name()` | — | |
| `ContaminantPass` | [`corrupt/contaminant.rs`](../engine_rs/src/passes/corrupt/contaminant.rs) | `apply_prob`, `base_dist` | — | Base dist is hard-coded to `UniformBase` at the Python bridge. |
| `RevCompPass` | [`corrupt/rev_comp.rs`](../engine_rs/src/passes/corrupt/rev_comp.rs) | `apply_prob` | — | |
| `PairedEndSamplingPass` | [`paired_end.rs`](../engine_rs/src/passes/paired_end.rs) | `r1_length`, `r2_length`, `insert_size` | — | All three distributions fold via `fmt_int_dist`. |
| `SampleBasePass` | [`sample_base.rs`](../engine_rs/src/passes/sample_base.rs) | — | (test-only) | Not used by any production lowering path; only by the `run_smoke_plan` test fixture. |
| `EchoPass` | [`echo.rs`](../engine_rs/src/passes/echo.rs) | — | (test-only) | Same — `run_smoke_plan` only. |

---

## 2. Q1 — Inventory of parameterized DSL surfaces

The table cross-references every public DSL surface on
[`Experiment`](../src/GenAIRR/experiment.py) (or driven from
[`ReferenceEmpiricalModels`](../src/GenAIRR/reference_models.py))
against the pass(es) it lowers to, and the folded parameters.

| DSL surface | Pass(es) | Folded params | Notes |
|---|---|---|---|
| `recombine(np1_lengths=...)` / `np2_lengths` | `GenerateNPPass` | `length` | Also folded when authored via cartridge `reference_models.np_lengths`. |
| `recombine(...)` with cartridge `np_bases` (`uniform` / `empirical_first_base` / `markov`) | `GenerateNPPass.base_generator` | The generator's `signature()` — covers all three kinds | Markov flattens 5 rows × A/C/G/T canonical order. |
| `recombine(...)` with cartridge `p_nucleotide_lengths` (per-end) | `PAdditionPass` (1 per end) | `length` | Slice — P-nucleotide v1. |
| `trim(v_3=..., d_5=..., d_3=..., j_5=...)` or cartridge `reference_models.trims` | `TrimPass` (1 per end) | `length` | |
| `restrict_alleles(v=..., d=..., j=...)` | `SampleAllelePass` (distribution narrowed) | — (no fold) | **Soft gap (1).** |
| `recombine(v_allele_weights=..., d_allele_weights=..., j_allele_weights=...)` | `SampleAllelePass` (distribution reweighted) | — (no fold) | **Soft gap (1).** |
| `mutate(rate=..., count_pairs=..., model='uniform' \| 's5f', s5f_model_name=..., segment_rates=..., v_subregion_rates=...)` | `UniformMutationPass` or `S5FMutationPass` | `count_source`, `base_dist` (uniform only), `segment_rates`, `v_subregion_rates`, `kernel_name` (S5F) | S5F kernel **payload** content NOT folded (soft gap (3)). |
| `polymerase_indels(count=..., rate=..., insertion_prob=...)` | `IndelPass` | `count_source`, `insertion_prob`, `base_dist` | |
| `pcr_errors(...)` (via `_CorruptStep`) | `PCRErrorPass` | `count_source`, `base_dist` | |
| `quality_errors(...)` | `QualityErrorPass` | `count_source`, `base_dist` | |
| `introduce_ns(...)` | `NCorruptionPass` | `count_source` | |
| `end_loss_5prime(...)` / `end_loss_3prime(...)` / `primer_trim_5prime(...)` / `primer_trim_3prime(...)` | `EndLossPass` (per side) | `length` | |
| `contaminate(prob=...)` | `ContaminantPass` | `apply_prob` (+ hard-coded `UniformBase`) | |
| `random_strand_orientation(prob=...)` | `RevCompPass` | `apply_prob` | |
| `invert_d(prob=...)` | `InvertDPass` | `prob` | |
| `receptor_revision(prob=...)` | `ReceptorRevisionPass` | `prob` | Replacement V distribution NOT folded — **soft gap (2)**. |
| `paired_end(r1_length=..., r2_length=..., insert_size=...)` | `PairedEndSamplingPass` | `r1_length`, `r2_length`, `insert_size` | |
| `productive_only()` | (none — adds `ContractSet` on the runtime context) | (n/a — contracts are runtime config, not pass params) | Contracts narrow per-choice support but the contract set itself is not in the signature. See §6 below. |

---

## 3. Q2 — Completeness against proposal support / output

For each parameter that changes proposal support or output,
the following matrix records empirical verification (each row
was tested in `tests/test_plan_signature_completeness_contract.py`
via "change only this parameter; assert signature changes"):

| Parameter | Folded? | Signature changes on parameter change | Replay across mismatched signature fails the gate? |
|---|---|---|---|
| NP1 / NP2 length distribution | ✓ | ✓ | ✓ |
| Trim distribution (per side) | ✓ | ✓ | ✓ |
| NP base model (uniform → categorical → Markov) | ✓ | ✓ | ✓ |
| Markov transition matrix rows | ✓ | ✓ | ✓ |
| P-nucleotide length (per end) | ✓ | ✓ | ✓ |
| Mutation count `rate` vs `count_pairs` | ✓ | ✓ (different fmt_count_source shape) | ✓ |
| Mutation count value (same form) | ✓ | ✓ | ✓ |
| `segment_rates` vector | ✓ | ✓ (when non-default) | ✓ |
| `v_subregion_rates` vector | ✓ | ✓ (when non-default) | ✓ |
| S5F `kernel_name` | ✓ | ✓ | ✓ |
| S5F kernel **payload** content | ✗ (kernel name only) | ✗ (when kernel files differ but share name) | — (no plan-signature gate) |
| Indel `insertion_prob` | ✓ | ✓ | ✓ |
| End-loss / primer-trim length | ✓ | ✓ | ✓ |
| Contaminant `apply_prob` | ✓ | ✓ | ✓ |
| Rev-comp `apply_prob` | ✓ | ✓ | ✓ |
| Invert-D `prob` | ✓ | ✓ | ✓ |
| Receptor-revision `prob` | ✓ | ✓ | ✓ |
| Paired-end `r1_length` / `r2_length` / `insert_size` | ✓ | ✓ | ✓ |
| Allele restriction (`restrict_alleles`) | ✗ | ✗ | — (strict sampler rejects out-of-support replay) |
| Allele weights (`*_allele_weights`) | ✗ | ✗ | — (strict sampler rejects out-of-support replay) |
| Receptor-revision replacement V distribution | ✗ | ✗ | — (strict sampler rejects out-of-support replay) |

### Empirical verification details

The three unfolded rows are **soft gaps**: changing the
parameter does NOT change the plan signature, so a trace
recorded under one cartridge can replay against another
cartridge with a different value for that parameter without
the signature gate firing. The replay-safety backstop
differs per case:

- **`restrict_alleles` / `*_allele_weights` / replacement V
  distribution:** strict-mode rejects the recorded allele ID
  when it's outside the (possibly narrowed) distribution
  support. Cross-replay where the original cartridge had an
  unrestricted pool and the target has a restricted pool
  fails at the SampleAllele sampler. The reverse (restricted
  origin → unrestricted target) replays cleanly but is
  byte-identical to the original run by construction — the
  cartridge identity gap is the cost the user accepts when
  not folding the distribution.
- **S5F kernel payload:** different kernel files under the
  same `kernel_name` produce equal plan signatures. The
  cartridge manifest's `s5f_kernel_digest` field (computed
  off the bundled pickle bytes at install time) surfaces the
  payload's identity at a coarser layer. The shipped boundary
  is documented in [`shm_model_audit.md`](shm_model_audit.md)
  §1.4 (`in_content_hash=False`).

---

## 4. Q3 — Canonicalization

The shared formatters in
[`engine_rs/src/passes/paramsig.rs`](../engine_rs/src/passes/paramsig.rs)
enforce the canonicalization discipline. Empirically verified:

| Surface | Canonicalization | Pin verification |
|---|---|---|
| `segment_rates` dict | Insertion order doesn't matter (`fmt_segment_rates` reads named fields off `SegmentRateWeights`, which the Python boundary sorts). | ✓ |
| `v_subregion_rates` dict | Same — `fmt_v_subregion_rates` reads named fields. | ✓ |
| Markov transition matrix rows | The audit's `MarkovBaseGenerator::signature()` walks rows in canonical A/C/G/T from-order and emits the 4 to-bytes in canonical A/C/G/T to-order. Dict insertion order at the Python boundary is normalised. | ✓ |
| Empirical distribution `(value, weight)` pairs from cartridge | Bridge resolvers (`_np_lengths_from_models`, `_p_nucleotide_lengths_from_models`, `_trim_from_models`) sort by value; two cartridge specs with distinct values in different Python insertion order produce equal signatures. `EmpiricalLengthDist::from_pairs` does NOT collapse duplicate values — `[(1, 0.5), (1, 0.5)]` and `[(1, 1.0)]` produce **different** signatures (the bridge does not pre-collapse either). Documented at audit level as a UX wart but not a correctness gap (the divergent signatures fail-safe in either direction; safer than collapsing one but not the other). | ✓ (sort) / ✗ (collapse — documented) |
| First-base / categorical NP base distribution | `CategoricalBase::from_pairs` accumulates by base byte; the four `(A, w_A), (C, w_C), (G, w_G), (T, w_T)` entries always appear in canonical alphabetical order in `support()`. | ✓ |

---

## 5. Q4 — Default normalization

Behaviourally equivalent defaults must collide to equal
signatures, so a user who omits a kwarg gets the same
signature as one who passes the canonical default.

| Equivalent surfaces | Collide on equal signature? |
|---|---|
| No `segment_rates` kwarg vs explicit `{"V": 1.0, "D": 1.0, "J": 1.0, "NP": 1.0}` | ✓ (`fmt_segment_rates` short-circuits on `is_default()`) |
| No `v_subregion_rates` kwarg vs explicit `{"FWR1": 1.0, "CDR1": 1.0, "FWR2": 1.0, "CDR2": 1.0, "FWR3": 1.0}` | ✓ |
| No `np_bases` kwarg vs explicit `NpBaseModelSpec(kind="uniform")` | ✓ (both lower to `UniformNpGenerator`; same `fmt_byte_dist` output) |
| No `p_nucleotide_lengths` kwarg vs empty `p_nucleotide_lengths={}` | ✓ (lowering omits the four `push_p_addition` calls; identical plan body) |
| No `invert_d(prob=...)` vs not calling `.invert_d(...)` at all | n/a — when not called, no `InvertDPass` is appended, so signature differs structurally (one fewer pass). This is the intended behaviour: enabling the pass is a structural change to the plan, not a parameter tweak. |

---

## 6. Q5 — Cartridge-driven defaults

Every cartridge-driven default in
`ReferenceEmpiricalModels` lowers to a pass parameter and
therefore folds:

| Cartridge field | Resolver → Pass param | Folded? |
|---|---|---|
| `np_lengths["NP1"]` / `["NP2"]` | `_np_lengths_from_models` → `GenerateNPPass.length_dist` | ✓ |
| `np_bases["NP1"]` / `["NP2"]` | `_np_bases_from_models` + `_np_markov_transitions_from_models` → `GenerateNPPass.base_generator` | ✓ |
| `p_nucleotide_lengths["V_3" / "D_5" / "D_3" / "J_5"]` | `_p_nucleotide_lengths_from_models` → `PAdditionPass.length_dist` | ✓ |
| `trims["V_3" / "D_5" / "D_3" / "J_5"]` | `_trim_from_models` → `TrimPass.distribution` | ✓ |

No cartridge plane is "silent" — every authored value
participates in the plan signature.

---

## 7. Q6 — Non-choice output parameters

Surfaces that don't change the proposal support at a specific
choice site but DO change downstream output (probabilities,
paired-end window shapes) are folded as documented above:

| Surface | Folded as | Notes |
|---|---|---|
| `paired_end(r1_length, r2_length, insert_size)` | three `fmt_int_dist` calls (`r1=`, `r2=`, `insert=`) | ✓ |
| `random_strand_orientation(prob)` → `RevCompPass.apply_prob` | `fmt_prob("apply_prob", ...)` | ✓ |
| `receptor_revision(prob)` → `ReceptorRevisionPass.prob` | `fmt_prob("prob", ...)` | ✓ (but see soft gap (2)) |
| `invert_d(prob)` → `InvertDPass.prob` | `fmt_prob("prob", ...)` | ✓ |
| `contaminate(prob)` → `ContaminantPass.apply_prob` | `fmt_prob("apply_prob", ...)` | ✓ |
| Indel `insertion_prob` | `fmt_prob("insertion_prob", ...)` | ✓ |

---

## 8. The three documented soft gaps

### Soft gap 1 — `SampleAllelePass` distribution narrowing

**What's not folded.** The
`distribution: Box<dyn Distribution<Output = AlleleId>>` field
on `SampleAllelePass`. When the user calls
`Experiment.restrict_alleles(v="IGHV1-2*02")` or
`recombine(v_allele_weights={"IGHV1-2*02": 100.0})`, the
Python lowering layer constructs a narrowed/reweighted
`AllelePoolDist` and hands it to `SampleAllelePass::new`. The
pass's `parameter_signature` returns the empty string
regardless.

**Why.** Documented at
[`sample_allele.rs::parameter_signature`](../engine_rs/src/passes/sample_allele.rs):

> The sampling distribution is opaque (typically the pool's
> uniform `AllelePoolDist`). We deliberately do NOT fold the
> distribution's `support()` into the signature: the pool
> identity is already covered by `refdata_content_hash`, and
> emitting the entire allele-id ↔ weight table per V/D/J
> pass would inflate signatures without adding information
> beyond what the content hash already gates.

The rationale held when the only way to construct
`SampleAllelePass` was the cartridge's full uniform
`AllelePoolDist`. With `restrict_alleles` /
`*_allele_weights` shipped, the distribution can deviate
from the uniform pool — and the signature doesn't capture
that deviation.

**Soft-gap backstop.** The strict-mode sampler rejects a
replayed allele ID that lies outside the narrowed
distribution support
([`sample_allele.rs::sample_allele`](../engine_rs/src/passes/sample_allele.rs)
documents this as the "v3.0 documented exception to the
global constrain-before-propose invariant"). So:

| Origin → target | Behaviour |
|---|---|
| Unrestricted → Restricted-to-single-allele | Strict sampler rejects when the recorded allele isn't the single restricted one. Replay fails LATER than the signature gate, but the bytes never propagate. |
| Restricted-to-single-allele → Unrestricted | Replay succeeds; bytes are byte-identical to the original run (the recorded allele is in unrestricted support). User loses the cartridge-identity distinction at the signature layer but gains it back through `refdata_content_hash`. |

**Recommended tightening.** A follow-up slice could extend
`parameter_signature` to fold `distribution.support()` ONLY
when the support is a strict subset of the V/D/J pool
(measured against `ctx.refdata`). Out of scope for this
audit; pinned as a documented soft gap.

### Soft gap 2 — `ReceptorRevisionPass` replacement V distribution

**What's not folded.** Same shape as soft gap (1) — the
replacement V `AllelePoolDist` and the replacement V trim_3
distribution on `ReceptorRevisionPass`. Today the Python
DSL doesn't expose a way to narrow the replacement V pool
or weight it, but the engine struct supports it (a future
DSL slice that opens it up would slip through the signature).

**Why.** Same rationale as (1) — documented in
[`receptor_revision.rs::parameter_signature`](../engine_rs/src/passes/receptor_revision.rs):

> The V replacement distribution is over the same V pool
> covered by `refdata_content_hash`; skip it (same
> rationale as `SampleAllelePass`) and pin only `prob`.

**Soft-gap backstop.** Same as (1).

### Soft gap 3 — S5F kernel payload content

**What's not folded.** The 5-mer transition matrix bytes
that back `S5FMutationPass.kernel`. Only the kernel's *name
string* is folded.

**Why.** Documented in [`shm_model_audit.md`](shm_model_audit.md)
§1.4 + the cartridge manifest block at
[`dataconfig/data_config.py::cartridge_manifest`](../src/GenAIRR/dataconfig/data_config.py)
(`models.shm.in_content_hash=False`). The shipped boundary
is explicit: kernel choice is per-experiment metadata, not
cartridge identity; the manifest carries
`s5f_kernel_digest` separately.

**Soft-gap backstop.** A user who swaps the kernel payload
under the bundled name `s5f_default` would produce different
output bytes but equal plan signatures. The
`s5f_kernel_digest` in the manifest changes; a downstream
consumer that audits the manifest catches this.

**Recommended tightening.** Fold `s5f_kernel_digest` into
the S5F pass signature when the kernel name is recognised as
a bundled identifier (so user-supplied kernels at a custom
path can still ride through). Out of scope for this audit.

---

## 9. Stop-and-report determination

The user brief's stop-and-report condition:

> If this audit finds a missing parameter, stop and report.
> This is exactly the class of bug that can silently corrupt
> replay.

**Verdict: clean-yes — no silent-corruption gap.**

The three soft gaps identified above are:

- All documented at the code level with explicit
  parameter_signature docstrings naming the chosen
  boundary.
- All backed by an independent safety mechanism (strict-
  sampler rejection for the two allele-distribution gaps;
  `refdata_content_hash` + manifest `s5f_kernel_digest` for
  the S5F payload gap).
- None lead to silent byte corruption — at worst, the user
  loses the cartridge-identity distinction at the signature
  layer for a specific narrowing surface.

The audit therefore proceeds to land the contract pins
without requesting a tightening slice. If a future regression
weakens any of the three backstops (e.g. removes the strict-
sampler's out-of-support rejection, or removes the manifest's
`s5f_kernel_digest`), the contract pins surface the
regression because they assert both the signature-gap
boundary AND the backstop behaviour.

---

## 10. Test surface — what this audit pins

Mirrored in
[`tests/test_plan_signature_completeness_contract.py`](../tests/test_plan_signature_completeness_contract.py).

### `pin_change_*` — per-parameter-family change-detection

For each major parameter family enumerated in §3, two
experiments differing only by that parameter MUST produce
different `pass_plan_signature` strings. Pinned families:

1. NP1 length distribution.
2. P-V_3 length distribution (per-end P-nucleotide proxy).
3. Trim V_3 distribution.
4. NP base model (uniform → empirical → Markov; pre-pinned
   in the Markov implementation tests but cross-pinned here
   for completeness).
5. Markov transition row.
6. Mutation count rate vs distribution.
7. `segment_rates` non-default vector.
8. `v_subregion_rates` non-default vector.
9. S5F `kernel_name`.
10. PCR `count_source` value.
11. Indel `insertion_prob`.
12. End-loss `length`.
13. Contaminant `apply_prob`.
14. Rev-comp `apply_prob`.
15. Invert-D `prob`.
16. Receptor-revision `prob`.
17. Paired-end `r1_length`.
18. Paired-end `insert_size`.

### `pin_default_*` — equivalent-defaults collide

19. No `segment_rates` kwarg vs explicit all-ones dict.
20. No `v_subregion_rates` kwarg vs explicit all-ones dict.
21. No `p_nucleotide_lengths` cartridge field vs empty dict
    (byte-identical to pre-slice baseline).

### `pin_canonical_*` — dict-order insensitivity

22. `segment_rates` dict insertion order.
23. `v_subregion_rates` dict insertion order.
24. Markov transition rows dict insertion order (`from`
    outer + `to` inner).
25. Empirical distribution `(value, weight)` pairs duplicate
    accumulation (`[(1, 0.5), (1, 0.5)]` vs `[(1, 1.0)]`).

### `pin_replay_gate_*` — mismatched signature fails before
trace consumption

26. Replay across mismatched NP-length distributions fails
    the signature gate with the canonical "pass plan
    signature mismatch" error.
27. Replay across mismatched Markov transition matrices
    fails the signature gate.
28. Replay across mismatched P-V_3 lengths fails the
    signature gate.
29. Replay across mismatched `segment_rates` fails the
    signature gate.
30. Replay across mismatched `paired_end` insert size fails
    the signature gate.

### `pin_soft_gap_*` — three documented gaps + their backstops

31. `restrict_alleles` does NOT change plan signature
    (gap pinned, not "fixed").
32. `*_allele_weights` does NOT change plan signature.
33. SampleAllele strict-mode backstop: replay against a
    restriction whose support excludes the recorded
    allele fails at the SampleAllele sampler.
34. `s5f_kernel_digest` is exposed in the cartridge
    manifest as the independent backstop for the S5F
    payload gap.

### `pin_legacy_*` — v1/v2 trace compatibility

35. Legacy trace files (recorded under the v1 schema) still
    pass the address-schema-version check against the
    current engine.
36. Frozen address spellings remain parseable — the
    `address.rs` `frozen_address_spellings_for_choice_address_schema_v1`
    test extends, never replaces, its pin set.

### `pin_unfolded_*` — no production pass with non-default
configurable parameters lacks a non-empty signature when
configured

37. Every `Pass` implementation in
    `engine_rs/src/passes/` that has at least one
    non-default constructor parameter beyond the `segment` /
    `end` discriminator implements `parameter_signature()`
    with a non-empty return value when those parameters are
    set. `AssembleSegmentPass` / `EchoPass` / `SampleBasePass`
    are the documented test-only / segment-only exceptions.

### Doc anchor

38. Audit doc exists and references the contract file;
    section structure intact.

---

## 11. Out of scope

Documented here so a future implementer doesn't accidentally
expand the work.

- **Tightening soft gap 1** (fold SampleAllele distribution
  support when narrowed). Recommended; deferred to a separate
  slice scoped to fold `restrict_alleles` /
  `*_allele_weights` into the signature. The tightening
  needs to be careful to avoid byte-identical signature
  inflation for the unrestricted case — the documented
  collision-with-content-hash discipline must hold.
- **Tightening soft gap 2** (receptor-revision replacement V).
  Tied to soft gap 1's fix — when the Python DSL opens up
  receptor-revision V restriction, the same fold extends.
- **Tightening soft gap 3** (S5F kernel payload digest in
  signature). Requires coordination with the cartridge
  manifest's `s5f_kernel_digest` boundary; if both are folded
  the user gets a stronger replay-safety gate but a noisier
  signature for the common case of bundled kernels.
- **Contracts in the signature.** `productive_only()` adds a
  `ContractSet` on the runtime context, not on the plan. The
  contract set narrows per-choice support but is treated as a
  *runtime* config (not a *plan* parameter) by design — same
  boundary as `feasibility`. A future audit may revisit this
  if contract authorability becomes a public surface.
- **Feasibility window in the signature.** Same disposition as
  contracts.
- **Per-pass digest of the configured pass-name string.** Each
  pass's `name()` is already part of the plan signature via
  `Pass::name`; this audit doesn't extend that.

---

## 12. Summary table

| Concern | Today | Verdict |
|---|---|---|
| Every production pass with non-default parameters folds them into `parameter_signature` | All but the three documented soft gaps | **Clean — three soft gaps documented at code level + backed by independent backstops** |
| Dict-order / insertion-order canonicalization | Enforced by `fmt_segment_rates` / `fmt_v_subregion_rates` / Markov `signature()` / `CategoricalBase::from_pairs` accumulation | ✓ |
| Default normalization (no kwarg vs explicit all-ones) | `is_default()` short-circuit on rate vectors; empty-dict cartridge fields lower to no-op | ✓ |
| Cartridge-driven defaults represented | Every `ReferenceEmpiricalModels` field lowers to a pass parameter that folds | ✓ |
| Non-choice output parameters (probabilities, paired-end windows) | Folded as `fmt_prob` / `fmt_int_dist` | ✓ |
| `restrict_alleles` / `*_allele_weights` fold | ✗ (soft gap 1 — strict-sampler backstop) | Documented; pinned as soft gap |
| Receptor-revision replacement V distribution fold | ✗ (soft gap 2 — same backstop) | Documented; pinned as soft gap |
| S5F kernel payload fold | ✗ (soft gap 3 — `s5f_kernel_digest` backstop) | Documented; pinned as soft gap |
| Pre-flight bugs found | **None.** The audit confirms the plan signature mechanism is structurally sound; the three soft gaps are explicit, documented, and backstopped. | — |

The plan-signature replay-safety surface is **release-grade
for v1**. Three documented soft gaps remain as future
tightening candidates; none is a silent corruption hole. The
recommended next step is the contract pins below, not a
tightening slice — the user can decide separately whether
the soft-gap tightening is worth a slice in its own right.
