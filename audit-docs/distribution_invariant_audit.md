# Distribution-invariant audit

**Status:** audit + golden tests (statistical). Drift items in §6.

This audit proves the engine isn't merely producing **legal** sequences
(the [productive_failure_mode](productive_failure_mode_audit.md) and
[productive_stress_matrix](../tests/test_productive_stress_matrix.py)
audits cover that), but sampling legal sequences with the **intended
conditional distribution**.

The unifying rule:

> For every constrained sampler, the probability of drawing value `x`
> is
>
> `Pr[x] ∝ natural_weight(x) * 1[contracts.admit(x)]`
>
> with the proportionality constant determined by inverse-CDF
> renormalization against the surviving mass.

That is — narrow the natural distribution to its admissible subset,
then **renormalize**. Never fall back to a different distribution,
never rejection-sample then uniform.

---

## 1. The central helper

[`sample_filtered_result`](../engine_rs/src/dist/filtered.rs#L65)
implements the rule. The algorithm (file:line 84–98):

1. Walk the natural distribution's enumerable support.
2. Apply the admissibility predicate to each candidate; collect the
   survivors with their **natural weights** (not 1.0).
3. Sum the surviving weights → `total`.
4. If `total <= 0` or non-finite, return `EmptyAdmissibleSupport` /
   `InvalidFilteredSupport`.
5. Otherwise inverse-CDF sample on the surviving (value, weight)
   list, with the implicit total as the normalizer.

The inverse-CDF step is the renormalization — there's no separate
divide-by-total operation, just `r = rng.next_f64() * total` and a
running cumulative sum. Mathematically equivalent to drawing from
the natural distribution conditioned on admissibility.

Two specialised paths use the same algorithm:

- [`sample_base_with_admit_mask`](../engine_rs/src/dist/filtered.rs#L220)
  — 4-bit base mask fast path for NP base sampling. Avoids `Vec`
  allocation by iterating the support twice in-place.
- [`sample_admissible_tuple`](../engine_rs/src/passes/corrupt/indel/tuple_sampler.rs)
  — multi-event indel sampler with **continuation weighting** (§3).

---

## 2. Per-mechanism implementation manifest

For each sampling slot in the engine, the (natural, admissibility,
implementation) triple:

| Mechanism                          | Natural distribution                                  | Admissibility filter                                            | Implementation                                                            |
|------------------------------------|--------------------------------------------------------|------------------------------------------------------------------|---------------------------------------------------------------------------|
| NP length (`generate_np`)          | `length_dist` (DSL-supplied empirical)                 | `contracts.admits_typed(Int(candidate))` (frame mod-3)           | Central `sample_filtered_with_policy`                                      |
| NP base — fast path                | `base_dist` (uniform A/C/G/T by default)               | 4-bit admit mask from `ProductiveAdmitMaskObserver`             | `sample_base_with_admit_mask` (in-place inverse CDF)                       |
| NP base — slow path                | `base_dist`                                            | `contracts.admits_typed(Base(candidate))`                        | Central `sample_filtered_with_policy`                                      |
| Trim length (V/D/J)                | `distribution` (DSL-supplied)                          | `LengthSupport` upper bound + feasibility filter                 | Central `sample_filtered_with_policy`                                      |
| SHM site + base (S5F)              | mutability × substitution-row weight per (site, base)  | per-site 4-bit admit mask                                        | Bespoke `sample_weighted_event` (same inverse-CDF algorithm)               |
| PCR error site                     | uniform on `[0, pool_len)`                             | per-site mask via `MutationTransaction` (slot rejected if mask empty) | Routed through `substitute_base_admitting`                                 |
| PCR error base                     | uniform A/C/G/T                                        | per-site mask                                                    | Routed through `substitute_base_admitting`                                 |
| Sequencing error site              | uniform on `[0, pool_len)`                             | per-site mask                                                    | Mechanically identical to PCR (lowercased before write)                    |
| Sequencing error base              | uniform A/C/G/T                                        | per-site mask                                                    | Same                                                                       |
| N-corruption site                  | uniform on `[0, pool_len)`                             | per-site mask with fixed base `N`                                | Per-site mask query; permissive skip on rejection                          |
| Indel kind/site/base tuple         | per-step independent (insertion_prob, uniform site, base_dist) | mod-3 DP for tuple-level frame closure + post-event contract verify | `sample_admissible_tuple` with **continuation weighting** (§3)             |
| Allele recombination               | `AllelePoolDist` per segment                            | NONE — documented v3.0 exception                                 | Central helper with always-true predicate                                  |

**Audit note on allele recombination**: this is the only sampling
slot in the engine that doesn't apply contract narrowing. The
documented rationale (see
[`sample_allele.rs`](../engine_rs/src/passes/sample_allele.rs)
near line 84) is that allele selection happens before V/J anchor
positions are known to the contract bundle, so per-allele admissibility
would be vacuous. Downstream passes (NP, trim, mutation) carry the
contract enforcement.

---

## 3. Continuation weighting (indels)

Indels are the unique mechanism where the constraint is **tuple-level**
(net frame change ≡ 0 mod 3 across all events), not per-event. A naive
"sample each event independently, then check post-state" would
oversample tuples that can never close the frame.

The
[`sample_admissible_tuple`](../engine_rs/src/passes/corrupt/indel/tuple_sampler.rs)
function implements **continuation weighting**: at step `k` of a
`count`-event tuple, the effective weight of candidate `c` is

```
effective_weight(c) = natural_weight(c) * dp[count - k - 1][required_residual_delta(c)]
```

where `dp[m][r]` is the number of `m`-step continuations whose
delta sum ≡ `r` mod 3. The DP is reseeded each step from the
current sim's per-step delta-mass profile.

Without this DP multiplier, the per-step distribution would be
biased — e.g., for count=2 under productive, naïvely sampling step
0 from the unconstrained distribution would yield
`Pr[step_0_delta = +1] ≈ 2/9 ≈ 0.2222` (since 2 of 9 sites are
FrameDelta(+1) on the baseline fixture). The corrected
continuation-weighted distribution yields `≈ 0.1353` — a
distinguishing 4σ gap at N=4000. The Rust test
[`indel_pass_count_2_matches_exact_continuation_weighting`](../engine_rs/src/passes/corrupt/indel/tests/constraints.rs)
pins this precisely.

This audit's load-bearing Python tests reproduce the indel count=1
and count=2 distribution checks, plus add coverage for the other
mechanisms.

---

## 4. The audit's positive + negative controls

For each mechanism, we want **both**:

- **Positive control**: with admissible support `S` and natural
  weights `w(x)`, the empirical distribution matches
  `w(x) / sum_{y∈S} w(y)`.
- **Negative control**: a buggy implementation (filter without
  renormalize, or rejection-sample then uniform) would produce a
  *different* observable distribution, and our test would fail.

A clean way to construct a negative control is via **non-uniform
natural distributions**: if natural is `[(0, 2.0), (3, 1.0)]` and
both are admissible, the correct conditional is `[2/3, 1/3]` while
a "filter then uniform" bug would give `[1/2, 1/2]` — distinguishable
at high confidence with modest N.

---

## 5. Per-mechanism distribution invariants (pinned by tests)

| # | Invariant                                                                                                                | Method                  | Pinned by                                                                                                   |
|---|--------------------------------------------------------------------------------------------------------------------------|-------------------------|-------------------------------------------------------------------------------------------------------------|
| D1| **NP length narrowing**: uniform natural over `{0..5}`, admissible `{0, 3}` → empirical `Pr[0] ≈ Pr[3] ≈ 0.5`.            | MC, 5σ                  | `test_np_length_narrowing_under_productive_is_uniform_on_admissible`                                        |
| D2| **NP length renormalization**: natural `[(0, 2.0), (3, 1.0)]`, both admissible → empirical `Pr[0] = 2/3, Pr[3] = 1/3`.   | MC, 5σ (negative control) | `test_np_length_renormalization_preserves_natural_weights`                                                  |
| D3| **PCR site distribution**: unconstrained → uniform over `[0, pool_len)`.                                                 | MC, 5σ                  | `test_pcr_site_distribution_is_uniform`                                                                     |
| D4| **Indel kind under productive (count=1)**: empirical `Pr[insertion] ≈ 0.5517` (Rust oracle).                              | MC, 5σ                  | `test_indel_count_1_kind_distribution_matches_rust_oracle`                                                  |
| D5| **Indel continuation weighting (count=2)**: empirical `Pr[step_0_delta = +1] ≈ 0.1353` (Rust oracle).                     | MC, 5σ (negative control: would be ≈ 0.2222 without DP) | `test_indel_count_2_continuation_weighting_matches_rust_oracle`                                              |
| D6| **SHM site narrowing under productive**: positions inside V/J anchor codons have measurably lower mutation rate than unconstrained. | MC qualitative          | `test_shm_under_productive_narrows_v_anchor_codon_substitutions`                                            |
| D7| **PCR base distribution**: unconstrained substitution to uniform A/C/G/T at a single site.                              | MC, 5σ                  | `test_pcr_substitution_base_distribution_is_uniform`                                                        |
| D8| **No drift on existing oracle tests**: Rust-side `indel_pass_count_1_..._matches_exact_enumeration` and `..._continuation_weighting` continue to pass. | Existing tests          | (Rust suite, run by CI)                                                                                     |
| D9| **D inversion is a Bernoulli draw**: at `prob=0.25` the empirical inversion rate over `N=4000` seeds lands within 5σ of `0.25`. Negative-controls the always-true / always-false / sign-swap regressions in `InvertDPass.execute_with_sampling_mode`. | MC, 5σ                  | `test_invert_d_bernoulli_draw_matches_prob_within_5_sigma`                                                  |

5σ for a Bernoulli proportion at p, N: `5 * sqrt(p * (1-p) / N)`.
For p ≈ 0.5, N = 4000: 5σ ≈ ±0.0395 (±4%). For p ≈ 0.55, N = 4000:
5σ ≈ ±0.0393. The Rust oracle tests use the same threshold.

---

## 6. Drift identified

### 6.1 No central distribution-invariant manifest

The architectural rule "natural distribution restricted to admissible
support, renormalized via inverse-CDF" is implemented in three
distinct code paths:

1. Central `sample_filtered_result` (used by NP length, NP base slow,
   trim, allele).
2. Fast path `sample_base_with_admit_mask` (NP base fast).
3. Bespoke `sample_weighted_event` for S5F.
4. Bespoke `sample_admissible_tuple` with DP continuation for indels.

A future refactor that changes one path's implementation could silently
diverge from the others without a test failing. The Rust-side oracle
tests catch divergences from the **theoretical** distribution but
not divergences *between* implementations.

**Possible fix:** add a unit-level invariant test that asserts every
in-tree filtered-sample callsite produces the same conditional
distribution as the central helper would, given identical inputs.
Out of scope for this slice — the existing tests catch theoretical
drift in each path independently.

### 6.2 No `cdr3`-relative substitution rate

S5F's natural weight per (site, base) is `mutability × substitution_row_weight`.
Under productive, the per-site mask varies by junction frame and
position. A potential consumer of "expected SHM rate inside CDR3"
must integrate mutability over admitted sites — there's no
denormalised "effective mutation rate" field today.

**Possible fix:** expose an `effective_mutation_rate` per record
(integrated mutability × admit-mass over the sequence). Deferred —
the AIRR record's `mutation_rate` already exposes the empirical
post-hoc rate.

### 6.3 No Python-side replication of the Rust DP continuation test

The gold-standard continuation-weighting test
(`indel_pass_count_2_matches_exact_continuation_weighting`) lives
only in Rust. Python users who modify the engine and want to verify
the DP wasn't broken would need to run cargo. The Python audit
test (D5) reproduces the *empirical* check (Pr ≈ 0.1353 ± 5σ) but
doesn't re-derive the oracle from first principles in Python.

**Possible fix:** none needed — running both cargo and pytest is
the standard CI path. Flagged for awareness.

---

## 7. Test coverage in this slice

[`tests/test_distribution_invariants.py`](../tests/test_distribution_invariants.py)
pins the §5 invariants:

- **D1**: NP length narrowing (uniform input → uniform on admissible
  subset).
- **D2**: NP length renormalization (non-uniform input → weights
  preserved on admissible subset).
- **D3**: PCR site distribution uniform.
- **D4**: Indel count=1 kind distribution matches Rust oracle.
- **D5**: Indel count=2 step-0 delta distribution matches Rust oracle
  (the continuation-weighting negative control).
- **D6**: SHM under productive narrows V_anchor_codon substitutions
  (qualitative).
- **D7**: PCR base distribution uniform under no contracts.
- **D9**: D-inversion Bernoulli draw matches its declared `prob`
  (Slice E consolidation; pins the `InvertDPass` RNG path).

All tests use Monte Carlo with N=4000 and 5σ tolerance, matching
the Rust oracle tests' precedent. Seeds are deterministic
(`range(4000)`).
