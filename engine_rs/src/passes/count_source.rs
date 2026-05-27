//! `CountSource` — the count abstraction used by
//! [`super::UniformMutationPass`] and [`super::S5FMutationPass`].
//!
//! Two variants:
//!
//! - **`Distribution`** — the classic shape: an explicit
//!   `Distribution<Output = i64>` (typically empirical or a fixed
//!   value) sampled once per pass execution. The pool length is
//!   ignored. This is the path the v1 `mutate(count=...)` DSL took
//!   and the only path before v2.0.
//! - **`Rate`** — a per-base mutation rate. Per execution, the
//!   pass samples `count ~ Poisson(rate * pool_len)`. This is the
//!   biologically natural way to specify SHM intensity: a rate of
//!   `0.03` reads as "3 % of bases mutated, on average," which
//!   matches how immunologists report SHM in the literature.
//!
//! The rate-mode samples per-record against the *current pool
//! length*, not against a refdata-mean length, so trimmed records,
//! VJ vs VDJ chains, and per-record stochastic length differences
//! all see a rate that's a true fraction of their own length. This
//! is intentional — the architecture review explicitly rejected
//! compile-time conversion via mean refdata length on the grounds
//! that it "silently breaks for non-IGH and trimmed records."

use crate::address::ChoiceAddress;
use crate::dist::Distribution;
use crate::pass::{PassContext, PassError};
use crate::rng::Rng;
use crate::trace::ChoiceValue;

pub enum CountSource {
    /// Empirical / explicit count distribution. Sampled once per
    /// pass execution; the result is the literal mutation count.
    Distribution(Box<dyn Distribution<Output = i64>>),
    /// Per-base mutation rate (e.g. `0.03` for 3 % SHM). Per pass
    /// execution, the count is drawn from `Poisson(rate * pool_len)`
    /// against the current pool length.
    Rate(f64),
}

impl CountSource {
    /// Sample the mutation count for one execution of the pass.
    ///
    /// `pool_len` is consulted only in the `Rate` variant. Returns
    /// an `i64` to match the existing `Distribution<Output = i64>`
    /// shape — the caller does the same negative / overflow checks
    /// as before.
    pub fn sample(&self, rng: &mut Rng, pool_len: u32) -> i64 {
        match self {
            Self::Distribution(d) => d.sample(rng),
            Self::Rate(rate) => {
                let lambda = (*rate) * (pool_len as f64);
                poisson_sample(rng, lambda) as i64
            }
        }
    }

    /// True iff this source samples against the pool length. Used
    /// by the trace / describe layers to label the count event.
    #[allow(dead_code)]
    pub fn is_rate(&self) -> bool {
        matches!(self, Self::Rate(_))
    }
}

/// Sample a mutation count for a count-driven pass, validate it,
/// and record it to the trace at `count_address`.
///
/// Shared prelude for the count-driven passes: uniform / S5F / PCR /
/// quality / N-corrupt. Each pass continues with its own pool-length
/// / zero-count short-circuit, transaction open, and per-step loop
/// after this returns.
///
/// Behavior:
/// - **Replay mode** (`ctx.replay_cursor.is_some()`): consumes the
///   recorded `Int` at `count_address` from the cursor instead of
///   drawing from `count_source`. The downstream validation path
///   (strict-mode range checks, debug asserts, trace record) runs
///   unchanged so a bad recorded count surfaces as
///   `InvalidDistributionOutput` rather than corrupting the IR.
/// - **Fresh-RNG mode** (`ctx.replay_cursor.is_none()`): draws
///   `count_raw: i64` from `count_source` (`Rate` mode consults
///   `pool_len`; `Distribution` mode ignores it).
/// - **Strict mode**: returns `PassError::InvalidDistributionOutput`
///   if `count_raw` is negative or exceeds `u32::MAX`.
/// - **Permissive mode**: the same out-of-range values trip
///   debug asserts (existing per-pass behavior).
/// - Records `ChoiceValue::Int(count_raw)` at `count_address`.
/// - Returns `count_raw as u32` on success.
pub(crate) fn sample_validated_count(
    count_source: &CountSource,
    ctx: &mut PassContext,
    pool_len: u32,
    pass_name: &str,
    count_address: ChoiceAddress,
    strict: bool,
) -> Result<u32, PassError> {
    // Trace-injected replay (Tier 2): consume the recorded count
    // from the cursor instead of drawing from `count_source`. The
    // post-consume validation path is identical so bad recorded
    // values still surface as `InvalidDistributionOutput`.
    let count_raw = if let Some(cursor) = ctx.replay_cursor.as_deref_mut() {
        cursor
            .expect_int(count_address)
            .map_err(|reason| PassError::replay(pass_name, reason))?
    } else {
        count_source.sample(ctx.rng, pool_len)
    };
    if strict && count_raw < 0 {
        return Err(PassError::invalid_distribution_output(
            pass_name,
            count_address.to_string(),
            count_raw,
            "negative_count",
        ));
    }
    if strict && count_raw > u32::MAX as i64 {
        return Err(PassError::invalid_distribution_output(
            pass_name,
            count_address.to_string(),
            count_raw,
            "count_exceeds_u32",
        ));
    }
    assert!(
        count_raw >= 0,
        "{}: count distribution returned negative {}",
        pass_name,
        count_raw
    );
    assert!(
        count_raw <= u32::MAX as i64,
        "{}: count distribution returned {} > u32::MAX",
        pass_name,
        count_raw
    );
    ctx.trace
        .record_choice(count_address, ChoiceValue::Int(count_raw));
    Ok(count_raw as u32)
}

/// Sample a `Poisson(lambda)`-distributed value using Knuth's
/// algorithm. Numerically stable for the regime expected by this
/// engine (`lambda` between 0 and ~50, since SHM rates are <10%
/// and sequence lengths are <500bp). For `lambda <= 0`, returns 0.
///
/// Reference: Knuth, *The Art of Computer Programming*, vol 2,
/// 3.4.1. Mean / variance are both `lambda`. Runtime is `O(lambda)`
/// expected, so the algorithm is unsuitable for large `lambda` —
/// but for the genuinely-immunological regime where rate × len is
/// at most a few dozen, it's the right pick.
fn poisson_sample(rng: &mut Rng, lambda: f64) -> u32 {
    if !lambda.is_finite() || lambda <= 0.0 {
        return 0;
    }
    let exp_neg_lambda = (-lambda).exp();
    let mut k: u32 = 0;
    let mut p: f64 = 1.0;
    loop {
        k += 1;
        p *= rng.next_f64();
        if p <= exp_neg_lambda {
            return k - 1;
        }
        // Safety net for pathological `lambda` (large): cap at a
        // very generous bound so we never spin forever. The real
        // engine guard is the upstream `count_raw > u32::MAX` check
        // in the pass; here we just stop the loop.
        if k > 1_000_000 {
            return k;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::EmpiricalLengthDist;

    fn fixed(n: i64) -> CountSource {
        CountSource::Distribution(Box::new(EmpiricalLengthDist::from_pairs(vec![(n, 1.0)])))
    }

    #[test]
    fn distribution_variant_ignores_pool_len() {
        let mut rng = Rng::new(42);
        let src = fixed(7);
        assert_eq!(src.sample(&mut rng, 1), 7);
        assert_eq!(src.sample(&mut rng, 1_000_000), 7);
    }

    #[test]
    fn rate_variant_returns_zero_on_zero_length() {
        let mut rng = Rng::new(42);
        let src = CountSource::Rate(0.05);
        assert_eq!(src.sample(&mut rng, 0), 0);
    }

    #[test]
    fn rate_variant_returns_zero_for_zero_rate() {
        let mut rng = Rng::new(42);
        let src = CountSource::Rate(0.0);
        assert_eq!(src.sample(&mut rng, 300), 0);
    }

    #[test]
    fn rate_variant_mean_matches_lambda_in_expectation() {
        // Run many samples and check the empirical mean is close
        // to lambda. lambda = 0.03 * 300 = 9 → expect 1000-sample
        // mean within ~10 % of 9.
        let mut rng = Rng::new(42);
        let src = CountSource::Rate(0.03);
        let pool_len: u32 = 300;
        let n = 1000usize;
        let total: i64 = (0..n).map(|_| src.sample(&mut rng, pool_len)).sum();
        let mean = total as f64 / n as f64;
        let expected = 0.03 * (pool_len as f64);
        let tolerance = expected * 0.15; // 15 % envelope on the Monte Carlo estimate
        assert!(
            (mean - expected).abs() < tolerance,
            "empirical mean {mean} differs from expected {expected} by more than {tolerance}"
        );
    }

    #[test]
    fn rate_variant_is_deterministic_under_same_seed() {
        let src = CountSource::Rate(0.05);
        let mut a = Rng::new(123);
        let mut b = Rng::new(123);
        for _ in 0..50 {
            assert_eq!(src.sample(&mut a, 250), src.sample(&mut b, 250));
        }
    }

    #[test]
    fn is_rate_discriminates_variants() {
        assert!(CountSource::Rate(0.05).is_rate());
        assert!(!fixed(8).is_rate());
    }

    // ── Trace-injected replay (Tier 2, count-only) ────────────────
    //
    // Direct unit tests of `sample_validated_count` so the count
    // consume path is verified without relying on full-plan
    // execution (the strict-positional cursor blocks count-pass
    // e2e replay during staged migration when unmigrated indexed
    // children come between the count and any earlier migrated
    // record). Tests construct a `PassContext` with a cursor in
    // isolation and assert:
    //  1. Replay consumes the recorded count and bypasses the RNG.
    //  2. Bad recorded counts still surface as
    //     `InvalidDistributionOutput` — proving the shared
    //     validation path runs.
    //  3. Wrong-kind / exhausted cursors surface as
    //     `PassError::Replay` with the expected `ReplayError`
    //     reason.

    use crate::pass::{PassContext, PassError};
    use crate::replay::{ReplayError, TraceCursor};
    use crate::trace::{ChoiceRecord, ChoiceValue, Trace};

    fn run_with_replay(
        records: Vec<ChoiceRecord>,
        source: CountSource,
        pool_len: u32,
        addr: ChoiceAddress,
        strict: bool,
    ) -> Result<(u32, Trace, u64), PassError> {
        // Use a deterministic seed so we can detect whether the
        // helper drew from the RNG (it shouldn't, in replay mode).
        let mut cursor = TraceCursor::from_owned(records);
        let mut trace = Trace::new();
        let mut rng = Rng::new(0xc0ff_ee);
        let result;
        let rng_words_after;
        {
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index: 0,
                refdata: None,
                contracts: None,
                feasibility: None,
                reference_index: None,
                replay_cursor: Some(&mut cursor),
                event_log_sink: None,
            };
            result = sample_validated_count(&source, &mut ctx, pool_len, "test.pass", addr, strict);
        }
        rng_words_after = rng.words_consumed();
        result.map(|count| (count, trace, rng_words_after))
    }

    fn rec(addr: ChoiceAddress, value: ChoiceValue) -> ChoiceRecord {
        ChoiceRecord::new(addr.to_string(), value)
    }

    #[test]
    fn helper_replay_consumes_recorded_count() {
        let addr = ChoiceAddress::MutateUniformCount;
        // `fixed(99)` would draw 99 under RNG mode; the cursor
        // overrides with 17. The returned count must be 17.
        let (count, trace, rng_words) =
            run_with_replay(vec![rec(addr, ChoiceValue::Int(17))], fixed(99), 0, addr, true)
                .unwrap();
        assert_eq!(count, 17);
        // The helper still records to the output trace.
        let rec = trace.find("mutate.uniform.count").unwrap();
        assert_eq!(rec.value, ChoiceValue::Int(17));
        // The RNG was not touched — replay bypassed `count_source.sample`.
        assert_eq!(rng_words, 0);
    }

    #[test]
    fn helper_replay_consumes_zero_count_record() {
        let addr = ChoiceAddress::CorruptPcrCount;
        let (count, _trace, rng_words) =
            run_with_replay(vec![rec(addr, ChoiceValue::Int(0))], fixed(5), 0, addr, true)
                .unwrap();
        assert_eq!(count, 0);
        assert_eq!(rng_words, 0);
    }

    #[test]
    fn helper_replay_negative_recorded_count_surfaces_validation_error() {
        // The invariant: replay consumes the recorded value, then
        // runs the same validation path as fresh sampling. A
        // recorded negative count must trip strict-mode
        // `InvalidDistributionOutput`, not corrupt the engine.
        let addr = ChoiceAddress::MutateS5fCount;
        let err = run_with_replay(
            vec![rec(addr, ChoiceValue::Int(-3))],
            fixed(0),
            0,
            addr,
            true,
        )
        .unwrap_err();
        match err {
            PassError::InvalidDistributionOutput { value, reason, .. } => {
                assert_eq!(value, -3);
                assert_eq!(reason, "negative_count");
            }
            other => panic!("expected InvalidDistributionOutput, got {other:?}"),
        }
    }

    #[test]
    fn helper_replay_oversized_recorded_count_surfaces_validation_error() {
        let addr = ChoiceAddress::CorruptQualityCount;
        let too_big: i64 = (u32::MAX as i64) + 1;
        let err = run_with_replay(
            vec![rec(addr, ChoiceValue::Int(too_big))],
            fixed(0),
            0,
            addr,
            true,
        )
        .unwrap_err();
        match err {
            PassError::InvalidDistributionOutput { value, reason, .. } => {
                assert_eq!(value, too_big);
                assert_eq!(reason, "count_exceeds_u32");
            }
            other => panic!("expected InvalidDistributionOutput, got {other:?}"),
        }
    }

    #[test]
    fn helper_replay_wrong_value_kind_surfaces_replay_error() {
        // Cursor has Base instead of Int — kind mismatch.
        let addr = ChoiceAddress::MutateUniformCount;
        let err = run_with_replay(
            vec![rec(addr, ChoiceValue::Base(b'A'))],
            fixed(0),
            0,
            addr,
            true,
        )
        .unwrap_err();
        match err {
            PassError::Replay { pass_name, reason } => {
                assert_eq!(pass_name, "test.pass");
                assert!(matches!(reason, ReplayError::ValueKindMismatch { .. }));
            }
            other => panic!("expected PassError::Replay, got {other:?}"),
        }
    }

    #[test]
    fn helper_replay_exhausted_cursor_surfaces_replay_error() {
        let addr = ChoiceAddress::MutateUniformCount;
        let err =
            run_with_replay(Vec::new(), fixed(5), 0, addr, true).unwrap_err();
        match err {
            PassError::Replay { reason, .. } => {
                assert!(matches!(reason, ReplayError::Exhausted { .. }));
            }
            other => panic!("expected PassError::Replay, got {other:?}"),
        }
    }

    #[test]
    fn helper_replay_address_mismatch_surfaces_replay_error() {
        // Cursor has the wrong typed address for what the pass asks.
        let cursor_addr = ChoiceAddress::CorruptPcrCount;
        let expected = ChoiceAddress::MutateUniformCount;
        let err = run_with_replay(
            vec![rec(cursor_addr, ChoiceValue::Int(3))],
            fixed(0),
            0,
            expected,
            true,
        )
        .unwrap_err();
        match err {
            PassError::Replay { reason, .. } => {
                assert!(matches!(reason, ReplayError::AddressMismatch { .. }));
            }
            other => panic!("expected PassError::Replay, got {other:?}"),
        }
    }

    #[test]
    fn helper_fresh_rng_path_still_works_when_cursor_absent() {
        // Direct call without replay_cursor (parallel test path).
        let addr = ChoiceAddress::MutateUniformCount;
        let source = fixed(11);
        let mut trace = Trace::new();
        let mut rng = Rng::new(0);
        let count = {
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index: 0,
                refdata: None,
                contracts: None,
                feasibility: None,
                reference_index: None,
                replay_cursor: None,
                event_log_sink: None,
            };
            sample_validated_count(&source, &mut ctx, 0, "test.pass", addr, true).unwrap()
        };
        assert_eq!(count, 11);
        assert_eq!(
            trace.find("mutate.uniform.count").unwrap().value,
            ChoiceValue::Int(11)
        );
    }
}
