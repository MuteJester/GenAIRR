use super::*;

#[test]
fn indel_pass_productive_filters_insertion_to_structurally_safe_site() {
    // Under v3.0 constrain-before-propose, a single-event indel
    // tuple under productive() must land on a FrameNeutral site
    // outside both anchor codons. With fixture V=TAC (anchor
    // codon at pool [0..3)) + J=TGG (anchor codon at [3..6)):
    //   - Sites 0..3: in V region and ≤ V_anchor_end=3 → Anchor V
    //     marks Forbidden.
    //   - Sites 3..6: in J region and ≤ J_anchor_end=6 → Anchor J
    //     marks Forbidden.
    //   - Site 6: outside both regions → FrameNeutral.
    // So the surviving sample space is {site=6} × 4 bases.
    let (cfg, sim) = make_substitution_productive_vj_fixture();
    let contracts = productive();

    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(1),
        1.0,
        Box::new(UniformBase),
    )));

    let outcome = PassRuntime::execute_with_context(&plan, sim, 0, Some(&cfg), Some(&contracts));

    // An insertion was committed (kind=true) at the only
    // bundle-admissible site for this fixture.
    assert_eq!(
        outcome.trace.find("corrupt.indel.kind[0]").unwrap().value,
        ChoiceValue::Bool(true)
    );
    assert_eq!(
        outcome.trace.find("corrupt.indel.site[0]").unwrap().value,
        ChoiceValue::Int(6)
    );
    // The bundle holds after the tuple commits — frame preserved
    // by the mod-3 DP, no stop introduced (post-event validator).
    assert!(contracts
        .verify(outcome.final_simulation(), Some(&cfg))
        .is_ok());
}

#[test]
fn indel_pass_permissive_reduces_to_noop_when_deletions_cannot_satisfy_contracts() {
    // When every (kind, site) candidate is rejected by the active
    // contract bundle, permissive mode honors the contract by
    // emitting `IndelEvent::NoOp` instead of falling through to
    // unconstrained sampling. The trace still records the slot
    // (kind=false, site=-1) so consumers can audit the reduction,
    // and the final simulation satisfies the contract bundle.
    let (cfg, sim) = make_substitution_productive_vj_fixture();
    let contracts = productive();

    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(1),
        0.0,
        Box::new(UniformBase),
    )));

    let outcome = PassRuntime::execute_with_context(&plan, sim, 0, Some(&cfg), Some(&contracts));

    // The slot was consumed as a NoOp — kind=false (matches the
    // deletion sentinel) and site=-1 (the "no-site" marker).
    assert_eq!(
        outcome.trace.find("corrupt.indel.kind[0]").unwrap().value,
        ChoiceValue::Bool(false)
    );
    assert_eq!(
        outcome.trace.find("corrupt.indel.site[0]").unwrap().value,
        ChoiceValue::Int(-1)
    );
    // The contract bundle holds after permissive reduce-and-skip.
    assert!(
        contracts
            .verify(outcome.final_simulation(), Some(&cfg))
            .is_ok(),
        "permissive NoOp must leave the final simulation contract-admissible"
    );
}

#[test]
fn indel_pass_tuple_sampler_balances_frame_under_productive_count_2() {
    // count=2 with insertion_prob=0.5 gives the mod-3 DP a real
    // tuple to balance. Insertions in junction-shifting positions
    // are FrameDelta(+1); deletions there are FrameDelta(-1).
    // The DP picks tuples whose net delta ≡ 0 (mod 3) — for
    // count=2 that means either two FrameNeutral or a +1/-1 pair.
    // Use a larger fixture so there are enough non-Forbidden
    // sites for the +1/-1 pairing to be realizable.
    let (cfg, sim) = make_productive_indel_balance_fixture();
    let contracts = productive();

    // Run the pass under productive() at a handful of seeds to
    // pin the architectural property: under v3.0 every tuple
    // commits leaves the bundle satisfied.
    for seed in 0..16u64 {
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(2),
            0.5,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute_with_context(
            &plan,
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        );
        assert!(
            contracts
                .verify(outcome.final_simulation(), Some(&cfg))
                .is_ok(),
            "seed {} produced a contract violation",
            seed
        );
    }
}

#[test]
fn indel_pass_tuple_sampler_balances_frame_under_productive_count_3() {
    // count=3 lets the DP exercise a 3-event balance: e.g. three
    // FrameNeutral events, or +1+1+1 (sum 3 ≡ 0 mod 3), or
    // -1-1-1 / +1-1+0 / +1+0-1 etc. Verify the post-state always
    // satisfies the bundle.
    let (cfg, sim) = make_productive_indel_balance_fixture();
    let contracts = productive();
    for seed in 0..16u64 {
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(3),
            0.5,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute_with_context(
            &plan,
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        );
        assert!(
            contracts
                .verify(outcome.final_simulation(), Some(&cfg))
                .is_ok(),
            "count=3 seed {} produced a contract violation",
            seed
        );
    }
}

/// Distribution-invariant test: exact enumeration oracle vs
/// empirical sample frequencies on a small fixture.
///
/// Fixture: V allele "TGTGGG" (anchor 0 → codon TGT = Cys), J
/// allele "TGGAAA" (anchor 0 → codon TGG = Trp). Pool = V[0..6)
/// + J[6..12), junction = [0..9) (length 9, in-frame). Under
/// productive():
///   - Sites 0,1,2: V anchor codon Forbidden.
///   - Sites 3,4,5: FrameDelta(±1) (in V, past anchor).
///   - Sites 6,7,8: J anchor codon Forbidden.
///   - Sites 9,10,11: FrameNeutral (in J, past anchor).
///   - Site 12: FrameNeutral (post-pool insertion).
///
/// For count=1 and insertion_prob=0.5, only delta=0 candidates
/// survive the DP filter (no remaining steps, target=0).
/// Surviving candidates and their natural weights:
///   - 16 FrameNeutral insertions (sites {9,10,11,12} × 4 bases),
///     each weight = 0.5 / (pool_len+1=13) / 4 ≈ 0.009615.
///   - 3 FrameNeutral deletions (sites {9,10,11}),
///     each weight = 0.5 / pool_len=12 ≈ 0.041667.
///   - Total insertion mass: 16 × 0.009615 ≈ 0.15385.
///   - Total deletion mass: 3 × 0.041667 = 0.125.
///   - Total mass: 0.27885.
///   - P(insertion | count=1) = 0.15385 / 0.27885 ≈ 0.5517.
///   - P(deletion | count=1)  = 0.125  / 0.27885 ≈ 0.4483.
///
/// The test runs the sampler over `n` independent seeds and
/// asserts the empirical insertion/deletion ratio matches the
/// theoretical to within a tolerance derived from the normal
/// approximation of the binomial.
#[test]
fn indel_pass_count_1_under_productive_matches_exact_enumeration() {
    let (cfg, sim) = balanced_dp_distribution_fixture();
    let contracts = productive();

    const N: u32 = 4000;
    let mut insertions = 0u32;
    let mut deletions = 0u32;
    for seed in 0..N as u64 {
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(1),
            0.5,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute_with_context(
            &plan,
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        );
        assert!(
            contracts
                .verify(outcome.final_simulation(), Some(&cfg))
                .is_ok(),
            "seed {} violated contracts",
            seed
        );
        match outcome.trace.find("corrupt.indel.kind[0]").unwrap().value {
            ChoiceValue::Bool(true) => insertions += 1,
            ChoiceValue::Bool(false) => deletions += 1,
            _ => panic!("unexpected kind variant"),
        }
    }

    // Theoretical: P(insertion) = 0.5517, P(deletion) = 0.4483.
    // Binomial std for N=4000: sqrt(N * p * (1-p)) ≈ sqrt(4000 *
    // 0.5517 * 0.4483) ≈ 31. ±5 std for a tight pin.
    let expected_insertions = (N as f64) * 0.5517;
    let std = ((N as f64) * 0.5517 * 0.4483).sqrt();
    let tolerance = 5.0 * std;
    let delta = (insertions as f64 - expected_insertions).abs();
    assert!(
        delta < tolerance,
        "empirical insertion count {} deviates from theoretical {} by {} > 5σ ({}); deletions={}",
        insertions,
        expected_insertions,
        delta,
        tolerance,
        deletions
    );
}

/// Distribution-invariant test for count=2: every accepted tuple
/// must have net frame delta ≡ 0 mod 3 (the mod-3 DP's core
/// guarantee). Verified by sampling many tuples and asserting
/// the post-tuple junction length is divisible by 3 against the
/// pre-tuple length.
#[test]
fn indel_pass_count_2_under_productive_preserves_frame_invariant() {
    use crate::junction::compute_junction;

    let (cfg, sim) = balanced_dp_distribution_fixture();
    let contracts = productive();
    let pre_junction = compute_junction(&sim, &cfg).expect("junction defined");
    assert!(pre_junction.is_in_frame());

    for seed in 0..256u64 {
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(2),
            0.5,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute_with_context(
            &plan,
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        );
        let post = outcome.final_simulation();
        // Either the tuple was accepted (junction in frame) or
        // every retry rejected and we collapsed to no-op
        // (junction unchanged from the in-frame input). Both
        // outcomes are admissible under v3.0; what's NOT allowed
        // is an out-of-frame post-state.
        match compute_junction(post, &cfg) {
            Some(j) => assert!(
                j.is_in_frame(),
                "seed {} produced out-of-frame junction (length {})",
                seed,
                j.length
            ),
            None => panic!("seed {} lost the junction window", seed),
        }
        assert!(
            contracts.verify(post, Some(&cfg)).is_ok(),
            "seed {} violated bundle",
            seed
        );
    }
}

/// Distribution-invariant oracle for count=2 — pins
/// **continuation weighting** (not just support filtering).
///
/// On the [`balanced_dp_distribution_fixture`] (V[0..6) +
/// J[6..12) under productive()), the per-step support has:
///   - 12 FrameDelta(+1) insertion candidates @ 1/104 each →
///     mass[+1] = 12/104 = 0.115385
///   - 3 FrameDelta(−1) deletion candidates @ 1/24 each →
///     mass[−1] = 3/24 = 0.125
///   - 16 FrameNeutral insertions + 3 FrameNeutral deletions →
///     mass[0] = 16/104 + 3/24 = 0.278846
///
/// For count=2 the step-0 continuation weight for a candidate
/// with delta `d` is `dp[1][(−d) mod 3]`:
///   - `dp[1][0] = mass[0]` (need next step delta=0)
///   - `dp[1][1] = mass[+1]` (next step must be +1)
///   - `dp[1][2] = mass[−1]` (next step must be −1)
///
/// Step-0 effective mass per delta bucket:
///   - delta=0: mass[0] × dp[1][0] = 0.278846² ≈ 0.07776
///   - delta=+1: mass[+1] × dp[1][2] = 0.115385 × 0.125 ≈ 0.01442
///   - delta=−1: mass[−1] × dp[1][1] = 0.125 × 0.115385 ≈ 0.01442
///   - total = 0.10660
///
/// P(step_0 delta=+1 | count=2) ≈ 0.01442 / 0.10660 ≈ 0.1353.
/// Without continuation weighting, the same probability would be
/// mass[+1] / (mass[0]+mass[+1]+mass[−1]) ≈ 0.2222 — a
/// distinguishing gap. The test asserts 0.1353 ± 5σ across N
/// seeds, so a regression that breaks continuation weighting
/// (e.g. dropping the DP factor) would land at ~0.22 and fail.
#[test]
fn indel_pass_count_2_matches_exact_continuation_weighting() {
    use crate::contract::{IndelEventClass, IndelKindHint};

    let (cfg, sim) = balanced_dp_distribution_fixture();
    let contracts = productive();

    const N: u32 = 4000;
    let mut delta_plus = 0u32;
    let mut delta_zero = 0u32;
    let mut delta_minus = 0u32;
    for seed in 0..N as u64 {
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(2),
            0.5,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute_with_context(
            &plan,
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        );

        // Extract step-0 (kind, site) from the trace. Classify
        // delta against the *initial* sim — that's the working
        // sim seen by step 0.
        let kind_0 = match outcome.trace.find("corrupt.indel.kind[0]").unwrap().value {
            ChoiceValue::Bool(b) => b,
            _ => panic!("unexpected kind variant"),
        };
        let site_0 = match outcome.trace.find("corrupt.indel.site[0]").unwrap().value {
            ChoiceValue::Int(s) => s as u32,
            _ => panic!("unexpected site variant"),
        };
        // Per the IndelEvent trace convention: kind=true is
        // insertion, kind=false is deletion. NoOp records site=−1,
        // which means permissive empty-support — skip it.
        if site_0 as i64 == -1 {
            continue;
        }
        let kind_hint = if kind_0 {
            IndelKindHint::Insertion
        } else {
            IndelKindHint::Deletion
        };
        let class = contracts.admissible_indel_class_at(&sim, Some(&cfg), site_0, kind_hint);
        match class {
            IndelEventClass::FrameNeutral => delta_zero += 1,
            IndelEventClass::FrameDelta(1) => delta_plus += 1,
            IndelEventClass::FrameDelta(-1) => delta_minus += 1,
            other => panic!(
                "unexpected delta class {:?} for sampled step-0 event",
                other
            ),
        }
    }
    let total = delta_plus + delta_zero + delta_minus;
    let p_plus = delta_plus as f64 / total as f64;
    let p_minus = delta_minus as f64 / total as f64;
    let p_zero = delta_zero as f64 / total as f64;

    // Theoretical step-0 marginals with continuation weighting:
    //   P(δ=+1) = 0.135338, P(δ=−1) = 0.135338, P(δ=0) = 0.729323.
    let expected_plus = 0.135338;
    let expected_minus = 0.135338;
    let expected_zero = 0.729323;
    let std_plus = ((total as f64) * expected_plus * (1.0 - expected_plus)).sqrt() / (total as f64);
    let std_zero = ((total as f64) * expected_zero * (1.0 - expected_zero)).sqrt() / (total as f64);
    let tol_plus = 5.0 * std_plus;
    let tol_zero = 5.0 * std_zero;

    assert!(
        (p_plus - expected_plus).abs() < tol_plus,
        "p(δ=+1) = {:.4}, expected {:.4} ± {:.4} (5σ); counts: + {} 0 {} − {}",
        p_plus,
        expected_plus,
        tol_plus,
        delta_plus,
        delta_zero,
        delta_minus
    );
    assert!(
        (p_minus - expected_minus).abs() < tol_plus,
        "p(δ=−1) = {:.4}, expected {:.4} ± {:.4} (5σ)",
        p_minus,
        expected_minus,
        tol_plus
    );
    assert!(
        (p_zero - expected_zero).abs() < tol_zero,
        "p(δ=0) = {:.4}, expected {:.4} ± {:.4} (5σ)",
        p_zero,
        expected_zero,
        tol_zero
    );
}

#[test]
fn indel_pass_strict_errors_when_deletion_filter_empty() {
    let (cfg, sim) = make_substitution_productive_vj_fixture();
    let contracts = productive();

    let mut plan = PassPlan::new();
    plan.push(Box::new(IndelPass::new(
        fixed_count(1),
        0.0,
        Box::new(UniformBase),
    )));

    let err = PassRuntime::execute_strict_with_context(&plan, sim, 0, Some(&cfg), Some(&contracts))
        .unwrap_err();

    assert_eq!(err.pass_name(), "corrupt.indel");
    assert_eq!(err.address(), "corrupt.indel.site[0]");
    assert_eq!(
        err.constraint_reason(),
        Some(FilteredSampleError::EmptyAdmissibleSupport)
    );
}
