//! `S5FMutationPass` — context-dependent SHM model (Phase E.3).

use crate::contract::ChoiceContext;
use crate::dist::{Distribution, FilteredSampleError};
use crate::ir::{NucHandle, NucleotidePool, Simulation};
use crate::pass::{Pass, PassContext, PassEffect, PassError};
use crate::s5f::S5FKernel;
use crate::trace::ChoiceValue;

/// Context-dependent SHM mutation pass using the S5F kernel
/// (Yaari et al. 2013).
///
/// **Algorithm** (per Yaari et al., per-mutation iterative form):
///
/// 1. Sample mutation count N from `count_dist`.
/// 2. For each of N iterations:
///    a. Walk the *current* pool, building a list of
///       `(position, mutability)` pairs for every position whose
///       5-mer context (positions `[pos-2, pos+2]`) is fully
///       defined in A/C/G/T and has non-zero kernel mutability.
///    b. Sample a position from the list, weighted by mutability.
///    c. Build the 5-mer context at the chosen position; look up
///       `kernel.substitution_row(context)`; sample a destination
///       base weighted by those four probabilities.
///    d. Apply the mutation via `sim.with_base_changed`. The
///       affected region's codon rail auto-refreshes (post-D
///       audit fix).
///
/// **Why iterative recomputation:** mutating a base changes the
/// 5-mer contexts of its neighbors (positions `[pos-2, pos+2]`),
/// which changes their mutabilities and substitution distributions.
/// The per-iteration recompute reflects the actual state. Cost is
/// O(N × pool_len) for N mutations on a pool of length pool_len —
/// for typical SHM (N≈30, pool_len≈400) that's 12,000 ops, well
/// within the per-simulation budget.
///
/// **Edge cases:**
/// - Pool shorter than 5 bases: no valid 5-mer contexts; pass is a
///   no-op (count is recorded, no mutations emitted).
/// - All contexts have mutability 0 in this pool: pass stops
///   early at the first iteration that finds an empty profile.
/// - Substitution row sums to 0 for the chosen context: skip this
///   iteration (unmutable destination).
///
/// **Trace addresses (D3):**
/// - `mutate.s5f.count` — sampled mutation count
/// - `mutate.s5f.site[i]` — pool position of the i-th mutation
/// - `mutate.s5f.base[i]` — destination base of the i-th mutation
pub struct S5FMutationPass {
    kernel: S5FKernel,
    count_dist: Box<dyn Distribution<Output = i64>>,
}

impl S5FMutationPass {
    pub fn new(kernel: S5FKernel, count_dist: Box<dyn Distribution<Output = i64>>) -> Self {
        Self { kernel, count_dist }
    }

    /// Build a 5-mer at `pos` from the current pool, encoded as a
    /// context index. Returns `None` if `pos` is too close to the
    /// pool boundary or if any base in the 5-mer is non-A/C/G/T.
    fn context_at(pool: &NucleotidePool, pos: u32) -> Option<u16> {
        if pos < 2 || pos + 2 >= pool.len() as u32 {
            return None;
        }
        let b1 = pool.get(NucHandle::new(pos - 2))?.base;
        let b2 = pool.get(NucHandle::new(pos - 1))?.base;
        let b3 = pool.get(NucHandle::new(pos))?.base;
        let b4 = pool.get(NucHandle::new(pos + 1))?.base;
        let b5 = pool.get(NucHandle::new(pos + 2))?.base;
        S5FKernel::encode_context(b1, b2, b3, b4, b5)
    }

    /// Walk the pool and build the per-position mutability profile —
    /// (position, mutability) pairs for every position with a valid
    /// non-zero-mutability context.
    fn build_profile(&self, pool: &NucleotidePool) -> Vec<(u32, f64)> {
        let n = pool.len() as u32;
        if n < 5 {
            return Vec::new();
        }
        let mut profile = Vec::with_capacity((n - 4) as usize);
        for pos in 2..n - 2 {
            if let Some(ctx) = Self::context_at(pool, pos) {
                let mu = self.kernel.mutability(ctx);
                if mu > 0.0 {
                    profile.push((pos, mu));
                }
            }
        }
        profile
    }

    fn positive_row_candidates(row: [f64; 4]) -> Vec<(u8, f64)> {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        BASES
            .iter()
            .copied()
            .zip(row)
            .filter(|(_, weight)| *weight > 0.0)
            .collect()
    }

    fn sample_weighted_base(
        rng: &mut crate::rng::Rng,
        candidates: &[(u8, f64)],
    ) -> Result<Option<u8>, FilteredSampleError> {
        if candidates.is_empty() {
            return Ok(None);
        }

        let total: f64 = candidates.iter().map(|(_, weight)| weight).sum();
        if !total.is_finite() || total <= 0.0 {
            return Err(FilteredSampleError::InvalidFilteredSupport);
        }

        let r = rng.next_f64() * total;
        let mut cum = 0.0;
        for &(base, weight) in candidates {
            cum += weight;
            if r < cum {
                return Ok(Some(base));
            }
        }

        Ok(Some(candidates.last().expect("non-empty candidates").0))
    }

    fn constraint_sampling_error(
        &self,
        address: impl Into<String>,
        reason: FilteredSampleError,
    ) -> PassError {
        PassError::constraint_sampling(self.name(), address, reason)
    }

    fn sample_base(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        address: &str,
        index: u32,
        count: u32,
        site: NucHandle,
        row: [f64; 4],
        strict: bool,
    ) -> Result<Option<u8>, PassError> {
        let candidates = Self::positive_row_candidates(row);
        if candidates.is_empty() {
            return Ok(None);
        }

        let refdata = ctx.refdata;
        let contracts = ctx.contracts;

        if let Some(contracts) = contracts {
            let context = ChoiceContext::indexed_target(index, count, site);
            let filtered: Vec<(u8, f64)> = candidates
                .iter()
                .copied()
                .filter(|(candidate, _)| {
                    contracts
                        .admits_with_context(
                            sim,
                            refdata,
                            address,
                            &ChoiceValue::Base(*candidate),
                            context,
                        )
                        .is_ok()
                })
                .collect();

            if filtered.is_empty() {
                if strict {
                    return Err(self.constraint_sampling_error(
                        address,
                        FilteredSampleError::EmptyAdmissibleSupport,
                    ));
                }
                return Self::sample_weighted_base(ctx.rng, &candidates)
                    .map_err(|reason| self.constraint_sampling_error(address, reason));
            }

            return Self::sample_weighted_base(ctx.rng, &filtered)
                .map_err(|reason| self.constraint_sampling_error(address, reason));
        }

        Self::sample_weighted_base(ctx.rng, &candidates)
            .map_err(|reason| self.constraint_sampling_error(address, reason))
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        // 1. Sample mutation count.
        let count_raw = self.count_dist.sample(ctx.rng);
        if strict && count_raw < 0 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                "mutate.s5f.count",
                count_raw,
                "negative_count",
            ));
        }
        if strict && count_raw > u32::MAX as i64 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                "mutate.s5f.count",
                count_raw,
                "count_exceeds_u32",
            ));
        }
        assert!(
            count_raw >= 0,
            "S5FMutationPass: count distribution returned negative {}",
            count_raw
        );
        assert!(
            count_raw <= u32::MAX as i64,
            "S5FMutationPass: count distribution returned {} > u32::MAX",
            count_raw
        );
        ctx.trace
            .record("mutate.s5f.count", ChoiceValue::Int(count_raw));

        let count = count_raw as u32;
        if count == 0 || sim.pool.len() < 5 {
            return Ok(sim.clone());
        }

        let mut current = sim.clone();

        // 2. Iteratively pick (position, base) and apply.
        for i in 0..count {
            let profile = self.build_profile(&current.pool);
            if profile.is_empty() {
                // No mutable positions in this state — stop early.
                break;
            }

            // Sample position weighted by mutability.
            let total: f64 = profile.iter().map(|(_, m)| m).sum();
            if total <= 0.0 || !total.is_finite() {
                break;
            }
            let r = ctx.rng.next_f64() * total;
            let mut cum = 0.0;
            let mut chosen_pos = profile[0].0;
            for &(pos, mu) in &profile {
                cum += mu;
                if r < cum {
                    chosen_pos = pos;
                    break;
                }
            }

            // Look up substitution distribution at chosen position's
            // current context. Skip this iteration if the row sums to 0.
            let context = match Self::context_at(&current.pool, chosen_pos) {
                Some(c) => c,
                None => continue,
            };
            let row = self.kernel.substitution_row(context);
            let base_address = format!("mutate.s5f.base[{}]", i);
            let chosen_base = match self.sample_base(
                &current,
                ctx,
                &base_address,
                i,
                count,
                NucHandle::new(chosen_pos),
                row,
                strict,
            )? {
                Some(base) => base,
                None => continue,
            };

            ctx.trace.record(
                format!("mutate.s5f.site[{}]", i),
                ChoiceValue::Int(chosen_pos as i64),
            );
            ctx.trace
                .record(base_address, ChoiceValue::Base(chosen_base));

            current = current.with_base_changed(NucHandle::new(chosen_pos), chosen_base);
        }

        Ok(current)
    }
}

impl Pass for S5FMutationPass {
    fn name(&self) -> &str {
        "mutate.s5f"
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx, false)
            .expect("S5FMutationPass permissive execution must not return PassError")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.execute_with_sampling_mode(sim, ctx, true)
    }

    fn declared_choices(&self) -> Vec<String> {
        vec![
            "mutate.s5f.count".to_string(),
            "mutate.s5f.site[0..n]".to_string(),
            "mutate.s5f.base[0..n]".to_string(),
        ]
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::EditBases]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assignment::AlleleInstance;
    use crate::contract::productive;
    use crate::dist::EmpiricalLengthDist;
    use crate::ir::{Nucleotide, Region, Segment};
    use crate::pass::{PassPlan, PassRuntime};
    use crate::refdata::{Allele, AlleleId, ChainType, RefDataConfig};
    use crate::s5f::{S5F_NUM_CONTEXTS, S5F_SUBSTITUTION_LEN};

    /// Helper: build a kernel where every context has uniform
    /// mutability and uniform substitution. Useful as a stress-test
    /// kernel where every position is mutable and any of A/C/G/T
    /// can be the destination.
    fn s5f_uniform_kernel() -> S5FKernel {
        S5FKernel::new(
            vec![1.0; S5F_NUM_CONTEXTS],
            vec![0.25; S5F_SUBSTITUTION_LEN],
        )
    }

    /// Helper: build a kernel where ALL mutabilities are 0 — no
    /// position is mutable. Used to verify the early-stop path.
    fn s5f_zero_kernel() -> S5FKernel {
        S5FKernel::new(
            vec![0.0; S5F_NUM_CONTEXTS],
            vec![0.25; S5F_SUBSTITUTION_LEN],
        )
    }

    fn s5f_stop_filter_kernel(include_safe_base: bool) -> S5FKernel {
        // In the fixture below, only pool position 2 has context
        // TACAA. Mutating that central C to A creates TAA, while C
        // remains productive and proves contract filtering can pick
        // a safe S5F destination from the same substitution row.
        let context =
            S5FKernel::encode_context(b'T', b'A', b'C', b'A', b'A').expect("canonical context");
        let mut mutability = vec![0.0; S5F_NUM_CONTEXTS];
        mutability[context as usize] = 1.0;

        let mut substitution = vec![0.0; S5F_SUBSTITUTION_LEN];
        let offset = context as usize * 4;
        substitution[offset] = 1.0; // A: would create TAA.
        if include_safe_base {
            substitution[offset + 1] = 1.0; // C: safe candidate.
        }

        S5FKernel::new(mutability, substitution)
    }

    fn s5f_single_mutation_plan(kernel: S5FKernel) -> PassPlan {
        let mut plan = PassPlan::new();
        plan.push(Box::new(S5FMutationPass::new(
            kernel,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
        )));
        plan
    }

    fn s5f_productive_vj_fixture() -> (RefDataConfig, Simulation) {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(Allele {
            name: "v_s5f*01".into(),
            gene: "v_s5f".into(),
            seq: b"TACAAA".to_vec(),
            segment: Segment::V,
            anchor: Some(0),
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j_s5f*01".into(),
            gene: "j_s5f".into(),
            seq: b"TGG".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });

        let mut sim = Simulation::new();
        for (i, &b) in b"TACAAA".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
            sim = next;
        }
        let v_region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(6))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(v_region);

        for (i, &b) in b"TGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::J));
            sim = next;
        }
        let j_region = Region::new(Segment::J, NucHandle::new(6), NucHandle::new(9))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(j_region);

        sim = sim
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

        (cfg, sim)
    }

    fn find_seed_for_s5f_unconstrained_base(
        kernel: &S5FKernel,
        sim: &Simulation,
        target_base: u8,
    ) -> u64 {
        for seed in 0..512u64 {
            let outcome =
                PassRuntime::execute(&s5f_single_mutation_plan(kernel.clone()), sim.clone(), seed);
            if let Some(rec) = outcome.trace.find("mutate.s5f.base[0]") {
                if rec.value == ChoiceValue::Base(target_base) {
                    return seed;
                }
            }
        }
        panic!(
            "no seed in search range produced S5F base {}",
            target_base as char
        );
    }

    fn s5f_test_sim() -> Simulation {
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGGTTTACGTACGT".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(20))
            .with_codon_rail_recomputed(&sim.pool);
        sim.with_region_added(region)
    }

    #[test]
    fn s5f_mutation_pass_zero_count_is_noop() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(S5FMutationPass::new(
            s5f_uniform_kernel(),
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
        )));

        let sim = s5f_test_sim();
        let outcome = PassRuntime::execute(&plan, sim.clone(), 0);

        assert_eq!(outcome.final_simulation().pool.len(), sim.pool.len());
        for i in 0..sim.pool.len() {
            assert_eq!(
                outcome
                    .final_simulation()
                    .pool
                    .get(NucHandle::new(i as u32))
                    .unwrap()
                    .base,
                sim.pool.get(NucHandle::new(i as u32)).unwrap().base
            );
        }
        assert_eq!(outcome.trace.len(), 1);
        assert_eq!(
            outcome.trace.find("mutate.s5f.count").unwrap().value,
            ChoiceValue::Int(0)
        );
    }

    #[test]
    fn s5f_mutation_pass_short_pool_emits_no_mutations() {
        // Pool of 4 bases — too short for any 5-mer context.
        let mut sim = Simulation::new();
        for (i, b) in b"AAAA".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }

        let mut plan = PassPlan::new();
        plan.push(Box::new(S5FMutationPass::new(
            s5f_uniform_kernel(),
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        )));
        let outcome = PassRuntime::execute(&plan, sim, 0);

        // Count is sampled but no actual mutations applied.
        assert_eq!(
            outcome.trace.find("mutate.s5f.count").unwrap().value,
            ChoiceValue::Int(5)
        );
        // No site/base entries.
        for i in 0..5 {
            assert!(outcome
                .trace
                .find(&format!("mutate.s5f.site[{}]", i))
                .is_none());
        }
    }

    #[test]
    fn s5f_mutation_pass_zero_mutability_kernel_emits_no_mutations() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(S5FMutationPass::new(
            s5f_zero_kernel(),
            Box::new(EmpiricalLengthDist::from_pairs(vec![(10, 1.0)])),
        )));
        let outcome = PassRuntime::execute(&plan, s5f_test_sim(), 0);

        assert_eq!(
            outcome.trace.find("mutate.s5f.count").unwrap().value,
            ChoiceValue::Int(10)
        );
        // No mutations recorded — early stop on empty profile.
        for i in 0..10 {
            assert!(
                outcome
                    .trace
                    .find(&format!("mutate.s5f.site[{}]", i))
                    .is_none(),
                "expected no mutation at index {}",
                i
            );
        }
    }

    #[test]
    fn s5f_mutation_pass_applies_n_mutations_with_uniform_kernel() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(S5FMutationPass::new(
            s5f_uniform_kernel(),
            Box::new(EmpiricalLengthDist::from_pairs(vec![(7, 1.0)])),
        )));
        let outcome = PassRuntime::execute(&plan, s5f_test_sim(), 1234);

        // 1 count + 7 sites + 7 bases = 15 trace records.
        assert_eq!(outcome.trace.len(), 15);
        for i in 0..7 {
            let site = outcome
                .trace
                .find(&format!("mutate.s5f.site[{}]", i))
                .unwrap();
            let base = outcome
                .trace
                .find(&format!("mutate.s5f.base[{}]", i))
                .unwrap();
            // Site is in the valid 5-mer range [2, pool_len - 2).
            match site.value {
                ChoiceValue::Int(s) => {
                    assert!(s >= 2 && s < 18, "site {} out of range [2, 18)", s);
                }
                _ => panic!("wrong variant"),
            }
            match base.value {
                ChoiceValue::Base(b) => assert!(matches!(b, b'A' | b'C' | b'G' | b'T')),
                _ => panic!("wrong variant"),
            }
        }
    }

    #[test]
    fn s5f_mutation_pass_pool_reflects_recorded_mutations() {
        // Faithfulness: the post-mutation pool's base at each site
        // (last write wins) matches the recorded base.
        let mut plan = PassPlan::new();
        plan.push(Box::new(S5FMutationPass::new(
            s5f_uniform_kernel(),
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        )));
        let outcome = PassRuntime::execute(&plan, s5f_test_sim(), 99);
        let final_sim = outcome.final_simulation();

        let mut last_at_site: std::collections::HashMap<u32, u8> = std::collections::HashMap::new();
        for i in 0..5 {
            let s = match outcome
                .trace
                .find(&format!("mutate.s5f.site[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Int(n) => n as u32,
                _ => unreachable!(),
            };
            let b = match outcome
                .trace
                .find(&format!("mutate.s5f.base[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Base(b) => b,
                _ => unreachable!(),
            };
            last_at_site.insert(s, b);
        }
        for (&site, &expected_base) in last_at_site.iter() {
            let actual = final_sim.pool.get(NucHandle::new(site)).unwrap().base;
            assert_eq!(
                actual, expected_base,
                "trace says site {} got base {}, but pool has {}",
                site, expected_base as char, actual as char
            );
        }
    }

    #[test]
    fn s5f_mutation_pass_refreshes_codon_rail() {
        // The post-D fix: every with_base_changed call refreshes
        // the affected region's codon rail. Verify that after S5F
        // runs, the stored amino_acids match a fresh recomputation.
        let mut plan = PassPlan::new();
        plan.push(Box::new(S5FMutationPass::new(
            s5f_uniform_kernel(),
            Box::new(EmpiricalLengthDist::from_pairs(vec![(8, 1.0)])),
        )));
        let outcome = PassRuntime::execute(&plan, s5f_test_sim(), 42);
        let final_sim = outcome.final_simulation();

        let stored_aa = &final_sim.sequence.regions[0].amino_acids;
        let fresh = final_sim.sequence.regions[0].with_codon_rail_recomputed(&final_sim.pool);
        assert_eq!(stored_aa, &fresh.amino_acids);
        assert_eq!(
            final_sim.sequence.regions[0].stop_codon_positions,
            fresh.stop_codon_positions
        );
    }

    #[test]
    fn s5f_mutation_pass_is_deterministic_under_same_seed() {
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(S5FMutationPass::new(
                s5f_uniform_kernel(),
                Box::new(EmpiricalLengthDist::from_pairs(vec![(6, 1.0)])),
            )));
            p
        };

        let oa = PassRuntime::execute(&plan(), s5f_test_sim(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&plan(), s5f_test_sim(), 0xc0ff_ee);
        assert_eq!(oa.trace.choices(), ob.trace.choices());
        for i in 0..oa.final_simulation().pool.len() {
            let h = NucHandle::new(i as u32);
            assert_eq!(
                oa.final_simulation().pool.get(h).unwrap().base,
                ob.final_simulation().pool.get(h).unwrap().base
            );
        }
    }

    #[test]
    fn s5f_mutation_pass_targets_mutability_hotspot() {
        // Build a kernel where context AAAAA (index 0) has
        // mutability 1.0 and ALL OTHER contexts have mutability 0.
        // Substitution: AAAAA always mutates to T.
        // Pool: AAAAAAAAAAAAAAAA (16 A's). Every internal position's
        // 5-mer is AAAAA → all mutable. Every mutation should be A→T.
        let mut mu = vec![0.0; S5F_NUM_CONTEXTS];
        mu[0] = 1.0; // AAAAA only
        let mut sub = vec![0.0; S5F_SUBSTITUTION_LEN];
        // Context 0: dest A=0, C=0, G=0, T=1.
        sub[3] = 1.0;
        let kernel = S5FKernel::new(mu, sub);

        let mut sim = Simulation::new();
        for i in 0..16 {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b'A', i as u16, Segment::V));
            sim = next;
        }

        let mut plan = PassPlan::new();
        plan.push(Box::new(S5FMutationPass::new(
            kernel,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        )));
        let outcome = PassRuntime::execute(&plan, sim, 0);

        // Every mutation should be base T.
        for i in 0..3 {
            let base_addr = format!("mutate.s5f.base[{}]", i);
            if let Some(rec) = outcome.trace.find(&base_addr) {
                match rec.value {
                    ChoiceValue::Base(b) => assert_eq!(
                        b, b'T',
                        "mutation {} produced base {} (expected T)",
                        i, b as char
                    ),
                    _ => panic!("wrong variant"),
                }
            }
        }

        // After 3 mutations, the pool should have at most a few A→T
        // changes. Since mutating A→T at position p changes the
        // contexts at positions p±1 and p±2 (some becoming non-AAAAA),
        // some later iterations may find an empty profile and stop.
    }

    #[test]
    fn s5f_mutation_pass_declared_choices() {
        let pass = S5FMutationPass::new(
            s5f_uniform_kernel(),
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        );
        let declared = pass.declared_choices();
        assert_eq!(declared.len(), 3);
        assert!(declared.contains(&"mutate.s5f.count".to_string()));
        assert!(declared.contains(&"mutate.s5f.site[0..n]".to_string()));
        assert!(declared.contains(&"mutate.s5f.base[0..n]".to_string()));
    }

    #[test]
    fn s5f_mutation_productive_filters_base_that_would_create_stop() {
        let kernel = s5f_stop_filter_kernel(true);
        let (cfg, sim) = s5f_productive_vj_fixture();
        let seed = find_seed_for_s5f_unconstrained_base(&kernel, &sim, b'A');
        let contracts = productive();

        let constrained = PassRuntime::execute_with_context(
            &s5f_single_mutation_plan(kernel.clone()),
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        );

        assert_eq!(
            constrained.trace.find("mutate.s5f.site[0]").unwrap().value,
            ChoiceValue::Int(2)
        );
        assert_eq!(
            constrained.trace.find("mutate.s5f.base[0]").unwrap().value,
            ChoiceValue::Base(b'C')
        );
        assert!(contracts
            .verify(constrained.final_simulation(), Some(&cfg))
            .is_ok());

        let unconstrained = PassRuntime::execute_with_context(
            &s5f_single_mutation_plan(kernel),
            sim,
            seed,
            Some(&cfg),
            None,
        );
        assert_eq!(
            unconstrained
                .trace
                .find("mutate.s5f.base[0]")
                .unwrap()
                .value,
            ChoiceValue::Base(b'A')
        );
        let violations = productive()
            .verify(unconstrained.final_simulation(), Some(&cfg))
            .unwrap_err();
        assert!(violations
            .iter()
            .any(|v| v.contract_name == "no_stop_codon_in_junction"));
    }

    #[test]
    fn s5f_mutation_strict_errors_when_base_filter_empty() {
        let (cfg, sim) = s5f_productive_vj_fixture();
        let contracts = productive();

        let err = PassRuntime::execute_strict_with_context(
            &s5f_single_mutation_plan(s5f_stop_filter_kernel(false)),
            sim,
            0,
            Some(&cfg),
            Some(&contracts),
        )
        .unwrap_err();

        assert_eq!(err.pass_name(), "mutate.s5f");
        assert_eq!(err.address(), "mutate.s5f.base[0]");
        assert_eq!(
            err.constraint_reason(),
            Some(FilteredSampleError::EmptyAdmissibleSupport)
        );
    }
}
