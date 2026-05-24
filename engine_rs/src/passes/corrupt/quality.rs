//! `QualityErrorPass` — sequencing-error model (E.5).

use crate::address;
use crate::dist::Distribution;
use crate::ir::{NucHandle, Simulation};
use crate::pass::{Pass, PassContext, PassEffect, PassError};
use crate::passes::mutation_transaction::MutationTransaction;
use crate::trace::ChoiceValue;

fn lowercase_base(base: u8) -> u8 {
    base.to_ascii_lowercase()
}

/// Models sequencing errors as random base substitutions written
/// in **lowercase** to mark the position as corrupted.
///
/// **Biological convention:** uppercase bases are germline-derived;
/// lowercase bases were mutated or corrupted at some point (SHM,
/// NP-derived, sequencing error, etc.). This pass writes lowercase
/// substitutions so downstream queries can distinguish "this
/// position was hit by a sequencing error" from "this position was
/// germline."
///
/// Mechanically identical to `PCRErrorPass` and
/// `UniformMutationPass` (count + per-error site + per-error
/// base + `with_base_changed`), with two distinctions:
/// - The destination base is lowercased before being written.
/// - Trace addresses use the `corrupt.quality.*` prefix.
/// When contracts are active, filtering evaluates the lowercased
/// value that will actually be written, while still sampling from the
/// configured base distribution.
///
/// Pass `name()` is `"corrupt.quality"`. Codon-rail consistency
/// is preserved by `with_base_changed` (the codon translator is
/// case-insensitive, so lowercase bases translate to the same
/// amino acid as their uppercase counterparts).
///
/// Trace addresses (D3):
/// - `corrupt.quality.count` — total errors applied
/// - `corrupt.quality.error_site[i]` — pool position of the i-th error
/// - `corrupt.quality.error_base[i]` — *lowercase* destination base
pub struct QualityErrorPass {
    count_dist: Box<dyn Distribution<Output = i64>>,
    base_dist: Box<dyn Distribution<Output = u8>>,
}

impl QualityErrorPass {
    pub fn new(
        count_dist: Box<dyn Distribution<Output = i64>>,
        base_dist: Box<dyn Distribution<Output = u8>>,
    ) -> Self {
        Self {
            count_dist,
            base_dist,
        }
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        let count_raw = self.count_dist.sample(ctx.rng);
        if strict && count_raw < 0 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::CORRUPT_QUALITY_COUNT,
                count_raw,
                "negative_count",
            ));
        }
        if strict && count_raw > u32::MAX as i64 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::CORRUPT_QUALITY_COUNT,
                count_raw,
                "count_exceeds_u32",
            ));
        }
        assert!(
            count_raw >= 0,
            "QualityErrorPass: count distribution returned negative {}",
            count_raw
        );
        assert!(
            count_raw <= u32::MAX as i64,
            "QualityErrorPass: count distribution returned {} > u32::MAX",
            count_raw
        );
        let count = count_raw as u32;
        ctx.trace
            .record(address::CORRUPT_QUALITY_COUNT, ChoiceValue::Int(count_raw));

        let pool_len = sim.pool.len() as u32;
        if pool_len == 0 || count == 0 {
            return Ok(sim.clone());
        }

        // Quality substitutions go through MutationTransaction with
        // a value-transform that lowercases the chosen base — the
        // sequencing-error convention. The TX handles observer
        // attach, contract filtering, trace recording, and the seal
        // protocol; this loop only samples sites and submits.
        let mut tx = MutationTransaction::open(sim, ctx, self.name(), strict);

        for i in 0..count {
            let site = tx.rng().range_u32(pool_len);
            let site_handle = NucHandle::new(site);
            tx.trace().record(
                address::corrupt_quality_site(i),
                ChoiceValue::Int(site as i64),
            );
            tx.substitute_base(
                site_handle,
                self.base_dist.as_ref(),
                &address::corrupt_quality_base(i),
                i,
                count,
                Some(lowercase_base),
            )?;
        }

        tx.commit()
    }
}

impl Pass for QualityErrorPass {
    fn name(&self) -> &str {
        address::CORRUPT_QUALITY
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx, false)
            .expect("QualityErrorPass permissive execution must not return PassError")
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
            address::CORRUPT_QUALITY_COUNT.to_string(),
            address::CORRUPT_QUALITY_SITE_PATTERN.to_string(),
            address::CORRUPT_QUALITY_BASE_PATTERN.to_string(),
        ]
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::EditBases]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::contract::productive;
    use crate::dist::{EmpiricalLengthDist, FilteredSampleError, UniformBase};
    use crate::ir::{compute_codon_rail, Nucleotide, Region, Segment};
    use crate::pass::PassPlan;
    use crate::pass::testing::PassRuntime;
    use crate::passes::test_support::{
        make_substitution_productive_vj_fixture, StopOnlyMutationBaseDist,
        StopThenSafeMutationBaseDist,
    };

    fn quality_test_sim() -> Simulation {
        // All-uppercase germline so the lowercase-after-error
        // distinction is visible.
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGGTTT".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        sim
    }

    fn quality_single_error_plan(base_dist: Box<dyn Distribution<Output = u8>>) -> PassPlan {
        let mut plan = PassPlan::new();
        plan.push(Box::new(QualityErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
            base_dist,
        )));
        plan
    }

    fn quality_error_site(seed: u64, sim: Simulation) -> u32 {
        let outcome = PassRuntime::execute(
            &quality_single_error_plan(Box::new(StopOnlyMutationBaseDist)),
            sim,
            seed,
        );
        match outcome
            .trace
            .find("corrupt.quality.error_site[0]")
            .unwrap()
            .value
        {
            ChoiceValue::Int(site) => site as u32,
            _ => panic!("wrong variant"),
        }
    }

    fn find_seed_for_quality_error_site(sim: &Simulation, target_site: u32) -> u64 {
        for seed in 0..512u64 {
            if quality_error_site(seed, sim.clone()) == target_site {
                return seed;
            }
        }
        panic!(
            "no seed in search range targeted quality site {}",
            target_site
        );
    }

    #[test]
    fn quality_error_pass_writes_lowercase_bases() {
        // The biological convention: every position hit by a
        // quality error should be lowercase in the post-pass pool.
        let mut plan = PassPlan::new();
        plan.push(Box::new(QualityErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, quality_test_sim(), 13);
        let final_sim = outcome.final_simulation();

        // Collect the (site, recorded_base) pairs, deduping for
        // last-write-wins semantics.
        let mut last_at_site: std::collections::HashMap<u32, u8> = std::collections::HashMap::new();
        for i in 0..5 {
            let s = match outcome
                .trace
                .find(&format!("corrupt.quality.error_site[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Int(n) => n as u32,
                _ => unreachable!(),
            };
            let b = match outcome
                .trace
                .find(&format!("corrupt.quality.error_base[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Base(b) => b,
                _ => unreachable!(),
            };
            // The TRACE value should be lowercase too — faithfulness.
            assert!(
                b.is_ascii_lowercase(),
                "trace base at error_base[{}] = {} is not lowercase",
                i,
                b as char
            );
            last_at_site.insert(s, b);
        }
        // Pool reflects the recorded lowercase bases.
        for (&site, &expected) in last_at_site.iter() {
            let actual = final_sim.pool.get(NucHandle::new(site)).unwrap().base;
            assert_eq!(actual, expected);
            assert!(actual.is_ascii_lowercase());
        }
    }

    #[test]
    fn quality_error_pass_codon_rail_unaffected_by_case() {
        // The codon translator is case-insensitive, so even though
        // bases are now lowercase the amino_acids should match
        // what an all-uppercase recompute would produce.
        let mut sim = Simulation::new();
        for (i, b) in b"ATGGGGGGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9));
        sim = sim.with_region_added(region);
        // Original codon rail: ATG GGG GGG → M G G.
        assert_eq!(compute_codon_rail(&sim.sequence.regions[0], &sim.pool).amino_acids, b"MGG");

        let mut plan = PassPlan::new();
        plan.push(Box::new(QualityErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, sim, 0);
        let final_sim = outcome.final_simulation();

        // After the pass, the codon rail still translates correctly
        // even though some bases are lowercase. Computed on demand.
        let fresh = compute_codon_rail(&final_sim.sequence.regions[0], &final_sim.pool);
        // No 'X' (ambiguous) amino acids — lowercase still translates.
        let stored = &fresh.amino_acids;
        // No 'X' (ambiguous) amino acids — lowercase still translates.
        for &aa in stored {
            assert_ne!(
                aa, b'X',
                "lowercase bases should still translate cleanly, got X in {:?}",
                stored
            );
        }
    }

    #[test]
    fn quality_error_pass_zero_count_is_noop() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(QualityErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, quality_test_sim(), 0);
        assert_eq!(outcome.trace.len(), 1);
        // No lowercase bases anywhere in the pool.
        for i in 0..outcome.final_simulation().pool.len() {
            let b = outcome
                .final_simulation()
                .pool
                .get(NucHandle::new(i as u32))
                .unwrap()
                .base;
            assert!(b.is_ascii_uppercase());
        }
    }

    #[test]
    fn quality_error_pass_is_deterministic() {
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(QualityErrorPass::new(
                Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
                Box::new(UniformBase),
            )));
            p
        };
        let oa = PassRuntime::execute(&plan(), quality_test_sim(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&plan(), quality_test_sim(), 0xc0ff_ee);
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
    fn quality_error_pass_declared_choices() {
        let pass = QualityErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
            Box::new(UniformBase),
        );
        let declared = pass.declared_choices();
        assert!(declared.contains(&"corrupt.quality.count".to_string()));
        assert!(declared.contains(&"corrupt.quality.error_site[0..n]".to_string()));
        assert!(declared.contains(&"corrupt.quality.error_base[0..n]".to_string()));
    }

    #[test]
    fn quality_error_productive_filters_lowercase_base_that_would_create_stop() {
        let (cfg, sim) = make_substitution_productive_vj_fixture();
        let seed = find_seed_for_quality_error_site(&sim, 2);
        let contracts = productive();

        let constrained = PassRuntime::execute_with_context(
            &quality_single_error_plan(Box::new(StopThenSafeMutationBaseDist)),
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        );

        assert_eq!(
            constrained
                .trace
                .find("corrupt.quality.error_site[0]")
                .unwrap()
                .value,
            ChoiceValue::Int(2)
        );
        assert_eq!(
            constrained
                .trace
                .find("corrupt.quality.error_base[0]")
                .unwrap()
                .value,
            ChoiceValue::Base(b'c')
        );
        assert!(contracts
            .verify(constrained.final_simulation(), Some(&cfg))
            .is_ok());

        let unconstrained = PassRuntime::execute_with_context(
            &quality_single_error_plan(Box::new(StopThenSafeMutationBaseDist)),
            sim,
            seed,
            Some(&cfg),
            None,
        );
        assert_eq!(
            unconstrained
                .trace
                .find("corrupt.quality.error_base[0]")
                .unwrap()
                .value,
            ChoiceValue::Base(b'a')
        );
        let violations = productive()
            .verify(unconstrained.final_simulation(), Some(&cfg))
            .unwrap_err();
        assert!(violations
            .iter()
            .any(|v| v.contract_name == "no_stop_codon_in_junction"));
    }

    #[test]
    fn quality_error_strict_errors_when_base_filter_empty() {
        let (cfg, sim) = make_substitution_productive_vj_fixture();
        let seed = find_seed_for_quality_error_site(&sim, 2);
        let contracts = productive();

        let err = PassRuntime::execute_strict_with_context(
            &quality_single_error_plan(Box::new(StopOnlyMutationBaseDist)),
            sim,
            seed,
            Some(&cfg),
            Some(&contracts),
        )
        .unwrap_err();

        assert_eq!(err.pass_name(), "corrupt.quality");
        assert_eq!(err.address(), "corrupt.quality.error_base[0]");
        assert_eq!(
            err.constraint_reason(),
            Some(FilteredSampleError::EmptyAdmissibleSupport)
        );
    }
}
