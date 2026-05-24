//! `S5FMutationPass` — context-dependent SHM model.

use crate::address;
use crate::dist::Distribution;
use crate::ir::Simulation;
use crate::pass::{Pass, PassContext, PassEffect, PassError};
use crate::s5f::S5FKernel;

mod event;
mod execution;
mod sampling;

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
///    d. Apply the mutation via `sim.with_base_changed`. Codon-rail
///       data is not stored on `Region`; callers that need
///       translation call [`crate::ir::compute_codon_rail`] on
///       demand.
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
}

impl Pass for S5FMutationPass {
    fn name(&self) -> &str {
        address::MUTATE_S5F
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
            address::MUTATE_S5F_COUNT.to_string(),
            address::MUTATE_S5F_SITE_PATTERN.to_string(),
            address::MUTATE_S5F_BASE_PATTERN.to_string(),
        ]
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::EditBases]
    }
}

#[cfg(test)]
mod tests;
