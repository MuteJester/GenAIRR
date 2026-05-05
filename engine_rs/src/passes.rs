//! Concrete pass implementations.
//!
//! Phase B/C/E gradually populates this module with the actual
//! biology — `SampleAlleleV`, `Trim`, `GenerateNP`, `MutateS5F`,
//! and so on. Phase B introduces just two reference passes that
//! anchor the two pass shapes:
//!
//! - **Transform passes** ([`EchoPass`], B.4) — deterministic
//!   modifications of the IR. No RNG draws, no trace records.
//!   Reference implementation for any pass that simply applies a
//!   known transformation.
//!
//! - **Sampling passes** (`SampleBasePass`, B.5 — landing next) —
//!   random draws from a `Distribution` that get recorded to the
//!   trace at an addressed key. Reference implementation for every
//!   pass that consumes RNG.
//!
//! Both reference passes are fully wired through the public API so
//! external code (tests, future integration) can build and run real
//! plans without touching internal machinery.

use crate::assignment::{AlleleInstance, TrimEnd};
use crate::contract::ChoiceContext;
use crate::dist::{sample_filtered_result, Distribution, FilteredSampleError};
use crate::ir::{
    flag, NucFlags, NucHandle, Nucleotide, NucleotidePool, Region, Segment, Simulation,
};
use crate::pass::{Pass, PassContext, PassError};
use crate::refdata::AlleleId;
use crate::s5f::S5FKernel;
use crate::trace::ChoiceValue;

// ──────────────────────────────────────────────────────────────────
// EchoPass — the deterministic transform reference
// ──────────────────────────────────────────────────────────────────

/// A deterministic transform pass that appends one configured
/// nucleotide to the simulation's pool when executed.
///
/// `EchoPass` makes no random choices — it consumes no RNG words
/// and writes nothing to the trace. The IR revision it produces
/// differs from the input only by one extra nucleotide at the end
/// of the pool.
///
/// Useful for:
/// - Building synthetic test plans where the input is fully
///   determined.
/// - Anchoring the "transform pass" pattern as a code reference
///   for future deterministic passes.
/// - Smoke-testing the runtime end to end without involving the
///   RNG or the contract system.
#[derive(Clone, Debug)]
pub struct EchoPass {
    base: u8,
    germline_pos: u16,
    segment: Segment,
}

impl EchoPass {
    /// Construct an `EchoPass` that, when executed, pushes a single
    /// `Nucleotide::germline(base, germline_pos, segment)` onto the
    /// simulation pool.
    pub fn new(base: u8, germline_pos: u16, segment: Segment) -> Self {
        Self {
            base,
            germline_pos,
            segment,
        }
    }

    /// Inspect the configured base.
    pub fn base(&self) -> u8 {
        self.base
    }

    /// Inspect the configured germline position.
    pub fn germline_pos(&self) -> u16 {
        self.germline_pos
    }

    /// Inspect the configured segment.
    pub fn segment(&self) -> Segment {
        self.segment
    }
}

impl Pass for EchoPass {
    fn name(&self) -> &str {
        "echo"
    }

    fn execute(&self, sim: &Simulation, _ctx: &mut PassContext) -> Simulation {
        let (next, _h) = sim.with_nucleotide_pushed(Nucleotide::germline(
            self.base,
            self.germline_pos,
            self.segment,
        ));
        next
    }
}

// ──────────────────────────────────────────────────────────────────
// SampleBasePass — the sampling reference
// ──────────────────────────────────────────────────────────────────

/// A sampling pass that draws one base byte from a `Distribution`,
/// records the choice to the trace at a configured address, and
/// appends a synthetic nucleotide carrying that base.
///
/// `SampleBasePass` is the reference implementation for any pass
/// that consumes RNG and must record its choice for the addressed-
/// trace contract (D3 + D9). Future biology-specific passes — N-nuc
/// generation, S5F substitution, PCR error injection — all follow
/// this shape: sample → record → mutate IR.
///
/// **Address discipline:** the address string is stored on the pass
/// instance and passed by reference into `trace.record`. For passes
/// that emit multiple choices per execution (e.g., a future
/// `GenerateNPBases` pass that draws N bases), the address would be
/// constructed per-draw via `format!("...[{}]", i)`. For
/// `SampleBasePass` (one draw per execution) the address is fixed.
pub struct SampleBasePass {
    address: String,
    distribution: Box<dyn Distribution<Output = u8>>,
    segment: Segment,
    flags: NucFlags,
}

impl SampleBasePass {
    /// Construct a sampling pass that draws from `distribution` and
    /// records to `address` in the trace. The pushed nucleotide will
    /// carry the given `segment` and `flags` (typically `flag::N_NUC`
    /// for TdT-like samples, `flag::P_NUC` for P-nucleotides, or
    /// `NucFlags::empty()` for plain test-only synthetic bases).
    pub fn new(
        address: impl Into<String>,
        distribution: Box<dyn Distribution<Output = u8>>,
        segment: Segment,
        flags: NucFlags,
    ) -> Self {
        Self {
            address: address.into(),
            distribution,
            segment,
            flags,
        }
    }

    /// The configured trace address.
    pub fn address(&self) -> &str {
        &self.address
    }

    /// The segment the pushed nucleotide will be tagged with.
    pub fn segment(&self) -> Segment {
        self.segment
    }
}

impl Pass for SampleBasePass {
    fn name(&self) -> &str {
        "sample_base"
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        let base = self.distribution.sample(ctx.rng);
        ctx.trace.record(&self.address[..], ChoiceValue::Base(base));
        let (next, _h) =
            sim.with_nucleotide_pushed(Nucleotide::synthetic(base, self.segment, self.flags));
        next
    }

    fn declared_choices(&self) -> Vec<String> {
        // SampleBasePass makes exactly one draw per execution at its
        // configured address. Phase D's upstream-bound propagation
        // and build-time validator both consume this list.
        vec![self.address.clone()]
    }
}

// ──────────────────────────────────────────────────────────────────
// SampleAllelePass — recombination-stage allele sampling (C.5)
// ──────────────────────────────────────────────────────────────────

/// Sample one allele from a pool distribution and assign it to the
/// simulation's V / D / J slot.
///
/// The pass is parameterized by:
/// - **Segment** — must be `V`, `D`, or `J` (NP segments and `C`
///   are not allowed; assemble passes don't run for NP, and
///   constant-region sampling is a future concern).
/// - **Distribution** — any `Box<dyn Distribution<Output = AlleleId>>`,
///   most commonly an `AllelePoolDist` (C.3) constructed against
///   the segment's pool in the active `RefDataConfig`. Construction
///   discipline guarantees every sampled `AlleleId` is in-bounds
///   for that pool.
///
/// On execute the pass:
/// 1. Draws one `AlleleId` from the distribution via `ctx.rng`.
/// 2. Records `ChoiceValue::AlleleId(id.index())` to the trace at
///    address `"sample_allele.{segment}"`.
/// 3. Constructs a fresh `AlleleInstance` (zero trims) and assigns
///    it to the simulation's slot for `segment`.
///
/// The pass does *not* read allele bases — that work belongs to
/// the assembly pass (C.8). Sampling here only chooses the id.
pub struct SampleAllelePass {
    segment: Segment,
    distribution: Box<dyn Distribution<Output = AlleleId>>,
}

impl SampleAllelePass {
    /// Construct a sampling pass for the given segment.
    ///
    /// Panics if `segment` is anything other than V, D, or J.
    /// The constructor catches the misuse at plan-build time so
    /// the pass-execute path is panic-free for valid plans.
    pub fn new(segment: Segment, distribution: Box<dyn Distribution<Output = AlleleId>>) -> Self {
        match segment {
            Segment::V | Segment::D | Segment::J => {}
            _ => panic!(
                "SampleAllelePass: segment must be V, D, or J — got {:?}",
                segment
            ),
        }
        Self {
            segment,
            distribution,
        }
    }

    /// Inspect the configured segment.
    pub fn segment(&self) -> Segment {
        self.segment
    }

    /// The hierarchical-string address (D3) at which this pass
    /// records its choice. Same string as `name()` since the pass
    /// makes exactly one choice per execution.
    fn address(&self) -> &'static str {
        match self.segment {
            Segment::V => "sample_allele.v",
            Segment::D => "sample_allele.d",
            Segment::J => "sample_allele.j",
            // Unreachable due to constructor validation; if this
            // fires it's a code defect, not a user error.
            _ => unreachable!("SampleAllelePass with non-V/D/J segment"),
        }
    }
}

impl Pass for SampleAllelePass {
    fn name(&self) -> &str {
        self.address()
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        let id = self.distribution.sample(ctx.rng);
        ctx.trace
            .record(self.address(), ChoiceValue::AlleleId(id.index()));
        sim.with_allele_assigned(self.segment, AlleleInstance::new(id))
    }

    fn declared_choices(&self) -> Vec<String> {
        vec![self.address().to_string()]
    }
}

// ──────────────────────────────────────────────────────────────────
// UniformMutationPass — simplest SHM model (Phase E.1)
// ──────────────────────────────────────────────────────────────────

/// The simplest mutation pass: pick `N` positions uniformly across
/// the assembled pool and substitute each with a base drawn from
/// `base_dist`.
///
/// Models position-independent point mutations. Real biology uses
/// the context-dependent S5F model (Phase E.3); this pass is the
/// architectural reference for any SHM-like pass — establishes the
/// trace-address shape, the per-mutation IR-revision flow, and the
/// integration with constraint-aware sampling.
///
/// **Determinism:** consumes RNG words deterministically — one for
/// the count, then two per mutation (site index + base). Same seed
/// → same mutations.
///
/// **Codon-rail consistency:** every `with_base_changed` call
/// auto-refreshes the affected region's codon rail (per the
/// post-D audit fix). Stop codons that get introduced will be
/// visible in `Region.amino_acids` immediately.
///
/// **Constraint awareness:** every candidate base is filtered
/// against the active `ContractSet` with the chosen target site in
/// [`ChoiceContext`]. This lets contracts such as
/// `NoStopCodonInJunction` reject substitutions that would leave the
/// junction non-productive before the mutation is committed. The
/// mutation count and site choice remain unconstrained for now; the
/// base draw is the first biologically meaningful hardening point.
///
/// Trace addresses (D3 hierarchical strings):
/// - `mutate.uniform.count` — total mutations applied
/// - `mutate.uniform.site[i]` — position of the i-th mutation
/// - `mutate.uniform.base[i]` — new base at the i-th mutation
pub struct UniformMutationPass {
    count_dist: Box<dyn Distribution<Output = i64>>,
    base_dist: Box<dyn Distribution<Output = u8>>,
}

impl UniformMutationPass {
    pub fn new(
        count_dist: Box<dyn Distribution<Output = i64>>,
        base_dist: Box<dyn Distribution<Output = u8>>,
    ) -> Self {
        Self {
            count_dist,
            base_dist,
        }
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
        strict: bool,
    ) -> Result<u8, PassError> {
        let refdata = ctx.refdata;
        let contracts = ctx.contracts;

        if let Some(contracts) = contracts {
            let context = ChoiceContext::indexed_target(index, count, site);
            match sample_filtered_result(ctx.rng, self.base_dist.as_ref(), |candidate: &u8| {
                contracts
                    .admits_with_context(
                        sim,
                        refdata,
                        address,
                        &ChoiceValue::Base(*candidate),
                        context,
                    )
                    .is_ok()
            }) {
                Ok(value) => return Ok(value),
                Err(reason) if strict => {
                    return Err(self.constraint_sampling_error(address, reason));
                }
                Err(_) => {
                    // Permissive legacy path: if filtering cannot
                    // produce a value, preserve the old unconstrained
                    // mutation behaviour.
                }
            }
        }

        Ok(self.base_dist.sample(ctx.rng))
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        // 1. Sample the number of mutations to apply.
        let count_raw = self.count_dist.sample(ctx.rng);
        assert!(
            count_raw >= 0,
            "UniformMutationPass: count distribution returned negative {}",
            count_raw
        );
        let count = count_raw as u32;
        ctx.trace
            .record("mutate.uniform.count", ChoiceValue::Int(count_raw));

        // No-op if the pool is empty — nothing to mutate.
        let pool_len = sim.pool.len() as u32;
        if pool_len == 0 || count == 0 {
            return Ok(sim.clone());
        }

        // 2. Apply `count` mutations sequentially. Each mutation
        //    samples a site (uniform in [0, pool_len)) and then draws
        //    an admissible base for that target when contracts are
        //    active. The IR evolves one mutation at a time through
        //    the persistent API.
        let mut current = sim.clone();
        for i in 0..count {
            let site = ctx.rng.range_u32(pool_len);
            let site_handle = NucHandle::new(site);
            let base_address = format!("mutate.uniform.base[{}]", i);

            ctx.trace.record(
                format!("mutate.uniform.site[{}]", i),
                ChoiceValue::Int(site as i64),
            );
            let new_base =
                self.sample_base(&current, ctx, &base_address, i, count, site_handle, strict)?;
            ctx.trace.record(base_address, ChoiceValue::Base(new_base));

            // `with_base_changed` auto-refreshes the codon rail of
            // any region containing `site` (post-D audit fix).
            current = current.with_base_changed(site_handle, new_base);
        }

        Ok(current)
    }
}

impl Pass for UniformMutationPass {
    fn name(&self) -> &str {
        "mutate.uniform"
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx, false)
            .expect("UniformMutationPass permissive execution must not return PassError")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.execute_with_sampling_mode(sim, ctx, true)
    }

    fn declared_choices(&self) -> Vec<String> {
        // Count is fixed-address; site and base are variable per
        // mutation. Use the [0..n] expansion convention from D3.
        vec![
            "mutate.uniform.count".to_string(),
            "mutate.uniform.site[0..n]".to_string(),
            "mutate.uniform.base[0..n]".to_string(),
        ]
    }
}

// ──────────────────────────────────────────────────────────────────
// S5FMutationPass — context-dependent SHM model (Phase E.3)
// ──────────────────────────────────────────────────────────────────

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
        assert!(
            count_raw >= 0,
            "S5FMutationPass: count distribution returned negative {}",
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
}

// ──────────────────────────────────────────────────────────────────
// PCRErrorPass — observation-stage PCR amplification errors (E.4)
// ──────────────────────────────────────────────────────────────────

/// Models PCR amplification errors as a small number of random
/// base substitutions across the sequence. Conceptually applied at
/// observation stage (after biology) but mechanically identical to
/// `UniformMutationPass` — a count + per-mutation (site, base) draw
/// applied through `with_base_changed`.
///
/// The biological distinction (PCR errors vs SHM) is preserved
/// through:
/// - Pass `name()` is `"corrupt.pcr"` (vs `"mutate.uniform"`)
/// - Trace addresses use `corrupt.pcr.*` prefix
/// - The count distribution typically encodes a much lower rate
///   (~10⁻⁵ per base per cycle in real PCR, but the rate is
///   user-supplied so it could be anything)
///
/// Like all substitution passes, every `with_base_changed` call
/// triggers automatic codon-rail refresh (post-D audit fix).
/// When contracts are active, replacement bases are sampled from
/// the admissible subset for the chosen target site, so observation
/// errors cannot silently break enforced contracts when a safe
/// replacement exists.
///
/// Trace addresses (D3):
/// - `corrupt.pcr.count` — total PCR errors applied
/// - `corrupt.pcr.error_site[i]` — pool position of the i-th error
/// - `corrupt.pcr.error_base[i]` — replacement base at the i-th error
pub struct PCRErrorPass {
    count_dist: Box<dyn Distribution<Output = i64>>,
    base_dist: Box<dyn Distribution<Output = u8>>,
}

impl PCRErrorPass {
    pub fn new(
        count_dist: Box<dyn Distribution<Output = i64>>,
        base_dist: Box<dyn Distribution<Output = u8>>,
    ) -> Self {
        Self {
            count_dist,
            base_dist,
        }
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
        strict: bool,
    ) -> Result<u8, PassError> {
        let refdata = ctx.refdata;
        let contracts = ctx.contracts;

        if let Some(contracts) = contracts {
            let context = ChoiceContext::indexed_target(index, count, site);
            match sample_filtered_result(ctx.rng, self.base_dist.as_ref(), |candidate: &u8| {
                contracts
                    .admits_with_context(
                        sim,
                        refdata,
                        address,
                        &ChoiceValue::Base(*candidate),
                        context,
                    )
                    .is_ok()
            }) {
                Ok(value) => return Ok(value),
                Err(reason) if strict => {
                    return Err(self.constraint_sampling_error(address, reason));
                }
                Err(_) => {
                    // Permissive legacy path: if support is missing
                    // or filtered empty, fall back to unconstrained
                    // PCR error sampling.
                }
            }
        }

        Ok(self.base_dist.sample(ctx.rng))
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        let count_raw = self.count_dist.sample(ctx.rng);
        assert!(
            count_raw >= 0,
            "PCRErrorPass: count distribution returned negative {}",
            count_raw
        );
        let count = count_raw as u32;
        ctx.trace
            .record("corrupt.pcr.count", ChoiceValue::Int(count_raw));

        let pool_len = sim.pool.len() as u32;
        if pool_len == 0 || count == 0 {
            return Ok(sim.clone());
        }

        let mut current = sim.clone();
        for i in 0..count {
            let site = ctx.rng.range_u32(pool_len);
            let site_handle = NucHandle::new(site);
            let base_address = format!("corrupt.pcr.error_base[{}]", i);

            ctx.trace.record(
                format!("corrupt.pcr.error_site[{}]", i),
                ChoiceValue::Int(site as i64),
            );
            let new_base =
                self.sample_base(&current, ctx, &base_address, i, count, site_handle, strict)?;
            ctx.trace.record(base_address, ChoiceValue::Base(new_base));

            current = current.with_base_changed(site_handle, new_base);
        }

        Ok(current)
    }
}

impl Pass for PCRErrorPass {
    fn name(&self) -> &str {
        "corrupt.pcr"
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx, false)
            .expect("PCRErrorPass permissive execution must not return PassError")
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
            "corrupt.pcr.count".to_string(),
            "corrupt.pcr.error_site[0..n]".to_string(),
            "corrupt.pcr.error_base[0..n]".to_string(),
        ]
    }
}

// ──────────────────────────────────────────────────────────────────
// QualityErrorPass — sequencing-error model (E.5)
// ──────────────────────────────────────────────────────────────────

/// Models sequencing errors as random base substitutions written
/// in **lowercase** to mark the position as corrupted.
///
/// **Biological convention preserved from V5:** uppercase bases
/// are germline-derived; lowercase bases were mutated or corrupted
/// at some point (SHM, NP-derived, sequencing error, etc.). This
/// pass writes lowercase substitutions so downstream queries can
/// distinguish "this position was hit by a sequencing error" from
/// "this position was germline."
///
/// Mechanically identical to `PCRErrorPass` and
/// `UniformMutationPass` (count + per-error site + per-error
/// base + `with_base_changed`), with two distinctions:
/// - The destination base is lowercased before being written.
/// - Trace addresses use the `corrupt.quality.*` prefix.
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
}

impl Pass for QualityErrorPass {
    fn name(&self) -> &str {
        "corrupt.quality"
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        let count_raw = self.count_dist.sample(ctx.rng);
        assert!(
            count_raw >= 0,
            "QualityErrorPass: count distribution returned negative {}",
            count_raw
        );
        let count = count_raw as u32;
        ctx.trace
            .record("corrupt.quality.count", ChoiceValue::Int(count_raw));

        let pool_len = sim.pool.len() as u32;
        if pool_len == 0 || count == 0 {
            return sim.clone();
        }

        let mut current = sim.clone();
        for i in 0..count {
            let site = ctx.rng.range_u32(pool_len);
            // Sample then lowercase. This is the biological marker
            // — the recorded trace value matches what's actually
            // written to the pool, preserving faithfulness.
            let new_base = self.base_dist.sample(ctx.rng).to_ascii_lowercase();

            ctx.trace.record(
                format!("corrupt.quality.error_site[{}]", i),
                ChoiceValue::Int(site as i64),
            );
            ctx.trace.record(
                format!("corrupt.quality.error_base[{}]", i),
                ChoiceValue::Base(new_base),
            );

            current = current.with_base_changed(NucHandle::new(site), new_base);
        }

        current
    }

    fn declared_choices(&self) -> Vec<String> {
        vec![
            "corrupt.quality.count".to_string(),
            "corrupt.quality.error_site[0..n]".to_string(),
            "corrupt.quality.error_base[0..n]".to_string(),
        ]
    }
}

// ──────────────────────────────────────────────────────────────────
// ContaminantPass — wholesale sequence replacement (E.6)
// ──────────────────────────────────────────────────────────────────

/// Models read contamination: with probability `apply_prob` the
/// entire assembled pool is overwritten with bases drawn from a
/// contaminant distribution. Used to simulate primer dimers,
/// bacterial DNA, or any non-receptor sequence that ends up in a
/// receptor-sequencing library.
///
/// **Architectural shape vs other corruption passes:**
/// - PCR / quality errors are *count-driven* — sample N positions,
///   substitute each.
/// - Contaminant is *probability-driven at the read level* — one
///   coin flip decides whether the *entire* read is wiped, then if
///   yes, every base gets a contaminant draw.
///
/// This means the trace begins with a single Boolean choice:
/// `corrupt.contaminant.applied`. When that's `Bool(true)`, the
/// trace continues with one base entry per pool position; when it's
/// `Bool(false)`, no further records are emitted (the pool is
/// returned unchanged).
///
/// Codon rail consistency is preserved per-base by
/// `with_base_changed`. After a contamination event, the affected
/// region's `amino_acids` reflects the contaminant content (not the
/// original germline) — which is the desired post-contamination
/// semantics. When contracts are active, each replacement base is
/// filtered against the current intermediate IR and the target site,
/// so a contamination event cannot transiently violate enforced
/// contracts when an admissible replacement exists.
///
/// Trace addresses (D3):
/// - `corrupt.contaminant.applied` — `Bool(true)` if contamination
///   was applied, `Bool(false)` otherwise.
/// - `corrupt.contaminant.bases[i]` — i-th replacement base, only
///   present when `applied = true`.
pub struct ContaminantPass {
    apply_prob: f64,
    base_dist: Box<dyn Distribution<Output = u8>>,
}

impl ContaminantPass {
    /// Construct a contaminant pass.
    ///
    /// Panics if `apply_prob` is not in `[0.0, 1.0]` or is non-finite.
    pub fn new(apply_prob: f64, base_dist: Box<dyn Distribution<Output = u8>>) -> Self {
        assert!(
            apply_prob.is_finite() && (0.0..=1.0).contains(&apply_prob),
            "ContaminantPass: apply_prob must be in [0.0, 1.0], got {}",
            apply_prob
        );
        Self {
            apply_prob,
            base_dist,
        }
    }

    pub fn apply_prob(&self) -> f64 {
        self.apply_prob
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
        strict: bool,
    ) -> Result<u8, PassError> {
        let refdata = ctx.refdata;
        let contracts = ctx.contracts;

        if let Some(contracts) = contracts {
            let context = ChoiceContext::indexed_target(index, count, site);
            match sample_filtered_result(ctx.rng, self.base_dist.as_ref(), |candidate: &u8| {
                contracts
                    .admits_with_context(
                        sim,
                        refdata,
                        address,
                        &ChoiceValue::Base(*candidate),
                        context,
                    )
                    .is_ok()
            }) {
                Ok(value) => return Ok(value),
                Err(reason) if strict => {
                    return Err(self.constraint_sampling_error(address, reason));
                }
                Err(_) => {
                    // Permissive legacy path: if filtering cannot
                    // produce a value, preserve unconstrained
                    // contaminant replacement.
                }
            }
        }

        Ok(self.base_dist.sample(ctx.rng))
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        // 1. Coin flip: is this read contaminated?
        let coin = ctx.rng.next_f64();
        let applied = coin < self.apply_prob;
        ctx.trace
            .record("corrupt.contaminant.applied", ChoiceValue::Bool(applied));

        if !applied {
            return Ok(sim.clone());
        }

        let pool_len = sim.pool.len() as u32;
        if pool_len == 0 {
            return Ok(sim.clone());
        }

        // 2. Replace every base in the pool with a contaminant draw.
        let mut current = sim.clone();
        for i in 0..pool_len {
            let site = NucHandle::new(i);
            let address = format!("corrupt.contaminant.bases[{}]", i);
            let new_base = self.sample_base(&current, ctx, &address, i, pool_len, site, strict)?;
            ctx.trace.record(address, ChoiceValue::Base(new_base));
            current = current.with_base_changed(site, new_base);
        }
        Ok(current)
    }
}

impl Pass for ContaminantPass {
    fn name(&self) -> &str {
        "corrupt.contaminant"
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx, false)
            .expect("ContaminantPass permissive execution must not return PassError")
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
            "corrupt.contaminant.applied".to_string(),
            "corrupt.contaminant.bases[0..n]".to_string(),
        ]
    }
}

// ──────────────────────────────────────────────────────────────────
// IndelPass — insertions + deletions at the observation stage (E.7)
// ──────────────────────────────────────────────────────────────────

/// Models PCR / sequencing indels: a small number of insertions
/// and deletions sprinkled across the assembled sequence.
///
/// **Architectural distinction from substitution passes:** indels
/// change pool length and shift downstream handle positions. The
/// underlying primitives (`Simulation::with_indel_inserted` /
/// `with_indel_deleted`) handle the position-shifting and
/// region-range adjustments. Codon rails refresh automatically
/// for any region whose range changed.
///
/// **Per-indel semantics:**
/// - Each indel is independently chosen to be an insertion or
///   deletion (50/50 by default — `insertion_prob` field gives
///   user control).
/// - The position is sampled uniformly within the current pool
///   (which means later indels may target positions that didn't
///   exist before earlier indels — that's fine, the pass walks
///   the pool's *current* state at each step).
/// - For insertions: a base is sampled from `base_dist`. Inserted
///   nucleotides are tagged with `flag::INDEL_INSERTED` and
///   carry no germline provenance (`germline_pos = NO_GERMLINE_POS`).
///
/// **Trace addresses (D3):**
/// - `corrupt.indel.count` — total indel events
/// - `corrupt.indel.kind[i]` — `Bool(true)` for insertion,
///   `Bool(false)` for deletion
/// - `corrupt.indel.site[i]` — pool position
/// - `corrupt.indel.base[i]` — only for insertions; the inserted base
///
/// **Edge cases handled:**
/// - Empty pool: deletion is a no-op for that step (insertion can
///   still happen at position 0).
/// - Deletion when pool length = 0: the kind is recorded but the
///   IR isn't modified (no-op).
pub struct IndelPass {
    count_dist: Box<dyn Distribution<Output = i64>>,
    insertion_prob: f64,
    base_dist: Box<dyn Distribution<Output = u8>>,
}

impl IndelPass {
    /// Construct an indel pass.
    ///
    /// Panics if `insertion_prob` is not in `[0.0, 1.0]` or is
    /// non-finite.
    pub fn new(
        count_dist: Box<dyn Distribution<Output = i64>>,
        insertion_prob: f64,
        base_dist: Box<dyn Distribution<Output = u8>>,
    ) -> Self {
        assert!(
            insertion_prob.is_finite() && (0.0..=1.0).contains(&insertion_prob),
            "IndelPass: insertion_prob must be in [0.0, 1.0], got {}",
            insertion_prob
        );
        Self {
            count_dist,
            insertion_prob,
            base_dist,
        }
    }

    /// Best-effort segment provenance for an inserted observation base.
    ///
    /// Insertions are synthetic, but downstream metadata still needs the
    /// nucleotide's segment to agree with the sequence context it lands in.
    /// This mirrors `Sequence::with_indel_adjusted`: if the insertion point
    /// is inside a region, or exactly at the start of a following region,
    /// the inserted nucleotide belongs to that region. Outside all regions,
    /// fall back to the nearest nucleotide's segment.
    fn insertion_segment(sim: &Simulation, at: u32) -> Segment {
        for region in &sim.sequence.regions {
            let start = region.start.index();
            let end = region.end.index();
            if start <= at && at < end {
                return region.segment;
            }
        }

        if at < sim.pool.len() as u32 {
            return sim
                .pool
                .get(NucHandle::new(at))
                .map(|n| n.segment)
                .unwrap_or(Segment::V);
        }

        if at > 0 {
            return sim
                .pool
                .get(NucHandle::new(at - 1))
                .map(|n| n.segment)
                .unwrap_or(Segment::V);
        }

        Segment::V
    }
}

impl Pass for IndelPass {
    fn name(&self) -> &str {
        "corrupt.indel"
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        let count_raw = self.count_dist.sample(ctx.rng);
        assert!(
            count_raw >= 0,
            "IndelPass: count distribution returned negative {}",
            count_raw
        );
        let count = count_raw as u32;
        ctx.trace
            .record("corrupt.indel.count", ChoiceValue::Int(count_raw));

        if count == 0 {
            return sim.clone();
        }

        let mut current = sim.clone();
        for i in 0..count {
            // Decide insertion vs deletion.
            let insertion = ctx.rng.next_f64() < self.insertion_prob;
            ctx.trace.record(
                format!("corrupt.indel.kind[{}]", i),
                ChoiceValue::Bool(insertion),
            );

            let pool_len = current.pool.len() as u32;

            if insertion {
                // Pool position can be in [0, pool_len] inclusive
                // (insertion at end is allowed).
                let site = if pool_len == 0 {
                    0
                } else {
                    ctx.rng.range_u32(pool_len + 1)
                };
                ctx.trace.record(
                    format!("corrupt.indel.site[{}]", i),
                    ChoiceValue::Int(site as i64),
                );

                let base = self.base_dist.sample(ctx.rng);
                ctx.trace.record(
                    format!("corrupt.indel.base[{}]", i),
                    ChoiceValue::Base(base),
                );

                let segment = Self::insertion_segment(&current, site);
                let new_nuc = Nucleotide::synthetic(base, segment, flag::INDEL_INSERTED);
                current = current.with_indel_inserted(site, new_nuc);
            } else {
                // Deletion: skip if pool empty.
                if pool_len == 0 {
                    // Record a sentinel position so the trace stays
                    // structurally consistent; downstream queries
                    // can detect the no-op via this value.
                    ctx.trace
                        .record(format!("corrupt.indel.site[{}]", i), ChoiceValue::Int(-1));
                    continue;
                }
                let site = ctx.rng.range_u32(pool_len);
                ctx.trace.record(
                    format!("corrupt.indel.site[{}]", i),
                    ChoiceValue::Int(site as i64),
                );
                current = current.with_indel_deleted(site);
            }
        }

        current
    }

    fn declared_choices(&self) -> Vec<String> {
        vec![
            "corrupt.indel.count".to_string(),
            "corrupt.indel.kind[0..n]".to_string(),
            "corrupt.indel.site[0..n]".to_string(),
            "corrupt.indel.base[0..n]".to_string(),
        ]
    }
}

// ──────────────────────────────────────────────────────────────────
// AssembleSegmentPass — copy a germline allele slice into the pool (C.8)
// ──────────────────────────────────────────────────────────────────

/// Assemble one germline segment (V, D, or J) from its assigned
/// `AlleleInstance` into the simulation pool.
///
/// Reads the allele's reference bases from
/// `ctx.refdata.expect(...)`, slices off `trim_5` from the start
/// and `trim_3` from the end, and pushes the retained bases as
/// germline-derived nucleotides into the pool. Then constructs a
/// `Region` for the segment with `frame_phase` chained from the
/// cumulative length of all prior regions, and recomputes its
/// codon rail.
///
/// **Determinism:** the pass makes no random choices.
/// `declared_choices()` returns the empty vector — there is
/// nothing to record to the trace.
///
/// **Pre-conditions:**
/// - An allele must be assigned to `segment` (panics otherwise —
///   propagated from `Simulation::assignments.get`).
/// - `ctx.refdata` must be `Some(...)` (panics otherwise —
///   `PassContext::refdata` is `Option` so this pass can be
///   constructed in any plan, but it requires the data at
///   execute time).
/// - `trim_5 + trim_3 <= allele.len()` (panics otherwise — caller
///   bug, asserted defensively).
///
/// **Frame phase:** computed as `(cumulative prior region length)
/// % 3`. This anchors the codon frame at the start of the *first*
/// region. Real biology anchors at the V Cys; that anchor-aware
/// offset is a future refinement (Phase D / E). For C.8 the raw
/// cumulative-length-mod-3 is correct as long as the user isn't
/// reading frame-relative properties yet.
pub struct AssembleSegmentPass {
    segment: Segment,
}

impl AssembleSegmentPass {
    /// Construct an assembly pass for the given segment.
    /// Panics if `segment` is not V, D, or J.
    pub fn new(segment: Segment) -> Self {
        match segment {
            Segment::V | Segment::D | Segment::J => {}
            _ => panic!(
                "AssembleSegmentPass: segment must be V, D, or J — got {:?}",
                segment
            ),
        }
        Self { segment }
    }

    pub fn segment(&self) -> Segment {
        self.segment
    }
}

impl Pass for AssembleSegmentPass {
    fn name(&self) -> &str {
        match self.segment {
            Segment::V => "assemble.v",
            Segment::D => "assemble.d",
            Segment::J => "assemble.j",
            _ => unreachable!("AssembleSegmentPass with non-V/D/J segment"),
        }
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        let refdata = ctx.refdata.unwrap_or_else(|| {
            panic!(
                "AssembleSegmentPass({:?}): PassContext.refdata is None — \
                 use PassRuntime::execute_with_refdata for plans containing \
                 assembly passes",
                self.segment
            )
        });

        let inst = sim
            .assignments
            .get(self.segment)
            .copied()
            .unwrap_or_else(|| {
                panic!(
                    "AssembleSegmentPass({:?}): no allele assigned — \
                 SampleAllelePass for this segment must run before assembly",
                    self.segment
                )
            });

        let allele = refdata
            .get(self.segment, inst.allele_id)
            .unwrap_or_else(|| {
                panic!(
                    "AssembleSegmentPass({:?}): allele_id {:?} out of bounds \
                 for refdata pool",
                    self.segment, inst.allele_id
                )
            });

        let trim_5 = inst.trim_5 as u32;
        let trim_3 = inst.trim_3 as u32;
        let allele_len = allele.len();

        assert!(
            trim_5 + trim_3 <= allele_len,
            "AssembleSegmentPass({:?}): trim_5 ({}) + trim_3 ({}) exceeds \
             allele length ({}) for {}",
            self.segment,
            trim_5,
            trim_3,
            allele_len,
            allele.name
        );

        let slice_start = trim_5;
        let slice_end = allele_len - trim_3;
        let slice_len = slice_end - slice_start;

        // Frame phase: cumulative length of prior regions, mod 3.
        // The frame is implicitly anchored at the start of the first
        // assembled region. Anchor-aware framing arrives in Phase D/E.
        let cumulative_len: u32 = sim.sequence.regions.iter().map(|r| r.len()).sum();
        let frame_phase = (cumulative_len % 3) as u8;

        let region_start = NucHandle::new(sim.pool.len() as u32);

        // Push the post-trim slice into the pool as germline nucleotides.
        // Each base carries its position-in-original-allele as `germline_pos`.
        let mut current = sim.clone();
        for i in 0..slice_len {
            let allele_pos = slice_start + i;
            let base = allele.seq[allele_pos as usize];
            let (next, _h) = current.with_nucleotide_pushed(Nucleotide::germline(
                base,
                allele_pos as u16,
                self.segment,
            ));
            current = next;
        }

        let region_end = NucHandle::new(current.pool.len() as u32);
        let region = Region::new(self.segment, region_start, region_end)
            .with_frame_phase(frame_phase)
            .with_codon_rail_recomputed(&current.pool);

        current.with_region_added(region)
    }

    fn declared_choices(&self) -> Vec<String> {
        // Deterministic — no choices.
        Vec::new()
    }
}

// ──────────────────────────────────────────────────────────────────
// GenerateNPPass — TdT-like N-nucleotide region generation (C.7)
// ──────────────────────────────────────────────────────────────────

/// Generate one NP (non-template) region — a stretch of bases
/// interpolated between adjacent V/D/J segments by the TdT enzyme
/// during recombination.
///
/// Single pass produces single IR revision (D2 boundary), but
/// internally records *each* atomic random choice to the trace at
/// its own indexed address. Length, then N bases. The trace stays
/// faithful at the choice level even though the IR revision count
/// stays at the biological-event level.
///
/// Parameterized by:
/// - **np_segment** — `Segment::Np1` or `Segment::Np2`. Other
///   segments rejected at construction.
/// - **length_dist** — `Box<dyn Distribution<Output = i64>>`.
///   Empirical NP length distribution. Negative or oversized
///   lengths are caller bugs and panic at execute time.
/// - **base_dist** — `Box<dyn Distribution<Output = u8>>`.
///   Per-base distribution (typically `UniformBase` until empirical
///   TdT models arrive in Phase E). When contracts are active and the
///   distribution exposes finite support, each base draw is filtered
///   through `ContractSet::admits` before being recorded.
///
/// On execute:
/// 1. Sample length L from `length_dist`. Validate `0 <= L <= u32::MAX`.
/// 2. Record `ChoiceValue::Int(L)` at `"np.{np1|np2}.length"`.
/// 3. For each of L positions: sample a base, record at
///    `"np.{np1|np2}.bases[i]"`, push a synthetic nucleotide
///    (`Nucleotide::synthetic(base, np_segment, N_NUC)`) onto the
///    pool.
/// 4. Construct `Region(np_segment, [start, start+L))` with
///    `frame_phase = 0` (cross-region frame chaining is C.8's
///    responsibility) and recompute its codon rail against the
///    new pool.
/// 5. Append the region to the sequence.
pub struct GenerateNPPass {
    np_segment: Segment,
    length_dist: Box<dyn Distribution<Output = i64>>,
    base_dist: Box<dyn Distribution<Output = u8>>,
}

impl GenerateNPPass {
    /// Construct an NP-generation pass.
    ///
    /// Panics if `np_segment` is not `Segment::Np1` or
    /// `Segment::Np2`.
    pub fn new(
        np_segment: Segment,
        length_dist: Box<dyn Distribution<Output = i64>>,
        base_dist: Box<dyn Distribution<Output = u8>>,
    ) -> Self {
        match np_segment {
            Segment::Np1 | Segment::Np2 => {}
            _ => panic!(
                "GenerateNPPass: np_segment must be Np1 or Np2 — got {:?}",
                np_segment
            ),
        }
        Self {
            np_segment,
            length_dist,
            base_dist,
        }
    }

    pub fn np_segment(&self) -> Segment {
        self.np_segment
    }

    fn length_address(&self) -> &'static str {
        match self.np_segment {
            Segment::Np1 => "np.np1.length",
            Segment::Np2 => "np.np2.length",
            _ => unreachable!("GenerateNPPass with non-NP segment"),
        }
    }

    fn bases_prefix(&self) -> &'static str {
        match self.np_segment {
            Segment::Np1 => "np.np1.bases",
            Segment::Np2 => "np.np2.bases",
            _ => unreachable!("GenerateNPPass with non-NP segment"),
        }
    }

    fn pass_name(&self) -> &'static str {
        match self.np_segment {
            Segment::Np1 => "generate_np.np1",
            Segment::Np2 => "generate_np.np2",
            _ => unreachable!("GenerateNPPass with non-NP segment"),
        }
    }

    /// Constraint-aware length sample (D.6).
    ///
    /// When `ctx.contracts` is `Some`, we use `sample_filtered_result`
    /// to draw only from values the contract set admits at the
    /// length address. In permissive mode, falls back to plain
    /// `sample` when:
    /// - no contracts are active,
    /// - the distribution doesn't expose `support()` (continuous
    ///   or too-large categorical), or
    /// - the filtered support is empty (no admissible candidate).
    ///
    /// In strict mode, the last two cases return a structured
    /// `PassError` instead of silently sampling outside the contract.
    fn sample_length(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        address: &'static str,
        strict: bool,
    ) -> Result<i64, PassError> {
        // Capture borrows of fields we'll read in the closure so
        // the borrow checker doesn't see overlapping `&mut ctx`
        // and `&ctx.contracts`.
        let refdata = ctx.refdata;
        let contracts = ctx.contracts;

        if let Some(contracts) = contracts {
            match sample_filtered_result(ctx.rng, self.length_dist.as_ref(), |&candidate: &i64| {
                contracts
                    .admits(sim, refdata, address, &ChoiceValue::Int(candidate))
                    .is_ok()
            }) {
                Ok(value) => return Ok(value),
                Err(reason) if strict => {
                    return Err(self.constraint_sampling_error(address, reason));
                }
                Err(_) => {
                    // Permissive legacy path: empty/non-enumerable filters
                    // fall through to unconstrained sampling.
                }
            }
        }
        Ok(self.length_dist.sample(ctx.rng))
    }

    /// Constraint-aware base sample (D.6).
    ///
    /// This is the base-level counterpart to `sample_length`: contracts
    /// can reject a candidate base before it is written to the IR. The
    /// immediate hardening target is `NoStopCodonInJunction`, which
    /// rejects NP bases that would complete a junction stop codon.
    fn sample_base(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        address: &str,
        index: u32,
        total_len: u32,
        strict: bool,
    ) -> Result<u8, PassError> {
        let refdata = ctx.refdata;
        let contracts = ctx.contracts;

        if let Some(contracts) = contracts {
            let context = ChoiceContext::indexed(index, total_len);
            match sample_filtered_result(ctx.rng, self.base_dist.as_ref(), |candidate: &u8| {
                contracts
                    .admits_with_context(
                        sim,
                        refdata,
                        address,
                        &ChoiceValue::Base(*candidate),
                        context,
                    )
                    .is_ok()
            }) {
                Ok(value) => return Ok(value),
                Err(reason) if strict => {
                    return Err(self.constraint_sampling_error(address, reason));
                }
                Err(_) => {}
            }
        }

        Ok(self.base_dist.sample(ctx.rng))
    }

    fn constraint_sampling_error(
        &self,
        address: impl Into<String>,
        reason: FilteredSampleError,
    ) -> PassError {
        PassError::constraint_sampling(self.pass_name(), address, reason)
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        // 1. Sample length, possibly constrained by active contracts.
        //
        // **Constraint-aware sampling (D.6):** if `ctx.contracts`
        // is set, we filter the length distribution's support to
        // values that the contract set admits at this address. In
        // strict mode, inability to draw an admissible value is a
        // structured error. In permissive mode, we preserve the
        // original Phase D.6 fallback behavior.
        let address = self.length_address();
        let length = self.sample_length(sim, ctx, address, strict)?;
        assert!(
            length >= 0,
            "GenerateNPPass({}): length distribution returned negative value {}",
            address,
            length
        );
        assert!(
            length <= u32::MAX as i64,
            "GenerateNPPass({}): length distribution returned {} > u32::MAX",
            address,
            length
        );
        ctx.trace.record(address, ChoiceValue::Int(length));

        let length = length as u32;
        let region_start = NucHandle::new(sim.pool.len() as u32);

        // 2. Sample bases and push them into the pool. Each push
        //    produces a one-step persistent IR update; the loop
        //    accumulates into `current`.
        let mut current = sim.clone();
        for i in 0..length {
            let base_address = format!("{}[{}]", self.bases_prefix(), i);
            let base = self.sample_base(&current, ctx, &base_address, i, length, strict)?;
            ctx.trace.record(base_address, ChoiceValue::Base(base));
            let (next, _h) = current.with_nucleotide_pushed(Nucleotide::synthetic(
                base,
                self.np_segment,
                flag::N_NUC,
            ));
            current = next;
        }

        // 3. Build the region and append it. Codon rail recomputes
        //    against the updated pool. Frame phase chains from
        //    cumulative prior region length, same convention as
        //    AssembleSegmentPass — keeps the codon frame coherent
        //    across V/NP/D/NP/J boundaries.
        let cumulative_len: u32 = sim.sequence.regions.iter().map(|r| r.len()).sum();
        let frame_phase = (cumulative_len % 3) as u8;

        let region_end = NucHandle::new(current.pool.len() as u32);
        let region = Region::new(self.np_segment, region_start, region_end)
            .with_frame_phase(frame_phase)
            .with_codon_rail_recomputed(&current.pool);
        Ok(current.with_region_added(region))
    }
}

impl Pass for GenerateNPPass {
    fn name(&self) -> &str {
        self.pass_name()
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx, false)
            .expect("GenerateNPPass permissive execution must not return PassError")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.execute_with_sampling_mode(sim, ctx, true)
    }

    fn declared_choices(&self) -> Vec<String> {
        // Length is fixed-address. Bases are variable-count; per
        // D3 we use the literal `[0..n]` form to indicate runtime
        // expansion.
        vec![
            self.length_address().to_string(),
            format!("{}[0..n]", self.bases_prefix()),
        ]
    }
}

// ──────────────────────────────────────────────────────────────────
// TrimPass — recombination-stage trim sampling (C.6)
// ──────────────────────────────────────────────────────────────────

/// Sample a trim amount from a distribution and apply it to the
/// assigned allele on the given segment / end.
///
/// The pass is parameterized by:
/// - **Segment** — must be V, D, or J (NP and C are rejected at
///   construction).
/// - **End** — `TrimEnd::Five` or `TrimEnd::Three`. All six
///   `(segment, end)` combinations are syntactically valid; whether
///   biology uses them is up to the plan author.
/// - **Distribution** — any `Box<dyn Distribution<Output = i64>>`,
///   typically an `EmpiricalLengthDist` constructed at plan-build
///   time with the empirical trim distribution for that segment.
///
/// On execute the pass:
/// 1. Draws one `i64` trim value from the distribution via
///    `ctx.rng`.
/// 2. Validates `0 <= value <= u16::MAX` — distributions that
///    return negative or oversized trims are caller bugs and
///    panic loudly with the offending value.
/// 3. Records `ChoiceValue::Int(value)` to the trace at address
///    `"trim.{segment}_{end}"` (e.g. `"trim.v_3"`, `"trim.j_5"`).
/// 4. Returns `sim.with_trim(segment, end, value as u16)`.
///
/// **Pre-condition:** an allele must already be assigned to
/// `segment` (i.e., the matching `SampleAllelePass` ran earlier in
/// the plan). `Simulation::with_trim` panics if the slot is empty,
/// and so does this pass by extension.
pub struct TrimPass {
    segment: Segment,
    end: TrimEnd,
    distribution: Box<dyn Distribution<Output = i64>>,
}

impl TrimPass {
    /// Construct a trim pass.
    ///
    /// Panics if `segment` is anything other than V, D, or J.
    pub fn new(
        segment: Segment,
        end: TrimEnd,
        distribution: Box<dyn Distribution<Output = i64>>,
    ) -> Self {
        match segment {
            Segment::V | Segment::D | Segment::J => {}
            _ => panic!("TrimPass: segment must be V, D, or J — got {:?}", segment),
        }
        Self {
            segment,
            end,
            distribution,
        }
    }

    pub fn segment(&self) -> Segment {
        self.segment
    }

    pub fn end(&self) -> TrimEnd {
        self.end
    }

    /// The hierarchical-string address (D3) at which this pass
    /// records its choice. Same string as `name()`.
    fn address(&self) -> &'static str {
        match (self.segment, self.end) {
            (Segment::V, TrimEnd::Five) => "trim.v_5",
            (Segment::V, TrimEnd::Three) => "trim.v_3",
            (Segment::D, TrimEnd::Five) => "trim.d_5",
            (Segment::D, TrimEnd::Three) => "trim.d_3",
            (Segment::J, TrimEnd::Five) => "trim.j_5",
            (Segment::J, TrimEnd::Three) => "trim.j_3",
            _ => unreachable!("TrimPass with non-V/D/J segment"),
        }
    }
}

impl Pass for TrimPass {
    fn name(&self) -> &str {
        self.address()
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        let value = self.distribution.sample(ctx.rng);
        assert!(
            value >= 0,
            "TrimPass({}): distribution returned negative trim {}",
            self.address(),
            value
        );
        assert!(
            value <= u16::MAX as i64,
            "TrimPass({}): distribution returned trim {} > u16::MAX",
            self.address(),
            value
        );

        ctx.trace.record(self.address(), ChoiceValue::Int(value));
        sim.with_trim(self.segment, self.end, value as u16)
    }

    fn declared_choices(&self) -> Vec<String> {
        vec![self.address().to_string()]
    }
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{NucHandle, Simulation};
    use crate::pass::{PassPlan, PassRuntime};

    #[test]
    fn echo_pass_appends_configured_nucleotide() {
        let plan = {
            let mut p = PassPlan::new();
            p.push(Box::new(EchoPass::new(b'C', 42, Segment::V)));
            p
        };
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);

        assert_eq!(outcome.revisions.len(), 2);
        let final_sim = outcome.final_simulation();
        assert_eq!(final_sim.pool.len(), 1);

        let n = final_sim.pool.get(NucHandle::new(0)).unwrap();
        assert_eq!(n.base, b'C');
        assert_eq!(n.germline, b'C'); // germline constructor: base == germline
        assert_eq!(n.germline_pos, 42);
        assert_eq!(n.segment, Segment::V);
    }

    #[test]
    fn echo_pass_records_nothing_to_trace() {
        // Determinism contract for transform passes: no RNG draws,
        // no trace writes. Verifies the design distinction between
        // transform and sampling passes.
        let plan = {
            let mut p = PassPlan::new();
            p.push(Box::new(EchoPass::new(b'A', 0, Segment::V)));
            p.push(Box::new(EchoPass::new(b'C', 1, Segment::V)));
            p.push(Box::new(EchoPass::new(b'G', 2, Segment::V)));
            p
        };
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);

        assert!(
            outcome.trace.is_empty(),
            "EchoPass produced trace entries: {:?}",
            outcome.trace.choices()
        );
    }

    #[test]
    fn echo_pass_is_deterministic_under_any_seed() {
        // Two runs with different seeds should produce identical
        // outcomes when only EchoPass instances are in the plan,
        // because EchoPass doesn't consume the RNG.
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(EchoPass::new(b'T', 5, Segment::J)));
            p.push(Box::new(EchoPass::new(b'A', 6, Segment::J)));
            p
        };

        let o1 = PassRuntime::execute(&plan(), Simulation::new(), 1);
        let o2 = PassRuntime::execute(&plan(), Simulation::new(), 0xdead_beef);

        assert_eq!(o1.revisions.len(), o2.revisions.len());
        for i in 0..o1.revisions[2].pool.len() {
            let h = NucHandle::new(i as u32);
            assert_eq!(
                o1.revisions[2].pool.get(h).unwrap().base,
                o2.revisions[2].pool.get(h).unwrap().base
            );
        }
    }

    #[test]
    fn echo_pass_composes_with_pass_plan_in_order() {
        let plan = {
            let mut p = PassPlan::new();
            p.push(Box::new(EchoPass::new(b'A', 0, Segment::V)));
            p.push(Box::new(EchoPass::new(b'C', 1, Segment::V)));
            p.push(Box::new(EchoPass::new(b'G', 2, Segment::V)));
            p.push(Box::new(EchoPass::new(b'T', 3, Segment::V)));
            p
        };
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);

        // 5 revisions: initial + 4 echo passes.
        assert_eq!(outcome.revisions.len(), 5);
        assert_eq!(outcome.pass_names, vec!["echo"; 4]);

        // Each successive revision has one more nucleotide.
        for (i, rev) in outcome.revisions.iter().enumerate() {
            assert_eq!(rev.pool.len(), i);
        }

        // Bases pushed in order ACGT.
        let final_pool = &outcome.final_simulation().pool;
        for (i, &expected) in [b'A', b'C', b'G', b'T'].iter().enumerate() {
            assert_eq!(
                final_pool.get(NucHandle::new(i as u32)).unwrap().base,
                expected
            );
        }
    }

    #[test]
    fn echo_pass_accessors_round_trip() {
        let p = EchoPass::new(b'G', 17, Segment::Np1);
        assert_eq!(p.base(), b'G');
        assert_eq!(p.germline_pos(), 17);
        assert_eq!(p.segment(), Segment::Np1);
        assert_eq!(p.name(), "echo");
    }

    #[test]
    fn echo_pass_declares_no_choices() {
        // Transform passes (no RNG draws) inherit the default empty
        // declared_choices(). Phase D's validator + upstream-bound
        // propagation will skip them by virtue of getting an empty
        // address list.
        let p = EchoPass::new(b'A', 0, Segment::V);
        assert!(p.declared_choices().is_empty());
    }

    #[test]
    fn echo_pass_preserves_history_chain_under_persistent_ir() {
        // D1 contract at the integration level: every revision retains
        // its own pool size after subsequent passes have run.
        let plan = {
            let mut p = PassPlan::new();
            p.push(Box::new(EchoPass::new(b'A', 0, Segment::V)));
            p.push(Box::new(EchoPass::new(b'C', 1, Segment::V)));
            p.push(Box::new(EchoPass::new(b'G', 2, Segment::V)));
            p
        };
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);

        // Initial revision: empty pool, even after 3 passes ran.
        assert_eq!(outcome.revisions[0].pool.len(), 0);
        // After first echo: 1 nucleotide, never grows.
        assert_eq!(outcome.revisions[1].pool.len(), 1);
        assert_eq!(
            outcome.revisions[1]
                .pool
                .get(NucHandle::new(0))
                .unwrap()
                .base,
            b'A'
        );
        // After second echo: 2 nucleotides; the original 'A' is unchanged.
        assert_eq!(outcome.revisions[2].pool.len(), 2);
        assert_eq!(
            outcome.revisions[2]
                .pool
                .get(NucHandle::new(0))
                .unwrap()
                .base,
            b'A'
        );
        assert_eq!(
            outcome.revisions[2]
                .pool
                .get(NucHandle::new(1))
                .unwrap()
                .base,
            b'C'
        );
    }

    // ── SampleBasePass tests ───────────────────────────────────────

    use crate::dist::UniformBase;

    #[test]
    fn sample_base_pass_appends_one_synthetic_nucleotide() {
        let pass = SampleBasePass::new(
            "test.np1.bases[0]",
            Box::new(UniformBase),
            Segment::Np1,
            crate::ir::flag::N_NUC,
        );
        let mut plan = PassPlan::new();
        plan.push(Box::new(pass));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 7);

        let final_sim = outcome.final_simulation();
        assert_eq!(final_sim.pool.len(), 1);
        let n = final_sim.pool.get(NucHandle::new(0)).unwrap();

        // Synthetic constructor: no germline provenance.
        assert_eq!(n.germline_pos, Nucleotide::NO_GERMLINE_POS);
        assert_eq!(n.segment, Segment::Np1);
        assert!(n.flags.contains(crate::ir::flag::N_NUC));
        // Base is one of A/C/G/T (UniformBase).
        assert!(matches!(n.base, b'A' | b'C' | b'G' | b'T'));
    }

    #[test]
    fn sample_base_pass_records_choice_at_configured_address() {
        let pass = SampleBasePass::new(
            "test.np1.bases[0]",
            Box::new(UniformBase),
            Segment::Np1,
            crate::ir::flag::N_NUC,
        );
        let mut plan = PassPlan::new();
        plan.push(Box::new(pass));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 7);

        // Trace contains exactly one entry, at the configured address.
        assert_eq!(outcome.trace.len(), 1);
        let rec = outcome.trace.find("test.np1.bases[0]").unwrap();
        match &rec.value {
            ChoiceValue::Base(b) => {
                assert!(matches!(*b, b'A' | b'C' | b'G' | b'T'));
            }
            other => panic!("expected ChoiceValue::Base, got {:?}", other),
        }
    }

    #[test]
    fn sample_base_pass_declares_its_address() {
        // Sampling passes override declared_choices() to report the
        // address(es) they will draw at. The runtime + Phase D
        // validator depend on this introspection.
        let pass = SampleBasePass::new(
            "test.np1.bases[0]",
            Box::new(UniformBase),
            Segment::Np1,
            crate::ir::flag::N_NUC,
        );
        assert_eq!(
            pass.declared_choices(),
            vec!["test.np1.bases[0]".to_string()]
        );
    }

    #[test]
    fn sample_base_pass_recorded_value_matches_pushed_nucleotide() {
        // Faithfulness: the trace must record exactly what was written
        // to the IR. Replay against a divergent record would produce
        // a different sequence; this test pins the invariant.
        let mut plan = PassPlan::new();
        for i in 0..10 {
            plan.push(Box::new(SampleBasePass::new(
                format!("test.bases[{}]", i),
                Box::new(UniformBase),
                Segment::Np1,
                crate::ir::flag::N_NUC,
            )));
        }

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0xdead_beef);
        let final_sim = outcome.final_simulation();

        for i in 0..10 {
            let addr = format!("test.bases[{}]", i);
            let recorded = match outcome.trace.find(&addr).unwrap().value {
                ChoiceValue::Base(b) => b,
                _ => panic!("wrong variant"),
            };
            let pushed = final_sim.pool.get(NucHandle::new(i)).unwrap().base;
            assert_eq!(
                recorded, pushed,
                "trace entry at {} = {} but IR has {} at handle {}",
                addr, recorded as char, pushed as char, i
            );
        }
    }

    // ── Phase B milestone: replay determinism ──────────────────────

    /// Build a representative mixed plan: alternating EchoPass (no
    /// RNG) and SampleBasePass (RNG-consuming). Used by the replay
    /// determinism test below.
    fn mixed_plan() -> PassPlan {
        let mut plan = PassPlan::new();
        for i in 0..8 {
            plan.push(Box::new(EchoPass::new(b'A', i as u16, Segment::V)));
            plan.push(Box::new(SampleBasePass::new(
                format!("np.np1.bases[{}]", i),
                Box::new(UniformBase),
                Segment::Np1,
                crate::ir::flag::N_NUC,
            )));
        }
        plan
    }

    #[test]
    fn replay_determinism_same_seed_same_trace_same_ir() {
        // The Phase B success criterion: two independent runs of the
        // same plan with the same seed must produce identical traces
        // and identical final IR revisions, byte for byte.
        let oa = PassRuntime::execute(&mixed_plan(), Simulation::new(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&mixed_plan(), Simulation::new(), 0xc0ff_ee);

        // Trace equality: same addresses, same values, same order.
        assert_eq!(oa.trace.len(), ob.trace.len());
        assert_eq!(oa.trace.choices().len(), 8); // 8 SampleBase + 8 Echo (no record)
        for (a, b) in oa.trace.choices().iter().zip(ob.trace.choices().iter()) {
            assert_eq!(a.address, b.address);
            assert_eq!(a.value, b.value);
        }

        // Final IR equality: same pool size, identical bases at every
        // handle.
        let pa = &oa.final_simulation().pool;
        let pb = &ob.final_simulation().pool;
        assert_eq!(pa.len(), pb.len());
        for i in 0..pa.len() {
            let h = NucHandle::new(i as u32);
            let na = pa.get(h).unwrap();
            let nb = pb.get(h).unwrap();
            assert_eq!(na.base, nb.base);
            assert_eq!(na.germline, nb.germline);
            assert_eq!(na.germline_pos, nb.germline_pos);
            assert_eq!(na.segment, nb.segment);
            assert_eq!(na.flags, nb.flags);
        }

        // Pass-name list also identical.
        assert_eq!(oa.pass_names, ob.pass_names);
    }

    #[test]
    fn replay_determinism_different_seed_diverges() {
        let oa = PassRuntime::execute(&mixed_plan(), Simulation::new(), 1);
        let ob = PassRuntime::execute(&mixed_plan(), Simulation::new(), 2);

        // Trace addresses should still match (plan is identical).
        for (a, b) in oa.trace.choices().iter().zip(ob.trace.choices().iter()) {
            assert_eq!(a.address, b.address);
        }

        // At least one sampled base must differ between the two seeds —
        // 8 independent draws have ~vanishing probability of full agreement.
        let any_diff = oa
            .trace
            .choices()
            .iter()
            .zip(ob.trace.choices().iter())
            .any(|(a, b)| a.value != b.value);
        assert!(
            any_diff,
            "different seeds produced byte-identical traces — RNG plumbing broken"
        );
    }

    // ── SampleAllelePass tests (C.5) ───────────────────────────────

    use crate::dist::AllelePoolDist;
    use crate::refdata::{Allele, AllelePool};

    /// Build an allele pool of `n` named alleles for testing.
    fn make_test_pool(n: usize, segment: Segment) -> AllelePool {
        let mut p = AllelePool::new();
        for i in 0..n {
            let _ = p.push(Allele {
                name: format!("test_allele_{}*01", i),
                gene: format!("test_allele_{}", i),
                seq: vec![b'A'; 30],
                segment,
                anchor: Some(10),
            });
        }
        p
    }

    #[test]
    #[should_panic(expected = "segment must be V, D, or J")]
    fn sample_allele_pass_rejects_np1() {
        let pool = make_test_pool(1, Segment::V);
        let _ = SampleAllelePass::new(Segment::Np1, Box::new(AllelePoolDist::uniform(&pool)));
    }

    #[test]
    #[should_panic(expected = "segment must be V, D, or J")]
    fn sample_allele_pass_rejects_np2() {
        let pool = make_test_pool(1, Segment::V);
        let _ = SampleAllelePass::new(Segment::Np2, Box::new(AllelePoolDist::uniform(&pool)));
    }

    #[test]
    fn sample_allele_pass_assigns_to_correct_slot_for_v() {
        let pool = make_test_pool(1, Segment::V);
        let dist = Box::new(AllelePoolDist::uniform(&pool));
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(Segment::V, dist)));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 42);

        let final_sim = outcome.final_simulation();
        // Single-allele dist always returns AlleleId(0).
        assert_eq!(
            final_sim.assignments.v.unwrap().allele_id,
            crate::refdata::AlleleId::new(0)
        );
        assert!(final_sim.assignments.d.is_none());
        assert!(final_sim.assignments.j.is_none());
        // Default trims.
        assert_eq!(final_sim.assignments.v.unwrap().trim_5, 0);
        assert_eq!(final_sim.assignments.v.unwrap().trim_3, 0);
    }

    #[test]
    fn sample_allele_pass_records_to_trace_at_segment_address() {
        let v_pool = make_test_pool(1, Segment::V);
        let d_pool = make_test_pool(1, Segment::D);
        let j_pool = make_test_pool(1, Segment::J);

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::D,
            Box::new(AllelePoolDist::uniform(&d_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&j_pool)),
        )));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 7);

        // Three choices, one per pass, at the canonical addresses.
        assert_eq!(outcome.trace.len(), 3);
        for addr in ["sample_allele.v", "sample_allele.d", "sample_allele.j"] {
            let rec = outcome
                .trace
                .find(addr)
                .unwrap_or_else(|| panic!("missing {}", addr));
            match rec.value {
                ChoiceValue::AlleleId(id) => {
                    assert_eq!(id, 0); // single-allele pool
                }
                _ => panic!("wrong variant at {}", addr),
            }
        }
    }

    #[test]
    fn sample_allele_pass_is_deterministic_under_same_seed() {
        let pool = make_test_pool(10, Segment::V);
        let mut plan_a = PassPlan::new();
        plan_a.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&pool)),
        )));
        let mut plan_b = PassPlan::new();
        plan_b.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&pool)),
        )));

        let oa = PassRuntime::execute(&plan_a, Simulation::new(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&plan_b, Simulation::new(), 0xc0ff_ee);

        assert_eq!(
            oa.final_simulation().assignments.v.unwrap().allele_id,
            ob.final_simulation().assignments.v.unwrap().allele_id
        );
        assert_eq!(oa.trace.choices()[0].value, ob.trace.choices()[0].value);
    }

    #[test]
    fn sample_allele_pass_full_recombination_chain_for_vdj() {
        // V + D + J sampling for a heavy chain: all three
        // assignments should be populated after the plan runs.
        let v_pool = make_test_pool(5, Segment::V);
        let d_pool = make_test_pool(3, Segment::D);
        let j_pool = make_test_pool(2, Segment::J);

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::D,
            Box::new(AllelePoolDist::uniform(&d_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&j_pool)),
        )));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 42);

        let sim = outcome.final_simulation();
        assert!(sim.assignments.v.is_some());
        assert!(sim.assignments.d.is_some());
        assert!(sim.assignments.j.is_some());
        assert!(sim.assignments.c.is_none());

        // Sampled ids are in-bounds for their pools (D-binding from C.3).
        assert!(sim.assignments.v.unwrap().allele_id.as_usize() < 5);
        assert!(sim.assignments.d.unwrap().allele_id.as_usize() < 3);
        assert!(sim.assignments.j.unwrap().allele_id.as_usize() < 2);
    }

    #[test]
    fn sample_allele_pass_declared_choices_returns_address() {
        let pool = make_test_pool(1, Segment::V);
        let pass_v = SampleAllelePass::new(Segment::V, Box::new(AllelePoolDist::uniform(&pool)));
        assert_eq!(
            pass_v.declared_choices(),
            vec!["sample_allele.v".to_string()]
        );

        let d_pool = make_test_pool(1, Segment::D);
        let pass_d = SampleAllelePass::new(Segment::D, Box::new(AllelePoolDist::uniform(&d_pool)));
        assert_eq!(
            pass_d.declared_choices(),
            vec!["sample_allele.d".to_string()]
        );
    }

    #[test]
    fn sample_allele_pass_segment_accessor() {
        let pool = make_test_pool(1, Segment::J);
        let pass = SampleAllelePass::new(Segment::J, Box::new(AllelePoolDist::uniform(&pool)));
        assert_eq!(pass.segment(), Segment::J);
        assert_eq!(pass.name(), "sample_allele.j");
    }

    // ── TrimPass tests (C.6) ───────────────────────────────────────

    use crate::dist::EmpiricalLengthDist;

    #[test]
    #[should_panic(expected = "segment must be V, D, or J")]
    fn trim_pass_rejects_np1() {
        let _ = TrimPass::new(
            Segment::Np1,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
        );
    }

    #[test]
    fn trim_pass_addresses_cover_all_six_combinations() {
        let dist = || Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)]));

        let cases = [
            (Segment::V, TrimEnd::Five, "trim.v_5"),
            (Segment::V, TrimEnd::Three, "trim.v_3"),
            (Segment::D, TrimEnd::Five, "trim.d_5"),
            (Segment::D, TrimEnd::Three, "trim.d_3"),
            (Segment::J, TrimEnd::Five, "trim.j_5"),
            (Segment::J, TrimEnd::Three, "trim.j_3"),
        ];
        for (seg, end, expected_addr) in cases {
            let pass = TrimPass::new(seg, end, dist());
            assert_eq!(pass.name(), expected_addr);
            assert_eq!(pass.declared_choices(), vec![expected_addr.to_string()]);
        }
    }

    #[test]
    fn trim_pass_applies_trim_and_records_choice() {
        // V allele assigned, then trim_3 sampled.
        let v_pool = make_test_pool(1, Segment::V);
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            // Single-value dist: always returns 4.
            Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
        )));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 7);

        // Trim recorded at canonical address.
        let rec = outcome.trace.find("trim.v_3").expect("trim recorded");
        assert_eq!(rec.value, ChoiceValue::Int(4));

        // Trim applied to the V allele instance.
        let v = outcome.final_simulation().assignments.v.unwrap();
        assert_eq!(v.trim_3, 4);
        assert_eq!(v.trim_5, 0);
    }

    #[test]
    #[should_panic(expected = "no instance assigned to segment")]
    fn trim_pass_panics_when_no_allele_assigned() {
        // The trim pass relies on `Simulation::with_trim`, which
        // panics when no instance is assigned to the segment. This
        // test pins the propagation.
        let mut plan = PassPlan::new();
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
        )));
        let _ = PassRuntime::execute(&plan, Simulation::new(), 0);
    }

    #[test]
    #[should_panic(expected = "negative trim")]
    fn trim_pass_panics_on_negative_distribution_output() {
        // Runtime defensive check. UniformInt::new(-5, -1) always
        // returns a negative trim, which TrimPass must reject loudly.
        use crate::dist::UniformInt;
        let v_pool = make_test_pool(1, Segment::V);
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(UniformInt::new(-5, -1)),
        )));
        let _ = PassRuntime::execute(&plan, Simulation::new(), 1);
    }

    #[test]
    fn trim_pass_full_vdj_chain_records_six_choices() {
        // Heavy chain shape: sample V/D/J, then trim each side.
        let v_pool = make_test_pool(1, Segment::V);
        let d_pool = make_test_pool(1, Segment::D);
        let j_pool = make_test_pool(1, Segment::J);
        let dist = || Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)]));

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::D,
            Box::new(AllelePoolDist::uniform(&d_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&j_pool)),
        )));
        plan.push(Box::new(TrimPass::new(Segment::V, TrimEnd::Three, dist())));
        plan.push(Box::new(TrimPass::new(Segment::D, TrimEnd::Five, dist())));
        plan.push(Box::new(TrimPass::new(Segment::D, TrimEnd::Three, dist())));
        plan.push(Box::new(TrimPass::new(Segment::J, TrimEnd::Five, dist())));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 42);

        // 3 sampling + 4 trim = 7 trace entries.
        assert_eq!(outcome.trace.len(), 7);

        // All four trim addresses present with value 2.
        for addr in ["trim.v_3", "trim.d_5", "trim.d_3", "trim.j_5"] {
            assert_eq!(outcome.trace.find(addr).unwrap().value, ChoiceValue::Int(2));
        }

        // Trims applied to the assignments.
        let sim = outcome.final_simulation();
        assert_eq!(sim.assignments.v.unwrap().trim_3, 2);
        assert_eq!(sim.assignments.d.unwrap().trim_5, 2);
        assert_eq!(sim.assignments.d.unwrap().trim_3, 2);
        assert_eq!(sim.assignments.j.unwrap().trim_5, 2);
    }

    #[test]
    fn trim_pass_is_deterministic_under_same_seed() {
        let v_pool = make_test_pool(1, Segment::V);
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(SampleAllelePass::new(
                Segment::V,
                Box::new(AllelePoolDist::uniform(&v_pool)),
            )));
            p.push(Box::new(TrimPass::new(
                Segment::V,
                TrimEnd::Three,
                Box::new(EmpiricalLengthDist::from_pairs(vec![
                    (0, 1.0),
                    (1, 2.0),
                    (2, 3.0),
                    (5, 1.0),
                ])),
            )));
            p
        };

        let oa = PassRuntime::execute(&plan(), Simulation::new(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&plan(), Simulation::new(), 0xc0ff_ee);
        assert_eq!(
            oa.final_simulation().assignments.v.unwrap().trim_3,
            ob.final_simulation().assignments.v.unwrap().trim_3
        );
        assert_eq!(oa.trace.choices(), ob.trace.choices());
    }

    // ── UniformMutationPass tests (E.1) ────────────────────────────

    #[derive(Clone, Debug)]
    struct StopThenSafeMutationBaseDist;

    impl Distribution for StopThenSafeMutationBaseDist {
        type Output = u8;

        fn sample(&self, _rng: &mut crate::rng::Rng) -> u8 {
            b'A'
        }

        fn support(&self) -> Option<Vec<(u8, f64)>> {
            Some(vec![(b'A', 1.0), (b'C', 1.0)])
        }
    }

    #[derive(Clone, Debug)]
    struct StopOnlyMutationBaseDist;

    impl Distribution for StopOnlyMutationBaseDist {
        type Output = u8;

        fn sample(&self, _rng: &mut crate::rng::Rng) -> u8 {
            b'A'
        }

        fn support(&self) -> Option<Vec<(u8, f64)>> {
            Some(vec![(b'A', 1.0)])
        }
    }

    fn uniform_mutation_plan(base_dist: Box<dyn Distribution<Output = u8>>) -> PassPlan {
        let mut plan = PassPlan::new();
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
            base_dist,
        )));
        plan
    }

    fn make_uniform_mutation_productive_vj_fixture() -> (RefDataConfig, Simulation) {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(Allele {
            name: "v_mut*01".into(),
            gene: "v_mut".into(),
            seq: b"TAC".to_vec(),
            segment: Segment::V,
            anchor: Some(0),
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j_mut*01".into(),
            gene: "j_mut".into(),
            seq: b"TGG".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });

        let mut sim = Simulation::new();
        for (i, &b) in b"TAC".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
            sim = next;
        }
        let v_region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(3))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(v_region);

        for (i, &b) in b"TGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::J));
            sim = next;
        }
        let j_region = Region::new(Segment::J, NucHandle::new(3), NucHandle::new(6))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(j_region);

        sim = sim
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

        (cfg, sim)
    }

    fn uniform_mutation_site(seed: u64, sim: Simulation) -> u32 {
        let outcome = PassRuntime::execute(
            &uniform_mutation_plan(Box::new(StopOnlyMutationBaseDist)),
            sim,
            seed,
        );
        match outcome.trace.find("mutate.uniform.site[0]").unwrap().value {
            ChoiceValue::Int(site) => site as u32,
            _ => panic!("wrong variant"),
        }
    }

    fn find_seed_for_uniform_mutation_site(sim: &Simulation, target_site: u32) -> u64 {
        for seed in 0..512u64 {
            if uniform_mutation_site(seed, sim.clone()) == target_site {
                return seed;
            }
        }
        panic!("no seed in search range targeted site {}", target_site);
    }

    #[test]
    fn uniform_mutation_pass_zero_count_is_noop() {
        // Build a sim with a few nucleotides + region, run uniform
        // mutation with count=0, verify nothing changed.
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region);

        let mut plan = PassPlan::new();
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
            Box::new(UniformBase),
        )));

        let outcome = PassRuntime::execute(&plan, sim.clone(), 42);
        let final_sim = outcome.final_simulation();

        // Pool unchanged.
        for i in 0..9 {
            assert_eq!(
                final_sim.pool.get(NucHandle::new(i)).unwrap().base,
                sim.pool.get(NucHandle::new(i)).unwrap().base
            );
        }
        // Codon rail unchanged.
        assert_eq!(
            final_sim.sequence.regions[0].amino_acids,
            sim.sequence.regions[0].amino_acids
        );
        // Trace recorded count=0, no site/base entries.
        assert_eq!(outcome.trace.len(), 1);
        assert_eq!(
            outcome.trace.find("mutate.uniform.count").unwrap().value,
            ChoiceValue::Int(0)
        );
    }

    #[test]
    fn uniform_mutation_pass_applies_n_mutations_with_traced_addresses() {
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGGTTT".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(12))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region);

        let mut plan = PassPlan::new();
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
            Box::new(UniformBase),
        )));

        let outcome = PassRuntime::execute(&plan, sim, 7);

        // Trace: 1 count + 5 site + 5 base = 11 records.
        assert_eq!(outcome.trace.len(), 11);
        assert_eq!(
            outcome.trace.find("mutate.uniform.count").unwrap().value,
            ChoiceValue::Int(5)
        );
        for i in 0..5 {
            let site_addr = format!("mutate.uniform.site[{}]", i);
            let base_addr = format!("mutate.uniform.base[{}]", i);
            assert!(outcome.trace.find(&site_addr).is_some());
            assert!(outcome.trace.find(&base_addr).is_some());

            // Site is in [0, 12).
            match outcome.trace.find(&site_addr).unwrap().value {
                ChoiceValue::Int(s) => assert!(s >= 0 && s < 12),
                _ => panic!("wrong variant"),
            }
            // Base is one of A/C/G/T.
            match outcome.trace.find(&base_addr).unwrap().value {
                ChoiceValue::Base(b) => assert!(matches!(b, b'A' | b'C' | b'G' | b'T')),
                _ => panic!("wrong variant"),
            }
        }
    }

    #[test]
    fn uniform_mutation_pass_pool_reflects_recorded_mutations() {
        // Faithfulness: the pool's base at each recorded site equals
        // the recorded new base (trace is honest).
        let mut sim = Simulation::new();
        for (i, b) in b"AAAAAAAAAA".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }

        let mut plan = PassPlan::new();
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(7, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, sim, 0xc0ff_ee);
        let final_sim = outcome.final_simulation();

        // Walk the trace: at each (site[i], base[i]) the pool should
        // hold base[i]. (NOTE: later mutations could overwrite earlier
        // ones at the same site, so we can only check the LAST
        // mutation per site.)
        let mut last_at_site = std::collections::HashMap::new();
        for i in 0..7 {
            let s = match outcome
                .trace
                .find(&format!("mutate.uniform.site[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Int(n) => n as u32,
                _ => unreachable!(),
            };
            let b = match outcome
                .trace
                .find(&format!("mutate.uniform.base[{}]", i))
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
    fn uniform_mutation_pass_refreshes_codon_rail_through_persistent_api() {
        // The post-D audit fix in action: every `with_base_changed`
        // call inside the pass updates Region.amino_acids
        // automatically. After a mutation pass runs, the codon rail
        // should still match a fresh recomputation against the pool.
        let mut sim = Simulation::new();
        for (i, b) in b"ATGGGGAAATTT".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(12))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region);

        let mut plan = PassPlan::new();
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(8, 1.0)])),
            Box::new(UniformBase),
        )));

        let outcome = PassRuntime::execute(&plan, sim, 99);
        let final_sim = outcome.final_simulation();

        // Stored codon rail equals a fresh recomputation against
        // the post-mutation pool. This is the staleness invariant.
        let stored_aa = &final_sim.sequence.regions[0].amino_acids;
        let fresh_region =
            final_sim.sequence.regions[0].with_codon_rail_recomputed(&final_sim.pool);
        assert_eq!(stored_aa, &fresh_region.amino_acids);
        assert_eq!(
            final_sim.sequence.regions[0].stop_codon_positions,
            fresh_region.stop_codon_positions
        );
    }

    #[test]
    fn uniform_mutation_pass_is_deterministic_under_same_seed() {
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(UniformMutationPass::new(
                Box::new(EmpiricalLengthDist::from_pairs(vec![
                    (3, 1.0),
                    (5, 2.0),
                    (7, 1.0),
                ])),
                Box::new(UniformBase),
            )));
            p
        };
        let build_sim = || {
            let mut s = Simulation::new();
            for (i, b) in b"AAACCCGGGTTTAAA".iter().enumerate() {
                let (next, _) =
                    s.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
                s = next;
            }
            s
        };

        let oa = PassRuntime::execute(&plan(), build_sim(), 0xfeed);
        let ob = PassRuntime::execute(&plan(), build_sim(), 0xfeed);

        assert_eq!(oa.trace.choices(), ob.trace.choices());
        // Pool byte-identical.
        for i in 0..15 {
            assert_eq!(
                oa.final_simulation()
                    .pool
                    .get(NucHandle::new(i))
                    .unwrap()
                    .base,
                ob.final_simulation()
                    .pool
                    .get(NucHandle::new(i))
                    .unwrap()
                    .base
            );
        }
    }

    #[test]
    fn uniform_mutation_pass_preserves_persistent_ir() {
        // Original sim must be untouched after the mutation pass
        // runs against a clone of it.
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCC".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let pre_pool: Vec<u8> = (0..6)
            .map(|i| sim.pool.get(NucHandle::new(i)).unwrap().base)
            .collect();

        let mut plan = PassPlan::new();
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(20, 1.0)])),
            Box::new(UniformBase),
        )));
        let _ = PassRuntime::execute(&plan, sim.clone(), 5);

        // Verify original sim unchanged.
        let post_pool: Vec<u8> = (0..6)
            .map(|i| sim.pool.get(NucHandle::new(i)).unwrap().base)
            .collect();
        assert_eq!(pre_pool, post_pool);
    }

    #[test]
    fn uniform_mutation_pass_declared_choices_uses_indexed_pattern() {
        let pass = UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        );
        let declared = pass.declared_choices();
        assert_eq!(declared.len(), 3);
        assert!(declared.contains(&"mutate.uniform.count".to_string()));
        assert!(declared.contains(&"mutate.uniform.site[0..n]".to_string()));
        assert!(declared.contains(&"mutate.uniform.base[0..n]".to_string()));
    }

    #[test]
    fn uniform_mutation_productive_filters_base_that_would_create_stop() {
        let (cfg, sim) = make_uniform_mutation_productive_vj_fixture();
        let seed = find_seed_for_uniform_mutation_site(&sim, 2);
        let contracts = productive();

        let constrained = PassRuntime::execute_with_context(
            &uniform_mutation_plan(Box::new(StopThenSafeMutationBaseDist)),
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        );

        assert_eq!(
            constrained
                .trace
                .find("mutate.uniform.site[0]")
                .unwrap()
                .value,
            ChoiceValue::Int(2)
        );
        assert_eq!(
            constrained
                .trace
                .find("mutate.uniform.base[0]")
                .unwrap()
                .value,
            ChoiceValue::Base(b'C')
        );
        assert!(contracts
            .verify(constrained.final_simulation(), Some(&cfg))
            .is_ok());

        let unconstrained = PassRuntime::execute_with_context(
            &uniform_mutation_plan(Box::new(StopThenSafeMutationBaseDist)),
            sim,
            seed,
            Some(&cfg),
            None,
        );
        assert_eq!(
            unconstrained
                .trace
                .find("mutate.uniform.base[0]")
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
    fn uniform_mutation_strict_errors_when_base_filter_empty() {
        let (cfg, sim) = make_uniform_mutation_productive_vj_fixture();
        let seed = find_seed_for_uniform_mutation_site(&sim, 2);
        let contracts = productive();

        let err = PassRuntime::execute_strict_with_context(
            &uniform_mutation_plan(Box::new(StopOnlyMutationBaseDist)),
            sim,
            seed,
            Some(&cfg),
            Some(&contracts),
        )
        .unwrap_err();

        assert_eq!(err.pass_name(), "mutate.uniform");
        assert_eq!(err.address(), "mutate.uniform.base[0]");
        assert_eq!(
            err.constraint_reason(),
            Some(FilteredSampleError::EmptyAdmissibleSupport)
        );
    }

    // ── ContaminantPass tests (E.6) ────────────────────────────────

    fn contaminant_test_sim() -> Simulation {
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        sim
    }

    #[derive(Clone, Debug)]
    struct StopOnlyContaminantBaseDist;

    impl Distribution for StopOnlyContaminantBaseDist {
        type Output = u8;

        fn sample(&self, _rng: &mut crate::rng::Rng) -> u8 {
            b'T'
        }

        fn support(&self) -> Option<Vec<(u8, f64)>> {
            Some(vec![(b'T', 1.0)])
        }
    }

    fn contaminant_plan(base_dist: Box<dyn Distribution<Output = u8>>) -> PassPlan {
        let mut plan = PassPlan::new();
        plan.push(Box::new(ContaminantPass::new(1.0, base_dist)));
        plan
    }

    fn make_contaminant_productive_vj_fixture() -> (RefDataConfig, Simulation) {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(Allele {
            name: "v_contam*01".into(),
            gene: "v_contam".into(),
            seq: b"AAA".to_vec(),
            segment: Segment::V,
            anchor: Some(0),
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j_contam*01".into(),
            gene: "j_contam".into(),
            seq: b"TGG".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });

        let mut sim = Simulation::new();
        for (i, &b) in b"AAA".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
            sim = next;
        }
        let v_region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(3))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(v_region);

        for (i, &b) in b"TGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::J));
            sim = next;
        }
        let j_region = Region::new(Segment::J, NucHandle::new(3), NucHandle::new(6))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(j_region);

        sim = sim
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

        (cfg, sim)
    }

    fn find_seed_for_unconstrained_contaminant_prefix(sim: &Simulation, expected: &[u8]) -> u64 {
        for seed in 0..4096u64 {
            let outcome =
                PassRuntime::execute(&contaminant_plan(Box::new(UniformBase)), sim.clone(), seed);
            let matches = expected.iter().enumerate().all(|(i, &base)| {
                outcome
                    .trace
                    .find(&format!("corrupt.contaminant.bases[{}]", i))
                    .map(|rec| rec.value == ChoiceValue::Base(base))
                    .unwrap_or(false)
            });
            if matches {
                return seed;
            }
        }
        panic!(
            "no seed in search range produced contaminant prefix {:?}",
            expected
        );
    }

    #[test]
    #[should_panic(expected = "apply_prob must be in [0.0, 1.0]")]
    fn contaminant_pass_rejects_negative_probability() {
        let _ = ContaminantPass::new(-0.1, Box::new(UniformBase));
    }

    #[test]
    #[should_panic(expected = "apply_prob must be in [0.0, 1.0]")]
    fn contaminant_pass_rejects_probability_above_one() {
        let _ = ContaminantPass::new(1.5, Box::new(UniformBase));
    }

    #[test]
    #[should_panic(expected = "apply_prob must be in [0.0, 1.0]")]
    fn contaminant_pass_rejects_nan_probability() {
        let _ = ContaminantPass::new(f64::NAN, Box::new(UniformBase));
    }

    #[test]
    fn contaminant_pass_zero_probability_never_applies() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(ContaminantPass::new(0.0, Box::new(UniformBase))));

        for seed in 0..20u64 {
            let outcome = PassRuntime::execute(&plan, contaminant_test_sim(), seed);
            // Trace records `applied: Bool(false)` and nothing else.
            assert_eq!(outcome.trace.len(), 1);
            assert_eq!(
                outcome
                    .trace
                    .find("corrupt.contaminant.applied")
                    .unwrap()
                    .value,
                ChoiceValue::Bool(false)
            );
            // Pool unchanged.
            for i in 0..9 {
                let b = outcome
                    .final_simulation()
                    .pool
                    .get(NucHandle::new(i as u32))
                    .unwrap()
                    .base;
                assert!(matches!(b, b'A' | b'C' | b'G'));
            }
        }
    }

    #[test]
    fn contaminant_pass_one_probability_always_applies() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(ContaminantPass::new(1.0, Box::new(UniformBase))));

        for seed in 0..20u64 {
            let outcome = PassRuntime::execute(&plan, contaminant_test_sim(), seed);
            // Trace records `applied: Bool(true)` followed by 9 base records.
            assert_eq!(outcome.trace.len(), 10);
            assert_eq!(
                outcome
                    .trace
                    .find("corrupt.contaminant.applied")
                    .unwrap()
                    .value,
                ChoiceValue::Bool(true)
            );
            for i in 0..9 {
                assert!(outcome
                    .trace
                    .find(&format!("corrupt.contaminant.bases[{}]", i))
                    .is_some());
            }
        }
    }

    #[test]
    fn contaminant_pass_one_probability_pool_reflects_recorded_bases() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(ContaminantPass::new(1.0, Box::new(UniformBase))));
        let outcome = PassRuntime::execute(&plan, contaminant_test_sim(), 42);
        let final_sim = outcome.final_simulation();

        for i in 0..9u32 {
            let recorded = match outcome
                .trace
                .find(&format!("corrupt.contaminant.bases[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Base(b) => b,
                _ => unreachable!(),
            };
            let pool_base = final_sim.pool.get(NucHandle::new(i)).unwrap().base;
            assert_eq!(recorded, pool_base, "trace lies at position {}", i);
        }
    }

    #[test]
    fn contaminant_pass_half_probability_mixed_outcomes() {
        // Across many seeds, p=0.5 should produce both applied and
        // not-applied outcomes. Negative-control proof that the
        // coin flip is actually firing.
        let mut plan = PassPlan::new();
        plan.push(Box::new(ContaminantPass::new(0.5, Box::new(UniformBase))));

        let mut applied_count = 0;
        let mut not_applied_count = 0;
        for seed in 0..50u64 {
            let outcome = PassRuntime::execute(&plan, contaminant_test_sim(), seed);
            match outcome
                .trace
                .find("corrupt.contaminant.applied")
                .unwrap()
                .value
            {
                ChoiceValue::Bool(true) => applied_count += 1,
                ChoiceValue::Bool(false) => not_applied_count += 1,
                _ => unreachable!(),
            }
        }
        assert!(applied_count > 0, "expected at least one applied outcome");
        assert!(
            not_applied_count > 0,
            "expected at least one not-applied outcome"
        );
    }

    #[test]
    fn contaminant_pass_codon_rail_refresh_after_application() {
        // When contamination applies, the assembled region's codon
        // rail must reflect the new (contaminant) bases.
        let mut sim = Simulation::new();
        for (i, b) in b"ATGGGGGGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region);
        // Before contamination: M G G.
        assert_eq!(sim.sequence.regions[0].amino_acids, b"MGG");

        let mut plan = PassPlan::new();
        plan.push(Box::new(ContaminantPass::new(1.0, Box::new(UniformBase))));
        let outcome = PassRuntime::execute(&plan, sim, 0);
        let final_sim = outcome.final_simulation();

        // After contamination: codon rail recomputed, may differ.
        let stored = &final_sim.sequence.regions[0].amino_acids;
        let fresh = final_sim.sequence.regions[0].with_codon_rail_recomputed(&final_sim.pool);
        assert_eq!(stored, &fresh.amino_acids);
    }

    #[test]
    fn contaminant_pass_is_deterministic() {
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(ContaminantPass::new(0.5, Box::new(UniformBase))));
            p
        };
        let oa = PassRuntime::execute(&plan(), contaminant_test_sim(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&plan(), contaminant_test_sim(), 0xc0ff_ee);
        assert_eq!(oa.trace.choices(), ob.trace.choices());
        for i in 0..9u32 {
            assert_eq!(
                oa.final_simulation()
                    .pool
                    .get(NucHandle::new(i))
                    .unwrap()
                    .base,
                ob.final_simulation()
                    .pool
                    .get(NucHandle::new(i))
                    .unwrap()
                    .base
            );
        }
    }

    #[test]
    fn contaminant_pass_empty_pool_skips_replacement() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(ContaminantPass::new(1.0, Box::new(UniformBase))));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);
        // Coin flip happened (Bool recorded) but no bases written.
        assert_eq!(outcome.trace.len(), 1);
        assert_eq!(
            outcome
                .trace
                .find("corrupt.contaminant.applied")
                .unwrap()
                .value,
            ChoiceValue::Bool(true)
        );
    }

    #[test]
    fn contaminant_pass_declared_choices() {
        let pass = ContaminantPass::new(0.1, Box::new(UniformBase));
        let declared = pass.declared_choices();
        assert!(declared.contains(&"corrupt.contaminant.applied".to_string()));
        assert!(declared.contains(&"corrupt.contaminant.bases[0..n]".to_string()));
    }

    #[test]
    fn contaminant_productive_filters_base_that_would_create_stop() {
        let (cfg, sim) = make_contaminant_productive_vj_fixture();
        let seed = find_seed_for_unconstrained_contaminant_prefix(&sim, b"TAA");
        let contracts = productive();

        let constrained = PassRuntime::execute_with_context(
            &contaminant_plan(Box::new(UniformBase)),
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        );

        assert_eq!(
            constrained
                .trace
                .find("corrupt.contaminant.applied")
                .unwrap()
                .value,
            ChoiceValue::Bool(true)
        );
        assert_ne!(
            constrained
                .trace
                .find("corrupt.contaminant.bases[0]")
                .unwrap()
                .value,
            ChoiceValue::Base(b'T')
        );
        assert!(contracts
            .verify(constrained.final_simulation(), Some(&cfg))
            .is_ok());

        let unconstrained = PassRuntime::execute_with_context(
            &contaminant_plan(Box::new(UniformBase)),
            sim,
            seed,
            Some(&cfg),
            None,
        );
        let first_three: Vec<u8> = (0..3)
            .map(|i| {
                match unconstrained
                    .trace
                    .find(&format!("corrupt.contaminant.bases[{}]", i))
                    .unwrap()
                    .value
                {
                    ChoiceValue::Base(b) => b,
                    _ => unreachable!(),
                }
            })
            .collect();
        assert_eq!(first_three, b"TAA");
        let violations = productive()
            .verify(unconstrained.final_simulation(), Some(&cfg))
            .unwrap_err();
        assert!(violations
            .iter()
            .any(|v| v.contract_name == "no_stop_codon_in_junction"));
    }

    #[test]
    fn contaminant_strict_errors_when_base_filter_empty() {
        let (cfg, sim) = make_contaminant_productive_vj_fixture();
        let contracts = productive();

        let err = PassRuntime::execute_strict_with_context(
            &contaminant_plan(Box::new(StopOnlyContaminantBaseDist)),
            sim,
            0,
            Some(&cfg),
            Some(&contracts),
        )
        .unwrap_err();

        assert_eq!(err.pass_name(), "corrupt.contaminant");
        assert_eq!(err.address(), "corrupt.contaminant.bases[0]");
        assert_eq!(
            err.constraint_reason(),
            Some(FilteredSampleError::EmptyAdmissibleSupport)
        );
    }

    // ── QualityErrorPass tests (E.5) ───────────────────────────────

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
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region);
        // Original codon rail: ATG GGG GGG → M G G.
        assert_eq!(sim.sequence.regions[0].amino_acids, b"MGG");

        let mut plan = PassPlan::new();
        plan.push(Box::new(QualityErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, sim, 0);
        let final_sim = outcome.final_simulation();

        // After the pass, the codon rail still translates correctly
        // even though some bases are lowercase.
        let stored = &final_sim.sequence.regions[0].amino_acids;
        let fresh = final_sim.sequence.regions[0].with_codon_rail_recomputed(&final_sim.pool);
        assert_eq!(stored, &fresh.amino_acids);
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

    // ── PCRErrorPass tests (E.4) ───────────────────────────────────

    fn pcr_test_sim() -> Simulation {
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGGTTTAAACCCGGG".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(21))
            .with_codon_rail_recomputed(&sim.pool);
        sim.with_region_added(region)
    }

    fn pcr_single_error_plan(base_dist: Box<dyn Distribution<Output = u8>>) -> PassPlan {
        let mut plan = PassPlan::new();
        plan.push(Box::new(PCRErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
            base_dist,
        )));
        plan
    }

    fn pcr_error_site(seed: u64, sim: Simulation) -> u32 {
        let outcome = PassRuntime::execute(
            &pcr_single_error_plan(Box::new(StopOnlyMutationBaseDist)),
            sim,
            seed,
        );
        match outcome
            .trace
            .find("corrupt.pcr.error_site[0]")
            .unwrap()
            .value
        {
            ChoiceValue::Int(site) => site as u32,
            _ => panic!("wrong variant"),
        }
    }

    fn find_seed_for_pcr_error_site(sim: &Simulation, target_site: u32) -> u64 {
        for seed in 0..512u64 {
            if pcr_error_site(seed, sim.clone()) == target_site {
                return seed;
            }
        }
        panic!("no seed in search range targeted PCR site {}", target_site);
    }

    #[test]
    fn pcr_error_pass_zero_count_is_noop() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(PCRErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, pcr_test_sim(), 0);
        assert_eq!(outcome.trace.len(), 1);
        assert_eq!(
            outcome.trace.find("corrupt.pcr.count").unwrap().value,
            ChoiceValue::Int(0)
        );
    }

    #[test]
    fn pcr_error_pass_applies_n_errors_at_canonical_addresses() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(PCRErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, pcr_test_sim(), 7);

        // 1 count + 4 error_site + 4 error_base = 9 records.
        assert_eq!(outcome.trace.len(), 9);
        for i in 0..4 {
            assert!(outcome
                .trace
                .find(&format!("corrupt.pcr.error_site[{}]", i))
                .is_some());
            assert!(outcome
                .trace
                .find(&format!("corrupt.pcr.error_base[{}]", i))
                .is_some());
        }
    }

    #[test]
    fn pcr_error_pass_pool_reflects_recorded_errors() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(PCRErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, pcr_test_sim(), 99);
        let final_sim = outcome.final_simulation();

        let mut last_at_site: std::collections::HashMap<u32, u8> = std::collections::HashMap::new();
        for i in 0..3 {
            let s = match outcome
                .trace
                .find(&format!("corrupt.pcr.error_site[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Int(n) => n as u32,
                _ => unreachable!(),
            };
            let b = match outcome
                .trace
                .find(&format!("corrupt.pcr.error_base[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Base(b) => b,
                _ => unreachable!(),
            };
            last_at_site.insert(s, b);
        }
        for (&site, &expected) in last_at_site.iter() {
            assert_eq!(
                final_sim.pool.get(NucHandle::new(site)).unwrap().base,
                expected
            );
        }
    }

    #[test]
    fn pcr_error_pass_codon_rail_stays_consistent() {
        // Same post-D fix invariant: stored amino_acids matches
        // a fresh recompute against the post-error pool.
        let mut plan = PassPlan::new();
        plan.push(Box::new(PCRErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, pcr_test_sim(), 1);
        let final_sim = outcome.final_simulation();

        let stored = &final_sim.sequence.regions[0].amino_acids;
        let fresh = final_sim.sequence.regions[0].with_codon_rail_recomputed(&final_sim.pool);
        assert_eq!(stored, &fresh.amino_acids);
    }

    #[test]
    fn pcr_error_pass_is_deterministic() {
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(PCRErrorPass::new(
                Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
                Box::new(UniformBase),
            )));
            p
        };
        let oa = PassRuntime::execute(&plan(), pcr_test_sim(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&plan(), pcr_test_sim(), 0xc0ff_ee);
        assert_eq!(oa.trace.choices(), ob.trace.choices());
    }

    #[test]
    fn pcr_error_pass_declared_choices() {
        let pass = PCRErrorPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
            Box::new(UniformBase),
        );
        let declared = pass.declared_choices();
        assert!(declared.contains(&"corrupt.pcr.count".to_string()));
        assert!(declared.contains(&"corrupt.pcr.error_site[0..n]".to_string()));
        assert!(declared.contains(&"corrupt.pcr.error_base[0..n]".to_string()));
    }

    #[test]
    fn pcr_error_productive_filters_base_that_would_create_stop() {
        let (cfg, sim) = make_uniform_mutation_productive_vj_fixture();
        let seed = find_seed_for_pcr_error_site(&sim, 2);
        let contracts = productive();

        let constrained = PassRuntime::execute_with_context(
            &pcr_single_error_plan(Box::new(StopThenSafeMutationBaseDist)),
            sim.clone(),
            seed,
            Some(&cfg),
            Some(&contracts),
        );

        assert_eq!(
            constrained
                .trace
                .find("corrupt.pcr.error_site[0]")
                .unwrap()
                .value,
            ChoiceValue::Int(2)
        );
        assert_eq!(
            constrained
                .trace
                .find("corrupt.pcr.error_base[0]")
                .unwrap()
                .value,
            ChoiceValue::Base(b'C')
        );
        assert!(contracts
            .verify(constrained.final_simulation(), Some(&cfg))
            .is_ok());

        let unconstrained = PassRuntime::execute_with_context(
            &pcr_single_error_plan(Box::new(StopThenSafeMutationBaseDist)),
            sim,
            seed,
            Some(&cfg),
            None,
        );
        assert_eq!(
            unconstrained
                .trace
                .find("corrupt.pcr.error_base[0]")
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
    fn pcr_error_strict_errors_when_base_filter_empty() {
        let (cfg, sim) = make_uniform_mutation_productive_vj_fixture();
        let seed = find_seed_for_pcr_error_site(&sim, 2);
        let contracts = productive();

        let err = PassRuntime::execute_strict_with_context(
            &pcr_single_error_plan(Box::new(StopOnlyMutationBaseDist)),
            sim,
            seed,
            Some(&cfg),
            Some(&contracts),
        )
        .unwrap_err();

        assert_eq!(err.pass_name(), "corrupt.pcr");
        assert_eq!(err.address(), "corrupt.pcr.error_base[0]");
        assert_eq!(
            err.constraint_reason(),
            Some(FilteredSampleError::EmptyAdmissibleSupport)
        );
    }

    // ── IndelPass tests (E.7) ──────────────────────────────────────

    /// Helper: build a sim with a 12-base germline V region. Used as
    /// the canonical input for indel-pass tests so length/region
    /// invariants are easy to assert against the post-pass state.
    fn indel_test_sim() -> Simulation {
        let mut sim = Simulation::new();
        for (i, b) in b"AAACCCGGGTTT".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(12))
            .with_codon_rail_recomputed(&sim.pool);
        sim.with_region_added(region)
    }

    /// Helper: build a contiguous multi-segment simulation so indel
    /// insertion provenance can be checked outside V.
    fn multi_segment_indel_context_sim() -> Simulation {
        let layout = [
            (Segment::V, b"AAA".as_slice()),
            (Segment::Np1, b"CC".as_slice()),
            (Segment::J, b"GGG".as_slice()),
        ];

        let mut sim = Simulation::new();
        let mut regions = Vec::new();
        for (segment, bases) in layout {
            let start = sim.pool.len() as u32;
            for (i, b) in bases.iter().enumerate() {
                let (next, _) =
                    sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, segment));
                sim = next;
            }
            let end = sim.pool.len() as u32;
            regions.push(Region::new(
                segment,
                NucHandle::new(start),
                NucHandle::new(end),
            ));
        }

        for region in regions {
            sim = sim.with_region_added(region.with_codon_rail_recomputed(&sim.pool));
        }
        sim
    }

    /// Helper: a count distribution that always returns the given
    /// value. Built on `EmpiricalLengthDist` for ergonomic test setup.
    fn fixed_count(n: i64) -> Box<dyn Distribution<Output = i64>> {
        Box::new(EmpiricalLengthDist::from_pairs(vec![(n, 1.0)]))
    }

    #[test]
    #[should_panic(expected = "insertion_prob must be in [0.0, 1.0]")]
    fn indel_pass_rejects_negative_insertion_prob() {
        let _ = IndelPass::new(fixed_count(0), -0.5, Box::new(UniformBase));
    }

    #[test]
    #[should_panic(expected = "insertion_prob must be in [0.0, 1.0]")]
    fn indel_pass_rejects_insertion_prob_above_one() {
        let _ = IndelPass::new(fixed_count(0), 1.5, Box::new(UniformBase));
    }

    #[test]
    #[should_panic(expected = "insertion_prob must be in [0.0, 1.0]")]
    fn indel_pass_rejects_nan_insertion_prob() {
        let _ = IndelPass::new(fixed_count(0), f64::NAN, Box::new(UniformBase));
    }

    #[test]
    fn indel_pass_zero_count_is_noop() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(0),
            0.5,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 0);

        // Only the count is recorded; no per-indel addresses.
        assert_eq!(outcome.trace.len(), 1);
        assert_eq!(
            outcome.trace.find("corrupt.indel.count").unwrap().value,
            ChoiceValue::Int(0)
        );

        // Pool length unchanged.
        let final_sim = outcome.final_simulation();
        assert_eq!(final_sim.pool.len(), 12);
    }

    #[test]
    #[should_panic(expected = "count distribution returned negative")]
    fn indel_pass_negative_count_panics() {
        // crate::dist::UniformInt::new(-3, -2) always emits -3 —
        // exercises the negative-count guard.
        use crate::dist::UniformInt;
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            Box::new(UniformInt::new(-3, -2)),
            0.5,
            Box::new(UniformBase),
        )));
        let _ = PassRuntime::execute(&plan, indel_test_sim(), 0);
    }

    #[test]
    fn indel_pass_insertion_prob_one_grows_pool() {
        // With p_ins = 1.0, every event is an insertion → pool grows
        // by exactly count.
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(4),
            1.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 7);
        assert_eq!(outcome.final_simulation().pool.len(), 12 + 4);

        // Every kind[i] = Bool(true).
        for i in 0..4 {
            assert_eq!(
                outcome
                    .trace
                    .find(&format!("corrupt.indel.kind[{}]", i))
                    .unwrap()
                    .value,
                ChoiceValue::Bool(true)
            );
        }
    }

    #[test]
    fn indel_pass_insertion_prob_zero_shrinks_pool() {
        // With p_ins = 0.0, every event is a deletion → pool shrinks
        // by exactly count (provided pool stays non-empty).
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(4),
            0.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 11);
        assert_eq!(outcome.final_simulation().pool.len(), 12 - 4);

        for i in 0..4 {
            assert_eq!(
                outcome
                    .trace
                    .find(&format!("corrupt.indel.kind[{}]", i))
                    .unwrap()
                    .value,
                ChoiceValue::Bool(false)
            );
        }
    }

    #[test]
    fn indel_pass_inserted_nucleotides_carry_indel_flag() {
        // Every nucleotide that landed via insertion must have the
        // INDEL_INSERTED flag set; existing germline nucleotides
        // must not.
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(3),
            1.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 19);
        let final_sim = outcome.final_simulation();

        let mut inserted = 0;
        let mut germline = 0;
        for i in 0..final_sim.pool.len() as u32 {
            let n = final_sim.pool.get(NucHandle::new(i)).unwrap();
            if n.flags.contains(flag::INDEL_INSERTED) {
                inserted += 1;
                assert_eq!(n.germline_pos, Nucleotide::NO_GERMLINE_POS);
            } else {
                germline += 1;
            }
        }
        assert_eq!(inserted, 3);
        assert_eq!(germline, 12);
    }

    #[test]
    fn indel_pass_insertion_segment_uses_sequence_context() {
        let sim = multi_segment_indel_context_sim();

        assert_eq!(IndelPass::insertion_segment(&sim, 1), Segment::V);
        assert_eq!(IndelPass::insertion_segment(&sim, 3), Segment::Np1);
        assert_eq!(IndelPass::insertion_segment(&sim, 5), Segment::J);
        assert_eq!(IndelPass::insertion_segment(&sim, 8), Segment::J);
        assert_eq!(
            IndelPass::insertion_segment(&Simulation::new(), 0),
            Segment::V
        );
    }

    #[test]
    fn indel_pass_inserted_nucleotide_inherits_non_v_segment() {
        let sim = multi_segment_indel_context_sim();
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(1),
            1.0,
            Box::new(UniformBase),
        )));

        for seed in 0..200u64 {
            let outcome = PassRuntime::execute(&plan, sim.clone(), seed);
            let site = match outcome.trace.find("corrupt.indel.site[0]").unwrap().value {
                ChoiceValue::Int(s) => s as u32,
                _ => unreachable!(),
            };

            if site >= 5 {
                let inserted = outcome
                    .final_simulation()
                    .pool
                    .get(NucHandle::new(site))
                    .unwrap();
                assert_eq!(inserted.segment, Segment::J);
                assert!(inserted.flags.contains(flag::INDEL_INSERTED));
                return;
            }
        }

        panic!("test setup did not produce an insertion in J context");
    }

    #[test]
    fn indel_pass_records_canonical_addresses() {
        // 1 count + 3 kind[i] + 3 site[i] + 3 base[i] = 10 records
        // (since p_ins=1.0 every event is an insertion and base[i]
        // is always recorded).
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(3),
            1.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 5);
        assert_eq!(outcome.trace.len(), 10);
        for i in 0..3 {
            assert!(outcome
                .trace
                .find(&format!("corrupt.indel.kind[{}]", i))
                .is_some());
            assert!(outcome
                .trace
                .find(&format!("corrupt.indel.site[{}]", i))
                .is_some());
            assert!(outcome
                .trace
                .find(&format!("corrupt.indel.base[{}]", i))
                .is_some());
        }
    }

    #[test]
    fn indel_pass_deletion_path_does_not_record_base() {
        // p_ins = 0.0 → every event is a deletion → no base[i]
        // records are emitted.
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(3),
            0.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 5);
        // 1 count + 3 kind + 3 site = 7. No base entries.
        assert_eq!(outcome.trace.len(), 7);
        for i in 0..3 {
            assert!(outcome
                .trace
                .find(&format!("corrupt.indel.base[{}]", i))
                .is_none());
        }
    }

    #[test]
    fn indel_pass_inserted_base_in_pool_matches_traced_base() {
        // After insertion, the freshly-inserted nucleotide at the
        // recorded site carries the recorded base. Subsequent indels
        // can shift earlier positions, so we replay the events
        // forward and track each insert's final position.
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(2),
            1.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 42);
        let final_sim = outcome.final_simulation();

        // Replay the recorded events to know each insert's final
        // resting position after subsequent shifts.
        let mut sites: Vec<u32> = (0..2)
            .map(|i| {
                match outcome
                    .trace
                    .find(&format!("corrupt.indel.site[{}]", i))
                    .unwrap()
                    .value
                {
                    ChoiceValue::Int(s) => s as u32,
                    _ => unreachable!(),
                }
            })
            .collect();
        let bases: Vec<u8> = (0..2)
            .map(|i| {
                match outcome
                    .trace
                    .find(&format!("corrupt.indel.base[{}]", i))
                    .unwrap()
                    .value
                {
                    ChoiceValue::Base(b) => b,
                    _ => unreachable!(),
                }
            })
            .collect();

        // Each later insertion at a position ≤ an earlier site shifts
        // the earlier site up by 1. Walk the events left to right,
        // updating the resting positions of all previously-inserted
        // nucleotides.
        let mut resting: Vec<u32> = Vec::new();
        for i in 0..2 {
            let s = sites[i];
            for r in resting.iter_mut() {
                if *r >= s {
                    *r += 1;
                }
            }
            resting.push(s);
        }
        sites = resting;

        for (s, expected) in sites.iter().zip(bases.iter()) {
            let pool_base = final_sim.pool.get(NucHandle::new(*s)).unwrap().base;
            assert_eq!(pool_base, *expected, "trace lies at site {}", s);
        }
    }

    #[test]
    fn indel_pass_insertion_grows_region_only_when_spanning() {
        // The pre-pass region is [0, 12). An insertion at site < 12
        // spans the region and grows `end` by 1; an insertion at
        // site == pool_len lands strictly after the region and
        // leaves it unchanged. Replay the recorded events to compute
        // the expected final end and assert against it.
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(3),
            1.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 0);
        let final_sim = outcome.final_simulation();

        // Replay: track the moving region end and shift it whenever
        // an insert site falls strictly inside the current region.
        let mut region_end: u32 = 12;
        for i in 0..3 {
            let site = match outcome
                .trace
                .find(&format!("corrupt.indel.site[{}]", i))
                .unwrap()
                .value
            {
                ChoiceValue::Int(s) => s as u32,
                _ => unreachable!(),
            };
            // Spanning rule: region.start (0) ≤ site < region.end.
            if site < region_end {
                region_end += 1;
            }
            // site == region_end is the "entirely before" case for
            // the region (region.end ≤ pos), so no shift.
        }

        assert_eq!(final_sim.sequence.regions.len(), 1);
        assert_eq!(final_sim.sequence.regions[0].start.index(), 0);
        assert_eq!(final_sim.sequence.regions[0].end.index(), region_end);
        assert_eq!(final_sim.pool.len(), 12 + 3);
    }

    #[test]
    fn indel_pass_codon_rail_refresh_after_insertion() {
        // Stored amino_acids on the spanning region must equal a
        // fresh recompute against the post-pass pool. Same invariant
        // proved for substitution passes; indel-aware version exercises
        // the rail-refresh path through `with_indel_inserted`.
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(2),
            1.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 1);
        let final_sim = outcome.final_simulation();

        let stored = &final_sim.sequence.regions[0].amino_acids;
        let fresh = final_sim.sequence.regions[0].with_codon_rail_recomputed(&final_sim.pool);
        assert_eq!(stored, &fresh.amino_acids);
    }

    #[test]
    fn indel_pass_codon_rail_refresh_after_deletion() {
        // Same invariant for the deletion path.
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(2),
            0.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, indel_test_sim(), 1);
        let final_sim = outcome.final_simulation();

        let stored = &final_sim.sequence.regions[0].amino_acids;
        let fresh = final_sim.sequence.regions[0].with_codon_rail_recomputed(&final_sim.pool);
        assert_eq!(stored, &fresh.amino_acids);
    }

    #[test]
    fn indel_pass_is_deterministic_under_same_seed() {
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(IndelPass::new(
                fixed_count(3),
                0.5,
                Box::new(UniformBase),
            )));
            p
        };
        let oa = PassRuntime::execute(&plan(), indel_test_sim(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&plan(), indel_test_sim(), 0xc0ff_ee);
        assert_eq!(oa.trace.choices(), ob.trace.choices());
        assert_eq!(
            oa.final_simulation().pool.len(),
            ob.final_simulation().pool.len()
        );
        for i in 0..oa.final_simulation().pool.len() as u32 {
            let h = NucHandle::new(i);
            assert_eq!(
                oa.final_simulation().pool.get(h).unwrap().base,
                ob.final_simulation().pool.get(h).unwrap().base
            );
        }
    }

    #[test]
    fn indel_pass_empty_pool_deletion_records_minus_one_site() {
        // Deletion against an empty pool: no IR change, but the trace
        // records site = -1 to keep address layout structurally
        // consistent with the non-empty case.
        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(2),
            0.0,
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);

        // Pool stays empty.
        assert_eq!(outcome.final_simulation().pool.len(), 0);

        // Both site[i] entries are -1.
        for i in 0..2 {
            assert_eq!(
                outcome
                    .trace
                    .find(&format!("corrupt.indel.site[{}]", i))
                    .unwrap()
                    .value,
                ChoiceValue::Int(-1)
            );
        }
    }

    #[test]
    fn indel_pass_half_probability_produces_mixed_kinds() {
        // p_ins = 0.5 should produce both insertions and deletions
        // across many seeds. Negative-control proof that the kind
        // coin actually fires.
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(IndelPass::new(
                fixed_count(4),
                0.5,
                Box::new(UniformBase),
            )));
            p
        };

        let mut saw_insertion = false;
        let mut saw_deletion = false;
        for seed in 0..50u64 {
            let outcome = PassRuntime::execute(&plan(), indel_test_sim(), seed);
            for i in 0..4 {
                match outcome
                    .trace
                    .find(&format!("corrupt.indel.kind[{}]", i))
                    .unwrap()
                    .value
                {
                    ChoiceValue::Bool(true) => saw_insertion = true,
                    ChoiceValue::Bool(false) => saw_deletion = true,
                    _ => unreachable!(),
                }
            }
            if saw_insertion && saw_deletion {
                break;
            }
        }
        assert!(saw_insertion, "expected at least one insertion");
        assert!(saw_deletion, "expected at least one deletion");
    }

    #[test]
    fn indel_pass_persistent_ir_preserves_input() {
        // Persistent IR contract: input simulation isn't mutated.
        let sim = indel_test_sim();
        let original_len = sim.pool.len();
        let original_region_end = sim.sequence.regions[0].end.index();

        let mut plan = PassPlan::new();
        plan.push(Box::new(IndelPass::new(
            fixed_count(3),
            1.0,
            Box::new(UniformBase),
        )));
        let _outcome = PassRuntime::execute(&plan, sim.clone(), 99);

        assert_eq!(sim.pool.len(), original_len);
        assert_eq!(sim.sequence.regions[0].end.index(), original_region_end);
    }

    #[test]
    fn indel_pass_declared_choices() {
        let pass = IndelPass::new(fixed_count(0), 0.5, Box::new(UniformBase));
        let declared = pass.declared_choices();
        assert!(declared.contains(&"corrupt.indel.count".to_string()));
        assert!(declared.contains(&"corrupt.indel.kind[0..n]".to_string()));
        assert!(declared.contains(&"corrupt.indel.site[0..n]".to_string()));
        assert!(declared.contains(&"corrupt.indel.base[0..n]".to_string()));
    }

    // ── S5FMutationPass tests (E.3) ────────────────────────────────

    use crate::s5f::{S5FKernel, S5F_NUM_CONTEXTS, S5F_SUBSTITUTION_LEN};

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
        sub[0 * 4 + 3] = 1.0;
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

    #[test]
    fn uniform_mutation_pass_works_on_empty_pool() {
        // Edge case: no nucleotides yet → mutation is a no-op
        // regardless of the count distribution.
        let mut plan = PassPlan::new();
        plan.push(Box::new(UniformMutationPass::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
            Box::new(UniformBase),
        )));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);
        assert_eq!(outcome.final_simulation().pool.len(), 0);
        // Trace still records the count (sample happened) but no
        // site/base entries because the loop didn't iterate.
        assert_eq!(outcome.trace.len(), 1);
        assert_eq!(
            outcome.trace.find("mutate.uniform.count").unwrap().value,
            ChoiceValue::Int(5)
        );
    }

    // ── AssembleSegmentPass tests (C.8) ────────────────────────────

    use crate::refdata::{ChainType, RefDataConfig};

    /// Build a tiny reference config with one V, D, and J allele
    /// each of known length / sequence for assembly tests.
    fn make_test_refdata() -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.v_pool.push(Allele {
            name: "v_test*01".to_string(),
            gene: "v_test".to_string(),
            seq: b"AAACCCGGG".to_vec(), // 9 bases
            segment: Segment::V,
            anchor: Some(6), // C of CCC at position 6
        });
        let _ = cfg.d_pool.push(Allele {
            name: "d_test*01".to_string(),
            gene: "d_test".to_string(),
            seq: b"TTTTTT".to_vec(), // 6 bases
            segment: Segment::D,
            anchor: None,
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j_test*01".to_string(),
            gene: "j_test".to_string(),
            seq: b"GGGCCC".to_vec(), // 6 bases
            segment: Segment::J,
            anchor: Some(0),
        });
        cfg
    }

    #[test]
    #[should_panic(expected = "segment must be V, D, or J")]
    fn assemble_segment_pass_rejects_np1() {
        let _ = AssembleSegmentPass::new(Segment::Np1);
    }

    #[test]
    #[should_panic(expected = "PassContext.refdata is None")]
    fn assemble_segment_pass_panics_without_refdata() {
        // Use plain `execute` (refdata = None) — should panic.
        let mut plan = PassPlan::new();
        let v_pool = make_test_pool(1, Segment::V);
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        let _ = PassRuntime::execute(&plan, Simulation::new(), 0);
    }

    #[test]
    #[should_panic(expected = "no allele assigned")]
    fn assemble_segment_pass_panics_without_allele_assignment() {
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        let _ = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
    }

    #[test]
    #[should_panic(expected = "exceeds allele length")]
    fn assemble_segment_pass_panics_when_trim_exceeds_length() {
        let refdata = make_test_refdata();
        // V is 9 bases. Trim 5 + 3 = 8, leaves 1 base — fine.
        // Trim 5' = 5, 3' = 5 → sum 10 > 9 → panic.
        let v_pool = make_test_pool(1, Segment::V);
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        // Note: SampleAllelePass uses make_test_pool, but AssembleSegmentPass
        // reads from refdata. They MUST be the same allele for the
        // test to succeed; we rely on test_refdata having allele 0
        // be the V whose length is 9.
        let _ = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
    }

    #[test]
    fn assemble_segment_pass_v_no_trim_pushes_full_allele() {
        let refdata = make_test_refdata();
        // Use the V pool from refdata directly so SampleAllele sees
        // the same alleles as Assemble. We sample id 0 (the only V).
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
        let sim = outcome.final_simulation();

        // V allele was 9 bases; no trim → 9 nucleotides in pool.
        assert_eq!(sim.pool.len(), 9);
        // One region (V) of length 9.
        assert_eq!(sim.sequence.region_count(), 1);
        let r = &sim.sequence.regions[0];
        assert_eq!(r.segment, Segment::V);
        assert_eq!(r.len(), 9);
        assert_eq!(r.frame_phase, 0); // first region

        // Bases match the allele sequence exactly, with germline_pos
        // matching position-in-allele.
        let expected = b"AAACCCGGG";
        for i in 0..9 {
            let n = sim.pool.get(NucHandle::new(i)).unwrap();
            assert_eq!(n.base, expected[i as usize]);
            assert_eq!(n.germline, expected[i as usize]);
            assert_eq!(n.germline_pos, i as u16);
            assert_eq!(n.segment, Segment::V);
        }
    }

    #[test]
    fn assemble_segment_pass_v_with_trim_pushes_post_trim_slice() {
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
        let sim = outcome.final_simulation();

        // Allele = "AAACCCGGG" (9 bases); trim_5=2, trim_3=3.
        // Slice = "ACCCG" (positions 2..6, length 4 — wait let me recompute).
        // Actually 9 - 2 - 3 = 4 → slice is positions [2, 6) = "ACCC".
        assert_eq!(sim.pool.len(), 4);
        let expected = b"ACCC";
        for i in 0..4 {
            let n = sim.pool.get(NucHandle::new(i)).unwrap();
            assert_eq!(n.base, expected[i as usize]);
            assert_eq!(n.germline_pos, (i + 2) as u16); // shifted by trim_5
        }
    }

    #[test]
    fn assemble_segment_pass_chains_frame_phase_across_regions() {
        // V (9 bases, frame_phase=0) + assemble J (6 bases,
        // frame_phase = 9 % 3 = 0). So both regions have phase 0.
        // Test that the chain is correct after assembly.
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&refdata.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
        let sim = outcome.final_simulation();

        assert_eq!(sim.pool.len(), 15); // 9 V + 6 J
        assert_eq!(sim.sequence.region_count(), 2);
        assert_eq!(sim.sequence.regions[0].segment, Segment::V);
        assert_eq!(sim.sequence.regions[0].frame_phase, 0);
        assert_eq!(sim.sequence.regions[1].segment, Segment::J);
        // J starts after 9 V bases; 9 % 3 = 0 → frame_phase 0.
        assert_eq!(sim.sequence.regions[1].frame_phase, 0);
    }

    #[test]
    fn assemble_segment_pass_chains_frame_phase_with_offset() {
        // Build a plan where V has 2 bases (after heavy trim), so
        // J's frame_phase = 2 % 3 = 2. Tests the non-zero phase chain.
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&refdata.j_pool)),
        )));
        // V is 9 bases. Trim 5'=4, 3'=3 → slice has 9-4-3 = 2 bases.
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
        let sim = outcome.final_simulation();

        assert_eq!(sim.sequence.regions[0].len(), 2);
        assert_eq!(sim.sequence.regions[0].frame_phase, 0);
        // J's first base sits at position 2 of the cumulative chain;
        // 2 % 3 = 2 → frame_phase 2.
        assert_eq!(sim.sequence.regions[1].frame_phase, 2);
    }

    #[test]
    fn assemble_segment_pass_does_not_record_to_trace() {
        // Determinism contract for transform passes: no trace
        // entries.
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 7, &refdata);

        // Trace contains only the SampleAllele entry, none from Assemble.
        assert_eq!(outcome.trace.len(), 1);
        assert_eq!(outcome.trace.choices()[0].address, "sample_allele.v");
    }

    #[test]
    fn assemble_segment_pass_declared_choices_is_empty() {
        let pass = AssembleSegmentPass::new(Segment::V);
        assert!(pass.declared_choices().is_empty());
    }

    #[test]
    fn assemble_segment_pass_codon_rail_computed_on_assembled_region() {
        // V allele "AAACCCGGG" (9 bases), no trim, frame_phase 0
        // → codons AAA, CCC, GGG → K, P, G.
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
        let sim = outcome.final_simulation();
        let r = &sim.sequence.regions[0];

        assert_eq!(r.amino_acids, b"KPG");
        assert!(r.stop_codon_positions.is_empty());
    }

    // ── GenerateNPPass tests (C.7) ─────────────────────────────────

    // (UniformBase already imported earlier in this test module.)

    #[test]
    #[should_panic(expected = "np_segment must be Np1 or Np2")]
    fn generate_np_pass_rejects_v_segment() {
        let _ = GenerateNPPass::new(
            Segment::V,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        );
    }

    #[test]
    #[should_panic(expected = "np_segment must be Np1 or Np2")]
    fn generate_np_pass_rejects_d_segment() {
        let _ = GenerateNPPass::new(
            Segment::D,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        );
    }

    #[test]
    fn generate_np_pass_zero_length_creates_empty_region() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);

        let sim = outcome.final_simulation();
        // Pool is unchanged.
        assert_eq!(sim.pool.len(), 0);
        // One region was added (NP1, empty).
        assert_eq!(sim.sequence.region_count(), 1);
        let r = &sim.sequence.regions[0];
        assert_eq!(r.segment, Segment::Np1);
        assert!(r.is_empty());
        // Trace recorded length, no bases.
        assert_eq!(outcome.trace.len(), 1);
        assert_eq!(
            outcome.trace.find("np.np1.length").unwrap().value,
            ChoiceValue::Int(0)
        );
    }

    #[test]
    fn generate_np_pass_pushes_n_bases_with_correct_metadata() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 1234);

        let sim = outcome.final_simulation();
        assert_eq!(sim.pool.len(), 5);
        assert_eq!(sim.sequence.region_count(), 1);

        // Region covers all 5 nucleotides.
        let r = &sim.sequence.regions[0];
        assert_eq!(r.segment, Segment::Np1);
        assert_eq!(r.start, NucHandle::new(0));
        assert_eq!(r.end, NucHandle::new(5));

        // Each pushed nucleotide is a synthetic NP1 N-nuc.
        for i in 0..5 {
            let n = sim.pool.get(NucHandle::new(i)).unwrap();
            assert_eq!(n.segment, Segment::Np1);
            assert_eq!(n.germline_pos, Nucleotide::NO_GERMLINE_POS);
            assert!(n.flags.contains(flag::N_NUC));
            assert!(matches!(n.base, b'A' | b'C' | b'G' | b'T'));
        }
    }

    #[test]
    fn generate_np_pass_records_length_and_indexed_bases() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np2,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 7);

        // Length entry + 3 base entries = 4 trace records.
        assert_eq!(outcome.trace.len(), 4);
        assert_eq!(
            outcome.trace.find("np.np2.length").unwrap().value,
            ChoiceValue::Int(3)
        );
        for i in 0..3 {
            let addr = format!("np.np2.bases[{}]", i);
            let rec = outcome
                .trace
                .find(&addr)
                .unwrap_or_else(|| panic!("missing {}", addr));
            assert!(matches!(rec.value, ChoiceValue::Base(_)));
        }
    }

    #[test]
    fn generate_np_pass_trace_value_matches_pushed_base_at_each_position() {
        // Faithfulness: trace[i] equals pool[i].base for each NP base.
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(8, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 99);
        let sim = outcome.final_simulation();

        for i in 0..8 {
            let addr = format!("np.np1.bases[{}]", i);
            let recorded = match outcome.trace.find(&addr).unwrap().value {
                ChoiceValue::Base(b) => b,
                _ => panic!("wrong variant"),
            };
            let pushed = sim.pool.get(NucHandle::new(i)).unwrap().base;
            assert_eq!(
                recorded, pushed,
                "trace at {} = {} but pool[{}] = {}",
                addr, recorded as char, i, pushed as char
            );
        }
    }

    #[test]
    fn generate_np_pass_is_deterministic_under_same_seed() {
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(GenerateNPPass::new(
                Segment::Np1,
                Box::new(EmpiricalLengthDist::from_pairs(vec![
                    (1, 1.0),
                    (3, 2.0),
                    (5, 1.0),
                ])),
                Box::new(UniformBase),
            )));
            p
        };
        let oa = PassRuntime::execute(&plan(), Simulation::new(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&plan(), Simulation::new(), 0xc0ff_ee);

        assert_eq!(
            oa.final_simulation().pool.len(),
            ob.final_simulation().pool.len()
        );
        assert_eq!(oa.trace.choices(), ob.trace.choices());
    }

    #[test]
    fn generate_np_pass_declared_choices_lists_length_and_bases_pattern() {
        let pass = GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        );
        let declared = pass.declared_choices();
        assert_eq!(declared.len(), 2);
        assert!(declared.contains(&"np.np1.length".to_string()));
        assert!(declared.contains(&"np.np1.bases[0..n]".to_string()));
    }

    #[test]
    fn generate_np_pass_appends_region_after_existing_pool_state() {
        // Pre-populate the pool with some nucleotides, then run NP
        // generation. The new region should start at the current
        // pool boundary, not at handle 0.
        let mut plan = PassPlan::new();
        // Three Echo passes push 3 nucleotides as germline.
        plan.push(Box::new(EchoPass::new(b'A', 0, Segment::V)));
        plan.push(Box::new(EchoPass::new(b'C', 1, Segment::V)));
        plan.push(Box::new(EchoPass::new(b'G', 2, Segment::V)));
        // Then NP1 of length 4.
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
            Box::new(UniformBase),
        )));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 5);
        let sim = outcome.final_simulation();

        // Total pool: 3 V + 4 NP1 = 7.
        assert_eq!(sim.pool.len(), 7);
        // One region (NP1) — Echo doesn't create regions.
        assert_eq!(sim.sequence.region_count(), 1);
        let r = &sim.sequence.regions[0];
        assert_eq!(r.start, NucHandle::new(3));
        assert_eq!(r.end, NucHandle::new(7));
        assert_eq!(r.segment, Segment::Np1);
    }

    #[test]
    #[should_panic(expected = "negative value")]
    fn generate_np_pass_panics_on_negative_length() {
        use crate::dist::UniformInt;
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(UniformInt::new(-3, -1)),
            Box::new(UniformBase),
        )));
        let _ = PassRuntime::execute(&plan, Simulation::new(), 0);
    }

    // ── Constraint-aware sampling tests (D.6) ──────────────────────

    use crate::contract::{
        productive, Contract, ContractSet, ContractViolation, ProductiveJunctionFrame,
    };

    /// Build a synthetic VJ refdata: V "AAACCCGGG" (9bp, anchor 6),
    /// J "TTTAAA" (6bp, anchor 0). Junction = V_anchor_to_end (3bp)
    /// + NP1 + J_anchor_to_W3 (3bp). For productive frame:
    /// (3 + NP1 + 3) % 3 == 0 → NP1 % 3 == 0.
    fn make_vj_refdata_for_filter() -> crate::refdata::RefDataConfig {
        let mut cfg = crate::refdata::RefDataConfig::empty(crate::refdata::ChainType::Vj);
        let _ = cfg.v_pool.push(crate::refdata::Allele {
            name: "v_test*01".into(),
            gene: "v_test".into(),
            seq: b"AAACCCGGG".to_vec(),
            segment: Segment::V,
            anchor: Some(6),
        });
        let _ = cfg.j_pool.push(crate::refdata::Allele {
            name: "j_test*01".into(),
            gene: "j_test".into(),
            seq: b"TTTAAA".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });
        cfg
    }

    #[test]
    fn productive_admits_filters_np1_for_vj_chain() {
        // Length distribution offers 0..6 (mix of in-frame and
        // out-of-frame). Constraint should filter to {0, 3} since
        // V_anchor_to_end is 3 (from anchor 6 to end of 9bp V),
        // and J anchor offset + 3 contributes 3 → fixed 6.
        // junction_length = 6 + NP1, so NP1 % 3 == 0.
        let cfg = make_vj_refdata_for_filter();
        let v_pool = &cfg.v_pool;
        let j_pool = &cfg.j_pool;
        let dist = EmpiricalLengthDist::from_pairs((0..7).map(|i| (i, 1.0)).collect::<Vec<_>>());

        // Build a plan: SampleAlleleV, SampleAlleleJ, AssembleV,
        // GenerateNP1 (constrained), then we inspect.
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(dist),
            Box::new(UniformBase),
        )));

        let contracts = ContractSet::new().with(Box::new(ProductiveJunctionFrame::new()));

        // Run multiple seeds; every NP1 sample must be 0, 3, or 6.
        for seed in 0..20u64 {
            let outcome = PassRuntime::execute_with_context(
                &plan,
                Simulation::new(),
                seed,
                Some(&cfg),
                Some(&contracts),
            );
            let np1_len = match outcome.trace.find("np.np1.length").unwrap().value {
                ChoiceValue::Int(n) => n,
                _ => panic!("wrong variant"),
            };
            assert!(
                np1_len % 3 == 0,
                "seed {} produced NP1 length {} (not divisible by 3)",
                seed,
                np1_len
            );
        }
    }

    #[test]
    fn productive_admits_unconstrained_without_contracts() {
        // Same plan, no contracts → NP1 length distributed as the
        // full empirical dist (we just verify lengths are in [0, 6]).
        let cfg = make_vj_refdata_for_filter();
        let dist = EmpiricalLengthDist::from_pairs((0..7).map(|i| (i, 1.0)).collect::<Vec<_>>());

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(dist),
            Box::new(UniformBase),
        )));

        // Without contracts: any value in [0, 6] is allowed.
        let mut seen_out_of_frame = false;
        for seed in 0..50u64 {
            let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), seed, &cfg);
            let np1_len = match outcome.trace.find("np.np1.length").unwrap().value {
                ChoiceValue::Int(n) => n,
                _ => panic!("wrong variant"),
            };
            assert!(np1_len >= 0 && np1_len <= 6);
            if np1_len % 3 != 0 {
                seen_out_of_frame = true;
            }
        }
        // With 50 seeds and 4/7 lengths out of frame, we should see
        // out-of-frame at least once. If not, RNG is suspect.
        assert!(
            seen_out_of_frame,
            "Without contracts, expected at least one out-of-frame NP1 sample"
        );
    }

    #[test]
    fn productive_admits_falls_back_when_filter_empty() {
        // Distribution offers ONLY out-of-frame values {1, 2, 4, 5}.
        // Filter is empty → fallback to unconstrained sampling. We
        // verify the pass doesn't panic and produces SOME value.
        let cfg = make_vj_refdata_for_filter();
        let dist = EmpiricalLengthDist::from_pairs(vec![(1, 1.0), (2, 1.0), (4, 1.0), (5, 1.0)]);

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(dist),
            Box::new(UniformBase),
        )));

        let contracts = ContractSet::new().with(Box::new(ProductiveJunctionFrame::new()));

        // Should not panic; produces some value from the original
        // (unfiltered) support.
        let outcome = PassRuntime::execute_with_context(
            &plan,
            Simulation::new(),
            0,
            Some(&cfg),
            Some(&contracts),
        );
        let np1_len = match outcome.trace.find("np.np1.length").unwrap().value {
            ChoiceValue::Int(n) => n,
            _ => panic!("wrong variant"),
        };
        assert!(matches!(np1_len, 1 | 2 | 4 | 5));
    }

    #[derive(Clone, Debug)]
    struct UnenumerableLengthDist;

    impl Distribution for UnenumerableLengthDist {
        type Output = i64;

        fn sample(&self, _rng: &mut crate::rng::Rng) -> i64 {
            0
        }
    }

    #[test]
    fn productive_strict_errors_when_length_filter_empty() {
        let cfg = make_vj_refdata_for_filter();
        let dist = EmpiricalLengthDist::from_pairs(vec![(1, 1.0), (2, 1.0), (4, 1.0), (5, 1.0)]);

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(dist),
            Box::new(UniformBase),
        )));

        let contracts = ContractSet::new().with(Box::new(ProductiveJunctionFrame::new()));
        let err = PassRuntime::execute_strict_with_context(
            &plan,
            Simulation::new(),
            0,
            Some(&cfg),
            Some(&contracts),
        )
        .unwrap_err();

        assert_eq!(err.pass_name(), "generate_np.np1");
        assert_eq!(err.address(), "np.np1.length");
        assert_eq!(
            err.constraint_reason(),
            Some(FilteredSampleError::EmptyAdmissibleSupport)
        );
        assert!(err.to_string().contains("no admissible candidates"));
    }

    #[test]
    fn productive_strict_errors_when_length_support_unavailable() {
        let cfg = make_vj_refdata_for_filter();

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(UnenumerableLengthDist),
            Box::new(UniformBase),
        )));

        let contracts = ContractSet::new().with(Box::new(ProductiveJunctionFrame::new()));
        let err = PassRuntime::execute_strict_with_context(
            &plan,
            Simulation::new(),
            0,
            Some(&cfg),
            Some(&contracts),
        )
        .unwrap_err();

        assert_eq!(err.pass_name(), "generate_np.np1");
        assert_eq!(err.address(), "np.np1.length");
        assert_eq!(
            err.constraint_reason(),
            Some(FilteredSampleError::SupportUnavailable)
        );
    }

    struct RejectNpBases;

    impl Contract for RejectNpBases {
        fn name(&self) -> &str {
            "reject_np_bases"
        }

        fn verify(
            &self,
            _sim: &Simulation,
            _refdata: Option<&crate::refdata::RefDataConfig>,
        ) -> Result<(), ContractViolation> {
            Ok(())
        }

        fn admits(
            &self,
            _sim: &Simulation,
            _refdata: Option<&crate::refdata::RefDataConfig>,
            address: &str,
            _candidate: &ChoiceValue,
        ) -> Result<(), ContractViolation> {
            if address.starts_with("np.np1.bases[") {
                return Err(ContractViolation::new(self.name(), "rejected by test"));
            }
            Ok(())
        }
    }

    #[test]
    fn productive_strict_errors_when_base_filter_empty() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
            Box::new(UniformBase),
        )));
        let contracts = ContractSet::new().with(Box::new(RejectNpBases));

        let err = PassRuntime::execute_strict_with_context(
            &plan,
            Simulation::new(),
            0,
            None,
            Some(&contracts),
        )
        .unwrap_err();

        assert_eq!(err.pass_name(), "generate_np.np1");
        assert_eq!(err.address(), "np.np1.bases[0]");
        assert_eq!(
            err.constraint_reason(),
            Some(FilteredSampleError::EmptyAdmissibleSupport)
        );
    }

    #[test]
    fn productive_strict_succeeds_when_admissible_length_exists() {
        let cfg = make_vj_refdata_for_filter();
        let dist = EmpiricalLengthDist::from_pairs((0..7).map(|i| (i, 1.0)).collect::<Vec<_>>());

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(dist),
            Box::new(UniformBase),
        )));

        let contracts = ContractSet::new().with(Box::new(ProductiveJunctionFrame::new()));
        let outcome = PassRuntime::execute_strict_with_context(
            &plan,
            Simulation::new(),
            0,
            Some(&cfg),
            Some(&contracts),
        )
        .unwrap();
        let np1_len = match outcome.trace.find("np.np1.length").unwrap().value {
            ChoiceValue::Int(n) => n,
            _ => panic!("wrong variant"),
        };

        assert_eq!(np1_len % 3, 0);
    }

    #[test]
    fn productive_admits_makes_full_pipeline_in_frame_for_vj() {
        // The architectural payoff: PRODUCTIVE_ONLY VJ produces
        // in-frame junctions BY CONSTRUCTION, with no retries needed.
        // Run many seeds; verify junction is in-frame every time.
        use crate::junction::compute_junction;

        let cfg = make_vj_refdata_for_filter();
        let dist = EmpiricalLengthDist::from_pairs((0..7).map(|i| (i, 1.0)).collect::<Vec<_>>());

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(dist),
            Box::new(UniformBase),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

        let contracts = ContractSet::new().with(Box::new(ProductiveJunctionFrame::new()));

        for seed in 0..30u64 {
            let outcome = PassRuntime::execute_with_context(
                &plan,
                Simulation::new(),
                seed,
                Some(&cfg),
                Some(&contracts),
            );
            let junction = compute_junction(outcome.final_simulation(), &cfg)
                .expect("junction should be defined");
            assert!(
                junction.is_in_frame(),
                "seed {} produced out-of-frame junction (length {})",
                seed,
                junction.length
            );
        }
    }

    #[test]
    fn productive_full_bundle_in_frame_and_admits_dispatch_works() {
        // productive() bundle includes ProductiveJunctionFrame with
        // its admits override. Verify the bundle works end to end.
        use crate::junction::compute_junction;

        let cfg = make_vj_refdata_for_filter();
        let dist = EmpiricalLengthDist::from_pairs((0..10).map(|i| (i, 1.0)).collect::<Vec<_>>());

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&cfg.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(dist),
            Box::new(UniformBase),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

        let contracts = productive();

        for seed in 0..20u64 {
            let outcome = PassRuntime::execute_with_context(
                &plan,
                Simulation::new(),
                seed,
                Some(&cfg),
                Some(&contracts),
            );
            let junction = compute_junction(outcome.final_simulation(), &cfg)
                .expect("junction should be defined");
            assert!(junction.is_in_frame());
        }
    }

    #[test]
    fn sample_allele_pass_replay_in_mixed_plan_with_echo_passes() {
        // Mixed plan: SampleAllele + Echo (transform) interleaved.
        // The trace records only the sampling choices; the IR
        // accumulates both the assignment and the echoed nucleotides.
        let pool = make_test_pool(3, Segment::V);
        let mut plan = PassPlan::new();
        plan.push(Box::new(EchoPass::new(b'A', 0, Segment::V)));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&pool)),
        )));
        plan.push(Box::new(EchoPass::new(b'C', 1, Segment::V)));

        let oa = PassRuntime::execute(&plan, Simulation::new(), 0xfeed);
        let ob = PassRuntime::execute(&plan, Simulation::new(), 0xfeed);

        // Determinism end to end.
        assert_eq!(
            oa.final_simulation().pool.len(),
            ob.final_simulation().pool.len()
        );
        assert_eq!(
            oa.final_simulation().assignments.v.unwrap().allele_id,
            ob.final_simulation().assignments.v.unwrap().allele_id
        );
        // Trace contains exactly one entry — only the sampling pass records.
        assert_eq!(oa.trace.len(), 1);
        assert_eq!(oa.trace.choices()[0].address, "sample_allele.v");
    }

    #[test]
    fn replay_determinism_holds_under_repeated_runs() {
        // Run the same plan with the same seed five times. Every
        // trace and every final IR must be identical to the first.
        let baseline = PassRuntime::execute(&mixed_plan(), Simulation::new(), 0xfeed_face);
        for _ in 0..5 {
            let other = PassRuntime::execute(&mixed_plan(), Simulation::new(), 0xfeed_face);
            assert_eq!(other.trace.len(), baseline.trace.len());
            for (a, b) in baseline
                .trace
                .choices()
                .iter()
                .zip(other.trace.choices().iter())
            {
                assert_eq!(a.address, b.address);
                assert_eq!(a.value, b.value);
            }
            assert_eq!(
                other.final_simulation().pool.len(),
                baseline.final_simulation().pool.len()
            );
        }
    }
}
