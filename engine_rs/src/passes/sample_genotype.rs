//! `SampleGenotypePass` — phased, genotype-aware V(D)J allele sampling.
//!
//! Replaces the three flat `SampleAllelePass` passes when a genotype is
//! attached. In ONE pass it:
//!   1. draws the chromosome (haplotype) once, restricted to **viable**
//!      haplotypes (those carrying a feasible allele for every required
//!      segment under the active contracts + feasibility), and
//!   2. for each of V/(D)/J, samples a gene present on that chromosome
//!      (usage-weighted) then the allele within the gene slot
//!      (single-copy slots are deterministic), assigning the slot.
//!
//! The chromosome is a local variable — there is no cross-pass state to
//! share (each pass gets a fresh per-pass trace), which is exactly why
//! haplotype + gene + allele all live in one pass. The canonical
//! `sample_allele.{seg}` choice + the slot assignment are emitted just
//! like the flat path, so every downstream pass (assemble/trim/AIRR/
//! replay) is unchanged. Gene + within-slot choices are recorded as
//! extra addresses for provenance and replay.
use std::sync::Arc;

use crate::address::{self, ChoiceAddress, ChoiceAddressPattern};
use crate::assignment::AlleleInstance;
use crate::contract::ChoiceContext;
use crate::dist::FilteredSampleError;
use crate::genotype::Genotype;
use crate::ir::{Segment, Simulation, SimulationBuilder};
use crate::pass::{AlleleIdSupport, Pass, PassCompileFact, PassContext, PassEffect, PassError};
use crate::refdata::{AlleleId, GeneId};
use crate::rng::Rng;
use crate::trace::ChoiceValue;

pub struct SampleGenotypePass {
    genotype: Arc<Genotype>,
    live_call_support: Arc<crate::live_call::LiveCallSupport>,
    d_required: bool,
    // Per-segment gene-usage weights (empty => uniform over present genes).
    usage_v: Vec<(GeneId, f64)>,
    usage_d: Vec<(GeneId, f64)>,
    usage_j: Vec<(GeneId, f64)>,
}

impl SampleGenotypePass {
    pub fn new(
        genotype: Arc<Genotype>,
        d_required: bool,
        usage_v: Vec<(GeneId, f64)>,
        usage_d: Vec<(GeneId, f64)>,
        usage_j: Vec<(GeneId, f64)>,
    ) -> Self {
        let live_call_support =
            Arc::new(crate::live_call::LiveCallSupport::from_genotype(&genotype));
        Self {
            genotype,
            live_call_support,
            d_required,
            usage_v,
            usage_d,
            usage_j,
        }
    }

    fn segments(&self) -> Vec<Segment> {
        if self.d_required {
            vec![Segment::V, Segment::D, Segment::J]
        } else {
            vec![Segment::V, Segment::J]
        }
    }

    fn usage(&self, seg: Segment) -> &[(GeneId, f64)] {
        match seg {
            Segment::V => &self.usage_v,
            Segment::D => &self.usage_d,
            Segment::J => &self.usage_j,
            _ => &[],
        }
    }

    fn usage_of(&self, seg: Segment, g: GeneId) -> f64 {
        let table = self.usage(seg);
        if table.is_empty() {
            1.0
        } else {
            table
                .iter()
                .find(|(gg, _)| *gg == g)
                .map(|(_, w)| *w)
                .unwrap_or(1.0)
        }
    }

    fn vseg(seg: Segment) -> address::VdjSegment {
        seg.try_into().expect("V/D/J segment")
    }

    fn allele_feasible(
        &self,
        c: usize,
        seg: Segment,
        id: AlleleId,
        sim: &Simulation,
        ctx: &PassContext,
    ) -> bool {
        let choice = ChoiceValue::AlleleId(id.index());
        let vseg = Self::vseg(seg);
        let contract_ok = ctx.contracts.map_or(true, |k| {
            k.admits_typed(
                sim,
                ctx.refdata,
                ChoiceContext::none().with_address(ChoiceAddress::SampleAllele(vseg)),
                &choice,
            )
            .is_ok()
        });
        let feasible_ok = ctx.feasibility.map_or(true, |f| {
            if f.has_p_productive() {
                // The compiled marginal support combines both chromosomes.
                // Productive completion must instead use future alleles that
                // can actually be sampled on the chromosome drawn here.
                let allowed = self.positive_allele_support(c);
                return allowed
                    .iter()
                    .any(|(s, ids)| *s == seg && ids.contains(&id))
                    && f.admits_with_allele_restrictions(
                        ctx.pass_index,
                        sim,
                        ctx.refdata,
                        address::sample_allele_vdj(seg),
                        &choice,
                        &allowed,
                    );
            }
            f.admits(
                ctx.pass_index,
                sim,
                ctx.refdata,
                address::sample_allele_vdj(seg),
                &choice,
            )
        });
        contract_ok && feasible_ok
    }

    fn positive_allele_support(&self, c: usize) -> Vec<(Segment, Vec<AlleleId>)> {
        let hap = self.genotype.haplotype(c);
        self.segments()
            .into_iter()
            .map(|seg| {
                let mut ids = Vec::new();
                for gene in hap.present_genes(seg) {
                    if self.usage_of(seg, gene) <= 0.0 {
                        continue;
                    }
                    for copy in hap.slot(seg, gene) {
                        if copy.copies > 0 && copy.weight > 0.0 && !ids.contains(&copy.allele) {
                            ids.push(copy.allele);
                        }
                    }
                }
                (seg, ids)
            })
            .collect()
    }

    /// Carried alleles on chromosome `c` for `seg`, optionally filtered
    /// by the active contracts + feasibility.
    fn admissible_alleles(
        &self,
        c: usize,
        seg: Segment,
        sim: &Simulation,
        ctx: &PassContext,
        filter: bool,
    ) -> Vec<AlleleId> {
        let carried = self.genotype.haplotype(c).carried_alleles(seg);
        if !filter {
            return carried;
        }
        carried
            .into_iter()
            .filter(|id| self.allele_feasible(c, seg, *id, sim, ctx))
            .collect()
    }

    fn is_viable(&self, c: usize, sim: &Simulation, ctx: &PassContext, filter: bool) -> bool {
        self.segments().iter().all(|seg| {
            !self
                .admissible_alleles(c, *seg, sim, ctx, filter)
                .is_empty()
        })
    }

    /// Viable chromosomes. The feasibility-filtered viable set is
    /// preferred; in permissive mode, only if NO chromosome is
    /// feasibility-viable do we fall back to presence-viability — exactly
    /// mirroring `SampleAllelePass` (filter first, fall back to the
    /// unconstrained draw only on empty admissible support). In strict
    /// mode the feasibility-viable set is authoritative (empty → the
    /// caller raises).
    fn viable_set(&self, sim: &Simulation, ctx: &PassContext, strict: bool) -> Vec<usize> {
        // A zero-weight chromosome is never expressed, so it can never make the
        // genotype viable — exclude it from the candidate set entirely (rather
        // than letting draw_haplotype fall back to a zero-weight chromosome).
        let w = self.genotype.chromosome_weights();
        let expressible = |c: usize| w[c] as f64 > 0.0;
        let feasible: Vec<usize> = (0..2)
            .filter(|&c| expressible(c) && self.is_viable(c, sim, ctx, true))
            .collect();
        if strict || !feasible.is_empty() {
            return feasible;
        }
        (0..2)
            .filter(|&c| expressible(c) && self.is_viable(c, sim, ctx, false))
            .collect()
    }

    fn draw_haplotype(&self, viable: &[usize], rng: &mut Rng) -> usize {
        let w = self.genotype.chromosome_weights();
        let total: f64 = viable.iter().map(|&c| w[c] as f64).sum();
        if total <= 0.0 {
            return *viable.first().expect("viable non-empty");
        }
        let mut x = rng.next_f64() * total;
        for &c in viable {
            x -= w[c] as f64;
            if x < 0.0 {
                return c;
            }
        }
        *viable.last().expect("viable non-empty")
    }

    fn weighted_pick<T: Copy>(items: &[(T, f64)], rng: &mut Rng) -> T {
        let total: f64 = items.iter().map(|(_, w)| *w).sum();
        if total <= 0.0 {
            return items.first().expect("non-empty").0;
        }
        let mut x = rng.next_f64() * total;
        for (t, w) in items {
            x -= *w;
            if x < 0.0 {
                return *t;
            }
        }
        items.last().expect("non-empty").0
    }

    fn infeasible_error(&self) -> PassError {
        PassError::constraint_sampling(
            self.name(),
            "sample_haplotype",
            FilteredSampleError::EmptyAdmissibleSupport,
        )
    }

    fn commit(
        &self,
        seg: Segment,
        sim: Simulation,
        id: AlleleId,
        c: usize,
        ctx: &mut PassContext,
    ) -> Simulation {
        let mut b = SimulationBuilder::from_simulation(sim);
        if ctx.event_log_sink.is_some() {
            b.attach_event_log_observer();
        }
        b.assign_allele(seg, AlleleInstance::new(id).with_haplotype(c as u8));
        if let Some(sink) = ctx.event_log_sink.as_deref_mut() {
            sink.extend(b.seal_event_log_observer());
        }
        b.seal()
    }

    /// Union of carried alleles across both haplotypes for a segment —
    /// the support advertised to the feasibility/schedule analyzer.
    fn union_support(&self, seg: Segment) -> Vec<(AlleleId, f64)> {
        let mut ids: Vec<AlleleId> = Vec::new();
        for c in 0..2 {
            for id in self.genotype.haplotype(c).carried_alleles(seg) {
                if !ids.contains(&id) {
                    ids.push(id);
                }
            }
        }
        ids.into_iter().map(|id| (id, 1.0)).collect()
    }

    /// Candidate `(allele, mass)` copies in a gene slot on chromosome
    /// `c`. `mass = weight * copies` (copy-number dosage). When `filter`
    /// is set, copies are restricted to those admissible under the active
    /// contracts + feasibility.
    fn slot_candidates(
        &self,
        seg: Segment,
        gene: GeneId,
        c: usize,
        sim: &Simulation,
        ctx: &PassContext,
        filter: bool,
    ) -> Vec<(AlleleId, f64)> {
        self.genotype
            .haplotype(c)
            .slot(seg, gene)
            .iter()
            .filter(|cp| !filter || self.allele_feasible(c, seg, cp.allele, sim, ctx))
            .map(|cp| (cp.allele, cp.weight as f64 * cp.copies as f64))
            .collect()
    }

    /// Live (fresh-RNG) sampling of one segment within chromosome `c`.
    /// Records in canonical order: `sample_gene` → `sample_allele_in_slot`
    /// (only when the gene slot has >1 raw copy) → `sample_allele`.
    fn sample_segment_live(
        &self,
        seg: Segment,
        c: usize,
        sim: Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        let hap = self.genotype.haplotype(c);
        let vseg = Self::vseg(seg);
        // Filter-then-fallback, mirroring SampleAllelePass: prefer
        // feasibility-admissible candidates; only when NONE are admissible
        // do we fall back to the unfiltered carried set (permissive), or
        // raise (strict). `filter` is true whenever feasible candidates
        // exist on this chromosome+segment.
        let feasible_exists = hap
            .present_genes(seg)
            .any(|g| !self.slot_candidates(seg, g, c, &sim, ctx, true).is_empty());
        if strict && !feasible_exists {
            return Err(PassError::constraint_sampling(
                self.name(),
                address::sample_allele_vdj(seg),
                FilteredSampleError::EmptyAdmissibleSupport,
            ));
        }
        let filter = strict || feasible_exists;
        // Genes present on this chromosome that have >=1 candidate copy,
        // weighted by gene usage * total copy dosage (so a duplicated
        // gene recombines more often).
        let mut genes: Vec<(GeneId, f64)> = Vec::new();
        for g in hap.present_genes(seg) {
            let cands = self.slot_candidates(seg, g, c, &sim, ctx, filter);
            if cands.is_empty() {
                continue;
            }
            let dosage: f64 = cands.iter().map(|(_, m)| *m).sum();
            genes.push((g, self.usage_of(seg, g) * dosage));
        }
        if genes.is_empty() {
            return Err(PassError::constraint_sampling(
                self.name(),
                address::sample_allele_vdj(seg),
                FilteredSampleError::EmptyAdmissibleSupport,
            ));
        }
        let gene = Self::weighted_pick(&genes, ctx.rng);
        let cands = self.slot_candidates(seg, gene, c, &sim, ctx, filter);
        let id = Self::weighted_pick(&cands, ctx.rng);

        // "Multi-copy slot" is decided by the RAW slot length (genotype
        // structure), independent of feasibility filtering, so replay can
        // reconstruct whether a within-slot record exists from the
        // genotype alone.
        let multi = hap.slot(seg, gene).len() > 1;
        ctx.trace.record_choice(
            ChoiceAddress::SampleGene(vseg),
            ChoiceValue::GeneId(gene.index()),
        );
        if multi {
            ctx.trace.record_choice(
                ChoiceAddress::SampleAlleleInSlot(vseg),
                ChoiceValue::AlleleId(id.index()),
            );
        }
        ctx.trace.record_choice(
            ChoiceAddress::SampleAllele(vseg),
            ChoiceValue::AlleleId(id.index()),
        );
        Ok(self.commit(seg, sim, id, c, ctx))
    }

    /// Replay (trace-injected) sampling of one segment within `c`.
    /// Consumes records in the same order live emits them — `sample_gene`,
    /// then `sample_allele_in_slot` (iff the raw slot has >1 copy), then
    /// the canonical `sample_allele` (the assigned id, source of truth).
    fn sample_segment_replay(
        &self,
        seg: Segment,
        c: usize,
        sim: Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        let vseg = Self::vseg(seg);
        let gene_idx = ctx
            .replay_cursor
            .as_deref_mut()
            .expect("replay cursor present")
            .expect_gene_id(ChoiceAddress::SampleGene(vseg))
            .map_err(|r| PassError::replay(self.name(), r))?;
        let gene = GeneId::new(gene_idx);
        let slot = self.genotype.haplotype(c).slot(seg, gene);
        // Replay validation: a recorded trace must be admissible against
        // this genotype — the recorded gene must be carried on the drawn
        // chromosome (non-empty slot). Mirrors SampleAllelePass's
        // "trace proposes, engine validates" contract.
        if slot.is_empty() {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::sample_allele_vdj(seg),
                gene_idx as i64,
                "genotype_gene_not_carried_on_haplotype",
            ));
        }
        let multi = slot.len() > 1;
        let slot_recorded = if multi {
            Some(
                ctx.replay_cursor
                    .as_deref_mut()
                    .expect("replay cursor present")
                    .expect_allele_id(ChoiceAddress::SampleAlleleInSlot(vseg))
                    .map_err(|r| PassError::replay(self.name(), r))?,
            )
        } else {
            None
        };
        let allele = ctx
            .replay_cursor
            .as_deref_mut()
            .expect("replay cursor present")
            .expect_allele_id(ChoiceAddress::SampleAllele(vseg))
            .map_err(|r| PassError::replay(self.name(), r))?;
        // The canonical allele must be one the gene slot actually carries.
        if !slot.iter().any(|cp| cp.allele.index() == allele) {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::sample_allele_vdj(seg),
                allele as i64,
                "genotype_allele_not_in_gene_slot",
            ));
        }
        // When a within-slot record exists it must agree with the
        // canonical allele (live writes the same id to both).
        if let Some(slot_id) = slot_recorded {
            if slot_id != allele {
                return Err(PassError::invalid_distribution_output(
                    self.name(),
                    address::sample_allele_vdj(seg),
                    slot_id as i64,
                    "genotype_slot_record_disagrees_with_canonical_allele",
                ));
            }
        }
        let id = AlleleId::new(allele);
        if ctx.feasibility.is_some_and(|f| f.has_p_productive())
            && !self.allele_feasible(c, seg, id, &sim, ctx)
        {
            return Err(PassError::constraint_sampling(
                self.name(),
                address::sample_allele_vdj(seg),
                FilteredSampleError::EmptyAdmissibleSupport,
            ));
        }

        ctx.trace.record_choice(
            ChoiceAddress::SampleGene(vseg),
            ChoiceValue::GeneId(gene_idx),
        );
        if let Some(slot_id) = slot_recorded {
            ctx.trace.record_choice(
                ChoiceAddress::SampleAlleleInSlot(vseg),
                ChoiceValue::AlleleId(slot_id),
            );
        }
        ctx.trace.record_choice(
            ChoiceAddress::SampleAllele(vseg),
            ChoiceValue::AlleleId(id.index()),
        );
        Ok(self.commit(seg, sim, id, c, ctx))
    }

    /// Shared execute body. `strict` selects feasibility-filtered
    /// viability/sampling (and structured errors) vs presence-based
    /// permissive sampling (feasibility advisory, never errors given a
    /// complete haplotype).
    fn run(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        let replaying = ctx.replay_cursor.is_some();
        let c = if replaying {
            let recorded = ctx
                .replay_cursor
                .as_deref_mut()
                .expect("replay cursor present")
                .expect_haplotype(ChoiceAddress::SampleHaplotype)
                .map_err(|r| PassError::replay(self.name(), r))?;
            let viable = self.viable_set(sim, ctx, strict);
            if !viable.contains(&(recorded as usize)) {
                return Err(self.infeasible_error());
            }
            ctx.trace.record_choice(
                ChoiceAddress::SampleHaplotype,
                ChoiceValue::Haplotype(recorded),
            );
            recorded as usize
        } else {
            let viable = self.viable_set(sim, ctx, strict);
            if viable.is_empty() {
                return Err(self.infeasible_error());
            }
            let c = self.draw_haplotype(&viable, ctx.rng);
            ctx.trace.record_choice(
                ChoiceAddress::SampleHaplotype,
                ChoiceValue::Haplotype(c as u8),
            );
            c
        };

        let mut current = sim.with_live_call_support(self.live_call_support.clone());
        for seg in self.segments() {
            current = if replaying {
                self.sample_segment_replay(seg, c, current, ctx)?
            } else {
                self.sample_segment_live(seg, c, current, ctx, strict)?
            };
        }
        Ok(current)
    }
}

impl Pass for SampleGenotypePass {
    fn name(&self) -> &str {
        "sample_genotype"
    }

    /// Encode the attached genotype so two different genotypes produce
    /// distinct plan signatures (replay-cache correctness): source
    /// cartridge hash, subject, chromosome weights, and the full
    /// per-haplotype/per-segment carried slots (allele id, copies,
    /// weight bits).
    fn parameter_signature(&self) -> String {
        use std::fmt::Write;
        let g = &self.genotype;
        let mut s = String::new();
        let _ = write!(
            s,
            "geno|hash={}|subj={}|cw={:?}|d={}",
            g.source_refdata_hash(),
            g.subject_id().unwrap_or(""),
            g.chromosome_weights(),
            self.d_required,
        );
        for c in 0..2usize {
            let hap = g.haplotype(c);
            for seg in self.segments() {
                for gene in hap.present_genes(seg) {
                    let mut copies: Vec<(u32, u8, u32)> = hap
                        .slot(seg, gene)
                        .iter()
                        .map(|cp| (cp.allele.index(), cp.copies, cp.weight.to_bits()))
                        .collect();
                    copies.sort_unstable();
                    let _ = write!(s, "|h{}.{:?}.g{}={:?}", c, seg, gene.index(), copies);
                }
            }
        }
        s
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        // Permissive: viability is presence-based (feasibility advisory),
        // so a genotype with >=1 complete haplotype (guaranteed by the
        // compile-time presence check) never errors here.
        self.run(sim, ctx, false)
            .expect("SampleGenotypePass permissive execution must not error")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.run(sim, ctx, true)
    }

    fn declared_choice_patterns(&self) -> Vec<ChoiceAddressPattern> {
        let mut patterns = vec![ChoiceAddressPattern::SampleHaplotype];
        for seg in self.segments() {
            let vseg = Self::vseg(seg);
            patterns.push(ChoiceAddressPattern::SampleGene(vseg));
            patterns.push(ChoiceAddressPattern::SampleAlleleInSlot(vseg));
            patterns.push(ChoiceAddressPattern::SampleAllele(vseg));
        }
        patterns
    }

    fn effects(&self) -> Vec<PassEffect> {
        self.segments()
            .into_iter()
            .map(PassEffect::AssignAllele)
            .collect()
    }

    fn compile_facts(&self) -> Vec<PassCompileFact> {
        self.segments()
            .into_iter()
            .map(|seg| PassCompileFact::AlleleSampleSupport {
                segment: seg,
                support: AlleleIdSupport::from_weighted_pairs(Some(self.union_support(seg))),
            })
            .collect()
    }
}

#[cfg(test)]
pub(crate) mod test_support {
    use super::*;
    use crate::genotype::{GeneCopy, Haplotype};

    pub fn copy(id: u32) -> GeneCopy {
        GeneCopy {
            allele: AlleleId::new(id),
            copies: 1,
            weight: 1.0,
        }
    }

    /// hap0 carries V (gene0 -> allele 0) + J (gene0 -> allele 100);
    /// hap1 carries V but no J → only hap0 viable.
    pub fn geno_chrom1_deletes_j() -> Genotype {
        let mut h0 = Haplotype::new();
        h0.set(Segment::V, GeneId::new(0), vec![copy(0)]);
        h0.set(Segment::J, GeneId::new(0), vec![copy(100)]);
        let mut h1 = Haplotype::new();
        h1.set(Segment::V, GeneId::new(0), vec![copy(1)]);
        Genotype::new(
            [h0, h1],
            [0.5, 0.5],
            Some("S1".into()),
            "sha256:test".into(),
        )
    }

    /// Neither haplotype carries a J gene → no viable haplotype.
    pub fn geno_both_delete_j() -> Genotype {
        let mut h0 = Haplotype::new();
        h0.set(Segment::V, GeneId::new(0), vec![copy(0)]);
        let mut h1 = Haplotype::new();
        h1.set(Segment::V, GeneId::new(0), vec![copy(1)]);
        Genotype::new(
            [h0, h1],
            [0.5, 0.5],
            Some("S1".into()),
            "sha256:test".into(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pass::testing::PassRuntime;
    use crate::pass::PassPlan;

    fn productive_p_fixture() -> (
        SampleGenotypePass,
        crate::refdata::RefDataConfig,
        crate::feasibility::FeasibilityContext,
    ) {
        use crate::address::PEnd;
        use crate::feasibility::{
            ChoiceDomain, FeasibilityContext, PProductiveFeasibility, SegmentCompletionDomain,
        };
        use crate::genotype::{Genotype, Haplotype};
        use crate::refdata::{Allele, ChainType, RefDataConfig};

        // Fixed V3 P=1 and N=0: only V0+J2 and V2+J0 are in frame.
        // hap0 carries V0+J0 (impossible); hap1 carries V0,V2+J2,
        // where only V0 is feasible despite V2's cross-haplotype support.
        let mut refdata = RefDataConfig::empty(ChainType::Vj);
        for (segment, gene, seq, anchor) in [
            (Segment::V, "V0", b"TGT".as_slice(), 0),
            (Segment::V, "V2", b"TGTCC".as_slice(), 0),
            (Segment::J, "J0", b"TTTGGG".as_slice(), 0),
            (Segment::J, "J2", b"CCTTTGGG".as_slice(), 2),
        ] {
            let allele = Allele {
                name: format!("{gene}*01"),
                gene: gene.into(),
                seq: seq.to_vec(),
                segment,
                anchor: Some(anchor),
                functional_status: None,
                subregions: Vec::new(),
            };
            match segment {
                Segment::V => {
                    let _ = refdata.v_pool.push(allele);
                }
                Segment::J => {
                    let _ = refdata.j_pool.push(allele);
                }
                _ => unreachable!(),
            }
        }
        let mut h0 = Haplotype::new();
        h0.set(Segment::V, GeneId::new(0), vec![test_support::copy(0)]);
        h0.set(Segment::J, GeneId::new(0), vec![test_support::copy(0)]);
        let mut h1 = Haplotype::new();
        h1.set(Segment::V, GeneId::new(0), vec![test_support::copy(0)]);
        h1.set(Segment::V, GeneId::new(1), vec![test_support::copy(1)]);
        h1.set(Segment::J, GeneId::new(1), vec![test_support::copy(1)]);
        let genotype = Arc::new(Genotype::new([h0, h1], [0.5, 0.5], None, "fixture".into()));
        let feasibility = FeasibilityContext::empty().with_p_productive(PProductiveFeasibility {
            processing: vec![],
            segments: [Segment::V, Segment::J]
                .into_iter()
                .map(|segment| SegmentCompletionDomain {
                    segment,
                    alleles: ChoiceDomain::new(Some(0), vec![AlleleId::new(0), AlleleId::new(1)]),
                    trim_5: ChoiceDomain::new(None, vec![0]),
                    trim_3: ChoiceDomain::new(None, vec![0]),
                })
                .collect(),
            p_lengths: vec![(PEnd::V3, ChoiceDomain::new(Some(1), vec![1]))],
            np: vec![],
            d_orientation: ChoiceDomain::new(None, vec![false]),
            check_frame: true,
            check_stops: true,
            check_coding: true,
            check_v_anchor: true,
            check_j_anchor: true,
        });
        (
            SampleGenotypePass::new(genotype, false, vec![], vec![], vec![]),
            refdata,
            feasibility,
        )
    }

    #[test]
    fn productive_p_completion_stays_on_the_drawn_haplotype() {
        use crate::trace::Trace;

        let (pass, refdata, feasibility) = productive_p_fixture();
        for seed in 0..32 {
            let mut trace = Trace::new();
            let mut rng = Rng::new(seed);
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index: 0,
                refdata: Some(&refdata),
                contracts: None,
                feasibility: Some(&feasibility),
                reference_index: None,
                replay_cursor: None,
                event_log_sink: None,
            };
            let sim = Simulation::new();
            // Marginals alone admit both bad decisions; chromosome-aware
            // completion must reject them before either is committed.
            assert!(feasibility.admits(
                0,
                &sim,
                Some(&refdata),
                "sample_allele.v",
                &ChoiceValue::AlleleId(0),
            ));
            assert!(!pass.allele_feasible(0, Segment::V, AlleleId::new(0), &sim, &ctx));
            assert!(feasibility.admits(
                0,
                &sim,
                Some(&refdata),
                "sample_allele.v",
                &ChoiceValue::AlleleId(1),
            ));
            assert!(!pass.allele_feasible(1, Segment::V, AlleleId::new(1), &sim, &ctx));
            assert_eq!(pass.viable_set(&sim, &ctx, true), vec![1]);
            let result = pass
                .run(&sim, &mut ctx, true)
                .expect("hap1 has a productive completion");
            let v = result.assignments.get(Segment::V).unwrap();
            let j = result.assignments.get(Segment::J).unwrap();
            assert_eq!((v.haplotype, v.allele_id), (Some(1), AlleleId::new(0)));
            assert_eq!((j.haplotype, j.allele_id), (Some(1), AlleleId::new(1)));
        }
    }

    #[test]
    fn productive_p_completion_does_not_use_zero_usage_future_alleles() {
        use crate::trace::Trace;

        let (mut pass, refdata, feasibility) = productive_p_fixture();
        pass.usage_v = vec![(GeneId::new(0), 0.0), (GeneId::new(1), 1.0)];
        let mut trace = Trace::new();
        let mut rng = Rng::new(0);
        let ctx = PassContext {
            trace: &mut trace,
            rng: &mut rng,
            pass_index: 0,
            refdata: Some(&refdata),
            contracts: None,
            feasibility: Some(&feasibility),
            reference_index: None,
            replay_cursor: None,
            event_log_sink: None,
        };
        assert!(pass.viable_set(&Simulation::new(), &ctx, true).is_empty());
    }

    #[test]
    fn productive_p_replay_rejects_an_allele_only_completable_on_another_haplotype() {
        use crate::replay::TraceCursor;
        use crate::trace::Trace;

        let (pass, refdata, feasibility) = productive_p_fixture();
        let mut trace = Trace::new();
        let mut rng = Rng::new(0);
        let sim = Simulation::new();
        {
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index: 0,
                refdata: Some(&refdata),
                contracts: None,
                feasibility: Some(&feasibility),
                reference_index: None,
                replay_cursor: None,
                event_log_sink: None,
            };
            pass.run(&sim, &mut ctx, true).expect("live trace");
        }
        let mut records = trace.choices().to_vec();
        for record in &mut records {
            match record.address.as_str() {
                "sample_gene.v" => record.value = ChoiceValue::GeneId(1),
                "sample_allele.v" => record.value = ChoiceValue::AlleleId(1),
                _ => {}
            }
        }
        let mut cursor = TraceCursor::from_owned(records);
        let mut replay_trace = Trace::new();
        let mut ctx = PassContext {
            trace: &mut replay_trace,
            rng: &mut rng,
            pass_index: 0,
            refdata: Some(&refdata),
            contracts: None,
            feasibility: Some(&feasibility),
            reference_index: None,
            replay_cursor: Some(&mut cursor),
            event_log_sink: None,
        };
        assert!(pass.run(&sim, &mut ctx, true).is_err());
    }

    #[test]
    fn draws_only_the_viable_haplotype_and_assigns_carried_alleles() {
        let g = Arc::new(test_support::geno_chrom1_deletes_j());
        let pass = SampleGenotypePass::new(g, false, vec![], vec![], vec![]);
        let mut plan = PassPlan::new();
        plan.push(Box::new(pass));
        for seed in 0..30u64 {
            let outcome = PassRuntime::execute(&plan, Simulation::new(), seed);
            match outcome.trace.find("sample_haplotype").unwrap().value {
                ChoiceValue::Haplotype(c) => assert_eq!(c, 0, "seed {seed}"),
                _ => panic!("wrong variant"),
            }
            let sim = outcome.final_simulation();
            assert_eq!(
                sim.assignments.get(Segment::V).unwrap().allele_id,
                AlleleId::new(0)
            );
            assert_eq!(
                sim.assignments.get(Segment::J).unwrap().allele_id,
                AlleleId::new(100)
            );
            // single-copy slots => no within-slot record
            assert!(outcome.trace.find("sample_allele_in_slot.v").is_none());
        }
    }

    #[test]
    fn commit_stamps_drawn_haplotype_on_assignments() {
        // Only hap0 is viable (hap1 lacks J), so the drawn chromosome is
        // always 0; the V and J assignments must carry haplotype Some(0).
        let g = Arc::new(test_support::geno_chrom1_deletes_j());
        let pass = SampleGenotypePass::new(g, false, vec![], vec![], vec![]);
        let mut plan = PassPlan::new();
        plan.push(Box::new(pass));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 7);
        let drawn = match outcome.trace.find("sample_haplotype").unwrap().value {
            ChoiceValue::Haplotype(h) => h,
            _ => panic!("expected Haplotype record"),
        };
        assert_eq!(drawn, 0);
        let sim = outcome.final_simulation();
        for seg in [Segment::V, Segment::J] {
            assert_eq!(
                sim.assignments.get(seg).unwrap().haplotype,
                Some(drawn),
                "{seg:?} assignment must carry the drawn haplotype",
            );
        }
    }

    #[test]
    fn records_canonical_sample_allele_and_gene_addresses() {
        let g = Arc::new(test_support::geno_chrom1_deletes_j());
        let pass = SampleGenotypePass::new(g, false, vec![], vec![], vec![]);
        let mut plan = PassPlan::new();
        plan.push(Box::new(pass));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 7);
        match outcome.trace.find("sample_allele.v").unwrap().value {
            ChoiceValue::AlleleId(id) => assert_eq!(id, 0),
            _ => panic!("expected AlleleId at sample_allele.v"),
        }
        match outcome.trace.find("sample_gene.v").unwrap().value {
            ChoiceValue::GeneId(g) => assert_eq!(g, 0),
            _ => panic!("expected GeneId at sample_gene.v"),
        }
    }

    #[test]
    fn errors_when_no_haplotype_viable() {
        let g = Arc::new(test_support::geno_both_delete_j());
        let pass = SampleGenotypePass::new(g, false, vec![], vec![], vec![]);
        let mut plan = PassPlan::new();
        plan.push(Box::new(pass));
        let result =
            PassRuntime::execute_strict_with_context(&plan, Simulation::new(), 0, None, None);
        assert!(result.is_err(), "expected genotype-infeasibility error");
    }

    #[test]
    fn replay_reproduces_live_phased_choices() {
        use crate::pass::PassContext;
        use crate::replay::TraceCursor;
        use crate::rng::Rng;
        use crate::trace::Trace;

        let g = Arc::new(test_support::geno_chrom1_deletes_j());
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleGenotypePass::new(
            g.clone(),
            false,
            vec![],
            vec![],
            vec![],
        )));
        // Live run — capture the full trace (haplotype + gene + allele).
        let live = PassRuntime::execute(&plan, Simulation::new(), 13);
        let records: Vec<_> = live.trace.choices().to_vec();
        let live_sim = live.final_simulation();
        let live_v = live_sim.assignments.get(Segment::V).unwrap().allele_id;
        let live_j = live_sim.assignments.get(Segment::J).unwrap().allele_id;

        // Replay the captured trace through a fresh pass instance.
        let replay_pass = SampleGenotypePass::new(g, false, vec![], vec![], vec![]);
        let mut cursor = TraceCursor::from_owned(records);
        let mut trace = Trace::new();
        let mut rng = Rng::new(999);
        let sim = Simulation::new();
        let result;
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
            result = replay_pass.run(&sim, &mut ctx, true);
        }
        let replayed = result.expect("genotype replay must succeed");
        assert_eq!(
            replayed.assignments.get(Segment::V).unwrap().allele_id,
            live_v
        );
        assert_eq!(
            replayed.assignments.get(Segment::J).unwrap().allele_id,
            live_j
        );
    }

    #[test]
    fn replay_rejects_allele_not_carried_in_genotype_slot() {
        use crate::pass::PassContext;
        use crate::replay::TraceCursor;
        use crate::rng::Rng;
        use crate::trace::Trace;

        let g = Arc::new(test_support::geno_chrom1_deletes_j());
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleGenotypePass::new(
            g.clone(),
            false,
            vec![],
            vec![],
            vec![],
        )));
        let live = PassRuntime::execute(&plan, Simulation::new(), 13);
        // Tamper: rewrite the canonical V allele to one not in the slot.
        let mut records: Vec<_> = live.trace.choices().to_vec();
        for r in records.iter_mut() {
            if r.address == "sample_allele.v" {
                r.value = ChoiceValue::AlleleId(999);
            }
        }
        let replay_pass = SampleGenotypePass::new(g, false, vec![], vec![], vec![]);
        let mut cursor = TraceCursor::from_owned(records);
        let mut trace = Trace::new();
        let mut rng = Rng::new(1);
        let sim = Simulation::new();
        let result;
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
            result = replay_pass.run(&sim, &mut ctx, true);
        }
        assert!(
            result.is_err(),
            "replay must reject an allele not carried in the genotype slot"
        );
    }
}
