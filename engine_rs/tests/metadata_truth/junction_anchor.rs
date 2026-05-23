//! Junction anchor + frame + productive-flag invariants.
//!
//! The junction is the V→...→J region centered on the V's Cys and
//! J's Phe/Trp anchors. Its boundaries (`junction_start`,
//! `junction_end`), translated amino acids (`junction_aa`), in-frame
//! flag (`vj_in_frame`), stop codon flag (`stop_codon`), and
//! productive flag must remain self-consistent under every kind of
//! corruption the engine supports.
//!
//! ## Invariants this module should test
//!
//! - **Anchor preserved under synonymous mutation**: SHM that flips
//!   the third codon position to another base producing the same
//!   amino acid keeps the productive flag set.
//! - **Anchor broken under nonsynonymous mutation at the anchor
//!   codon**: SHM that changes the C/W amino acid resets productive
//!   to false. Live call set may or may not narrow (a different
//!   allele might naturally code that mutation).
//! - **Junction frame shifts under structural indel**: a 1- or 2-base
//!   indel before the J anchor shifts `vj_in_frame` to false.
//!   A 3-base indel preserves frame.
//! - **Stop codon in junction**: SHM that introduces a stop codon
//!   inside the junction sets `stop_codon = true` and `productive
//!   = false` regardless of frame.
//! - **End-loss removing the V anchor**: the pass strips bases such
//!   that the V Cys is gone. Productive must be false. v_call may
//!   widen (less informative bases left) or fall back to truth.
//! - **End-loss removing the J anchor**: similar.
//! - **Anchor lookup is germline_pos-based, not pool-offset based**:
//!   verify by inserting bases before the anchor and confirming
//!   `junction_start` shifts in pool coords but `junction_aa`
//!   recovers the same amino acids.
//!
//! ## Adjacent in-tree coverage
//!
//! - `engine_rs/src/airr_record/tests/anchors.rs`
//! - `engine_rs/src/contract/anchor_preserved/tests/`
//! - `engine_rs/src/compiled/tests/curated.rs::junction_locates_anchors_via_germline_pos`
//! - `engine_rs/src/compiled/tests/curated.rs::productive_uses_germline_pos_anchor_after_indents`

use super::common;
use genairr_engine::airr_record::{build_airr_record, AirrRecord};
use genairr_engine::compiled::{CompiledSimulator, ExecutionPolicy};
use genairr_engine::dist::{AllelePoolDist, EmpiricalLengthDist, UniformBase};
use genairr_engine::ir::{flag, NucHandle, Nucleotide, Segment, Simulation};
use genairr_engine::pass::{Pass, PassContext, PassEffect, PassPlan};
use genairr_engine::passes::{
    AssembleSegmentPass, EndLossPass, GenerateNPPass, LossEnd, SampleAllelePass, TrimPass,
};
use genairr_engine::refdata::{Allele, AlleleId, ChainType, RefDataConfig};

// ──────────────────────────────────────────────────────────────
// Test-local helpers.
// ──────────────────────────────────────────────────────────────

#[derive(Clone, Debug)]
struct EditBaseAtPass {
    handle: NucHandle,
    new_base: u8,
}

impl EditBaseAtPass {
    fn new(pool_index: u32, new_base: u8) -> Self {
        Self {
            handle: NucHandle::new(pool_index),
            new_base,
        }
    }
}

impl Pass for EditBaseAtPass {
    fn name(&self) -> &str {
        "test.edit_base_at"
    }
    fn execute(&self, sim: &Simulation, _ctx: &mut PassContext) -> Simulation {
        sim.with_base_changed(self.handle, self.new_base)
    }
    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::EditBases]
    }
}

#[derive(Clone, Debug)]
struct InsertAtPass {
    at: u32,
    base: u8,
    segment: Segment,
}

impl InsertAtPass {
    fn new(at: u32, base: u8, segment: Segment) -> Self {
        Self { at, base, segment }
    }
}

impl Pass for InsertAtPass {
    fn name(&self) -> &str {
        "test.insert_at"
    }
    fn execute(&self, sim: &Simulation, _ctx: &mut PassContext) -> Simulation {
        let nuc = Nucleotide::synthetic(self.base, self.segment, flag::INDEL_INSERTED);
        sim.with_indel_inserted(self.at, nuc)
    }
    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::StructuralIndel]
    }
}

/// VJ refdata with biologically-real anchors:
/// - V (12 bp, anchor=6): `AAACCCTGTAAA` — Cys (TGT) codon at pos
///   6-8, distinguishing tail `AAA` at pos 9-11.
/// - V2 (12 bp, anchor=6): `AAACCCTGTGGG` — same Cys codon,
///   distinguishing tail `GGG`.
/// - J (9 bp, anchor=3): `AAATGGACG` — Trp (TGG) codon at pos 3-5.
/// Single J entry — keeps the J call deterministic.
///
/// Baseline assembly (no trim, no NP): V + J = 21 bases.
///   pool[0..12] = V; pool[12..21] = J.
///   V anchor pool=6 (TGT). J anchor pool=15 (TGG, since J starts
///   at pool 12 + anchor=3).
///   Junction = pool[6..18] = TGT AAA AAA TGG  (3+3+3+3 = 12 bytes).
///   Codons: TGT AAA AAA TGG = C K K W → productive.
fn productive_vj_refdata() -> RefDataConfig {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let _ = cfg.v_pool.push(Allele {
        name: "Vp*01".into(),
        gene: "Vp".into(),
        seq: b"AAACCCTGTAAA".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
    });
    let _ = cfg.v_pool.push(Allele {
        name: "Vp*02".into(),
        gene: "Vp".into(),
        seq: b"AAACCCTGTGGG".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "Jp*01".into(),
        gene: "Jp".into(),
        seq: b"AAATGGACG".to_vec(),
        segment: Segment::J,
        anchor: Some(3),
    });
    cfg
}

fn vj_baseline_plan(cfg: &RefDataConfig, v_id: AlleleId, j_id: AlleleId) -> PassPlan {
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.v_pool, vec![v_id])),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::restricted_uniform(&cfg.j_pool, vec![j_id])),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
    plan
}

fn run_plan_with_record(
    cfg: &RefDataConfig,
    plan: PassPlan,
    seed: u64,
) -> (Simulation, AirrRecord) {
    let compiled = CompiledSimulator::compile(&plan, Some(cfg), None, ExecutionPolicy::Permissive)
        .expect("plan should compile");
    let outcome = compiled.run_one(seed).expect("plan should run");
    let rec = build_airr_record(&outcome, cfg, "test");
    (outcome.final_simulation().clone(), rec)
}

// ──────────────────────────────────────────────────────────────
// Tests.
// ──────────────────────────────────────────────────────────────

#[test]
fn baseline_no_corruption_is_productive() {
    // Sanity: with no corruption, V + J assemble produces a 12-byte
    // junction (CKKW = C K K W) — productive, in-frame, no stop.
    let cfg = productive_vj_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "Vp*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "Jp*01");
    let plan = vj_baseline_plan(&cfg, v_id, j_id);
    let (_sim, rec) = run_plan_with_record(&cfg, plan, 0);

    assert_eq!(rec.junction_start, Some(6));
    assert_eq!(rec.junction_end, Some(18));
    assert_eq!(rec.junction_length, Some(12));
    assert_eq!(rec.junction, "TGTAAAAAATGG");
    assert_eq!(rec.junction_aa, "CKKW");
    assert_eq!(rec.productive, Some(true));
    assert_eq!(rec.vj_in_frame, Some(true));
    assert_eq!(rec.stop_codon, Some(false));
}

#[test]
fn synonymous_mutation_at_v_anchor_third_position_preserves_productive() {
    // Vp's anchor codon is TGT (Cys) at pool 6-8. The third position
    // (pool 8) can flip T → C without changing the amino acid: TGC
    // also encodes Cys. The productive flag must remain true.
    let cfg = productive_vj_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "Vp*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "Jp*01");
    let mut plan = vj_baseline_plan(&cfg, v_id, j_id);
    plan.push(Box::new(EditBaseAtPass::new(8, b'C')));
    let (_sim, rec) = run_plan_with_record(&cfg, plan, 0);

    assert!(
        rec.junction_aa.starts_with('C'),
        "synonymous edit must preserve V Cys anchor; junction_aa = {:?}",
        rec.junction_aa,
    );
    assert_eq!(
        rec.productive,
        Some(true),
        "synonymous V-anchor edit must preserve productive flag (junction_aa={:?})",
        rec.junction_aa,
    );
    assert_eq!(rec.vj_in_frame, Some(true));
    assert_eq!(rec.stop_codon, Some(false));
}

#[test]
fn nonsynonymous_mutation_at_v_anchor_flips_productive_to_false() {
    // Vp anchor TGT (Cys) at pool 6-8. Flip pool 6 T → A: anchor
    // codon becomes AGT (Ser) → no longer Cys → productive=false.
    let cfg = productive_vj_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "Vp*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "Jp*01");
    let mut plan = vj_baseline_plan(&cfg, v_id, j_id);
    plan.push(Box::new(EditBaseAtPass::new(6, b'A')));
    let (_sim, rec) = run_plan_with_record(&cfg, plan, 0);

    // Anchor codon is now AGT (Ser).
    assert!(
        !rec.junction_aa.starts_with('C'),
        "nonsynonymous edit must destroy Cys; junction_aa = {:?}",
        rec.junction_aa,
    );
    assert_eq!(
        rec.productive,
        Some(false),
        "nonsynonymous V-anchor edit must flip productive to false (junction_aa={:?})",
        rec.junction_aa,
    );
}

#[test]
fn stop_codon_introduced_in_junction_flips_productive_and_stop_codon() {
    // Insert a stop codon (TAA) inside the junction by editing
    // pool 9, 10, 11 (the AAA triplet after the V anchor).
    // AAA → TAA: that's a stop. Productive must flip to false and
    // stop_codon to true regardless of frame.
    let cfg = productive_vj_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "Vp*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "Jp*01");
    let mut plan = vj_baseline_plan(&cfg, v_id, j_id);
    plan.push(Box::new(EditBaseAtPass::new(9, b'T')));
    let (_sim, rec) = run_plan_with_record(&cfg, plan, 0);

    assert!(
        rec.junction_aa.contains('*') || rec.stop_codon == Some(true),
        "stop codon must be recorded in junction (junction_aa={:?}, stop_codon={:?})",
        rec.junction_aa,
        rec.stop_codon,
    );
    assert_eq!(
        rec.stop_codon,
        Some(true),
        "stop_codon flag must be true; junction_aa={:?}",
        rec.junction_aa,
    );
    assert_eq!(
        rec.productive,
        Some(false),
        "productive must flip to false when stop codon is in junction",
    );
}

#[test]
fn one_base_indel_before_j_anchor_flips_vj_in_frame_false() {
    // Insert a 1-base 'C' inside V (pool 10, inside the Vp distinguishing
    // tail). The junction length grows by 1 (was 12, becomes 13) →
    // 13 % 3 != 0 → vj_in_frame = false.
    let cfg = productive_vj_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "Vp*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "Jp*01");
    let mut plan = vj_baseline_plan(&cfg, v_id, j_id);
    plan.push(Box::new(InsertAtPass::new(10, b'C', Segment::V)));
    let (_sim, rec) = run_plan_with_record(&cfg, plan, 0);

    assert_eq!(rec.junction_length, Some(13));
    assert_eq!(
        rec.vj_in_frame,
        Some(false),
        "single-base insertion before J anchor must flip vj_in_frame to false",
    );
    assert_eq!(
        rec.productive,
        Some(false),
        "out-of-frame junction must flip productive to false",
    );
}

#[test]
fn two_base_indel_before_j_anchor_also_flips_frame() {
    // Same fixture, but two synthetic insertions inside V → junction
    // length = 14 (still not divisible by 3) → vj_in_frame=false.
    let cfg = productive_vj_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "Vp*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "Jp*01");
    let mut plan = vj_baseline_plan(&cfg, v_id, j_id);
    plan.push(Box::new(InsertAtPass::new(10, b'C', Segment::V)));
    plan.push(Box::new(InsertAtPass::new(10, b'C', Segment::V)));
    let (_sim, rec) = run_plan_with_record(&cfg, plan, 0);

    assert_eq!(rec.junction_length, Some(14));
    assert_eq!(rec.vj_in_frame, Some(false));
}

#[test]
fn end_loss_removing_v_anchor_makes_productive_false() {
    // The V anchor (Cys, pool 6-8 baseline) is removed by a 5'
    // end-loss of 7 bases: pool[0..7] gets stripped, including
    // anchor pool 6. After the loss, the junction can't form
    // around a Cys anchor → productive=false (or junction
    // coordinates absent).
    let cfg = productive_vj_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "Vp*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "Jp*01");
    let mut plan = vj_baseline_plan(&cfg, v_id, j_id);
    plan.push(Box::new(EndLossPass::new(
        LossEnd::Five,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(7, 1.0)])),
    )));
    let (_sim, rec) = run_plan_with_record(&cfg, plan, 0);

    // Either productive is explicitly false, or junction coords
    // are absent — both are valid expressions of "V anchor is
    // gone".
    let productive_lost = rec.productive == Some(false) || rec.productive.is_none();
    assert!(
        productive_lost,
        "removing V anchor must flip productive away from true; got {:?}",
        rec.productive,
    );
}

#[test]
fn junction_frame_correctness_under_v_trim() {
    // V_3 trim shrinks the V region but does NOT shift the V
    // anchor (anchor germline_pos=6 sits within the surviving
    // 5' bases). Junction span shrinks accordingly — still
    // in-frame as long as the trim is a multiple of 3.
    let cfg = productive_vj_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "Vp*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "Jp*01");

    let make_plan = |v_trim_3: i64| -> PassPlan {
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::restricted_uniform(&cfg.v_pool, vec![v_id])),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::restricted_uniform(&cfg.j_pool, vec![j_id])),
        )));
        if v_trim_3 > 0 {
            plan.push(Box::new(TrimPass::new(
                Segment::V,
                genairr_engine::assignment::TrimEnd::Three,
                Box::new(EmpiricalLengthDist::from_pairs(vec![(v_trim_3, 1.0)])),
            )));
        }
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
            Box::new(UniformBase),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));
        plan
    };

    // 3-byte trim — keeps frame (3 % 3 = 0).
    let (_s3, rec3) = run_plan_with_record(&cfg, make_plan(3), 0);
    assert_eq!(
        rec3.vj_in_frame,
        Some(true),
        "3-byte V_3 trim preserves frame (got {:?})",
        rec3.junction_length,
    );
    // 1-byte trim — out of frame (1 % 3 != 0).
    let (_s1, rec1) = run_plan_with_record(&cfg, make_plan(1), 0);
    // Verify the junction-length parity reflects the trim.
    if let (Some(len3), Some(len1)) = (rec3.junction_length, rec1.junction_length) {
        assert!(
            len3 % 3 == 0,
            "3-byte trim must produce in-frame junction len {} (3 | len)",
            len3,
        );
        assert!(
            len1 % 3 != 0,
            "1-byte trim must produce out-of-frame junction len {} (3 ∤ len)",
            len1,
        );
        assert_eq!(rec1.vj_in_frame, Some(false));
    }
}

#[test]
fn synthetic_inserted_byte_inside_junction_appears_in_junction_aa() {
    // The productive VJ fixture has anchors:
    //   V Vp*01 = AAACCCTGTAAA  (anchor=6, Cys codon TGT at pos 6-8)
    //   J Jp*01 = AAATGGACG     (anchor=3, Trp codon TGG at pos 3-5)
    // Baseline junction (no NP, no trim): pool[6..18] = "TGTAAAAAATGG"
    // → junction_aa = "CKKW".
    //
    // Insert a synthetic 3-base run inside V at pool 9 (right after
    // the Cys anchor codon, inside the junction window). The 3 bases
    // preserve frame (junction grows from 12 to 15 — still divisible
    // by 3). With three inserts of 'G' at the same position, the
    // junction nt becomes "TGT GGG AAAAAA TGG" — junction_aa
    // includes a 'G' codon (GGG = Glycine).
    //
    // Synthetic bytes have `germline_pos == NONE` (no allele evidence
    // contribution) but DO appear in the junction string and in
    // junction_aa.
    let cfg = productive_vj_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "Vp*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "Jp*01");
    let mut plan = vj_baseline_plan(&cfg, v_id, j_id);
    // Insert 3 synthetic 'G' bytes at pool 9 (inside the junction,
    // after the Cys anchor). Total junction grows by 3 — stays
    // in-frame.
    plan.push(Box::new(InsertAtPass::new(9, b'G', Segment::V)));
    plan.push(Box::new(InsertAtPass::new(9, b'G', Segment::V)));
    plan.push(Box::new(InsertAtPass::new(9, b'G', Segment::V)));
    let (_sim, rec) = run_plan_with_record(&cfg, plan, 0);

    // Junction length grew from 12 to 15 (still divisible by 3).
    assert_eq!(rec.junction_length, Some(15));
    assert_eq!(rec.vj_in_frame, Some(true));
    // The synthetic 'G' bytes appear in the junction nt string.
    assert!(
        rec.junction.contains("GGG"),
        "junction nt must contain the synthetic GGG triplet; got {:?}",
        rec.junction,
    );
    // The junction_aa must include the synthetic G codon (Glycine).
    assert!(
        rec.junction_aa.contains('G'),
        "junction_aa must include Glycine from synthetic GGG codon; got {:?}",
        rec.junction_aa,
    );
    // The Cys anchor at the start is still preserved (synthetic
    // bytes were inserted AFTER the anchor codon).
    assert!(
        rec.junction_aa.starts_with('C'),
        "junction_aa must still start with Cys; got {:?}",
        rec.junction_aa,
    );
}

#[test]
fn anchor_lookup_via_germline_pos_survives_insertion_before_anchor() {
    // Insert a synthetic 'C' at pool position 3 (inside V, BEFORE
    // the Cys anchor at pool 6). The anchor's germline_pos=6 now
    // sits at pool index 7. The anchor lookup must follow:
    // `junction_start` shifts from 6 → 7, but `junction_aa` still
    // starts with 'C' (anchor codon recovered via germline_pos).
    let cfg = productive_vj_refdata();
    let v_id = common::allele_id_by_name(&cfg, Segment::V, "Vp*01");
    let j_id = common::allele_id_by_name(&cfg, Segment::J, "Jp*01");
    let mut plan = vj_baseline_plan(&cfg, v_id, j_id);
    plan.push(Box::new(InsertAtPass::new(3, b'C', Segment::V)));
    let (_sim, rec) = run_plan_with_record(&cfg, plan, 0);

    assert_eq!(
        rec.junction_start,
        Some(7),
        "anchor germline_pos lookup must shift junction_start by the insertion offset",
    );
    assert!(
        rec.junction.starts_with("TGT"),
        "junction must still start with the V Cys codon, got {:?}",
        rec.junction,
    );
    assert!(
        rec.junction_aa.starts_with('C'),
        "junction_aa must still begin with Cys after pre-anchor insertion; junction_aa={:?}",
        rec.junction_aa,
    );
}
