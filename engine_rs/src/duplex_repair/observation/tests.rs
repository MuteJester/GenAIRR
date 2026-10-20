use super::*;
use crate::duplex_repair::{
    bindings::DuplexReferenceBindings, repair_coupled_annotated, seed::seed_from_coupled, JoinName,
    JoinScript, RepairOperation, Side,
};
use crate::ir::NucHandle;
use crate::refdata::{Allele, ChainType};

fn seed() -> (Simulation, RefDataConfig) {
    let mut r = RefDataConfig::empty(ChainType::Vdj);
    for (s, name, seq, a) in [
        (Segment::V, "V*01", "AAATGT", Some(3)),
        (Segment::D, "D*01", "CCC", None),
        (Segment::J, "J*01", "TTTAAA", Some(0)),
    ] {
        let allele = Allele {
            name: name.into(),
            gene: name.split('*').next().unwrap().into(),
            seq: seq.as_bytes().to_vec(),
            segment: s,
            anchor: a,
            functional_status: None,
            subregions: vec![],
        };
        match s {
            Segment::V => r.v_pool.push(allele),
            Segment::D => r.d_pool.push(allele),
            Segment::J => r.j_pool.push(allele),
            _ => unreachable!(),
        };
    }
    let scripts: Vec<_> = [(JoinName::VD, 6), (JoinName::DJ, 3)]
        .into_iter()
        .map(|(join, right_start)| JoinScript {
            join,
            operations: vec![
                RepairOperation::Open {
                    side: Side::Left,
                    offset: 0,
                },
                RepairOperation::Open {
                    side: Side::Right,
                    offset: 0,
                },
                RepairOperation::Synapse { right_start },
                RepairOperation::Ligate {},
            ],
        })
        .collect();
    let b = DuplexReferenceBindings::new(
        &r,
        ["V*01", "D*01", "J*01"],
        ["AAATGT", "CCC", "TTTAAA"],
        false,
        true,
    )
    .unwrap();
    let a =
        repair_coupled_annotated("AAATGT", "CCC", "TTTAAA", false, &scripts, 256, 1000).unwrap();
    (seed_from_coupled(a, b, &r).unwrap(), r)
}

#[test]
fn founder_is_immutable_and_current_substitution_is_not_repair_synthesis() {
    let (parent, r) = seed();
    let child = parent.with_base_changed(NucHandle::new(6), b'a');
    let projected = project_duplex_simulation(&child, &r).unwrap();
    assert_eq!(projected["record"]["sequence"], "AAATGTaCCTTTAAA");
    assert_eq!(
        projected["record"]["duplex_repair"]["product"]["sequence"],
        "AAATGTCCCTTTAAA"
    );
    assert_eq!(parent.pool.as_slice()[6].base, b'C');
    assert!(projected["record"]["n_mutations"].is_null());
    assert_eq!(projected["record"]["d_identity"], 2.0 / 3.0);
}

#[test]
fn complete_structural_d_loss_keeps_deletion_path_and_np_boundaries() {
    let (mut sim, r) = seed();
    for _ in 0..3 {
        sim = sim.with_indel_deleted(6);
    }
    let x = project_duplex_simulation(&sim, &r).unwrap();
    assert_eq!(x["record"]["sequence"], "AAATGTTTTAAA");
    assert_eq!(x["record"]["d_cigar"], "3D");
    assert_eq!(x["record"]["d_germline_start"], 0);
    assert_eq!(x["record"]["d_germline_end"], 3);
    assert!(x["record"]["d_sequence_start"].is_null());
    assert_eq!(x["paths"]["d"]["query_boundary"], 6);
    assert_eq!(x["record"]["np1_length"], 0);
    assert_eq!(x["record"]["np2_length"], 0);
    assert_eq!(
        x["record"]["duplex_repair"]["counts"]["founder_query"]["columns"],
        15
    );
}

#[test]
fn low_level_cache_tampering_cannot_be_exported_as_valid_duplex_metadata() {
    let (sim, r) = seed();
    let mut n = sim.pool.as_slice()[3];
    n.germline_pos = crate::ir::GermlinePos::NONE;
    let forged = sim.with_nucleotide_changed(NucHandle::new(3), n);
    assert!(project_duplex_simulation(&forged, &r).is_err());
}
