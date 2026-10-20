use super::*;
use crate::refdata::{Allele, ChainType};

fn data(v: &str, d: &str, j: &str, anchor: Option<u16>) -> RefDataConfig {
    let mut r = RefDataConfig::empty(ChainType::Vdj);
    for (s, name, seq, a) in [
        (Segment::V, "V*01", v, anchor),
        (Segment::D, "D*01", d, None),
        (Segment::J, "J*01", j, Some(0)),
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
    r
}
fn open(join: JoinName) -> JoinScript {
    JoinScript {
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
        ],
    }
}
fn blunt(join: JoinName, at: i64) -> JoinScript {
    let mut j = open(join);
    j.operations.extend([
        RepairOperation::Synapse { right_start: at },
        RepairOperation::Ligate {},
    ]);
    j
}
fn project(r: &RefDataConfig, joins: &[JoinScript], inverse: bool) -> Value {
    let seqs: [&str; 3] = std::array::from_fn(|i| {
        std::str::from_utf8(
            &r.get(
                [Segment::V, Segment::D, Segment::J][i],
                crate::refdata::AlleleId::new(0),
            )
            .unwrap()
            .seq,
        )
        .unwrap()
    });
    let b = DuplexReferenceBindings::new(r, ["V*01", "D*01", "J*01"], seqs, inverse, true).unwrap();
    let p = repair_coupled_annotated(seqs[0], seqs[1], seqs[2], inverse, joins, 256, 1000).unwrap();
    project_coupled(&p, &b, r).unwrap()
}

#[test]
fn restored_anchor_is_placed_without_becoming_an_original_nucleotide() {
    let r = data("AAATGT", "CCC", "TTTAAA", Some(3));
    let mut vd = open(JoinName::VD);
    vd.operations.extend([
        RepairOperation::Resect {
            side: Side::Left,
            strand: Strand::Top,
            count: 1,
        },
        RepairOperation::Synapse { right_start: 6 },
        RepairOperation::Fill {
            strand: Strand::Top,
        },
        RepairOperation::Ligate {},
    ]);
    let x = project(&r, &[vd, blunt(JoinName::DJ, 3)], false);
    let rec = &x["record"];
    assert_eq!(rec["sequence"], "AAATGTCCCTTTAAA");
    assert_eq!(rec["v_cigar"], "6M");
    assert_eq!(rec["productive"], true);
    assert_eq!(rec["junction"], "TGTCCCTTT");
    assert_eq!(rec["junction_aa"], "CPF");
    assert_eq!(rec["v_trim_3"], 0);
    assert_eq!(
        rec["duplex_repair"]["product"]["columns"][5]["top"]["kind"],
        "templated_synthesis"
    );
    assert!(rec["p_v_3_length"].is_null());
}

#[test]
fn true_body_overlap_preserves_both_alignments_and_omits_combined_alignment() {
    let r = data("ACG", "GTA", "CCC", None);
    let mut vd = open(JoinName::VD);
    vd.operations.extend([
        RepairOperation::Resect {
            side: Side::Left,
            strand: Strand::Bottom,
            count: 1,
        },
        RepairOperation::Resect {
            side: Side::Right,
            strand: Strand::Top,
            count: 1,
        },
        RepairOperation::Synapse { right_start: 2 },
        RepairOperation::Ligate {},
    ]);
    let x = project(&r, &[vd, blunt(JoinName::DJ, 3)], false);
    assert_eq!(x["record"]["sequence"], "ACGTACCC");
    assert_eq!(x["record"]["v_sequence_end"], 3);
    assert_eq!(x["record"]["d_sequence_start"], 2);
    assert_eq!(x["record"]["v_cigar"], "3M");
    assert_eq!(x["record"]["d_cigar"], "3M");
    assert!(x["record"]["sequence_alignment"].is_null());
    assert_eq!(x["record"]["np1_length"], 0);
    assert_eq!(x["capabilities"]["simulation_seed"], false);
    assert_eq!(
        x["record"]["duplex_repair"]["overlaps"][0]["query_position"],
        2
    );
}

#[test]
fn snapshot_and_certificate_tampering_are_errors_not_partial_annotations() {
    let r = data("AAATGT", "CCC", "TTTAAA", Some(3));
    let b = DuplexReferenceBindings::new(
        &r,
        ["V*01", "D*01", "J*01"],
        ["AAATGT", "CCC", "TTTAAA"],
        false,
        true,
    )
    .unwrap();
    let mut p = repair_coupled_annotated(
        "AAATGT",
        "CCC",
        "TTTAAA",
        false,
        &[blunt(JoinName::VD, 6), blunt(JoinName::DJ, 3)],
        256,
        1000,
    )
    .unwrap();
    p.evidence.placements.get_mut(&0).unwrap().body = None;
    assert!(project_coupled(&p, &b, &r)
        .unwrap_err()
        .contains("exact replay"));
}
