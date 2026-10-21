use super::*;
use std::collections::BTreeSet;

fn join(join: JoinName, offset: usize, right_start: i64) -> JoinScript {
    JoinScript {
        join,
        operations: vec![
            RepairOperation::Open {
                side: Side::Left,
                offset,
            },
            RepairOperation::Open {
                side: Side::Right,
                offset,
            },
            RepairOperation::Synapse { right_start },
            RepairOperation::Ligate {},
        ],
    }
}

fn run(v: &str, d: &str, j: &str, inverted: bool, joins: &[JoinScript]) -> CoupledProduct {
    repair_coupled(v, d, j, inverted, joins, 64, 4096).unwrap()
}

fn check_ledger(product: &CoupledProduct) {
    let mut surviving = BTreeSet::new();
    for column in &product.columns {
        assert_eq!(complement(column.top.base), column.bottom.base);
        for base in [&column.top, &column.bottom] {
            assert!(surviving.insert(base.id), "duplicated physical nucleotide");
            assert_eq!(&product.ledger[base.id as usize], base);
        }
    }
    for base in &product.removed {
        assert!(
            surviving.insert(base.id),
            "removed nucleotide remains in product"
        );
        assert_eq!(&product.ledger[base.id as usize], base);
    }
    assert_eq!(surviving, (0..product.ledger.len() as u64).collect());
    for base in &product.ledger {
        if let Some(parent) = base.template_id {
            assert!(parent < base.id);
            assert_eq!(complement(product.ledger[parent as usize].base), base.base);
        }
    }
}

#[test]
fn blunt_both_orders_and_inversion_preserve_canonical_ids() {
    for inverted in [false, true] {
        let vd = join(JoinName::VD, 0, 3);
        let dj = join(JoinName::DJ, 0, 3);
        let a = run("AAC", "GAT", "TCC", inverted, &[vd.clone(), dj.clone()]);
        let b = run("AAC", "GAT", "TCC", inverted, &[dj, vd]);
        assert_eq!(a.sequence, if inverted { "AACATCTCC" } else { "AACGATTCC" });
        assert_eq!(a.columns, b.columns);
        assert_eq!(a.ledger, b.ledger);
        assert_eq!(
            a.origin_positions,
            BTreeMap::from([(Substrate::V, 0), (Substrate::D, 3), (Substrate::J, 6),])
        );
        for p in 0..3 {
            let base = &a.columns[3 + p].top;
            let original = base.original.as_ref().unwrap();
            assert_eq!(original.segment, Substrate::D);
            assert_eq!(
                original.reference_position,
                if inverted { 2 - p } else { p }
            );
            assert_eq!(
                original.strand,
                if inverted {
                    Strand::Bottom
                } else {
                    Strand::Top
                }
            );
            let expected_id = 6 + 2 * original.reference_position + usize::from(inverted);
            assert_eq!(base.id, expected_id as u64);
        }
        check_ledger(&a);
        check_ledger(&b);
    }
}

#[test]
fn successive_short_d_opening_uses_actual_component_stem() {
    let vd = join(JoinName::VD, 1, 3);
    let dj = join(JoinName::DJ, 1, 2);
    for history in [vec![vd.clone(), dj.clone()], vec![dj, vd.clone()]] {
        let product = run("AAG", "CA", "TAA", false, &history);
        assert_eq!(product.sequence, "AAGCATAA");
        assert_eq!(product.ledger.len(), 16);
        assert!(product.removed.is_empty());
        check_ledger(&product);
    }
    // The second opening may extend beyond the original D-length heuristic:
    // its actual paired stem includes material contributed at the first join.
    let product = run("AAG", "CA", "TGA", false, &[vd, join(JoinName::DJ, 2, 2)]);
    assert_eq!(product.sequence, "AAGCATGA");
    assert_eq!(product.ligation_seams[0].top, 4);
    assert_eq!(product.ligation_seams[0].bottom, 2);
    check_ledger(&product);
}

fn restored_d_then_vd() -> Vec<JoinScript> {
    let mut dj = join(JoinName::DJ, 0, 3);
    dj.operations.insert(
        2,
        RepairOperation::Resect {
            side: Side::Left,
            strand: Strand::Top,
            count: 2,
        },
    );
    dj.operations.insert(
        4,
        RepairOperation::Fill {
            strand: Strand::Top,
        },
    );
    vec![dj, join(JoinName::VD, 2, 3)]
}

#[test]
fn prior_synthesis_moves_across_second_hairpin_with_template_links_intact() {
    let product = run("ATG", "CAA", "TAA", false, &restored_d_then_vd());
    assert_eq!(product.sequence, "ATGCAATAA");
    assert_eq!(product.ledger.len(), 20);
    assert_eq!(
        product.removed.iter().map(|b| b.id).collect::<Vec<_>>(),
        vec![10, 8]
    );
    let moved = &product.columns[1].bottom;
    assert_eq!(moved.id, 18);
    assert_eq!(moved.kind, "templated_synthesis");
    assert_eq!(moved.template_id, Some(9));
    assert_eq!(moved.original, None);
    assert_eq!(product.columns[5].top.id, 19);
    check_ledger(&product);
}

#[test]
fn global_ledger_retains_a_removed_template_of_prior_and_successive_synthesis() {
    let mut history = restored_d_then_vd();
    let mut vd = join(JoinName::VD, 0, 3);
    vd.operations.insert(
        2,
        RepairOperation::Resect {
            side: Side::Right,
            strand: Strand::Bottom,
            count: 2,
        },
    );
    vd.operations.insert(
        4,
        RepairOperation::Fill {
            strand: Strand::Bottom,
        },
    );
    history[1] = vd;
    let product = run("ACG", "CAA", "TAA", false, &history);
    assert_eq!(product.sequence, "ACGCAATAA");
    assert_eq!(product.ledger.len(), 22);
    assert_eq!(product.removed.len(), 4);
    assert_eq!(product.columns[4].top.id, 18);
    assert_eq!(product.columns[4].top.template_id, Some(9));
    assert!(product.removed.iter().any(|b| b.id == 9));
    assert_eq!(product.columns[4].bottom.template_id, Some(18));
    check_ledger(&product);
}

#[test]
fn translated_first_join_seams_stay_anchored_after_dj_then_vd() {
    let product = run(
        "AAA",
        "CCC",
        "GGG",
        false,
        &[join(JoinName::DJ, 0, 3), join(JoinName::VD, 0, 3)],
    );
    assert_eq!(
        product.ligation_seams,
        vec![
            LigationSeam {
                join: JoinName::DJ,
                top: 6,
                bottom: 6
            },
            LigationSeam {
                join: JoinName::VD,
                top: 3,
                bottom: 3
            },
        ]
    );
}

#[test]
fn second_resection_and_opening_cannot_break_a_prior_ligation_bond() {
    for (first, second, side, strand) in [
        (JoinName::VD, JoinName::DJ, Side::Left, Strand::Top),
        (JoinName::VD, JoinName::DJ, Side::Left, Strand::Bottom),
        (JoinName::DJ, JoinName::VD, Side::Right, Strand::Top),
        (JoinName::DJ, JoinName::VD, Side::Right, Strand::Bottom),
    ] {
        let mut later = join(second, 0, 3);
        later.operations.insert(
            2,
            RepairOperation::Resect {
                side,
                strand,
                count: 3,
            },
        );
        let error = repair_coupled(
            "AAA",
            "CCC",
            "GGG",
            false,
            &[join(first, 0, 3), later],
            64,
            4096,
        )
        .unwrap_err();
        assert!(error.contains("previous ligation seam"), "{error}");
    }
    // After a one-base opening of D's left end, opening the other at 3 would
    // cut the completed bottom-strand seam rather than an unprocessed hairpin.
    let error = repair_coupled(
        "AAG",
        "CA",
        "TGA",
        false,
        &[join(JoinName::VD, 1, 3), join(JoinName::DJ, 3, 2)],
        64,
        4096,
    )
    .unwrap_err();
    assert!(error.contains("previous ligation seam"), "{error}");
}

#[test]
fn global_budgets_count_original_d_once_and_never_refund_resection() {
    let history = [join(JoinName::VD, 0, 3), join(JoinName::DJ, 0, 3)];
    assert!(repair_coupled("AAA", "CCC", "GGG", false, &history, 8, 18).is_ok());
    assert!(repair_coupled("AAA", "CCC", "GGG", false, &history, 7, 18)
        .unwrap_err()
        .contains("operation budget"));
    assert!(repair_coupled("AAA", "CCC", "GGG", false, &history, 8, 17)
        .unwrap_err()
        .contains("base budget"));
    let history = restored_d_then_vd();
    assert!(repair_coupled("ATG", "CAA", "TAA", false, &history, 10, 19)
        .unwrap_err()
        .contains("base budget"));
    assert!(repair_coupled("ATG", "CAA", "TAA", false, &history, 10, 20).is_ok());
}

#[test]
fn script_replay_is_exact_and_invalid_join_sets_are_rejected() {
    let history = restored_d_then_vd();
    let json = serde_json::to_string(&history).unwrap();
    let replay: Vec<JoinScript> = serde_json::from_str(&json).unwrap();
    assert_eq!(
        run("ATG", "CAA", "TAA", false, &history),
        run("ATG", "CAA", "TAA", false, &replay)
    );
    assert!(json.contains("\"join\":\"dj\""));
    for history in [
        vec![],
        vec![join(JoinName::VD, 0, 3)],
        vec![join(JoinName::DJ, 0, 3), join(JoinName::DJ, 0, 3)],
    ] {
        assert!(
            repair_coupled("AAA", "CCC", "GGG", false, &history, 64, 4096)
                .unwrap_err()
                .contains("exactly one VD and one DJ")
        );
    }
    assert!(
        serde_json::from_str::<JoinScript>(r#"{"join":"vd","operations":[],"invented":true}"#)
            .is_err()
    );
}

#[test]
fn huge_opening_or_synapsis_coordinates_are_bounded_errors() {
    for history in [
        vec![join(JoinName::VD, usize::MAX, 3), join(JoinName::DJ, 0, 3)],
        vec![join(JoinName::VD, 0, i64::MAX), join(JoinName::DJ, 0, 3)],
        vec![join(JoinName::VD, 0, 3), join(JoinName::DJ, 0, i64::MIN)],
    ] {
        assert!(repair_coupled("AAA", "CCC", "GGG", false, &history, 64, 4096).is_err());
    }
}
