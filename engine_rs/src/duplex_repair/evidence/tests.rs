use super::*;

fn opened(join: JoinName) -> JoinScript {
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
    let mut s = opened(join);
    s.operations.extend([
        RepairOperation::Synapse { right_start: at },
        RepairOperation::Ligate {},
    ]);
    s
}

#[test]
fn evidence_does_not_change_physical_products_in_either_order_or_orientation() {
    for inverse in [false, true] {
        for reverse in [false, true] {
            let mut joins = vec![blunt(JoinName::VD, 3), blunt(JoinName::DJ, 3)];
            if reverse {
                joins.reverse();
            }
            let p =
                repair_coupled_annotated("ACG", "AAT", "TCC", inverse, &joins, 256, 1000).unwrap();
            assert_eq!(
                p.product,
                repair_coupled("ACG", "AAT", "TCC", inverse, &joins, 256, 1000).unwrap()
            );
            assert_eq!(p.evidence.placements.len(), 18);
            for column in &p.product.columns {
                for (strand, base) in [(Strand::Top, &column.top), (Strand::Bottom, &column.bottom)]
                {
                    let placement = &p.evidence.placements[&base.id];
                    assert_eq!(
                        placement.location,
                        Some(PhysicalLocation {
                            frame: Substrate::V,
                            strand,
                            position: column.position
                        })
                    );
                    assert_eq!(placement.body, base.original);
                }
            }
        }
    }
}

#[test]
fn restoration_and_resection_keep_different_chemical_and_body_facts() {
    let mut vd = opened(JoinName::VD);
    vd.operations.extend([
        RepairOperation::Resect {
            side: Side::Right,
            strand: Strand::Top,
            count: 1,
        },
        RepairOperation::Synapse { right_start: 3 },
        RepairOperation::Fill {
            strand: Strand::Top,
        },
        RepairOperation::Ligate {},
    ]);
    let p = repair_coupled_annotated(
        "ACG",
        "CAT",
        "TAA",
        false,
        &[vd, blunt(JoinName::DJ, 3)],
        256,
        1000,
    )
    .unwrap();
    let restored = &p.product.columns[3].top;
    assert_eq!(restored.kind, "templated_synthesis");
    let placement = &p.evidence.placements[&restored.id];
    assert_eq!(placement.role, PlacementRole::RestoredBody);
    assert_eq!(placement.body.as_ref().unwrap().segment, Substrate::D);
    assert_eq!(placement.body.as_ref().unwrap().reference_position, 0);
    assert_eq!(p.product.removed.len(), 1);
    let removed = &p.evidence.placements[&p.product.removed[0].id];
    assert!(removed.body.is_none() && removed.location.is_none());
    assert!(removed.removal_operation.is_some());
}

#[test]
fn hairpin_transfer_revokes_body_without_losing_original_identity() {
    let scripts: Vec<_> = [(JoinName::VD, 3), (JoinName::DJ, 2)]
        .into_iter()
        .map(|(join, at)| JoinScript {
            join,
            operations: vec![
                RepairOperation::Open {
                    side: Side::Left,
                    offset: 1,
                },
                RepairOperation::Open {
                    side: Side::Right,
                    offset: 1,
                },
                RepairOperation::Synapse { right_start: at },
                RepairOperation::Ligate {},
            ],
        })
        .collect();
    let p = repair_coupled_annotated("AAG", "CA", "TAA", false, &scripts, 256, 1000).unwrap();
    assert_eq!(p.product.sequence, "AAGCATAA");
    let top = &p.product.columns[3].top;
    assert_eq!(top.original.as_ref().unwrap().segment, Substrate::V);
    assert_eq!(
        p.evidence.placements[&top.id].role,
        PlacementRole::HairpinTransfer
    );
    assert!(p.evidence.placements[&top.id].body.is_none());
    let opposite = &p.product.columns[3].bottom;
    assert_eq!(
        p.evidence.placements[&opposite.id]
            .body
            .as_ref()
            .unwrap()
            .segment,
        Substrate::D
    );
    assert_eq!(
        p.evidence
            .events
            .iter()
            .filter(|e| matches!(e, RepairEvidenceEvent::Transferred { .. }))
            .count(),
        4
    );
}
