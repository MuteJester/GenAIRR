use super::*;
use crate::duplex_repair::{JoinName, RepairOperation, Side};
use crate::refdata::{Allele, ChainType};

pub(crate) fn fixture() -> (RefDataConfig, FiniteDuplexModel) {
    let mut refdata = RefDataConfig::empty(ChainType::Vdj);
    for (segment, name, sequence, anchor) in [
        (Segment::V, "V*01", "AAATGT", Some(3)),
        (Segment::D, "D*01", "ACT", None),
        (Segment::J, "J*01", "TTTAAA", Some(0)),
    ] {
        let allele = Allele {
            name: name.into(),
            gene: name.split('*').next().unwrap().into(),
            seq: sequence.as_bytes().to_vec(),
            segment,
            anchor,
            functional_status: None,
            subregions: vec![],
        };
        match segment {
            Segment::V => refdata.v_pool.push(allele),
            Segment::D => refdata.d_pool.push(allele),
            Segment::J => refdata.j_pool.push(allele),
            _ => unreachable!(),
        };
    }
    let script = CoupledDuplexScript {
        model: COUPLED_MODEL.into(),
        v: "AAATGT".into(),
        d: "ACT".into(),
        j: "TTTAAA".into(),
        d_inverted: false,
        max_operations: 16,
        max_bases: 100,
        joins: [(JoinName::VD, 6), (JoinName::DJ, 3)]
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
            .collect(),
    };
    let rows = [("a", 1.0), ("b", 3.0), ("disabled", 0.0)]
        .into_iter()
        .map(|(id, weight)| DuplexHistoryRow {
            id: id.into(),
            weight,
            v_call: "V*01".into(),
            d_call: "D*01".into(),
            j_call: "J*01".into(),
            script: script.clone(),
        })
        .collect();
    let model = FiniteDuplexModel {
        model: FINITE_DUPLEX_MODEL.into(),
        reference_content_hash: crate::trace_file::refdata_content_hash(&refdata),
        reference_ends_are_coding_ends: true,
        max_preflight_operations: 1000,
        max_preflight_bases: 1000,
        rows,
    };
    (refdata, model)
}

#[test]
fn model_round_trip_and_digest_bind_all_authored_history_inputs() {
    let (_, model) = fixture();
    model.validate_schema().unwrap();
    let json = serde_json::to_string(&model).unwrap();
    let round_trip: FiniteDuplexModel = serde_json::from_str(&json).unwrap();
    assert_eq!(model, round_trip);
    assert_eq!(model.digest(), round_trip.digest());
    let mut variants = Vec::new();
    let mut v = model.clone();
    v.rows.swap(0, 1);
    variants.push(v);
    let mut v = model.clone();
    v.rows[0].id = "different".into();
    variants.push(v);
    let mut v = model.clone();
    v.rows[0].weight = 2.0;
    variants.push(v);
    let mut v = model.clone();
    v.rows[0].script.max_operations += 1;
    variants.push(v);
    let mut v = model.clone();
    v.rows[0].script.joins.reverse();
    variants.push(v);
    let mut v = model.clone();
    v.max_preflight_bases += 1;
    variants.push(v);
    for variant in variants {
        assert_ne!(model.digest(), variant.digest());
    }
    let mut value = serde_json::to_value(&model).unwrap();
    value["rows"][0]["unexpected"] = true.into();
    assert!(serde_json::from_value::<FiniteDuplexModel>(value).is_err());
}

#[test]
fn schema_rejects_invalid_weights_ids_hash_and_version() {
    let (_, model) = fixture();
    for weight in [-1.0, f64::NAN, f64::INFINITY] {
        let mut invalid = model.clone();
        invalid.rows[0].weight = weight;
        assert!(invalid.validate_schema().is_err());
    }
    let mut invalid = model.clone();
    invalid.rows[1].id = invalid.rows[0].id.clone();
    assert!(invalid.validate_schema().is_err());
    let mut invalid = model.clone();
    invalid.reference_content_hash = "sha256:fake".into();
    assert!(invalid.validate_schema().is_err());
    let mut invalid = model.clone();
    invalid.model = "unknown".into();
    assert!(invalid.validate_schema().is_err());
}

#[test]
fn normalization_never_silently_erases_positive_authored_rows() {
    let (_, mut model) = fixture();
    model.rows[0].weight = f64::from_bits(1);
    model.rows[1].weight = 1e308;
    assert!(model
        .validate_weight_support(&[0, 1])
        .unwrap_err()
        .contains("numerical"));
    // A sole tiny row is still fully representable after support filtering.
    model.validate_weight_support(&[0]).unwrap();
    model.rows[0].weight = 1e308;
    model.validate_weight_support(&[0, 1]).unwrap();
    model.rows[0].weight = 1.0;
    model.rows[1].weight = 1e-20;
    assert!(model
        .validate_weight_support(&[0, 1])
        .unwrap_err()
        .contains("numerical"));
    model.rows[1].weight = 3.0;
    model.validate_weight_support(&[0, 1]).unwrap();
}

#[test]
fn cached_digest_preparation_preserves_the_verified_public_result() {
    let (refdata, model) = fixture();
    let digest = model.digest();
    for index in [0, 1] {
        let direct = model.prepare_row(index, &refdata).unwrap();
        let cached = model
            .prepare_row_with_digest(index, &refdata, &digest)
            .unwrap();
        assert_eq!(direct.context, cached.context);
        assert_eq!(direct.seed.pool.as_slice(), cached.seed.pool.as_slice());
        assert_eq!(direct.seed.duplex_provenance, cached.seed.duplex_provenance);
    }
}

#[test]
fn preparation_preserves_distinct_history_identity_and_rejects_reference_mismatch() {
    let (refdata, model) = fixture();
    let a = model.prepare_row(0, &refdata).unwrap();
    let b = model.prepare_row(1, &refdata).unwrap();
    assert_eq!(a.seed.pool.as_slice(), b.seed.pool.as_slice());
    let pa = a.seed.duplex_provenance.as_ref().unwrap();
    let pb = b.seed.duplex_provenance.as_ref().unwrap();
    assert_eq!(pa.founder.annotated, pb.founder.annotated);
    assert_ne!(pa.founder.selection, pb.founder.selection);
    assert_eq!(
        serde_json::to_value(&pa.founder.selection).unwrap()["natural_weight"],
        1.0
    );
    let mut wrong = refdata.clone();
    wrong.identity.name = Some("different cartridge identity".into());
    assert!(model.prepare_row(0, &wrong).is_err());
    let mut bad = model.clone();
    bad.rows[0].script.v = "AAATGC".into();
    assert!(bad.prepare_row(0, &refdata).is_err());
}
