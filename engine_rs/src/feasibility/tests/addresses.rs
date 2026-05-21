use super::*;

#[test]
fn is_vj_productive_choice_accepts_the_six_known_addresses() {
    for addr in [
        "sample_allele.v",
        "sample_allele.j",
        "trim.v_5",
        "trim.v_3",
        "trim.j_5",
        "trim.j_3",
    ] {
        assert!(is_vj_productive_choice(addr), "expected {addr:?} to count");
    }
}

#[test]
fn is_vj_productive_choice_rejects_unrelated_addresses() {
    for addr in [
        "sample_allele.d",
        "trim.d_5",
        "trim.d_3",
        "np.np1.length",
        "mutate.s5f.site[0]",
        "",
        "sample_allele.v.something",
    ] {
        assert!(
            !is_vj_productive_choice(addr),
            "expected {addr:?} to NOT count"
        );
    }
}

#[test]
fn trim_address_maps_each_supported_segment_end_pair() {
    assert_eq!(trim_address(Segment::V, TrimEnd::Five), "trim.v_5");
    assert_eq!(trim_address(Segment::V, TrimEnd::Three), "trim.v_3");
    assert_eq!(trim_address(Segment::J, TrimEnd::Five), "trim.j_5");
    assert_eq!(trim_address(Segment::J, TrimEnd::Three), "trim.j_3");
}

#[test]
fn trim_address_returns_unsupported_marker_for_other_segments() {
    assert_eq!(
        trim_address(Segment::D, TrimEnd::Five),
        "trim.<unsupported>"
    );
    assert_eq!(
        trim_address(Segment::Np1, TrimEnd::Three),
        "trim.<unsupported>"
    );
}

#[test]
fn sample_allele_address_maps_v_and_j_only() {
    assert_eq!(sample_allele_address(Segment::V), "sample_allele.v");
    assert_eq!(sample_allele_address(Segment::J), "sample_allele.j");
    assert_eq!(
        sample_allele_address(Segment::D),
        "sample_allele.<unsupported>"
    );
}
