//! Canonical trace/pass address strings.
//!
//! Trace addresses are part of the engine's internal contract: passes
//! write them, contracts and feasibility filters inspect them, AIRR
//! projection reads them, and tests assert them. Keep the spelling here
//! as the single source of truth.

use crate::assignment::TrimEnd;
use crate::ir::Segment;

pub const SAMPLE_ALLELE_V: &str = "sample_allele.v";
pub const SAMPLE_ALLELE_D: &str = "sample_allele.d";
pub const SAMPLE_ALLELE_J: &str = "sample_allele.j";
pub const SAMPLE_ALLELE_INVALID: &str = "sample_allele.<invalid>";
pub const SAMPLE_ALLELE_UNSUPPORTED: &str = "sample_allele.<unsupported>";

pub const TRIM_V_5: &str = "trim.v_5";
pub const TRIM_V_3: &str = "trim.v_3";
pub const TRIM_D_5: &str = "trim.d_5";
pub const TRIM_D_3: &str = "trim.d_3";
pub const TRIM_J_5: &str = "trim.j_5";
pub const TRIM_J_3: &str = "trim.j_3";
pub const TRIM_INVALID: &str = "trim.<invalid>";
pub const TRIM_UNSUPPORTED: &str = "trim.<unsupported>";

pub const ASSEMBLE_V: &str = "assemble.v";
pub const ASSEMBLE_D: &str = "assemble.d";
pub const ASSEMBLE_J: &str = "assemble.j";

pub const GENERATE_NP1: &str = "generate_np.np1";
pub const GENERATE_NP2: &str = "generate_np.np2";

pub const NP1_LENGTH: &str = "np.np1.length";
pub const NP2_LENGTH: &str = "np.np2.length";
pub const NP_INVALID_LENGTH: &str = "np.<invalid>.length";
pub const NP1_BASES_PREFIX: &str = "np.np1.bases";
pub const NP2_BASES_PREFIX: &str = "np.np2.bases";
pub const NP1_BASES_INDEX_PREFIX: &str = "np.np1.bases[";
pub const NP2_BASES_INDEX_PREFIX: &str = "np.np2.bases[";

pub const MUTATE_UNIFORM: &str = "mutate.uniform";
pub const MUTATE_UNIFORM_COUNT: &str = "mutate.uniform.count";
pub const MUTATE_UNIFORM_SITE_PATTERN: &str = "mutate.uniform.site[0..n]";
pub const MUTATE_UNIFORM_BASE_PATTERN: &str = "mutate.uniform.base[0..n]";
pub const MUTATE_UNIFORM_SITE_PREFIX: &str = "mutate.uniform.site[";
pub const MUTATE_UNIFORM_BASE_PREFIX: &str = "mutate.uniform.base[";

pub const MUTATE_S5F: &str = "mutate.s5f";
pub const MUTATE_S5F_COUNT: &str = "mutate.s5f.count";
pub const MUTATE_S5F_SITE_PATTERN: &str = "mutate.s5f.site[0..n]";
pub const MUTATE_S5F_BASE_PATTERN: &str = "mutate.s5f.base[0..n]";
pub const MUTATE_S5F_SITE_PREFIX: &str = "mutate.s5f.site[";
pub const MUTATE_S5F_BASE_PREFIX: &str = "mutate.s5f.base[";

pub const CORRUPT_PCR: &str = "corrupt.pcr";
pub const CORRUPT_PCR_COUNT: &str = "corrupt.pcr.count";
pub const CORRUPT_PCR_SITE_PATTERN: &str = "corrupt.pcr.error_site[0..n]";
pub const CORRUPT_PCR_BASE_PATTERN: &str = "corrupt.pcr.error_base[0..n]";
pub const CORRUPT_PCR_SITE_PREFIX: &str = "corrupt.pcr.error_site[";
pub const CORRUPT_PCR_BASE_PREFIX: &str = "corrupt.pcr.error_base[";

pub const CORRUPT_QUALITY: &str = "corrupt.quality";
pub const CORRUPT_QUALITY_COUNT: &str = "corrupt.quality.count";
pub const CORRUPT_QUALITY_SITE_PATTERN: &str = "corrupt.quality.error_site[0..n]";
pub const CORRUPT_QUALITY_BASE_PATTERN: &str = "corrupt.quality.error_base[0..n]";
pub const CORRUPT_QUALITY_SITE_PREFIX: &str = "corrupt.quality.error_site[";
pub const CORRUPT_QUALITY_BASE_PREFIX: &str = "corrupt.quality.error_base[";

pub const CORRUPT_CONTAMINANT: &str = "corrupt.contaminant";
pub const CORRUPT_CONTAMINANT_APPLIED: &str = "corrupt.contaminant.applied";
pub const CORRUPT_CONTAMINANT_BASES_PATTERN: &str = "corrupt.contaminant.bases[0..n]";
pub const CORRUPT_CONTAMINANT_BASES_PREFIX: &str = "corrupt.contaminant.bases[";

pub const CORRUPT_INDEL: &str = "corrupt.indel";
pub const CORRUPT_INDEL_COUNT: &str = "corrupt.indel.count";
pub const CORRUPT_INDEL_KIND_PATTERN: &str = "corrupt.indel.kind[0..n]";
pub const CORRUPT_INDEL_SITE_PATTERN: &str = "corrupt.indel.site[0..n]";
pub const CORRUPT_INDEL_BASE_PATTERN: &str = "corrupt.indel.base[0..n]";
pub const CORRUPT_INDEL_KIND_PREFIX: &str = "corrupt.indel.kind[";
pub const CORRUPT_INDEL_SITE_PREFIX: &str = "corrupt.indel.site[";
pub const CORRUPT_INDEL_BASE_PREFIX: &str = "corrupt.indel.base[";

pub const CORRUPT_NS: &str = "corrupt.ns";
pub const CORRUPT_NS_COUNT: &str = "corrupt.ns.count";
pub const CORRUPT_NS_SITE_PATTERN: &str = "corrupt.ns.site[0..n]";
pub const CORRUPT_NS_SITE_PREFIX: &str = "corrupt.ns.site[";

pub const CORRUPT_END_LOSS_5: &str = "corrupt.end_loss.5";
pub const CORRUPT_END_LOSS_3: &str = "corrupt.end_loss.3";

pub const CORRUPT_REV_COMP: &str = "corrupt.rev_comp";
pub const CORRUPT_REV_COMP_APPLIED: &str = "corrupt.rev_comp.applied";

pub fn sample_allele(segment: Segment) -> Option<&'static str> {
    match segment {
        Segment::V => Some(SAMPLE_ALLELE_V),
        Segment::D => Some(SAMPLE_ALLELE_D),
        Segment::J => Some(SAMPLE_ALLELE_J),
        Segment::Np1 | Segment::Np2 => None,
    }
}

pub fn sample_allele_vdj(segment: Segment) -> &'static str {
    sample_allele(segment).expect("sample allele address is only defined for V/D/J")
}

pub fn trim(segment: Segment, end: TrimEnd) -> Option<&'static str> {
    match (segment, end) {
        (Segment::V, TrimEnd::Five) => Some(TRIM_V_5),
        (Segment::V, TrimEnd::Three) => Some(TRIM_V_3),
        (Segment::D, TrimEnd::Five) => Some(TRIM_D_5),
        (Segment::D, TrimEnd::Three) => Some(TRIM_D_3),
        (Segment::J, TrimEnd::Five) => Some(TRIM_J_5),
        (Segment::J, TrimEnd::Three) => Some(TRIM_J_3),
        _ => None,
    }
}

pub fn trim_vdj(segment: Segment, end: TrimEnd) -> &'static str {
    trim(segment, end).expect("trim address is only defined for V/D/J")
}

pub fn assemble(segment: Segment) -> Option<&'static str> {
    match segment {
        Segment::V => Some(ASSEMBLE_V),
        Segment::D => Some(ASSEMBLE_D),
        Segment::J => Some(ASSEMBLE_J),
        Segment::Np1 | Segment::Np2 => None,
    }
}

pub fn assemble_vdj(segment: Segment) -> &'static str {
    assemble(segment).expect("assemble address is only defined for V/D/J")
}

pub fn generate_np(segment: Segment) -> Option<&'static str> {
    match segment {
        Segment::Np1 => Some(GENERATE_NP1),
        Segment::Np2 => Some(GENERATE_NP2),
        _ => None,
    }
}

pub fn generate_np_region(segment: Segment) -> &'static str {
    generate_np(segment).expect("generate_np address is only defined for NP1/NP2")
}

pub fn np_length(segment: Segment) -> Option<&'static str> {
    match segment {
        Segment::Np1 => Some(NP1_LENGTH),
        Segment::Np2 => Some(NP2_LENGTH),
        _ => None,
    }
}

pub fn np_length_region(segment: Segment) -> &'static str {
    np_length(segment).expect("NP length address is only defined for NP1/NP2")
}

pub fn np_bases_prefix(segment: Segment) -> Option<&'static str> {
    match segment {
        Segment::Np1 => Some(NP1_BASES_PREFIX),
        Segment::Np2 => Some(NP2_BASES_PREFIX),
        _ => None,
    }
}

pub fn np_bases_region_prefix(segment: Segment) -> &'static str {
    np_bases_prefix(segment).expect("NP bases address is only defined for NP1/NP2")
}

pub fn np_base(segment: Segment, index: u32) -> Option<String> {
    Some(format!("{}[{}]", np_bases_prefix(segment)?, index))
}

pub fn np_bases_pattern(segment: Segment) -> Option<String> {
    Some(format!("{}[0..n]", np_bases_prefix(segment)?))
}

pub fn mutate_uniform_site(index: u32) -> String {
    indexed(MUTATE_UNIFORM_SITE_PREFIX, index)
}

pub fn mutate_uniform_base(index: u32) -> String {
    indexed(MUTATE_UNIFORM_BASE_PREFIX, index)
}

pub fn mutate_s5f_site(index: u32) -> String {
    indexed(MUTATE_S5F_SITE_PREFIX, index)
}

pub fn mutate_s5f_base(index: u32) -> String {
    indexed(MUTATE_S5F_BASE_PREFIX, index)
}

pub fn corrupt_pcr_site(index: u32) -> String {
    indexed(CORRUPT_PCR_SITE_PREFIX, index)
}

pub fn corrupt_pcr_base(index: u32) -> String {
    indexed(CORRUPT_PCR_BASE_PREFIX, index)
}

pub fn corrupt_quality_site(index: u32) -> String {
    indexed(CORRUPT_QUALITY_SITE_PREFIX, index)
}

pub fn corrupt_quality_base(index: u32) -> String {
    indexed(CORRUPT_QUALITY_BASE_PREFIX, index)
}

pub fn corrupt_contaminant_base(index: u32) -> String {
    indexed(CORRUPT_CONTAMINANT_BASES_PREFIX, index)
}

pub fn corrupt_indel_kind(index: u32) -> String {
    indexed(CORRUPT_INDEL_KIND_PREFIX, index)
}

pub fn corrupt_indel_site(index: u32) -> String {
    indexed(CORRUPT_INDEL_SITE_PREFIX, index)
}

pub fn corrupt_indel_base(index: u32) -> String {
    indexed(CORRUPT_INDEL_BASE_PREFIX, index)
}

pub fn corrupt_ns_site(index: u32) -> String {
    indexed(CORRUPT_NS_SITE_PREFIX, index)
}

pub fn parse_np_base(address: &str) -> Option<(Segment, u32)> {
    if let Some(index) = parse_indexed(address, NP1_BASES_INDEX_PREFIX) {
        Some((Segment::Np1, index))
    } else {
        parse_indexed(address, NP2_BASES_INDEX_PREFIX).map(|index| (Segment::Np2, index))
    }
}

pub fn parse_np_length(address: &str) -> Option<Segment> {
    match address {
        NP1_LENGTH => Some(Segment::Np1),
        NP2_LENGTH => Some(Segment::Np2),
        _ => None,
    }
}

pub fn parse_indexed(address: &str, prefix: &str) -> Option<u32> {
    let rest = address.strip_prefix(prefix)?;
    rest.strip_suffix(']')?.parse::<u32>().ok()
}

fn indexed(prefix: &str, index: u32) -> String {
    format!("{prefix}{index}]")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn segment_helpers_emit_existing_recombination_addresses() {
        assert_eq!(sample_allele_vdj(Segment::V), "sample_allele.v");
        assert_eq!(sample_allele_vdj(Segment::D), "sample_allele.d");
        assert_eq!(sample_allele_vdj(Segment::J), "sample_allele.j");
        assert_eq!(trim_vdj(Segment::V, TrimEnd::Five), "trim.v_5");
        assert_eq!(trim_vdj(Segment::J, TrimEnd::Three), "trim.j_3");
        assert_eq!(assemble_vdj(Segment::D), "assemble.d");
        assert_eq!(np_length_region(Segment::Np1), "np.np1.length");
        assert_eq!(np_base(Segment::Np2, 3).unwrap(), "np.np2.bases[3]");
    }

    #[test]
    fn parsers_accept_existing_indexed_addresses() {
        assert_eq!(parse_np_base("np.np1.bases[17]"), Some((Segment::Np1, 17)));
        assert_eq!(parse_np_base("np.np2.bases[4]"), Some((Segment::Np2, 4)));
        assert_eq!(parse_np_base("np.np2.bases[x]"), None);
        assert_eq!(parse_np_length("np.np1.length"), Some(Segment::Np1));
        assert_eq!(
            parse_indexed("mutate.s5f.base[2]", MUTATE_S5F_BASE_PREFIX),
            Some(2)
        );
    }
}
