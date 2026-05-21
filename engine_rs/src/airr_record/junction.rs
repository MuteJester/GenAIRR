use crate::codon::translate_codon_slice;
use crate::ir::{NucHandle, Region, Segment, Simulation};
use crate::refdata::RefDataConfig;

use super::projection::lookup_allele;

/// locate the pool position where an allele's anchor
/// codon resides in the assembled sequence. Scans the segment's
/// structural region for the first node whose `germline_pos`
/// matches the anchor's allele coordinate. Returns `None` when
/// the anchor was trimmed off, deleted by a structural-indel
/// pass, or otherwise not present in the live data.
///
/// Biologically correct under indels: each surviving germline
/// node carries its original `germline_pos`, so the anchor codon
/// can be found wherever the actual base ended up - even after
/// insertions or deletions shifted things around inside the
/// segment.
pub(super) fn anchor_pool_position(
    sim: &Simulation,
    region: &Region,
    anchor_ref: u32,
) -> Option<u32> {
    let r_start = region.start.index();
    let r_end = region.end.index();
    let pool = sim.pool.as_slice();
    for i in r_start..r_end {
        let nuc = &pool[i as usize];
        if nuc.germline_pos.get().map(|p| p as u32) == Some(anchor_ref) {
            return Some(i);
        }
    }
    None
}

pub(super) fn anchor_amino_acid_preserved(
    sim: &Simulation,
    refdata: &RefDataConfig,
    segment: Segment,
    region: &Region,
    allele_id: Option<crate::refdata::AlleleId>,
    _trim_5: i64,
) -> bool {
    let Some(allele) = lookup_allele(refdata, segment, allele_id) else {
        return true;
    };
    let Some(anchor) = allele.anchor else {
        return true;
    };

    let anchor = anchor as i64;
    if anchor < 0 || anchor + 3 > allele.seq.len() as i64 {
        return false;
    }

    // locate the anchor by `germline_pos`, not by a
    // structural offset. Same rationale as the junction-window
    // scanner - a structural offset miscomputes the anchor pool
    // position when indels disturb V/J's coding region.
    let Some(pool_start) = anchor_pool_position(sim, region, anchor as u32) else {
        return false;
    };
    if pool_start + 3 > region.end.index() {
        return false;
    }

    let mut live = [b'N'; 3];
    for offset in 0..3 {
        let Some(nuc) = sim.pool.get(NucHandle::new(pool_start + offset)) else {
            return false;
        };
        live[offset as usize] = nuc.base;
    }

    let ref_start = anchor as usize;
    let reference = [
        allele.seq[ref_start],
        allele.seq[ref_start + 1],
        allele.seq[ref_start + 2],
    ];
    translate_codon_slice(&live) == translate_codon_slice(&reference)
}

pub(super) fn aa_slice_for_region(
    region_start: usize,
    region_end: usize,
    frame_offset: usize,
    sequence_aa: &str,
) -> String {
    if region_end <= region_start {
        return String::new();
    }
    // First codon-aligned position at or after `region_start`. Mirror
    // the Python `(frame_offset - region_start) % 3` semantics with
    // `rem_euclid` so a negative arithmetic difference still produces
    // a non-negative remainder in `[0, 3)`.
    let delta = (frame_offset as i64 - region_start as i64).rem_euclid(3) as usize;
    let codon_start = region_start + delta;
    if codon_start >= region_end {
        return String::new();
    }
    let n_codons = (region_end - codon_start) / 3;
    if n_codons == 0 {
        return String::new();
    }
    let aa_index = (codon_start - frame_offset) / 3;
    let aa_bytes = sequence_aa.as_bytes();
    if aa_index >= aa_bytes.len() {
        return String::new();
    }
    let end = (aa_index + n_codons).min(aa_bytes.len());
    String::from_utf8_lossy(&aa_bytes[aa_index..end]).into_owned()
}

pub(super) fn junction_has_stop(seq: &str) -> bool {
    let bytes = seq.as_bytes();
    let n = bytes.len() / 3;
    for i in 0..n {
        let start = i * 3;
        let mut codon = [0u8; 3];
        for j in 0..3 {
            codon[j] = match bytes[start + j] {
                b'U' | b'u' => b'T',
                c => c.to_ascii_uppercase(),
            };
        }
        if (codon == *b"TAA") || (codon == *b"TAG") || (codon == *b"TGA") {
            return true;
        }
    }
    false
}
