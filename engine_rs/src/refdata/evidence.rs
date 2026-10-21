//! Versioned, sequence-bound evidence; absence denotes a legacy annotation.
use super::{Allele, AlleleId, RefDataConfig};
use crate::ir::Segment;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AnchorEvidenceStatus {
    Verified,
    Inferred,
    Unknown,
    Invalid,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CodingBoundaryStatus {
    Complete,
    Partial,
    Unknown,
}

/// Claims are supplied by the author, not certified merely by passing validation.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AlleleEvidence {
    pub schema_version: u32,
    pub sequence_sha256: String,
    pub coding_phase: Option<u8>,
    pub anchor_status: AnchorEvidenceStatus,
    pub boundary_5: CodingBoundaryStatus,
    pub boundary_3: CodingBoundaryStatus,
    pub source: Option<String>,
    pub method: String,
    #[serde(default)]
    pub notes: Vec<String>,
}

pub fn sequence_sha256(seq: &[u8]) -> String {
    format!("sha256:{:x}", Sha256::digest(seq))
}

impl AlleleEvidence {
    pub fn problems(&self, allele: &Allele) -> Vec<String> {
        let mut out = Vec::new();
        if self.schema_version != 1 {
            out.push("unsupported evidence schema_version".into());
        }
        if self.sequence_sha256 != sequence_sha256(&allele.seq) {
            out.push("evidence sequence checksum does not match allele".into());
        }
        if self.method.trim().is_empty() {
            out.push("evidence method must be nonempty".into());
        }
        if self.anchor_status == AnchorEvidenceStatus::Verified
            && self.source.as_deref().is_none_or(|s| s.trim().is_empty())
        {
            out.push("verified landmark evidence requires a source identifier".into());
        }
        if self.coding_phase.is_some_and(|p| p > 2) {
            out.push("coding_phase must be 0, 1 or 2".into());
        }
        if matches!(
            self.anchor_status,
            AnchorEvidenceStatus::Verified | AnchorEvidenceStatus::Inferred
        ) {
            if matches!(allele.segment, Segment::V | Segment::J) && allele.anchor.is_none() {
                out.push("usable anchor evidence requires an anchor coordinate".into());
            }
            if let (Some(phase), Some(anchor)) = (self.coding_phase, allele.anchor) {
                if u16::from(phase) != anchor % 3 {
                    out.push("coding_phase disagrees with anchor coordinate".into());
                }
            }
        }
        out
    }
}

impl RefDataConfig {
    pub fn allele_evidence(&self, segment: Segment, id: AlleleId) -> Option<&AlleleEvidence> {
        if self.evidence.is_empty() {
            return None;
        }
        let allele = self.get(segment, id)?;
        self.evidence.get(&(segment, allele.name.clone()))
    }
    /// Unknown boundaries may be covered by the explicit caller assertion.
    /// A known partial junction-facing end contradicts original-end processing.
    pub fn coding_ends_are_not_known_partial(&self, segment: Segment, id: AlleleId) -> bool {
        let Some(e) = self.allele_evidence(segment, id) else {
            return true;
        };
        match segment {
            Segment::V => e.boundary_3 != CodingBoundaryStatus::Partial,
            Segment::J => e.boundary_5 != CodingBoundaryStatus::Partial,
            Segment::D => {
                e.boundary_5 != CodingBoundaryStatus::Partial
                    && e.boundary_3 != CodingBoundaryStatus::Partial
            }
            _ => false,
        }
    }

    pub fn set_allele_evidence(
        &mut self,
        segment: Segment,
        id: AlleleId,
        evidence: AlleleEvidence,
    ) -> Result<(), String> {
        let allele = self.get(segment, id).ok_or("unknown allele id")?;
        let problems = evidence.problems(allele);
        if !problems.is_empty() {
            return Err(problems.join("; "));
        }
        self.evidence
            .insert((segment, allele.name.clone()), evidence);
        Ok(())
    }
}

/// Conservative authoring resolver. It never equates an arbitrary terminal
/// cysteine with the IMGT conserved V landmark. A unique J F/W-G-X-G motif
/// is an explicitly *inferred* candidate, not independently verified evidence.
#[derive(Clone, Debug, Serialize)]
pub struct AnchorResolution {
    pub position: Option<u16>,
    pub status: AnchorEvidenceStatus,
    pub method: String,
    pub reason: String,
    pub candidates: Vec<u16>,
}

pub fn resolve_anchor(
    seq: &[u8],
    segment: Segment,
    locus: &str,
    gapped: Option<&[u8]>,
    imgt_gapped: bool,
    explicit: Option<u16>,
) -> AnchorResolution {
    let expected = |aa: char| match segment {
        Segment::V => aa == 'C',
        Segment::J => {
            if locus.eq_ignore_ascii_case("IGH") {
                aa == 'W'
            } else if ["IGK", "IGL", "TRA", "TRB", "TRG", "TRD"]
                .iter()
                .any(|x| x.eq_ignore_ascii_case(locus))
            {
                aa == 'F'
            } else {
                aa == 'F' || aa == 'W'
            }
        }
        _ => false,
    };
    let valid = |pos: usize| {
        seq.get(pos..pos + 3)
            .map(|c| expected(crate::codon::translate_codon_slice(c)))
            .unwrap_or(false)
    };
    let result = |position, status, method: &str, reason: &str, candidates| AnchorResolution {
        position,
        status,
        method: method.into(),
        reason: reason.into(),
        candidates,
    };
    if let Some(pos) = explicit {
        return if valid(pos as usize) {
            result(
                Some(pos),
                AnchorEvidenceStatus::Inferred,
                "explicit_coordinate",
                "coordinate supplied by author; codon validated",
                vec![pos],
            )
        } else {
            result(
                None,
                AnchorEvidenceStatus::Invalid,
                "explicit_coordinate",
                "supplied coordinate does not encode the expected landmark",
                vec![pos],
            )
        };
    }
    if segment == Segment::V && imgt_gapped {
        if let Some(g) = gapped {
            let ungapped: Vec<_> = g
                .iter()
                .copied()
                .filter(|b| *b != b'.' && *b != b'-')
                .collect();
            if !ungapped.eq_ignore_ascii_case(seq) {
                return result(
                    None,
                    AnchorEvidenceStatus::Invalid,
                    "imgt_gapped_104",
                    "gapped and ungapped sequences disagree",
                    vec![],
                );
            }
            // IMGT conserved Cys104 starts at nucleotide 309 (zero-based).
            if g.len() >= 312 && !g[309..312].iter().any(|b| *b == b'.' || *b == b'-') {
                let pos = g[..309]
                    .iter()
                    .filter(|b| **b != b'.' && **b != b'-')
                    .count();
                if pos <= u16::MAX as usize && valid(pos) {
                    return result(
                        Some(pos as u16),
                        AnchorEvidenceStatus::Inferred,
                        "imgt_gapped_104",
                        "declared IMGT alignment; conserved codon validated",
                        vec![pos as u16],
                    );
                }
            }
            return result(
                None,
                AnchorEvidenceStatus::Invalid,
                "imgt_gapped_104",
                "declared IMGT landmark is absent or does not encode Cys",
                vec![],
            );
        }
    }
    if segment == Segment::J {
        let mut candidates = Vec::new();
        for p in 0..seq.len().saturating_sub(11) {
            if !valid(p) {
                continue;
            }
            let aa = |q| crate::codon::translate_codon_slice(&seq[q..q + 3]);
            if aa(p + 3) == 'G' && aa(p + 9) == 'G' && p <= u16::MAX as usize {
                candidates.push(p as u16);
            }
        }
        return if candidates.len() == 1 {
            result(
                Some(candidates[0]),
                AnchorEvidenceStatus::Inferred,
                "unique_j_fwgxg",
                "unique locus-compatible F/W-G-X-G candidate; not source-verified",
                candidates,
            )
        } else {
            result(
                None,
                AnchorEvidenceStatus::Unknown,
                "unique_j_fwgxg",
                "zero or multiple candidate landmarks; source annotation required",
                candidates,
            )
        };
    }
    result(
        None,
        AnchorEvidenceStatus::Unknown,
        "unresolved",
        "no independently positioned landmark supplied",
        vec![],
    )
}
