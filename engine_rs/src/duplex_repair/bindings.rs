//! Verified whole-allele bindings for standalone duplex projection.
//!
//! Source identity is chosen explicitly from reference data, never inferred
//! from an identical sequence. The checksum preserves the exact stored bytes;
//! substrate comparison permits only the kernel's ASCII case normalization.
use super::{CoupledProduct, Substrate};
use crate::ir::Segment;
use crate::refdata::evidence::sequence_sha256;
use crate::refdata::{
    AlleleEvidence, AlleleId, CodingBoundaryStatus, RefDataConfig, RefDataValidationMode,
};
use serde::Serialize;
use serde_json::{json, Value};

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct DuplexAlleleBinding {
    pub source_instance: String,
    pub segment: Substrate,
    pub allele_id: u32,
    pub allele_name: String,
    /// Exact canonical bytes as stored in the reference catalogue.
    pub sequence: String,
    pub sequence_sha256: String,
    pub canonical_start: usize,
    pub canonical_end: usize,
    pub anchor: Option<usize>,
    pub landmark_ready: bool,
    pub coding_phase: Option<usize>,
    pub evidence: Option<AlleleEvidence>,
    pub boundary_5: CodingBoundaryStatus,
    pub boundary_3: CodingBoundaryStatus,
    /// Complete evidence is an author declaration, not independent certification.
    pub boundary_contract: String,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct DuplexReferenceBindings {
    pub schema_version: u32,
    pub reference_content_hash: String,
    pub reference_identity: Value,
    pub substrate_normalization: String,
    pub v: DuplexAlleleBinding,
    pub d: DuplexAlleleBinding,
    pub j: DuplexAlleleBinding,
    pub d_inverted: bool,
    pub reference_ends_are_coding_ends: bool,
}

impl DuplexReferenceBindings {
    /// Bind the exact V/D/J names and complete canonical allele substrates.
    /// Unknown junction-facing boundaries require an explicit caller assertion.
    /// A known partial boundary always contradicts that assertion.
    pub fn new(
        refdata: &RefDataConfig,
        calls: [&str; 3],
        substrates: [&str; 3],
        d_inverted: bool,
        reference_ends_are_coding_ends: bool,
    ) -> Result<Self, String> {
        if !refdata.chain_type.has_d() {
            return Err("duplex projection requires VDJ reference data".into());
        }
        // Curatable landmark uncertainty belongs in a partial projection.
        // Structural corruption and invalid sequence-bound evidence do not.
        refdata
            .validate_with_mode(RefDataValidationMode::AllowCuratable)
            .map_err(|e| format!("invalid duplex reference data: {e}"))?;
        let v = bind_allele(
            refdata,
            Substrate::V,
            calls[0],
            substrates[0],
            reference_ends_are_coding_ends,
        )?;
        let d = bind_allele(
            refdata,
            Substrate::D,
            calls[1],
            substrates[1],
            reference_ends_are_coding_ends,
        )?;
        let j = bind_allele(
            refdata,
            Substrate::J,
            calls[2],
            substrates[2],
            reference_ends_are_coding_ends,
        )?;
        let identity = &refdata.identity;
        Ok(Self {
            schema_version: 1,
            reference_content_hash: crate::trace_file::refdata_content_hash(refdata),
            reference_identity: json!({
                "species": identity.species,
                "locus": identity.locus,
                "reference_set": identity.reference_set,
                "name": identity.name,
                "source": identity.source,
            }),
            substrate_normalization: "ascii_uppercase_v1".into(),
            v,
            d,
            j,
            d_inverted,
            reference_ends_are_coding_ends,
        })
    }

    pub fn for_substrate(&self, substrate: Substrate) -> &DuplexAlleleBinding {
        match substrate {
            Substrate::V => &self.v,
            Substrate::D => &self.d,
            Substrate::J => &self.j,
        }
    }

    pub fn for_segment(&self, segment: Segment) -> Option<&DuplexAlleleBinding> {
        match segment {
            Segment::V => Some(&self.v),
            Segment::D => Some(&self.d),
            Segment::J => Some(&self.j),
            _ => None,
        }
    }

    pub fn allele_id(&self, substrate: Substrate) -> AlleleId {
        AlleleId::new(self.for_substrate(substrate).allele_id)
    }

    /// Rebuild all evidence from the supplied reference and reject edited
    /// identities, readiness, coverage, orientation or substrate information.
    pub fn validate(
        &self,
        product: &CoupledProduct,
        refdata: &RefDataConfig,
    ) -> Result<(), String> {
        if product.d_inverted != self.d_inverted {
            return Err("duplex binding D orientation does not match the repair product".into());
        }
        let expected = Self::new(
            refdata,
            [
                &self.v.allele_name,
                &self.d.allele_name,
                &self.j.allele_name,
            ],
            [
                &product.v_reference,
                &product.d_reference,
                &product.j_reference,
            ],
            product.d_inverted,
            self.reference_ends_are_coding_ends,
        )?;
        if self != &expected {
            return Err("duplex binding identity or evidence does not match reference data".into());
        }
        Ok(())
    }
}

fn bind_allele(
    refdata: &RefDataConfig,
    substrate: Substrate,
    call: &str,
    supplied: &str,
    asserted: bool,
) -> Result<DuplexAlleleBinding, String> {
    let (segment, label) = match substrate {
        Substrate::V => (Segment::V, "v"),
        Substrate::D => (Segment::D, "d"),
        Substrate::J => (Segment::J, "j"),
    };
    let pool = refdata.pool_for(segment).ok_or("missing reference pool")?;
    let mut matches = pool.iter().filter(|(_, allele)| allele.name == call);
    let (id, allele) = matches
        .next()
        .ok_or_else(|| format!("unknown generating {label} allele {call:?}"))?;
    if matches.next().is_some() {
        return Err(format!("ambiguous generating {label} allele name {call:?}"));
    }
    if allele.segment != segment {
        return Err(format!(
            "bound {label} allele has an inconsistent segment role"
        ));
    }
    if allele.seq.is_empty()
        || !allele
            .seq
            .iter()
            .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
    {
        return Err(format!(
            "bound {label} allele must contain nonempty resolved A/C/G/T"
        ));
    }
    if !supplied.as_bytes().eq_ignore_ascii_case(&allele.seq) {
        return Err(format!(
            "{label} substrate must match the entire canonical generating allele {call:?}"
        ));
    }
    let evidence = refdata.allele_evidence(segment, id).cloned();
    let (boundary_5, boundary_3) = evidence
        .as_ref()
        .map(|e| (e.boundary_5, e.boundary_3))
        .unwrap_or((CodingBoundaryStatus::Unknown, CodingBoundaryStatus::Unknown));
    let ends: &[CodingBoundaryStatus] = match substrate {
        Substrate::V => &[boundary_3],
        Substrate::D => &[boundary_5, boundary_3],
        Substrate::J => &[boundary_5],
    };
    if ends.contains(&CodingBoundaryStatus::Partial) {
        return Err(format!(
            "bound {label} allele has a known partial coding end; a caller assertion cannot override it"
        ));
    }
    let complete = ends.iter().all(|b| *b == CodingBoundaryStatus::Complete);
    if !complete && !asserted {
        return Err(format!(
            "bound {label} allele has an unknown coding end; explicitly set reference_ends_are_coding_ends=True"
        ));
    }
    Ok(DuplexAlleleBinding {
        source_instance: format!("{label}:0"),
        segment: substrate,
        allele_id: id.index(),
        allele_name: allele.name.clone(),
        sequence: String::from_utf8(allele.seq.clone()).expect("resolved ASCII was checked"),
        sequence_sha256: sequence_sha256(&allele.seq),
        canonical_start: 0,
        canonical_end: allele.seq.len(),
        anchor: allele.anchor.map(usize::from),
        landmark_ready: crate::productivity::reference_ready(refdata, segment, id),
        coding_phase: crate::productivity::reference_phase_for(refdata, segment, id),
        evidence,
        boundary_5,
        boundary_3,
        boundary_contract: if complete {
            "declared_complete"
        } else {
            "caller_asserted"
        }
        .into(),
    })
}

#[cfg(test)]
mod tests;
