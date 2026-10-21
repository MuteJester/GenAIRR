//! Portable effective native reference; ordered pools preserve allele IDs.
use crate::ir::Segment;
use crate::refdata::*;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SegmentSpec {
    V,
    D,
    J,
    Np1,
    Np2,
}
impl SegmentSpec {
    pub fn native(self) -> Segment {
        match self {
            Self::V => Segment::V,
            Self::D => Segment::D,
            Self::J => Segment::J,
            Self::Np1 => Segment::Np1,
            Self::Np2 => Segment::Np2,
        }
    }
}
impl From<Segment> for SegmentSpec {
    fn from(s: Segment) -> Self {
        match s {
            Segment::V => Self::V,
            Segment::D => Self::D,
            Segment::J => Self::J,
            Segment::Np1 => Self::Np1,
            Segment::Np2 => Self::Np2,
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IdentitySpec {
    pub species: Option<String>,
    pub locus: Option<String>,
    pub reference_set: Option<String>,
    pub name: Option<String>,
    pub source: Option<String>,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SeveritySpec {
    Fatal,
    Curatable,
}
impl SeveritySpec {
    fn native(self) -> RefDataIssueSeverity {
        match self {
            Self::Fatal => RefDataIssueSeverity::Fatal,
            Self::Curatable => RefDataIssueSeverity::Curatable,
        }
    }
    fn capture(s: RefDataIssueSeverity) -> Self {
        match s {
            RefDataIssueSeverity::Fatal => Self::Fatal,
            RefDataIssueSeverity::Curatable => Self::Curatable,
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AnchorSpec {
    pub required: bool,
    pub expected_amino_acids: Vec<char>,
    pub missing_severity: SeveritySpec,
    pub mismatch_severity: SeveritySpec,
}
impl AnchorSpec {
    fn capture(r: &AnchorRule) -> Self {
        Self {
            required: r.required,
            expected_amino_acids: r.expected_amino_acids.clone(),
            missing_severity: SeveritySpec::capture(r.missing_severity),
            mismatch_severity: SeveritySpec::capture(r.mismatch_severity),
        }
    }
    fn native(&self) -> AnchorRule {
        AnchorRule {
            required: self.required,
            expected_amino_acids: self.expected_amino_acids.clone(),
            missing_severity: self.missing_severity.native(),
            mismatch_severity: self.mismatch_severity.native(),
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RulesSpec {
    pub alphabet: Vec<u8>,
    pub v_anchor: AnchorSpec,
    pub j_anchor: AnchorSpec,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FunctionSpec {
    Functional,
    Orf,
    Pseudogene,
    Unknown,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SubregionSpec {
    pub label: String,
    pub start: u16,
    pub end: u16,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AlleleSpec {
    pub name: String,
    pub gene: String,
    pub sequence: String,
    pub segment: SegmentSpec,
    pub anchor: Option<u16>,
    pub functionality: Option<FunctionSpec>,
    pub subregions: Vec<SubregionSpec>,
    pub evidence: Option<AlleleEvidence>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferenceSpecV1 {
    pub schema_version: u32,
    pub chain_type: String,
    pub identity: IdentitySpec,
    pub rules: RulesSpec,
    pub v: Vec<AlleleSpec>,
    pub d: Vec<AlleleSpec>,
    pub j: Vec<AlleleSpec>,
    pub c: Vec<AlleleSpec>,
}
impl ReferenceSpecV1 {
    pub fn capture(r: &RefDataConfig) -> Result<Self, String> {
        let pools = [&r.v_pool, &r.d_pool, &r.j_pool, &r.c_pool];
        let mut keys = std::collections::HashSet::new();
        for pool in pools {
            for (_, allele) in pool.iter() {
                crate::family_replay::primitive_work(1);
                crate::family_replay::primitive_bases(allele.seq.len() as u64);
                keys.insert((allele.segment, allele.name.as_str()));
            }
        }
        for key in r.evidence.keys() {
            crate::family_replay::primitive_work(1);
            if !keys.contains(&(key.0, key.1.as_str())) {
                return Err("orphan reference evidence cannot be encoded".into());
            }
        }
        let pool = |p: &AllelePool| -> Result<Vec<AlleleSpec>, String> {
            p.iter()
                .map(|(_, a)| {
                    Ok(AlleleSpec {
                        name: a.name.clone(),
                        gene: a.gene.clone(),
                        sequence: String::from_utf8(a.seq.clone())
                            .map_err(|_| "reference sequence is not UTF-8")?,
                        segment: a.segment.into(),
                        anchor: a.anchor,
                        functionality: a.functional_status.map(|s| match s {
                            FunctionalStatus::Functional => FunctionSpec::Functional,
                            FunctionalStatus::Orf => FunctionSpec::Orf,
                            FunctionalStatus::Pseudogene => FunctionSpec::Pseudogene,
                            FunctionalStatus::Unknown => FunctionSpec::Unknown,
                        }),
                        subregions: a
                            .subregions
                            .iter()
                            .map(|s| SubregionSpec {
                                label: s.label.as_str().into(),
                                start: s.start,
                                end: s.end,
                            })
                            .collect(),
                        evidence: r.evidence.get(&(a.segment, a.name.clone())).cloned(),
                    })
                })
                .collect()
        };
        let id = &r.identity;
        Ok(Self {
            schema_version: 1,
            chain_type: match r.chain_type {
                ChainType::Vj => "vj",
                ChainType::Vdj => "vdj",
            }
            .into(),
            identity: IdentitySpec {
                species: id.species.clone(),
                locus: id.locus.clone(),
                reference_set: id.reference_set.clone(),
                name: id.name.clone(),
                source: id.source.clone(),
            },
            rules: RulesSpec {
                alphabet: r.rules.alphabet.allowed.clone(),
                v_anchor: AnchorSpec::capture(&r.rules.v_anchor),
                j_anchor: AnchorSpec::capture(&r.rules.j_anchor),
            },
            v: pool(&r.v_pool)?,
            d: pool(&r.d_pool)?,
            j: pool(&r.j_pool)?,
            c: pool(&r.c_pool)?,
        })
    }
    pub fn reconstruct(&self) -> Result<RefDataConfig, String> {
        if self.schema_version != 1 {
            return Err("unsupported reference schema".into());
        }
        let ct = match self.chain_type.as_str() {
            "vj" => ChainType::Vj,
            "vdj" => ChainType::Vdj,
            _ => return Err("unknown reference chain type".into()),
        };
        let mut r = RefDataConfig::empty(ct);
        let id = &self.identity;
        r.identity = ReferenceIdentity {
            species: id.species.clone(),
            locus: id.locus.clone(),
            reference_set: id.reference_set.clone(),
            name: id.name.clone(),
            source: id.source.clone(),
        };
        r.rules = ReferenceRules {
            alphabet: ReferenceAlphabet {
                allowed: self.rules.alphabet.clone(),
            },
            v_anchor: self.rules.v_anchor.native(),
            j_anchor: self.rules.j_anchor.native(),
        };
        for (label, values) in [
            ("v", &self.v),
            ("d", &self.d),
            ("j", &self.j),
            ("c", &self.c),
        ] {
            let mut alleles = Vec::with_capacity(values.len());
            for a in values {
                if !a.sequence.is_ascii() {
                    return Err("reference bases must be ASCII".into());
                }
                let segment = a.segment.native();
                let expected = match label {
                    "v" => Some(Segment::V),
                    "d" => Some(Segment::D),
                    "j" => Some(Segment::J),
                    _ => None,
                };
                if expected.is_some_and(|s| s != segment) {
                    return Err("reference pool/segment mismatch".into());
                }
                let subregions = a
                    .subregions
                    .iter()
                    .map(|s| {
                        Ok(VSubregion {
                            label: VSubregionLabel::from_str(&s.label)
                                .ok_or("unknown V subregion")?,
                            start: s.start,
                            end: s.end,
                        })
                    })
                    .collect::<Result<Vec<_>, String>>()?;
                let allele = Allele {
                    name: a.name.clone(),
                    gene: a.gene.clone(),
                    seq: a.sequence.as_bytes().to_vec(),
                    segment,
                    anchor: a.anchor,
                    functional_status: a.functionality.map(|s| match s {
                        FunctionSpec::Functional => FunctionalStatus::Functional,
                        FunctionSpec::Orf => FunctionalStatus::Orf,
                        FunctionSpec::Pseudogene => FunctionalStatus::Pseudogene,
                        FunctionSpec::Unknown => FunctionalStatus::Unknown,
                    }),
                    subregions,
                };
                if let Some(e) = &a.evidence {
                    if !e.problems(&allele).is_empty() {
                        return Err("invalid sequence-bound reference evidence".into());
                    }
                    if r.evidence
                        .insert((segment, a.name.clone()), e.clone())
                        .is_some()
                    {
                        return Err("duplicate reference evidence key".into());
                    }
                }
                alleles.push(allele);
            }
            let p = AllelePool::from_vec(alleles);
            match label {
                "v" => r.v_pool = p,
                "d" => r.d_pool = p,
                "j" => r.j_pool = p,
                _ => r.c_pool = p,
            }
        }
        let fatal: Vec<_> = r
            .validate()
            .into_iter()
            .filter(|i| i.severity() == RefDataIssueSeverity::Fatal)
            .collect();
        if !fatal.is_empty() {
            return Err(format!("invalid reference: {fatal:?}"));
        }
        Ok(r)
    }
    pub fn digest(&self) -> Result<String, String> {
        crate::family_replay::canonical_digest(
            "reference_spec_v1",
            &serde_json::to_value(self).map_err(|e| e.to_string())?,
        )
    }
}
