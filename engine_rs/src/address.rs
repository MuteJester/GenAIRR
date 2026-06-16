//! Canonical trace/pass address strings.
//!
//! Trace addresses are part of the engine's internal contract: passes
//! write them, contracts and feasibility filters inspect them, AIRR
//! projection reads them, and tests assert them. Keep the spelling here
//! as the single source of truth.
//!
//! # Address schema versioning
//!
//! The strings produced by [`ChoiceAddress::Display`] and consumed by
//! [`ChoiceAddress::parse`] form the **persisted vocabulary** of every
//! trace file written by the engine. Once a trace file has been emitted,
//! that vocabulary must remain parseable forever or replay of the file
//! will fail at the cursor's address-match step.
//!
//! [`ADDRESS_SCHEMA_VERSION`] tracks that vocabulary as a single
//! integer carried inside every [`crate::trace_file::TraceFile`].
//!
//! ## What triggers a bump
//!
//! Any change that breaks `Display ↔ parse` for an existing typed
//! variant — renaming a constant (`"np.np1.length"` → `"np.1.length"`),
//! changing an indexed prefix, dropping a variant, changing how a
//! variant is parameterised. The compile-fence
//! `frozen_address_spellings` test below pins one representative of
//! every variant so a code change that drifts the on-disk vocabulary
//! has to come with an explicit version bump.
//!
//! ## What does *not* trigger a bump
//!
//! Adding **new** variants whose Display strings don't collide with
//! existing prefixes. Old traces won't reference the new addresses;
//! new traces are parseable by future engines.

use crate::assignment::TrimEnd;
use crate::ir::Segment;
use std::{fmt, str::FromStr};

/// Current revision of the persisted trace-address vocabulary. Bumped
/// when any existing [`ChoiceAddress`] variant changes its
/// `Display` spelling or `parse` shape. See module docs for the
/// bump policy.
pub const ADDRESS_SCHEMA_VERSION: u32 = 1;

const SAMPLE_ALLELE_V: &str = "sample_allele.v";
const SAMPLE_ALLELE_D: &str = "sample_allele.d";
const SAMPLE_ALLELE_J: &str = "sample_allele.j";
/// Bool — `true` when the sampled D segment is to be assembled in
/// reverse-complement (V(D)J inversion event). Recorded by the
/// `InvertDPass` (see [`crate::passes::InvertDPass`]) so trace
/// replay can re-fire the same orientation decision deterministically.
/// Sub-namespace of `sample_allele.d` so an existing
/// `sample_allele.*`-aware consumer keeps working.
const SAMPLE_ALLELE_D_INVERTED: &str = "sample_allele.d.inverted";
pub const SAMPLE_ALLELE_INVALID: &str = "sample_allele.<invalid>";
pub const SAMPLE_ALLELE_UNSUPPORTED: &str = "sample_allele.<unsupported>";
/// Per-rearrangement chromosome choice for a phased genotype.
const SAMPLE_HAPLOTYPE: &str = "sample_haplotype";

/// Pass name for `InvertDPass`. Used as the `name()` return value
/// and as the pass-plan signature token, so external consumers
/// (trace replay, MCP audits, schedule reports) keep a stable key
/// once Slice C lands the pass.
pub const INVERT_D: &str = "invert_d";

const TRIM_V_5: &str = "trim.v_5";
const TRIM_V_3: &str = "trim.v_3";
const TRIM_D_5: &str = "trim.d_5";
const TRIM_D_3: &str = "trim.d_3";
const TRIM_J_5: &str = "trim.j_5";
const TRIM_J_3: &str = "trim.j_3";
pub const TRIM_INVALID: &str = "trim.<invalid>";
pub const TRIM_UNSUPPORTED: &str = "trim.<unsupported>";

const ASSEMBLE_V: &str = "assemble.v";
const ASSEMBLE_D: &str = "assemble.d";
const ASSEMBLE_J: &str = "assemble.j";

const GENERATE_NP1: &str = "generate_np.np1";
const GENERATE_NP2: &str = "generate_np.np2";

const NP1_LENGTH: &str = "np.np1.length";
const NP2_LENGTH: &str = "np.np2.length";
pub const NP_INVALID_LENGTH: &str = "np.<invalid>.length";
const NP1_BASES_INDEX_PREFIX: &str = "np.np1.bases[";
const NP2_BASES_INDEX_PREFIX: &str = "np.np2.bases[";
const NP1_BASES_PATTERN: &str = "np.np1.bases[0..n]";
const NP2_BASES_PATTERN: &str = "np.np2.bases[0..n]";

// P-nucleotide (palindromic addition) length addresses. Per-end:
// each `PAdditionPass` records exactly one `Int(length)` choice
// at the canonical spelling below. Bases derive deterministically
// from `(assigned allele, trim, orientation, length)` via
// `complement_base` — no per-base trace records exist.
//
// Pass-name spellings follow the per-end suffix convention used
// by `trim.{v_3,d_5,d_3,j_5}` so a pass-name-aware consumer
// keeps working.
const P_V3_LENGTH: &str = "p.v_3.length";
const P_D5_LENGTH: &str = "p.d_5.length";
const P_D3_LENGTH: &str = "p.d_3.length";
const P_J5_LENGTH: &str = "p.j_5.length";
pub const P_ADDITION_V_3: &str = "p_addition.v_3";
pub const P_ADDITION_D_5: &str = "p_addition.d_5";
pub const P_ADDITION_D_3: &str = "p_addition.d_3";
pub const P_ADDITION_J_5: &str = "p_addition.j_5";

pub const MUTATE_UNIFORM: &str = "mutate.uniform";
const MUTATE_UNIFORM_COUNT: &str = "mutate.uniform.count";
const MUTATE_UNIFORM_SITE_PATTERN: &str = "mutate.uniform.site[0..n]";
const MUTATE_UNIFORM_BASE_PATTERN: &str = "mutate.uniform.base[0..n]";
const MUTATE_UNIFORM_SITE_PREFIX: &str = "mutate.uniform.site[";
const MUTATE_UNIFORM_BASE_PREFIX: &str = "mutate.uniform.base[";

pub const MUTATE_S5F: &str = "mutate.s5f";
const MUTATE_S5F_COUNT: &str = "mutate.s5f.count";
const MUTATE_S5F_SITE_PATTERN: &str = "mutate.s5f.site[0..n]";
const MUTATE_S5F_BASE_PATTERN: &str = "mutate.s5f.base[0..n]";
const MUTATE_S5F_SITE_PREFIX: &str = "mutate.s5f.site[";
const MUTATE_S5F_BASE_PREFIX: &str = "mutate.s5f.base[";

pub const CORRUPT_PCR: &str = "corrupt.pcr";
const CORRUPT_PCR_COUNT: &str = "corrupt.pcr.count";
const CORRUPT_PCR_SITE_PATTERN: &str = "corrupt.pcr.error_site[0..n]";
const CORRUPT_PCR_BASE_PATTERN: &str = "corrupt.pcr.error_base[0..n]";
const CORRUPT_PCR_SITE_PREFIX: &str = "corrupt.pcr.error_site[";
const CORRUPT_PCR_BASE_PREFIX: &str = "corrupt.pcr.error_base[";

pub const CORRUPT_QUALITY: &str = "corrupt.quality";
const CORRUPT_QUALITY_COUNT: &str = "corrupt.quality.count";
const CORRUPT_QUALITY_SITE_PATTERN: &str = "corrupt.quality.error_site[0..n]";
const CORRUPT_QUALITY_BASE_PATTERN: &str = "corrupt.quality.error_base[0..n]";
const CORRUPT_QUALITY_SITE_PREFIX: &str = "corrupt.quality.error_site[";
const CORRUPT_QUALITY_BASE_PREFIX: &str = "corrupt.quality.error_base[";

pub const CORRUPT_CONTAMINANT: &str = "corrupt.contaminant";
const CORRUPT_CONTAMINANT_APPLIED: &str = "corrupt.contaminant.applied";
const CORRUPT_CONTAMINANT_BASES_PATTERN: &str = "corrupt.contaminant.bases[0..n]";
const CORRUPT_CONTAMINANT_BASES_PREFIX: &str = "corrupt.contaminant.bases[";

pub const CORRUPT_INDEL: &str = "corrupt.indel";
const CORRUPT_INDEL_COUNT: &str = "corrupt.indel.count";
const CORRUPT_INDEL_KIND_PATTERN: &str = "corrupt.indel.kind[0..n]";
const CORRUPT_INDEL_SITE_PATTERN: &str = "corrupt.indel.site[0..n]";
const CORRUPT_INDEL_BASE_PATTERN: &str = "corrupt.indel.base[0..n]";
const CORRUPT_INDEL_KIND_PREFIX: &str = "corrupt.indel.kind[";
const CORRUPT_INDEL_SITE_PREFIX: &str = "corrupt.indel.site[";
const CORRUPT_INDEL_BASE_PREFIX: &str = "corrupt.indel.base[";

pub const CORRUPT_NS: &str = "corrupt.ns";
const CORRUPT_NS_COUNT: &str = "corrupt.ns.count";
const CORRUPT_NS_SITE_PATTERN: &str = "corrupt.ns.site[0..n]";
const CORRUPT_NS_SITE_PREFIX: &str = "corrupt.ns.site[";

pub const CORRUPT_END_LOSS_5: &str = "corrupt.end_loss.5";
pub const CORRUPT_END_LOSS_3: &str = "corrupt.end_loss.3";

pub const CORRUPT_REV_COMP: &str = "corrupt.rev_comp";
const CORRUPT_REV_COMP_APPLIED: &str = "corrupt.rev_comp.applied";

/// Pass name for `ReceptorRevisionPass`. Stable identifier referenced
/// by `pass_plan_signature` and trace consumers.
pub const RECEPTOR_REVISION: &str = "receptor_revision";

/// Bool — `true` when the receptor-revision pass replaces the V
/// segment for this simulation. One record per simulation, regardless
/// of outcome.
const RECEPTOR_REVISION_APPLIED: &str = "receptor_revision.applied";
/// AlleleId — the replacement V allele chosen on `applied=true`.
/// Absent on `applied=false`.
const RECEPTOR_REVISION_V_ALLELE: &str = "receptor_revision.v_allele";
/// Int — the 3' trim sampled for the replacement V allele on
/// `applied=true`. Absent on `applied=false`. 5' trim stays at 0
/// for receptor revision v1 (single-length-preserving constraint).
const RECEPTOR_REVISION_V_TRIM_3: &str = "receptor_revision.v_trim_3";

/// Pass name for `PairedEndSamplingPass`. Stable identifier
/// referenced by `pass_plan_signature` and trace consumers.
pub const PAIRED_END: &str = "paired_end";

/// Int — the R1 read length sampled for this simulation. Recorded
/// once per simulation by `PairedEndSamplingPass`; absent when no
/// paired-end layout is in the plan.
const PAIRED_END_R1_LENGTH: &str = "paired_end.r1_length";
/// Int — the R2 read length sampled for this simulation.
const PAIRED_END_R2_LENGTH: &str = "paired_end.r2_length";
/// Int — the fragment insert size sampled for this simulation.
const PAIRED_END_INSERT_SIZE: &str = "paired_end.insert_size";

pub fn sample_allele_vdj(segment: Segment) -> &'static str {
    match segment {
        Segment::V => SAMPLE_ALLELE_V,
        Segment::D => SAMPLE_ALLELE_D,
        Segment::J => SAMPLE_ALLELE_J,
        Segment::Np1 | Segment::Np2 => {
            panic!("sample allele address is only defined for V/D/J")
        }
    }
}

pub fn trim_vdj(segment: Segment, end: TrimEnd) -> &'static str {
    match segment {
        Segment::V => match end {
            TrimEnd::Five => TRIM_V_5,
            TrimEnd::Three => TRIM_V_3,
        },
        Segment::D => match end {
            TrimEnd::Five => TRIM_D_5,
            TrimEnd::Three => TRIM_D_3,
        },
        Segment::J => match end {
            TrimEnd::Five => TRIM_J_5,
            TrimEnd::Three => TRIM_J_3,
        },
        Segment::Np1 | Segment::Np2 => panic!("trim address is only defined for V/D/J"),
    }
}

pub fn assemble_vdj(segment: Segment) -> &'static str {
    match segment {
        Segment::V => ASSEMBLE_V,
        Segment::D => ASSEMBLE_D,
        Segment::J => ASSEMBLE_J,
        Segment::Np1 | Segment::Np2 => panic!("assemble address is only defined for V/D/J"),
    }
}

pub fn generate_np_region(segment: Segment) -> &'static str {
    match segment {
        Segment::Np1 => GENERATE_NP1,
        Segment::Np2 => GENERATE_NP2,
        Segment::V | Segment::D | Segment::J => {
            panic!("generate_np address is only defined for NP1/NP2")
        }
    }
}

pub fn np_length_region(segment: Segment) -> &'static str {
    match segment {
        Segment::Np1 => NP1_LENGTH,
        Segment::Np2 => NP2_LENGTH,
        Segment::V | Segment::D | Segment::J => {
            panic!("NP length address is only defined for NP1/NP2")
        }
    }
}

fn parse_indexed(address: &str, prefix: &str) -> Option<u32> {
    let rest = address.strip_prefix(prefix)?;
    rest.strip_suffix(']')?.parse::<u32>().ok()
}

/// V/D/J-only segment used by typed trace addresses.
///
/// This avoids constructing invalid addresses like
/// `sample_allele.np1` or `trim.np2_3` while still converting back to
/// the engine-wide [`Segment`] enum at the boundary.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum VdjSegment {
    V,
    D,
    J,
}

impl VdjSegment {
    const fn suffix(self) -> &'static str {
        match self {
            Self::V => "v",
            Self::D => "d",
            Self::J => "j",
        }
    }
}

impl TryFrom<Segment> for VdjSegment {
    type Error = ();

    fn try_from(value: Segment) -> Result<Self, Self::Error> {
        match value {
            Segment::V => Ok(Self::V),
            Segment::D => Ok(Self::D),
            Segment::J => Ok(Self::J),
            Segment::Np1 | Segment::Np2 => Err(()),
        }
    }
}

impl From<VdjSegment> for Segment {
    fn from(value: VdjSegment) -> Self {
        match value {
            VdjSegment::V => Segment::V,
            VdjSegment::D => Segment::D,
            VdjSegment::J => Segment::J,
        }
    }
}

/// NP1/NP2-only segment used by typed trace addresses.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum NpSegment {
    Np1,
    Np2,
}

impl NpSegment {
    const fn suffix(self) -> &'static str {
        match self {
            Self::Np1 => "np1",
            Self::Np2 => "np2",
        }
    }
}

impl TryFrom<Segment> for NpSegment {
    type Error = ();

    fn try_from(value: Segment) -> Result<Self, Self::Error> {
        match value {
            Segment::Np1 => Ok(Self::Np1),
            Segment::Np2 => Ok(Self::Np2),
            Segment::V | Segment::D | Segment::J => Err(()),
        }
    }
}

impl From<NpSegment> for Segment {
    fn from(value: NpSegment) -> Self {
        match value {
            NpSegment::Np1 => Segment::Np1,
            NpSegment::Np2 => Segment::Np2,
        }
    }
}

/// Generic 5'/3' end for non-allele trace choices such as end-loss.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum PrimeEnd {
    Five,
    Three,
}

/// The four V(D)J coding-end junction sides at which a
/// `PAdditionPass` can extend the assembled sequence with
/// templated palindromic (P-)nucleotides. Order in the enum
/// matches biological position along the V→J axis so iteration
/// produces the natural left-to-right surface for AIRR
/// projection and lowering.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum PEnd {
    /// V's 3' coding end — P bytes sit between V and NP1.
    V3,
    /// D's 5' coding end — P bytes sit between NP1 and D.
    D5,
    /// D's 3' coding end — P bytes sit between D and NP2.
    D3,
    /// J's 5' coding end — P bytes sit between NP2 and J.
    J5,
}

impl PEnd {
    /// On-disk suffix used in trace addresses and pass names.
    /// Matches the existing `trim.{v_3,d_5,d_3,j_5}` convention.
    pub const fn suffix(self) -> &'static str {
        match self {
            Self::V3 => "v_3",
            Self::D5 => "d_5",
            Self::D3 => "d_3",
            Self::J5 => "j_5",
        }
    }

    /// The V/D/J segment from whose coding flank the P-bytes
    /// derive (palindrome source). V3 → V, D5/D3 → D, J5 → J.
    pub const fn source_segment(self) -> Segment {
        match self {
            Self::V3 => Segment::V,
            Self::D5 | Self::D3 => Segment::D,
            Self::J5 => Segment::J,
        }
    }

    /// All four ends in canonical V→J biological order.
    pub const fn all() -> [PEnd; 4] {
        [Self::V3, Self::D5, Self::D3, Self::J5]
    }
}

impl From<TrimEnd> for PrimeEnd {
    fn from(value: TrimEnd) -> Self {
        match value {
            TrimEnd::Five => Self::Five,
            TrimEnd::Three => Self::Three,
        }
    }
}

/// Typed form of every built-in concrete stochastic-choice address.
///
/// Persisted traces still store strings; this enum gives internal
/// code a typed surface for new call sites and a single round-trip
/// parser for replay / inspection paths. Pattern specs such as
/// `mutate.s5f.site[0..n]` are intentionally not represented here:
/// they describe declared choice families, not concrete choices made
/// during one run.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ChoiceAddress {
    SampleAllele(VdjSegment),
    Trim { segment: VdjSegment, end: TrimEnd },
    NpLength(NpSegment),
    NpBase { segment: NpSegment, index: u32 },
    MutateUniformCount,
    MutateUniformSite(u32),
    MutateUniformBase(u32),
    MutateS5fCount,
    MutateS5fSite(u32),
    MutateS5fBase(u32),
    CorruptPcrCount,
    CorruptPcrSite(u32),
    CorruptPcrBase(u32),
    CorruptQualityCount,
    CorruptQualitySite(u32),
    CorruptQualityBase(u32),
    CorruptContaminantApplied,
    CorruptContaminantBase(u32),
    CorruptIndelCount,
    CorruptIndelKind(u32),
    CorruptIndelSite(u32),
    CorruptIndelBase(u32),
    CorruptNsCount,
    CorruptNsSite(u32),
    CorruptEndLoss(PrimeEnd),
    CorruptRevCompApplied,
    /// Per-simulation Bool: `true` iff the D segment is to be
    /// assembled in reverse-complement orientation. Singleton choice
    /// (one record per simulation). See [`SAMPLE_ALLELE_D_INVERTED`]
    /// for the on-disk spelling.
    SampleAlleleDInverted,
    /// Per-simulation Bool: `true` iff `ReceptorRevisionPass` replaces
    /// the V segment for this simulation. Always recorded — even on
    /// `false`, so the trace carries the outcome unambiguously.
    /// See [`RECEPTOR_REVISION_APPLIED`].
    ReceptorRevisionApplied,
    /// AlleleId: the replacement V allele chosen on `applied=true`.
    /// Recorded only on `applied=true`. See
    /// [`RECEPTOR_REVISION_V_ALLELE`].
    ReceptorRevisionVAllele,
    /// Int: the 3' trim applied to the replacement V allele on
    /// `applied=true`. Recorded only on `applied=true`. See
    /// [`RECEPTOR_REVISION_V_TRIM_3`].
    ReceptorRevisionVTrim3,
    /// Int: the R1 read length for the paired-end projection.
    /// Recorded once per simulation by `PairedEndSamplingPass`.
    PairedEndR1Length,
    /// Int: the R2 read length for the paired-end projection.
    PairedEndR2Length,
    /// Int: the fragment insert size for the paired-end projection.
    PairedEndInsertSize,
    /// Int: per-end P-nucleotide length sampled by a
    /// `PAdditionPass`. P-bases themselves are deterministic
    /// from `(allele, trim, orientation, length)` so only the
    /// length needs a trace address.
    PLength { end: PEnd },
    /// Haplotype (0/1): the chromosome drawn once per rearrangement by
    /// `SampleHaplotypePass` for a phased genotype. V/D/J read it back.
    SampleHaplotype,
    /// GeneId: the gene chosen within the drawn chromosome for a segment
    /// by `SampleGeneAllelePass`.
    SampleGene(VdjSegment),
    /// AlleleId: the within-slot allele draw, recorded only when a gene
    /// slot carries more than one copy (single-copy slots are
    /// deterministic and record nothing here).
    SampleAlleleInSlot(VdjSegment),
}

impl ChoiceAddress {
    pub fn parse(address: &str) -> Option<Self> {
        address.parse().ok()
    }
}

impl fmt::Display for ChoiceAddress {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::SampleAllele(segment) => {
                write!(f, "sample_allele.{}", segment.suffix())
            }
            Self::Trim { segment, end } => {
                let end = match end {
                    TrimEnd::Five => "5",
                    TrimEnd::Three => "3",
                };
                write!(f, "trim.{}_{}", segment.suffix(), end)
            }
            Self::NpLength(segment) => write!(f, "np.{}.length", segment.suffix()),
            Self::NpBase { segment, index } => {
                write!(f, "np.{}.bases[{}]", segment.suffix(), index)
            }
            Self::MutateUniformCount => f.write_str(MUTATE_UNIFORM_COUNT),
            Self::MutateUniformSite(index) => {
                write!(f, "{}{}]", MUTATE_UNIFORM_SITE_PREFIX, index)
            }
            Self::MutateUniformBase(index) => {
                write!(f, "{}{}]", MUTATE_UNIFORM_BASE_PREFIX, index)
            }
            Self::MutateS5fCount => f.write_str(MUTATE_S5F_COUNT),
            Self::MutateS5fSite(index) => {
                write!(f, "{}{}]", MUTATE_S5F_SITE_PREFIX, index)
            }
            Self::MutateS5fBase(index) => {
                write!(f, "{}{}]", MUTATE_S5F_BASE_PREFIX, index)
            }
            Self::CorruptPcrCount => f.write_str(CORRUPT_PCR_COUNT),
            Self::CorruptPcrSite(index) => write!(f, "{}{}]", CORRUPT_PCR_SITE_PREFIX, index),
            Self::CorruptPcrBase(index) => write!(f, "{}{}]", CORRUPT_PCR_BASE_PREFIX, index),
            Self::CorruptQualityCount => f.write_str(CORRUPT_QUALITY_COUNT),
            Self::CorruptQualitySite(index) => {
                write!(f, "{}{}]", CORRUPT_QUALITY_SITE_PREFIX, index)
            }
            Self::CorruptQualityBase(index) => {
                write!(f, "{}{}]", CORRUPT_QUALITY_BASE_PREFIX, index)
            }
            Self::CorruptContaminantApplied => f.write_str(CORRUPT_CONTAMINANT_APPLIED),
            Self::CorruptContaminantBase(index) => {
                write!(f, "{}{}]", CORRUPT_CONTAMINANT_BASES_PREFIX, index)
            }
            Self::CorruptIndelCount => f.write_str(CORRUPT_INDEL_COUNT),
            Self::CorruptIndelKind(index) => {
                write!(f, "{}{}]", CORRUPT_INDEL_KIND_PREFIX, index)
            }
            Self::CorruptIndelSite(index) => {
                write!(f, "{}{}]", CORRUPT_INDEL_SITE_PREFIX, index)
            }
            Self::CorruptIndelBase(index) => {
                write!(f, "{}{}]", CORRUPT_INDEL_BASE_PREFIX, index)
            }
            Self::CorruptNsCount => f.write_str(CORRUPT_NS_COUNT),
            Self::CorruptNsSite(index) => write!(f, "{}{}]", CORRUPT_NS_SITE_PREFIX, index),
            Self::CorruptEndLoss(PrimeEnd::Five) => f.write_str(CORRUPT_END_LOSS_5),
            Self::CorruptEndLoss(PrimeEnd::Three) => f.write_str(CORRUPT_END_LOSS_3),
            Self::CorruptRevCompApplied => f.write_str(CORRUPT_REV_COMP_APPLIED),
            Self::SampleAlleleDInverted => f.write_str(SAMPLE_ALLELE_D_INVERTED),
            Self::ReceptorRevisionApplied => f.write_str(RECEPTOR_REVISION_APPLIED),
            Self::ReceptorRevisionVAllele => f.write_str(RECEPTOR_REVISION_V_ALLELE),
            Self::ReceptorRevisionVTrim3 => f.write_str(RECEPTOR_REVISION_V_TRIM_3),
            Self::PairedEndR1Length => f.write_str(PAIRED_END_R1_LENGTH),
            Self::PairedEndR2Length => f.write_str(PAIRED_END_R2_LENGTH),
            Self::PairedEndInsertSize => f.write_str(PAIRED_END_INSERT_SIZE),
            Self::PLength { end } => write!(f, "p.{}.length", end.suffix()),
            Self::SampleHaplotype => f.write_str(SAMPLE_HAPLOTYPE),
            Self::SampleGene(segment) => write!(f, "sample_gene.{}", segment.suffix()),
            Self::SampleAlleleInSlot(segment) => {
                write!(f, "sample_allele_in_slot.{}", segment.suffix())
            }
        }
    }
}

impl From<ChoiceAddress> for String {
    fn from(value: ChoiceAddress) -> Self {
        value.to_string()
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct ChoiceAddressParseError;

impl fmt::Display for ChoiceAddressParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("unrecognized choice address")
    }
}

impl std::error::Error for ChoiceAddressParseError {}

impl FromStr for ChoiceAddress {
    type Err = ChoiceAddressParseError;

    fn from_str(address: &str) -> Result<Self, Self::Err> {
        parse_choice_address(address).ok_or(ChoiceAddressParseError)
    }
}

fn parse_choice_address(address: &str) -> Option<ChoiceAddress> {
    let exact = match address {
        SAMPLE_ALLELE_V => Some(ChoiceAddress::SampleAllele(VdjSegment::V)),
        SAMPLE_ALLELE_D => Some(ChoiceAddress::SampleAllele(VdjSegment::D)),
        SAMPLE_ALLELE_J => Some(ChoiceAddress::SampleAllele(VdjSegment::J)),
        TRIM_V_5 => Some(ChoiceAddress::Trim {
            segment: VdjSegment::V,
            end: TrimEnd::Five,
        }),
        TRIM_V_3 => Some(ChoiceAddress::Trim {
            segment: VdjSegment::V,
            end: TrimEnd::Three,
        }),
        TRIM_D_5 => Some(ChoiceAddress::Trim {
            segment: VdjSegment::D,
            end: TrimEnd::Five,
        }),
        TRIM_D_3 => Some(ChoiceAddress::Trim {
            segment: VdjSegment::D,
            end: TrimEnd::Three,
        }),
        TRIM_J_5 => Some(ChoiceAddress::Trim {
            segment: VdjSegment::J,
            end: TrimEnd::Five,
        }),
        TRIM_J_3 => Some(ChoiceAddress::Trim {
            segment: VdjSegment::J,
            end: TrimEnd::Three,
        }),
        NP1_LENGTH => Some(ChoiceAddress::NpLength(NpSegment::Np1)),
        NP2_LENGTH => Some(ChoiceAddress::NpLength(NpSegment::Np2)),
        MUTATE_UNIFORM_COUNT => Some(ChoiceAddress::MutateUniformCount),
        MUTATE_S5F_COUNT => Some(ChoiceAddress::MutateS5fCount),
        CORRUPT_PCR_COUNT => Some(ChoiceAddress::CorruptPcrCount),
        CORRUPT_QUALITY_COUNT => Some(ChoiceAddress::CorruptQualityCount),
        CORRUPT_CONTAMINANT_APPLIED => Some(ChoiceAddress::CorruptContaminantApplied),
        CORRUPT_INDEL_COUNT => Some(ChoiceAddress::CorruptIndelCount),
        CORRUPT_NS_COUNT => Some(ChoiceAddress::CorruptNsCount),
        CORRUPT_END_LOSS_5 => Some(ChoiceAddress::CorruptEndLoss(PrimeEnd::Five)),
        CORRUPT_END_LOSS_3 => Some(ChoiceAddress::CorruptEndLoss(PrimeEnd::Three)),
        CORRUPT_REV_COMP_APPLIED => Some(ChoiceAddress::CorruptRevCompApplied),
        SAMPLE_ALLELE_D_INVERTED => Some(ChoiceAddress::SampleAlleleDInverted),
        RECEPTOR_REVISION_APPLIED => Some(ChoiceAddress::ReceptorRevisionApplied),
        RECEPTOR_REVISION_V_ALLELE => Some(ChoiceAddress::ReceptorRevisionVAllele),
        RECEPTOR_REVISION_V_TRIM_3 => Some(ChoiceAddress::ReceptorRevisionVTrim3),
        PAIRED_END_R1_LENGTH => Some(ChoiceAddress::PairedEndR1Length),
        PAIRED_END_R2_LENGTH => Some(ChoiceAddress::PairedEndR2Length),
        PAIRED_END_INSERT_SIZE => Some(ChoiceAddress::PairedEndInsertSize),
        P_V3_LENGTH => Some(ChoiceAddress::PLength { end: PEnd::V3 }),
        P_D5_LENGTH => Some(ChoiceAddress::PLength { end: PEnd::D5 }),
        P_D3_LENGTH => Some(ChoiceAddress::PLength { end: PEnd::D3 }),
        P_J5_LENGTH => Some(ChoiceAddress::PLength { end: PEnd::J5 }),
        SAMPLE_HAPLOTYPE => Some(ChoiceAddress::SampleHaplotype),
        "sample_gene.v" => Some(ChoiceAddress::SampleGene(VdjSegment::V)),
        "sample_gene.d" => Some(ChoiceAddress::SampleGene(VdjSegment::D)),
        "sample_gene.j" => Some(ChoiceAddress::SampleGene(VdjSegment::J)),
        "sample_allele_in_slot.v" => Some(ChoiceAddress::SampleAlleleInSlot(VdjSegment::V)),
        "sample_allele_in_slot.d" => Some(ChoiceAddress::SampleAlleleInSlot(VdjSegment::D)),
        "sample_allele_in_slot.j" => Some(ChoiceAddress::SampleAlleleInSlot(VdjSegment::J)),
        _ => None,
    };
    if exact.is_some() {
        return exact;
    }

    if let Some(index) = parse_indexed(address, NP1_BASES_INDEX_PREFIX) {
        return Some(ChoiceAddress::NpBase {
            segment: NpSegment::Np1,
            index,
        });
    }
    if let Some(index) = parse_indexed(address, NP2_BASES_INDEX_PREFIX) {
        return Some(ChoiceAddress::NpBase {
            segment: NpSegment::Np2,
            index,
        });
    }
    if let Some(index) = parse_indexed(address, MUTATE_UNIFORM_SITE_PREFIX) {
        return Some(ChoiceAddress::MutateUniformSite(index));
    }
    if let Some(index) = parse_indexed(address, MUTATE_UNIFORM_BASE_PREFIX) {
        return Some(ChoiceAddress::MutateUniformBase(index));
    }
    if let Some(index) = parse_indexed(address, MUTATE_S5F_SITE_PREFIX) {
        return Some(ChoiceAddress::MutateS5fSite(index));
    }
    if let Some(index) = parse_indexed(address, MUTATE_S5F_BASE_PREFIX) {
        return Some(ChoiceAddress::MutateS5fBase(index));
    }
    if let Some(index) = parse_indexed(address, CORRUPT_PCR_SITE_PREFIX) {
        return Some(ChoiceAddress::CorruptPcrSite(index));
    }
    if let Some(index) = parse_indexed(address, CORRUPT_PCR_BASE_PREFIX) {
        return Some(ChoiceAddress::CorruptPcrBase(index));
    }
    if let Some(index) = parse_indexed(address, CORRUPT_QUALITY_SITE_PREFIX) {
        return Some(ChoiceAddress::CorruptQualitySite(index));
    }
    if let Some(index) = parse_indexed(address, CORRUPT_QUALITY_BASE_PREFIX) {
        return Some(ChoiceAddress::CorruptQualityBase(index));
    }
    if let Some(index) = parse_indexed(address, CORRUPT_CONTAMINANT_BASES_PREFIX) {
        return Some(ChoiceAddress::CorruptContaminantBase(index));
    }
    if let Some(index) = parse_indexed(address, CORRUPT_INDEL_KIND_PREFIX) {
        return Some(ChoiceAddress::CorruptIndelKind(index));
    }
    if let Some(index) = parse_indexed(address, CORRUPT_INDEL_SITE_PREFIX) {
        return Some(ChoiceAddress::CorruptIndelSite(index));
    }
    if let Some(index) = parse_indexed(address, CORRUPT_INDEL_BASE_PREFIX) {
        return Some(ChoiceAddress::CorruptIndelBase(index));
    }
    if let Some(index) = parse_indexed(address, CORRUPT_NS_SITE_PREFIX) {
        return Some(ChoiceAddress::CorruptNsSite(index));
    }

    None
}

/// Typed form of declared stochastic-choice address families.
///
/// A [`ChoiceAddress`] names one concrete draw from one run, while
/// `ChoiceAddressPattern` names the family a pass declares at
/// compile/report time. Singleton choices (for example
/// `mutate.s5f.count`) appear as one-member families; indexed choices
/// display with the existing `[0..n]` pattern string.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ChoiceAddressPattern {
    SampleAllele(VdjSegment),
    Trim { segment: VdjSegment, end: TrimEnd },
    NpLength(NpSegment),
    NpBase(NpSegment),
    MutateUniformCount,
    MutateUniformSite,
    MutateUniformBase,
    MutateS5fCount,
    MutateS5fSite,
    MutateS5fBase,
    CorruptPcrCount,
    CorruptPcrSite,
    CorruptPcrBase,
    CorruptQualityCount,
    CorruptQualitySite,
    CorruptQualityBase,
    CorruptContaminantApplied,
    CorruptContaminantBase,
    CorruptIndelCount,
    CorruptIndelKind,
    CorruptIndelSite,
    CorruptIndelBase,
    CorruptNsCount,
    CorruptNsSite,
    CorruptEndLoss(PrimeEnd),
    CorruptRevCompApplied,
    /// Singleton-family mirror of
    /// [`ChoiceAddress::SampleAlleleDInverted`]. The
    /// `InvertDPass` declares this pattern from
    /// [`crate::pass::Pass::declared_choice_patterns`] so the
    /// schedule analyser and trace-replay validator recognise the
    /// address.
    SampleAlleleDInverted,
    /// Singleton-family mirror of
    /// [`ChoiceAddress::ReceptorRevisionApplied`]. Declared by
    /// `ReceptorRevisionPass` regardless of whether replacement
    /// fires, since the Bool is always recorded.
    ReceptorRevisionApplied,
    /// Singleton-family mirror of
    /// [`ChoiceAddress::ReceptorRevisionVAllele`]. Declared even
    /// though the record is conditional on `applied=true`; the
    /// schedule analyser treats it as a *potential* draw the pass
    /// may make.
    ReceptorRevisionVAllele,
    /// Singleton-family mirror of
    /// [`ChoiceAddress::ReceptorRevisionVTrim3`].
    ReceptorRevisionVTrim3,
    /// Singleton-family mirror of
    /// [`ChoiceAddress::PairedEndR1Length`]. Declared by
    /// `PairedEndSamplingPass` so the schedule analyser and
    /// trace-replay validator recognise the address.
    PairedEndR1Length,
    /// Singleton-family mirror of
    /// [`ChoiceAddress::PairedEndR2Length`].
    PairedEndR2Length,
    /// Singleton-family mirror of
    /// [`ChoiceAddress::PairedEndInsertSize`].
    PairedEndInsertSize,
    /// Singleton-family mirror of
    /// [`ChoiceAddress::PLength`]. One pattern instance per
    /// `PEnd` — declared by `PAdditionPass`.
    PLength { end: PEnd },
    /// Singleton-family mirror of [`ChoiceAddress::SampleHaplotype`].
    /// Declared by `SampleHaplotypePass`.
    SampleHaplotype,
    /// Family mirror of [`ChoiceAddress::SampleGene`]. Declared by
    /// `SampleGeneAllelePass`.
    SampleGene(VdjSegment),
    /// Family mirror of [`ChoiceAddress::SampleAlleleInSlot`]. Declared
    /// by `SampleGeneAllelePass` (a potential draw on multi-copy slots).
    SampleAlleleInSlot(VdjSegment),
}

impl ChoiceAddressPattern {
    pub fn parse(address: &str) -> Option<Self> {
        address.parse().ok()
    }
}

impl fmt::Display for ChoiceAddressPattern {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::SampleAllele(segment) => ChoiceAddress::SampleAllele(segment).fmt(f),
            Self::Trim { segment, end } => ChoiceAddress::Trim { segment, end }.fmt(f),
            Self::NpLength(segment) => ChoiceAddress::NpLength(segment).fmt(f),
            Self::NpBase(NpSegment::Np1) => f.write_str(NP1_BASES_PATTERN),
            Self::NpBase(NpSegment::Np2) => f.write_str(NP2_BASES_PATTERN),
            Self::MutateUniformCount => f.write_str(MUTATE_UNIFORM_COUNT),
            Self::MutateUniformSite => f.write_str(MUTATE_UNIFORM_SITE_PATTERN),
            Self::MutateUniformBase => f.write_str(MUTATE_UNIFORM_BASE_PATTERN),
            Self::MutateS5fCount => f.write_str(MUTATE_S5F_COUNT),
            Self::MutateS5fSite => f.write_str(MUTATE_S5F_SITE_PATTERN),
            Self::MutateS5fBase => f.write_str(MUTATE_S5F_BASE_PATTERN),
            Self::CorruptPcrCount => f.write_str(CORRUPT_PCR_COUNT),
            Self::CorruptPcrSite => f.write_str(CORRUPT_PCR_SITE_PATTERN),
            Self::CorruptPcrBase => f.write_str(CORRUPT_PCR_BASE_PATTERN),
            Self::CorruptQualityCount => f.write_str(CORRUPT_QUALITY_COUNT),
            Self::CorruptQualitySite => f.write_str(CORRUPT_QUALITY_SITE_PATTERN),
            Self::CorruptQualityBase => f.write_str(CORRUPT_QUALITY_BASE_PATTERN),
            Self::CorruptContaminantApplied => f.write_str(CORRUPT_CONTAMINANT_APPLIED),
            Self::CorruptContaminantBase => f.write_str(CORRUPT_CONTAMINANT_BASES_PATTERN),
            Self::CorruptIndelCount => f.write_str(CORRUPT_INDEL_COUNT),
            Self::CorruptIndelKind => f.write_str(CORRUPT_INDEL_KIND_PATTERN),
            Self::CorruptIndelSite => f.write_str(CORRUPT_INDEL_SITE_PATTERN),
            Self::CorruptIndelBase => f.write_str(CORRUPT_INDEL_BASE_PATTERN),
            Self::CorruptNsCount => f.write_str(CORRUPT_NS_COUNT),
            Self::CorruptNsSite => f.write_str(CORRUPT_NS_SITE_PATTERN),
            Self::CorruptEndLoss(end) => ChoiceAddress::CorruptEndLoss(end).fmt(f),
            Self::CorruptRevCompApplied => f.write_str(CORRUPT_REV_COMP_APPLIED),
            Self::SampleAlleleDInverted => f.write_str(SAMPLE_ALLELE_D_INVERTED),
            Self::ReceptorRevisionApplied => f.write_str(RECEPTOR_REVISION_APPLIED),
            Self::ReceptorRevisionVAllele => f.write_str(RECEPTOR_REVISION_V_ALLELE),
            Self::ReceptorRevisionVTrim3 => f.write_str(RECEPTOR_REVISION_V_TRIM_3),
            Self::PairedEndR1Length => f.write_str(PAIRED_END_R1_LENGTH),
            Self::PairedEndR2Length => f.write_str(PAIRED_END_R2_LENGTH),
            Self::PairedEndInsertSize => f.write_str(PAIRED_END_INSERT_SIZE),
            Self::PLength { end } => ChoiceAddress::PLength { end }.fmt(f),
            Self::SampleHaplotype => f.write_str(SAMPLE_HAPLOTYPE),
            Self::SampleGene(segment) => ChoiceAddress::SampleGene(segment).fmt(f),
            Self::SampleAlleleInSlot(segment) => {
                ChoiceAddress::SampleAlleleInSlot(segment).fmt(f)
            }
        }
    }
}

impl From<ChoiceAddressPattern> for String {
    fn from(value: ChoiceAddressPattern) -> Self {
        value.to_string()
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct ChoiceAddressPatternParseError;

impl fmt::Display for ChoiceAddressPatternParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("unrecognized choice address pattern")
    }
}

impl std::error::Error for ChoiceAddressPatternParseError {}

impl FromStr for ChoiceAddressPattern {
    type Err = ChoiceAddressPatternParseError;

    fn from_str(address: &str) -> Result<Self, Self::Err> {
        parse_choice_address_pattern(address).ok_or(ChoiceAddressPatternParseError)
    }
}

fn parse_choice_address_pattern(address: &str) -> Option<ChoiceAddressPattern> {
    let exact = match address {
        SAMPLE_ALLELE_V => Some(ChoiceAddressPattern::SampleAllele(VdjSegment::V)),
        SAMPLE_ALLELE_D => Some(ChoiceAddressPattern::SampleAllele(VdjSegment::D)),
        SAMPLE_ALLELE_J => Some(ChoiceAddressPattern::SampleAllele(VdjSegment::J)),
        TRIM_V_5 => Some(ChoiceAddressPattern::Trim {
            segment: VdjSegment::V,
            end: TrimEnd::Five,
        }),
        TRIM_V_3 => Some(ChoiceAddressPattern::Trim {
            segment: VdjSegment::V,
            end: TrimEnd::Three,
        }),
        TRIM_D_5 => Some(ChoiceAddressPattern::Trim {
            segment: VdjSegment::D,
            end: TrimEnd::Five,
        }),
        TRIM_D_3 => Some(ChoiceAddressPattern::Trim {
            segment: VdjSegment::D,
            end: TrimEnd::Three,
        }),
        TRIM_J_5 => Some(ChoiceAddressPattern::Trim {
            segment: VdjSegment::J,
            end: TrimEnd::Five,
        }),
        TRIM_J_3 => Some(ChoiceAddressPattern::Trim {
            segment: VdjSegment::J,
            end: TrimEnd::Three,
        }),
        NP1_LENGTH => Some(ChoiceAddressPattern::NpLength(NpSegment::Np1)),
        NP2_LENGTH => Some(ChoiceAddressPattern::NpLength(NpSegment::Np2)),
        NP1_BASES_PATTERN => Some(ChoiceAddressPattern::NpBase(NpSegment::Np1)),
        NP2_BASES_PATTERN => Some(ChoiceAddressPattern::NpBase(NpSegment::Np2)),
        MUTATE_UNIFORM_COUNT => Some(ChoiceAddressPattern::MutateUniformCount),
        MUTATE_UNIFORM_SITE_PATTERN => Some(ChoiceAddressPattern::MutateUniformSite),
        MUTATE_UNIFORM_BASE_PATTERN => Some(ChoiceAddressPattern::MutateUniformBase),
        MUTATE_S5F_COUNT => Some(ChoiceAddressPattern::MutateS5fCount),
        MUTATE_S5F_SITE_PATTERN => Some(ChoiceAddressPattern::MutateS5fSite),
        MUTATE_S5F_BASE_PATTERN => Some(ChoiceAddressPattern::MutateS5fBase),
        CORRUPT_PCR_COUNT => Some(ChoiceAddressPattern::CorruptPcrCount),
        CORRUPT_PCR_SITE_PATTERN => Some(ChoiceAddressPattern::CorruptPcrSite),
        CORRUPT_PCR_BASE_PATTERN => Some(ChoiceAddressPattern::CorruptPcrBase),
        CORRUPT_QUALITY_COUNT => Some(ChoiceAddressPattern::CorruptQualityCount),
        CORRUPT_QUALITY_SITE_PATTERN => Some(ChoiceAddressPattern::CorruptQualitySite),
        CORRUPT_QUALITY_BASE_PATTERN => Some(ChoiceAddressPattern::CorruptQualityBase),
        CORRUPT_CONTAMINANT_APPLIED => Some(ChoiceAddressPattern::CorruptContaminantApplied),
        CORRUPT_CONTAMINANT_BASES_PATTERN => Some(ChoiceAddressPattern::CorruptContaminantBase),
        CORRUPT_INDEL_COUNT => Some(ChoiceAddressPattern::CorruptIndelCount),
        CORRUPT_INDEL_KIND_PATTERN => Some(ChoiceAddressPattern::CorruptIndelKind),
        CORRUPT_INDEL_SITE_PATTERN => Some(ChoiceAddressPattern::CorruptIndelSite),
        CORRUPT_INDEL_BASE_PATTERN => Some(ChoiceAddressPattern::CorruptIndelBase),
        CORRUPT_NS_COUNT => Some(ChoiceAddressPattern::CorruptNsCount),
        CORRUPT_NS_SITE_PATTERN => Some(ChoiceAddressPattern::CorruptNsSite),
        CORRUPT_END_LOSS_5 => Some(ChoiceAddressPattern::CorruptEndLoss(PrimeEnd::Five)),
        CORRUPT_END_LOSS_3 => Some(ChoiceAddressPattern::CorruptEndLoss(PrimeEnd::Three)),
        CORRUPT_REV_COMP_APPLIED => Some(ChoiceAddressPattern::CorruptRevCompApplied),
        SAMPLE_ALLELE_D_INVERTED => Some(ChoiceAddressPattern::SampleAlleleDInverted),
        RECEPTOR_REVISION_APPLIED => Some(ChoiceAddressPattern::ReceptorRevisionApplied),
        RECEPTOR_REVISION_V_ALLELE => Some(ChoiceAddressPattern::ReceptorRevisionVAllele),
        RECEPTOR_REVISION_V_TRIM_3 => Some(ChoiceAddressPattern::ReceptorRevisionVTrim3),
        PAIRED_END_R1_LENGTH => Some(ChoiceAddressPattern::PairedEndR1Length),
        PAIRED_END_R2_LENGTH => Some(ChoiceAddressPattern::PairedEndR2Length),
        PAIRED_END_INSERT_SIZE => Some(ChoiceAddressPattern::PairedEndInsertSize),
        P_V3_LENGTH => Some(ChoiceAddressPattern::PLength { end: PEnd::V3 }),
        P_D5_LENGTH => Some(ChoiceAddressPattern::PLength { end: PEnd::D5 }),
        P_D3_LENGTH => Some(ChoiceAddressPattern::PLength { end: PEnd::D3 }),
        P_J5_LENGTH => Some(ChoiceAddressPattern::PLength { end: PEnd::J5 }),
        SAMPLE_HAPLOTYPE => Some(ChoiceAddressPattern::SampleHaplotype),
        "sample_gene.v" => Some(ChoiceAddressPattern::SampleGene(VdjSegment::V)),
        "sample_gene.d" => Some(ChoiceAddressPattern::SampleGene(VdjSegment::D)),
        "sample_gene.j" => Some(ChoiceAddressPattern::SampleGene(VdjSegment::J)),
        "sample_allele_in_slot.v" => Some(ChoiceAddressPattern::SampleAlleleInSlot(VdjSegment::V)),
        "sample_allele_in_slot.d" => Some(ChoiceAddressPattern::SampleAlleleInSlot(VdjSegment::D)),
        "sample_allele_in_slot.j" => Some(ChoiceAddressPattern::SampleAlleleInSlot(VdjSegment::J)),
        _ => None,
    };

    exact
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
    }

    #[test]
    fn parsers_accept_existing_indexed_addresses() {
        assert_eq!(
            parse_indexed("mutate.s5f.base[2]", MUTATE_S5F_BASE_PREFIX),
            Some(2)
        );
        assert_eq!(
            ChoiceAddress::parse("np.np1.bases[17]"),
            Some(ChoiceAddress::NpBase {
                segment: NpSegment::Np1,
                index: 17,
            })
        );
        assert_eq!(
            ChoiceAddress::parse("np.np2.bases[4]"),
            Some(ChoiceAddress::NpBase {
                segment: NpSegment::Np2,
                index: 4,
            })
        );
        assert_eq!(ChoiceAddress::parse("np.np2.bases[x]"), None);
        assert_eq!(
            ChoiceAddress::parse("np.np1.length"),
            Some(ChoiceAddress::NpLength(NpSegment::Np1))
        );
    }

    #[test]
    fn typed_choice_addresses_round_trip_persisted_strings() {
        let cases = [
            (ChoiceAddress::SampleAllele(VdjSegment::V), SAMPLE_ALLELE_V),
            (ChoiceAddress::SampleAllele(VdjSegment::D), SAMPLE_ALLELE_D),
            (ChoiceAddress::SampleAllele(VdjSegment::J), SAMPLE_ALLELE_J),
            (
                ChoiceAddress::Trim {
                    segment: VdjSegment::V,
                    end: TrimEnd::Three,
                },
                TRIM_V_3,
            ),
            (
                ChoiceAddress::Trim {
                    segment: VdjSegment::D,
                    end: TrimEnd::Five,
                },
                TRIM_D_5,
            ),
            (
                ChoiceAddress::Trim {
                    segment: VdjSegment::J,
                    end: TrimEnd::Five,
                },
                TRIM_J_5,
            ),
            (ChoiceAddress::NpLength(NpSegment::Np1), NP1_LENGTH),
            (
                ChoiceAddress::NpBase {
                    segment: NpSegment::Np2,
                    index: 17,
                },
                "np.np2.bases[17]",
            ),
            (ChoiceAddress::MutateUniformCount, MUTATE_UNIFORM_COUNT),
            (
                ChoiceAddress::MutateUniformSite(2),
                "mutate.uniform.site[2]",
            ),
            (
                ChoiceAddress::MutateUniformBase(3),
                "mutate.uniform.base[3]",
            ),
            (ChoiceAddress::MutateS5fCount, MUTATE_S5F_COUNT),
            (ChoiceAddress::MutateS5fSite(4), "mutate.s5f.site[4]"),
            (ChoiceAddress::MutateS5fBase(5), "mutate.s5f.base[5]"),
            (ChoiceAddress::CorruptPcrCount, CORRUPT_PCR_COUNT),
            (
                ChoiceAddress::CorruptPcrSite(6),
                "corrupt.pcr.error_site[6]",
            ),
            (
                ChoiceAddress::CorruptPcrBase(7),
                "corrupt.pcr.error_base[7]",
            ),
            (ChoiceAddress::CorruptQualityCount, CORRUPT_QUALITY_COUNT),
            (
                ChoiceAddress::CorruptQualitySite(8),
                "corrupt.quality.error_site[8]",
            ),
            (
                ChoiceAddress::CorruptQualityBase(9),
                "corrupt.quality.error_base[9]",
            ),
            (
                ChoiceAddress::CorruptContaminantApplied,
                CORRUPT_CONTAMINANT_APPLIED,
            ),
            (
                ChoiceAddress::CorruptContaminantBase(10),
                "corrupt.contaminant.bases[10]",
            ),
            (ChoiceAddress::CorruptIndelCount, CORRUPT_INDEL_COUNT),
            (
                ChoiceAddress::CorruptIndelKind(11),
                "corrupt.indel.kind[11]",
            ),
            (
                ChoiceAddress::CorruptIndelSite(12),
                "corrupt.indel.site[12]",
            ),
            (
                ChoiceAddress::CorruptIndelBase(13),
                "corrupt.indel.base[13]",
            ),
            (ChoiceAddress::CorruptNsCount, CORRUPT_NS_COUNT),
            (ChoiceAddress::CorruptNsSite(14), "corrupt.ns.site[14]"),
            (
                ChoiceAddress::CorruptEndLoss(PrimeEnd::Five),
                CORRUPT_END_LOSS_5,
            ),
            (
                ChoiceAddress::CorruptEndLoss(PrimeEnd::Three),
                CORRUPT_END_LOSS_3,
            ),
            (
                ChoiceAddress::CorruptRevCompApplied,
                CORRUPT_REV_COMP_APPLIED,
            ),
            (
                ChoiceAddress::SampleAlleleDInverted,
                SAMPLE_ALLELE_D_INVERTED,
            ),
            (
                ChoiceAddress::ReceptorRevisionApplied,
                RECEPTOR_REVISION_APPLIED,
            ),
            (
                ChoiceAddress::ReceptorRevisionVAllele,
                RECEPTOR_REVISION_V_ALLELE,
            ),
            (
                ChoiceAddress::ReceptorRevisionVTrim3,
                RECEPTOR_REVISION_V_TRIM_3,
            ),
            (ChoiceAddress::PairedEndR1Length, PAIRED_END_R1_LENGTH),
            (ChoiceAddress::PairedEndR2Length, PAIRED_END_R2_LENGTH),
            (
                ChoiceAddress::PairedEndInsertSize,
                PAIRED_END_INSERT_SIZE,
            ),
        ];

        for (typed, raw) in cases {
            assert_eq!(typed.to_string(), raw);
            assert_eq!(ChoiceAddress::parse(raw), Some(typed));
            assert_eq!(raw.parse::<ChoiceAddress>().unwrap(), typed);
        }
    }

    #[test]
    fn typed_choice_address_rejects_patterns_and_unknown_strings() {
        for raw in [
            "mutate.s5f.site[0..n]",
            "np.np1.bases[x]",
            "assemble.v",
            "sample_allele.np1",
            "trim.np1_3",
            "custom.choice",
            "",
        ] {
            assert_eq!(ChoiceAddress::parse(raw), None, "raw={raw:?}");
            assert!(raw.parse::<ChoiceAddress>().is_err(), "raw={raw:?}");
        }
    }

    #[test]
    fn typed_segment_conversions_reject_wrong_segment_family() {
        assert_eq!(VdjSegment::try_from(Segment::V), Ok(VdjSegment::V));
        assert_eq!(VdjSegment::try_from(Segment::Np1), Err(()));
        assert_eq!(NpSegment::try_from(Segment::Np2), Ok(NpSegment::Np2));
        assert_eq!(NpSegment::try_from(Segment::J), Err(()));
    }

    #[test]
    fn typed_choice_address_patterns_round_trip_report_strings() {
        let cases = [
            (
                ChoiceAddressPattern::SampleAllele(VdjSegment::V),
                SAMPLE_ALLELE_V,
            ),
            (
                ChoiceAddressPattern::SampleAllele(VdjSegment::D),
                SAMPLE_ALLELE_D,
            ),
            (
                ChoiceAddressPattern::SampleAllele(VdjSegment::J),
                SAMPLE_ALLELE_J,
            ),
            (
                ChoiceAddressPattern::Trim {
                    segment: VdjSegment::V,
                    end: TrimEnd::Five,
                },
                TRIM_V_5,
            ),
            (
                ChoiceAddressPattern::Trim {
                    segment: VdjSegment::D,
                    end: TrimEnd::Three,
                },
                TRIM_D_3,
            ),
            (
                ChoiceAddressPattern::Trim {
                    segment: VdjSegment::J,
                    end: TrimEnd::Five,
                },
                TRIM_J_5,
            ),
            (ChoiceAddressPattern::NpLength(NpSegment::Np1), NP1_LENGTH),
            (ChoiceAddressPattern::NpLength(NpSegment::Np2), NP2_LENGTH),
            (
                ChoiceAddressPattern::NpBase(NpSegment::Np1),
                NP1_BASES_PATTERN,
            ),
            (
                ChoiceAddressPattern::NpBase(NpSegment::Np2),
                NP2_BASES_PATTERN,
            ),
            (
                ChoiceAddressPattern::MutateUniformCount,
                MUTATE_UNIFORM_COUNT,
            ),
            (
                ChoiceAddressPattern::MutateUniformSite,
                MUTATE_UNIFORM_SITE_PATTERN,
            ),
            (
                ChoiceAddressPattern::MutateUniformBase,
                MUTATE_UNIFORM_BASE_PATTERN,
            ),
            (ChoiceAddressPattern::MutateS5fCount, MUTATE_S5F_COUNT),
            (ChoiceAddressPattern::MutateS5fSite, MUTATE_S5F_SITE_PATTERN),
            (ChoiceAddressPattern::MutateS5fBase, MUTATE_S5F_BASE_PATTERN),
            (ChoiceAddressPattern::CorruptPcrCount, CORRUPT_PCR_COUNT),
            (
                ChoiceAddressPattern::CorruptPcrSite,
                CORRUPT_PCR_SITE_PATTERN,
            ),
            (
                ChoiceAddressPattern::CorruptPcrBase,
                CORRUPT_PCR_BASE_PATTERN,
            ),
            (
                ChoiceAddressPattern::CorruptQualityCount,
                CORRUPT_QUALITY_COUNT,
            ),
            (
                ChoiceAddressPattern::CorruptQualitySite,
                CORRUPT_QUALITY_SITE_PATTERN,
            ),
            (
                ChoiceAddressPattern::CorruptQualityBase,
                CORRUPT_QUALITY_BASE_PATTERN,
            ),
            (
                ChoiceAddressPattern::CorruptContaminantApplied,
                CORRUPT_CONTAMINANT_APPLIED,
            ),
            (
                ChoiceAddressPattern::CorruptContaminantBase,
                CORRUPT_CONTAMINANT_BASES_PATTERN,
            ),
            (ChoiceAddressPattern::CorruptIndelCount, CORRUPT_INDEL_COUNT),
            (
                ChoiceAddressPattern::CorruptIndelKind,
                CORRUPT_INDEL_KIND_PATTERN,
            ),
            (
                ChoiceAddressPattern::CorruptIndelSite,
                CORRUPT_INDEL_SITE_PATTERN,
            ),
            (
                ChoiceAddressPattern::CorruptIndelBase,
                CORRUPT_INDEL_BASE_PATTERN,
            ),
            (ChoiceAddressPattern::CorruptNsCount, CORRUPT_NS_COUNT),
            (ChoiceAddressPattern::CorruptNsSite, CORRUPT_NS_SITE_PATTERN),
            (
                ChoiceAddressPattern::CorruptEndLoss(PrimeEnd::Five),
                CORRUPT_END_LOSS_5,
            ),
            (
                ChoiceAddressPattern::CorruptEndLoss(PrimeEnd::Three),
                CORRUPT_END_LOSS_3,
            ),
            (
                ChoiceAddressPattern::CorruptRevCompApplied,
                CORRUPT_REV_COMP_APPLIED,
            ),
            (
                ChoiceAddressPattern::SampleAlleleDInverted,
                SAMPLE_ALLELE_D_INVERTED,
            ),
            (
                ChoiceAddressPattern::ReceptorRevisionApplied,
                RECEPTOR_REVISION_APPLIED,
            ),
            (
                ChoiceAddressPattern::ReceptorRevisionVAllele,
                RECEPTOR_REVISION_V_ALLELE,
            ),
            (
                ChoiceAddressPattern::ReceptorRevisionVTrim3,
                RECEPTOR_REVISION_V_TRIM_3,
            ),
            (
                ChoiceAddressPattern::PairedEndR1Length,
                PAIRED_END_R1_LENGTH,
            ),
            (
                ChoiceAddressPattern::PairedEndR2Length,
                PAIRED_END_R2_LENGTH,
            ),
            (
                ChoiceAddressPattern::PairedEndInsertSize,
                PAIRED_END_INSERT_SIZE,
            ),
        ];

        for (typed, raw) in cases {
            assert_eq!(typed.to_string(), raw);
            assert_eq!(ChoiceAddressPattern::parse(raw), Some(typed));
            assert_eq!(raw.parse::<ChoiceAddressPattern>().unwrap(), typed);
        }
    }

    #[test]
    fn typed_choice_address_pattern_rejects_concrete_indexes_and_unknown_strings() {
        for raw in [
            "mutate.s5f.site[2]",
            "np.np1.bases[0]",
            "corrupt.indel.kind[7]",
            "mutate.s5f.site[x]",
            "assemble.v",
            "custom.choice",
            "",
        ] {
            assert_eq!(ChoiceAddressPattern::parse(raw), None, "raw={raw:?}");
            assert!(raw.parse::<ChoiceAddressPattern>().is_err(), "raw={raw:?}");
        }
    }

    // ── Frozen address vocabulary (compile-fence) ─────────────────
    //
    // Every persisted trace file carries an `address_schema_version`.
    // Bumping the constant is the explicit signal that the on-disk
    // vocabulary has changed; this test makes accidental drift loud.
    //
    // For each `ChoiceAddress` variant we pin:
    //   (a) the exact `Display` string (the on-disk spelling), and
    //   (b) that `Display → parse → Display` round-trips to the same
    //       string (the bidirectional contract).
    //
    // A code change that touches the Display or parse paths must
    // either preserve the pinned strings or bump
    // `ADDRESS_SCHEMA_VERSION` and re-pin the new spellings here.

    fn assert_pinned(addr: ChoiceAddress, expected: &str) {
        let s = addr.to_string();
        assert_eq!(
            s, expected,
            "ChoiceAddress::Display drift detected; if intentional, bump ADDRESS_SCHEMA_VERSION",
        );
        let parsed = ChoiceAddress::parse(expected)
            .unwrap_or_else(|| panic!("frozen address string fails to parse: {expected:?}"));
        assert_eq!(
            parsed, addr,
            "ChoiceAddress::parse drift; round-trip would break replay of committed traces",
        );
        let round = parsed.to_string();
        assert_eq!(round, expected, "Display ∘ parse must equal Display");
    }

    #[test]
    fn frozen_address_spellings_for_choice_address_schema_v1() {
        assert_eq!(
            ADDRESS_SCHEMA_VERSION, 1,
            "if you bumped ADDRESS_SCHEMA_VERSION, also re-pin this test",
        );

        // SampleAllele.{v,d,j}
        assert_pinned(
            ChoiceAddress::SampleAllele(VdjSegment::V),
            "sample_allele.v",
        );
        assert_pinned(
            ChoiceAddress::SampleAllele(VdjSegment::D),
            "sample_allele.d",
        );
        assert_pinned(
            ChoiceAddress::SampleAllele(VdjSegment::J),
            "sample_allele.j",
        );

        // Trim.{segment}_{end}
        for (seg, seg_str) in [
            (VdjSegment::V, "v"),
            (VdjSegment::D, "d"),
            (VdjSegment::J, "j"),
        ] {
            for (end, end_str) in [(TrimEnd::Five, "5"), (TrimEnd::Three, "3")] {
                let expected = format!("trim.{seg_str}_{end_str}");
                assert_pinned(ChoiceAddress::Trim { segment: seg, end }, &expected);
            }
        }

        // NpLength.{np1,np2}
        assert_pinned(
            ChoiceAddress::NpLength(NpSegment::Np1),
            "np.np1.length",
        );
        assert_pinned(
            ChoiceAddress::NpLength(NpSegment::Np2),
            "np.np2.length",
        );

        // NpBase[i] — pin a few representative indices.
        for i in [0u32, 1, 5, 42] {
            assert_pinned(
                ChoiceAddress::NpBase {
                    segment: NpSegment::Np1,
                    index: i,
                },
                &format!("np.np1.bases[{i}]"),
            );
            assert_pinned(
                ChoiceAddress::NpBase {
                    segment: NpSegment::Np2,
                    index: i,
                },
                &format!("np.np2.bases[{i}]"),
            );
        }

        // Mutate kernels.
        assert_pinned(ChoiceAddress::MutateUniformCount, "mutate.uniform.count");
        for i in [0u32, 7, 99] {
            assert_pinned(
                ChoiceAddress::MutateUniformSite(i),
                &format!("mutate.uniform.site[{i}]"),
            );
            assert_pinned(
                ChoiceAddress::MutateUniformBase(i),
                &format!("mutate.uniform.base[{i}]"),
            );
        }
        assert_pinned(ChoiceAddress::MutateS5fCount, "mutate.s5f.count");
        for i in [0u32, 3] {
            assert_pinned(
                ChoiceAddress::MutateS5fSite(i),
                &format!("mutate.s5f.site[{i}]"),
            );
            assert_pinned(
                ChoiceAddress::MutateS5fBase(i),
                &format!("mutate.s5f.base[{i}]"),
            );
        }

        // Corruption passes — PCR, quality, contaminant, ns, end_loss,
        // rev_comp, indel.
        assert_pinned(ChoiceAddress::CorruptPcrCount, "corrupt.pcr.count");
        for i in [0u32, 2] {
            assert_pinned(
                ChoiceAddress::CorruptPcrSite(i),
                &format!("corrupt.pcr.error_site[{i}]"),
            );
            assert_pinned(
                ChoiceAddress::CorruptPcrBase(i),
                &format!("corrupt.pcr.error_base[{i}]"),
            );
        }
        assert_pinned(ChoiceAddress::CorruptQualityCount, "corrupt.quality.count");
        for i in [0u32, 2] {
            assert_pinned(
                ChoiceAddress::CorruptQualitySite(i),
                &format!("corrupt.quality.error_site[{i}]"),
            );
            assert_pinned(
                ChoiceAddress::CorruptQualityBase(i),
                &format!("corrupt.quality.error_base[{i}]"),
            );
        }
        assert_pinned(
            ChoiceAddress::CorruptContaminantApplied,
            "corrupt.contaminant.applied",
        );
        for i in [0u32, 5] {
            assert_pinned(
                ChoiceAddress::CorruptContaminantBase(i),
                &format!("corrupt.contaminant.bases[{i}]"),
            );
        }
        assert_pinned(ChoiceAddress::CorruptIndelCount, "corrupt.indel.count");
        for i in [0u32, 4] {
            assert_pinned(
                ChoiceAddress::CorruptIndelKind(i),
                &format!("corrupt.indel.kind[{i}]"),
            );
            assert_pinned(
                ChoiceAddress::CorruptIndelSite(i),
                &format!("corrupt.indel.site[{i}]"),
            );
            assert_pinned(
                ChoiceAddress::CorruptIndelBase(i),
                &format!("corrupt.indel.base[{i}]"),
            );
        }
        assert_pinned(ChoiceAddress::CorruptNsCount, "corrupt.ns.count");
        for i in [0u32, 3] {
            assert_pinned(
                ChoiceAddress::CorruptNsSite(i),
                &format!("corrupt.ns.site[{i}]"),
            );
        }
        assert_pinned(
            ChoiceAddress::CorruptEndLoss(PrimeEnd::Five),
            "corrupt.end_loss.5",
        );
        assert_pinned(
            ChoiceAddress::CorruptEndLoss(PrimeEnd::Three),
            "corrupt.end_loss.3",
        );
        assert_pinned(
            ChoiceAddress::CorruptRevCompApplied,
            "corrupt.rev_comp.applied",
        );

        // D inversion (Slice C). The on-disk spelling sits under the
        // existing `sample_allele.d` namespace so the persisted
        // address-schema-version stays at v1 — no migration burden
        // for traces emitted by pre-Slice-C engines that simply
        // never wrote this address.
        assert_pinned(
            ChoiceAddress::SampleAlleleDInverted,
            "sample_allele.d.inverted",
        );

        // Receptor revision (Slice C of the receptor-revision roadmap).
        // New top-level `receptor_revision.*` namespace; additive
        // under the v1 vocabulary policy — old traces don't reference
        // these strings, new traces parse on engines that postdate
        // this slice. No ADDRESS_SCHEMA_VERSION bump required.
        assert_pinned(
            ChoiceAddress::ReceptorRevisionApplied,
            "receptor_revision.applied",
        );
        assert_pinned(
            ChoiceAddress::ReceptorRevisionVAllele,
            "receptor_revision.v_allele",
        );
        assert_pinned(
            ChoiceAddress::ReceptorRevisionVTrim3,
            "receptor_revision.v_trim_3",
        );

        // Paired-end / read layout (Slice C of the paired-end
        // roadmap). New top-level `paired_end.*` namespace; same
        // additive policy as receptor revision — old traces don't
        // reference these strings, no ADDRESS_SCHEMA_VERSION bump.
        assert_pinned(
            ChoiceAddress::PairedEndR1Length,
            "paired_end.r1_length",
        );
        assert_pinned(
            ChoiceAddress::PairedEndR2Length,
            "paired_end.r2_length",
        );
        assert_pinned(
            ChoiceAddress::PairedEndInsertSize,
            "paired_end.insert_size",
        );

        // P-nucleotide length per end (palindromic addition
        // slice). New top-level `p.*` namespace; same additive
        // policy as receptor revision / paired-end — old traces
        // don't reference these strings, no
        // ADDRESS_SCHEMA_VERSION bump.
        assert_pinned(ChoiceAddress::PLength { end: PEnd::V3 }, "p.v_3.length");
        assert_pinned(ChoiceAddress::PLength { end: PEnd::D5 }, "p.d_5.length");
        assert_pinned(ChoiceAddress::PLength { end: PEnd::D3 }, "p.d_3.length");
        assert_pinned(ChoiceAddress::PLength { end: PEnd::J5 }, "p.j_5.length");

        // Phased genotype (genotype-modeling PR1). New top-level
        // `sample_haplotype` + `sample_gene.*` + `sample_allele_in_slot.*`
        // namespaces; same additive policy as receptor revision /
        // paired-end — old traces don't reference these strings, no
        // ADDRESS_SCHEMA_VERSION bump.
        assert_pinned(ChoiceAddress::SampleHaplotype, "sample_haplotype");
        assert_pinned(ChoiceAddress::SampleGene(VdjSegment::V), "sample_gene.v");
        assert_pinned(ChoiceAddress::SampleGene(VdjSegment::D), "sample_gene.d");
        assert_pinned(ChoiceAddress::SampleGene(VdjSegment::J), "sample_gene.j");
        assert_pinned(
            ChoiceAddress::SampleAlleleInSlot(VdjSegment::V),
            "sample_allele_in_slot.v",
        );
        assert_pinned(
            ChoiceAddress::SampleAlleleInSlot(VdjSegment::D),
            "sample_allele_in_slot.d",
        );
        assert_pinned(
            ChoiceAddress::SampleAlleleInSlot(VdjSegment::J),
            "sample_allele_in_slot.j",
        );
    }
}
