//! Reference-data validator.
//!
//! Rejects malformed refdata before simulation runs, instead of
//! letting bad anchors / empty pools / invalid bytes surface later
//! as confusing contract failures, projection mismatches, or
//! cache-parity divergences.
//!
//! ## Usage
//!
//! ```ignore
//! let cfg: RefDataConfig = load_v6dat("human_igh.v6dat")?;
//! cfg.validate_strict()?;  // fail fast at compile time
//! // ... or, to collect every issue rather than short-circuit ...
//! for issue in cfg.validate() { eprintln!("warn: {issue}"); }
//! ```
//!
//! ## Scope of this slice
//!
//! Read-only: the validator NEVER mutates the config. It runs at
//! O(N * L) over the alleles. It is NOT yet wired into compile —
//! that lives in the next slice. For now the validator is exposed
//! and tested but the compile path is unchanged so existing
//! consumers don't observe a behavior change.
//!
//! ## J anchor convention
//!
//! Expected amino acid at the J anchor is locus-driven, matching
//! the bundled refdata's actual conventions:
//!
//! - IGH → W (Trp)
//! - IGK → F (Phe)
//! - IGL → F (Phe)
//! - TRA / TRB / TRG / TRD → F (TCR convention)
//! - Unknown / no recognised locus prefix → accept either W or F
//!
//! This is deliberately not a biology rule — it mirrors what the
//! bundled `*.v6dat` files actually contain. Tighten per-locus
//! later if needed.

use super::{Allele, AlleleId, AllelePool, ChainType, RefDataConfig};
use crate::codon::translate_codon;
use crate::ir::Segment;
use std::collections::HashSet;

/// Classification of how an issue should gate compilation.
///
/// - [`Fatal`](Self::Fatal): structural problems the engine cannot
///   work around — empty required pools, duplicate allele names,
///   invalid sequence bytes, anchor positions that step past the
///   end of the allele sequence. These reject compile in every mode.
///
/// - [`Curatable`](Self::Curatable): biologically real but
///   functionally non-canonical entries. Real reference catalogues
///   (IMGT, the bundled mouse_igh and human_tcrb data) include
///   pseudogenes and ORF alleles whose anchor codons don't translate
///   to the conserved Cys / W / F residue, or where the anchor is
///   absent. These reject compile under [`RefDataValidationMode::Strict`]
///   but pass under [`RefDataValidationMode::AllowCuratable`]. The
///   `Curatable` label is deliberate: long-term these alleles should
///   be filtered by an explicit curation policy (a future
///   `filter_functional_alleles()` step), not by a blanket "skip
///   validation" escape.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum RefDataIssueSeverity {
    /// Cannot be opted out of; rejects compile in every mode.
    Fatal,
    /// Reflects a pseudogene/ORF / non-canonical allele. Rejects
    /// compile under strict validation; user opts in via
    /// [`RefDataValidationMode::AllowCuratable`] to accept the
    /// catalogue as-is, or filters before compile.
    Curatable,
}

// ──────────────────────────────────────────────────────────────────
// Reference rules — the programmable interpretation layer
// ──────────────────────────────────────────────────────────────────

/// Allowed nucleotide alphabet for an allele sequence. Case-folded;
/// `is_allowed` accepts upper or lower case. Default: `A/C/G/T/N`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ReferenceAlphabet {
    /// Uppercase letters considered valid. Mixed-case input is
    /// folded before comparison so callers don't need to think
    /// about case.
    pub allowed: Vec<u8>,
}

impl ReferenceAlphabet {
    /// Default `A/C/G/T/N`.
    pub fn standard_dna_with_n() -> Self {
        Self {
            allowed: vec![b'A', b'C', b'G', b'T', b'N'],
        }
    }

    /// Returns whether `byte` is in the allowed set (case-insensitive).
    pub fn is_allowed(&self, byte: u8) -> bool {
        let upper = byte.to_ascii_uppercase();
        self.allowed.iter().any(|&a| a.to_ascii_uppercase() == upper)
    }
}

impl Default for ReferenceAlphabet {
    fn default() -> Self {
        Self::standard_dna_with_n()
    }
}

/// Anchor rule for one segment (V or J). Drives both whether an
/// anchor is required and how anchor-codon mismatches are classified.
///
/// - `required = true` + anchor missing → emits `MissingAnchor`
///   tagged with `missing_severity`.
/// - anchor codon AA outside `expected_amino_acids` → emits the
///   appropriate `VAnchorNotCys` / `JAnchorUnexpectedAa` variant
///   tagged with `mismatch_severity`.
///
/// Default severities (`Curatable`) preserve the brief's intent:
/// pseudogene-shape anomalies surface but don't gate strict-mode
/// compile unnecessarily for catalogues that include them.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AnchorRule {
    /// When `true`, alleles without an anchor produce a
    /// `MissingAnchor` issue.
    pub required: bool,
    /// Amino acids the anchor codon is allowed to translate to.
    /// Anchor codons translating to any other AA produce a mismatch
    /// issue. The set itself is locus-specific (V = `['C']` always;
    /// J = `['W']` for IGH, `['F']` for IGK/IGL/TR*); construct
    /// configs with the rule appropriate for the catalogue's locus.
    pub expected_amino_acids: Vec<char>,
    /// Severity carried on `MissingAnchor` issues this rule emits.
    pub missing_severity: RefDataIssueSeverity,
    /// Severity carried on anchor-codon mismatch issues this rule
    /// emits.
    pub mismatch_severity: RefDataIssueSeverity,
}

impl AnchorRule {
    /// Default V rule: Cys-only, anchor required, both severities
    /// Curatable. Real V catalogues with pseudogenes still load
    /// under `AllowCuratable` mode.
    pub fn cys_required_curatable() -> Self {
        Self {
            required: true,
            expected_amino_acids: vec!['C'],
            missing_severity: RefDataIssueSeverity::Curatable,
            mismatch_severity: RefDataIssueSeverity::Curatable,
        }
    }

    /// Default J rule: accepts `W` or `F` (lenient — preserves the
    /// previous "unknown locus" behaviour for synthetic test
    /// fixtures whose allele names don't match an AIRR locus
    /// prefix). Bundled loaders narrow this to the locus-specific
    /// set (`['W']` for IGH, `['F']` for IGK/IGL/TR*).
    pub fn w_or_f_required_curatable() -> Self {
        Self {
            required: true,
            expected_amino_acids: vec!['W', 'F'],
            missing_severity: RefDataIssueSeverity::Curatable,
            mismatch_severity: RefDataIssueSeverity::Curatable,
        }
    }
}

/// The programmable rules layer of a reference cartridge.
///
/// Holds the rules the validator and projection layers consult when
/// interpreting a catalogue. Today this carries the anchor rules and
/// the allowed alphabet; later slices will extend it (junction
/// convention, productivity convention, etc.). Two `RefDataConfig`s
/// with identical catalogues but different `ReferenceRules` validate
/// differently — see [`RefDataConfig::validate_with_mode`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReferenceRules {
    pub alphabet: ReferenceAlphabet,
    pub v_anchor: AnchorRule,
    pub j_anchor: AnchorRule,
}

impl Default for ReferenceRules {
    fn default() -> Self {
        Self {
            alphabet: ReferenceAlphabet::default(),
            v_anchor: AnchorRule::cys_required_curatable(),
            j_anchor: AnchorRule::w_or_f_required_curatable(),
        }
    }
}

impl ReferenceRules {
    /// Return the rules appropriate for an AIRR locus prefix
    /// (`"IGH"`, `"IGK"`, `"IGL"`, `"TRA"`, `"TRB"`, `"TRG"`,
    /// `"TRD"`). Unknown prefixes fall back to `Default` (J accepts
    /// `W` or `F`).
    ///
    /// Used by the Python `dataconfig_to_refdata` loader to stamp
    /// locus-appropriate defaults onto bundled refdata. Synthetic
    /// test refdata that doesn't go through that loader keeps the
    /// lenient default — explicit setters override either path.
    pub fn for_locus(locus_prefix: &str) -> Self {
        let mut rules = Self::default();
        let expected_j: &[char] = match locus_prefix.to_ascii_uppercase().as_str() {
            "IGH" => &['W'],
            "IGK" | "IGL" => &['F'],
            "TRA" | "TRB" | "TRG" | "TRD" => &['F'],
            _ => return rules,
        };
        rules.j_anchor.expected_amino_acids = expected_j.to_vec();
        rules
    }
}

/// How the compile gate treats curatable validation issues.
///
/// Fatal issues are always rejected regardless of mode — they
/// represent structural corruption that the engine cannot work
/// around. The mode toggle only controls whether curatable issues
/// (pseudogene-shaped allele anomalies) gate compile.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum RefDataValidationMode {
    /// Reject every validation issue, including curatable ones.
    /// Default for production compile paths — high-confidence
    /// functional simulations should not consume pseudogenes
    /// silently.
    Strict,
    /// Accept curatable issues; still reject Fatal ones. Use when
    /// the simulation explicitly intends to sample from the raw
    /// catalogue (including pseudogenes/ORFs), with the
    /// understanding that productive contracts may reject more
    /// records at runtime.
    AllowCuratable,
}

/// One validation finding produced by [`RefDataConfig::validate`].
///
/// Each variant carries enough structured context that a downstream
/// consumer (compile error, MCP response, validator UI) can render
/// it precisely without re-parsing a free-form message.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RefDataValidationIssue {
    /// A pool that is required by the chain type is empty. VJ
    /// requires V and J; VDJ requires V, D, and J.
    EmptyRequiredPool { segment: Segment },
    /// Two alleles in the same pool share a name. Allele names are
    /// the user-visible identifier and must round-trip uniquely.
    DuplicateAlleleName { segment: Segment, name: String },
    /// Allele sequence byte is not in the allowed alphabet
    /// (A/C/G/T/N, case-insensitive). Gaps (`.`) and IUPAC codes
    /// other than `N` must be resolved before alleles enter
    /// `RefDataConfig`.
    InvalidAlleleByte {
        segment: Segment,
        allele_id: AlleleId,
        pos: u32,
        byte: u8,
    },
    /// Anchor position would step past the end of the allele
    /// sequence — there's no full codon at `[anchor, anchor+3)`.
    AnchorOutOfBounds {
        segment: Segment,
        allele_id: AlleleId,
        anchor: u16,
        len: u32,
    },
    /// V anchor codon translates to an amino acid outside the
    /// allowed set for the V anchor rule. Severity is set by the
    /// rule (default `Curatable`).
    VAnchorNotCys {
        allele_id: AlleleId,
        codon: [u8; 3],
        aa: char,
        severity: RefDataIssueSeverity,
    },
    /// J anchor codon translates to an amino acid outside the
    /// allowed set for the J anchor rule. Severity is set by the
    /// rule (default `Curatable`).
    JAnchorUnexpectedAa {
        allele_id: AlleleId,
        codon: [u8; 3],
        aa: char,
        expected: Vec<char>,
        severity: RefDataIssueSeverity,
    },
    /// V or J allele has no anchor. Emitted only when the rule for
    /// that segment has `required = true`; severity is set by the
    /// rule's `missing_severity` (default `Curatable`).
    MissingAnchor {
        segment: Segment,
        allele_id: AlleleId,
        severity: RefDataIssueSeverity,
    },
    /// Declared identity locus disagrees with the cartridge's chain
    /// topology. `IGH`/`TRB`/`TRD` require `ChainType::Vdj`;
    /// `IGK`/`IGL`/`TRA`/`TRG` require `ChainType::Vj`. Always
    /// `Fatal` — chain topology drives recombination shape, and a
    /// mismatch would mis-wire the assembly pipeline. Unknown
    /// locus prefixes don't produce this issue.
    LocusChainTypeMismatch {
        locus: String,
        chain_type: ChainType,
    },
}

impl RefDataValidationIssue {
    /// Severity classification — see [`RefDataIssueSeverity`].
    ///
    /// Structural problems (empty pools, duplicate names, invalid
    /// bytes, anchor out of bounds) are always `Fatal` — the engine
    /// cannot work around them. Rule-controlled variants
    /// (`VAnchorNotCys`, `JAnchorUnexpectedAa`, `MissingAnchor`)
    /// carry the severity assigned by the active [`AnchorRule`] at
    /// the moment they were emitted; default is `Curatable`.
    pub fn severity(&self) -> RefDataIssueSeverity {
        use RefDataIssueSeverity::*;
        use RefDataValidationIssue::*;
        match self {
            EmptyRequiredPool { .. }
            | DuplicateAlleleName { .. }
            | InvalidAlleleByte { .. }
            | AnchorOutOfBounds { .. }
            | LocusChainTypeMismatch { .. } => Fatal,
            VAnchorNotCys { severity, .. }
            | JAnchorUnexpectedAa { severity, .. }
            | MissingAnchor { severity, .. } => *severity,
        }
    }
}

impl std::fmt::Display for RefDataValidationIssue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use RefDataValidationIssue::*;
        match self {
            EmptyRequiredPool { segment } => write!(
                f,
                "{segment:?} pool is empty but required for this chain type"
            ),
            DuplicateAlleleName { segment, name } => {
                write!(f, "duplicate allele name '{name}' in {segment:?} pool")
            }
            InvalidAlleleByte {
                segment,
                allele_id,
                pos,
                byte,
            } => write!(
                f,
                "{segment:?} allele {} byte at position {pos} is 0x{byte:02x} ('{}'); allowed: A/C/G/T/N (case-insensitive)",
                allele_id.index(),
                *byte as char,
            ),
            AnchorOutOfBounds {
                segment,
                allele_id,
                anchor,
                len,
            } => write!(
                f,
                "{segment:?} allele {} anchor at {anchor} but full codon needs <= {} (len={len})",
                allele_id.index(),
                len.saturating_sub(3),
            ),
            VAnchorNotCys {
                allele_id,
                codon,
                aa,
                ..
            } => write!(
                f,
                "V allele {} anchor codon '{}' translates to '{aa}' (expected C)",
                allele_id.index(),
                String::from_utf8_lossy(codon),
            ),
            JAnchorUnexpectedAa {
                allele_id,
                codon,
                aa,
                expected,
                ..
            } => {
                let expected_str: String = expected.iter().collect();
                write!(
                    f,
                    "J allele {} anchor codon '{}' translates to '{aa}' (expected one of {expected_str:?})",
                    allele_id.index(),
                    String::from_utf8_lossy(codon),
                )
            }
            MissingAnchor { segment, allele_id, .. } => write!(
                f,
                "{segment:?} allele {} has no anchor",
                allele_id.index()
            ),
            LocusChainTypeMismatch { locus, chain_type } => write!(
                f,
                "identity locus '{locus}' is incompatible with chain_type \
                 {chain_type:?}; IGH/TRB/TRD require Vdj, IGK/IGL/TRA/TRG \
                 require Vj",
            ),
        }
    }
}

/// Aggregated error result with structured issue list. Returned by
/// [`RefDataConfig::validate_strict`] when one or more issues exist.
///
/// Implements [`std::error::Error`] so it slots into the
/// compile path's error pipeline once that wiring lands.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RefDataValidationErrors {
    pub issues: Vec<RefDataValidationIssue>,
    /// Mode under which validation was run. Determines whether the
    /// trailing remediation hint mentions `allow_curatable_refdata`
    /// (Strict — some issues are curatable) or just the catalogue
    /// fix-up path (AllowCuratable — only Fatal issues remained).
    pub mode: RefDataValidationMode,
}

impl RefDataValidationErrors {
    /// Number of issues by severity (`(fatal_count, curatable_count)`).
    pub fn severity_counts(&self) -> (usize, usize) {
        let mut fatal = 0;
        let mut curatable = 0;
        for i in &self.issues {
            match i.severity() {
                RefDataIssueSeverity::Fatal => fatal += 1,
                RefDataIssueSeverity::Curatable => curatable += 1,
            }
        }
        (fatal, curatable)
    }
}

impl std::fmt::Display for RefDataValidationErrors {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (fatal, curatable) = self.severity_counts();
        writeln!(
            f,
            "{} refdata validation issue(s) ({} fatal, {} curatable):",
            self.issues.len(),
            fatal,
            curatable,
        )?;
        for issue in &self.issues {
            let tag = match issue.severity() {
                RefDataIssueSeverity::Fatal => "fatal",
                RefDataIssueSeverity::Curatable => "curatable",
            };
            writeln!(f, "  - [{tag}] {issue}")?;
        }
        // Remediation hint: if any curatable issues were the only
        // blockers (Strict mode with no Fatal), the user can opt in
        // via `allow_curatable_refdata` or filter their catalogue.
        // If Fatal issues exist, neither opt-in helps and the
        // catalogue must be fixed.
        if curatable > 0 && fatal == 0 {
            writeln!(
                f,
                "These may represent pseudogene/ORF alleles. Use \
                 allow_curatable_refdata() or filter reference alleles \
                 before simulation."
            )?;
        }
        Ok(())
    }
}

impl std::error::Error for RefDataValidationErrors {}

impl RefDataConfig {
    /// Return the full list of refdata validation issues. Empty
    /// list = passes every gate. Read-only; never panics.
    pub fn validate(&self) -> Vec<RefDataValidationIssue> {
        let mut issues = Vec::new();

        // Required pools per chain type.
        if self.v_pool.is_empty() {
            issues.push(RefDataValidationIssue::EmptyRequiredPool {
                segment: Segment::V,
            });
        }
        if self.j_pool.is_empty() {
            issues.push(RefDataValidationIssue::EmptyRequiredPool {
                segment: Segment::J,
            });
        }
        if self.chain_type.has_d() && self.d_pool.is_empty() {
            issues.push(RefDataValidationIssue::EmptyRequiredPool {
                segment: Segment::D,
            });
        }

        // Identity ↔ chain_type consistency. If the cartridge
        // declares a locus, it must match the chain topology: IGH /
        // TRB / TRD are VDJ; IGK / IGL / TRA / TRG are VJ. Unknown
        // loci are silent — the engine doesn't enforce expectations
        // about catalogues whose origin it can't recognise.
        if let Some(locus) = self.identity.locus.as_deref() {
            let upper = locus.to_ascii_uppercase();
            let expected_has_d = match upper.as_str() {
                "IGH" | "TRB" | "TRD" => Some(true),
                "IGK" | "IGL" | "TRA" | "TRG" => Some(false),
                _ => None,
            };
            if let Some(expected) = expected_has_d {
                if expected != self.chain_type.has_d() {
                    issues.push(RefDataValidationIssue::LocusChainTypeMismatch {
                        locus: upper,
                        chain_type: self.chain_type,
                    });
                }
            }
        }

        validate_pool(&self.v_pool, Segment::V, &self.rules, &mut issues);
        validate_pool(&self.d_pool, Segment::D, &self.rules, &mut issues);
        validate_pool(&self.j_pool, Segment::J, &self.rules, &mut issues);

        issues
    }

    /// Strict mode: returns `Err` listing every issue, or `Ok(())`
    /// when validation passes. Equivalent to
    /// `validate_with_mode(RefDataValidationMode::Strict)`. Designed
    /// for compile-time gating; the compile path can convert this
    /// into a structured `CompileError` once wired in.
    pub fn validate_strict(&self) -> Result<(), RefDataValidationErrors> {
        self.validate_with_mode(RefDataValidationMode::Strict)
    }

    /// Mode-aware validation gate.
    ///
    /// Under [`RefDataValidationMode::Strict`], any issue rejects.
    /// Under [`RefDataValidationMode::AllowCuratable`], issues
    /// classified as [`RefDataIssueSeverity::Curatable`] are filtered
    /// out — Fatal issues still reject. The returned error always
    /// preserves the **full** issue list (no surprise dropping); the
    /// mode only controls whether the result is `Ok` or `Err`.
    pub fn validate_with_mode(
        &self,
        mode: RefDataValidationMode,
    ) -> Result<(), RefDataValidationErrors> {
        let issues = self.validate();
        if issues.is_empty() {
            return Ok(());
        }
        let any_blocking = issues.iter().any(|i| match mode {
            RefDataValidationMode::Strict => true,
            RefDataValidationMode::AllowCuratable => {
                i.severity() == RefDataIssueSeverity::Fatal
            }
        });
        if any_blocking {
            Err(RefDataValidationErrors { issues, mode })
        } else {
            Ok(())
        }
    }
}

fn validate_pool(
    pool: &AllelePool,
    segment: Segment,
    rules: &ReferenceRules,
    issues: &mut Vec<RefDataValidationIssue>,
) {
    let mut seen_names: HashSet<&str> = HashSet::with_capacity(pool.len());
    for (id, allele) in pool.iter() {
        // Identity: duplicate names.
        if !seen_names.insert(allele.name.as_str()) {
            issues.push(RefDataValidationIssue::DuplicateAlleleName {
                segment,
                name: allele.name.clone(),
            });
        }
        // Sequence byte alphabet — driven by `rules.alphabet`.
        for (pos, &byte) in allele.seq.iter().enumerate() {
            if !rules.alphabet.is_allowed(byte) {
                issues.push(RefDataValidationIssue::InvalidAlleleByte {
                    segment,
                    allele_id: id,
                    pos: pos as u32,
                    byte,
                });
            }
        }
        // Anchor checks only apply to V/J — D segments are typically
        // anchorless in real reference data, and the validator
        // shouldn't invent expectations the bundled data doesn't meet.
        if matches!(segment, Segment::V | Segment::J) {
            validate_anchor(allele, id, segment, rules, issues);
        }
    }
}

fn validate_anchor(
    allele: &Allele,
    id: AlleleId,
    segment: Segment,
    rules: &ReferenceRules,
    issues: &mut Vec<RefDataValidationIssue>,
) {
    let rule = match segment {
        Segment::V => &rules.v_anchor,
        Segment::J => &rules.j_anchor,
        _ => unreachable!("validate_anchor only called for V/J"),
    };
    let Some(anchor) = allele.anchor else {
        if rule.required {
            issues.push(RefDataValidationIssue::MissingAnchor {
                segment,
                allele_id: id,
                severity: rule.missing_severity,
            });
        }
        return;
    };
    let len = allele.seq.len() as u32;
    let anchor_u32 = anchor as u32;
    if anchor_u32 + 3 > len {
        // AnchorOutOfBounds is always structural Fatal — the engine
        // cannot index into a codon that isn't there. No rule
        // controls this.
        issues.push(RefDataValidationIssue::AnchorOutOfBounds {
            segment,
            allele_id: id,
            anchor,
            len,
        });
        return;
    }
    let a = anchor as usize;
    let codon = [allele.seq[a], allele.seq[a + 1], allele.seq[a + 2]];
    let aa = translate_codon(codon[0], codon[1], codon[2]) as char;
    if rule.expected_amino_acids.contains(&aa) {
        return;
    }
    match segment {
        Segment::V => issues.push(RefDataValidationIssue::VAnchorNotCys {
            allele_id: id,
            codon,
            aa,
            severity: rule.mismatch_severity,
        }),
        Segment::J => issues.push(RefDataValidationIssue::JAnchorUnexpectedAa {
            allele_id: id,
            codon,
            aa,
            expected: rule.expected_amino_acids.clone(),
            severity: rule.mismatch_severity,
        }),
        _ => unreachable!(),
    }
}

#[allow(dead_code)]
fn _locus_prefix(name: &str) -> String {
    name.chars().take(3).flat_map(|c| c.to_uppercase()).collect()
}

/// Discourage `ChainType` being silently swapped — bring it into
/// scope here so the validator module compiles cleanly even when
/// the only consumer of the import is a doctest example.
#[allow(dead_code)]
fn _chain_type_imported(ct: ChainType) -> bool {
    ct.has_d()
}

// ──────────────────────────────────────────────────────────────────
// tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::refdata::{Allele, AlleleId, ChainType, RefDataConfig};

    fn v_allele(name: &str, seq: &[u8], anchor: Option<u16>) -> Allele {
        Allele {
            name: name.to_string(),
            gene: name.split('*').next().unwrap_or(name).to_string(),
            seq: seq.to_vec(),
            segment: Segment::V,
            anchor,
            functional_status: None,
            subregions: Vec::new(),
        }
    }
    fn d_allele(name: &str, seq: &[u8]) -> Allele {
        Allele {
            name: name.to_string(),
            gene: name.split('*').next().unwrap_or(name).to_string(),
            seq: seq.to_vec(),
            segment: Segment::D,
            anchor: None,
            functional_status: None,
            subregions: Vec::new(),
        }
    }
    fn j_allele(name: &str, seq: &[u8], anchor: Option<u16>) -> Allele {
        Allele {
            name: name.to_string(),
            gene: name.split('*').next().unwrap_or(name).to_string(),
            seq: seq.to_vec(),
            segment: Segment::J,
            anchor,
            functional_status: None,
            subregions: Vec::new(),
        }
    }

    /// Build a minimal valid VDJ refdata so each test only varies
    /// the one field under examination.
    fn minimal_vdj() -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        // V allele: anchor codon TGT (Cys) at position 0.
        let _ = cfg.v_pool.push(v_allele("IGHV1-1*01", b"TGTAAACCC", Some(0)));
        // D allele: no anchor.
        let _ = cfg.d_pool.push(d_allele("IGHD1-1*01", b"GGGCCCAAA"));
        // J allele: anchor codon TGG (Trp = W) at position 0.
        let _ = cfg.j_pool.push(j_allele("IGHJ1*01", b"TGGAAACCC", Some(0)));
        cfg
    }

    fn minimal_vj() -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(v_allele("IGKV1-1*01", b"TGTAAACCC", Some(0)));
        // IGK J anchor codon TTC (Phe = F).
        let _ = cfg.j_pool.push(j_allele("IGKJ1*01", b"TTCAAACCC", Some(0)));
        cfg
    }

    #[test]
    fn minimal_vdj_validates_clean() {
        assert!(minimal_vdj().validate().is_empty());
        assert!(minimal_vdj().validate_strict().is_ok());
    }

    #[test]
    fn minimal_vj_validates_clean() {
        assert!(minimal_vj().validate().is_empty());
    }

    #[test]
    fn empty_v_pool_fails_vj() {
        let mut cfg = minimal_vj();
        cfg.v_pool = super::AllelePool::new();
        let issues = cfg.validate();
        assert!(issues.contains(&RefDataValidationIssue::EmptyRequiredPool {
            segment: Segment::V,
        }));
    }

    #[test]
    fn empty_v_pool_fails_vdj() {
        let mut cfg = minimal_vdj();
        cfg.v_pool = super::AllelePool::new();
        let issues = cfg.validate();
        assert!(issues.contains(&RefDataValidationIssue::EmptyRequiredPool {
            segment: Segment::V,
        }));
    }

    #[test]
    fn empty_j_pool_fails_vj() {
        let mut cfg = minimal_vj();
        cfg.j_pool = super::AllelePool::new();
        let issues = cfg.validate();
        assert!(issues.contains(&RefDataValidationIssue::EmptyRequiredPool {
            segment: Segment::J,
        }));
    }

    #[test]
    fn empty_d_pool_allowed_vj() {
        // VJ chains have empty D pools by construction. This must
        // not flag an EmptyRequiredPool for D.
        let cfg = minimal_vj();
        let issues = cfg.validate();
        assert!(!issues
            .iter()
            .any(|i| matches!(i, RefDataValidationIssue::EmptyRequiredPool { segment } if *segment == Segment::D)));
    }

    #[test]
    fn empty_d_pool_rejected_vdj() {
        let mut cfg = minimal_vdj();
        cfg.d_pool = super::AllelePool::new();
        let issues = cfg.validate();
        assert!(issues.contains(&RefDataValidationIssue::EmptyRequiredPool {
            segment: Segment::D,
        }));
    }

    #[test]
    fn duplicate_v_allele_names_fail() {
        let mut cfg = minimal_vdj();
        let _ = cfg.v_pool.push(v_allele("IGHV1-1*01", b"TGTAAA", Some(0)));
        let issues = cfg.validate();
        assert!(issues.contains(&RefDataValidationIssue::DuplicateAlleleName {
            segment: Segment::V,
            name: "IGHV1-1*01".to_string(),
        }));
    }

    #[test]
    fn duplicate_j_allele_names_fail() {
        let mut cfg = minimal_vdj();
        let _ = cfg.j_pool.push(j_allele("IGHJ1*01", b"TGGAAA", Some(0)));
        let issues = cfg.validate();
        assert!(issues.contains(&RefDataValidationIssue::DuplicateAlleleName {
            segment: Segment::J,
            name: "IGHJ1*01".to_string(),
        }));
    }

    #[test]
    fn invalid_byte_in_v_seq_fails_with_exact_position() {
        let mut cfg = minimal_vdj();
        let _ = cfg.v_pool.push(v_allele("IGHV2*01", b"TGT.AAA", Some(0)));
        let issues = cfg.validate();
        let bad_v_id = AlleleId::new(1);
        let expected = RefDataValidationIssue::InvalidAlleleByte {
            segment: Segment::V,
            allele_id: bad_v_id,
            pos: 3,
            byte: b'.',
        };
        assert!(
            issues.contains(&expected),
            "expected invalid-byte issue at pos 3, got: {issues:?}"
        );
    }

    #[test]
    fn lowercase_bases_accepted() {
        let mut cfg = minimal_vdj();
        let _ = cfg.v_pool.push(v_allele("IGHV3*01", b"tgtaaaccc", Some(0)));
        // Lowercase A/C/G/T/N is allowed. The anchor codon `tgt` still
        // translates to Cys, so V anchor check also passes.
        let issues = cfg.validate();
        assert!(
            issues.is_empty(),
            "lowercase canonical bases should validate, got: {issues:?}"
        );
    }

    #[test]
    fn iupac_ambiguity_other_than_n_fails() {
        let mut cfg = minimal_vdj();
        let _ = cfg.v_pool.push(v_allele("IGHV4*01", b"TGTRAACCC", Some(0)));
        let issues = cfg.validate();
        let v_id = AlleleId::new(1);
        assert!(issues.contains(&RefDataValidationIssue::InvalidAlleleByte {
            segment: Segment::V,
            allele_id: v_id,
            pos: 3,
            byte: b'R',
        }));
    }

    #[test]
    fn v_anchor_out_of_bounds_fails() {
        let mut cfg = minimal_vdj();
        let _ = cfg.v_pool.push(v_allele("IGHV5*01", b"TGT", Some(5)));
        let issues = cfg.validate();
        let v_id = AlleleId::new(1);
        assert!(issues.contains(&RefDataValidationIssue::AnchorOutOfBounds {
            segment: Segment::V,
            allele_id: v_id,
            anchor: 5,
            len: 3,
        }));
    }

    #[test]
    fn v_anchor_at_boundary_with_no_full_codon_fails() {
        // seq.len() == 5, anchor == 3 → codon would span [3, 6) but
        // len is 5 → out of bounds.
        let mut cfg = minimal_vdj();
        let _ = cfg.v_pool.push(v_allele("IGHV6*01", b"AAATG", Some(3)));
        let issues = cfg.validate();
        let v_id = AlleleId::new(1);
        assert!(issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::AnchorOutOfBounds { segment, allele_id, .. }
                if *segment == Segment::V && *allele_id == v_id
        )));
    }

    #[test]
    fn v_anchor_non_cys_fails() {
        let mut cfg = minimal_vdj();
        // GGG = Gly, not Cys.
        let _ = cfg.v_pool.push(v_allele("IGHV7*01", b"GGGAAACCC", Some(0)));
        let issues = cfg.validate();
        let v_id = AlleleId::new(1);
        assert!(issues.contains(&RefDataValidationIssue::VAnchorNotCys {
            allele_id: v_id,
            codon: [b'G', b'G', b'G'],
            aa: 'G',
            severity: RefDataIssueSeverity::Curatable,
        }));
    }

    #[test]
    fn v_anchor_tgc_is_cys_and_passes() {
        // TGC is also Cys (degenerate; both TGT and TGC code for C).
        let mut cfg = minimal_vdj();
        let _ = cfg.v_pool.push(v_allele("IGHV8*01", b"TGCAAACCC", Some(0)));
        let issues = cfg.validate();
        assert!(
            issues.is_empty(),
            "TGC anchor codon should pass V Cys check, got: {issues:?}"
        );
    }

    #[test]
    fn j_anchor_unexpected_aa_under_igh_rule_fails() {
        // The catalogue claims IGH (rules.j_anchor.expected = ['W']);
        // a J allele with anchor codon TTC (Phe = F) is flagged.
        // Allele names are no longer used as the source of truth for
        // expected J anchors — the rule on the config is.
        let mut cfg = minimal_vdj();
        cfg.rules = ReferenceRules::for_locus("IGH");
        let _ = cfg.j_pool.push(j_allele("IGHJ2*01", b"TTCAAACCC", Some(0)));
        let issues = cfg.validate();
        let j_id = AlleleId::new(1);
        assert!(issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::JAnchorUnexpectedAa { allele_id, aa, .. }
                if *allele_id == j_id && *aa == 'F'
        )));
    }

    #[test]
    fn j_anchor_accepts_f_under_igl_rule() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules = ReferenceRules::for_locus("IGL");
        let _ = cfg.v_pool.push(v_allele("IGLV1*01", b"TGTAAA", Some(0)));
        let _ = cfg.j_pool.push(j_allele("IGLJ1*01", b"TTCAAACCC", Some(0)));
        assert!(cfg.validate().is_empty());
    }

    #[test]
    fn j_anchor_w_under_igl_rule_is_flagged() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules = ReferenceRules::for_locus("IGL");
        let _ = cfg.v_pool.push(v_allele("IGLV1*01", b"TGTAAA", Some(0)));
        // IGL rule expects F — a W-coded J anchor is flagged.
        let _ = cfg.j_pool.push(j_allele("IGLJ2*01", b"TGGAAACCC", Some(0)));
        let issues = cfg.validate();
        assert!(issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::JAnchorUnexpectedAa { aa, .. } if *aa == 'W'
        )));
    }

    #[test]
    fn default_rule_accepts_w_or_f_for_unconfigured_locus() {
        // `RefDataConfig::empty` initialises with the lenient default
        // (J expects W or F). Any allele-name pattern that the loader
        // doesn't recognise lands here.
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(v_allele("MYV1*01", b"TGTAAA", Some(0)));
        let _ = cfg.j_pool.push(j_allele("MYJ1*01", b"TGGAAACCC", Some(0)));
        assert!(
            cfg.validate().is_empty(),
            "default-locus J anchor 'W' should be accepted"
        );

        let mut cfg2 = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg2.v_pool.push(v_allele("MYV1*01", b"TGTAAA", Some(0)));
        let _ = cfg2.j_pool.push(j_allele("MYJ2*01", b"TTCAAACCC", Some(0)));
        assert!(
            cfg2.validate().is_empty(),
            "default-locus J anchor 'F' should be accepted"
        );
    }

    #[test]
    fn missing_anchor_on_v_is_reported() {
        let mut cfg = minimal_vdj();
        let _ = cfg.v_pool.push(v_allele("IGHV9*01", b"AAACCCGGG", None));
        let issues = cfg.validate();
        let v_id = AlleleId::new(1);
        assert!(issues.contains(&RefDataValidationIssue::MissingAnchor {
            segment: Segment::V,
            allele_id: v_id,
            severity: RefDataIssueSeverity::Curatable,
        }));
    }

    #[test]
    fn missing_anchor_on_j_is_reported() {
        let mut cfg = minimal_vdj();
        let _ = cfg.j_pool.push(j_allele("IGHJ-orphan*01", b"GGGAAA", None));
        let issues = cfg.validate();
        let j_id = AlleleId::new(1);
        assert!(issues.contains(&RefDataValidationIssue::MissingAnchor {
            segment: Segment::J,
            allele_id: j_id,
            severity: RefDataIssueSeverity::Curatable,
        }));
    }

    #[test]
    fn d_anchor_absence_is_not_reported() {
        // The default minimal_vdj D allele has no anchor; this must
        // not flag a MissingAnchor issue (D anchors are not part of
        // the contract).
        let cfg = minimal_vdj();
        let issues = cfg.validate();
        assert!(!issues
            .iter()
            .any(|i| matches!(i, RefDataValidationIssue::MissingAnchor { segment, .. } if *segment == Segment::D)));
    }

    #[test]
    fn validate_strict_packages_errors() {
        let mut cfg = minimal_vdj();
        let _ = cfg.v_pool.push(v_allele("IGHV1-1*01", b"TGTAAA", Some(0)));
        let err = cfg.validate_strict().expect_err("duplicate must fail strict");
        assert_eq!(err.issues.len(), 1);
        assert!(format!("{err}").contains("duplicate allele name"));
    }

    #[test]
    fn validate_strict_collects_all_issues_not_just_first() {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        // Multiple problems at once: empty V, empty J, empty D.
        let err = cfg.validate_strict().expect_err("all empty should fail");
        let segs: Vec<_> = err
            .issues
            .iter()
            .filter_map(|i| match i {
                RefDataValidationIssue::EmptyRequiredPool { segment } => Some(*segment),
                _ => None,
            })
            .collect();
        assert!(segs.contains(&Segment::V));
        assert!(segs.contains(&Segment::J));
        assert!(segs.contains(&Segment::D));

        // Plug V, leave D and J empty: still has two issues.
        let _ = cfg.v_pool.push(v_allele("IGHV1*01", b"TGTAAA", Some(0)));
        let err2 = cfg.validate_strict().expect_err("two still missing");
        assert_eq!(err2.issues.len(), 2);
    }

    // ── Severity classification + mode-aware gating ───────────────

    /// VJ refdata with only curatable issues (Gly V anchor + no J
    /// anchor). Useful to verify mode toggling without bleed-over
    /// from Fatal issues.
    fn vj_curatable_only() -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(v_allele("IGKV1-1*01", b"AAACCCGGG", Some(6)));
        let _ = cfg.j_pool.push(Allele {
            name: "IGKJ-orphan*01".into(),
            gene: "IGKJ-orphan".into(),
            seq: b"TTCAAACCC".to_vec(),
            segment: Segment::J,
            anchor: None,
            functional_status: None,
            subregions: Vec::new(),
        });
        cfg
    }

    #[test]
    fn validate_with_strict_mode_rejects_curatable_issues() {
        let cfg = vj_curatable_only();
        let err = cfg
            .validate_with_mode(RefDataValidationMode::Strict)
            .expect_err("strict must reject curatable issues");
        let (fatal, curatable) = err.severity_counts();
        assert_eq!(fatal, 0);
        assert!(curatable >= 2);
        assert_eq!(err.mode, RefDataValidationMode::Strict);
    }

    #[test]
    fn validate_with_allow_curatable_accepts_curatable_only() {
        let cfg = vj_curatable_only();
        cfg.validate_with_mode(RefDataValidationMode::AllowCuratable)
            .expect("allow_curatable must accept curatable-only fixtures");
    }

    #[test]
    fn validate_with_allow_curatable_still_rejects_fatal() {
        let mut cfg = vj_curatable_only();
        // Add an invalid byte allele (Fatal).
        let _ = cfg.v_pool.push(v_allele("IGKV2*01", b"TGT.AAACC", Some(0)));
        let err = cfg
            .validate_with_mode(RefDataValidationMode::AllowCuratable)
            .expect_err("allow_curatable must still reject Fatal issues");
        let (fatal, curatable) = err.severity_counts();
        assert!(fatal >= 1, "expected ≥1 fatal, got: {:?}", err.issues);
        assert!(curatable >= 1, "expected ≥1 curatable preserved, got: {:?}", err.issues);
        assert_eq!(err.mode, RefDataValidationMode::AllowCuratable);
    }

    #[test]
    fn display_remediation_hint_when_only_curatable_issues_remain() {
        let cfg = vj_curatable_only();
        let err = cfg.validate_strict().expect_err("must fail");
        let msg = format!("{err}");
        assert!(
            msg.contains("allow_curatable_refdata"),
            "msg must suggest allow_curatable_refdata when only curatable: {msg}"
        );
        assert!(
            msg.contains("pseudogene/ORF"),
            "msg must mention pseudogene/ORF: {msg}"
        );
    }

    #[test]
    fn display_omits_remediation_hint_when_fatal_present() {
        let mut cfg = vj_curatable_only();
        let _ = cfg.v_pool.push(v_allele("IGKV2*01", b"TGT.AAACC", Some(0)));
        let err = cfg.validate_strict().expect_err("must fail");
        let msg = format!("{err}");
        // Fatal issues exist — opting into allow_curatable_refdata
        // wouldn't help, so the hint must NOT appear.
        assert!(
            !msg.contains("allow_curatable_refdata"),
            "msg must NOT suggest allow_curatable_refdata when Fatal issues remain: {msg}"
        );
    }

    #[test]
    fn each_issue_variant_carries_documented_severity() {
        use RefDataIssueSeverity::*;
        let aid = AlleleId::new(0);
        // Fatal set.
        assert_eq!(
            RefDataValidationIssue::EmptyRequiredPool { segment: Segment::V }.severity(),
            Fatal
        );
        assert_eq!(
            RefDataValidationIssue::DuplicateAlleleName {
                segment: Segment::J,
                name: "x".into()
            }
            .severity(),
            Fatal
        );
        assert_eq!(
            RefDataValidationIssue::InvalidAlleleByte {
                segment: Segment::V,
                allele_id: aid,
                pos: 0,
                byte: b'.'
            }
            .severity(),
            Fatal
        );
        assert_eq!(
            RefDataValidationIssue::AnchorOutOfBounds {
                segment: Segment::V,
                allele_id: aid,
                anchor: 99,
                len: 5
            }
            .severity(),
            Fatal
        );
        // Rule-controlled set — severity is whatever the construction
        // site stamped on the issue (`severity()` just returns it).
        assert_eq!(
            RefDataValidationIssue::VAnchorNotCys {
                allele_id: aid,
                codon: [b'G', b'G', b'G'],
                aa: 'G',
                severity: Curatable,
            }
            .severity(),
            Curatable
        );
        assert_eq!(
            RefDataValidationIssue::JAnchorUnexpectedAa {
                allele_id: aid,
                codon: [b'T', b'T', b'A'],
                aa: 'L',
                expected: vec!['W'],
                severity: Curatable,
            }
            .severity(),
            Curatable
        );
        assert_eq!(
            RefDataValidationIssue::MissingAnchor {
                segment: Segment::V,
                allele_id: aid,
                severity: Curatable,
            }
            .severity(),
            Curatable
        );
        // The rule-controlled variants honour the severity carried
        // on the issue itself — a Fatal anchor mismatch (custom rule)
        // reports as Fatal.
        assert_eq!(
            RefDataValidationIssue::MissingAnchor {
                segment: Segment::V,
                allele_id: aid,
                severity: Fatal,
            }
            .severity(),
            Fatal
        );
    }

    // ── ReferenceRules v1 — configurable anchor + alphabet ──────

    #[test]
    fn default_rules_match_legacy_behavior() {
        // V → ['C'], J → ['W', 'F'], alphabet ACGTN, all Curatable.
        let rules = ReferenceRules::default();
        assert_eq!(rules.v_anchor.expected_amino_acids, vec!['C']);
        assert_eq!(rules.j_anchor.expected_amino_acids, vec!['W', 'F']);
        assert!(rules.v_anchor.required);
        assert!(rules.j_anchor.required);
        for &b in b"acgtnACGTN" {
            assert!(rules.alphabet.is_allowed(b), "should allow {}", b as char);
        }
        assert!(!rules.alphabet.is_allowed(b'.'));
        assert!(!rules.alphabet.is_allowed(b'R'));
    }

    #[test]
    fn for_locus_returns_locus_appropriate_j_rule() {
        assert_eq!(
            ReferenceRules::for_locus("IGH").j_anchor.expected_amino_acids,
            vec!['W']
        );
        assert_eq!(
            ReferenceRules::for_locus("IGK").j_anchor.expected_amino_acids,
            vec!['F']
        );
        assert_eq!(
            ReferenceRules::for_locus("IGL").j_anchor.expected_amino_acids,
            vec!['F']
        );
        for tr in ["TRA", "TRB", "TRG", "TRD"] {
            assert_eq!(
                ReferenceRules::for_locus(tr).j_anchor.expected_amino_acids,
                vec!['F'],
                "{tr} J should expect F"
            );
        }
        // Unknown prefix → default lenient set.
        assert_eq!(
            ReferenceRules::for_locus("XXX").j_anchor.expected_amino_acids,
            vec!['W', 'F']
        );
    }

    #[test]
    fn custom_j_rule_y_accepts_tat_and_tac() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules.j_anchor.expected_amino_acids = vec!['Y'];
        let _ = cfg.v_pool.push(v_allele("MYV1*01", b"TGTAAA", Some(0)));
        // TAT → Y.
        let _ = cfg.j_pool.push(j_allele("MYJ1*01", b"TATAAACCC", Some(0)));
        assert!(cfg.validate().is_empty());
        // TAC → Y.
        let mut cfg2 = cfg.clone();
        cfg2.j_pool = crate::refdata::AllelePool::new();
        let _ = cfg2.j_pool.push(j_allele("MYJ2*01", b"TACAAACCC", Some(0)));
        assert!(cfg2.validate().is_empty());
    }

    #[test]
    fn custom_j_rule_y_rejects_tgg() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules.j_anchor.expected_amino_acids = vec!['Y'];
        let _ = cfg.v_pool.push(v_allele("MYV1*01", b"TGTAAA", Some(0)));
        // TGG → W. Under a Y-only J rule, this is flagged.
        let _ = cfg.j_pool.push(j_allele("MYJ-bad*01", b"TGGAAACCC", Some(0)));
        let issues = cfg.validate();
        assert!(issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::JAnchorUnexpectedAa { aa, .. } if *aa == 'W'
        )));
    }

    #[test]
    fn anchor_required_false_suppresses_missing_anchor_issue() {
        // Build a config where the J rule says anchor is optional.
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules.j_anchor.required = false;
        let _ = cfg.v_pool.push(v_allele("MYV1*01", b"TGTAAA", Some(0)));
        // J allele with no anchor → would normally emit MissingAnchor.
        let _ = cfg.j_pool.push(j_allele("MYJ-orphan*01", b"GGG", None));
        assert!(
            cfg.validate().is_empty(),
            "rule.required=false must suppress MissingAnchor"
        );
    }

    #[test]
    fn missing_severity_follows_rule() {
        // Bump missing_severity to Fatal — the resulting issue's
        // severity() reflects that, and AllowCuratable mode still
        // rejects.
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules.j_anchor.missing_severity = RefDataIssueSeverity::Fatal;
        let _ = cfg.v_pool.push(v_allele("MYV1*01", b"TGTAAA", Some(0)));
        let _ = cfg.j_pool.push(j_allele("MYJ*01", b"TGGAAA", None));
        let issues = cfg.validate();
        let missing = issues
            .iter()
            .find(|i| matches!(i, RefDataValidationIssue::MissingAnchor { .. }))
            .expect("must emit MissingAnchor");
        assert_eq!(missing.severity(), RefDataIssueSeverity::Fatal);

        // AllowCuratable still rejects.
        cfg.validate_with_mode(RefDataValidationMode::AllowCuratable)
            .expect_err("Fatal missing-anchor must not be opt-outable");
    }

    #[test]
    fn mismatch_severity_follows_rule() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules.v_anchor.mismatch_severity = RefDataIssueSeverity::Fatal;
        // GGG (Gly) at V anchor.
        let _ = cfg.v_pool.push(v_allele("MYV-bad*01", b"GGGAAACCC", Some(0)));
        let _ = cfg.j_pool.push(j_allele("MYJ*01", b"TGGAAA", Some(0)));
        let issues = cfg.validate();
        let mismatch = issues
            .iter()
            .find(|i| matches!(i, RefDataValidationIssue::VAnchorNotCys { .. }))
            .expect("must emit VAnchorNotCys");
        assert_eq!(mismatch.severity(), RefDataIssueSeverity::Fatal);

        // AllowCuratable still rejects this catalogue.
        cfg.validate_with_mode(RefDataValidationMode::AllowCuratable)
            .expect_err("Fatal V anchor mismatch must not be opt-outable");
    }

    #[test]
    fn invalid_byte_remains_fatal_regardless_of_rule() {
        // Even with an extremely permissive anchor rule, an invalid
        // sequence byte stays Fatal — structural problems aren't
        // rule-controlled.
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules.v_anchor.expected_amino_acids = vec!['A', 'C', 'D', 'E', 'F'];
        cfg.rules.j_anchor.expected_amino_acids = vec!['A', 'C', 'D', 'E', 'F'];
        let _ = cfg.v_pool.push(v_allele("MYV-bad*01", b"TGT.AAACC", Some(0)));
        let _ = cfg.j_pool.push(j_allele("MYJ*01", b"TGGAAA", Some(0)));
        let issues = cfg.validate();
        let bad = issues
            .iter()
            .find(|i| matches!(i, RefDataValidationIssue::InvalidAlleleByte { .. }))
            .expect("must emit InvalidAlleleByte");
        assert_eq!(bad.severity(), RefDataIssueSeverity::Fatal);
    }

    #[test]
    fn custom_alphabet_extends_allowed_set() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules.alphabet = ReferenceAlphabet {
            allowed: vec![b'A', b'C', b'G', b'T', b'N', b'R'], // R = puRine ambig
        };
        let _ = cfg.v_pool.push(v_allele("MYV*01", b"TGTRAA", Some(0)));
        let _ = cfg.j_pool.push(j_allele("MYJ*01", b"TGGAAA", Some(0)));
        let issues = cfg.validate();
        // R is now in-alphabet → no InvalidAlleleByte issue.
        assert!(
            !issues
                .iter()
                .any(|i| matches!(i, RefDataValidationIssue::InvalidAlleleByte { .. })),
            "R should be allowed under the extended alphabet; got: {issues:?}"
        );
    }

    #[test]
    fn reference_alphabet_is_case_insensitive() {
        let alphabet = ReferenceAlphabet::default();
        assert!(alphabet.is_allowed(b'a'));
        assert!(alphabet.is_allowed(b'A'));
        assert!(alphabet.is_allowed(b'n'));
        assert!(!alphabet.is_allowed(b'.'));
    }

    // ── Identity ↔ chain-type mismatch (Reference Identity slice) ──

    fn vdj_cfg_with_identity_locus(locus: &str) -> RefDataConfig {
        let mut cfg = minimal_vdj();
        cfg.identity.locus = Some(locus.to_string());
        cfg
    }

    fn vj_cfg_with_identity_locus(locus: &str) -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(v_allele("IGKV1-1*01", b"TGTAAACCC", Some(0)));
        let _ = cfg.j_pool.push(j_allele("IGKJ1*01", b"TTCAAACCC", Some(0)));
        cfg.identity.locus = Some(locus.to_string());
        cfg
    }

    #[test]
    fn vj_cartridge_with_igh_locus_is_flagged() {
        // VJ topology declaring an IGH (VDJ) locus is structurally
        // wrong — recombination shape and locus disagree.
        let cfg = vj_cfg_with_identity_locus("IGH");
        let issues = cfg.validate();
        assert!(issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::LocusChainTypeMismatch { locus, .. } if locus == "IGH"
        )));
    }

    #[test]
    fn vj_cartridge_with_igk_locus_passes() {
        // Locus matches topology — no LocusChainTypeMismatch.
        let cfg = vj_cfg_with_identity_locus("IGK");
        let issues = cfg.validate();
        assert!(!issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::LocusChainTypeMismatch { .. }
        )));
    }

    #[test]
    fn vdj_cartridge_with_igk_locus_is_flagged() {
        let cfg = vdj_cfg_with_identity_locus("IGK");
        assert!(cfg.validate().iter().any(|i| matches!(
            i,
            RefDataValidationIssue::LocusChainTypeMismatch { locus, .. } if locus == "IGK"
        )));
    }

    #[test]
    fn locus_chain_type_mismatch_is_fatal() {
        let cfg = vj_cfg_with_identity_locus("IGH");
        let issue = cfg
            .validate()
            .into_iter()
            .find(|i| matches!(i, RefDataValidationIssue::LocusChainTypeMismatch { .. }))
            .expect("must emit LocusChainTypeMismatch");
        assert_eq!(issue.severity(), RefDataIssueSeverity::Fatal);
        // AllowCuratable mode cannot opt this out.
        cfg.validate_with_mode(RefDataValidationMode::AllowCuratable)
            .expect_err("Fatal mismatch must not be opt-outable");
    }

    #[test]
    fn unknown_locus_does_not_fail_validation() {
        let cfg = vj_cfg_with_identity_locus("XYZ");
        assert!(!cfg.validate().iter().any(|i| matches!(
            i,
            RefDataValidationIssue::LocusChainTypeMismatch { .. }
        )));
    }

    #[test]
    fn empty_identity_does_not_fail_validation() {
        // No locus declared → no mismatch issue regardless of chain.
        let cfg_vj = RefDataConfig::empty(ChainType::Vj);
        assert!(!cfg_vj.validate().iter().any(|i| matches!(
            i,
            RefDataValidationIssue::LocusChainTypeMismatch { .. }
        )));
        let cfg_vdj = RefDataConfig::empty(ChainType::Vdj);
        assert!(!cfg_vdj.validate().iter().any(|i| matches!(
            i,
            RefDataValidationIssue::LocusChainTypeMismatch { .. }
        )));
    }

    #[test]
    fn locus_is_case_insensitive() {
        // Lowercase / mixed-case locus still flags correctly.
        let cfg = vj_cfg_with_identity_locus("igh");
        let issues = cfg.validate();
        assert!(issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::LocusChainTypeMismatch { locus, .. } if locus == "IGH"
        )));
    }
}
