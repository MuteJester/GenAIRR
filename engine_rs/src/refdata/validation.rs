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
mod tests;
