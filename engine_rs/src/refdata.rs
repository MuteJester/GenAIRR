//! Reference data — immutable shared inputs to the simulation.
//!
//! ## What lives here
//!
//! Allele sequences, gene metadata, chain configuration. Loaded once
//! at simulator construction (from `*.v6dat` files per D10), shared
//! across every simulation that uses it. *Never mutated* by the
//! simulator; the per-simulation state lives on `AlleleInstance`.
//!
//! Empirical distributions over alleles (frequencies) and trim/NP
//! distributions live in `dist.rs` rather than here — they are sampled
//! against the reference data, not part of it.
//!
//! ## Scope
//!
//! Just the typed data shapes plus construction / lookup helpers.
//! No biology, no PyO3, no serde. `serde` derives let the types
//! round-trip through `bincode` (D10).

use crate::ir::Segment;

mod validation;

pub use validation::{
    AnchorRule, RefDataIssueSeverity, RefDataValidationErrors, RefDataValidationIssue,
    RefDataValidationMode, ReferenceAlphabet, ReferenceRules,
};

// ──────────────────────────────────────────────────────────────────
// ReferenceIdentity — self-describing cartridge metadata
// ──────────────────────────────────────────────────────────────────

/// Self-describing identity for a reference cartridge.
///
/// Every field is `Option<String>` because identity is opt-in:
/// synthetic test cartridges and user-built configs may carry none
/// of it, while bundled `Experiment.on("human_igh")` data has all
/// five fields populated by the loader. Identity participates in the
/// content hash — two cartridges with identical catalogues + rules
/// but different declared species/locus/source produce different
/// hashes — so trace-replay integrity and debugging both benefit.
///
/// The fields are deliberately strings, not enums: an enum would
/// pin biology into engine code, exactly the implicit-knowledge
/// problem the cartridge architecture exists to eliminate. The
/// Python loader's enum (`Species`, `ChainType`) is stringified at
/// the bridge.
#[derive(Clone, Debug, Default, Eq, PartialEq, Hash)]
pub struct ReferenceIdentity {
    /// Common species label (e.g., `"Human"`, `"Mouse"`,
    /// `"Cynomolgus Macaque"`). Free-form by design; whoever ships
    /// the cartridge defines the vocabulary.
    pub species: Option<String>,
    /// AIRR locus prefix (`"IGH"`, `"IGK"`, `"IGL"`, `"TRA"`,
    /// `"TRB"`, `"TRG"`, `"TRD"`). When set, the validator
    /// cross-checks it against `chain_type` and feeds
    /// [`ReferenceRules::for_locus`] when no explicit rules are
    /// attached.
    pub locus: Option<String>,
    /// Reference catalogue family — `"IMGT"`, `"OGRDB"`,
    /// `"AIRR-C"`, `"custom"`, etc. Useful for cross-source
    /// provenance.
    pub reference_set: Option<String>,
    /// Human-readable cartridge name (e.g., `"human_igh"`,
    /// `"my_custom_kappa"`). Often the path or alias the cartridge
    /// was loaded from.
    pub name: Option<String>,
    /// Loader / file / user that produced this cartridge
    /// (`"DataConfig"`, `"v6dat"`, `"user_python"`).
    pub source: Option<String>,
}

// ──────────────────────────────────────────────────────────────────
// AlleleId — typed u32 newtype, distinct from NucHandle / RegionHandle
// ──────────────────────────────────────────────────────────────────

/// Index into an `AllelePool`. Stable for the lifetime of the
/// `RefDataConfig` it was issued for.
///
/// Distinct type from `NucHandle` and `RegionHandle` so the compiler
/// catches handle-confusion at call sites — this is the same
/// discipline as the IR handles in §3 of the design doc.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct AlleleId(u32);

impl AlleleId {
    pub const fn new(idx: u32) -> Self {
        Self(idx)
    }
    pub const fn index(self) -> u32 {
        self.0
    }
    pub const fn as_usize(self) -> usize {
        self.0 as usize
    }
}

/// Stable, refdata-local identifier for a gene within one segment's
/// pool. Assigned in first-appearance order over the pool's alleles, so
/// it is deterministic for a fixed cartridge (and therefore safe to
/// record in the trace for replay, gated by `refdata_content_hash`).
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct GeneId(u32);

impl GeneId {
    pub const fn new(idx: u32) -> Self {
        Self(idx)
    }
    pub const fn index(self) -> u32 {
        self.0
    }
    pub const fn as_usize(self) -> usize {
        self.0 as usize
    }
}

/// Gene-level grouping over a single `AllelePool`, derived from each
/// allele's `gene` string. Built on demand; the pool itself is
/// unchanged, so `Allele` and `refdata_content_hash` are untouched.
#[derive(Clone, Debug)]
pub struct GeneIndex {
    names: Vec<String>,                                   // GeneId.index() -> gene name
    by_name: std::collections::HashMap<String, GeneId>,   // gene name -> GeneId
    alleles: Vec<Vec<AlleleId>>,                          // GeneId.index() -> alleles, pool order
    gene_of: Vec<GeneId>,                                 // AlleleId.index() -> GeneId
}

impl GeneIndex {
    /// Build the gene grouping from a pool. Genes are numbered in the
    /// order their first allele appears in the pool.
    pub fn build(pool: &AllelePool) -> Self {
        let mut names: Vec<String> = Vec::new();
        let mut by_name: std::collections::HashMap<String, GeneId> =
            std::collections::HashMap::new();
        let mut alleles: Vec<Vec<AlleleId>> = Vec::new();
        let mut gene_of: Vec<GeneId> = Vec::with_capacity(pool.len());
        for (id, allele) in pool.iter() {
            let gid = *by_name.entry(allele.gene.clone()).or_insert_with(|| {
                let g = GeneId::new(names.len() as u32);
                names.push(allele.gene.clone());
                alleles.push(Vec::new());
                g
            });
            alleles[gid.as_usize()].push(id);
            gene_of.push(gid);
        }
        Self {
            names,
            by_name,
            alleles,
            gene_of,
        }
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }
    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }
    pub fn gene_id(&self, name: &str) -> Option<GeneId> {
        self.by_name.get(name).copied()
    }
    pub fn gene_name(&self, g: GeneId) -> &str {
        &self.names[g.as_usize()]
    }
    pub fn alleles_of(&self, g: GeneId) -> &[AlleleId] {
        &self.alleles[g.as_usize()]
    }
    pub fn gene_of(&self, a: AlleleId) -> GeneId {
        self.gene_of[a.as_usize()]
    }
    pub fn genes(&self) -> impl Iterator<Item = (GeneId, &str)> {
        self.names
            .iter()
            .enumerate()
            .map(|(i, n)| (GeneId::new(i as u32), n.as_str()))
    }
}

// ──────────────────────────────────────────────────────────────────
// ChainType — VJ (light) vs VDJ (heavy)
// ──────────────────────────────────────────────────────────────────

/// The biological chain configuration.
///
/// `Vj` chains (light: kappa, lambda, TCR alpha/gamma) recombine V
/// and J only — there is no D segment, NP1 spans the V→J interval,
/// and NP2 / D pools are unused.
///
/// `Vdj` chains (heavy: IGH, TCR beta/delta) recombine V, D, and J
/// with NP1 between V and D and NP2 between D and J.
///
/// Other axes (allele frequency models, isotype constants, etc.)
/// are independent of `ChainType` and live elsewhere in the config.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum ChainType {
    Vj,
    Vdj,
}

impl ChainType {
    /// Whether this chain has a D segment (and therefore an NP2 region).
    pub const fn has_d(self) -> bool {
        matches!(self, ChainType::Vdj)
    }
}

// ──────────────────────────────────────────────────────────────────
// Allele — one germline allele entry
// ──────────────────────────────────────────────────────────────────

/// IMGT functional classification carried per-allele.
///
/// `Functional`, `Orf`, `Pseudogene` mirror IMGT's `F` / `ORF` / `P`
/// labels. `Unknown` is the explicit "annotation present but value
/// unrecognised" placeholder — distinct from absence.
///
/// At the `Allele` level the field is `Option<FunctionalStatus>`:
/// - `None` means *the cartridge did not provide this annotation*,
///   so curation policies that filter on status either keep or drop
///   the entry according to their `keep_unannotated` flag.
/// - `Some(Unknown)` means *the cartridge explicitly says unknown*,
///   so the entry is filterable as Unknown — different semantics from
///   no-annotation.
///
/// The enum is deliberately closed: bundled catalogues that ship with
/// any future label (e.g. `[F]` shadowed) round-trip into `Unknown`,
/// preserving the audit trail without leaking unknown labels into the
/// engine.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum FunctionalStatus {
    Functional,
    Orf,
    Pseudogene,
    Unknown,
}

impl FunctionalStatus {
    /// Stable lowercase identifier used in the content hash and
    /// curation tag strings. Stable across engine versions so trace
    /// files round-trip.
    pub fn as_str(self) -> &'static str {
        match self {
            FunctionalStatus::Functional => "functional",
            FunctionalStatus::Orf => "orf",
            FunctionalStatus::Pseudogene => "pseudogene",
            FunctionalStatus::Unknown => "unknown",
        }
    }
}

/// One germline allele in a reference data set.
///
/// Immutable once constructed. Per-simulation state (which allele was
/// sampled, what trim was applied, what ambiguity set the post-trim
/// retained bases project to) lives on `AlleleInstance`.
///
/// **Field discipline:**
/// - `seq` is uppercase (`b'A'`, `b'C'`, `b'G'`, `b'T'`). Mixed case
///   is not allowed — case carries semantic meaning later in the
///   simulation (lowercase = SHM-mutated, etc.) and should not leak
///   into reference data.
/// - `anchor` is `Some(pos)` when the allele has a conserved codon
///   (Cys for V, W/F for J) at position `pos` (0-indexed within
///   `seq`). `None` for anchorless / partial alleles. The
///   `AnchorPreserved` contract reads this field.
/// - `name` is the canonical allele identifier (e.g.,
///   `"IGHV1-2*01"`). `gene` is the truncation to the gene level
///   (e.g., `"IGHV1-2"`).
/// - `functional_status` is `Some(...)` when the cartridge carries
///   an IMGT-style annotation, `None` when it doesn't. The
///   distinction matters for `RefDataCurationPolicy::FunctionalStatus`
///   which can opt to keep or drop unannotated entries.
/// - `subregions` holds IMGT-derived V-region substructure
///   (FWR1 / CDR1 / FWR2 / CDR2 / FWR3) when the cartridge carries
///   gapped sequence data. Empty `Vec` for D / J / C alleles and
///   for V alleles without gapped-sequence metadata. See
///   `docs/v_region_substructure_audit.md`.
#[derive(Clone, Debug)]
pub struct Allele {
    pub name: String,
    pub gene: String,
    pub seq: Vec<u8>,
    pub segment: Segment,
    pub anchor: Option<u16>,
    pub functional_status: Option<FunctionalStatus>,
    pub subregions: Vec<VSubregion>,
}

/// Canonical V-region substructure labels (IMGT convention).
///
/// Discriminants are explicit so the type can be used as a
/// `PerLabel`-style array index and serialised stably to wire
/// formats. Order matches biological 5'→3' order on the V allele.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
#[repr(u8)]
pub enum VSubregionLabel {
    Fwr1 = 0,
    Cdr1 = 1,
    Fwr2 = 2,
    Cdr2 = 3,
    Fwr3 = 4,
}

impl VSubregionLabel {
    /// Canonical string name (uppercase IMGT label). Used by the
    /// content hash, PyO3 serialisation, and the manifest layer.
    pub const fn as_str(&self) -> &'static str {
        match self {
            VSubregionLabel::Fwr1 => "FWR1",
            VSubregionLabel::Cdr1 => "CDR1",
            VSubregionLabel::Fwr2 => "FWR2",
            VSubregionLabel::Cdr2 => "CDR2",
            VSubregionLabel::Fwr3 => "FWR3",
        }
    }

    /// Parse a canonical string name back into the enum. Returns
    /// `None` for unknown labels (case-sensitive — the bridge
    /// normalises to uppercase before reaching this point).
    pub fn from_str(label: &str) -> Option<Self> {
        match label {
            "FWR1" => Some(VSubregionLabel::Fwr1),
            "CDR1" => Some(VSubregionLabel::Cdr1),
            "FWR2" => Some(VSubregionLabel::Fwr2),
            "CDR2" => Some(VSubregionLabel::Cdr2),
            "FWR3" => Some(VSubregionLabel::Fwr3),
            _ => None,
        }
    }
}

/// One V-region substructure annotation. Coordinates are ungapped,
/// allele-relative, half-open. The `RefDataConfig::validate` pass
/// rejects malformed annotations (start >= end, out-of-bounds end,
/// overlap with another annotation on the same allele, duplicate
/// label).
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct VSubregion {
    pub label: VSubregionLabel,
    pub start: u16,
    pub end: u16,
}

impl Allele {
    /// Length of the allele sequence in nucleotides.
    pub fn len(&self) -> u32 {
        self.seq.len() as u32
    }

    /// Whether the allele has a known anchor position.
    pub fn has_anchor(&self) -> bool {
        self.anchor.is_some()
    }
}

// ──────────────────────────────────────────────────────────────────
// AllelePool — flat indexed collection of alleles for one segment
// ──────────────────────────────────────────────────────────────────

/// All alleles available for one segment role (V, D, J, or C).
///
/// Indexed by `AlleleId`. The pool is the reference data; *which*
/// allele a particular simulation sampled is recorded by an
/// `AlleleInstance` referring to a stable `AlleleId`.
#[derive(Clone, Debug, Default)]
pub struct AllelePool {
    alleles: Vec<Allele>,
}

impl AllelePool {
    pub fn new() -> Self {
        Self::default()
    }

    /// Construct from an existing `Vec<Allele>`. Caller is responsible
    /// for ensuring all alleles have the right segment role for the
    /// pool being built.
    pub fn from_vec(alleles: Vec<Allele>) -> Self {
        Self { alleles }
    }

    /// Number of alleles in the pool.
    pub fn len(&self) -> usize {
        self.alleles.len()
    }

    /// Whether the pool contains zero alleles.
    pub fn is_empty(&self) -> bool {
        self.alleles.is_empty()
    }

    /// Append an allele. Returns the issued `AlleleId` so callers can
    /// reference the just-added allele in subsequent setup. This is
    /// the construction-time API; reference data is sealed before any
    /// simulation runs against it.
    #[must_use = "AllelePool::push returns the issued AlleleId; use \
                  `_ = pool.push(...)` if it isn't needed"]
    pub fn push(&mut self, allele: Allele) -> AlleleId {
        let id = AlleleId::new(self.alleles.len() as u32);
        self.alleles.push(allele);
        id
    }

    /// Look up an allele by id. Returns `None` if `id` is out of
    /// bounds — defensive to allow `RefDataConfig` to expose
    /// fallible lookup at the boundary even though correct callers
    /// always have valid handles.
    pub fn get(&self, id: AlleleId) -> Option<&Allele> {
        self.alleles.get(id.as_usize())
    }

    /// Read-only slice of all alleles in pool order. Iterated by
    /// `AllelePoolDist` to build cumulative frequency tables.
    pub fn as_slice(&self) -> &[Allele] {
        &self.alleles
    }

    /// Iterator over `(AlleleId, &Allele)` pairs in pool order.
    pub fn iter(&self) -> impl Iterator<Item = (AlleleId, &Allele)> {
        self.alleles
            .iter()
            .enumerate()
            .map(|(i, a)| (AlleleId::new(i as u32), a))
    }

    /// Find the first allele whose name matches exactly. O(N).
    /// Used in tests and tooling; production sampling goes through
    /// `AllelePoolDist` (C.3) and never name-resolves here.
    pub fn find_by_name(&self, name: &str) -> Option<(AlleleId, &Allele)> {
        self.iter().find(|(_, a)| a.name == name)
    }
}

// ──────────────────────────────────────────────────────────────────
// RefDataConfig — top-level container for one species/chain
// ──────────────────────────────────────────────────────────────────

/// Top-level immutable reference configuration for one simulation
/// target (e.g., human IGH).
///
/// In production, this is loaded from a `.v6dat` file (D10) at
/// `Experiment.compile()` time. In tests, build directly via
/// `RefDataConfig::builder()` or by constructing the fields
/// manually.
///
/// **Invariant:** for `chain_type == ChainType::Vj`, the `d_pool` is
/// expected to be empty. Construction does not enforce this — callers
/// of the builder API in tests are expected to honor it; the
/// assembly pass in C.8 will explicitly skip D for VJ chains.
#[derive(Clone, Debug)]
pub struct RefDataConfig {
    pub chain_type: ChainType,
    /// Self-describing identity (species, locus, reference set,
    /// name, source). Default is empty; the bundled-data loader and
    /// user code populate it. See [`ReferenceIdentity`].
    pub identity: ReferenceIdentity,
    /// Programmable interpretation rules — anchor expectations,
    /// allowed alphabet. Drives both the validator and downstream
    /// projection layers. Default constructs a lenient set; bundled
    /// loaders override with locus-appropriate rules. See
    /// [`ReferenceRules`].
    pub rules: ReferenceRules,
    pub v_pool: AllelePool,
    pub d_pool: AllelePool,
    pub j_pool: AllelePool,
    pub c_pool: AllelePool,
}

impl RefDataConfig {
    /// Empty config for the given chain type. Use the builder /
    /// direct field assignment to populate the pools. `identity`
    /// initialises to [`ReferenceIdentity::default`] (all fields
    /// `None`); `rules` initialises to [`ReferenceRules::default`]
    /// (lenient alphabet, Cys V, W-or-F J).
    pub fn empty(chain_type: ChainType) -> Self {
        Self {
            chain_type,
            identity: ReferenceIdentity::default(),
            rules: ReferenceRules::default(),
            v_pool: AllelePool::new(),
            d_pool: AllelePool::new(),
            j_pool: AllelePool::new(),
            c_pool: AllelePool::new(),
        }
    }

    /// Resolve an allele by its segment role and id. Returns `None`
    /// if the segment isn't in this config's pools (e.g., asking for
    /// a D allele on a VJ chain) or if the id is out of bounds.
    pub fn get(&self, segment: Segment, id: AlleleId) -> Option<&Allele> {
        match segment {
            Segment::V => self.v_pool.get(id),
            Segment::D => self.d_pool.get(id),
            Segment::J => self.j_pool.get(id),
            // NP regions don't have alleles; return None defensively
            // rather than panicking — a caller asking for an NP
            // allele is misusing the API.
            Segment::Np1 | Segment::Np2 => None,
        }
    }

    /// Pool for the given segment, or `None` for NP segments.
    pub fn pool_for(&self, segment: Segment) -> Option<&AllelePool> {
        match segment {
            Segment::V => Some(&self.v_pool),
            Segment::D => Some(&self.d_pool),
            Segment::J => Some(&self.j_pool),
            Segment::Np1 | Segment::Np2 => None,
        }
    }

    /// Return a curated copy of this cartridge under `policy`.
    ///
    /// Curation selects which alleles participate in simulation; it
    /// does NOT fix structural problems. Specifically:
    ///
    /// - [`RefDataCurationPolicy::Raw`] is identity — returns a clone.
    /// - [`RefDataCurationPolicy::FunctionalAnchorsOnly`] removes V
    ///   and J alleles that fail the active [`AnchorRule`] for their
    ///   segment: missing anchor, anchor codon out of bounds, or
    ///   anchor codon translating to an amino acid outside
    ///   `expected_amino_acids`. D and C pools pass through unchanged.
    ///
    /// Fatal structural issues (duplicate names, invalid sequence
    /// bytes, `LocusChainTypeMismatch`) are NOT filtered — they
    /// surface from [`Self::validate`] both before and after
    /// curation. Curation must never silently hide corruption.
    ///
    /// Identity carries a curation tag: `identity.source` is
    /// extended with `|curated:<policy>` (or set to that string if
    /// previously empty). This guarantees a curated cartridge's
    /// content hash differs from the raw one and that downstream
    /// trace files attribute outputs to the right artefact.
    ///
    /// Curation can leave a required pool empty — strict compile
    /// will then fail with `EmptyRequiredPool`, which is the
    /// intended diagnostic (the catalogue couldn't support a
    /// functional simulation under these rules).
    pub fn curated(&self, policy: RefDataCurationPolicy) -> RefDataConfig {
        let mut out = self.clone();
        match &policy {
            RefDataCurationPolicy::Raw => {
                // No filtering, no provenance tag. A no-op pass; the
                // user explicitly chose to bypass curation.
                out
            }
            RefDataCurationPolicy::FunctionalAnchorsOnly => {
                out.v_pool = filter_pool_by_anchor_rule(&self.v_pool, &self.rules.v_anchor);
                out.j_pool = filter_pool_by_anchor_rule(&self.j_pool, &self.rules.j_anchor);
                // D + C pools pass through unchanged.
                tag_identity_source(&mut out.identity, &policy.tag());
                out
            }
            RefDataCurationPolicy::FunctionalStatus { allowed, keep_unannotated } => {
                out.v_pool =
                    filter_pool_by_functional_status(&self.v_pool, allowed, *keep_unannotated);
                out.d_pool =
                    filter_pool_by_functional_status(&self.d_pool, allowed, *keep_unannotated);
                out.j_pool =
                    filter_pool_by_functional_status(&self.j_pool, allowed, *keep_unannotated);
                // C pool passes through unchanged.
                tag_identity_source(&mut out.identity, &policy.tag());
                out
            }
        }
    }
}

/// Curation policy — which subset of the catalogue participates in
/// simulation.
///
/// Curation is distinct from validation: validation **describes** the
/// catalogue, curation **selects** from it. A pseudogene-bearing
/// catalogue can fail strict validation but pass a curated
/// simulation; conversely, structural corruption (duplicate names,
/// invalid bytes) is not fixable by curation and continues to fail
/// validation either way.
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub enum RefDataCurationPolicy {
    /// Identity policy — keep every allele. Existing behaviour.
    Raw,
    /// Drop V/J alleles whose anchor doesn't satisfy the active
    /// [`AnchorRule`]: missing anchor, anchor out of bounds, or
    /// anchor codon AA outside `expected_amino_acids`. D and C
    /// pools pass through unchanged.
    FunctionalAnchorsOnly,
    /// Filter V/D/J pools by per-allele [`FunctionalStatus`].
    ///
    /// `allowed` is the set of statuses to keep (e.g.
    /// `[Functional]`, `[Functional, Orf]`).
    ///
    /// `keep_unannotated` controls what happens to alleles whose
    /// `functional_status == None`. Bundled `.pkl` cartridges leave
    /// status unannotated for now; setting `keep_unannotated=true`
    /// preserves backward compatibility with those catalogues, while
    /// `keep_unannotated=false` lets a curated cartridge enforce
    /// "annotation required".
    ///
    /// C pool passes through unchanged.
    FunctionalStatus {
        allowed: Vec<FunctionalStatus>,
        keep_unannotated: bool,
    },
}

impl RefDataCurationPolicy {
    /// Short identifier used in [`ReferenceIdentity::source`]
    /// curation tags. Stable strings so trace files written with one
    /// engine version replay under another.
    ///
    /// For [`Self::FunctionalStatus`] the tag is the parametrised
    /// form `functional_status:functional,orf|keep_unannotated=true`
    /// (allowed list joined by `,` in canonical lowercase order,
    /// followed by the `keep_unannotated` flag). Two inputs that
    /// produce identical curated catalogues produce identical tags.
    pub fn tag(&self) -> String {
        match self {
            RefDataCurationPolicy::Raw => "raw".to_string(),
            RefDataCurationPolicy::FunctionalAnchorsOnly => {
                "functional_anchors_only".to_string()
            }
            RefDataCurationPolicy::FunctionalStatus { allowed, keep_unannotated } => {
                let mut canonical: Vec<&'static str> =
                    allowed.iter().map(|s| s.as_str()).collect();
                canonical.sort_unstable();
                canonical.dedup();
                format!(
                    "functional_status:{}|keep_unannotated={}",
                    canonical.join(","),
                    keep_unannotated,
                )
            }
        }
    }
}

/// Apply `rule` against each allele in `pool`; keep alleles whose
/// anchor passes the rule.
///
/// "Passes" means:
/// - if the rule has `required=false`, anchorless alleles are kept;
/// - the anchor position (if present) must leave room for a full
///   codon (`anchor + 3 <= seq.len()`);
/// - the anchor codon must translate to an AA in
///   `expected_amino_acids`.
///
/// Alleles whose sequence is shorter than 3 bytes or whose codon
/// can't translate (synthetic non-DNA bytes) fail the rule
/// regardless — they can't participate in a functional projection.
fn filter_pool_by_anchor_rule(pool: &AllelePool, rule: &validation::AnchorRule) -> AllelePool {
    let kept: Vec<Allele> = pool
        .iter()
        .filter(|(_, allele)| anchor_passes_rule(allele, rule))
        .map(|(_, a)| a.clone())
        .collect();
    AllelePool::from_vec(kept)
}

/// Apply functional-status filtering to a pool.
///
/// An allele is kept iff:
/// - `allele.functional_status == Some(s)` and `s` is in `allowed`; OR
/// - `allele.functional_status.is_none()` and `keep_unannotated` is true.
///
/// Order is preserved so trace-replay determinism is unaffected by
/// the filter step.
fn filter_pool_by_functional_status(
    pool: &AllelePool,
    allowed: &[FunctionalStatus],
    keep_unannotated: bool,
) -> AllelePool {
    let kept: Vec<Allele> = pool
        .iter()
        .filter(|(_, allele)| match allele.functional_status {
            Some(s) => allowed.contains(&s),
            None => keep_unannotated,
        })
        .map(|(_, a)| a.clone())
        .collect();
    AllelePool::from_vec(kept)
}

fn anchor_passes_rule(allele: &Allele, rule: &validation::AnchorRule) -> bool {
    use crate::codon::translate_codon;
    let Some(anchor) = allele.anchor else {
        // No anchor: only keep when the rule treats anchors as optional.
        return !rule.required;
    };
    let len = allele.seq.len() as u32;
    let anchor_u32 = anchor as u32;
    if anchor_u32 + 3 > len {
        return false;
    }
    let a = anchor as usize;
    let codon = [allele.seq[a], allele.seq[a + 1], allele.seq[a + 2]];
    let aa = translate_codon(codon[0], codon[1], codon[2]) as char;
    rule.expected_amino_acids.contains(&aa)
}

/// Append `|curated:<tag>` to `identity.source`, or set it to
/// `curated:<tag>` if no source was previously declared. Repeated
/// curation appends each pass so trace consumers can read the full
/// curation chain.
fn tag_identity_source(identity: &mut ReferenceIdentity, tag: &str) {
    let suffix = format!("curated:{tag}");
    identity.source = Some(match identity.source.as_deref() {
        Some(prev) if !prev.is_empty() => format!("{prev}|{suffix}"),
        _ => suffix,
    });
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::mem::size_of;

    /// Helper: a tiny synthetic V allele for tests.
    fn make_v(name: &str, gene: &str, seq: &[u8], anchor: Option<u16>) -> Allele {
        Allele {
            name: name.to_string(),
            gene: gene.to_string(),
            seq: seq.to_vec(),
            segment: Segment::V,
            anchor,
            functional_status: None,
            subregions: Vec::new(),
        }
    }

    #[test]
    fn allele_id_is_zero_cost_newtype() {
        assert_eq!(size_of::<AlleleId>(), size_of::<u32>());
    }

    #[test]
    fn allele_id_round_trip() {
        let id = AlleleId::new(7);
        assert_eq!(id.index(), 7);
        assert_eq!(id.as_usize(), 7);
    }

    #[test]
    fn chain_type_has_d_distinction() {
        assert!(!ChainType::Vj.has_d());
        assert!(ChainType::Vdj.has_d());
    }

    #[test]
    fn allele_basic_accessors() {
        let a = make_v("IGHV1-2*01", "IGHV1-2", b"ACGTACGT", Some(3));
        assert_eq!(a.len(), 8);
        assert!(a.has_anchor());
        assert_eq!(a.anchor, Some(3));
        assert_eq!(a.segment, Segment::V);
    }

    #[test]
    fn allele_anchorless_round_trip() {
        let a = make_v("IGHV-pseudo*01", "IGHV-pseudo", b"ACGT", None);
        assert!(!a.has_anchor());
        assert_eq!(a.anchor, None);
    }

    #[test]
    fn allele_pool_starts_empty() {
        let p = AllelePool::new();
        assert_eq!(p.len(), 0);
        assert!(p.is_empty());
        assert!(p.get(AlleleId::new(0)).is_none());
        assert!(p.find_by_name("nonexistent").is_none());
    }

    #[test]
    fn allele_pool_push_returns_sequential_ids() {
        let mut p = AllelePool::new();
        let id0 = p.push(make_v("a*01", "a", b"AA", Some(0)));
        let id1 = p.push(make_v("b*01", "b", b"CC", Some(0)));
        let id2 = p.push(make_v("c*01", "c", b"GG", None));

        assert_eq!(id0.index(), 0);
        assert_eq!(id1.index(), 1);
        assert_eq!(id2.index(), 2);
        assert_eq!(p.len(), 3);
    }

    #[test]
    fn allele_pool_get_returns_stored_allele() {
        let mut p = AllelePool::new();
        let id = p.push(make_v("x*01", "x", b"ATGC", Some(2)));
        let got = p.get(id).expect("just pushed");
        assert_eq!(got.name, "x*01");
        assert_eq!(got.seq, b"ATGC");
        assert_eq!(got.anchor, Some(2));
    }

    #[test]
    fn allele_pool_get_out_of_bounds_returns_none() {
        let p = AllelePool::new();
        assert!(p.get(AlleleId::new(99)).is_none());
    }

    #[test]
    fn allele_pool_iter_yields_id_allele_pairs_in_order() {
        let mut p = AllelePool::new();
        let _ = p.push(make_v("a*01", "a", b"AA", None));
        let _ = p.push(make_v("b*01", "b", b"CC", None));
        let _ = p.push(make_v("c*01", "c", b"GG", None));

        let collected: Vec<(u32, String)> = p
            .iter()
            .map(|(id, a)| (id.index(), a.name.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![
                (0, "a*01".to_string()),
                (1, "b*01".to_string()),
                (2, "c*01".to_string()),
            ]
        );
    }

    #[test]
    fn allele_pool_find_by_name_locates_exact_match() {
        let mut p = AllelePool::new();
        let _ = p.push(make_v("IGHV1-2*01", "IGHV1-2", b"AA", None));
        let target_id = p.push(make_v("IGHV1-2*02", "IGHV1-2", b"AC", None));
        let _ = p.push(make_v("IGHV3-23*01", "IGHV3-23", b"GG", None));

        let (id, allele) = p.find_by_name("IGHV1-2*02").expect("name should exist");
        assert_eq!(id, target_id);
        assert_eq!(allele.seq, b"AC");

        // Partial / similar names should not match.
        assert!(p.find_by_name("IGHV1-2").is_none());
        assert!(p.find_by_name("IGHV1-2*03").is_none());
    }

    #[test]
    fn ref_data_config_empty_for_chain_type() {
        let cfg = RefDataConfig::empty(ChainType::Vdj);
        assert_eq!(cfg.chain_type, ChainType::Vdj);
        assert!(cfg.v_pool.is_empty());
        assert!(cfg.d_pool.is_empty());
        assert!(cfg.j_pool.is_empty());
        assert!(cfg.c_pool.is_empty());
    }

    #[test]
    fn ref_data_config_pool_for_segment_routes_correctly() {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.v_pool.push(make_v("v*01", "v", b"AA", None));
        let _ = cfg.d_pool.push(Allele {
            name: "d*01".into(),
            gene: "d".into(),
            seq: b"GG".to_vec(),
            segment: Segment::D,
            anchor: None,
            functional_status: None,
            subregions: Vec::new(),
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j*01".into(),
            gene: "j".into(),
            seq: b"TT".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
            functional_status: None,
            subregions: Vec::new(),
        });

        assert_eq!(cfg.pool_for(Segment::V).unwrap().len(), 1);
        assert_eq!(cfg.pool_for(Segment::D).unwrap().len(), 1);
        assert_eq!(cfg.pool_for(Segment::J).unwrap().len(), 1);
        assert!(cfg.pool_for(Segment::Np1).is_none());
        assert!(cfg.pool_for(Segment::Np2).is_none());
    }

    #[test]
    fn ref_data_config_get_resolves_segment_id_pair() {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let v_id = cfg.v_pool.push(make_v("v*01", "v", b"AAAT", Some(1)));

        let v = cfg.get(Segment::V, v_id).expect("v*01 should resolve");
        assert_eq!(v.name, "v*01");
        assert_eq!(v.anchor, Some(1));

        // Wrong segment -> None.
        assert!(cfg.get(Segment::J, v_id).is_none());
        // NP segment -> None defensively.
        assert!(cfg.get(Segment::Np1, v_id).is_none());
    }

    // ── Curation policy ───────────────────────────────────────────

    fn allele_with_anchor(name: &str, seq: &[u8], anchor: Option<u16>, segment: Segment) -> Allele {
        Allele {
            name: name.to_string(),
            gene: name.split('*').next().unwrap_or(name).to_string(),
            seq: seq.to_vec(),
            segment,
            anchor,
            functional_status: None,
            subregions: Vec::new(),
        }
    }

    fn vj_cfg_for_curation() -> RefDataConfig {
        // Two V alleles: one Cys-anchored, one Gly-anchored.
        // Two J alleles: one W-anchored, one anchorless.
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg
            .v_pool
            .push(allele_with_anchor("MYV-good*01", b"TGTAAACCC", Some(0), Segment::V));
        let _ = cfg
            .v_pool
            .push(allele_with_anchor("MYV-gly*01", b"GGGAAACCC", Some(0), Segment::V));
        let _ = cfg
            .j_pool
            .push(allele_with_anchor("MYJ-good*01", b"TGGAAACCC", Some(0), Segment::J));
        let _ = cfg
            .j_pool
            .push(allele_with_anchor("MYJ-orphan*01", b"GGG", None, Segment::J));
        cfg
    }

    #[test]
    fn curation_raw_is_identity() {
        let cfg = vj_cfg_for_curation();
        let curated = cfg.curated(RefDataCurationPolicy::Raw);
        assert_eq!(curated.v_pool.len(), 2);
        assert_eq!(curated.j_pool.len(), 2);
        // Raw policy does NOT touch identity.
        assert_eq!(curated.identity.source, None);
    }

    #[test]
    fn curation_functional_anchors_only_filters_v_and_j() {
        let cfg = vj_cfg_for_curation();
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        assert_eq!(curated.v_pool.len(), 1);
        assert_eq!(curated.v_pool.iter().next().unwrap().1.name, "MYV-good*01");
        assert_eq!(curated.j_pool.len(), 1);
        assert_eq!(curated.j_pool.iter().next().unwrap().1.name, "MYJ-good*01");
    }

    #[test]
    fn curation_tags_identity_source() {
        let cfg = vj_cfg_for_curation();
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        assert_eq!(
            curated.identity.source.as_deref(),
            Some("curated:functional_anchors_only"),
        );
    }

    #[test]
    fn curation_extends_existing_identity_source() {
        let mut cfg = vj_cfg_for_curation();
        cfg.identity.source = Some("DataConfig".to_string());
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        assert_eq!(
            curated.identity.source.as_deref(),
            Some("DataConfig|curated:functional_anchors_only"),
        );
    }

    #[test]
    fn curation_preserves_other_identity_fields() {
        let mut cfg = vj_cfg_for_curation();
        cfg.identity.species = Some("Human".into());
        cfg.identity.locus = Some("IGK".into());
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        assert_eq!(curated.identity.species.as_deref(), Some("Human"));
        assert_eq!(curated.identity.locus.as_deref(), Some("IGK"));
    }

    #[test]
    fn curation_does_not_silently_fix_duplicate_names() {
        // Add a duplicate V allele AFTER the good one. Both are
        // Cys-anchored — they survive curation, then validation
        // surfaces the duplicate as a Fatal issue. Curation never
        // removes a structurally-corrupt entry to make a catalogue
        // look healthier.
        let mut cfg = vj_cfg_for_curation();
        let _ = cfg.v_pool.push(allele_with_anchor(
            "MYV-good*01",
            b"TGTAAACCC",
            Some(0),
            Segment::V,
        ));
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        let issues = curated.validate();
        assert!(issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::DuplicateAlleleName { .. }
        )));
    }

    #[test]
    fn curation_does_not_silently_fix_invalid_bytes() {
        // Add an allele with a gap byte ('.') and a valid Cys anchor.
        // The invalid byte is Fatal and must still surface post-curation.
        let mut cfg = vj_cfg_for_curation();
        let _ = cfg.v_pool.push(allele_with_anchor(
            "MYV-badbyte*01",
            b"TGT.AAACC",
            Some(0),
            Segment::V,
        ));
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        assert!(curated.validate().iter().any(|i| matches!(
            i,
            RefDataValidationIssue::InvalidAlleleByte { .. }
        )));
    }

    #[test]
    fn curation_to_empty_pool_surfaces_empty_required_pool() {
        // Every V allele is non-Cys — curation empties V.
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg
            .v_pool
            .push(allele_with_anchor("MYV-gly*01", b"GGGAAACCC", Some(0), Segment::V));
        let _ = cfg
            .j_pool
            .push(allele_with_anchor("MYJ-good*01", b"TGGAAACCC", Some(0), Segment::J));
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        assert_eq!(curated.v_pool.len(), 0);
        let issues = curated.validate();
        assert!(issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::EmptyRequiredPool { segment: Segment::V }
        )));
    }

    #[test]
    fn curation_d_pool_passes_through_unchanged() {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg
            .v_pool
            .push(allele_with_anchor("MYV*01", b"TGTAAACCC", Some(0), Segment::V));
        let _ = cfg.d_pool.push(allele_with_anchor(
            "MYD*01",
            b"GGGCCCAAA",
            None,
            Segment::D,
        ));
        let _ = cfg
            .j_pool
            .push(allele_with_anchor("MYJ*01", b"TGGAAACCC", Some(0), Segment::J));
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        // Anchorless D survives — D pool is never filtered.
        assert_eq!(curated.d_pool.len(), 1);
    }

    #[test]
    fn curation_respects_custom_j_anchor_rule() {
        // Rule expects Y at J anchor. TGG (W) and TTC (F) are dropped;
        // TAT (Y) is kept.
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules.j_anchor.expected_amino_acids = vec!['Y'];
        let _ = cfg
            .v_pool
            .push(allele_with_anchor("MYV*01", b"TGTAAACCC", Some(0), Segment::V));
        let _ = cfg
            .j_pool
            .push(allele_with_anchor("MYJ-w*01", b"TGGAAACCC", Some(0), Segment::J));
        let _ = cfg
            .j_pool
            .push(allele_with_anchor("MYJ-y*01", b"TATAAACCC", Some(0), Segment::J));
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        assert_eq!(curated.j_pool.len(), 1);
        assert_eq!(curated.j_pool.iter().next().unwrap().1.name, "MYJ-y*01");
    }

    #[test]
    fn curation_anchor_required_false_keeps_anchorless_alleles() {
        // If the rule says anchor is optional, anchorless alleles
        // are kept under FunctionalAnchorsOnly.
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        cfg.rules.j_anchor.required = false;
        let _ = cfg
            .v_pool
            .push(allele_with_anchor("MYV*01", b"TGTAAACCC", Some(0), Segment::V));
        let _ = cfg
            .j_pool
            .push(allele_with_anchor("MYJ-orphan*01", b"GGG", None, Segment::J));
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalAnchorsOnly);
        assert_eq!(curated.j_pool.len(), 1);
    }

    // ── Functional-status curation ────────────────────────────────

    fn allele_with_status(
        name: &str,
        seq: &[u8],
        segment: Segment,
        anchor: Option<u16>,
        status: Option<FunctionalStatus>,
    ) -> Allele {
        Allele {
            name: name.to_string(),
            gene: name.split('*').next().unwrap_or(name).to_string(),
            seq: seq.to_vec(),
            segment,
            anchor,
            functional_status: status,
            subregions: Vec::new(),
        }
    }

    fn mixed_status_cfg() -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.v_pool.push(allele_with_status(
            "v-f*01", b"TGTAAACCC", Segment::V, Some(0),
            Some(FunctionalStatus::Functional),
        ));
        let _ = cfg.v_pool.push(allele_with_status(
            "v-o*01", b"TGTAAACCC", Segment::V, Some(0),
            Some(FunctionalStatus::Orf),
        ));
        let _ = cfg.v_pool.push(allele_with_status(
            "v-p*01", b"TGTAAACCC", Segment::V, Some(0),
            Some(FunctionalStatus::Pseudogene),
        ));
        let _ = cfg.v_pool.push(allele_with_status(
            "v-na*01", b"TGTAAACCC", Segment::V, Some(0), None,
        ));
        let _ = cfg.d_pool.push(allele_with_status(
            "d-f*01", b"GGG", Segment::D, None,
            Some(FunctionalStatus::Functional),
        ));
        let _ = cfg.d_pool.push(allele_with_status(
            "d-p*01", b"GGG", Segment::D, None,
            Some(FunctionalStatus::Pseudogene),
        ));
        let _ = cfg.j_pool.push(allele_with_status(
            "j-f*01", b"TGGAAACCC", Segment::J, Some(0),
            Some(FunctionalStatus::Functional),
        ));
        cfg
    }

    #[test]
    fn curation_functional_status_filters_v_d_j() {
        let cfg = mixed_status_cfg();
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalStatus {
            allowed: vec![FunctionalStatus::Functional],
            keep_unannotated: false,
        });
        let v_names: Vec<&str> =
            curated.v_pool.iter().map(|(_, a)| a.name.as_str()).collect();
        let d_names: Vec<&str> =
            curated.d_pool.iter().map(|(_, a)| a.name.as_str()).collect();
        let j_names: Vec<&str> =
            curated.j_pool.iter().map(|(_, a)| a.name.as_str()).collect();
        assert_eq!(v_names, vec!["v-f*01"]);
        assert_eq!(d_names, vec!["d-f*01"]);
        assert_eq!(j_names, vec!["j-f*01"]);
    }

    #[test]
    fn curation_functional_status_keeps_unannotated_when_flagged() {
        let cfg = mixed_status_cfg();
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalStatus {
            allowed: vec![FunctionalStatus::Functional],
            keep_unannotated: true,
        });
        let v_names: Vec<&str> =
            curated.v_pool.iter().map(|(_, a)| a.name.as_str()).collect();
        assert_eq!(v_names, vec!["v-f*01", "v-na*01"]);
    }

    #[test]
    fn curation_functional_status_tag_is_canonical() {
        let p = RefDataCurationPolicy::FunctionalStatus {
            allowed: vec![FunctionalStatus::Orf, FunctionalStatus::Functional],
            keep_unannotated: false,
        };
        // allowed list is sorted into canonical lexical order; the
        // policy tag is the source of identity-source provenance, so
        // two policies producing identical curated catalogues must
        // produce identical tags regardless of input order.
        assert_eq!(
            p.tag(),
            "functional_status:functional,orf|keep_unannotated=false",
        );
    }

    #[test]
    fn curation_functional_status_dedupes_allowed_in_tag() {
        let p = RefDataCurationPolicy::FunctionalStatus {
            allowed: vec![
                FunctionalStatus::Functional,
                FunctionalStatus::Functional,
            ],
            keep_unannotated: true,
        };
        assert_eq!(p.tag(), "functional_status:functional|keep_unannotated=true");
    }

    #[test]
    fn curation_functional_status_empty_v_surfaces_empty_required_pool() {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(allele_with_status(
            "v-p*01", b"TGTAAACCC", Segment::V, Some(0),
            Some(FunctionalStatus::Pseudogene),
        ));
        let _ = cfg.j_pool.push(allele_with_status(
            "j-f*01", b"TGGAAACCC", Segment::J, Some(0),
            Some(FunctionalStatus::Functional),
        ));
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalStatus {
            allowed: vec![FunctionalStatus::Functional],
            keep_unannotated: false,
        });
        assert!(curated.v_pool.is_empty());
        let issues = curated.validate();
        assert!(issues.iter().any(|i| matches!(
            i,
            RefDataValidationIssue::EmptyRequiredPool { segment: Segment::V }
        )));
    }

    #[test]
    fn curation_functional_status_tags_identity_source() {
        let cfg = mixed_status_cfg();
        let curated = cfg.curated(RefDataCurationPolicy::FunctionalStatus {
            allowed: vec![FunctionalStatus::Functional],
            keep_unannotated: true,
        });
        let src = curated.identity.source.as_deref().unwrap_or("");
        assert!(src.contains("curated:functional_status:functional|keep_unannotated=true"));
    }

    #[test]
    fn ref_data_config_supports_vj_chain_with_empty_d_pool() {
        // VJ chains: d_pool is conventionally empty. Construction
        // does not enforce this; assembly (C.8) handles VJ vs VDJ
        // explicitly via chain_type.has_d().
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(make_v("v*01", "v", b"AA", Some(0)));
        let _ = cfg.j_pool.push(Allele {
            name: "j*01".into(),
            gene: "j".into(),
            seq: b"TT".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
            functional_status: None,
            subregions: Vec::new(),
        });

        assert!(!cfg.chain_type.has_d());
        assert!(cfg.d_pool.is_empty());
        assert_eq!(cfg.v_pool.len(), 1);
        assert_eq!(cfg.j_pool.len(), 1);
    }
}

#[cfg(test)]
mod gene_index_tests {
    use super::*;

    fn pool_with(genes: &[(&str, &str)]) -> AllelePool {
        // genes: (allele_name, gene_name)
        let mut p = AllelePool::new();
        for (name, gene) in genes {
            let _ = p.push(Allele {
                name: (*name).to_string(),
                gene: (*gene).to_string(),
                seq: vec![b'A'; 10],
                segment: Segment::V,
                anchor: Some(3),
                functional_status: None,
                subregions: Vec::new(),
            });
        }
        p
    }

    #[test]
    fn gene_index_groups_alleles_by_gene_in_first_appearance_order() {
        let pool = pool_with(&[
            ("IGHV1-2*01", "IGHV1-2"),
            ("IGHV1-2*02", "IGHV1-2"),
            ("IGHV3-23*01", "IGHV3-23"),
        ]);
        let idx = GeneIndex::build(&pool);

        assert_eq!(idx.len(), 2);
        let g12 = idx.gene_id("IGHV1-2").unwrap();
        let g323 = idx.gene_id("IGHV3-23").unwrap();
        assert_eq!(g12.index(), 0); // first appearance
        assert_eq!(g323.index(), 1);
        assert_eq!(idx.gene_name(g12), "IGHV1-2");
        assert_eq!(idx.alleles_of(g12), &[AlleleId::new(0), AlleleId::new(1)]);
        assert_eq!(idx.alleles_of(g323), &[AlleleId::new(2)]);
        assert_eq!(idx.gene_of(AlleleId::new(1)), g12);
        assert!(idx.gene_id("nope").is_none());
    }
}
