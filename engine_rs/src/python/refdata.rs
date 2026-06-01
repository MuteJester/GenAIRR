//! `PyAllele` + `PyRefDataConfig` — the Python-facing reference data.
//!
//! `PyRefDataConfig` is a *mutable builder* on the Python side: you
//! construct it for a chain type (`vj` or `vdj`) and grow each
//! segment pool by adding alleles. Once built, it's passed to a
//! runner (e.g. `run_vj_recombination`) which clones the inner
//! Rust `RefDataConfig` for the duration of the simulation.
//!
//! `PyAllele` is the corresponding read-only view returned by
//! pool accessors. F.3 keeps the Python surface narrow: enough
//! to drive recombination and inspect what was supplied, no more.

use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyBytes;

use crate::ir::Segment;
use crate::refdata::{
    Allele, AlleleId, AllelePool, ChainType, FunctionalStatus, RefDataConfig,
    RefDataCurationPolicy, RefDataIssueSeverity, RefDataValidationIssue, RefDataValidationMode,
    ReferenceIdentity,
};

/// A read-only view of one V / D / J allele.
#[pyclass(name = "Allele", module = "GenAIRR._engine", frozen)]
pub struct PyAllele {
    pub(crate) inner: Allele,
}

impl PyAllele {
    pub(crate) fn new(inner: Allele) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyAllele {
    /// IMGT allele identifier (e.g. `"IGHV1-2*01"`).
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    /// IMGT gene name (e.g. `"IGHV1-2"`).
    #[getter]
    fn gene(&self) -> &str {
        &self.inner.gene
    }

    /// Reference nucleotide sequence as ASCII bytes.
    fn seq<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new_bound(py, &self.inner.seq)
    }

    /// Length of the reference sequence in bases.
    fn __len__(&self) -> usize {
        self.inner.len() as usize
    }

    /// Segment role (`"V"`, `"D"`, or `"J"`).
    #[getter]
    fn segment(&self) -> &'static str {
        match self.inner.segment {
            Segment::V => "V",
            Segment::D => "D",
            Segment::J => "J",
            Segment::Np1 | Segment::Np2 => unreachable!("alleles are only V/D/J"),
        }
    }

    /// Position of the conserved anchor codon (V Cys / J W/F) in
    /// the allele coordinate system, or `None` for anchorless
    /// pseudogenes / partial alleles.
    #[getter]
    fn anchor(&self) -> Option<u16> {
        self.inner.anchor
    }

    /// IMGT-style functional classification carried by the cartridge,
    /// as a lowercase string (`"functional"` / `"orf"` / `"pseudogene"`
    /// / `"unknown"`), or `None` if the cartridge did not annotate this
    /// allele. Bundled `.pkl` catalogues currently leave this `None`.
    #[getter]
    fn functional_status(&self) -> Option<&'static str> {
        self.inner.functional_status.map(functional_status_to_str)
    }

    /// V-region substructure annotations as a list of
    /// ``(label, start, end)`` tuples. Empty for D/J alleles and
    /// for V alleles without IMGT subregion metadata. Coordinates
    /// are ungapped, allele-relative, half-open. Labels are the
    /// canonical IMGT uppercase strings (``"FWR1"``, ``"CDR1"``,
    /// ``"FWR2"``, ``"CDR2"``, ``"FWR3"``). See
    /// ``docs/v_region_substructure_audit.md``.
    #[getter]
    fn subregions(&self) -> Vec<(String, u16, u16)> {
        self.inner
            .subregions
            .iter()
            .map(|s| (s.label.as_str().to_string(), s.start, s.end))
            .collect()
    }

    fn __repr__(&self) -> String {
        match self.inner.anchor {
            Some(a) => format!(
                "<Allele {} {} {}bp anchor={}>",
                self.segment(),
                self.inner.name,
                self.inner.len(),
                a
            ),
            None => format!(
                "<Allele {} {} {}bp anchor=None>",
                self.segment(),
                self.inner.name,
                self.inner.len()
            ),
        }
    }
}

/// Mutable builder for the V/D/J reference pools used by every
/// recombination plan.
///
/// Construction:
/// ```python
/// cfg = RefDataConfig.vj()
/// cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
/// cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
/// ```
#[pyclass(name = "RefDataConfig", module = "GenAIRR._engine")]
pub struct PyRefDataConfig {
    pub(crate) inner: RefDataConfig,
}

impl PyRefDataConfig {
    /// Borrowed view of the inner Rust config — used by runners that
    /// need to thread it into a `PassRuntime`. Internal API only.
    pub(crate) fn inner(&self) -> &RefDataConfig {
        &self.inner
    }
}

#[pymethods]
impl PyRefDataConfig {
    /// Construct an empty refdata for the given chain type.
    /// `chain_type` is `"vj"` (light chain) or `"vdj"` (heavy chain),
    /// case-insensitive.
    #[new]
    fn new(chain_type: &str) -> PyResult<Self> {
        let ct = parse_chain_type(chain_type)?;
        Ok(Self {
            inner: RefDataConfig::empty(ct),
        })
    }

    /// Construct an empty VJ-chain (light) refdata.
    #[staticmethod]
    fn vj() -> Self {
        Self {
            inner: RefDataConfig::empty(ChainType::Vj),
        }
    }

    /// Construct an empty VDJ-chain (heavy) refdata.
    #[staticmethod]
    fn vdj() -> Self {
        Self {
            inner: RefDataConfig::empty(ChainType::Vdj),
        }
    }

    /// Chain type as a lowercase string (`"vj"` or `"vdj"`).
    #[getter]
    fn chain_type(&self) -> &'static str {
        match self.inner.chain_type {
            ChainType::Vj => "vj",
            ChainType::Vdj => "vdj",
        }
    }

    /// Whether the chain type carries a D segment (i.e. is heavy).
    fn has_d(&self) -> bool {
        self.inner.chain_type.has_d()
    }

    // ──────────────────────────────────────────────────────────────
    // Reference rules — the programmable interpretation layer.
    // ──────────────────────────────────────────────────────────────

    /// Configure the V anchor rule.
    ///
    /// ``expected_aa`` is a list of single-character amino-acid
    /// strings the V anchor codon may translate to (default
    /// ``['C']``). ``required`` controls whether anchorless V alleles
    /// emit a `MissingAnchor` issue. ``missing_severity`` and
    /// ``mismatch_severity`` are ``"fatal"`` or ``"curatable"`` —
    /// they tag the issues this rule emits, controlling whether
    /// `AllowCuratable` mode lets them through.
    #[pyo3(signature = (*, expected_aa, required = true, missing_severity = "curatable", mismatch_severity = "curatable"))]
    fn set_v_anchor_rule(
        &mut self,
        expected_aa: Vec<String>,
        required: bool,
        missing_severity: &str,
        mismatch_severity: &str,
    ) -> PyResult<()> {
        self.inner.rules.v_anchor = build_anchor_rule(
            "V",
            expected_aa,
            required,
            missing_severity,
            mismatch_severity,
        )?;
        Ok(())
    }

    /// Configure the J anchor rule. Same options as
    /// [`Self::set_v_anchor_rule`]; default expected AAs are
    /// ``['W', 'F']`` (lenient — pre-load default).
    #[pyo3(signature = (*, expected_aa, required = true, missing_severity = "curatable", mismatch_severity = "curatable"))]
    fn set_j_anchor_rule(
        &mut self,
        expected_aa: Vec<String>,
        required: bool,
        missing_severity: &str,
        mismatch_severity: &str,
    ) -> PyResult<()> {
        self.inner.rules.j_anchor = build_anchor_rule(
            "J",
            expected_aa,
            required,
            missing_severity,
            mismatch_severity,
        )?;
        Ok(())
    }

    /// Configure the allowed nucleotide alphabet for allele
    /// sequences. ``bases`` is a list of single-character strings;
    /// case-folded comparison is used. Default ``["A", "C", "G",
    /// "T", "N"]``.
    fn set_allowed_bases(&mut self, bases: Vec<String>) -> PyResult<()> {
        let mut allowed: Vec<u8> = Vec::with_capacity(bases.len());
        for b in &bases {
            let bytes = b.as_bytes();
            if bytes.len() != 1 {
                return Err(PyValueError::new_err(format!(
                    "alphabet entries must be single characters, got {b:?}"
                )));
            }
            allowed.push(bytes[0].to_ascii_uppercase());
        }
        self.inner.rules.alphabet =
            crate::refdata::ReferenceAlphabet { allowed };
        Ok(())
    }

    /// Read the V anchor rule as a dict: ``{"expected_aa": [...],
    /// "required": bool, "missing_severity": str, "mismatch_severity":
    /// str}``.
    fn v_anchor_rule<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        anchor_rule_to_pydict(py, &self.inner.rules.v_anchor)
    }

    /// Read the J anchor rule. Same shape as [`Self::v_anchor_rule`].
    fn j_anchor_rule<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        anchor_rule_to_pydict(py, &self.inner.rules.j_anchor)
    }

    /// Read the allowed alphabet as a list of single-character strings
    /// (upper-case canonical form).
    fn allowed_bases(&self) -> Vec<String> {
        self.inner
            .rules
            .alphabet
            .allowed
            .iter()
            .map(|&b| (b as char).to_string())
            .collect()
    }

    /// Deterministic content hash of this refdata. Two configs with
    /// identical catalogues AND identical [`ReferenceRules`] AND
    /// identical [`ReferenceIdentity`] produce the same hash; any
    /// difference produces a different hash. Format:
    /// ``"sha256:{hex}"``.
    fn content_hash(&self) -> String {
        crate::trace_file::refdata_content_hash(&self.inner)
    }

    // ──────────────────────────────────────────────────────────────
    // Reference identity — self-describing cartridge metadata.
    // ──────────────────────────────────────────────────────────────

    /// Set one or more identity fields. Every field is keyword-only
    /// and optional; passing ``None`` (or omitting) leaves the
    /// existing value unchanged. Pass an empty string to clear a
    /// field (rare — useful for tests). See [`ReferenceIdentity`].
    #[pyo3(signature = (
        *,
        species = None,
        locus = None,
        reference_set = None,
        name = None,
        source = None,
    ))]
    fn set_identity(
        &mut self,
        species: Option<String>,
        locus: Option<String>,
        reference_set: Option<String>,
        name: Option<String>,
        source: Option<String>,
    ) {
        let id = &mut self.inner.identity;
        if let Some(v) = species {
            id.species = into_some_or_none(v);
        }
        if let Some(v) = locus {
            id.locus = into_some_or_none(v);
        }
        if let Some(v) = reference_set {
            id.reference_set = into_some_or_none(v);
        }
        if let Some(v) = name {
            id.name = into_some_or_none(v);
        }
        if let Some(v) = source {
            id.source = into_some_or_none(v);
        }
    }

    /// Return a curated copy of this cartridge under ``policy``.
    ///
    /// ``policy`` is one of:
    ///   - ``"raw"`` — identity policy; returns an unmodified clone.
    ///   - ``"functional_anchors_only"`` — drop V/J alleles whose
    ///     anchor doesn't satisfy the active anchor rule (missing,
    ///     out of bounds, or codon AA outside `expected_amino_acids`).
    ///     D and C pools pass through unchanged.
    ///   - ``"functional_status"`` — filter V/D/J pools by IMGT
    ///     functional status. ``allowed`` is the list of statuses to
    ///     keep (case-insensitive strings, default
    ///     ``["functional"]``). ``keep_unannotated`` (default
    ///     ``True``) controls whether alleles with no annotation
    ///     survive the filter — bundled ``.pkl`` catalogues currently
    ///     leave the field unannotated, so the default preserves
    ///     them. C pool passes through unchanged.
    ///
    /// Curation never fixes structural problems (duplicate names,
    /// invalid bytes, locus/chain mismatch) — those continue to
    /// surface from ``validate()`` on the curated cartridge. The
    /// curated cartridge's ``identity.source`` is tagged
    /// ``|curated:<policy>`` so its content hash differs from the
    /// raw and trace files attribute outputs correctly.
    #[pyo3(signature = (policy, *, allowed=None, keep_unannotated=true))]
    fn curated(
        &self,
        policy: &str,
        allowed: Option<Vec<String>>,
        keep_unannotated: bool,
    ) -> PyResult<Self> {
        let policy = build_curation_policy(policy, allowed, keep_unannotated)?;
        Ok(Self {
            inner: self.inner.curated(policy),
        })
    }

    /// In-place sibling of [`Self::curated`] — applies the policy
    /// and mutates this cartridge. Convenient for fluent builder
    /// code that doesn't want to thread a rebound name through.
    #[pyo3(signature = (policy, *, allowed=None, keep_unannotated=true))]
    fn curate(
        &mut self,
        policy: &str,
        allowed: Option<Vec<String>>,
        keep_unannotated: bool,
    ) -> PyResult<()> {
        let policy = build_curation_policy(policy, allowed, keep_unannotated)?;
        self.inner = self.inner.curated(policy);
        Ok(())
    }

    /// Read the identity as a dict. Missing fields are present with
    /// value ``None``; consumers can use ``dict.get(...)`` without
    /// guarding.
    fn identity<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        use pyo3::types::PyDict;
        let d = PyDict::new_bound(py);
        let id = &self.inner.identity;
        d.set_item("species", id.species.as_deref())?;
        d.set_item("locus", id.locus.as_deref())?;
        d.set_item("reference_set", id.reference_set.as_deref())?;
        d.set_item("name", id.name.as_deref())?;
        d.set_item("source", id.source.as_deref())?;
        Ok(d)
    }

    /// Append a V allele. Returns the new allele id (`u32`).
    ///
    /// ``functional_status`` is an optional IMGT-style classification
    /// string — ``"functional"`` / ``"orf"`` / ``"pseudogene"`` /
    /// ``"unknown"`` — case-insensitive. ``None`` means the cartridge
    /// did not provide this annotation (the default for legacy
    /// bundled catalogues). Invalid strings raise ``ValueError``.
    #[pyo3(signature = (name, gene, seq, *, anchor=None, functional_status=None, subregions=None))]
    fn add_v_allele(
        &mut self,
        name: &str,
        gene: &str,
        seq: &[u8],
        anchor: Option<u16>,
        functional_status: Option<&str>,
        subregions: Option<Vec<(String, u16, u16)>>,
    ) -> PyResult<u32> {
        let status = parse_functional_status(functional_status)?;
        let subregions =
            parse_subregions(subregions.as_deref(), seq.len() as u16)?;
        Ok(push_allele(
            &mut self.inner.v_pool,
            Segment::V,
            name,
            gene,
            seq,
            anchor,
            status,
            subregions,
        ))
    }

    /// Append a D allele. Returns the new allele id (`u32`). Panics
    /// if called on a VJ chain — D pools are only valid for heavy chains.
    ///
    /// See [`Self::add_v_allele`] for ``functional_status`` semantics.
    #[pyo3(signature = (name, gene, seq, *, anchor=None, functional_status=None))]
    fn add_d_allele(
        &mut self,
        name: &str,
        gene: &str,
        seq: &[u8],
        anchor: Option<u16>,
        functional_status: Option<&str>,
    ) -> PyResult<u32> {
        if !self.inner.chain_type.has_d() {
            return Err(PyValueError::new_err(
                "cannot add D allele to a VJ-chain RefDataConfig",
            ));
        }
        let status = parse_functional_status(functional_status)?;
        Ok(push_allele(
            &mut self.inner.d_pool,
            Segment::D,
            name,
            gene,
            seq,
            anchor,
            status,
            Vec::new(), // D alleles never carry V-region subregions
        ))
    }

    /// Append a J allele. Returns the new allele id (`u32`).
    ///
    /// See [`Self::add_v_allele`] for ``functional_status`` semantics.
    #[pyo3(signature = (name, gene, seq, *, anchor=None, functional_status=None))]
    fn add_j_allele(
        &mut self,
        name: &str,
        gene: &str,
        seq: &[u8],
        anchor: Option<u16>,
        functional_status: Option<&str>,
    ) -> PyResult<u32> {
        let status = parse_functional_status(functional_status)?;
        Ok(push_allele(
            &mut self.inner.j_pool,
            Segment::J,
            name,
            gene,
            seq,
            anchor,
            status,
            Vec::new(), // J alleles never carry V-region subregions
        ))
    }

    /// Number of V alleles in the pool.
    fn v_pool_size(&self) -> usize {
        self.inner.v_pool.len()
    }

    /// Number of D alleles in the pool. Always `0` for VJ chains.
    fn d_pool_size(&self) -> usize {
        self.inner.d_pool.len()
    }

    /// Number of J alleles in the pool.
    fn j_pool_size(&self) -> usize {
        self.inner.j_pool.len()
    }

    /// V allele at the given id, or raises `IndexError`.
    fn v_allele(&self, id: u32) -> PyResult<PyAllele> {
        self.inner
            .v_pool
            .get(AlleleId::new(id))
            .cloned()
            .map(PyAllele::new)
            .ok_or_else(|| PyIndexError::new_err(format!("V allele id {} out of range", id)))
    }

    /// D allele at the given id, or raises `IndexError`.
    fn d_allele(&self, id: u32) -> PyResult<PyAllele> {
        self.inner
            .d_pool
            .get(AlleleId::new(id))
            .cloned()
            .map(PyAllele::new)
            .ok_or_else(|| PyIndexError::new_err(format!("D allele id {} out of range", id)))
    }

    /// J allele at the given id, or raises `IndexError`.
    fn j_allele(&self, id: u32) -> PyResult<PyAllele> {
        self.inner
            .j_pool
            .get(AlleleId::new(id))
            .cloned()
            .map(PyAllele::new)
            .ok_or_else(|| PyIndexError::new_err(format!("J allele id {} out of range", id)))
    }

    fn __repr__(&self) -> String {
        format!(
            "<RefDataConfig {} V={} D={} J={}>",
            self.chain_type(),
            self.v_pool_size(),
            self.d_pool_size(),
            self.j_pool_size(),
        )
    }

    /// Run the structural validator on this refdata. Returns a list
    /// of issue dicts; empty list means the refdata passes every
    /// gate.
    ///
    /// Each issue dict carries:
    ///   - `kind`: variant name (e.g. `"EmptyRequiredPool"`,
    ///     `"VAnchorNotCys"`, `"JAnchorUnexpectedAa"`).
    ///   - `severity`: `"fatal"` or `"curatable"` — fatal issues
    ///     gate compile in every mode; curatable issues reflect
    ///     pseudogene/ORF alleles and gate compile only in the
    ///     default strict mode.
    ///   - plus variant-specific fields (`segment`, `allele_id`,
    ///     `pos`, `byte`, `anchor`, `len`, `codon`, `aa`,
    ///     `expected`, `name`).
    ///
    /// Read-only; never raises. To force errors instead, call
    /// [`Self::validate_strict`].
    fn validate<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Vec<Bound<'py, pyo3::types::PyDict>>> {
        let issues = self.inner.validate();
        issues
            .into_iter()
            .map(|issue| issue_to_pydict(py, &issue))
            .collect()
    }

    /// Strict variant: raises `ValueError` with a multi-line summary
    /// of every issue, tagged by severity. Returns `None` on success.
    fn validate_strict(&self) -> PyResult<()> {
        match self.inner.validate_strict() {
            Ok(()) => Ok(()),
            Err(errs) => Err(PyValueError::new_err(format!("{errs}"))),
        }
    }

    /// Mode-aware validation gate. ``mode`` is ``"strict"`` (default)
    /// or ``"allow_curatable"``. Strict rejects every validation
    /// issue. AllowCuratable accepts Curatable issues
    /// (pseudogene-shape anomalies) but still rejects Fatal ones
    /// (empty pool, duplicate name, invalid byte, anchor out of
    /// bounds). Raises ``ValueError`` on failure with the full
    /// severity-tagged issue list embedded in the message; returns
    /// ``None`` on success.
    #[pyo3(signature = (*, mode = "strict"))]
    fn validate_with_mode(&self, mode: &str) -> PyResult<()> {
        let mode = match mode {
            "strict" => RefDataValidationMode::Strict,
            "allow_curatable" => RefDataValidationMode::AllowCuratable,
            other => {
                return Err(PyValueError::new_err(format!(
                    "unknown validation mode {other:?}; expected 'strict' or 'allow_curatable'",
                )));
            }
        };
        match self.inner.validate_with_mode(mode) {
            Ok(()) => Ok(()),
            Err(errs) => Err(PyValueError::new_err(format!("{errs}"))),
        }
    }
}

/// Identity field policy: an empty string clears the slot, any
/// other value sets it. This lets ``set_identity(species="")``
/// reset a field while ``set_identity(species=None)`` leaves it
/// untouched.
fn into_some_or_none(s: String) -> Option<String> {
    if s.is_empty() {
        None
    } else {
        Some(s)
    }
}

#[allow(dead_code)]
fn _identity_marker(_: ReferenceIdentity) {}

/// Parse a Python ``functional_status`` argument into the Rust
/// enum. Accepted strings (case-insensitive):
/// - `"functional"`, `"f"` → [`FunctionalStatus::Functional`]
/// - `"orf"` → [`FunctionalStatus::Orf`]
/// - `"pseudogene"`, `"p"` → [`FunctionalStatus::Pseudogene`]
/// - `"unknown"` → [`FunctionalStatus::Unknown`]
///
/// `None` (or omitted at the call site) yields `Ok(None)`, meaning
/// the cartridge did not provide an annotation. Any other value
/// raises `ValueError`.
pub(crate) fn parse_functional_status(
    s: Option<&str>,
) -> PyResult<Option<FunctionalStatus>> {
    let Some(raw) = s else { return Ok(None) };
    match raw.to_ascii_lowercase().as_str() {
        "functional" | "f" => Ok(Some(FunctionalStatus::Functional)),
        "orf" => Ok(Some(FunctionalStatus::Orf)),
        "pseudogene" | "p" => Ok(Some(FunctionalStatus::Pseudogene)),
        "unknown" => Ok(Some(FunctionalStatus::Unknown)),
        other => Err(PyValueError::new_err(format!(
            "unknown functional_status {other:?}; expected one of \
             'functional' / 'orf' / 'pseudogene' / 'unknown' \
             (case-insensitive; 'F' / 'P' aliases also accepted)",
        ))),
    }
}

fn functional_status_to_str(s: FunctionalStatus) -> &'static str {
    s.as_str()
}

/// Translate the Python ``[(label, start, end), …]`` shape into a
/// validated ``Vec<VSubregion>``. Empty input → empty Vec; ``None`` →
/// empty Vec (the bridge default). Validation rules (audit §4):
///
/// - Each label must be one of ``"FWR1"``, ``"CDR1"``, ``"FWR2"``,
///   ``"CDR2"``, ``"FWR3"`` (case-sensitive — Python side normalises
///   to uppercase before calling).
/// - Each interval must satisfy ``start < end`` and
///   ``end <= seq_len``.
/// - No duplicate labels.
/// - No overlapping intervals.
///
/// Failure raises ``ValueError`` with a diagnostic naming the
/// offending label / interval. Designed so a custom cartridge that
/// hand-attaches a malformed annotation fails loudly at the bridge.
fn parse_subregions(
    entries: Option<&[(String, u16, u16)]>,
    seq_len: u16,
) -> PyResult<Vec<crate::refdata::VSubregion>> {
    let Some(entries) = entries else {
        return Ok(Vec::new());
    };
    if entries.is_empty() {
        return Ok(Vec::new());
    }
    let mut seen: std::collections::HashSet<crate::refdata::VSubregionLabel> =
        std::collections::HashSet::new();
    let mut out: Vec<crate::refdata::VSubregion> = Vec::with_capacity(entries.len());
    for (raw_label, start, end) in entries {
        let label = crate::refdata::VSubregionLabel::from_str(raw_label.as_str())
            .ok_or_else(|| {
                PyValueError::new_err(format!(
                    "subregions: unknown label {:?}; expected one of \
                     'FWR1' / 'CDR1' / 'FWR2' / 'CDR2' / 'FWR3' \
                     (case-sensitive)",
                    raw_label,
                ))
            })?;
        if !seen.insert(label) {
            return Err(PyValueError::new_err(format!(
                "subregions: duplicate label {:?}",
                label.as_str()
            )));
        }
        if start >= end {
            return Err(PyValueError::new_err(format!(
                "subregions[{:?}]: start ({}) must be strictly less than end ({})",
                label.as_str(),
                start,
                end,
            )));
        }
        if *end > seq_len {
            return Err(PyValueError::new_err(format!(
                "subregions[{:?}]: end ({}) exceeds allele length ({})",
                label.as_str(),
                end,
                seq_len,
            )));
        }
        out.push(crate::refdata::VSubregion {
            label,
            start: *start,
            end: *end,
        });
    }
    // Overlap check: sort by start, walk pairs.
    let mut sorted = out.clone();
    sorted.sort_by_key(|s| s.start);
    for pair in sorted.windows(2) {
        if pair[0].end > pair[1].start {
            return Err(PyValueError::new_err(format!(
                "subregions: intervals {:?} [{}..{}) and {:?} [{}..{}) overlap",
                pair[0].label.as_str(),
                pair[0].start,
                pair[0].end,
                pair[1].label.as_str(),
                pair[1].start,
                pair[1].end,
            )));
        }
    }
    Ok(out)
}

fn parse_curation_policy(name: &str) -> PyResult<RefDataCurationPolicy> {
    match name {
        "raw" => Ok(RefDataCurationPolicy::Raw),
        "functional_anchors_only" => Ok(RefDataCurationPolicy::FunctionalAnchorsOnly),
        other => Err(PyValueError::new_err(format!(
            "unknown curation policy {other:?}; expected 'raw', \
             'functional_anchors_only', or 'functional_status'",
        ))),
    }
}

/// Build a [`RefDataCurationPolicy`] from a Python policy name plus
/// the optional kwargs used by [`RefDataCurationPolicy::FunctionalStatus`].
///
/// - `"raw"` and `"functional_anchors_only"` ignore the kwargs.
/// - `"functional_status"` uses `allowed` (default `["functional"]`)
///   and `keep_unannotated` (default `true`). An empty `allowed` is
///   accepted (it filters every annotated allele); a non-empty list
///   with an unknown string raises `ValueError` via
///   [`parse_functional_status`].
fn build_curation_policy(
    name: &str,
    allowed: Option<Vec<String>>,
    keep_unannotated: bool,
) -> PyResult<RefDataCurationPolicy> {
    match name {
        "raw" | "functional_anchors_only" => {
            if allowed.is_some() {
                return Err(PyValueError::new_err(format!(
                    "curation policy {name:?} does not accept an 'allowed' \
                     argument; remove it or switch to 'functional_status'",
                )));
            }
            parse_curation_policy(name)
        }
        "functional_status" => {
            let raw = allowed.unwrap_or_else(|| vec!["functional".to_string()]);
            let mut parsed = Vec::with_capacity(raw.len());
            for s in raw {
                let normalised = parse_functional_status(Some(s.as_str()))?
                    .expect("Some(string) produces Some(status)");
                if !parsed.contains(&normalised) {
                    parsed.push(normalised);
                }
            }
            Ok(RefDataCurationPolicy::FunctionalStatus {
                allowed: parsed,
                keep_unannotated,
            })
        }
        other => Err(PyValueError::new_err(format!(
            "unknown curation policy {other:?}; expected 'raw', \
             'functional_anchors_only', or 'functional_status'",
        ))),
    }
}

fn parse_severity(name: &str, field: &str) -> PyResult<RefDataIssueSeverity> {
    match name {
        "fatal" => Ok(RefDataIssueSeverity::Fatal),
        "curatable" => Ok(RefDataIssueSeverity::Curatable),
        other => Err(PyValueError::new_err(format!(
            "{field}={other:?} must be 'fatal' or 'curatable'"
        ))),
    }
}

fn parse_expected_aa(segment: &str, raw: &[String]) -> PyResult<Vec<char>> {
    let mut out = Vec::with_capacity(raw.len());
    for s in raw {
        let mut chars = s.chars();
        match (chars.next(), chars.next()) {
            (Some(c), None) => out.push(c.to_ascii_uppercase()),
            _ => {
                return Err(PyValueError::new_err(format!(
                    "{segment} anchor expected_aa entries must be single \
                     amino-acid letters, got {s:?}"
                )));
            }
        }
    }
    if out.is_empty() {
        return Err(PyValueError::new_err(format!(
            "{segment} anchor expected_aa must be a non-empty list"
        )));
    }
    Ok(out)
}

fn build_anchor_rule(
    segment: &str,
    expected_aa: Vec<String>,
    required: bool,
    missing_severity: &str,
    mismatch_severity: &str,
) -> PyResult<crate::refdata::AnchorRule> {
    Ok(crate::refdata::AnchorRule {
        required,
        expected_amino_acids: parse_expected_aa(segment, &expected_aa)?,
        missing_severity: parse_severity(missing_severity, "missing_severity")?,
        mismatch_severity: parse_severity(mismatch_severity, "mismatch_severity")?,
    })
}

fn anchor_rule_to_pydict<'py>(
    py: Python<'py>,
    rule: &crate::refdata::AnchorRule,
) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    use pyo3::types::PyDict;
    let d = PyDict::new_bound(py);
    let expected: Vec<String> = rule
        .expected_amino_acids
        .iter()
        .map(|c| c.to_string())
        .collect();
    d.set_item("expected_aa", expected)?;
    d.set_item("required", rule.required)?;
    d.set_item("missing_severity", severity_str(rule.missing_severity))?;
    d.set_item("mismatch_severity", severity_str(rule.mismatch_severity))?;
    Ok(d)
}

fn severity_str(s: RefDataIssueSeverity) -> &'static str {
    match s {
        RefDataIssueSeverity::Fatal => "fatal",
        RefDataIssueSeverity::Curatable => "curatable",
    }
}

fn segment_str(s: Segment) -> &'static str {
    match s {
        Segment::V => "V",
        Segment::D => "D",
        Segment::J => "J",
        Segment::Np1 => "Np1",
        Segment::Np2 => "Np2",
    }
}

fn issue_to_pydict<'py>(
    py: Python<'py>,
    issue: &RefDataValidationIssue,
) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    use pyo3::types::PyDict;
    let d = PyDict::new_bound(py);
    let severity = match issue.severity() {
        RefDataIssueSeverity::Fatal => "fatal",
        RefDataIssueSeverity::Curatable => "curatable",
    };
    d.set_item("severity", severity)?;
    match issue {
        RefDataValidationIssue::EmptyRequiredPool { segment } => {
            d.set_item("kind", "EmptyRequiredPool")?;
            d.set_item("segment", segment_str(*segment))?;
        }
        RefDataValidationIssue::DuplicateAlleleName { segment, name } => {
            d.set_item("kind", "DuplicateAlleleName")?;
            d.set_item("segment", segment_str(*segment))?;
            d.set_item("name", name)?;
        }
        RefDataValidationIssue::InvalidAlleleByte {
            segment,
            allele_id,
            pos,
            byte,
        } => {
            d.set_item("kind", "InvalidAlleleByte")?;
            d.set_item("segment", segment_str(*segment))?;
            d.set_item("allele_id", allele_id.index())?;
            d.set_item("pos", pos)?;
            d.set_item("byte", *byte as u32)?;
        }
        RefDataValidationIssue::AnchorOutOfBounds {
            segment,
            allele_id,
            anchor,
            len,
        } => {
            d.set_item("kind", "AnchorOutOfBounds")?;
            d.set_item("segment", segment_str(*segment))?;
            d.set_item("allele_id", allele_id.index())?;
            d.set_item("anchor", *anchor as u32)?;
            d.set_item("len", len)?;
        }
        RefDataValidationIssue::VAnchorNotCys {
            allele_id,
            codon,
            aa,
            severity: _,
        } => {
            d.set_item("kind", "VAnchorNotCys")?;
            d.set_item("allele_id", allele_id.index())?;
            d.set_item("codon", PyBytes::new_bound(py, codon))?;
            d.set_item("aa", aa.to_string())?;
        }
        RefDataValidationIssue::JAnchorUnexpectedAa {
            allele_id,
            codon,
            aa,
            expected,
            severity: _,
        } => {
            d.set_item("kind", "JAnchorUnexpectedAa")?;
            d.set_item("allele_id", allele_id.index())?;
            d.set_item("codon", PyBytes::new_bound(py, codon))?;
            d.set_item("aa", aa.to_string())?;
            d.set_item("expected", expected.iter().map(|c| c.to_string()).collect::<Vec<_>>())?;
        }
        RefDataValidationIssue::MissingAnchor {
            segment,
            allele_id,
            severity: _,
        } => {
            d.set_item("kind", "MissingAnchor")?;
            d.set_item("segment", segment_str(*segment))?;
            d.set_item("allele_id", allele_id.index())?;
        }
        RefDataValidationIssue::LocusChainTypeMismatch { locus, chain_type } => {
            d.set_item("kind", "LocusChainTypeMismatch")?;
            d.set_item("locus", locus)?;
            d.set_item(
                "chain_type",
                match chain_type {
                    crate::refdata::ChainType::Vj => "vj",
                    crate::refdata::ChainType::Vdj => "vdj",
                },
            )?;
        }
    }
    Ok(d.into())
}

fn parse_chain_type(s: &str) -> PyResult<ChainType> {
    match s.to_ascii_lowercase().as_str() {
        "vj" => Ok(ChainType::Vj),
        "vdj" => Ok(ChainType::Vdj),
        other => Err(PyValueError::new_err(format!(
            "unknown chain type {:?}; expected 'vj' or 'vdj'",
            other
        ))),
    }
}

fn push_allele(
    pool: &mut AllelePool,
    segment: Segment,
    name: &str,
    gene: &str,
    seq: &[u8],
    anchor: Option<u16>,
    functional_status: Option<FunctionalStatus>,
    subregions: Vec<crate::refdata::VSubregion>,
) -> u32 {
    let id = pool.push(Allele {
        name: name.to_string(),
        gene: gene.to_string(),
        seq: seq.to_vec(),
        segment,
        anchor,
        functional_status,
        subregions,
    });
    id.index()
}
