//! `PyOutcome` — read-only view of [`crate::pass::Outcome`].

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::address::PrimeEnd;
use crate::airr_record::{
    build_airr_record, AirrRecord, AlleleOrderReason, ProductiveDecidedBy, RecordValidationIssue,
};
use crate::ir::Segment;
use crate::pass::Outcome;

use super::event::PyEventRecord;
use super::refdata::PyRefDataConfig;
use super::simulation::PySimulation;
use super::trace::PyTrace;

/// The result of executing a `PassPlan`: the IR revision history,
/// the per-revision pass names, the addressed-choice trace, and the
/// committed event ledger.
#[pyclass(name = "Outcome", module = "GenAIRR._engine", frozen)]
pub struct PyOutcome {
    pub(crate) inner: Outcome,
}

impl PyOutcome {
    pub(crate) fn new(inner: Outcome) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyOutcome {
    /// Number of revisions (= initial + one per executed pass).
    fn revision_count(&self) -> usize {
        self.inner.revisions.len()
    }

    /// Revision at the given index. `0` is the initial simulation
    /// (input to the first pass); `i+1` is the output of the `i`-th
    /// pass. Raises `IndexError` when `index` is out of range.
    fn revision(&self, index: usize) -> PyResult<PySimulation> {
        self.inner
            .revisions
            .get(index)
            .cloned()
            .map(PySimulation::new)
            .ok_or_else(|| {
                PyIndexError::new_err(format!(
                    "revision index {} out of range (have {})",
                    index,
                    self.inner.revisions.len(),
                ))
            })
    }

    /// The final IR revision after every pass has run.
    fn final_simulation(&self) -> PySimulation {
        PySimulation::new(
            self.inner
                .revisions
                .last()
                .cloned()
                .expect("Outcome guarantees at least one revision"),
        )
    }

    /// Pass names in execution order. `pass_names()[i]` is the name
    /// of the pass that produced `revision(i + 1)`. Length equals the
    /// number of executed passes.
    fn pass_names(&self) -> Vec<String> {
        self.inner.pass_names.clone()
    }

    /// First revision produced by the pass with the given `name`, or
    /// `None` if no pass with that name ran.
    fn revision_after(&self, name: &str) -> Option<PySimulation> {
        self.inner
            .revision_after(name)
            .cloned()
            .map(PySimulation::new)
    }

    /// The addressed-choice trace recorded during the run.
    fn trace(&self) -> PyTrace {
        PyTrace::new(self.inner.trace.clone())
    }

    /// Number of committed event records.
    fn event_count(&self) -> usize {
        self.inner.events.len()
    }

    /// The committed event ledger, in execution order.
    fn events(&self) -> Vec<PyEventRecord> {
        self.inner
            .events
            .iter()
            .cloned()
            .map(PyEventRecord::new)
            .collect()
    }

    /// Build a fully-populated AIRR Rearrangement record dict from
    /// this outcome and the reference data it ran against.
    ///
    /// Returns a Python `dict` with the AIRR Rearrangement field
    /// shape — ~50 columns, 0-based half-open coordinates. The
    /// Python `SimulationResult` layer adds the AIRR-strict
    /// 1-based-inclusive transformation at TSV/CSV/DataFrame export
    /// time.
    ///
    /// `sequence_id` defaults to an empty string when omitted.
    #[pyo3(signature = (refdata, *, sequence_id = ""))]
    fn airr_record<'py>(
        &self,
        py: Python<'py>,
        refdata: &PyRefDataConfig,
        sequence_id: &str,
    ) -> PyResult<Bound<'py, PyDict>> {
        let rec = build_airr_record(&self.inner, refdata.inner(), sequence_id);
        airr_record_to_pydict(py, &rec)
    }

    /// Build the AIRR record and validate it against this outcome.
    /// Returns a list of issue dicts. Empty list means the record
    /// passed every check.
    ///
    /// Each issue is a dict with at least `kind` (stable identifier
    /// matching the Rust enum variant name). Variants that involve
    /// a segment include `segment` (`"V"`, `"D"`, or `"J"`). Variants
    /// that compare engine-reported vs validator-recomputed values
    /// include `reported` and `expected`. Variant-specific extras
    /// land under `details`.
    ///
    /// See `docs/airr_record_validator.md` for the check catalogue.
    #[pyo3(signature = (refdata, *, sequence_id = ""))]
    fn validate_record<'py>(
        &self,
        py: Python<'py>,
        refdata: &PyRefDataConfig,
        sequence_id: &str,
    ) -> PyResult<Vec<Bound<'py, PyDict>>> {
        use crate::airr_record::validate_airr_record;
        let rec = build_airr_record(&self.inner, refdata.inner(), sequence_id);
        let issues = validate_airr_record(&rec, &self.inner, refdata.inner());
        issues
            .into_iter()
            .map(|i| issue_to_pydict(py, i))
            .collect()
    }

    /// **Internal live-call cache correctness check.**
    ///
    /// Compare the cached `SegmentLiveCall` on the final simulation
    /// against a fresh from-scratch recomputation, per V/D/J. This is
    /// the engine-side guard: "does the cached state that feeds
    /// projection equal what a from-scratch walk would produce?"
    ///
    /// Returns a list of per-segment parity dicts:
    ///
    /// ```text
    /// {
    ///   "segment": "V" / "D" / "J",
    ///   "tie_set_matches": bool,
    ///   "cached_tie_set": [allele_id, ...],
    ///   "fresh_tie_set": [allele_id, ...],
    ///   "cached_present": bool,
    ///   "fresh_present": bool,
    ///   "hypothesis_bounds_match": bool | None,
    ///   "cached_hypothesis": {seq_start, seq_end, ref_start, ref_end} | None,
    ///   "fresh_hypothesis": {seq_start, seq_end, ref_start, ref_end} | None,
    /// }
    /// ```
    ///
    /// Use this in tests as an explicit cache-equivalence guard:
    ///
    /// ```text
    /// for p in outcome.check_live_call_cache_parity(refdata):
    ///     assert p["tie_set_matches"], p
    /// ```
    ///
    /// **Companion check** — for downstream AIRR-record correctness,
    /// see `SimulationResult.validate_records(refdata)` (the
    /// user-facing postcondition over projected output).
    ///
    /// **Troubleshooting rule** — if a CI run has both this parity
    /// harness AND `validate_records` failing on the same batch,
    /// fix the parity divergence FIRST: a stale cache leaks into
    /// projection and produces spurious validator failures
    /// downstream. Once parity is green, rerun the validator;
    /// remaining failures point at a real projection-layer bug.
    ///
    /// See `docs/airr_record_validator.md` §5.2 for the bug-class
    /// history that motivated this harness.
    fn check_live_call_cache_parity<'py>(
        &self,
        py: Python<'py>,
        refdata: &PyRefDataConfig,
    ) -> PyResult<Vec<Bound<'py, PyDict>>> {
        use crate::live_call::check_segment_calls_parity;
        let sim = self.inner.final_simulation();
        let results = check_segment_calls_parity(sim, refdata.inner());
        results
            .into_iter()
            .map(|r| parity_to_pydict(py, r))
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "<Outcome revisions={} passes={} trace_len={} events={}>",
            self.inner.revisions.len(),
            self.inner.pass_names.len(),
            self.inner.trace.len(),
            self.inner.events.len(),
        )
    }
}

/// Convert an `AirrRecord` into a Python `dict` matching the field
/// names + types the Python `_airr_record` builder used to emit.
///
/// `Option<i64>` becomes `int | None`, `Option<f64>` becomes `float
/// | None`, `Option<bool>` becomes `bool | None`. Strings and the
/// few `bool`/`i64`/`f64` non-optional fields go through directly.
fn airr_record_to_pydict<'py>(py: Python<'py>, rec: &AirrRecord) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new_bound(py);

    // AIRR metadata
    dict.set_item("sequence_id", &rec.sequence_id)?;
    dict.set_item("sequence", &rec.sequence)?;
    dict.set_item("sequence_aa", &rec.sequence_aa)?;
    dict.set_item("sequence_alignment", &rec.sequence_alignment)?;
    dict.set_item("germline_alignment", &rec.germline_alignment)?;
    dict.set_item("germline_alignment_d_mask", &rec.germline_alignment_d_mask)?;
    dict.set_item("sequence_length", rec.sequence_length)?;
    dict.set_item("rev_comp", rec.rev_comp)?;
    dict.set_item("locus", &rec.locus)?;

    // V
    dict.set_item("v_call", &rec.v_call)?;
    dict.set_item("v_cigar", &rec.v_cigar)?;
    set_opt_f64(&dict, "v_score", rec.v_score)?;
    set_opt_f64(&dict, "v_identity", rec.v_identity)?;
    set_opt_f64(&dict, "v_support", rec.v_support)?;
    set_opt_i64(&dict, "v_sequence_start", rec.v_sequence_start)?;
    set_opt_i64(&dict, "v_sequence_end", rec.v_sequence_end)?;
    set_opt_i64(&dict, "v_alignment_start", rec.v_alignment_start)?;
    set_opt_i64(&dict, "v_alignment_end", rec.v_alignment_end)?;
    set_opt_i64(&dict, "v_germline_start", rec.v_germline_start)?;
    set_opt_i64(&dict, "v_germline_end", rec.v_germline_end)?;
    dict.set_item("v_trim_5", rec.v_trim_5)?;
    dict.set_item("v_trim_3", rec.v_trim_3)?;

    // D
    dict.set_item("d_call", &rec.d_call)?;
    dict.set_item("d_cigar", &rec.d_cigar)?;
    set_opt_f64(&dict, "d_score", rec.d_score)?;
    set_opt_f64(&dict, "d_identity", rec.d_identity)?;
    set_opt_f64(&dict, "d_support", rec.d_support)?;
    set_opt_i64(&dict, "d_sequence_start", rec.d_sequence_start)?;
    set_opt_i64(&dict, "d_sequence_end", rec.d_sequence_end)?;
    set_opt_i64(&dict, "d_alignment_start", rec.d_alignment_start)?;
    set_opt_i64(&dict, "d_alignment_end", rec.d_alignment_end)?;
    set_opt_i64(&dict, "d_germline_start", rec.d_germline_start)?;
    set_opt_i64(&dict, "d_germline_end", rec.d_germline_end)?;
    dict.set_item("d_trim_5", rec.d_trim_5)?;
    dict.set_item("d_trim_3", rec.d_trim_3)?;

    // J
    dict.set_item("j_call", &rec.j_call)?;
    dict.set_item("j_cigar", &rec.j_cigar)?;
    set_opt_f64(&dict, "j_score", rec.j_score)?;
    set_opt_f64(&dict, "j_identity", rec.j_identity)?;
    set_opt_f64(&dict, "j_support", rec.j_support)?;
    set_opt_i64(&dict, "j_sequence_start", rec.j_sequence_start)?;
    set_opt_i64(&dict, "j_sequence_end", rec.j_sequence_end)?;
    set_opt_i64(&dict, "j_alignment_start", rec.j_alignment_start)?;
    set_opt_i64(&dict, "j_alignment_end", rec.j_alignment_end)?;
    set_opt_i64(&dict, "j_germline_start", rec.j_germline_start)?;
    set_opt_i64(&dict, "j_germline_end", rec.j_germline_end)?;
    dict.set_item("j_trim_5", rec.j_trim_5)?;
    dict.set_item("j_trim_3", rec.j_trim_3)?;

    // C
    dict.set_item("c_call", &rec.c_call)?;

    // Junction
    dict.set_item("junction", &rec.junction)?;
    dict.set_item("junction_aa", &rec.junction_aa)?;
    set_opt_i64(&dict, "junction_start", rec.junction_start)?;
    set_opt_i64(&dict, "junction_end", rec.junction_end)?;
    set_opt_i64(&dict, "junction_length", rec.junction_length)?;

    // NP regions
    dict.set_item("np1", &rec.np1)?;
    dict.set_item("np1_aa", &rec.np1_aa)?;
    dict.set_item("np1_length", rec.np1_length)?;
    dict.set_item("np2", &rec.np2)?;
    dict.set_item("np2_aa", &rec.np2_aa)?;
    dict.set_item("np2_length", rec.np2_length)?;

    // Functionality
    set_opt_bool(&dict, "productive", rec.productive)?;
    set_opt_bool(&dict, "vj_in_frame", rec.vj_in_frame)?;
    set_opt_bool(&dict, "stop_codon", rec.stop_codon)?;

    // SHM + corruption
    dict.set_item("n_mutations", rec.n_mutations)?;
    dict.set_item("mutation_rate", rec.mutation_rate)?;
    dict.set_item("n_pcr_errors", rec.n_pcr_errors)?;
    dict.set_item("n_quality_errors", rec.n_quality_errors)?;
    dict.set_item("n_indels", rec.n_indels)?;
    // Per-segment indel counters (audit §6.2 fix). Count
    // `IndelInserted` + `IndelDeleted` events attributed to V/D/J
    // by the event's `segment` field. NP1 / NP2 events count
    // toward `n_indels` but not these per-segment counters.
    dict.set_item("n_v_indels", rec.n_v_indels)?;
    dict.set_item("n_d_indels", rec.n_d_indels)?;
    dict.set_item("n_j_indels", rec.n_j_indels)?;
    // Observation-stage end-loss / primer-trim amounts (audit
    // §6.1 fix). Distinct from recombination-stage v_trim_*/
    // j_trim_*; defaults to 0 when no end-loss pass ran.
    dict.set_item("end_loss_5_length", rec.end_loss_5_length)?;
    dict.set_item("end_loss_3_length", rec.end_loss_3_length)?;
    dict.set_item("is_contaminant", rec.is_contaminant)?;

    Ok(dict)
}

#[inline]
fn set_opt_i64(dict: &Bound<'_, PyDict>, key: &str, val: Option<i64>) -> PyResult<()> {
    match val {
        Some(v) => dict.set_item(key, v),
        None => dict.set_item(key, dict.py().None()),
    }
}

#[inline]
fn set_opt_f64(dict: &Bound<'_, PyDict>, key: &str, val: Option<f64>) -> PyResult<()> {
    match val {
        Some(v) => dict.set_item(key, v),
        None => dict.set_item(key, dict.py().None()),
    }
}

#[inline]
fn set_opt_bool(dict: &Bound<'_, PyDict>, key: &str, val: Option<bool>) -> PyResult<()> {
    match val {
        Some(v) => dict.set_item(key, v),
        None => dict.set_item(key, dict.py().None()),
    }
}

// ──────────────────────────────────────────────────────────────────
// Validator issue serialization
//
// Each `RecordValidationIssue` variant becomes a Python dict shaped
// as `{ "kind": <PascalCase variant>, ... }`. Variants that involve
// a segment include `segment` ("V" / "D" / "J"). Variants that
// compare engine-reported vs validator-recomputed values surface
// `reported` / `expected` at the top level. Variant-specific extras
// land under `details`.
// ──────────────────────────────────────────────────────────────────

fn segment_str(s: Segment) -> &'static str {
    match s {
        Segment::V => "V",
        Segment::Np1 => "NP1",
        Segment::D => "D",
        Segment::Np2 => "NP2",
        Segment::J => "J",
    }
}

fn prime_end_str(e: PrimeEnd) -> &'static str {
    match e {
        PrimeEnd::Five => "5'",
        PrimeEnd::Three => "3'",
    }
}

fn productive_reason_str(r: &ProductiveDecidedBy) -> &'static str {
    match r {
        ProductiveDecidedBy::OutOfFrame => "out_of_frame",
        ProductiveDecidedBy::JunctionStopCodon => "junction_stop_codon",
        ProductiveDecidedBy::VAnchorAaChanged => "v_anchor_aa_changed",
        ProductiveDecidedBy::JAnchorAaChanged => "j_anchor_aa_changed",
        ProductiveDecidedBy::InFrameAndAnchorsPreserved => "in_frame_and_anchors_preserved",
    }
}

fn allele_order_reason_str(r: &AlleleOrderReason) -> &'static str {
    match r {
        AlleleOrderReason::TruthFirstIfInTieSet => "truth_first_if_in_tie_set",
        AlleleOrderReason::AscendingAlleleIdOtherwise => "ascending_allele_id_otherwise",
    }
}

fn issue_to_pydict<'py>(
    py: Python<'py>,
    issue: RecordValidationIssue,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new_bound(py);
    let details = PyDict::new_bound(py);
    match issue {
        // ── C1: Structural ───────────────────────────────────────
        RecordValidationIssue::SequenceLengthMismatch {
            reported,
            actual_bytes,
        } => {
            d.set_item("kind", "SequenceLengthMismatch")?;
            d.set_item("reported", reported)?;
            d.set_item("expected", actual_bytes as i64)?;
        }
        RecordValidationIssue::SequenceContentMismatch {
            reported_prefix,
            actual_prefix,
        } => {
            d.set_item("kind", "SequenceContentMismatch")?;
            d.set_item("reported", reported_prefix)?;
            d.set_item("expected", actual_prefix)?;
        }
        RecordValidationIssue::SegmentCoordinatesOutOfOrder { segment, start, end } => {
            d.set_item("kind", "SegmentCoordinatesOutOfOrder")?;
            d.set_item("segment", segment_str(segment))?;
            details.set_item("start", start)?;
            details.set_item("end", end)?;
        }
        RecordValidationIssue::GermlineCoordinatesOutOfOrder { segment, start, end } => {
            d.set_item("kind", "GermlineCoordinatesOutOfOrder")?;
            d.set_item("segment", segment_str(segment))?;
            details.set_item("start", start)?;
            details.set_item("end", end)?;
        }
        RecordValidationIssue::CigarReadsInvalid { segment, reason } => {
            d.set_item("kind", "CigarReadsInvalid")?;
            d.set_item("segment", segment_str(segment))?;
            details.set_item("reason", reason)?;
        }
        RecordValidationIssue::CigarSpanMismatch {
            segment,
            cigar_query_span,
            sequence_span,
        } => {
            d.set_item("kind", "CigarSpanMismatch")?;
            d.set_item("segment", segment_str(segment))?;
            d.set_item("reported", cigar_query_span as i64)?;
            d.set_item("expected", sequence_span as i64)?;
        }

        // ── C2: Counter provenance ──────────────────────────────
        RecordValidationIssue::NMutationsMismatch {
            reported,
            sim_count,
        } => {
            d.set_item("kind", "NMutationsMismatch")?;
            d.set_item("reported", reported)?;
            d.set_item("expected", sim_count)?;
            details.set_item("source", "sim.mutation_count")?;
        }
        RecordValidationIssue::NPcrErrorsMismatch {
            reported,
            trace_count,
        } => {
            d.set_item("kind", "NPcrErrorsMismatch")?;
            d.set_item("reported", reported)?;
            d.set_item("expected", trace_count)?;
            details.set_item("source", "trace:corrupt.pcr.count")?;
        }
        RecordValidationIssue::NQualityErrorsMismatch {
            reported,
            trace_count,
        } => {
            d.set_item("kind", "NQualityErrorsMismatch")?;
            d.set_item("reported", reported)?;
            d.set_item("expected", trace_count)?;
            details.set_item("source", "trace:corrupt.quality.count")?;
        }
        RecordValidationIssue::NIndelsMismatch {
            reported,
            event_count,
        } => {
            d.set_item("kind", "NIndelsMismatch")?;
            d.set_item("reported", reported)?;
            d.set_item("expected", event_count)?;
            details.set_item("source", "events:corrupt.indel")?;
        }
        RecordValidationIssue::NSegmentIndelsMismatch {
            segment,
            reported,
            event_count,
        } => {
            d.set_item("kind", "NSegmentIndelsMismatch")?;
            d.set_item("segment", segment_str(segment))?;
            d.set_item("reported", reported)?;
            d.set_item("expected", event_count)?;
            details.set_item("source", "events:corrupt.indel")?;
        }
        RecordValidationIssue::EndLossLengthMismatch {
            side,
            reported,
            trace_count,
        } => {
            d.set_item("kind", "EndLossLengthMismatch")?;
            d.set_item("reported", reported)?;
            d.set_item("expected", trace_count)?;
            details.set_item("side", prime_end_str(side))?;
        }

        // ── C3: Junction truth ──────────────────────────────────
        RecordValidationIssue::JunctionLengthMismatch {
            reported,
            recomputed,
        } => {
            d.set_item("kind", "JunctionLengthMismatch")?;
            set_opt_i64(&d, "reported", reported)?;
            d.set_item("expected", recomputed as i64)?;
        }
        RecordValidationIssue::JunctionContentMismatch {
            reported,
            recomputed,
        } => {
            d.set_item("kind", "JunctionContentMismatch")?;
            d.set_item("reported", reported)?;
            d.set_item("expected", recomputed)?;
        }
        RecordValidationIssue::JunctionAaMismatch {
            reported,
            recomputed,
        } => {
            d.set_item("kind", "JunctionAaMismatch")?;
            d.set_item("reported", reported)?;
            d.set_item("expected", recomputed)?;
        }
        RecordValidationIssue::VjInFrameMismatch {
            reported,
            recomputed,
        } => {
            d.set_item("kind", "VjInFrameMismatch")?;
            set_opt_bool(&d, "reported", reported)?;
            d.set_item("expected", recomputed)?;
        }
        RecordValidationIssue::StopCodonMismatch {
            reported,
            recomputed,
        } => {
            d.set_item("kind", "StopCodonMismatch")?;
            set_opt_bool(&d, "reported", reported)?;
            d.set_item("expected", recomputed)?;
        }
        RecordValidationIssue::ProductiveMismatch {
            reported,
            recomputed,
            reason,
        } => {
            d.set_item("kind", "ProductiveMismatch")?;
            set_opt_bool(&d, "reported", reported)?;
            d.set_item("expected", recomputed)?;
            details.set_item("reason", productive_reason_str(&reason))?;
        }

        // ── C4: Allele oracle ───────────────────────────────────
        RecordValidationIssue::AlleleCallTieSetMismatch {
            segment,
            reported,
            recomputed,
        } => {
            d.set_item("kind", "AlleleCallTieSetMismatch")?;
            d.set_item("segment", segment_str(segment))?;
            d.set_item("reported", reported)?;
            d.set_item("expected", recomputed)?;
        }
        RecordValidationIssue::AlleleCallOrderMismatch {
            segment,
            reported_first,
            expected_first,
            reason,
        } => {
            d.set_item("kind", "AlleleCallOrderMismatch")?;
            d.set_item("segment", segment_str(segment))?;
            d.set_item("reported", reported_first)?;
            d.set_item("expected", expected_first)?;
            details.set_item("reason", allele_order_reason_str(&reason))?;
        }

        // ── C5: Structural invariants ───────────────────────────
        RecordValidationIssue::MultipleRegionsForSegment { segment, count } => {
            d.set_item("kind", "MultipleRegionsForSegment")?;
            d.set_item("segment", segment_str(segment))?;
            details.set_item("count", count as i64)?;
        }
        RecordValidationIssue::MultipleHypothesesInLiveCall { segment, count } => {
            d.set_item("kind", "MultipleHypothesesInLiveCall")?;
            d.set_item("segment", segment_str(segment))?;
            details.set_item("count", count as i64)?;
        }
    }
    d.set_item("details", details)?;
    Ok(d)
}

// ──────────────────────────────────────────────────────────────────
// Live-call cache parity serialization
// ──────────────────────────────────────────────────────────────────

fn hypothesis_bounds_to_pydict<'py>(
    py: Python<'py>,
    bounds: crate::live_call::HypothesisBounds,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new_bound(py);
    d.set_item("seq_start", bounds.seq_start)?;
    d.set_item("seq_end", bounds.seq_end)?;
    d.set_item("ref_start", bounds.ref_start)?;
    d.set_item("ref_end", bounds.ref_end)?;
    Ok(d)
}

fn parity_to_pydict<'py>(
    py: Python<'py>,
    parity: crate::live_call::SegmentParity,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new_bound(py);
    d.set_item("segment", segment_str(parity.segment))?;
    d.set_item("tie_set_matches", parity.tie_set_matches)?;
    d.set_item("cached_present", parity.cached_present)?;
    d.set_item("fresh_present", parity.fresh_present)?;
    d.set_item(
        "cached_tie_set",
        parity
            .cached_tie_set
            .iter()
            .map(|id| id.as_usize() as i64)
            .collect::<Vec<i64>>(),
    )?;
    d.set_item(
        "fresh_tie_set",
        parity
            .fresh_tie_set
            .iter()
            .map(|id| id.as_usize() as i64)
            .collect::<Vec<i64>>(),
    )?;
    match parity.hypothesis_bounds_match {
        Some(v) => d.set_item("hypothesis_bounds_match", v)?,
        None => d.set_item("hypothesis_bounds_match", py.None())?,
    }
    match parity.cached_hypothesis {
        Some(b) => d.set_item("cached_hypothesis", hypothesis_bounds_to_pydict(py, b)?)?,
        None => d.set_item("cached_hypothesis", py.None())?,
    }
    match parity.fresh_hypothesis {
        Some(b) => d.set_item("fresh_hypothesis", hypothesis_bounds_to_pydict(py, b)?)?,
        None => d.set_item("fresh_hypothesis", py.None())?,
    }
    Ok(d)
}
