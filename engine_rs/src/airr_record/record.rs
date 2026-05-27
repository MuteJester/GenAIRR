/// One AIRR Rearrangement record. All ~50 fields populated as
/// ground truth from the IR, refdata, and trace - no aligner.
///
/// `Option<i64>` for coordinate fields lets us emit `null` /
/// missing values for VJ chains' D coords or pre-recombine sims.
/// Strings are owned `String`s; per-record allocation cost is
/// dwarfed by the actual sequence string lengths anyway.
#[derive(Debug, Clone, Default)]
pub struct AirrRecord {
    // AIRR metadata
    pub sequence_id: String,
    pub sequence: String,
    pub sequence_aa: String,
    pub sequence_alignment: String,
    pub germline_alignment: String,
    pub germline_alignment_d_mask: String,
    pub sequence_length: i64,
    pub rev_comp: bool,
    pub locus: String,

    // V
    pub v_call: String,
    pub v_cigar: String,
    pub v_score: Option<f64>,
    pub v_identity: Option<f64>,
    pub v_support: Option<f64>,
    pub v_sequence_start: Option<i64>,
    pub v_sequence_end: Option<i64>,
    pub v_alignment_start: Option<i64>,
    pub v_alignment_end: Option<i64>,
    pub v_germline_start: Option<i64>,
    pub v_germline_end: Option<i64>,
    pub v_trim_5: i64,
    pub v_trim_3: i64,

    // D
    pub d_call: String,
    pub d_cigar: String,
    pub d_score: Option<f64>,
    pub d_identity: Option<f64>,
    pub d_support: Option<f64>,
    pub d_sequence_start: Option<i64>,
    pub d_sequence_end: Option<i64>,
    pub d_alignment_start: Option<i64>,
    pub d_alignment_end: Option<i64>,
    pub d_germline_start: Option<i64>,
    pub d_germline_end: Option<i64>,
    pub d_trim_5: i64,
    pub d_trim_3: i64,

    // J
    pub j_call: String,
    pub j_cigar: String,
    pub j_score: Option<f64>,
    pub j_identity: Option<f64>,
    pub j_support: Option<f64>,
    pub j_sequence_start: Option<i64>,
    pub j_sequence_end: Option<i64>,
    pub j_alignment_start: Option<i64>,
    pub j_alignment_end: Option<i64>,
    pub j_germline_start: Option<i64>,
    pub j_germline_end: Option<i64>,
    pub j_trim_5: i64,
    pub j_trim_3: i64,

    // C
    pub c_call: String,

    // Junction
    pub junction: String,
    pub junction_aa: String,
    pub junction_start: Option<i64>,
    pub junction_end: Option<i64>,
    pub junction_length: Option<i64>,

    // NP regions
    pub np1: String,
    pub np1_aa: String,
    pub np1_length: i64,
    pub np2: String,
    pub np2_aa: String,
    pub np2_length: i64,

    // Functionality
    pub productive: Option<bool>,
    pub vj_in_frame: Option<bool>,
    pub stop_codon: Option<bool>,

    // SHM + corruption (non-AIRR; GenAIRR additions)
    pub n_mutations: i64,
    pub mutation_rate: f64,
    pub n_pcr_errors: i64,
    pub n_quality_errors: i64,
    pub n_indels: i64,
    /// Per-segment indel counters: the count of `IndelInserted` +
    /// `IndelDeleted` events emitted by the `corrupt.indel` pass
    /// whose `segment` attributes the event to V / D / J.
    /// Insertions attribute via [`insertion_segment`] (region-
    /// containment with end-of-region falling to the following
    /// region); deletions attribute to the deleted nucleotide's
    /// own segment. Events landing in NP1 / NP2 are excluded from
    /// these counters but still contribute to `n_indels`. See
    /// [`docs/indel_provenance_audit.md`](../../../docs/indel_provenance_audit.md)
    /// §6.2.
    pub n_v_indels: i64,
    pub n_d_indels: i64,
    pub n_j_indels: i64,
    /// Number of bases removed from the 5' end of the assembled
    /// sequence by the observation-stage `EndLossPass` (the same
    /// engine pass that backs `Experiment.primer_trim_5prime`).
    /// Distinct from `v_trim_5`, which is the recombination-stage
    /// exonuclease trim — see
    /// [`docs/primer_trim_end_loss_audit.md`](../../../docs/primer_trim_end_loss_audit.md)
    /// §6.1. Defaults to `0` when no end-loss pass ran.
    pub end_loss_5_length: i64,
    /// Number of bases removed from the 3' end. Mirror of
    /// `end_loss_5_length`.
    pub end_loss_3_length: i64,
    pub is_contaminant: bool,
}
