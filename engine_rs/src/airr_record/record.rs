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

    // Per-end P-nucleotide length counters (Slice —
    // P-nucleotide v1). Counts of templated palindromic
    // (P-)nucleotide bytes emitted by `PAdditionPass` at each
    // V(D)J coding-end junction side. Aggregated by walking
    // `outcome.events()` filtered to
    // `SimulationEvent::PRegionAdded { end, region }` and
    // summing `region.len()` per `end`. Bases themselves are
    // deterministic from `(allele, trim, orientation,
    // length)` — see [`docs/p_nucleotide_design.md`].
    //
    // VJ chains never carry D-end values (no D segment); v1
    // ships every field unconditionally so the record shape
    // stays uniform across chain types. Records from a VJ
    // simulation have `p_d_5_length == 0` and `p_d_3_length
    // == 0` by construction.
    //
    // The validator surfaces `PLengthMismatch { end }` per
    // end when a downstream record's reported value
    // disagrees with the event-ledger recompute.
    pub p_v_3_length: i64,
    pub p_d_5_length: i64,
    pub p_d_3_length: i64,
    pub p_j_5_length: i64,

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
    /// Per-segment biological-SHM mutation counters. Aggregated by
    /// walking `outcome.events()`, filtering to the SHM passes
    /// (`mutate.uniform`, `mutate.s5f`), and bucketing every
    /// `SimulationEvent::BaseChanged` by its carried `segment`.
    /// NP1 and NP2 events roll into [`Self::n_np_mutations`] — same
    /// rollup the per-segment SHM rates DSL (`segment_rates={"NP":
    /// …}`) uses. Corruption-stage `BaseChanged` events (PCR,
    /// quality, contaminant) are excluded by the pass-name filter
    /// so the four fields cleanly partition realised SHM. The sum
    /// across the four fields equals [`Self::n_mutations`] by
    /// construction; the validator surfaces
    /// [`RecordValidationIssue::MutationCountSumMismatch`] when a
    /// downstream record's reported values violate that invariant.
    /// See `docs/mutation_provenance_audit.md`.
    pub n_v_mutations: i64,
    pub n_d_mutations: i64,
    pub n_j_mutations: i64,
    pub n_np_mutations: i64,
    /// Per-V-subregion biological-SHM mutation counters. Together
    /// with [`Self::n_v_unannotated_mutations`] they form a
    /// **partition of [`Self::n_v_mutations`]**:
    ///
    /// ```text
    /// n_fwr1 + n_cdr1 + n_fwr2 + n_cdr2 + n_fwr3
    ///     + n_v_unannotated == n_v_mutations
    /// ```
    ///
    /// Aggregation walks `outcome.events()` with the same
    /// pass-name filter as the per-segment counters (`mutate.uniform`
    /// + `mutate.s5f` only — PCR / quality / contaminant / receptor
    /// revision base changes are excluded). For each V-segment
    /// `BaseChanged` event, the carried `germline_pos: Option<u16>`
    /// is matched against the assigned V allele's
    /// [`crate::refdata::Allele::subregions`] table; events landing
    /// inside one of the five canonical IMGT intervals route to the
    /// matching bucket, everything else routes to
    /// [`Self::n_v_unannotated_mutations`].
    ///
    /// The validator surfaces six per-field
    /// `N{Fwr1,Cdr1,Fwr2,Cdr2,Fwr3,VUnannotated}MutationsMismatch`
    /// issue kinds plus the cross-field
    /// `VSubregionMutationCountSumMismatch` invariant. See
    /// `docs/v_subregion_mutation_counters_audit.md`.
    pub n_fwr1_mutations: i64,
    pub n_cdr1_mutations: i64,
    pub n_fwr2_mutations: i64,
    pub n_cdr2_mutations: i64,
    pub n_fwr3_mutations: i64,
    /// V-segment SHM events that can't be attributed to one of
    /// the five canonical IMGT V subregions. Three cases:
    ///
    /// - **Unannotated V allele** — assigned V has empty
    ///   `subregions` (legacy cartridge / mixed user cartridge
    ///   with hand-authored alleles).
    /// - **V-side CDR3 contribution** — `germline_pos` is valid
    ///   but lies between `FWR3.end` and `len(allele.seq)`; the
    ///   five-label annotation set deliberately stops at FWR3.
    /// - **Indel-inserted V base** — `germline_pos = None`
    ///   (sentinel for nucleotides without a germline origin).
    ///   Reachable only when an SHM pass runs after an indel
    ///   pass (non-canonical pass order).
    ///
    /// On bundled human OGRDB cartridges (100% V-subregion
    /// coverage) this counter is `0` on every record under the
    /// canonical pass order. Non-zero values signal one of the
    /// three cases above; downstream consumers can audit by
    /// inspecting the assigned V allele and the cartridge
    /// manifest's `v_subregion_support.annotated_v_count`.
    pub n_v_unannotated_mutations: i64,
    /// Number of bases removed from the 5' end of the assembled
    /// sequence by the observation-stage `EndLossPass` (the engine
    /// pass that backs `Experiment.end_loss_5prime`, with
    /// `Experiment.primer_trim_5prime` surviving as a backwards-
    /// compatible alias). Distinct from `v_trim_5`, which is the
    /// recombination-stage exonuclease trim — see
    /// [`docs/primer_trim_end_loss_audit.md`](../../../docs/primer_trim_end_loss_audit.md)
    /// §6.1. Defaults to `0` when no end-loss pass ran.
    pub end_loss_5_length: i64,
    /// Number of bases removed from the 3' end. Mirror of
    /// `end_loss_5_length`.
    pub end_loss_3_length: i64,
    pub is_contaminant: bool,
    /// GenAIRR-specific provenance flag: `true` iff the simulation's
    /// D allele was committed in
    /// [`crate::assignment::SegmentOrientation::ReverseComplement`]
    /// orientation, i.e. an
    /// [`crate::passes::InvertDPass`] decision fired with `Bool(true)`
    /// and AssembleSegment(D) emitted the reverse-complemented D
    /// slice into the pool (Slice B). Defaults to `false` for:
    ///   - VJ chains (no D segment),
    ///   - VDJ chains with no `Experiment.invert_d(...)` step,
    ///   - VDJ chains where InvertDPass fired with `Bool(false)`.
    ///
    /// Sourced from the final IR (`Simulation.assignments`), not the
    /// trace — the post-pipeline orientation is the canonical ground
    /// truth. Trace replay reproduces the same value through the same
    /// commit path. See
    /// [`docs/d_inversion_design.md`](../../../docs/d_inversion_design.md)
    /// §6.3 for the rationale of `d_inverted: bool` over
    /// `d_orientation: String`.
    pub d_inverted: bool,
    /// GenAIRR-specific provenance flag: `true` iff an
    /// [`crate::passes::ReceptorRevisionPass`] fired with
    /// `applied=Bool(true)` during this simulation, replacing the V
    /// slot with a different germline V allele after initial
    /// recombination. Defaults to `false` for:
    ///   - VJ chains (the pass rejects them at the DSL boundary),
    ///   - VDJ chains with no `Experiment.receptor_revision(...)`
    ///     step,
    ///   - VDJ chains where the pass fired with `Bool(false)`.
    ///
    /// Sourced from the trace's `receptor_revision.applied` record.
    /// Unlike `d_inverted` (which reads sim-level orientation
    /// metadata that the final IR carries unambiguously), receptor
    /// revision is a per-simulation Bool the pass emits once; the
    /// trace is the canonical source.
    ///
    /// `v_call` remains the *post-revision* V call — the allele
    /// whose bytes are actually in the assembled pool. See
    /// [`Self::original_v_call`] for the pre-revision V identity.
    pub receptor_revision_applied: bool,
    /// GenAIRR-specific provenance string: the V allele name the
    /// recombine pass originally committed to the V slot, before
    /// `ReceptorRevisionPass` (if any) overwrote it. Empty when:
    ///   - no receptor-revision step ran, or
    ///   - the step ran but `receptor_revision_applied=false`.
    ///
    /// Sourced from the trace's first `sample_allele.v` record;
    /// `v_call` continues to report the live-evidence call (which,
    /// after a successful revision, points at the replacement
    /// allele). Holding both lets downstream consumers
    /// distinguish "not revised" from "revised but unrelated to
    /// the original allele" without a separate boolean diff.
    ///
    /// Design note: empty (not equal to `v_call`) on no-revision is
    /// intentional. Setting `original_v_call == v_call` would
    /// require the consumer to cross-check
    /// `receptor_revision_applied` to disambiguate; the empty
    /// sentinel reads cleaner and matches optional-provenance
    /// style. See `docs/receptor_revision_design.md` §7.
    pub original_v_call: String,

    // ── Paired-end / read layout (Slice A: schema only) ──
    //
    // Eight fields landed in Slice A of the paired-end roadmap
    // ([`docs/paired_end_design.md`](../../../docs/paired_end_design.md)
    // §10). Defaults match the additive-field precedent
    // (`end_loss_5/3_length`, `d_inverted`, `receptor_revision_applied`):
    // every existing record carries these at their `Default`
    // values until the projection layer (Slice B) populates them.
    //
    // Slice A enforces ONE invariant via the validator: if
    // `read_layout` is empty, every other paired-end field must
    // also be at its default. Geometry / window-byte checks land
    // in Slice B/C once the projection helper exists.

    /// String tag for the read-layout convention. v1 sanctioned
    /// values: `""` (legacy single-molecule output — the default),
    /// `"paired_end"` (Slice C+ — R1/R2 windows populated),
    /// `"single_end"` (reserved for an explicit single-read
    /// layout, §2 of the design doc).
    pub read_layout: String,

    /// Forward-strand R1 read substring carved out of `sequence`
    /// at `[r1_start, r1_end)`. Empty under the no-layout default.
    pub r1_sequence: String,

    /// Reverse-complemented R2 read substring carved out of
    /// `sequence` at `[r2_start, r2_end)`. Empty under the
    /// no-layout default.
    pub r2_sequence: String,

    /// 5' position of the R1 window in `sequence`. `None` under
    /// the no-layout default — a `0` sentinel would be ambiguous
    /// because `0` is a valid biological coordinate once the
    /// projection layer lands.
    pub r1_start: Option<i64>,
    /// Exclusive 3' position of the R1 window in `sequence`.
    /// `None` under the no-layout default.
    pub r1_end: Option<i64>,
    /// 5' position of the R2 window in `sequence`. `None` under
    /// the no-layout default.
    pub r2_start: Option<i64>,
    /// Exclusive 3' position of the R2 window in `sequence`.
    /// `None` under the no-layout default.
    pub r2_end: Option<i64>,

    /// Fragment insert size (the distance from R1's 5' to R2's
    /// 3' on the projected molecule). `0` under the no-layout
    /// default — distinct from `None` because the field is
    /// always populated when a layout has been set; the `0`
    /// sentinel reads as "no fragment geometry attached" without
    /// requiring an `Option`.
    pub insert_size: i64,
}
