pub(in crate::airr_record) struct AlignmentWalk {
    pub(in crate::airr_record) sa: String,
    pub(in crate::airr_record) galn: String,
    pub(in crate::airr_record) dmask: String,
    /// One CIGAR string per V/D/J segment (in that order).
    pub(in crate::airr_record) cigars: [String; 3],
    /// Per-segment alignment-string spans `(start, end)` (0-based
    /// half-open). `None` when the segment contributed no columns.
    pub(in crate::airr_record) align_ranges: [Option<(i64, i64)>; 3],
    /// Per-segment pool-position spans `(start, end)` (0-based
    /// half-open). Tracks the structural region plus any NP columns
    /// the column-walker claimed for the segment — i.e. the slice of
    /// the sequence that this segment "owns" after live-call overlap
    /// resolution. Used for `*_sequence_start/end`.
    pub(in crate::airr_record) seq_ranges: [Option<(i64, i64)>; 3],
    /// Per-segment allele-position spans `(start, end)` (0-based
    /// half-open). Mirrors `seq_ranges` but in reference space —
    /// the union of ref positions consumed by `M` and `D` ops on
    /// the segment. Used for `*_germline_start/end` so the AIRR
    /// record's germline span matches the CIGAR exactly:
    /// `germline_span == M + D` by construction.
    pub(in crate::airr_record) ref_ranges: [Option<(i64, i64)>; 3],
    /// Per-segment identity (matches / total), or `None` when no
    /// columns. Indexed V=0, D=1, J=2.
    pub(in crate::airr_record) identities: [Option<f64>; 3],
}
