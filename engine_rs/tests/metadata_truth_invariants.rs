//! Metadata-truth invariants: end-to-end tests that the simulator's
//! per-record metadata is an honest ground truth.
//!
//! GenAIRR's flagship guarantee is that every simulated sequence
//! carries metadata (allele calls, region bounds, junction
//! coordinates, productive flags, mutation counts, etc.) that exactly
//! reflects what biology / a perfect aligner would extract from the
//! sequence given only its bases. This suite exists to catch any
//! drift between the recorded metadata and what the post-pass IR
//! actually supports.
//!
//! See `tests/metadata_truth/README.md` for the test catalog,
//! category boundaries, and pointers to the in-tree tests that
//! cover adjacent ground.

mod metadata_truth;
