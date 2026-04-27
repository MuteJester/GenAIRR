/**
 * aseq.h — Annotated Sequence: the smart linked list.
 *
 * ASeq is the central data structure. It is a doubly-linked list of
 * Nuc nodes, backed by a contiguous arena pool for cache-friendly
 * allocation. Each nucleotide carries its own segment identity,
 * germline origin, and event flags.
 *
 * Codon Rail:
 *   After assembly, aseq_build_codon_rail() assigns frame_phase to
 *   every node and translates each codon. This enables:
 *     - O(1) amino acid update on point mutation (1 codon retranslation)
 *     - O(1) stop codon count maintenance
 *     - O(1) productivity query via cached aggregates
 *     - Fast codon-level traversal via codon_next skip pointers
 *
 * Segment boundaries are implicit (where tags change) and optionally
 * cached via seg_first[]/seg_last[] pointers for O(1) access.
 *
 * All position-dependent metadata (mutation positions, indel sites,
 * error locations) is intrinsic to the nodes — no external dicts
 * that must be kept in sync.
 *
 * Key operations:
 *   - aseq_append_segment()  — bulk-append a germline segment
 *   - aseq_insert_after()    — splice a single node (indel)
 *   - aseq_delete()          — unlink a node (indel / trim)
 *   - aseq_mutate()          — change a node's base with flag + codon update
 *   - aseq_build_codon_rail()— build frame phase and translate all codons
 *   - aseq_reverse_complement() — reverse list + complement
 *   - aseq_to_string()       — serialize to flat string
 *   - aseq_segment_*()       — segment boundary queries
 */

#ifndef GENAIRR_ASEQ_H
#define GENAIRR_ASEQ_H

#include "nuc.h"

typedef struct {
    /* Arena pool: contiguous memory for all nodes in this sequence. */
    Nuc          pool[GENAIRR_MAX_SEQ_LEN];
    int          pool_used;

    /* Doubly-linked list endpoints. */
    Nuc         *head;
    Nuc         *tail;
    int          length;        /* number of active (linked) nodes */

    /* Segment boundary cache: first/last node of each segment.
     * NULL if segment is not present in the sequence. */
    Nuc         *seg_first[SEG_COUNT];
    Nuc         *seg_last[SEG_COUNT];

    /* ── Codon rail state ────────────────────────────────────── */
    bool         codon_rail_valid;  /* true after aseq_build_codon_rail() */
    Nuc         *codon_head;        /* first phase-0 node (codon rail head)*/
    int          n_codons;          /* number of complete codons */
    int          n_stop_codons;     /* current stop codon count */

    /* ── Reactive productivity state ─────────────────────────── */
    Nuc         *v_anchor_node;     /* cached V anchor (conserved Cys)      */
    Nuc         *j_anchor_node;     /* cached J anchor (conserved W/F)      */
    bool         junction_in_frame; /* junction start/end/len all % 3 == 0  */
} ASeq;

/* ── Lifecycle ────────────────────────────────────────────────── */

/** Initialize an ASeq to empty state. Does not allocate. */
void  aseq_init(ASeq *seq);

/** Reset an ASeq for reuse (clears all nodes, resets pool). */
void  aseq_reset(ASeq *seq);

/* ── Node allocation ──────────────────────────────────────────── */

/**
 * Allocate a new Nuc from the arena pool.
 * Returns NULL if pool is exhausted.
 */
Nuc  *aseq_alloc_node(ASeq *seq);

/* ── Building: append germline segments ───────────────────────── */

/**
 * Append a germline segment to the tail of the sequence.
 *
 * @param seq       The annotated sequence.
 * @param bases     Nucleotide string (e.g. the trimmed V allele).
 * @param len       Length of bases.
 * @param segment   Segment tag for all appended nodes.
 * @param germ_offset  Germline position of the first base (e.g. trim_5).
 * @param anchor_pos   If >= 0, the node at this germline_pos gets FLAG_ANCHOR.
 *                     Pass -1 for no anchor.
 * @return Number of nodes appended, or -1 on pool exhaustion.
 */
int   aseq_append_segment(ASeq *seq, const char *bases, int len,
                          Segment segment, int germ_offset, int anchor_pos);

/**
 * Append NP-region nodes (no germline origin).
 *
 * @param seq       The annotated sequence.
 * @param bases     NP nucleotide string.
 * @param len       Length of bases.
 * @param segment   SEG_NP1 or SEG_NP2.
 * @param np_flags  Base flags for each node (e.g. NUC_FLAG_N_NUCLEOTIDE).
 * @return Number of nodes appended, or -1 on pool exhaustion.
 */
int   aseq_append_np(ASeq *seq, const char *bases, int len,
                     Segment segment, uint16_t np_flags);

/* ── Splicing operations ──────────────────────────────────────── */

/**
 * Insert a new node after the given position.
 * Updates segment cache. If codon rail is valid, propagates frame
 * shift and retranslates affected codons.
 *
 * @return The newly inserted node, or NULL on pool exhaustion.
 */
Nuc  *aseq_insert_after(ASeq *seq, Nuc *pos, char base,
                        Segment segment, uint16_t flags);

/**
 * Insert a new node before the given position.
 * Updates segment cache. If codon rail is valid, propagates frame
 * shift and retranslates affected codons.
 *
 * @return The newly inserted node, or NULL on pool exhaustion.
 */
Nuc  *aseq_insert_before(ASeq *seq, Nuc *pos, char base,
                         Segment segment, uint16_t flags);

/**
 * Delete (unlink) a node from the sequence.
 * Updates segment cache. If codon rail is valid, propagates frame
 * shift and retranslates affected codons.
 * Does NOT free memory (arena-managed).
 */
void  aseq_delete(ASeq *seq, Nuc *node);

/* ── Batch splicing operations ────────────────────────────────── */

/**
 * Delete N nodes from the head of the sequence in one batch.
 * Unlinks all N nodes, updates segment caches, then rebuilds the
 * codon rail ONCE (instead of N separate propagations).
 *
 * @param seq  The annotated sequence.
 * @param n    Number of nodes to delete from head.
 * @return     Number of nodes actually deleted (may be < n if seq too short).
 */
int   aseq_delete_head_n(ASeq *seq, int n);

/**
 * Delete N nodes from the tail of the sequence in one batch.
 * Same optimization as aseq_delete_head_n but from the 3' end.
 *
 * @param seq  The annotated sequence.
 * @param n    Number of nodes to delete from tail.
 * @return     Number of nodes actually deleted.
 */
int   aseq_delete_tail_n(ASeq *seq, int n);

/**
 * Prepend N bases before the head of the sequence in one batch.
 * Allocates N nodes, links them all, then rebuilds the codon rail ONCE.
 *
 * @param seq      The annotated sequence.
 * @param bases    Array of base characters to prepend (index 0 = new head).
 * @param n        Number of bases.
 * @param segment  Segment tag for all new nodes.
 * @param flags    Flags for all new nodes.
 * @return         Number of nodes actually inserted (may be < n on pool exhaustion).
 */
int   aseq_prepend_bases(ASeq *seq, const char *bases, int n,
                         Segment segment, uint16_t flags);

/**
 * Append N bases after the tail of the sequence in one batch.
 * Allocates N nodes, links them all, then rebuilds the codon rail ONCE.
 *
 * @param seq      The annotated sequence.
 * @param bases    Array of base characters to append.
 * @param n        Number of bases.
 * @param segment  Segment tag for all new nodes.
 * @param flags    Flags for all new nodes.
 * @return         Number of nodes actually inserted.
 */
int   aseq_append_bases(ASeq *seq, const char *bases, int n,
                        Segment segment, uint16_t flags);

/* ── Mutation ─────────────────────────────────────────────────── */

/**
 * Mutate a node: change its current base and set the given flag.
 * The germline base is preserved.
 *
 * If the codon rail is valid, retranslates the affected codon and
 * updates the stop codon count. This is O(1) — only the single
 * codon containing this node is retranslated.
 *
 * @param seq       The annotated sequence (needed for codon rail update).
 * @param node      The nucleotide to mutate.
 * @param new_base  The new base character.
 * @param flag      Flag(s) to set (e.g. NUC_FLAG_MUTATED).
 */
void  aseq_mutate(ASeq *seq, Nuc *node, char new_base, uint16_t flag);

/**
 * Revert a node to its germline base.
 * Clears the MUTATED flag. Updates codon rail if valid.
 * Used by selection pressure for reverting unfavored mutations.
 *
 * @param seq   The annotated sequence.
 * @param node  The nucleotide to revert. Must have a valid germline base.
 */
void  aseq_revert(ASeq *seq, Nuc *node);

/* ── Codon rail ───────────────────────────────────────────────── */

/**
 * Build the codon rail over the current sequence.
 *
 * Assigns frame_phase (0, 1, 2) to every node starting from head.
 * Translates each complete codon and stores the amino acid on the
 * phase-0 node. Links phase-0 nodes via codon_next skip pointers.
 * Counts stop codons and stores the aggregate.
 *
 * Call this once after assembly is complete. The rail is automatically
 * maintained by aseq_mutate() for point mutations. Structural ops
 * (insert/delete) propagate the frame shift and retranslate.
 *
 * Calling this on a sequence with an existing rail rebuilds it.
 */
void  aseq_build_codon_rail(ASeq *seq);

/**
 * Invalidate the codon rail.
 * Call this when bulk operations make incremental maintenance
 * impractical (e.g. reverse complement, spike contaminants).
 * After invalidation, codon rail functions are no-ops until
 * aseq_build_codon_rail() is called again.
 */
void  aseq_invalidate_codon_rail(ASeq *seq);

/**
 * Retranslate a single codon starting at a phase-0 node.
 * Updates the amino acid and the stop codon count.
 * Internal helper exposed for testing.
 */
void  aseq_retranslate_codon(ASeq *seq, Nuc *codon_head);

/* ── Reverse complement ───────────────────────────────────────── */

/**
 * Reverse-complement the entire sequence in-place.
 * Reverses the linked list and complements each current base.
 * Flags and germline annotations are preserved.
 * Invalidates the codon rail (must rebuild if needed).
 */
void  aseq_reverse_complement(ASeq *seq);

/* ── Queries ──────────────────────────────────────────────────── */

/** Get the number of active nodes. O(1). */
static inline int aseq_length(const ASeq *seq) {
    return seq->length;
}

/** Check if a segment is present. O(1). */
static inline bool aseq_has_segment(const ASeq *seq, Segment seg) {
    return seq->seg_first[seg] != NULL;
}

/**
 * O(1) check: does the sequence have any stop codons?
 * Only valid when codon_rail_valid is true.
 */
static inline bool aseq_has_stop_codons(const ASeq *seq) {
    return seq->codon_rail_valid && seq->n_stop_codons > 0;
}

/**
 * Reactive productivity query. O(L) with early exit on first non-productive.
 *
 * Each Nuc carries a `productive` flag that is auto-maintained by the
 * codon rail: stop codons and broken conserved anchors mark their
 * codon's 3 nodes as non-productive. The global `junction_in_frame`
 * flag covers the frame constraint (updated on indels).
 *
 * Returns true only if:
 *   - codon rail is valid
 *   - junction is in frame
 *   - every node in the sequence has productive == true
 */
bool aseq_is_productive(const ASeq *seq);

/**
 * Count the length of a specific segment. O(n) on segment length.
 */
int   aseq_segment_length(const ASeq *seq, Segment seg);

/**
 * Compute the absolute position (0-based from head) of a node.
 * O(n) — use only at serialization time, not in hot loops.
 */
int   aseq_position_of(const ASeq *seq, const Nuc *node);

/**
 * Find the first anchor node (FLAG_ANCHOR) in the given segment.
 * Returns NULL if no anchor found.
 */
Nuc  *aseq_find_anchor(const ASeq *seq, Segment seg);

/**
 * Tag every node in the junction span [V Cys, J W/F + 3) with
 * NUC_FLAG_JUNCTION. Returns the junction length on success, or 0
 * when V or J anchors are missing or the W/F codon does not fully
 * fit (no flag changes are made in that case).
 *
 * Called at end-of-assembly so the junction can be tracked by node
 * flags through the rest of the pipeline (mutations preserve flags;
 * corruption deletes nodes; indels auto-inherit). AIRR derive then
 * scans for first/last junction-flagged node to report the live span.
 */
int   aseq_mark_junction(ASeq *seq);

/* ── Serialization ────────────────────────────────────────────── */

/**
 * Write the current bases to a flat string buffer.
 * @param buf    Output buffer (must be at least seq->length + 1 bytes).
 * @param maxlen Size of buf.
 * @return Number of characters written (excluding NUL).
 */
int   aseq_to_string(const ASeq *seq, char *buf, int maxlen);

/**
 * Extract the subsequence for a specific segment into buf.
 * @return Length of the extracted segment.
 */
int   aseq_segment_to_string(const ASeq *seq, Segment seg,
                             char *buf, int maxlen);

/* ── Derived statistics ───────────────────────────────────────── */

typedef struct {
    int  total_germline;    /* number of V+D+J+C nodes */
    int  mutated;           /* nodes with FLAG_MUTATED  */
    int  seq_errors;        /* nodes with FLAG_SEQ_ERROR */
    int  pcr_errors;        /* nodes with FLAG_PCR_ERROR */
    int  insertions;        /* nodes with FLAG_INDEL_INS */
    int  n_bases;           /* nodes with FLAG_IS_N */
    double mutation_rate;   /* mutated / total_germline */
} ASeqStats;

/**
 * Compute derived statistics by walking the list once.
 */
ASeqStats aseq_stats(const ASeq *seq);

/* ── Debug ────────────────────────────────────────────────────── */

/**
 * Print a colored, annotated view of the sequence to stderr.
 * Each segment gets a different color. Mutated bases are highlighted.
 */
void  aseq_debug_print(const ASeq *seq);

#endif /* GENAIRR_ASEQ_H */
