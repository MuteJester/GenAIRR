/**
 * nuc.h — The annotated nucleotide node.
 *
 * Each Nuc is one nucleotide in the simulated sequence. It carries
 * its current base, germline base, segment identity, germline
 * position, and event flags. Nodes are doubly-linked.
 *
 * Codon rail: when built, each node carries a frame_phase (0, 1, 2)
 * indicating its position within the reading frame. Phase-0 nodes
 * are "codon heads" — they hold the translated amino_acid for
 * positions [phase-0, phase-1, phase-2] and a codon_next skip
 * pointer to the next phase-0 node, forming a fast codon-level
 * traversal rail over the nucleotide list.
 *
 * Nuc instances are always allocated from an ASeq's arena pool —
 * never individually malloc'd.
 */

#ifndef GENAIRR_NUC_H
#define GENAIRR_NUC_H

#include "types.h"

typedef struct Nuc {
    /* ── Primary nucleotide data ─────────────────────────────── */
    char            current;       /* current base (after mutations/errors) */
    char            germline;      /* original germline base ('\0' for NP)  */
    Segment         segment;       /* which segment this belongs to         */
    uint16_t        germline_pos;  /* position within germline allele       */
    uint16_t        flags;         /* NucFlags bitmask                      */

    /* ── Doubly-linked list ──────────────────────────────────── */
    struct Nuc     *next;
    struct Nuc     *prev;

    /* ── Codon rail (valid only when ASeq.codon_rail_valid) ── */
    uint8_t         frame_phase;   /* 0, 1, or 2 within reading frame       */
    char            amino_acid;    /* translated AA (phase-0 only, else 0)   */
    struct Nuc     *codon_next;    /* skip pointer: phase-0 → next phase-0   */

    /* ── Reactive productivity (valid when codon rail is valid) ─ */
    bool            productive;    /* false if this node is part of a stop   */
                                   /* codon or a broken conserved anchor     */
} Nuc;

/* ── Convenience flag checks ──────────────────────────────────── */

static inline bool nuc_is_mutated(const Nuc *n) {
    return (n->flags & NUC_FLAG_MUTATED) != 0;
}

static inline bool nuc_is_germline_segment(const Nuc *n) {
    return n->segment == SEG_V || n->segment == SEG_D ||
           n->segment == SEG_J || n->segment == SEG_C;
}

/* True if the node is part of the receptor's coding region (V/D/J
 * plus the NP segments between them). Used to scope the stop-codon
 * count so adapter / UMI / contaminant stops do not break
 * productivity post-corruption. */
static inline bool nuc_is_coding_segment(const Nuc *n) {
    return n->segment == SEG_V || n->segment == SEG_NP1 ||
           n->segment == SEG_D || n->segment == SEG_NP2 ||
           n->segment == SEG_J;
}

static inline bool nuc_is_anchor(const Nuc *n) {
    return (n->flags & NUC_FLAG_ANCHOR) != 0;
}

/* ── Codon rail convenience ───────────────────────────────────── */

/** True if this node is a codon head (first base of a codon). */
static inline bool nuc_is_codon_head(const Nuc *n) {
    return n->frame_phase == 0;
}

/**
 * Get the codon head for a given node.
 * If the node is phase-0, returns itself.
 * If phase-1, returns prev. If phase-2, returns prev->prev.
 * Returns NULL if codon rail is not built or boundary issues.
 */
static inline struct Nuc *nuc_codon_head_of(struct Nuc *n) {
    if (n->frame_phase == 0) return n;
    if (n->frame_phase == 1 && n->prev) return n->prev;
    if (n->frame_phase == 2 && n->prev && n->prev->prev) return n->prev->prev;
    return NULL;
}

#endif /* GENAIRR_NUC_H */
