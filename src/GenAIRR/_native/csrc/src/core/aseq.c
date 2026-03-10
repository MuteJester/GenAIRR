/**
 * aseq.c — Annotated Sequence implementation.
 *
 * The ASeq is the heart of GenAIRR-C. All sequence operations
 * (building, mutation, indels, corruption) work through this
 * interface. Positions and metadata are intrinsic to the nodes,
 * so operations never need external bookkeeping synchronization.
 *
 * Codon Rail:
 *   After assembly, aseq_build_codon_rail() assigns reading frame
 *   phases and translates all codons. Point mutations auto-update
 *   the affected codon in O(1). Structural changes (insert/delete)
 *   propagate the frame shift downstream.
 */

#include "genairr/aseq.h"
#include "genairr/codon.h"
#include <string.h>
#include <stdio.h>
#include <assert.h>

/* Forward declarations */
static void aseq_update_junction_in_frame(ASeq *seq);

/* ── Lifecycle ────────────────────────────────────────────────── */

void aseq_init(ASeq *seq) {
    memset(seq, 0, sizeof(*seq));
}

void aseq_reset(ASeq *seq) {
    aseq_init(seq);
}

/* ── Node allocation ──────────────────────────────────────────── */

Nuc *aseq_alloc_node(ASeq *seq) {
    if (seq->pool_used >= GENAIRR_MAX_SEQ_LEN) {
        return NULL;
    }
    Nuc *n = &seq->pool[seq->pool_used++];
    memset(n, 0, sizeof(*n));
    return n;
}

/* ── Internal: link a node at the tail ────────────────────────── */

static void link_at_tail(ASeq *seq, Nuc *node) {
    node->prev = seq->tail;
    node->next = NULL;

    if (seq->tail) {
        seq->tail->next = node;
    } else {
        seq->head = node;
    }
    seq->tail = node;
    seq->length++;
}

/* ── Internal: update segment cache after append ──────────────── */

static void update_seg_cache_append(ASeq *seq, Nuc *node, Segment seg) {
    if (seq->seg_first[seg] == NULL) {
        seq->seg_first[seg] = node;
    }
    seq->seg_last[seg] = node;
}

/* ── Internal: update segment cache after delete ──────────────── */

static void update_seg_cache_delete(ASeq *seq, Nuc *node) {
    Segment seg = node->segment;

    if (seq->seg_first[seg] == node) {
        /* Find next node with same segment, or NULL */
        Nuc *next = node->next;
        seq->seg_first[seg] = (next && next->segment == seg) ? next : NULL;
    }

    if (seq->seg_last[seg] == node) {
        Nuc *prev = node->prev;
        seq->seg_last[seg] = (prev && prev->segment == seg) ? prev : NULL;
    }

    /* If first became NULL but last didn't (or vice versa), clear both */
    if (seq->seg_first[seg] == NULL || seq->seg_last[seg] == NULL) {
        seq->seg_first[seg] = NULL;
        seq->seg_last[seg] = NULL;
    }
}

/* ── Building: append germline segments ───────────────────────── */

int aseq_append_segment(ASeq *seq, const char *bases, int len,
                        Segment segment, int germ_offset, int anchor_pos) {
    for (int i = 0; i < len; i++) {
        Nuc *n = aseq_alloc_node(seq);
        if (!n) return -1;

        n->current     = bases[i];
        n->germline    = bases[i];
        n->segment     = segment;
        n->germline_pos = (uint16_t)(germ_offset + i);
        n->flags       = NUC_FLAG_NONE;

        if (anchor_pos >= 0 && n->germline_pos == (uint16_t)anchor_pos) {
            n->flags |= NUC_FLAG_ANCHOR;
            /* Cache anchor node for reactive productivity */
            if (segment == SEG_V)  seq->v_anchor_node = n;
            if (segment == SEG_J)  seq->j_anchor_node = n;
        }

        n->productive = true;  /* default until codon rail sets otherwise */

        link_at_tail(seq, n);
        update_seg_cache_append(seq, n, segment);
    }
    return len;
}

int aseq_append_np(ASeq *seq, const char *bases, int len,
                   Segment segment, uint16_t np_flags) {
    for (int i = 0; i < len; i++) {
        Nuc *n = aseq_alloc_node(seq);
        if (!n) return -1;

        n->current     = bases[i];
        n->germline    = '\0';           /* NP has no germline */
        n->segment     = segment;
        n->germline_pos = 0xFFFF;        /* sentinel: not germline */
        n->flags       = np_flags;
        n->productive  = true;           /* default until codon rail */

        link_at_tail(seq, n);
        update_seg_cache_append(seq, n, segment);
    }
    return len;
}

/* ═══════════════════════════════════════════════════════════════
 * Codon Rail
 *
 * The codon rail overlays the nucleotide linked list with reading
 * frame information. Each node gets a frame_phase (0, 1, 2).
 * Phase-0 nodes are "codon heads" — they store the translated
 * amino acid and link to the next codon head via codon_next.
 *
 * This enables:
 *   - O(1) amino acid update on point mutation
 *   - O(1) stop codon count maintenance
 *   - Fast codon-level iteration for productivity checks
 * ═══════════════════════════════════════════════════════════════ */

void aseq_retranslate_codon(ASeq *seq, Nuc *ch) {
    if (!ch || ch->frame_phase != 0) return;

    char old_aa = ch->amino_acid;

    Nuc *n2 = ch->next;
    Nuc *n3 = n2 ? n2->next : NULL;

    if (n2 && n3) {
        ch->amino_acid = translate_codon(ch->current, n2->current, n3->current);
    } else {
        /* Incomplete codon at end of sequence */
        ch->amino_acid = '?';
    }

    /* Update stop codon aggregate */
    if (old_aa == '*' && ch->amino_acid != '*') {
        seq->n_stop_codons--;
    } else if (old_aa != '*' && ch->amino_acid == '*') {
        seq->n_stop_codons++;
    }

    /* ── Update per-nucleotide productive flags ──────────────── */
    /* Determine if this codon causes a productivity violation:
     *   - Stop codon at any position
     *   - Anchor codon (V=Cys, J=W/F) with wrong amino acid  */
    bool violation = false;

    if (ch->amino_acid == '*') {
        violation = true;
    }

    /* Check if this codon contains an anchor node */
    if (ch->flags & NUC_FLAG_ANCHOR) {
        if (ch->segment == SEG_V && ch->amino_acid != 'C')
            violation = true;
        else if (ch->segment == SEG_J &&
                 ch->amino_acid != 'W' && ch->amino_acid != 'F')
            violation = true;
    }

    ch->productive = !violation;
    if (n2) n2->productive = !violation;
    if (n3) n3->productive = !violation;
}

void aseq_build_codon_rail(ASeq *seq) {
    /* Reset codon state */
    seq->codon_head = NULL;
    seq->n_codons = 0;
    seq->n_stop_codons = 0;

    /* Assign frame phases: 0, 1, 2, 0, 1, 2, ... */
    int phase = 0;
    Nuc *last_codon_head = NULL;

    for (Nuc *n = seq->head; n; n = n->next) {
        n->frame_phase = (uint8_t)phase;
        n->amino_acid = '\0';
        n->codon_next = NULL;
        n->productive = true;  /* default; overwritten below for violations */

        if (phase == 0) {
            /* This is a codon head */
            if (!seq->codon_head) {
                seq->codon_head = n;
            }

            /* Link previous codon head to this one */
            if (last_codon_head) {
                last_codon_head->codon_next = n;
            }
            last_codon_head = n;

            /* Translate this codon */
            Nuc *n2 = n->next;
            Nuc *n3 = n2 ? n2->next : NULL;

            if (n2 && n3) {
                n->amino_acid = translate_codon(n->current, n2->current, n3->current);
                seq->n_codons++;

                /* Check for productivity violations */
                bool violation = false;
                if (n->amino_acid == '*') {
                    seq->n_stop_codons++;
                    violation = true;
                }
                /* Anchor codon checks */
                if (n->flags & NUC_FLAG_ANCHOR) {
                    if (n->segment == SEG_V && n->amino_acid != 'C')
                        violation = true;
                    else if (n->segment == SEG_J &&
                             n->amino_acid != 'W' && n->amino_acid != 'F')
                        violation = true;
                }
                if (violation) {
                    n->productive  = false;
                    n2->productive = false;
                    n3->productive = false;
                }
            } else {
                /* Incomplete trailing codon */
                n->amino_acid = '?';
            }
        }

        phase = (phase + 1) % 3;
    }

    seq->codon_rail_valid = true;
}

void aseq_invalidate_codon_rail(ASeq *seq) {
    seq->codon_rail_valid = false;
    seq->codon_head = NULL;
    seq->n_codons = 0;
    seq->n_stop_codons = 0;
}

/* ── Internal: propagate frame shift from a node downstream ───── */

/**
 * After an insert or delete, the reading frame shifts for all nodes
 * downstream of the affected position. This function walks from
 * `start` to the end of the list, reassigning frame phases,
 * retranslating affected codons, and relinking codon_next pointers.
 *
 * This is O(L) in the worst case, which is irreducible since the
 * reading frame genuinely shifts for every downstream nucleotide.
 */
static void codon_rail_propagate_from(ASeq *seq, Nuc *start) {
    if (!seq->codon_rail_valid || !start) return;

    /* Determine what phase `start` should have based on its predecessor */
    int phase;
    if (start->prev) {
        phase = (start->prev->frame_phase + 1) % 3;
    } else {
        phase = 0; /* Head of sequence */
    }

    /* Walk downstream, reassigning phases and retranslating */
    seq->n_codons = 0;
    seq->n_stop_codons = 0;
    seq->codon_head = NULL;
    Nuc *last_codon_head = NULL;

    /* We need to recount from the beginning for accurate aggregates.
     * Walk the part BEFORE start to count existing codons. */
    for (Nuc *n = seq->head; n && n != start; n = n->next) {
        if (n->frame_phase == 0) {
            if (!seq->codon_head) seq->codon_head = n;
            if (last_codon_head) last_codon_head->codon_next = n;
            last_codon_head = n;
            n->codon_next = NULL;

            if (n->amino_acid == '*') seq->n_stop_codons++;
            seq->n_codons++;
        }
    }

    /* Now walk from start to end, fixing phases and retranslating */
    for (Nuc *n = start; n; n = n->next) {
        n->frame_phase = (uint8_t)phase;
        n->amino_acid = '\0';
        n->codon_next = NULL;
        n->productive = true;  /* default; overwritten for violations */

        if (phase == 0) {
            if (!seq->codon_head) seq->codon_head = n;
            if (last_codon_head) last_codon_head->codon_next = n;
            last_codon_head = n;

            Nuc *n2 = n->next;
            Nuc *n3 = n2 ? n2->next : NULL;
            if (n2 && n3) {
                n->amino_acid = translate_codon(n->current, n2->current, n3->current);
                seq->n_codons++;

                bool violation = false;
                if (n->amino_acid == '*') {
                    seq->n_stop_codons++;
                    violation = true;
                }
                if (n->flags & NUC_FLAG_ANCHOR) {
                    if (n->segment == SEG_V && n->amino_acid != 'C')
                        violation = true;
                    else if (n->segment == SEG_J &&
                             n->amino_acid != 'W' && n->amino_acid != 'F')
                        violation = true;
                }
                if (violation) {
                    n->productive  = false;
                    n2->productive = false;
                    n3->productive = false;
                }
            } else {
                n->amino_acid = '?';
            }
        }

        phase = (phase + 1) % 3;
    }
}

/* ── Splicing operations ──────────────────────────────────────── */

Nuc *aseq_insert_after(ASeq *seq, Nuc *pos, char base,
                       Segment segment, uint16_t flags) {
    /* Save codon head of insertion point BEFORE splice (phases still valid) */
    Nuc *codon_start = NULL;
    if (seq->codon_rail_valid) {
        codon_start = nuc_codon_head_of(pos);
    }

    Nuc *n = aseq_alloc_node(seq);
    if (!n) return NULL;

    n->current     = base;
    n->germline    = '\0';
    n->segment     = segment;
    n->germline_pos = 0xFFFF;
    n->flags       = flags;

    /* Splice into doubly-linked list */
    n->prev = pos;
    n->next = pos->next;

    if (pos->next) {
        pos->next->prev = n;
    } else {
        seq->tail = n;
    }
    pos->next = n;
    seq->length++;

    /* Update segment cache */
    if (seq->seg_first[segment] == NULL) {
        seq->seg_first[segment] = n;
        seq->seg_last[segment] = n;
    } else if (pos == seq->seg_last[segment] && pos->segment == segment) {
        seq->seg_last[segment] = n;
    }

    /* Propagate codon rail from the codon head that contains the
     * insertion point, since that codon's bases have shifted. */
    if (seq->codon_rail_valid) {
        codon_rail_propagate_from(seq, codon_start ? codon_start : n);
    }

    aseq_update_junction_in_frame(seq);
    return n;
}

Nuc *aseq_insert_before(ASeq *seq, Nuc *pos, char base,
                        Segment segment, uint16_t flags) {
    /* Save codon head of pos BEFORE splice (prev pointers still valid) */
    Nuc *codon_start = NULL;
    if (seq->codon_rail_valid) {
        codon_start = nuc_codon_head_of(pos);
    }

    Nuc *n = aseq_alloc_node(seq);
    if (!n) return NULL;

    n->current     = base;
    n->germline    = '\0';
    n->segment     = segment;
    n->germline_pos = 0xFFFF;
    n->flags       = flags;

    /* Splice */
    n->next = pos;
    n->prev = pos->prev;

    if (pos->prev) {
        pos->prev->next = n;
    } else {
        seq->head = n;
    }
    pos->prev = n;
    seq->length++;

    /* Update segment cache */
    if (seq->seg_first[segment] == NULL) {
        seq->seg_first[segment] = n;
        seq->seg_last[segment] = n;
    } else if (pos == seq->seg_first[segment] && pos->segment == segment) {
        seq->seg_first[segment] = n;
    }

    /* Propagate codon rail from the codon head that contains the
     * insertion point, or from the new node if it's at the head. */
    if (seq->codon_rail_valid) {
        codon_rail_propagate_from(seq, codon_start ? codon_start : n);
    }

    aseq_update_junction_in_frame(seq);
    return n;
}

void aseq_delete(ASeq *seq, Nuc *node) {
    assert(node != NULL);

    /* Save codon head BEFORE unlinking (prev pointers still valid) */
    Nuc *codon_start = NULL;
    if (seq->codon_rail_valid) {
        codon_start = nuc_codon_head_of(node);
    }

    update_seg_cache_delete(seq, node);

    /* ── Clear cached anchor pointers if this node is an anchor ── */
    if (node == seq->v_anchor_node)
        seq->v_anchor_node = NULL;
    if (node == seq->j_anchor_node)
        seq->j_anchor_node = NULL;

    Nuc *after = node->next;

    if (node->prev) {
        node->prev->next = node->next;
    } else {
        seq->head = node->next;
    }

    if (node->next) {
        node->next->prev = node->prev;
    } else {
        seq->tail = node->prev;
    }

    node->next = NULL;
    node->prev = NULL;
    seq->length--;

    /* Propagate codon rail from the codon head that contained the
     * deleted node (that codon's bases may have changed). */
    if (seq->codon_rail_valid) {
        /* codon_start might be the deleted node itself (if it was phase 0).
         * In that case, start from `after` (or the node before, or rebuild). */
        if (codon_start == node) {
            codon_start = after;  /* deleted node was a codon head */
        }
        if (codon_start) {
            codon_rail_propagate_from(seq, codon_start);
        } else {
            /* No propagation start available — rebuild (rare: deleted tail) */
            aseq_build_codon_rail(seq);
        }
    }

    aseq_update_junction_in_frame(seq);
}

/* ── Batch splicing operations ────────────────────────────────── */

int aseq_delete_head_n(ASeq *seq, int n) {
    int deleted = 0;

    for (int i = 0; i < n && seq->head; i++) {
        Nuc *node = seq->head;

        /* Update segment cache */
        update_seg_cache_delete(seq, node);

        /* Clear cached anchor pointers */
        if (node == seq->v_anchor_node) seq->v_anchor_node = NULL;
        if (node == seq->j_anchor_node) seq->j_anchor_node = NULL;

        /* Unlink from head */
        seq->head = node->next;
        if (seq->head) {
            seq->head->prev = NULL;
        } else {
            seq->tail = NULL;
        }
        node->next = NULL;
        node->prev = NULL;
        seq->length--;
        deleted++;
    }

    /* Single codon rail rebuild for all deletions */
    if (deleted > 0 && seq->codon_rail_valid) {
        aseq_build_codon_rail(seq);
    }
    if (deleted > 0) {
        aseq_update_junction_in_frame(seq);
    }

    return deleted;
}

int aseq_delete_tail_n(ASeq *seq, int n) {
    int deleted = 0;

    for (int i = 0; i < n && seq->tail; i++) {
        Nuc *node = seq->tail;

        /* Update segment cache */
        update_seg_cache_delete(seq, node);

        /* Clear cached anchor pointers */
        if (node == seq->v_anchor_node) seq->v_anchor_node = NULL;
        if (node == seq->j_anchor_node) seq->j_anchor_node = NULL;

        /* Unlink from tail */
        seq->tail = node->prev;
        if (seq->tail) {
            seq->tail->next = NULL;
        } else {
            seq->head = NULL;
        }
        node->next = NULL;
        node->prev = NULL;
        seq->length--;
        deleted++;
    }

    /* Single codon rail rebuild for all deletions */
    if (deleted > 0 && seq->codon_rail_valid) {
        aseq_build_codon_rail(seq);
    }
    if (deleted > 0) {
        aseq_update_junction_in_frame(seq);
    }

    return deleted;
}

int aseq_prepend_bases(ASeq *seq, const char *bases, int n,
                       Segment segment, uint16_t flags) {
    int inserted = 0;

    /* Build nodes in reverse so bases[0] ends up at head */
    for (int i = n - 1; i >= 0; i--) {
        Nuc *node = aseq_alloc_node(seq);
        if (!node) break;

        node->current     = bases[i];
        node->germline    = '\0';
        node->segment     = segment;
        node->germline_pos = 0xFFFF;
        node->flags       = flags;

        /* Link before current head */
        node->next = seq->head;
        node->prev = NULL;
        if (seq->head) {
            seq->head->prev = node;
        } else {
            seq->tail = node;
        }
        seq->head = node;
        seq->length++;

        /* Update segment cache */
        if (seq->seg_first[segment] == NULL) {
            seq->seg_first[segment] = node;
            seq->seg_last[segment] = node;
        } else {
            /* New node is before any existing node of this segment */
            seq->seg_first[segment] = node;
        }

        inserted++;
    }

    /* Single codon rail rebuild */
    if (inserted > 0 && seq->codon_rail_valid) {
        aseq_build_codon_rail(seq);
    }
    if (inserted > 0) {
        aseq_update_junction_in_frame(seq);
    }

    return inserted;
}

int aseq_append_bases(ASeq *seq, const char *bases, int n,
                      Segment segment, uint16_t flags) {
    int inserted = 0;

    for (int i = 0; i < n; i++) {
        Nuc *node = aseq_alloc_node(seq);
        if (!node) break;

        node->current     = bases[i];
        node->germline    = '\0';
        node->segment     = segment;
        node->germline_pos = 0xFFFF;
        node->flags       = flags;

        /* Link after current tail */
        node->prev = seq->tail;
        node->next = NULL;
        if (seq->tail) {
            seq->tail->next = node;
        } else {
            seq->head = node;
        }
        seq->tail = node;
        seq->length++;

        /* Update segment cache */
        if (seq->seg_first[segment] == NULL) {
            seq->seg_first[segment] = node;
        }
        seq->seg_last[segment] = node;

        inserted++;
    }

    /* Single codon rail rebuild */
    if (inserted > 0 && seq->codon_rail_valid) {
        aseq_build_codon_rail(seq);
    }
    if (inserted > 0) {
        aseq_update_junction_in_frame(seq);
    }

    return inserted;
}

/* ── Mutation ─────────────────────────────────────────────────── */

void aseq_mutate(ASeq *seq, Nuc *node, char new_base, uint16_t flag) {
    /* If germline is unset (e.g. indel-inserted node), record the
     * pre-mutation base so annotations produce valid "X>Y" strings
     * instead of writing a null byte that truncates the buffer. */
    if (node->germline == '\0') {
        node->germline = node->current;
    }
    node->current = new_base;
    node->flags |= flag;

    /* If codon rail is active, retranslate the affected codon in O(1) */
    if (seq->codon_rail_valid) {
        Nuc *ch = nuc_codon_head_of(node);
        if (ch) {
            aseq_retranslate_codon(seq, ch);
        }
    }
}

void aseq_revert(ASeq *seq, Nuc *node) {
    if (node->germline == '\0') return;  /* NP node — cannot revert */

    node->current = node->germline;
    node->flags &= ~(uint16_t)NUC_FLAG_MUTATED;

    if (seq->codon_rail_valid) {
        Nuc *ch = nuc_codon_head_of(node);
        if (ch) {
            aseq_retranslate_codon(seq, ch);
        }
    }
}

/* ── Reverse complement ───────────────────────────────────────── */

void aseq_reverse_complement(ASeq *seq) {
    /* Reverse the linked list */
    Nuc *cur = seq->head;
    while (cur) {
        /* Swap prev and next */
        Nuc *tmp = cur->next;
        cur->next = cur->prev;
        cur->prev = tmp;

        /* Complement the current base */
        cur->current = complement(cur->current);

        cur = tmp;  /* advance (was next, now prev) */
    }

    /* Swap head and tail */
    Nuc *tmp = seq->head;
    seq->head = seq->tail;
    seq->tail = tmp;

    /* Swap segment cache first/last pointers */
    for (int s = 0; s < SEG_COUNT; s++) {
        Nuc *t = seq->seg_first[s];
        seq->seg_first[s] = seq->seg_last[s];
        seq->seg_last[s] = t;
    }

    /* Invalidate codon rail — reading frame is reversed */
    aseq_invalidate_codon_rail(seq);
}

/* ── Queries ──────────────────────────────────────────────────── */

int aseq_segment_length(const ASeq *seq, Segment seg) {
    if (!seq->seg_first[seg]) return 0;

    int count = 0;
    for (Nuc *n = seq->seg_first[seg]; n && n->segment == seg; n = n->next) {
        count++;
    }
    return count;
}

int aseq_position_of(const ASeq *seq, const Nuc *node) {
    int pos = 0;
    for (Nuc *n = seq->head; n; n = n->next) {
        if (n == node) return pos;
        pos++;
    }
    return -1;  /* not found */
}

Nuc *aseq_find_anchor(const ASeq *seq, Segment seg) {
    if (!seq->seg_first[seg]) return NULL;

    for (Nuc *n = seq->seg_first[seg]; n && n->segment == seg; n = n->next) {
        if (n->flags & NUC_FLAG_ANCHOR) return n;
    }
    return NULL;
}

/* ── Serialization ────────────────────────────────────────────── */

int aseq_to_string(const ASeq *seq, char *buf, int maxlen) {
    int i = 0;
    for (Nuc *n = seq->head; n && i < maxlen - 1; n = n->next) {
        buf[i++] = n->current;
    }
    buf[i] = '\0';
    return i;
}

int aseq_segment_to_string(const ASeq *seq, Segment seg,
                           char *buf, int maxlen) {
    if (!seq->seg_first[seg]) {
        buf[0] = '\0';
        return 0;
    }

    int i = 0;
    for (Nuc *n = seq->seg_first[seg]; n && n->segment == seg && i < maxlen - 1; n = n->next) {
        buf[i++] = n->current;
    }
    buf[i] = '\0';
    return i;
}

/* ── Reactive productivity ────────────────────────────────────── */

/**
 * Recompute junction_in_frame after structural changes (insert/delete).
 * Matches the FrameAlignmentRule from functionality.c:
 *   junction_start % 3 == 0 && junction_end % 3 == 0 && junction_len % 3 == 0
 */
static void aseq_update_junction_in_frame(ASeq *seq) {
    Nuc *v_anchor = seq->v_anchor_node;
    Nuc *j_anchor = seq->j_anchor_node;

    if (!v_anchor || !j_anchor) {
        seq->junction_in_frame = false;
        return;
    }

    int jstart = aseq_position_of(seq, v_anchor);
    int jend   = aseq_position_of(seq, j_anchor) + 3;
    int jlen   = jend - jstart;

    seq->junction_in_frame = (jstart % 3 == 0) &&
                             (jend   % 3 == 0) &&
                             (jlen   % 3 == 0);
}

bool aseq_is_productive(const ASeq *seq) {
    if (!seq->codon_rail_valid || !seq->junction_in_frame)
        return false;
    /* Both anchors must still exist (may be deleted by 5'/3' corruption) */
    if (!seq->v_anchor_node || !seq->j_anchor_node)
        return false;
    for (Nuc *n = seq->head; n; n = n->next) {
        if (!n->productive) return false;
    }
    return true;
}

/* ── Derived statistics ───────────────────────────────────────── */

ASeqStats aseq_stats(const ASeq *seq) {
    ASeqStats s = {0};

    for (Nuc *n = seq->head; n; n = n->next) {
        if (nuc_is_germline_segment(n)) {
            s.total_germline++;
            if (n->flags & NUC_FLAG_MUTATED)  s.mutated++;
        }
        if (n->flags & NUC_FLAG_SEQ_ERROR)  s.seq_errors++;
        if (n->flags & NUC_FLAG_PCR_ERROR)  s.pcr_errors++;
        if (n->flags & NUC_FLAG_INDEL_INS)  s.insertions++;
        if (n->flags & NUC_FLAG_IS_N)       s.n_bases++;
    }

    s.mutation_rate = s.total_germline > 0
        ? (double)s.mutated / s.total_germline
        : 0.0;

    return s;
}

/* ── Debug ────────────────────────────────────────────────────── */

/* ANSI color codes for segments */
static const char *seg_colors[] = {
    "\033[34m",   /* SEG_V:    blue     */
    "\033[33m",   /* SEG_NP1:  yellow   */
    "\033[31m",   /* SEG_D:    red      */
    "\033[35m",   /* SEG_NP2:  magenta  */
    "\033[32m",   /* SEG_J:    green    */
    "\033[36m",   /* SEG_C:    cyan     */
    "\033[37m",   /* SEG_UMI:  white    */
    "\033[90m",   /* SEG_ADAPTER: gray  */
};

static const char *seg_names[] = {
    "V", "NP1", "D", "NP2", "J", "C", "UMI", "ADP"
};

void aseq_debug_print(const ASeq *seq) {
    fprintf(stderr, "ASeq [%d nodes]:\n  ", seq->length);

    Segment last_seg = SEG_COUNT;
    for (Nuc *n = seq->head; n; n = n->next) {
        if (n->segment != last_seg) {
            if (last_seg != SEG_COUNT) fprintf(stderr, "\033[0m|");
            fprintf(stderr, "%s[%s:", seg_colors[n->segment], seg_names[n->segment]);
            last_seg = n->segment;
        }

        if (n->flags & NUC_FLAG_MUTATED) {
            fprintf(stderr, "\033[1;4m%c\033[22;24m", n->current);
        } else {
            fprintf(stderr, "%c", n->current);
        }
    }
    fprintf(stderr, "\033[0m]\n");

    /* Print segment summary */
    fprintf(stderr, "  Segments: ");
    for (int s = 0; s < SEG_COUNT; s++) {
        if (seq->seg_first[s]) {
            fprintf(stderr, "%s=%d ", seg_names[s], aseq_segment_length(seq, s));
        }
    }

    /* Print codon rail info if valid */
    if (seq->codon_rail_valid) {
        fprintf(stderr, "| Codons=%d Stops=%d", seq->n_codons, seq->n_stop_codons);
    }
    fprintf(stderr, "\n");
}
