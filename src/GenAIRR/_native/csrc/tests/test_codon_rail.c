/**
 * test_codon_rail.c — Tests for the codon rail overlay on ASeq.
 *
 * Tests:
 *   - Frame phase assignment (0, 1, 2, 0, 1, 2, ...)
 *   - Amino acid translation on phase-0 nodes
 *   - codon_next skip pointer chain
 *   - Stop codon counting
 *   - O(1) retranslation on point mutation
 *   - Frame propagation on insert/delete
 *   - Codon rail invalidation on reverse complement
 *   - aseq_revert with codon rail
 */

#include "genairr/aseq.h"
#include "genairr/codon.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

static int tests_passed = 0;
static int tests_total  = 0;

#define RUN_TEST(fn)  do {                                          \
    tests_total++;                                                  \
    printf("  %-55s ", #fn);                                        \
    fn();                                                           \
    printf("PASS\n");                                               \
    tests_passed++;                                                 \
} while (0)

/* ── Helpers ─────────────────────────────────────────────────── */

static void build_simple_seq(ASeq *seq, const char *bases, int len) {
    aseq_init(seq);
    aseq_append_segment(seq, bases, len, SEG_V, 0, -1);
}

/* ═══════════════════════════════════════════════════════════════
 * Frame phase assignment
 * ═══════════════════════════════════════════════════════════════ */

static void test_frame_phase_assignment(void) {
    ASeq seq;
    build_simple_seq(&seq, "ATGATGATG", 9);  /* 3 complete codons */

    aseq_build_codon_rail(&seq);
    assert(seq.codon_rail_valid);

    /* Check phases: 0,1,2,0,1,2,0,1,2 */
    int phase = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        assert(n->frame_phase == phase);
        phase = (phase + 1) % 3;
    }

    aseq_reset(&seq);
}

static void test_frame_phase_incomplete_codon(void) {
    ASeq seq;
    build_simple_seq(&seq, "ATGATGA", 7);  /* 2 complete codons + 1 leftover */

    aseq_build_codon_rail(&seq);
    assert(seq.n_codons == 2);  /* only 2 complete codons */

    /* The last base (phase 0) has amino_acid = '?' (incomplete) */
    Nuc *last_codon_head = NULL;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->frame_phase == 0) last_codon_head = n;
    }
    assert(last_codon_head != NULL);
    assert(last_codon_head->amino_acid == '?');

    aseq_reset(&seq);
}

/* ═══════════════════════════════════════════════════════════════
 * Amino acid translation
 * ═══════════════════════════════════════════════════════════════ */

static void test_amino_acid_translation(void) {
    /* ATG ATG TGT = M M C */
    ASeq seq;
    build_simple_seq(&seq, "ATGATGTGT", 9);

    aseq_build_codon_rail(&seq);

    Nuc *n = seq.codon_head;
    assert(n != NULL);
    assert(n->amino_acid == 'M');  /* ATG */

    n = n->codon_next;
    assert(n != NULL);
    assert(n->amino_acid == 'M');  /* ATG */

    n = n->codon_next;
    assert(n != NULL);
    assert(n->amino_acid == 'C');  /* TGT */

    assert(n->codon_next == NULL);  /* end of chain */

    aseq_reset(&seq);
}

static void test_codon_next_chain(void) {
    ASeq seq;
    build_simple_seq(&seq, "ATGATGTGTGGG", 12);  /* 4 codons */

    aseq_build_codon_rail(&seq);
    assert(seq.n_codons == 4);

    /* Walk codon_next chain and count */
    int count = 0;
    for (Nuc *n = seq.codon_head; n; n = n->codon_next) {
        count++;
        assert(n->frame_phase == 0);
    }
    assert(count == 4);

    aseq_reset(&seq);
}

/* ═══════════════════════════════════════════════════════════════
 * Stop codon counting
 * ═══════════════════════════════════════════════════════════════ */

static void test_stop_codon_count_none(void) {
    ASeq seq;
    build_simple_seq(&seq, "ATGATGTGT", 9);

    aseq_build_codon_rail(&seq);
    assert(seq.n_stop_codons == 0);

    aseq_reset(&seq);
}

static void test_stop_codon_count_one(void) {
    /* TAA = stop codon at pos 0 */
    ASeq seq;
    build_simple_seq(&seq, "TAAATGTGT", 9);

    aseq_build_codon_rail(&seq);
    assert(seq.n_stop_codons == 1);

    aseq_reset(&seq);
}

static void test_stop_codon_count_multiple(void) {
    /* TAA TAA TGA = 3 stop codons */
    ASeq seq;
    build_simple_seq(&seq, "TAATAATGA", 9);

    aseq_build_codon_rail(&seq);
    assert(seq.n_stop_codons == 3);

    aseq_reset(&seq);
}

/* ═══════════════════════════════════════════════════════════════
 * Point mutation retranslation
 * ═══════════════════════════════════════════════════════════════ */

static void test_mutate_retranslates_codon(void) {
    /* ATG ATG TGT → M M C */
    ASeq seq;
    build_simple_seq(&seq, "ATGATGTGT", 9);
    aseq_build_codon_rail(&seq);

    /* Verify initial state */
    assert(seq.codon_head->amino_acid == 'M');
    assert(seq.n_stop_codons == 0);

    /* Mutate pos 0: A→T → TTG ATG TGT → L M C */
    aseq_mutate(&seq, seq.head, 'T', NUC_FLAG_MUTATED);
    assert(seq.codon_head->amino_acid == 'L');  /* TTG = Leu */
    assert(seq.n_stop_codons == 0);

    aseq_reset(&seq);
}

static void test_mutate_creates_stop_codon(void) {
    /* ATG ATG TGT → M M C */
    ASeq seq;
    build_simple_seq(&seq, "ATGATGTGT", 9);
    aseq_build_codon_rail(&seq);
    assert(seq.n_stop_codons == 0);

    /* Mutate pos 2: G→A → ATA ATG TGT. No, that's I M C. */
    /* To create TAA: mutate pos 0→T, pos 1→A, pos 2→A */
    Nuc *n0 = seq.head;
    Nuc *n1 = n0->next;
    Nuc *n2 = n1->next;

    aseq_mutate(&seq, n0, 'T', NUC_FLAG_MUTATED);
    aseq_mutate(&seq, n1, 'A', NUC_FLAG_MUTATED);
    aseq_mutate(&seq, n2, 'A', NUC_FLAG_MUTATED);

    assert(seq.codon_head->amino_acid == '*');  /* TAA = stop */
    assert(seq.n_stop_codons == 1);

    aseq_reset(&seq);
}

static void test_mutate_removes_stop_codon(void) {
    /* TAA ATG TGT → stop M C */
    ASeq seq;
    build_simple_seq(&seq, "TAAATGTGT", 9);
    aseq_build_codon_rail(&seq);
    assert(seq.n_stop_codons == 1);

    /* Mutate pos 0: T→A → AAA ATG TGT → K M C */
    aseq_mutate(&seq, seq.head, 'A', NUC_FLAG_MUTATED);
    assert(seq.codon_head->amino_acid == 'K');  /* AAA = Lys */
    assert(seq.n_stop_codons == 0);

    aseq_reset(&seq);
}

static void test_mutate_phase1_retranslates_correct_codon(void) {
    /* ATG ATG TGT → M M C */
    ASeq seq;
    build_simple_seq(&seq, "ATGATGTGT", 9);
    aseq_build_codon_rail(&seq);

    /* Mutate pos 1 (phase 1): T→C → ACG ATG TGT → T M C */
    Nuc *n1 = seq.head->next;
    assert(n1->frame_phase == 1);

    aseq_mutate(&seq, n1, 'C', NUC_FLAG_MUTATED);
    assert(seq.codon_head->amino_acid == 'T');  /* ACG = Thr */

    aseq_reset(&seq);
}

static void test_mutate_phase2_retranslates_correct_codon(void) {
    /* ATG ATG TGT → M M C */
    ASeq seq;
    build_simple_seq(&seq, "ATGATGTGT", 9);
    aseq_build_codon_rail(&seq);

    /* Mutate pos 2 (phase 2): G→A → ATA ATG TGT → I M C */
    Nuc *n2 = seq.head->next->next;
    assert(n2->frame_phase == 2);

    aseq_mutate(&seq, n2, 'A', NUC_FLAG_MUTATED);
    assert(seq.codon_head->amino_acid == 'I');  /* ATA = Ile */

    aseq_reset(&seq);
}

/* ═══════════════════════════════════════════════════════════════
 * Revert
 * ═══════════════════════════════════════════════════════════════ */

static void test_revert_restores_amino_acid(void) {
    /* ATG → M. Mutate to TAA (stop). Revert → back to ATG → M */
    ASeq seq;
    build_simple_seq(&seq, "ATGATG", 6);
    aseq_build_codon_rail(&seq);

    assert(seq.codon_head->amino_acid == 'M');
    assert(seq.n_stop_codons == 0);

    /* Mutate pos 0: A→T → TTG = L */
    aseq_mutate(&seq, seq.head, 'T', NUC_FLAG_MUTATED);
    assert(seq.codon_head->amino_acid == 'L');

    /* Revert */
    aseq_revert(&seq, seq.head);
    assert(seq.head->current == 'A');
    assert(seq.codon_head->amino_acid == 'M');

    aseq_reset(&seq);
}

/* ═══════════════════════════════════════════════════════════════
 * Insert / delete (frame propagation)
 * ═══════════════════════════════════════════════════════════════ */

static void test_insert_shifts_frame(void) {
    /* ATG ATG → M M (2 codons) */
    ASeq seq;
    build_simple_seq(&seq, "ATGATG", 6);
    aseq_build_codon_rail(&seq);

    assert(seq.n_codons == 2);
    assert(seq.codon_head->amino_acid == 'M');

    /* Insert 'C' after pos 0 → A C T G A T G → 7 bases, 2 complete codons */
    /* New codons: ACT (T), GAT (D) + trailing G */
    Nuc *inserted = aseq_insert_after(&seq, seq.head, 'C', SEG_V, 0);
    assert(inserted != NULL);
    assert(seq.length == 7);
    assert(seq.codon_rail_valid);

    /* First codon should now be ACT = T */
    assert(seq.codon_head->amino_acid == 'T');

    aseq_reset(&seq);
}

static void test_delete_shifts_frame(void) {
    /* ATG ATG TGT → M M C (3 codons) */
    ASeq seq;
    build_simple_seq(&seq, "ATGATGTGT", 9);
    aseq_build_codon_rail(&seq);
    assert(seq.n_codons == 3);

    /* Delete first node (A) → TG ATG TGT → 8 bases, 2 complete codons */
    /* Codons: TGA (stop!), TGT (?? only 2 more = TG) */
    /* Actually: T G A T G T G T → TGA TGT + GT */
    Nuc *to_delete = seq.head;
    aseq_delete(&seq, to_delete);
    assert(seq.length == 8);
    assert(seq.codon_rail_valid);

    /* First codon = TGA = stop */
    assert(seq.codon_head->amino_acid == '*');
    assert(seq.n_stop_codons == 1);

    aseq_reset(&seq);
}

/* ═══════════════════════════════════════════════════════════════
 * Reverse complement invalidation
 * ═══════════════════════════════════════════════════════════════ */

static void test_reverse_complement_invalidates_rail(void) {
    ASeq seq;
    build_simple_seq(&seq, "ATGATG", 6);
    aseq_build_codon_rail(&seq);
    assert(seq.codon_rail_valid);

    aseq_reverse_complement(&seq);
    assert(!seq.codon_rail_valid);

    aseq_reset(&seq);
}

/* ═══════════════════════════════════════════════════════════════
 * Rebuild
 * ═══════════════════════════════════════════════════════════════ */

static void test_rebuild_codon_rail(void) {
    ASeq seq;
    build_simple_seq(&seq, "ATGATG", 6);
    aseq_build_codon_rail(&seq);
    assert(seq.codon_rail_valid);
    assert(seq.n_codons == 2);

    /* Invalidate */
    aseq_invalidate_codon_rail(&seq);
    assert(!seq.codon_rail_valid);

    /* Rebuild */
    aseq_build_codon_rail(&seq);
    assert(seq.codon_rail_valid);
    assert(seq.n_codons == 2);
    assert(seq.codon_head->amino_acid == 'M');

    aseq_reset(&seq);
}

/* ═══════════════════════════════════════════════════════════════
 * Multi-segment sequence
 * ═══════════════════════════════════════════════════════════════ */

static void test_codon_rail_multi_segment(void) {
    /* V(6bp) + NP1(3bp) + J(6bp) = 15bp = 5 codons */
    ASeq seq;
    aseq_init(&seq);
    aseq_append_segment(&seq, "ATGATG", 6, SEG_V, 0, 3);
    aseq_append_np(&seq, "GGG", 3, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "CCCTGG", 6, SEG_J, 0, 3);

    aseq_build_codon_rail(&seq);
    assert(seq.codon_rail_valid);
    assert(seq.n_codons == 5);

    /* Codons: ATG ATG GGG CCC TGG → M M G P W */
    Nuc *c = seq.codon_head;
    assert(c->amino_acid == 'M'); c = c->codon_next;
    assert(c->amino_acid == 'M'); c = c->codon_next;
    assert(c->amino_acid == 'G'); c = c->codon_next;
    assert(c->amino_acid == 'P'); c = c->codon_next;
    assert(c->amino_acid == 'W'); c = c->codon_next;
    assert(c == NULL);

    /* Phases cross segment boundaries correctly */
    int phase = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        assert(n->frame_phase == phase);
        phase = (phase + 1) % 3;
    }

    aseq_reset(&seq);
}

static void test_aseq_has_stop_codons(void) {
    ASeq seq;
    build_simple_seq(&seq, "ATGATG", 6);
    aseq_build_codon_rail(&seq);
    assert(!aseq_has_stop_codons(&seq));

    /* Mutate to create stop */
    Nuc *n0 = seq.head;
    Nuc *n1 = n0->next;
    Nuc *n2 = n1->next;
    aseq_mutate(&seq, n0, 'T', NUC_FLAG_MUTATED);
    aseq_mutate(&seq, n1, 'A', NUC_FLAG_MUTATED);
    aseq_mutate(&seq, n2, 'A', NUC_FLAG_MUTATED);
    assert(aseq_has_stop_codons(&seq));

    aseq_reset(&seq);
}

/* ═══════════════════════════════════════════════════════════════
 * V-anchored reading frame (T0-4 regression)
 *
 * After 5' corruption (delete bases from head, prepend adapter at
 * head, UMI prepend, etc.) the codon rail must continue to use V's
 * reading frame, not seq->head's. The rail derives initial phase
 * from V's first surviving node and its germline_pos.
 * ═══════════════════════════════════════════════════════════════ */

/* Build an assembled-style sequence with V (germline_pos 0..) and
 * J (anchor at known V-frame phase 0) — the minimal shape the rail
 * inspects. */
static void build_v_then_j(ASeq *seq, const char *v_bases, int v_len,
                           int v_anchor, const char *j_bases, int j_len,
                           int j_anchor) {
    aseq_init(seq);
    aseq_append_segment(seq, v_bases, v_len, SEG_V, 0, v_anchor);
    aseq_append_segment(seq, j_bases, j_len, SEG_J, 0, j_anchor);
}

/* Without head changes the V-anchored rail must match the legacy
 * head-anchored rail: when V starts at index 0 with germline_pos 0,
 * initial_phase = 0. */
static void test_v_anchored_no_head_change(void) {
    ASeq seq;
    /* V len 9, anchor at 6 (codon head, %3==0).
     * J len 6, anchor at 3 (codon head). Junction = anchor..end. */
    build_v_then_j(&seq, "ATGTGTTGT", 9, 6, "TGGGGC", 6, 3);

    aseq_build_codon_rail(&seq);

    /* Phase 0 starts at head. */
    assert(seq.head->frame_phase == 0);
    /* V anchor (position 6) is at phase 0 (codon head). */
    Nuc *va = aseq_find_anchor(&seq, SEG_V);
    assert(va && va->frame_phase == 0);
    /* J anchor (position 9+3=12) is at phase 0. */
    Nuc *ja = aseq_find_anchor(&seq, SEG_J);
    assert(ja && ja->frame_phase == 0);

    aseq_reset(&seq);
}

/* Delete N bases (N % 3 != 0) from the head: the rebuilt rail must
 * keep V_anchor at phase 0 (V's frame is intact, only the head
 * changed). Pre-T0-4 the rail would re-anchor to the new head and
 * V_anchor's phase would drift. */
static void test_v_anchored_after_5prime_delete(void) {
    ASeq seq;
    build_v_then_j(&seq, "ATGTGTTGT", 9, 6, "TGGGGC", 6, 3);
    aseq_build_codon_rail(&seq);
    seq.v_anchor_node = aseq_find_anchor(&seq, SEG_V);
    seq.j_anchor_node = aseq_find_anchor(&seq, SEG_J);

    /* Delete 2 bases (≡ 2 mod 3) from the head. V's first surviving
     * node now has germline_pos = 2; head's initial phase should
     * be (2 - 0) % 3 = 2 so that V's first node lands on phase 2,
     * and V[6] (the anchor) lands on phase (6-2) % 3 = 1... wait,
     * we want V_anchor at phase 0. Walking from head (phase 2),
     * positions 0,1,2,...: phases 2,0,1,2,0,1,2,0,1,...
     *   V's first surviving node (was V[2]) at seq pos 0 → phase 2 ✓
     *   V_anchor (was V[6]) at seq pos 4 → phase 0 ✓ */
    aseq_delete_head_n(&seq, 2);

    Nuc *va = aseq_find_anchor(&seq, SEG_V);
    Nuc *ja = aseq_find_anchor(&seq, SEG_J);
    assert(va && va->frame_phase == 0);
    assert(ja && ja->frame_phase == 0);
    assert(seq.junction_in_frame);

    aseq_reset(&seq);
}

/* Prepend N adapter bases (N % 3 != 0): rebuilt rail keeps V_anchor
 * at phase 0. */
static void test_v_anchored_after_5prime_prepend(void) {
    ASeq seq;
    build_v_then_j(&seq, "ATGTGTTGT", 9, 6, "TGGGGC", 6, 3);
    aseq_build_codon_rail(&seq);
    seq.v_anchor_node = aseq_find_anchor(&seq, SEG_V);
    seq.j_anchor_node = aseq_find_anchor(&seq, SEG_J);

    /* Prepend 2 adapter bases. V_anchor's seq position becomes 6+2=8;
     * head's initial phase must be such that V[0] (seq pos 2) lands
     * on phase 0 → head phase = (0 - 2) % 3 = 1. V_anchor at seq
     * pos 8 → phase (1 + 8) % 3 = 0. ✓ */
    aseq_prepend_bases(&seq, "AC", 2, SEG_ADAPTER, 0);

    Nuc *va = aseq_find_anchor(&seq, SEG_V);
    Nuc *ja = aseq_find_anchor(&seq, SEG_J);
    assert(va && va->frame_phase == 0);
    assert(ja && ja->frame_phase == 0);
    assert(seq.junction_in_frame);

    aseq_reset(&seq);
}

/* Stop codons in the 5' adapter must NOT be counted toward
 * n_stop_codons (they are sequencing artifacts, not biology). */
static void test_stop_codon_in_adapter_ignored(void) {
    ASeq seq;
    /* V starts with ATG (start codon, no stops). Anchor at 6. */
    build_v_then_j(&seq, "ATGGCAGCA", 9, 6, "TGGGGC", 6, 3);
    aseq_build_codon_rail(&seq);
    seq.v_anchor_node = aseq_find_anchor(&seq, SEG_V);
    seq.j_anchor_node = aseq_find_anchor(&seq, SEG_J);

    int stops_pre = seq.n_stop_codons;

    /* Prepend 3 adapter bases that translate to a stop codon (TAA).
     * Whatever frame the adapter ends up in, if a TAA-style codon
     * lands on a phase-0 SEG_ADAPTER node it must NOT bump
     * n_stop_codons. */
    aseq_prepend_bases(&seq, "TAA", 3, SEG_ADAPTER, 0);

    /* V's first node is now at seq pos 3, germline_pos 0 → phase 0.
     * Head's initial phase = (0 - 3) % 3 = 0. So adapter[0..2]
     * are at phases 0,1,2 — adapter[0] is a phase-0 codon head. */
    int stops_post = seq.n_stop_codons;
    assert(stops_post == stops_pre);  /* adapter stop ignored */

    aseq_reset(&seq);
}

/* When V is entirely missing (no SEG_V nodes), the rail falls back
 * to phase 0 from head and aseq_is_productive returns false because
 * v_anchor_node is NULL. */
static void test_no_v_falls_back_not_productive(void) {
    ASeq seq;
    aseq_init(&seq);
    /* Pure NP and J — no V at all. */
    aseq_append_segment(&seq, "TGGGGC", 6, SEG_J, 0, 3);
    aseq_build_codon_rail(&seq);

    /* No V → no anchor cached → not productive. */
    assert(seq.v_anchor_node == NULL);
    assert(!aseq_is_productive(&seq));

    aseq_reset(&seq);
}

/* ═══════════════════════════════════════════════════════════════
 * Main
 * ═══════════════════════════════════════════════════════════════ */

int main(void) {
    printf("test_codon_rail: Codon rail overlay tests\n");
    printf("═══════════════════════════════════════════════════════════\n\n");

    printf("Frame phase:\n");
    RUN_TEST(test_frame_phase_assignment);
    RUN_TEST(test_frame_phase_incomplete_codon);

    printf("\nAmino acid translation:\n");
    RUN_TEST(test_amino_acid_translation);
    RUN_TEST(test_codon_next_chain);

    printf("\nStop codon counting:\n");
    RUN_TEST(test_stop_codon_count_none);
    RUN_TEST(test_stop_codon_count_one);
    RUN_TEST(test_stop_codon_count_multiple);

    printf("\nPoint mutation retranslation:\n");
    RUN_TEST(test_mutate_retranslates_codon);
    RUN_TEST(test_mutate_creates_stop_codon);
    RUN_TEST(test_mutate_removes_stop_codon);
    RUN_TEST(test_mutate_phase1_retranslates_correct_codon);
    RUN_TEST(test_mutate_phase2_retranslates_correct_codon);

    printf("\nRevert:\n");
    RUN_TEST(test_revert_restores_amino_acid);

    printf("\nInsert / delete (frame propagation):\n");
    RUN_TEST(test_insert_shifts_frame);
    RUN_TEST(test_delete_shifts_frame);

    printf("\nReverse complement:\n");
    RUN_TEST(test_reverse_complement_invalidates_rail);

    printf("\nRebuild:\n");
    RUN_TEST(test_rebuild_codon_rail);

    printf("\nMulti-segment:\n");
    RUN_TEST(test_codon_rail_multi_segment);
    RUN_TEST(test_aseq_has_stop_codons);

    printf("\nV-anchored frame (T0-4):\n");
    RUN_TEST(test_v_anchored_no_head_change);
    RUN_TEST(test_v_anchored_after_5prime_delete);
    RUN_TEST(test_v_anchored_after_5prime_prepend);
    RUN_TEST(test_stop_codon_in_adapter_ignored);
    RUN_TEST(test_no_v_falls_back_not_productive);

    printf("\n%d/%d tests passed.\n", tests_passed, tests_total);
    return tests_passed == tests_total ? 0 : 1;
}
