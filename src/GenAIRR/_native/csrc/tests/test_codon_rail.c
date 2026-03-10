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

    printf("\n%d/%d tests passed.\n", tests_passed, tests_total);
    return tests_passed == tests_total ? 0 : 1;
}
