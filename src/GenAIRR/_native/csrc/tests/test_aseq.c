/**
 * test_aseq.c — Tests for the ASeq annotated sequence data structure.
 *
 * Minimal test harness (no external dependencies). Each test function
 * returns 0 on success, 1 on failure.
 */

#include "genairr/genairr.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>

#define TEST(name) \
    do { \
        printf("  %-50s", #name); \
        int _r = name(); \
        printf("%s\n", _r == 0 ? "PASS" : "FAIL"); \
        if (_r) failures++; \
        total++; \
    } while (0)

/* ── Test: basic init ─────────────────────────────────────────── */

static int test_init(void) {
    ASeq seq;
    aseq_init(&seq);
    if (seq.head != NULL) return 1;
    if (seq.tail != NULL) return 1;
    if (seq.length != 0)  return 1;
    if (seq.pool_used != 0) return 1;
    return 0;
}

/* ── Test: append a V segment ─────────────────────────────────── */

static int test_append_v_segment(void) {
    ASeq seq;
    aseq_init(&seq);

    const char *v_seq = "ATCGATCG";
    int len = strlen(v_seq);
    int ret = aseq_append_segment(&seq, v_seq, len, SEG_V, 0, -1);

    if (ret != len) return 1;
    if (seq.length != len) return 1;
    if (seq.head == NULL || seq.tail == NULL) return 1;

    /* Verify sequence content */
    char buf[64];
    aseq_to_string(&seq, buf, sizeof(buf));
    if (strcmp(buf, v_seq) != 0) return 1;

    /* Verify all nodes are tagged V */
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->segment != SEG_V) return 1;
    }

    /* Verify segment cache */
    if (seq.seg_first[SEG_V] != seq.head) return 1;
    if (seq.seg_last[SEG_V] != seq.tail) return 1;

    return 0;
}

/* ── Test: append multiple segments ───────────────────────────── */

static int test_append_multiple_segments(void) {
    ASeq seq;
    aseq_init(&seq);

    aseq_append_segment(&seq, "AAAA", 4, SEG_V, 0, -1);
    aseq_append_np(&seq, "CC", 2, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "GG", 2, SEG_D, 0, -1);
    aseq_append_np(&seq, "TT", 2, SEG_NP2, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "CCCC", 4, SEG_J, 0, -1);

    if (seq.length != 14) return 1;

    char buf[64];
    aseq_to_string(&seq, buf, sizeof(buf));
    if (strcmp(buf, "AAAACCGGTTCCCC") != 0) return 1;

    /* Verify segment lengths */
    if (aseq_segment_length(&seq, SEG_V) != 4) return 1;
    if (aseq_segment_length(&seq, SEG_NP1) != 2) return 1;
    if (aseq_segment_length(&seq, SEG_D) != 2) return 1;
    if (aseq_segment_length(&seq, SEG_NP2) != 2) return 1;
    if (aseq_segment_length(&seq, SEG_J) != 4) return 1;

    return 0;
}

/* ── Test: anchor detection ───────────────────────────────────── */

static int test_anchor(void) {
    ASeq seq;
    aseq_init(&seq);

    /* V allele with anchor at germline position 3 */
    aseq_append_segment(&seq, "ATCGATCG", 8, SEG_V, 0, 3);

    Nuc *anchor = aseq_find_anchor(&seq, SEG_V);
    if (anchor == NULL) return 1;
    if (anchor->germline_pos != 3) return 1;
    if (!(anchor->flags & NUC_FLAG_ANCHOR)) return 1;

    return 0;
}

/* ── Test: insert after ───────────────────────────────────────── */

static int test_insert_after(void) {
    ASeq seq;
    aseq_init(&seq);

    aseq_append_segment(&seq, "ACGT", 4, SEG_V, 0, -1);
    /* Insert 'X' after the second node (C) */
    Nuc *second = seq.head->next;
    Nuc *inserted = aseq_insert_after(&seq, second, 'N', SEG_V, NUC_FLAG_INDEL_INS);

    if (inserted == NULL) return 1;
    if (seq.length != 5) return 1;

    char buf[64];
    aseq_to_string(&seq, buf, sizeof(buf));
    if (strcmp(buf, "ACNGT") != 0) return 1;

    /* Verify the insertion's prev/next links */
    if (inserted->prev != second) return 1;
    if (inserted->next->current != 'G') return 1;

    return 0;
}

/* ── Test: delete node ────────────────────────────────────────── */

static int test_delete(void) {
    ASeq seq;
    aseq_init(&seq);

    aseq_append_segment(&seq, "ACGT", 4, SEG_V, 0, -1);
    /* Delete the third node (G) */
    Nuc *third = seq.head->next->next;
    aseq_delete(&seq, third);

    if (seq.length != 3) return 1;

    char buf[64];
    aseq_to_string(&seq, buf, sizeof(buf));
    if (strcmp(buf, "ACT") != 0) return 1;

    return 0;
}

/* ── Test: delete head ────────────────────────────────────────── */

static int test_delete_head(void) {
    ASeq seq;
    aseq_init(&seq);

    aseq_append_segment(&seq, "ACGT", 4, SEG_V, 0, -1);
    aseq_delete(&seq, seq.head);

    if (seq.length != 3) return 1;

    char buf[64];
    aseq_to_string(&seq, buf, sizeof(buf));
    if (strcmp(buf, "CGT") != 0) return 1;
    if (seq.head->prev != NULL) return 1;

    return 0;
}

/* ── Test: delete tail ────────────────────────────────────────── */

static int test_delete_tail(void) {
    ASeq seq;
    aseq_init(&seq);

    aseq_append_segment(&seq, "ACGT", 4, SEG_V, 0, -1);
    aseq_delete(&seq, seq.tail);

    if (seq.length != 3) return 1;

    char buf[64];
    aseq_to_string(&seq, buf, sizeof(buf));
    if (strcmp(buf, "ACG") != 0) return 1;
    if (seq.tail->next != NULL) return 1;

    return 0;
}

/* ── Test: mutation preserves germline ────────────────────────── */

static int test_mutation(void) {
    ASeq seq;
    aseq_init(&seq);

    aseq_append_segment(&seq, "AAAA", 4, SEG_V, 0, -1);

    /* Mutate second node A -> G */
    Nuc *n = seq.head->next;
    aseq_mutate(&seq, n, 'G', NUC_FLAG_MUTATED);

    if (n->current != 'G') return 1;
    if (n->germline != 'A') return 1;
    if (!(n->flags & NUC_FLAG_MUTATED)) return 1;

    char buf[64];
    aseq_to_string(&seq, buf, sizeof(buf));
    if (strcmp(buf, "AGAA") != 0) return 1;

    return 0;
}

/* ── Test: reverse complement ─────────────────────────────────── */

static int test_reverse_complement(void) {
    ASeq seq;
    aseq_init(&seq);

    aseq_append_segment(&seq, "ATCG", 4, SEG_V, 0, -1);
    aseq_reverse_complement(&seq);

    char buf[64];
    aseq_to_string(&seq, buf, sizeof(buf));
    /* RC of ATCG is CGAT */
    if (strcmp(buf, "CGAT") != 0) return 1;

    /* Verify head/tail swapped correctly */
    if (seq.head->current != 'C') return 1;
    if (seq.tail->current != 'T') return 1;
    /* Germline rail must mirror current after RC for aligned serialization */
    if (seq.head->germline != 'C') return 1;
    if (seq.tail->germline != 'T') return 1;

    return 0;
}

/* ── Test: stats computation ──────────────────────────────────── */

static int test_stats(void) {
    ASeq seq;
    aseq_init(&seq);

    aseq_append_segment(&seq, "AAAA", 4, SEG_V, 0, -1);
    aseq_append_np(&seq, "CC", 2, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "GGGG", 4, SEG_J, 0, -1);

    /* Mutate 2 V bases */
    aseq_mutate(&seq, seq.head, 'T', NUC_FLAG_MUTATED);
    aseq_mutate(&seq, seq.head->next, 'G', NUC_FLAG_MUTATED);

    ASeqStats stats = aseq_stats(&seq);

    if (stats.total_germline != 8) return 1;  /* 4 V + 4 J */
    if (stats.mutated != 2) return 1;
    if (stats.mutation_rate < 0.24 || stats.mutation_rate > 0.26) return 1; /* 2/8 = 0.25 */

    return 0;
}

/* ── Test: segment_to_string ──────────────────────────────────── */

static int test_segment_to_string(void) {
    ASeq seq;
    aseq_init(&seq);

    aseq_append_segment(&seq, "AAAA", 4, SEG_V, 0, -1);
    aseq_append_np(&seq, "CC", 2, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "TTTT", 4, SEG_J, 0, -1);

    char buf[64];

    aseq_segment_to_string(&seq, SEG_V, buf, sizeof(buf));
    if (strcmp(buf, "AAAA") != 0) return 1;

    aseq_segment_to_string(&seq, SEG_NP1, buf, sizeof(buf));
    if (strcmp(buf, "CC") != 0) return 1;

    aseq_segment_to_string(&seq, SEG_J, buf, sizeof(buf));
    if (strcmp(buf, "TTTT") != 0) return 1;

    /* Non-existent segment */
    aseq_segment_to_string(&seq, SEG_D, buf, sizeof(buf));
    if (strcmp(buf, "") != 0) return 1;

    return 0;
}

/* ── Test: reset and reuse ────────────────────────────────────── */

static int test_reset(void) {
    ASeq seq;
    aseq_init(&seq);

    aseq_append_segment(&seq, "ATCGATCG", 8, SEG_V, 0, -1);
    if (seq.length != 8) return 1;

    aseq_reset(&seq);
    if (seq.length != 0) return 1;
    if (seq.head != NULL) return 1;
    if (seq.pool_used != 0) return 1;

    /* Reuse after reset */
    aseq_append_segment(&seq, "GCTA", 4, SEG_J, 0, -1);
    if (seq.length != 4) return 1;

    char buf[64];
    aseq_to_string(&seq, buf, sizeof(buf));
    if (strcmp(buf, "GCTA") != 0) return 1;

    return 0;
}

/* ── Test: insert does not break positions ────────────────────── */

static int test_insert_preserves_structure(void) {
    ASeq seq;
    aseq_init(&seq);

    /* Build V + NP1 + J */
    aseq_append_segment(&seq, "AAAA", 4, SEG_V, 0, 2);
    aseq_append_np(&seq, "CC", 2, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "TTTT", 4, SEG_J, 0, 1);

    /* Insert an indel in V segment */
    Nuc *v_second = seq.head->next;
    aseq_insert_after(&seq, v_second, 'G', SEG_V, NUC_FLAG_INDEL_INS);

    if (seq.length != 11) return 1;

    /* V anchor should still be found and correct */
    Nuc *v_anchor = aseq_find_anchor(&seq, SEG_V);
    if (v_anchor == NULL) return 1;
    if (v_anchor->germline_pos != 2) return 1;

    /* J anchor should still be found */
    Nuc *j_anchor = aseq_find_anchor(&seq, SEG_J);
    if (j_anchor == NULL) return 1;
    if (j_anchor->germline_pos != 1) return 1;

    /* NP1 segment should be unchanged */
    if (aseq_segment_length(&seq, SEG_NP1) != 2) return 1;

    return 0;
}

/* ── Test: insert_after returns NULL on pool exhaustion (T0-6) ── */

/* When pool_used == GENAIRR_MAX_SEQ_LEN, no node can be allocated.
 * aseq_insert_after must return NULL, leaving the sequence
 * unchanged. Callers (insert_indels, long_read_errors) rely on this
 * to avoid over-counting n_insertions. */
static int test_insert_after_pool_exhausted(void) {
    ASeq seq;
    aseq_init(&seq);

    /* Fill the pool to exactly the cap. We append a large germline
     * segment at the limit. */
    char *big = malloc(GENAIRR_MAX_SEQ_LEN);
    if (!big) return 1;
    for (int i = 0; i < GENAIRR_MAX_SEQ_LEN; i++) big[i] = 'A';
    aseq_append_segment(&seq, big, GENAIRR_MAX_SEQ_LEN, SEG_V, 0, -1);
    free(big);

    if (seq.pool_used != GENAIRR_MAX_SEQ_LEN) return 1;
    int len_before = seq.length;

    /* Try to insert one more node — must fail. */
    Nuc *inserted = aseq_insert_after(&seq, seq.head, 'C', SEG_V,
                                       NUC_FLAG_INDEL_INS);
    if (inserted != NULL) return 1;
    if (seq.length != len_before) return 1;

    return 0;
}

/* ── Main ─────────────────────────────────────────────────────── */

int main(void) {
    int failures = 0;
    int total = 0;

    printf("=== ASeq Tests ===\n");

    TEST(test_init);
    TEST(test_append_v_segment);
    TEST(test_append_multiple_segments);
    TEST(test_anchor);
    TEST(test_insert_after);
    TEST(test_insert_after_pool_exhausted);
    TEST(test_delete);
    TEST(test_delete_head);
    TEST(test_delete_tail);
    TEST(test_mutation);
    TEST(test_reverse_complement);
    TEST(test_stats);
    TEST(test_segment_to_string);
    TEST(test_reset);
    TEST(test_insert_preserves_structure);

    printf("\n%d/%d tests passed.\n", total - failures, total);
    return failures > 0 ? 1 : 0;
}
