/**
 * test_bitmap.c — Tests for AlleleBitmap allele call derivation.
 *
 * Verifies that the bitmap-based allele call correctly identifies
 * equivalent alleles under trimming, mutation, and deletion.
 */

#include "genairr/genairr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ── Test framework ──────────────────────────────────────────── */

static int tests_run = 0;
static int tests_passed = 0;

#define ASSERT(cond, msg) do { \
    if (!(cond)) { \
        printf("  FAIL: %s (line %d)\n    %s\n", __func__, __LINE__, msg); \
        return 0; \
    } \
} while (0)

#define RUN(fn) do { \
    tests_run++; \
    printf("  %-52s", #fn); \
    if (fn()) { tests_passed++; printf("PASS\n"); } \
    else { printf("FAIL\n"); } \
} while (0)

/* ═══════════════════════════════════════════════════════════════
 *  Basic bitmap construction and lookup
 * ═══════════════════════════════════════════════════════════════ */

static int test_bitmap_build_basic(void) {
    AllelePool pool = allele_pool_create(4);

    Allele a1 = { .name = "V1*01", .seq = "ATGC", .length = 4,
                  .segment_type = SEG_V };
    Allele a2 = { .name = "V1*02", .seq = "ATGC", .length = 4,
                  .segment_type = SEG_V };
    Allele a3 = { .name = "V2*01", .seq = "TTGC", .length = 4,
                  .segment_type = SEG_V };

    allele_pool_add(&pool, &a1);
    allele_pool_add(&pool, &a2);
    allele_pool_add(&pool, &a3);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    ASSERT(bm.n_alleles == 3, "Should have 3 alleles");
    ASSERT(bm.max_pos == 4, "Max pos should be 4");
    ASSERT(bm.data != NULL, "Data should be allocated");

    /* Position 0, base 'A': alleles 0 and 1 have A, allele 2 has T */
    const uint64_t *bits = allele_bitmap_lookup(&bm, 0, 'A');
    ASSERT(bits != NULL, "Lookup should succeed");
    ASSERT((*bits & 1) != 0, "Allele 0 should match A at pos 0");
    ASSERT((*bits & 2) != 0, "Allele 1 should match A at pos 0");
    ASSERT((*bits & 4) == 0, "Allele 2 should NOT match A at pos 0");

    /* Position 0, base 'T': only allele 2 */
    bits = allele_bitmap_lookup(&bm, 0, 'T');
    ASSERT(bits != NULL, "Lookup should succeed");
    ASSERT((*bits & 4) != 0, "Allele 2 should match T at pos 0");
    ASSERT((*bits & 3) == 0, "Alleles 0,1 should NOT match T at pos 0");

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 1;
}

/* ═══════════════════════════════════════════════════════════════
 *  Allele call derivation — identical alleles
 * ═══════════════════════════════════════════════════════════════ */

static int test_derive_identical_alleles(void) {
    AllelePool pool = allele_pool_create(4);

    Allele a1 = { .name = "V1*01", .seq = "ATGCATGCATGC", .length = 12,
                  .anchor = 10, .segment_type = SEG_V };
    Allele a2 = { .name = "V1*02", .seq = "ATGCATGCATGC", .length = 12,
                  .anchor = 10, .segment_type = SEG_V };

    allele_pool_add(&pool, &a1);
    allele_pool_add(&pool, &a2);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    /* Build sequence from allele 0 */
    ASeq seq;
    aseq_init(&seq);
    aseq_append_segment(&seq, a1.seq, a1.length, SEG_V, 0, a1.anchor);

    AlleleCallResult result;
    allele_call_derive(&bm, &seq, SEG_V, &result, NULL);

    ASSERT(result.count == 2, "Both alleles should tie");
    ASSERT(result.best_score == 12, "All 12 positions match");

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 1;
}

/* ═══════════════════════════════════════════════════════════════
 *  Allele call — single SNP difference
 * ═══════════════════════════════════════════════════════════════ */

static int test_derive_snp_difference(void) {
    AllelePool pool = allele_pool_create(4);

    /*                                     v  (pos 3 differs) */
    Allele a1 = { .name = "V1*01", .seq = "ATGAATGCATGC", .length = 12,
                  .anchor = 10, .segment_type = SEG_V };
    Allele a2 = { .name = "V1*02", .seq = "ATGCATGCATGC", .length = 12,
                  .anchor = 10, .segment_type = SEG_V };

    allele_pool_add(&pool, &a1);
    allele_pool_add(&pool, &a2);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    /* Build sequence from allele 0 (has 'A' at pos 3) */
    ASeq seq;
    aseq_init(&seq);
    aseq_append_segment(&seq, a1.seq, a1.length, SEG_V, 0, a1.anchor);

    AlleleCallResult result;
    allele_call_derive(&bm, &seq, SEG_V, &result, NULL);

    /* Allele 0 should score 12 (all match), allele 1 scores 11 (pos 3 differs) */
    ASSERT(result.count == 1, "Only allele 0 should win");
    ASSERT(result.indices[0] == 0, "Allele 0 should be the best");
    ASSERT(result.best_score == 12, "Best score should be 12");

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 1;
}

/* ═══════════════════════════════════════════════════════════════
 *  Allele call — trimming makes alleles equivalent
 * ═══════════════════════════════════════════════════════════════ */

static int test_derive_trimming_creates_equivalence(void) {
    AllelePool pool = allele_pool_create(4);

    /* Alleles differ only at position 11 (the last position) */
    Allele a1 = { .name = "V1*01", .seq = "ATGCATGCATGA", .length = 12,
                  .anchor = 8, .segment_type = SEG_V };
    Allele a2 = { .name = "V1*02", .seq = "ATGCATGCATGC", .length = 12,
                  .anchor = 8, .segment_type = SEG_V };

    allele_pool_add(&pool, &a1);
    allele_pool_add(&pool, &a2);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    /* Build sequence from allele 0, but with 3' trim of 2
     * (positions 10,11 trimmed — the differing position is gone) */
    ASeq seq;
    aseq_init(&seq);
    aseq_append_segment(&seq, a1.seq, 10, SEG_V, 0, a1.anchor);

    AlleleCallResult result;
    allele_call_derive(&bm, &seq, SEG_V, &result, NULL);

    /* Both alleles should match on the retained 10 positions */
    ASSERT(result.count == 2, "Both alleles should tie after trim");
    ASSERT(result.best_score == 10, "Score = 10 retained positions");

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 1;
}

/* ═══════════════════════════════════════════════════════════════
 *  Allele call — mutations CHANGE the call (god-aligner semantics)
 * ═══════════════════════════════════════════════════════════════ */

static int test_derive_mutations_change_call(void) {
    AllelePool pool = allele_pool_create(4);

    /* Alleles differ only at position 0: a1 has 'A', a2 has 'T' */
    Allele a1 = { .name = "V1*01", .seq = "ATGCATGCATGC", .length = 12,
                  .anchor = 10, .segment_type = SEG_V };
    Allele a2 = { .name = "V2*01", .seq = "TTGCATGCATGC", .length = 12,
                  .anchor = 10, .segment_type = SEG_V };

    allele_pool_add(&pool, &a1);
    allele_pool_add(&pool, &a2);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    /* Build from allele 0, then mutate position 0 from A→T
     * This makes the sequence match allele 2 perfectly */
    ASeq seq;
    aseq_init(&seq);
    aseq_append_segment(&seq, a1.seq, a1.length, SEG_V, 0, a1.anchor);

    /* Mutate A→T at pos 0: now current = 'T', germline = 'A' */
    aseq_mutate(&seq, seq.head, 'T', NUC_FLAG_MUTATED);

    AlleleCallResult result;
    allele_call_derive(&bm, &seq, SEG_V, &result, NULL);

    /* God-aligner: current[0] = T, which matches allele 2 at pos 0.
     * Allele 2 scores 12 (all positions match), allele 1 scores 11.
     * The mutation CHANGED the call from allele 1 to allele 2. */
    ASSERT(result.count == 1, "Only allele 2 (T at pos 0) should win");
    ASSERT(result.indices[0] == 1, "Allele 2 is the best match");
    ASSERT(result.best_score == 12, "Score = 12 (all current bases match)");

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 1;
}

/* ═══════════════════════════════════════════════════════════════
 *  Allele call — mutation to non-allele base creates tie
 * ═══════════════════════════════════════════════════════════════ */

static int test_derive_mutation_creates_tie(void) {
    AllelePool pool = allele_pool_create(4);

    /* Alleles differ at pos 0: a1='A', a2='T' */
    Allele a1 = { .name = "V1*01", .seq = "ATGCATGCATGC", .length = 12,
                  .anchor = 10, .segment_type = SEG_V };
    Allele a2 = { .name = "V2*01", .seq = "TTGCATGCATGC", .length = 12,
                  .anchor = 10, .segment_type = SEG_V };

    allele_pool_add(&pool, &a1);
    allele_pool_add(&pool, &a2);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    ASeq seq;
    aseq_init(&seq);
    aseq_append_segment(&seq, a1.seq, a1.length, SEG_V, 0, a1.anchor);

    /* Mutate A→G at pos 0: 'G' matches neither allele at pos 0 */
    aseq_mutate(&seq, seq.head, 'G', NUC_FLAG_MUTATED);

    AlleleCallResult result;
    allele_call_derive(&bm, &seq, SEG_V, &result, NULL);

    /* God-aligner: current[0] = G, matches neither allele at pos 0.
     * Both alleles now score 11 (positions 1-11 are identical).
     * The distinguishing evidence at pos 0 is lost. */
    ASSERT(result.count == 2, "Both alleles should tie");
    ASSERT(result.best_score == 11, "Score = 11 (pos 0 lost)");

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 1;
}

/* ═══════════════════════════════════════════════════════════════
 *  Allele call — NP nodes don't contribute
 * ═══════════════════════════════════════════════════════════════ */

static int test_derive_skips_np_nodes(void) {
    AllelePool pool = allele_pool_create(4);

    Allele a1 = { .name = "V1*01", .seq = "ATGC", .length = 4,
                  .anchor = 2, .segment_type = SEG_V };
    allele_pool_add(&pool, &a1);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    ASeq seq;
    aseq_init(&seq);
    aseq_append_segment(&seq, a1.seq, a1.length, SEG_V, 0, a1.anchor);
    /* Add NP1 — should not affect V call */
    aseq_append_np(&seq, "AAAA", 4, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);

    AlleleCallResult result;
    allele_call_derive(&bm, &seq, SEG_V, &result, NULL);

    ASSERT(result.total_positions == 4, "Only 4 V positions counted");
    ASSERT(result.count == 1, "One allele matches");

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 1;
}

/* ═══════════════════════════════════════════════════════════════
 *  Allele call — deletion removes that position's vote
 * ═══════════════════════════════════════════════════════════════ */

static int test_derive_deletion_removes_vote(void) {
    AllelePool pool = allele_pool_create(4);

    /* Differ at pos 0: a1 has A, a2 has T */
    Allele a1 = { .name = "V1*01", .seq = "ATGCATGC", .length = 8,
                  .anchor = 6, .segment_type = SEG_V };
    Allele a2 = { .name = "V1*02", .seq = "TTGCATGC", .length = 8,
                  .anchor = 6, .segment_type = SEG_V };

    allele_pool_add(&pool, &a1);
    allele_pool_add(&pool, &a2);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    /* Build from allele 0 */
    ASeq seq;
    aseq_init(&seq);
    aseq_append_segment(&seq, a1.seq, a1.length, SEG_V, 0, a1.anchor);

    /* Before deletion: allele 0 wins with 8, allele 1 has 7 */
    AlleleCallResult result;
    allele_call_derive(&bm, &seq, SEG_V, &result, NULL);
    ASSERT(result.count == 1, "Before deletion: allele 0 wins");
    ASSERT(result.indices[0] == 0, "Allele 0 is the best");

    /* Delete the head node (pos 0, base A) — the distinguishing position */
    aseq_delete(&seq, seq.head);

    allele_call_derive(&bm, &seq, SEG_V, &result, NULL);

    /* Now both alleles match on positions 1-7 (7 positions each) */
    ASSERT(result.count == 2, "After deleting distinguishing position, both tie");
    ASSERT(result.best_score == 7, "Score = 7 remaining positions");

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 1;
}

/* ═══════════════════════════════════════════════════════════════
 *  Short-D detection
 * ═══════════════════════════════════════════════════════════════ */

static int test_short_d_detection(void) {
    AllelePool pool = allele_pool_create(4);

    Allele d1 = { .name = "D1*01", .seq = "GGTATAACTGGAAC", .length = 14,
                  .segment_type = SEG_D };
    Allele d2 = { .name = "D2*01", .seq = "AGGATATTGTACT", .length = 13,
                  .segment_type = SEG_D };

    allele_pool_add(&pool, &d1);
    allele_pool_add(&pool, &d2);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    /* Build a very short D segment (2 bases) */
    ASeq seq;
    aseq_init(&seq);
    aseq_append_segment(&seq, "AT", 2, SEG_D, 5, -1);

    AlleleCallResult result;
    allele_call_derive(&bm, &seq, SEG_D, &result, NULL);

    ASSERT(allele_call_is_short_d(&result, &bm),
           "2-base D should be Short-D");

    /* Build a long D segment */
    aseq_reset(&seq);
    aseq_append_segment(&seq, d1.seq, d1.length, SEG_D, 0, -1);

    allele_call_derive(&bm, &seq, SEG_D, &result, NULL);

    ASSERT(!allele_call_is_short_d(&result, &bm),
           "Full-length D should NOT be Short-D");

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 1;
}

/* ═══════════════════════════════════════════════════════════════
 *  Format allele call string
 * ═══════════════════════════════════════════════════════════════ */

static int test_format_call_string(void) {
    AllelePool pool = allele_pool_create(4);

    Allele a1 = { .name = "V1*01", .seq = "ATGC", .length = 4,
                  .segment_type = SEG_V };
    Allele a2 = { .name = "V1*02", .seq = "ATGC", .length = 4,
                  .segment_type = SEG_V };
    Allele a3 = { .name = "V2*01", .seq = "TTGC", .length = 4,
                  .segment_type = SEG_V };

    allele_pool_add(&pool, &a1);
    allele_pool_add(&pool, &a2);
    allele_pool_add(&pool, &a3);

    AlleleCallResult result;
    result.count = 2;
    result.indices[0] = 0;
    result.indices[1] = 1;
    result.best_score = 4;

    char buf[256];
    allele_call_format(&result, &pool, 0, buf, sizeof(buf));

    ASSERT(strcmp(buf, "V1*01,V1*02") == 0, "True allele first, then ties");

    /* Format with true_idx = 1 */
    allele_call_format(&result, &pool, 1, buf, sizeof(buf));

    ASSERT(strcmp(buf, "V1*02,V1*01") == 0, "True allele first when idx=1");

    allele_pool_destroy(&pool);
    return 1;
}

/* ═══════════════════════════════════════════════════════════════
 *  AlleleCorrectionSet lifecycle
 * ═══════════════════════════════════════════════════════════════ */

static int test_correction_set_lifecycle(void) {
    SimConfig cfg;
    sim_config_init(&cfg, CHAIN_IGH);

    Allele v = { .name = "V1*01", .seq = "ATGCATGC", .length = 8,
                 .segment_type = SEG_V };
    Allele d = { .name = "D1*01", .seq = "GGTA", .length = 4,
                 .segment_type = SEG_D };
    Allele j = { .name = "J1*01", .seq = "TGGGGC", .length = 6,
                 .segment_type = SEG_J };

    allele_pool_add(&cfg.v_alleles, &v);
    allele_pool_add(&cfg.d_alleles, &d);
    allele_pool_add(&cfg.j_alleles, &j);

    AlleleCorrectionSet set = allele_correction_set_build(&cfg);

    ASSERT(set.v_bitmap.n_alleles == 1, "V bitmap has 1 allele");
    ASSERT(set.d_bitmap.n_alleles == 1, "D bitmap has 1 allele");
    ASSERT(set.j_bitmap.n_alleles == 1, "J bitmap has 1 allele");

    allele_correction_set_destroy(&set);
    sim_config_destroy(&cfg);
    return 1;
}

/* ═══════════════════════════════════════════════════════════════
 *  Real-world scenario: V trimming with multiple similar alleles
 * ═══════════════════════════════════════════════════════════════ */

static int test_realistic_v_trimming(void) {
    AllelePool pool = allele_pool_create(4);

    /* Three V alleles: differ at pos 95 and 98 */
    char seq1[128] = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAAAGCAAGC";
    char seq2[128] = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAAAGCATGC";
    char seq3[128] = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCTTGGCATGC";

    Allele a1 = { .name = "V1*01", .gene = "V1", .length = 100,
                  .anchor = 90, .segment_type = SEG_V };
    memcpy(a1.seq, seq1, 100);

    Allele a2 = { .name = "V1*02", .gene = "V1", .length = 100,
                  .anchor = 90, .segment_type = SEG_V };
    memcpy(a2.seq, seq2, 100);

    Allele a3 = { .name = "V2*01", .gene = "V2", .length = 100,
                  .anchor = 90, .segment_type = SEG_V };
    memcpy(a3.seq, seq3, 100);

    allele_pool_add(&pool, &a1);
    allele_pool_add(&pool, &a2);
    allele_pool_add(&pool, &a3);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    /* Simulate: use allele 0, trim 3 from 3' end (remove positions 97-99)
     * After trim, diff at pos 98 is gone. Only pos 95 differs between a1/a3 */
    ASeq seq;
    aseq_init(&seq);
    aseq_append_segment(&seq, a1.seq, 97, SEG_V, 0, a1.anchor);

    AlleleCallResult result;
    allele_call_derive(&bm, &seq, SEG_V, &result, NULL);

    /* a1 and a2 differ at pos 98 (trimmed) → both score 97
     * a3 differs at pos 95 (A vs T) and pos 96 (A vs G) → scores 95
     * So a1 and a2 should tie at 97 */
    ASSERT(result.count == 2, "V1*01 and V1*02 should tie after trim");
    ASSERT(result.best_score == 97, "Best score should be 97");

    /* Verify a3 is NOT in the result */
    bool a3_found = false;
    for (int i = 0; i < result.count; i++) {
        if (result.indices[i] == 2) a3_found = true;
    }
    ASSERT(!a3_found, "V2*01 should NOT be in the call (differs at pos 95)");

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 1;
}

/* ═══════════════════════════════════════════════════════════════
 *  Mutation + trimming combined
 * ═══════════════════════════════════════════════════════════════ */

static int test_mutation_and_trim_combined(void) {
    AllelePool pool = allele_pool_create(4);

    /* Two alleles differ at pos 2 */
    Allele a1 = { .name = "V1*01", .seq = "AAGCATGC", .length = 8,
                  .anchor = 6, .segment_type = SEG_V };
    Allele a2 = { .name = "V1*02", .seq = "AATCATGC", .length = 8,
                  .anchor = 6, .segment_type = SEG_V };

    allele_pool_add(&pool, &a1);
    allele_pool_add(&pool, &a2);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    /* Build from allele 0 with 5' trim of 3 (removes positions 0,1,2) */
    ASeq seq;
    aseq_init(&seq);
    aseq_append_segment(&seq, a1.seq + 3, 5, SEG_V, 3, a1.anchor);

    /* Also mutate position 4 (germline_pos=4) from A→G */
    Nuc *n = seq.head;
    n = n->next;  /* pos 4 (0-indexed from trim start at pos 3) */
    aseq_mutate(&seq, n, 'G', NUC_FLAG_MUTATED);

    AlleleCallResult result;
    allele_call_derive(&bm, &seq, SEG_V, &result, NULL);

    /* Retained positions: 3,4,5,6,7 (5 positions)
     * Differing position 2 is trimmed → both alleles match on retained.
     * GOD-ALIGNER: Mutation at pos 4 (A→G) uses current base 'G'.
     * Both alleles have 'A' at pos 4, so neither matches there.
     * Both score 4 (positions 3,5,6,7 match, pos 4 lost). */
    ASSERT(result.count == 2, "Both alleles should tie (diff pos trimmed, mutated pos lost)");
    ASSERT(result.best_score == 4, "Score = 4 (pos 4 mutated, no allele has G there)");

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 1;
}

/* ═══════════════════════════════════════════════════════════════ */

int main(void) {
    printf("test_bitmap: AlleleBitmap allele call derivation\n");
    printf("═══════════════════════════════════════════════════\n\n");

    printf("Construction & lookup:\n");
    RUN(test_bitmap_build_basic);

    printf("\nAllele call derivation:\n");
    RUN(test_derive_identical_alleles);
    RUN(test_derive_snp_difference);
    RUN(test_derive_trimming_creates_equivalence);
    RUN(test_derive_mutations_change_call);
    RUN(test_derive_mutation_creates_tie);
    RUN(test_derive_skips_np_nodes);
    RUN(test_derive_deletion_removes_vote);

    printf("\nShort-D:\n");
    RUN(test_short_d_detection);

    printf("\nFormatting:\n");
    RUN(test_format_call_string);

    printf("\nLifecycle:\n");
    RUN(test_correction_set_lifecycle);

    printf("\nRealistic scenarios:\n");
    RUN(test_realistic_v_trimming);
    RUN(test_mutation_and_trim_combined);

    printf("\n%d/%d tests passed.\n", tests_passed, tests_run);
    return tests_passed == tests_run ? 0 : 1;
}
