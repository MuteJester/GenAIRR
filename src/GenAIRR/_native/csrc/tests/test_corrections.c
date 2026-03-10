/**
 * test_corrections.c — Tests for bitmap-based allele call correction.
 *
 * Validates that AlleleBitmap correctly derives allele calls from
 * ASeq segment state, replacing the old TrimCorrectionMap approach.
 */

#include "genairr/genairr.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define TEST(name) \
    do { \
        printf("  %-55s", #name); \
        int _r = name(); \
        printf("%s\n", _r == 0 ? "PASS" : "FAIL"); \
        if (_r) failures++; \
        total++; \
    } while (0)

/* ── Helper: create allele with known sequence ─────────────────── */

static Allele make_allele(const char *name, const char *seq,
                          int anchor, Segment seg) {
    Allele a = {0};
    strncpy(a.name, name, sizeof(a.name) - 1);
    strncpy(a.gene, name, sizeof(a.gene) - 1);
    strncpy(a.family, name, sizeof(a.family) - 1);
    strncpy(a.seq, seq, sizeof(a.seq) - 1);
    a.length = strlen(seq);
    a.anchor = anchor;
    a.segment_type = seg;
    return a;
}

/* ══════════════════════════════════════════════════════════════════
 * V Allele Correction Tests (bitmap-based)
 * ══════════════════════════════════════════════════════════════════ */

/**
 * Test: Two V alleles differ at positions 9,10,11 (last 3 bases).
 * Full sequence → only true allele. Trimmed by 3 → both tie.
 */
static int test_v_correction_distinct(void) {
    AllelePool pool = allele_pool_create(4);
    Allele v1 = make_allele("V1*01", "ATCGATCGATCG", 7, SEG_V);
    Allele v2 = make_allele("V2*01", "ATCGATCGACCC", 7, SEG_V);
    allele_pool_add(&pool, &v1);
    allele_pool_add(&pool, &v2);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    /* Full V (12bp): alleles differ → only true matches best */
    {
        ASeq seq;
        aseq_init(&seq);
        aseq_append_segment(&seq, "ATCGATCGATCG", 12, SEG_V, 0, 7);

        AlleleCallResult r;
        allele_call_derive(&bm, &seq, SEG_V, &r, NULL);

        /* V1 should score 12 (all positions match), V2 should score 9 */
        if (r.count != 1) {
            fprintf(stderr, "Expected 1 best at full, got %d\n", r.count);
            allele_bitmap_destroy(&bm);
            allele_pool_destroy(&pool);
            return 1;
        }
        if (r.indices[0] != 0) {  /* V1*01 is index 0 */
            fprintf(stderr, "Expected V1*01 (idx 0), got idx %d\n", r.indices[0]);
            allele_bitmap_destroy(&bm);
            allele_pool_destroy(&pool);
            return 1;
        }
    }

    /* Trimmed by 3 (only first 9bp remain): alleles are identical → tie */
    {
        ASeq seq;
        aseq_init(&seq);
        aseq_append_segment(&seq, "ATCGATCGA", 9, SEG_V, 0, 7);

        AlleleCallResult r;
        allele_call_derive(&bm, &seq, SEG_V, &r, NULL);

        if (r.count != 2) {
            fprintf(stderr, "Expected 2 tied at trim_3=3, got %d\n", r.count);
            allele_bitmap_destroy(&bm);
            allele_pool_destroy(&pool);
            return 1;
        }
    }

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 0;
}

/**
 * Test: Three identical V alleles → always tie at any trimming.
 */
static int test_v_correction_identical(void) {
    AllelePool pool = allele_pool_create(4);
    Allele v1 = make_allele("V1*01", "ATCGATCGATCG", 7, SEG_V);
    Allele v2 = make_allele("V1*02", "ATCGATCGATCG", 7, SEG_V);
    Allele v3 = make_allele("V1*03", "ATCGATCGATCG", 7, SEG_V);
    allele_pool_add(&pool, &v1);
    allele_pool_add(&pool, &v2);
    allele_pool_add(&pool, &v3);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    ASeq seq;
    aseq_init(&seq);
    aseq_append_segment(&seq, "ATCGATCGATCG", 12, SEG_V, 0, 7);

    AlleleCallResult r;
    allele_call_derive(&bm, &seq, SEG_V, &r, NULL);

    if (r.count != 3) {
        fprintf(stderr, "Expected 3 tied for identical alleles, got %d\n", r.count);
        allele_bitmap_destroy(&bm);
        allele_pool_destroy(&pool);
        return 1;
    }

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 0;
}

/* ══════════════════════════════════════════════════════════════════
 * J Allele Correction Tests
 * ══════════════════════════════════════════════════════════════════ */

/**
 * Test: Two J alleles differ at positions 0,1.
 * Full J → distinct. Trimmed 5' by 2 → tie.
 */
static int test_j_correction_distinct(void) {
    AllelePool pool = allele_pool_create(4);
    Allele j1 = make_allele("J1*01", "AATTTGGCCAAGGG", 2, SEG_J);
    Allele j2 = make_allele("J2*01", "CCTTTGGCCAAGGG", 2, SEG_J);
    allele_pool_add(&pool, &j1);
    allele_pool_add(&pool, &j2);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    /* Full J → differ at pos 0,1 → only true matches best */
    {
        ASeq seq;
        aseq_init(&seq);
        aseq_append_segment(&seq, "AATTTGGCCAAGGG", 14, SEG_J, 0, 2);

        AlleleCallResult r;
        allele_call_derive(&bm, &seq, SEG_J, &r, NULL);
        if (r.count != 1) {
            fprintf(stderr, "Expected 1 at full J, got %d\n", r.count);
            allele_bitmap_destroy(&bm);
            allele_pool_destroy(&pool);
            return 1;
        }
    }

    /* Trimmed 5' by 2: remaining "TTTGGCCAAGGG" is identical → tie */
    {
        ASeq seq;
        aseq_init(&seq);
        aseq_append_segment(&seq, "TTTGGCCAAGGG", 12, SEG_J, 2, 2);

        AlleleCallResult r;
        allele_call_derive(&bm, &seq, SEG_J, &r, NULL);
        if (r.count != 2) {
            fprintf(stderr, "Expected 2 tied at trim_5=2, got %d\n", r.count);
            allele_bitmap_destroy(&bm);
            allele_pool_destroy(&pool);
            return 1;
        }
    }

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 0;
}

/* ══════════════════════════════════════════════════════════════════
 * D Allele Correction Tests
 * ══════════════════════════════════════════════════════════════════ */

/**
 * Test: Two D alleles share 10bp core "GTATTACTGT", flanked differently.
 * Full D → distinct. Trimmed 2+2 → tie.
 */
static int test_d_correction(void) {
    AllelePool pool = allele_pool_create(4);
    Allele d1 = make_allele("D1*01", "AAGTATTACTGTAA", 0, SEG_D);
    Allele d2 = make_allele("D2*01", "CCGTATTACTGTCC", 0, SEG_D);
    allele_pool_add(&pool, &d1);
    allele_pool_add(&pool, &d2);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    /* Full D → distinct */
    {
        ASeq seq;
        aseq_init(&seq);
        aseq_append_segment(&seq, "AAGTATTACTGTAA", 14, SEG_D, 0, -1);

        AlleleCallResult r;
        allele_call_derive(&bm, &seq, SEG_D, &r, NULL);
        if (r.count != 1) {
            fprintf(stderr, "Expected 1 at full D, got %d\n", r.count);
            allele_bitmap_destroy(&bm);
            allele_pool_destroy(&pool);
            return 1;
        }
    }

    /* Trimmed 2+2: only core "GTATTACTGT" remains → tie */
    {
        ASeq seq;
        aseq_init(&seq);
        aseq_append_segment(&seq, "GTATTACTGT", 10, SEG_D, 2, -1);

        AlleleCallResult r;
        allele_call_derive(&bm, &seq, SEG_D, &r, NULL);
        if (r.count != 2) {
            fprintf(stderr, "Expected 2 at D trim 2+2, got %d\n", r.count);
            allele_bitmap_destroy(&bm);
            allele_pool_destroy(&pool);
            return 1;
        }
    }

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 0;
}

/**
 * Test: Asymmetric D trimming.
 */
static int test_d_correction_asymmetric(void) {
    AllelePool pool = allele_pool_create(4);
    Allele d1 = make_allele("D1*01", "AABCDEFGHAA", 0, SEG_D);
    Allele d2 = make_allele("D2*01", "CCBCDEFGHCC", 0, SEG_D);
    allele_pool_add(&pool, &d1);
    allele_pool_add(&pool, &d2);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    /* trim_5=2, trim_3=0: "BCDEFGHAA" vs "BCDEFGHCC" → differ at end */
    {
        ASeq seq;
        aseq_init(&seq);
        aseq_append_segment(&seq, "BCDEFGHAA", 9, SEG_D, 2, -1);

        AlleleCallResult r;
        allele_call_derive(&bm, &seq, SEG_D, &r, NULL);
        if (r.count != 1) {
            fprintf(stderr, "Expected 1 at D (2,0), got %d\n", r.count);
            allele_bitmap_destroy(&bm);
            allele_pool_destroy(&pool);
            return 1;
        }
    }

    /* trim_5=2, trim_3=2: "BCDEFGH" → identical → tie */
    {
        ASeq seq;
        aseq_init(&seq);
        aseq_append_segment(&seq, "BCDEFGH", 7, SEG_D, 2, -1);

        AlleleCallResult r;
        allele_call_derive(&bm, &seq, SEG_D, &r, NULL);
        if (r.count != 2) {
            fprintf(stderr, "Expected 2 at D (2,2), got %d\n", r.count);
            allele_bitmap_destroy(&bm);
            allele_pool_destroy(&pool);
            return 1;
        }
    }

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 0;
}

/* ══════════════════════════════════════════════════════════════════
 * Short-D Tests (bitmap-based)
 * ══════════════════════════════════════════════════════════════════ */

/**
 * Test: Short-D detection via bitmap.
 */
static int test_short_d_detection(void) {
    AllelePool pool = allele_pool_create(4);
    Allele d1 = make_allele("D1*01", "AATTGGCCAATT", 0, SEG_D);
    Allele d2 = make_allele("D2*01", "AATTGGCCTTAA", 0, SEG_D);
    allele_pool_add(&pool, &d1);
    allele_pool_add(&pool, &d2);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    /* 1bp D → should be short */
    {
        ASeq seq;
        aseq_init(&seq);
        aseq_append_segment(&seq, "A", 1, SEG_D, 5, -1);

        AlleleCallResult r;
        allele_call_derive(&bm, &seq, SEG_D, &r, NULL);

        if (!allele_call_is_short_d(&r, &bm)) {
            fprintf(stderr, "1bp D should be Short-D (threshold=%d)\n",
                    bm.short_d_threshold);
            allele_bitmap_destroy(&bm);
            allele_pool_destroy(&pool);
            return 1;
        }
    }

    /* Full D → should NOT be short */
    {
        ASeq seq;
        aseq_init(&seq);
        aseq_append_segment(&seq, "AATTGGCCAATT", 12, SEG_D, 0, -1);

        AlleleCallResult r;
        allele_call_derive(&bm, &seq, SEG_D, &r, NULL);

        if (allele_call_is_short_d(&r, &bm)) {
            fprintf(stderr, "Full D should NOT be Short-D\n");
            allele_bitmap_destroy(&bm);
            allele_pool_destroy(&pool);
            return 1;
        }
    }

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 0;
}

/* ══════════════════════════════════════════════════════════════════
 * Mutation Tests (god-aligner: uses current base)
 * ══════════════════════════════════════════════════════════════════ */

/**
 * Test: Mutations at distinguishing position change the allele call.
 * GOD-ALIGNER: bitmap uses n->current, so the call reflects the
 * observable sequence.
 */
static int test_mutations_change_call(void) {
    AllelePool pool = allele_pool_create(4);
    Allele v1 = make_allele("V1*01", "ATCGATCG", 5, SEG_V);
    Allele v2 = make_allele("V2*01", "ATCGAACG", 5, SEG_V);  /* differs at pos 5 */
    allele_pool_add(&pool, &v1);
    allele_pool_add(&pool, &v2);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    /* Build V1 sequence then mutate position 5 (the distinguishing one) */
    ASeq seq;
    aseq_init(&seq);
    aseq_append_segment(&seq, "ATCGATCG", 8, SEG_V, 0, 5);

    /* Mutate pos 5 from T→A (now current matches V2 at the diff pos) */
    Nuc *n = seq.head;
    for (int i = 0; i < 5; i++) n = n->next;
    aseq_mutate(&seq, n, 'A', NUC_FLAG_MUTATED);

    AlleleCallResult r;
    allele_call_derive(&bm, &seq, SEG_V, &r, NULL);

    /* GOD-ALIGNER: V2 should win (8/8) since current at pos 5 is 'A'
     * which matches V2. V1 only scores 7/8. */
    if (r.count != 1 || r.indices[0] != 1) {
        fprintf(stderr, "Expected V2*01 after mutation at diff pos, "
                "got count=%d idx=%d\n", r.count,
                r.count > 0 ? r.indices[0] : -1);
        allele_bitmap_destroy(&bm);
        allele_pool_destroy(&pool);
        return 1;
    }

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 0;
}

/**
 * Test: Trimming away the distinguishing position creates a tie
 * (equivalent to "trimmed region excludes diff" from old reassess test).
 */
static int test_trimming_creates_equivalence(void) {
    AllelePool pool = allele_pool_create(4);
    Allele v1 = make_allele("V1*01", "ATCGATCG", 3, SEG_V);
    Allele v2 = make_allele("V2*01", "ATCGAACG", 3, SEG_V);  /* differs at pos 5 */
    allele_pool_add(&pool, &v1);
    allele_pool_add(&pool, &v2);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    /* Only include first 5 positions (0..4) → pos 5 is trimmed away */
    ASeq seq;
    aseq_init(&seq);
    aseq_append_segment(&seq, "ATCGA", 5, SEG_V, 0, 3);

    AlleleCallResult r;
    allele_call_derive(&bm, &seq, SEG_V, &r, NULL);

    /* Both alleles are identical in positions 0..4 → tie */
    if (r.count != 2) {
        fprintf(stderr, "Expected 2 tied (trimmed diff away), got %d\n", r.count);
        allele_bitmap_destroy(&bm);
        allele_pool_destroy(&pool);
        return 1;
    }

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 0;
}

/* ══════════════════════════════════════════════════════════════════
 * Integration: full pipeline → bitmap correction
 * ══════════════════════════════════════════════════════════════════ */

static int test_integration_pipeline_corrections(void) {
    srand(42);

    SimConfig cfg;
    sim_config_init(&cfg, CHAIN_IGH);

    cfg.v_alleles = allele_pool_create(4);
    cfg.d_alleles = allele_pool_create(4);
    cfg.j_alleles = allele_pool_create(4);

    Allele v1 = make_allele("IGHV1*01", "ATCGATCGATCGATCGATCGATCGATCGAT", 24, SEG_V);
    Allele v2 = make_allele("IGHV1*02", "ATCGATCGATCGATCGATCGATCGATCXXX", 24, SEG_V);
    allele_pool_add(&cfg.v_alleles, &v1);
    allele_pool_add(&cfg.v_alleles, &v2);

    Allele d1 = make_allele("IGHD1*01", "GTATTACTGTGC", 0, SEG_D);
    Allele d2 = make_allele("IGHD2*01", "GTATTACTAAAA", 0, SEG_D);
    allele_pool_add(&cfg.d_alleles, &d1);
    allele_pool_add(&cfg.d_alleles, &d2);

    Allele j1 = make_allele("IGHJ1*01", "ACTACTTTGACTACTGGGGC", 5, SEG_J);
    allele_pool_add(&cfg.j_alleles, &j1);

    AlleleCorrectionSet corr = allele_correction_set_build(&cfg);

    Pipeline p = pipeline_build(&cfg);
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    pipeline_execute(&p, &cfg, &seq, &rec, 0, NULL, NULL);

    /* Derive V call */
    AlleleCallResult vr;
    allele_call_derive(&corr.v_bitmap, &seq, SEG_V, &vr, NULL);
    char v_call[256];
    int v_idx = allele_pool_find_index(&cfg.v_alleles, rec.v_allele->name);
    allele_call_format(&vr, &cfg.v_alleles, v_idx, v_call, sizeof(v_call));

    printf("\n    V=%s → call=%s (score=%d/%d)",
           rec.v_allele->name, v_call, vr.best_score, vr.total_positions);

    /* Should find at least the true allele */
    if (vr.count < 1) {
        fprintf(stderr, "V call should have at least 1 match\n");
        allele_correction_set_destroy(&corr);
        sim_config_destroy(&cfg);
        return 1;
    }

    printf("\n    ");

    allele_correction_set_destroy(&corr);
    sim_config_destroy(&cfg);
    return 0;
}

/**
 * Test: allele_call_format puts true allele first.
 */
static int test_format_true_first(void) {
    AllelePool pool = allele_pool_create(4);
    Allele v1 = make_allele("V1*01", "AAAA", 2, SEG_V);
    Allele v2 = make_allele("V2*01", "AAAA", 2, SEG_V);
    allele_pool_add(&pool, &v1);
    allele_pool_add(&pool, &v2);

    AlleleBitmap bm = allele_bitmap_build(&pool);

    ASeq seq;
    aseq_init(&seq);
    aseq_append_segment(&seq, "AAAA", 4, SEG_V, 0, 2);

    AlleleCallResult r;
    allele_call_derive(&bm, &seq, SEG_V, &r, NULL);

    /* Format with true_idx=1 (V2*01): V2*01 should be listed first */
    char buf[256];
    allele_call_format(&r, &pool, 1, buf, sizeof(buf));

    if (strncmp(buf, "V2*01", 5) != 0) {
        fprintf(stderr, "Expected V2*01 first, got '%s'\n", buf);
        allele_bitmap_destroy(&bm);
        allele_pool_destroy(&pool);
        return 1;
    }
    if (strstr(buf, "V1*01") == NULL) {
        fprintf(stderr, "Expected V1*01 in call string: '%s'\n", buf);
        allele_bitmap_destroy(&bm);
        allele_pool_destroy(&pool);
        return 1;
    }

    allele_bitmap_destroy(&bm);
    allele_pool_destroy(&pool);
    return 0;
}

/* ── Main ─────────────────────────────────────────────────────── */

int main(void) {
    int failures = 0;
    int total = 0;

    printf("=== Bitmap Correction Tests ===\n");

    printf("\nV corrections:\n");
    TEST(test_v_correction_distinct);
    TEST(test_v_correction_identical);

    printf("\nJ corrections:\n");
    TEST(test_j_correction_distinct);

    printf("\nD corrections:\n");
    TEST(test_d_correction);
    TEST(test_d_correction_asymmetric);
    TEST(test_short_d_detection);

    printf("\nGod-aligner semantics:\n");
    TEST(test_mutations_change_call);
    TEST(test_trimming_creates_equivalence);

    printf("\nIntegration:\n");
    TEST(test_integration_pipeline_corrections);
    TEST(test_format_true_first);

    printf("\n%d/%d tests passed.\n", total - failures, total);
    return failures > 0 ? 1 : 0;
}
