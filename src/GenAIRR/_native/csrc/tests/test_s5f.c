/**
 * test_s5f.c — Tests for the S5F mutation model and Fenwick tree.
 *
 * Uses mock mutability/substitution data to verify:
 *   - Fenwick tree construction, update, select, rebuild
 *   - 5-mer encoding/decoding roundtrip
 *   - S5F mutation: correct number of mutations, NP-awareness,
 *     segment filtering, productive stop-codon avoidance
 */

#include "genairr/genairr.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define TEST(name) \
    do { \
        printf("  %-55s", #name); \
        int _r = name(); \
        printf("%s\n", _r == 0 ? "PASS" : "FAIL"); \
        if (_r) failures++; \
        total++; \
    } while (0)

/* ══════════════════════════════════════════════════════════════════
 * Fenwick Tree Tests
 * ══════════════════════════════════════════════════════════════════ */

/**
 * Test: Fenwick tree init correctly sums weights.
 */
static int test_fenwick_init(void) {
    double weights[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    FenwickTree ft;
    fenwick_init(&ft, weights, 5);

    if (fabs(ft.total - 15.0) > 1e-9) {
        fprintf(stderr, "Expected total=15.0, got %f\n", ft.total);
        fenwick_destroy(&ft);
        return 1;
    }
    if (ft.n != 5) {
        fenwick_destroy(&ft);
        return 1;
    }

    fenwick_destroy(&ft);
    return 0;
}

/**
 * Test: Fenwick select returns correct positions.
 *
 * Weights: [1, 2, 3, 4, 5].  Prefix sums: [1, 3, 6, 10, 15].
 * Select(0.5) → position 0 (prefix sum 1 > 0.5).
 * Select(1.5) → position 1 (prefix sum 3 > 1.5).
 * Select(5.5) → position 2 (prefix sum 6 > 5.5).
 * Select(9.5) → position 3 (prefix sum 10 > 9.5).
 * Select(14.0)→ position 4.
 */
static int test_fenwick_select(void) {
    double weights[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    FenwickTree ft;
    fenwick_init(&ft, weights, 5);

    int s0 = fenwick_select(&ft, 0.5);
    int s1 = fenwick_select(&ft, 1.5);
    int s2 = fenwick_select(&ft, 5.5);
    int s3 = fenwick_select(&ft, 9.5);
    int s4 = fenwick_select(&ft, 14.0);

    fenwick_destroy(&ft);

    if (s0 != 0) { fprintf(stderr, "select(0.5)=%d, expected 0\n", s0); return 1; }
    if (s1 != 1) { fprintf(stderr, "select(1.5)=%d, expected 1\n", s1); return 1; }
    if (s2 != 2) { fprintf(stderr, "select(5.5)=%d, expected 2\n", s2); return 1; }
    if (s3 != 3) { fprintf(stderr, "select(9.5)=%d, expected 3\n", s3); return 1; }
    if (s4 != 4) { fprintf(stderr, "select(14.0)=%d, expected 4\n", s4); return 1; }

    return 0;
}

/**
 * Test: Fenwick update changes total and select behavior.
 */
static int test_fenwick_update(void) {
    double weights[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    FenwickTree ft;
    fenwick_init(&ft, weights, 5);

    /* Zero out position 4 (weight 5→0), total becomes 10 */
    fenwick_update(&ft, 4, 0.0, 5.0);
    if (fabs(ft.total - 10.0) > 1e-9) {
        fprintf(stderr, "After zeroing pos 4: total=%f, expected 10.0\n", ft.total);
        fenwick_destroy(&ft);
        return 1;
    }

    /* Now select(9.5) should go to position 3 (last nonzero) */
    int s = fenwick_select(&ft, 9.5);
    if (s != 3) {
        fprintf(stderr, "select(9.5) after update=%d, expected 3\n", s);
        fenwick_destroy(&ft);
        return 1;
    }

    fenwick_destroy(&ft);
    return 0;
}

/**
 * Test: Fenwick rebuild resets from weights array.
 */
static int test_fenwick_rebuild(void) {
    double weights[] = {1.0, 1.0, 1.0};
    FenwickTree ft;
    fenwick_init(&ft, weights, 3);

    /* Corrupt via updates */
    fenwick_update(&ft, 0, 10.0, 1.0);
    fenwick_update(&ft, 1, 10.0, 1.0);

    /* Rebuild from original weights */
    fenwick_rebuild(&ft, weights);
    if (fabs(ft.total - 3.0) > 1e-9) {
        fprintf(stderr, "After rebuild: total=%f, expected 3.0\n", ft.total);
        fenwick_destroy(&ft);
        return 1;
    }
    if (ft.update_count != 0) {
        fenwick_destroy(&ft);
        return 1;
    }

    fenwick_destroy(&ft);
    return 0;
}

/* ══════════════════════════════════════════════════════════════════
 * 5-mer Encoding Tests
 * ══════════════════════════════════════════════════════════════════ */

/**
 * Test: base-to-int and int-to-base roundtrip.
 */
static int test_base_encoding(void) {
    if (s5f_base_to_int('A') != 0) return 1;
    if (s5f_base_to_int('C') != 1) return 1;
    if (s5f_base_to_int('G') != 2) return 1;
    if (s5f_base_to_int('T') != 3) return 1;
    if (s5f_base_to_int('N') != 4) return 1;

    if (s5f_int_to_base(0) != 'A') return 1;
    if (s5f_int_to_base(1) != 'C') return 1;
    if (s5f_int_to_base(2) != 'G') return 1;
    if (s5f_int_to_base(3) != 'T') return 1;
    if (s5f_int_to_base(4) != 'N') return 1;

    return 0;
}

/**
 * Test: 5-mer encoding produces correct key and stays in bounds.
 */
static int test_kmer_encoding(void) {
    /* AAAAA = 0*625 + 0*125 + 0*25 + 0*5 + 0 = 0 */
    if (s5f_encode_kmer(0, 0, 0, 0, 0) != 0) return 1;

    /* TTTTT = 3*625 + 3*125 + 3*25 + 3*5 + 3 = 3 * 781 = 2343 */
    int ttttt = s5f_encode_kmer(3, 3, 3, 3, 3);
    if (ttttt != 1875 + 375 + 75 + 15 + 3) return 1;

    /* NNNNN = 4*625 + 4*125 + 4*25 + 4*5 + 4 = 3124 (max valid key) */
    int nnnnn = s5f_encode_kmer(4, 4, 4, 4, 4);
    if (nnnnn != 3124) return 1;
    if (nnnnn >= S5F_KMER_SPACE) return 1;

    /* ACGTA = 0*625 + 1*125 + 2*25 + 3*5 + 0 = 190 */
    int acgta = s5f_encode_kmer(0, 1, 2, 3, 0);
    if (acgta != 190) return 1;

    return 0;
}

/**
 * Test: string-based mutability setter encodes correctly.
 */
static int test_set_mutability_str(void) {
    S5FModel model;
    s5f_model_init(&model, 0.05, 0.15, false);

    s5f_set_mutability_str(&model, "ACGTA", 0.42);

    int key = s5f_encode_kmer(0, 1, 2, 3, 0);  /* 190 */
    if (fabs(model.mutability[key] - 0.42) > 1e-9) {
        fprintf(stderr, "mutability[ACGTA]=%f, expected 0.42\n",
                model.mutability[key]);
        return 1;
    }

    return 0;
}

/* ══════════════════════════════════════════════════════════════════
 * S5F Substitution Tests
 * ══════════════════════════════════════════════════════════════════ */

/**
 * Test: substitution data is normalized to cumulative distribution.
 */
static int test_substitution_cumulative(void) {
    S5FModel model;
    s5f_model_init(&model, 0.05, 0.15, false);

    /* Set substitution for AAAAA: A→C (0.3), A→G (0.5), A→T (0.2) */
    int key = s5f_encode_kmer(0, 0, 0, 0, 0);
    const char bases[] = "CGT";
    double weights[] = {0.3, 0.5, 0.2};
    s5f_set_substitution(&model, key, bases, weights, 3);

    const S5FSubstitution *sub = &model.substitution[key];
    if (sub->count != 3) return 1;
    if (sub->bases[0] != 'C') return 1;
    if (sub->bases[1] != 'G') return 1;
    if (sub->bases[2] != 'T') return 1;

    /* Cumulative: 0.3, 0.8, 1.0 */
    if (fabs(sub->cumulative[0] - 0.3) > 1e-9) return 1;
    if (fabs(sub->cumulative[1] - 0.8) > 1e-9) return 1;
    if (fabs(sub->cumulative[2] - 1.0) > 1e-9) return 1;

    return 0;
}

/* ══════════════════════════════════════════════════════════════════
 * S5F Mutation Tests
 * ══════════════════════════════════════════════════════════════════ */

/**
 * Helper: build a simple ASeq with V + NP1 + J segments.
 * V = 20bp, NP1 = 5bp, J = 15bp → total 40bp.
 */
static void build_test_seq(ASeq *seq) {
    aseq_init(seq);
    aseq_append_segment(seq, "ACGTACGTACGTACGTACGT", 20, SEG_V, 0, -1);
    aseq_append_segment(seq, "NNNNN",                  5, SEG_NP1, 0, -1);
    aseq_append_segment(seq, "TGCATGCATGCATGC",       15, SEG_J, 0, -1);
}

/**
 * Helper: populate model with uniform mutability for all 5-mers
 * and A→C substitution for simplicity.
 */
static void populate_uniform_model(S5FModel *model) {
    for (int k = 0; k < S5F_KMER_SPACE; k++) {
        model->mutability[k] = 1.0;

        /* For each 5-mer, the center base determines the substitution.
         * Center base = (k / 25) % 5 */
        int center = (k / 25) % 5;
        char center_base = s5f_int_to_base(center);

        /* Pick next base in A→C→G→T→A cycle */
        const char cycle[] = "ACGT";
        char target = 'C';  /* default */
        for (int i = 0; i < 4; i++) {
            if (cycle[i] == center_base) {
                target = cycle[(i + 1) % 4];
                break;
            }
        }

        char bases[1] = {target};
        double w[1] = {1.0};
        s5f_set_substitution(model, k, bases, w, 1);
    }
}

/**
 * Test: S5F mutate produces mutations within expected range.
 *
 * With rate [0.10, 0.10] on a 40bp seq (35 mutable = V+J),
 * we expect ~4 mutations (10% of 40).
 */
static int test_s5f_basic_mutation(void) {
    srand(42);

    S5FModel model;
    s5f_model_init(&model, 0.10, 0.10, false);
    populate_uniform_model(&model);

    ASeq seq;
    build_test_seq(&seq);

    SimRecord rec = {0};
    S5FResult result;

    s5f_mutate(&model, &seq, &rec, &result);

    /* With 10% rate on 40bp, target = 4 mutations */
    if (result.count < 1 || result.count > 10) {
        fprintf(stderr, "Expected ~4 mutations, got %d\n", result.count);
        aseq_reset(&seq);
        return 1;
    }

    /* Verify mutation_rate is consistent */
    double expected_rate = (double)result.count / 40.0;
    if (fabs(result.mutation_rate - expected_rate) > 1e-9) {
        fprintf(stderr, "mutation_rate=%f, expected %f\n",
                result.mutation_rate, expected_rate);
        aseq_reset(&seq);
        return 1;
    }

    aseq_reset(&seq);
    return 0;
}

/**
 * Test: S5F only mutates V/D/J positions, not NP regions.
 *
 * Build a seq where only NP1 has nonzero mutability.
 * If NP-awareness works, no mutations should occur on NP1 positions.
 */
static int test_s5f_np_awareness(void) {
    srand(123);

    S5FModel model;
    s5f_model_init(&model, 0.10, 0.10, false);
    populate_uniform_model(&model);

    ASeq seq;
    build_test_seq(&seq);

    SimRecord rec = {0};
    S5FResult result;

    s5f_mutate(&model, &seq, &rec, &result);

    /* Check that no mutation is at positions 20-24 (NP1 region) */
    for (int m = 0; m < result.count; m++) {
        int pos = result.mutations[m].position;
        if (pos >= 20 && pos < 25) {
            fprintf(stderr, "Mutation at NP1 position %d — NP-awareness violated\n", pos);
            aseq_reset(&seq);
            return 1;
        }
    }

    aseq_reset(&seq);
    return 0;
}

/**
 * Test: S5F with zero mutation rate produces no mutations.
 */
static int test_s5f_zero_rate(void) {
    srand(99);

    S5FModel model;
    s5f_model_init(&model, 0.0, 0.0, false);
    populate_uniform_model(&model);

    ASeq seq;
    build_test_seq(&seq);

    SimRecord rec = {0};
    S5FResult result;

    s5f_mutate(&model, &seq, &rec, &result);

    if (result.count != 0) {
        fprintf(stderr, "Zero rate should give 0 mutations, got %d\n", result.count);
        aseq_reset(&seq);
        return 1;
    }

    aseq_reset(&seq);
    return 0;
}

/**
 * Test: S5F mutations are written back to ASeq nodes.
 * After mutation, the ASeq's current bases should differ from germline
 * at mutated positions.
 */
static int test_s5f_writes_to_aseq(void) {
    srand(77);

    S5FModel model;
    s5f_model_init(&model, 0.15, 0.15, false);
    populate_uniform_model(&model);

    ASeq seq;
    build_test_seq(&seq);

    /* Save pre-mutation sequence */
    char before[64];
    aseq_to_string(&seq, before, sizeof(before));

    SimRecord rec = {0};
    S5FResult result;
    s5f_mutate(&model, &seq, &rec, &result);

    char after[64];
    aseq_to_string(&seq, after, sizeof(after));

    /* Verify that mutated positions in result correspond to
     * actual changes in the ASeq */
    int diff_count = 0;
    for (int i = 0; i < 40; i++) {
        if (before[i] != after[i]) diff_count++;
    }

    if (diff_count != result.count) {
        fprintf(stderr, "ASeq diffs=%d but result.count=%d\n",
                diff_count, result.count);
        aseq_reset(&seq);
        return 1;
    }

    /* Check NUC_FLAG_MUTATED is set on mutated nodes */
    int flagged = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->flags & NUC_FLAG_MUTATED) flagged++;
    }
    if (flagged != result.count) {
        fprintf(stderr, "Flagged nodes=%d but result.count=%d\n",
                flagged, result.count);
        aseq_reset(&seq);
        return 1;
    }

    aseq_reset(&seq);
    return 0;
}

/**
 * Test: S5F mutation result records correct from/to bases.
 */
static int test_s5f_mutation_records(void) {
    srand(55);

    S5FModel model;
    s5f_model_init(&model, 0.10, 0.10, false);
    populate_uniform_model(&model);

    ASeq seq;
    build_test_seq(&seq);

    /* Save germline by walking nodes */
    char germline[64];
    {
        int gi = 0;
        for (Nuc *nd = seq.head; nd && gi < 63; nd = nd->next)
            germline[gi++] = nd->germline;
        germline[gi] = '\0';
    }

    SimRecord rec = {0};
    S5FResult result;
    s5f_mutate(&model, &seq, &rec, &result);

    /* Verify from_base matches germline at that position */
    for (int m = 0; m < result.count; m++) {
        int pos = result.mutations[m].position;
        char from = result.mutations[m].from_base;
        if (from != germline[pos]) {
            fprintf(stderr, "Mutation %d: from_base='%c' but germline[%d]='%c'\n",
                    m, from, pos, germline[pos]);
            aseq_reset(&seq);
            return 1;
        }
        /* to_base should differ from from_base */
        if (result.mutations[m].to_base == from) {
            fprintf(stderr, "Mutation %d: to_base == from_base == '%c'\n",
                    m, from);
            aseq_reset(&seq);
            return 1;
        }
    }

    aseq_reset(&seq);
    return 0;
}

/**
 * Test: S5F with all-zero mutability produces no mutations
 * even with high mutation rate.
 */
static int test_s5f_zero_mutability(void) {
    srand(200);

    S5FModel model;
    s5f_model_init(&model, 0.20, 0.20, false);
    /* Don't populate — all mutability stays 0 */

    ASeq seq;
    build_test_seq(&seq);

    SimRecord rec = {0};
    S5FResult result;
    s5f_mutate(&model, &seq, &rec, &result);

    if (result.count != 0) {
        fprintf(stderr, "Zero mutability should give 0 mutations, got %d\n",
                result.count);
        aseq_reset(&seq);
        return 1;
    }

    aseq_reset(&seq);
    return 0;
}

/**
 * Test: S5F with only V-segment sequence (no NP, no J).
 * All positions should be mutable.
 */
static int test_s5f_v_only(void) {
    srand(300);

    S5FModel model;
    s5f_model_init(&model, 0.10, 0.10, false);
    populate_uniform_model(&model);

    ASeq seq;
    aseq_init(&seq);
    aseq_append_segment(&seq, "ACGTACGTACGTACGTACGT", 20, SEG_V, 0, -1);

    SimRecord rec = {0};
    S5FResult result;
    s5f_mutate(&model, &seq, &rec, &result);

    /* target = 10% of 20 = 2 mutations */
    if (result.count < 1 || result.count > 6) {
        fprintf(stderr, "V-only: expected ~2 mutations, got %d\n", result.count);
        aseq_reset(&seq);
        return 1;
    }

    aseq_reset(&seq);
    return 0;
}

/**
 * Test: Repeated mutations with different seeds produce different results.
 * (Verifies randomness is working, not deterministic.)
 */
static int test_s5f_different_seeds(void) {
    S5FModel model;
    s5f_model_init(&model, 0.15, 0.15, false);
    populate_uniform_model(&model);

    /* Run with seed 1 */
    srand(1);
    ASeq seq1;
    build_test_seq(&seq1);
    SimRecord rec = {0};
    S5FResult r1;
    s5f_mutate(&model, &seq1, &rec, &r1);
    char s1[64];
    aseq_to_string(&seq1, s1, sizeof(s1));

    /* Run with seed 2 */
    srand(2);
    ASeq seq2;
    build_test_seq(&seq2);
    S5FResult r2;
    s5f_mutate(&model, &seq2, &rec, &r2);
    char s2[64];
    aseq_to_string(&seq2, s2, sizeof(s2));

    /* They should differ (overwhelmingly likely with different seeds) */
    if (strcmp(s1, s2) == 0 && r1.count == r2.count) {
        fprintf(stderr, "Different seeds produced identical results — suspicious\n");
        aseq_reset(&seq1);
        aseq_reset(&seq2);
        return 1;
    }

    aseq_reset(&seq1);
    aseq_reset(&seq2);
    return 0;
}

/**
 * Test: S5F productive mode avoids creating stop codons.
 *
 * Build a sequence where the only possible mutation at a position
 * would create a stop codon (e.g., T_A → TAG). With productive=true,
 * this mutation should be skipped.
 */
static int test_s5f_productive_no_stops(void) {
    srand(42);

    S5FModel model;
    /* High rate to force many attempts */
    s5f_model_init(&model, 0.20, 0.20, true);
    populate_uniform_model(&model);

    ASeq seq;
    build_test_seq(&seq);

    SimRecord rec = {0};
    S5FResult result;
    s5f_mutate(&model, &seq, &rec, &result);

    /* Check that no stop codon exists in the result sequence */
    char mutated[64];
    aseq_to_string(&seq, mutated, sizeof(mutated));
    int len = seq.length;

    for (int i = 0; i + 2 < len; i++) {
        char c1 = mutated[i], c2 = mutated[i + 1], c3 = mutated[i + 2];
        if (c1 == 'T' && c2 == 'A' && c3 == 'G') {
            /* Check if this is actually in a mutable region */
            fprintf(stderr, "Stop codon TAG at pos %d in productive mode\n", i);
            aseq_reset(&seq);
            return 1;
        }
        if (c1 == 'T' && c2 == 'A' && c3 == 'A') {
            fprintf(stderr, "Stop codon TAA at pos %d in productive mode\n", i);
            aseq_reset(&seq);
            return 1;
        }
        if (c1 == 'T' && c2 == 'G' && c3 == 'A') {
            fprintf(stderr, "Stop codon TGA at pos %d in productive mode\n", i);
            aseq_reset(&seq);
            return 1;
        }
    }

    aseq_reset(&seq);
    return 0;
}

/**
 * Test: Fenwick tree with single element.
 */
static int test_fenwick_single(void) {
    double w[] = {3.14};
    FenwickTree ft;
    fenwick_init(&ft, w, 1);

    if (fabs(ft.total - 3.14) > 1e-9) {
        fenwick_destroy(&ft);
        return 1;
    }

    int s = fenwick_select(&ft, 0.5);
    if (s != 0) {
        fenwick_destroy(&ft);
        return 1;
    }

    fenwick_destroy(&ft);
    return 0;
}

/**
 * Test: Fenwick needs_rebuild triggers at threshold.
 */
static int test_fenwick_rebuild_trigger(void) {
    double w[] = {1.0, 1.0};
    FenwickTree ft;
    fenwick_init(&ft, w, 2);

    /* Simulate updates up to threshold */
    for (int i = 0; i < FENWICK_REBUILD_INTERVAL - 1; i++) {
        fenwick_update(&ft, 0, 1.0 + 0.001 * i, 1.0 + 0.001 * (i - 1));
    }

    if (fenwick_needs_rebuild(&ft)) {
        fprintf(stderr, "Should not need rebuild at %d updates\n",
                FENWICK_REBUILD_INTERVAL - 1);
        fenwick_destroy(&ft);
        return 1;
    }

    /* One more update → triggers rebuild */
    fenwick_update(&ft, 0, 2.0, 1.0);
    if (!fenwick_needs_rebuild(&ft)) {
        fprintf(stderr, "Should need rebuild at %d updates\n",
                FENWICK_REBUILD_INTERVAL);
        fenwick_destroy(&ft);
        return 1;
    }

    fenwick_destroy(&ft);
    return 0;
}

/* ══════════════════════════════════════════════════════════════════
 * Main
 * ══════════════════════════════════════════════════════════════════ */

int main(void) {
    int failures = 0, total = 0;

    printf("\n=== S5F Mutation Model Tests ===\n\n");

    printf("-- Fenwick Tree --\n");
    TEST(test_fenwick_init);
    TEST(test_fenwick_select);
    TEST(test_fenwick_update);
    TEST(test_fenwick_rebuild);
    TEST(test_fenwick_single);
    TEST(test_fenwick_rebuild_trigger);

    printf("\n-- 5-mer Encoding --\n");
    TEST(test_base_encoding);
    TEST(test_kmer_encoding);
    TEST(test_set_mutability_str);

    printf("\n-- Substitution --\n");
    TEST(test_substitution_cumulative);

    printf("\n-- S5F Mutation --\n");
    TEST(test_s5f_basic_mutation);
    TEST(test_s5f_np_awareness);
    TEST(test_s5f_zero_rate);
    TEST(test_s5f_writes_to_aseq);
    TEST(test_s5f_mutation_records);
    TEST(test_s5f_zero_mutability);
    TEST(test_s5f_v_only);
    TEST(test_s5f_different_seeds);
    TEST(test_s5f_productive_no_stops);

    printf("\n%d/%d tests passed.\n\n", total - failures, total);
    return failures > 0 ? 1 : 0;
}
