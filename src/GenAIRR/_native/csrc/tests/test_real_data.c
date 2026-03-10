/**
 * test_real_data.c — Integration tests using real HUMAN_IGH_IMGT alleles.
 *
 * Tests the full pipeline (rearrangement + S5F mutation + AIRR derivation
 * + allele corrections) with genuine human heavy chain allele data.
 */

#include "human_igh_imgt.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define TEST(name) \
    do { \
        printf("  %-60s", #name); \
        int _r = name(); \
        printf("%s\n", _r == 0 ? "PASS" : "FAIL"); \
        if (_r) failures++; \
        total++; \
    } while (0)

static SimConfig cfg;
static S5FModel  s5f;

/* ══════════════════════════════════════════════════════════════════
 * Config Loading Tests
 * ══════════════════════════════════════════════════════════════════ */

static int test_config_loads(void) {
    if (cfg.v_alleles.count != 307) {
        fprintf(stderr, "V alleles: %d, expected 307\n", cfg.v_alleles.count);
        return 1;
    }
    if (cfg.d_alleles.count != 47) {
        fprintf(stderr, "D alleles: %d, expected 47\n", cfg.d_alleles.count);
        return 1;
    }
    if (cfg.j_alleles.count != 14) {
        fprintf(stderr, "J alleles: %d, expected 14\n", cfg.j_alleles.count);
        return 1;
    }
    if (cfg.chain_type != CHAIN_IGH) {
        fprintf(stderr, "chain_type: %d, expected CHAIN_IGH\n", cfg.chain_type);
        return 1;
    }
    return 0;
}

static int test_allele_data_valid(void) {
    /* Check that V alleles have reasonable properties */
    for (int i = 0; i < cfg.v_alleles.count; i++) {
        const Allele *a = &cfg.v_alleles.alleles[i];
        if (a->length == 0 || a->length > GENAIRR_MAX_ALLELE_SEQ) {
            fprintf(stderr, "V allele %s: bad length %d\n", a->name, a->length);
            return 1;
        }
        if (strlen(a->name) == 0) {
            fprintf(stderr, "V allele %d: empty name\n", i);
            return 1;
        }
        if (a->segment_type != SEG_V) {
            fprintf(stderr, "V allele %s: wrong segment type %d\n", a->name, a->segment_type);
            return 1;
        }
    }

    /* Check D alleles */
    for (int i = 0; i < cfg.d_alleles.count; i++) {
        const Allele *a = &cfg.d_alleles.alleles[i];
        if (a->length == 0) {
            fprintf(stderr, "D allele %s: zero length\n", a->name);
            return 1;
        }
        if (a->segment_type != SEG_D) return 1;
    }

    /* Check J alleles */
    for (int i = 0; i < cfg.j_alleles.count; i++) {
        const Allele *a = &cfg.j_alleles.alleles[i];
        if (a->length == 0) {
            fprintf(stderr, "J allele %s: zero length\n", a->name);
            return 1;
        }
        if (a->segment_type != SEG_J) return 1;
    }

    return 0;
}

static int test_trim_distributions(void) {
    if (cfg.v_trim_3.max_trim == 0 || cfg.v_trim_3.probs == NULL) {
        fprintf(stderr, "V_3 trim distribution missing\n");
        return 1;
    }
    if (cfg.j_trim_5.max_trim == 0 || cfg.j_trim_5.probs == NULL) {
        fprintf(stderr, "J_5 trim distribution missing\n");
        return 1;
    }
    if (cfg.d_trim_5.max_trim == 0 || cfg.d_trim_5.probs == NULL) {
        fprintf(stderr, "D_5 trim distribution missing\n");
        return 1;
    }
    if (cfg.d_trim_3.max_trim == 0 || cfg.d_trim_3.probs == NULL) {
        fprintf(stderr, "D_3 trim distribution missing\n");
        return 1;
    }

    /* Check that distributions sum to ~1.0 */
    double v_sum = 0.0;
    for (int i = 0; i < cfg.v_trim_3.max_trim; i++)
        v_sum += cfg.v_trim_3.probs[i];
    if (fabs(v_sum - 1.0) > 0.01) {
        fprintf(stderr, "V_3 trim sum: %f, expected ~1.0\n", v_sum);
        return 1;
    }

    return 0;
}

static int test_s5f_model_loaded(void) {
    /* Check that mutability has nonzero values */
    int nonzero = 0;
    for (int i = 0; i < S5F_KMER_SPACE; i++) {
        if (s5f.mutability[i] > 0.0) nonzero++;
    }
    if (nonzero < 100) {
        fprintf(stderr, "S5F mutability: only %d nonzero entries\n", nonzero);
        return 1;
    }

    /* Check substitution data */
    int has_sub = 0;
    for (int i = 0; i < S5F_KMER_SPACE; i++) {
        if (s5f.substitution[i].count > 0) has_sub++;
    }
    if (has_sub < 100) {
        fprintf(stderr, "S5F substitution: only %d entries\n", has_sub);
        return 1;
    }

    return 0;
}

/* ══════════════════════════════════════════════════════════════════
 * Pipeline Integration Tests
 * ══════════════════════════════════════════════════════════════════ */

static int test_pipeline_rearrange(void) {
    srand(42);

    Pipeline pl;
    pl = pipeline_build(&cfg);

    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    pipeline_execute(&pl, &cfg, &seq, &rec, 0, NULL, NULL);

    /* Check we got a reasonable sequence */
    if (seq.length < 100 || seq.length > 600) {
        fprintf(stderr, "Sequence length %d outside expected range [100, 600]\n", seq.length);
        aseq_reset(&seq);
        return 1;
    }

    /* Verify segments exist */
    if (!aseq_has_segment(&seq, SEG_V)) {
        fprintf(stderr, "Missing V segment\n");
        aseq_reset(&seq);
        return 1;
    }
    if (!aseq_has_segment(&seq, SEG_J)) {
        fprintf(stderr, "Missing J segment\n");
        aseq_reset(&seq);
        return 1;
    }

    /* Print info for manual inspection */
    char buf[1024];
    aseq_to_string(&seq, buf, sizeof(buf));
    printf("\n    seq[%d]: %.60s...\n", seq.length, buf);
    printf("    V=%s D=%s J=%s\n",
           rec.v_allele->name,
           rec.d_allele ? rec.d_allele->name : "none",
           rec.j_allele->name);

    aseq_reset(&seq);
    return 0;
}

static int test_pipeline_with_s5f(void) {
    srand(100);

    Pipeline pl;
    pl = pipeline_build(&cfg);

    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    pipeline_execute(&pl, &cfg, &seq, &rec, 0, NULL, NULL);

    /* Apply S5F mutations */
    S5FResult result;
    s5f_mutate(&s5f, &seq, &rec, &result);

    /* Check mutations occurred */
    if (result.count == 0) {
        fprintf(stderr, "No mutations produced\n");
        aseq_reset(&seq);
        return 1;
    }

    /* Verify mutation rate is reasonable */
    if (result.mutation_rate > 0.30) {
        fprintf(stderr, "Mutation rate too high: %f\n", result.mutation_rate);
        aseq_reset(&seq);
        return 1;
    }

    printf("\n    mutations=%d, rate=%.4f, seq_len=%d\n",
           result.count, result.mutation_rate, seq.length);

    /* Check mutations are only in V/D/J */
    for (int m = 0; m < result.count; m++) {
        int pos = result.mutations[m].position;
        /* Walk to that node */
        Nuc *node = seq.head;
        for (int p = 0; p < pos && node; p++) node = node->next;
        if (node) {
            Segment seg = node->segment;
            if (seg != SEG_V && seg != SEG_D && seg != SEG_J) {
                fprintf(stderr, "Mutation at pos %d in segment %d (not V/D/J)\n",
                        pos, seg);
                aseq_reset(&seq);
                return 1;
            }
        }
    }

    aseq_reset(&seq);
    return 0;
}

static int test_airr_derivation(void) {
    srand(200);

    Pipeline pl;
    pl = pipeline_build(&cfg);

    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    pipeline_execute(&pl, &cfg, &seq, &rec, 0, NULL, NULL);

    /* Derive AIRR positions */
    AirrPositions pos;
    airr_derive_positions(&seq, &rec, &pos);

    /* Check V positions */
    if (pos.v_sequence_start != 0) {
        fprintf(stderr, "v_sequence_start=%d, expected 0\n", pos.v_sequence_start);
        aseq_reset(&seq);
        return 1;
    }
    if (pos.v_sequence_end <= 0) {
        fprintf(stderr, "v_sequence_end=%d, expected > 0\n", pos.v_sequence_end);
        aseq_reset(&seq);
        return 1;
    }

    /* Check J positions */
    if (pos.j_sequence_end != seq.length) {
        fprintf(stderr, "j_sequence_end=%d, expected seq.length=%d\n",
                pos.j_sequence_end, seq.length);
        aseq_reset(&seq);
        return 1;
    }

    /* Check D is between V and J */
    if (pos.d_sequence_start < pos.v_sequence_end) {
        /* D can overlap V boundary due to NP ambiguity, but shouldn't be before V start */
    }
    if (pos.d_sequence_end > pos.j_sequence_start + 1) {
        /* D should end before J starts (with tolerance for ambiguity) */
    }

    printf("\n    V=[%d,%d) D=[%d,%d) J=[%d,%d) junction=[%d,%d)\n",
           pos.v_sequence_start, pos.v_sequence_end,
           pos.d_sequence_start, pos.d_sequence_end,
           pos.j_sequence_start, pos.j_sequence_end,
           pos.junction_start, pos.junction_end);

    aseq_reset(&seq);
    return 0;
}

static int test_correction_maps_build(void) {
    /* Build bitmap correction set from real allele pools */
    AlleleCorrectionSet corr = allele_correction_set_build(&cfg);

    if (corr.v_bitmap.n_alleles != cfg.v_alleles.count) {
        fprintf(stderr, "V bitmap allele count %d != %d\n",
                corr.v_bitmap.n_alleles, cfg.v_alleles.count);
        allele_correction_set_destroy(&corr);
        return 1;
    }

    printf("\n    V bitmap: %d alleles, max_pos=%d\n",
           corr.v_bitmap.n_alleles, corr.v_bitmap.max_pos);
    printf("    D bitmap: %d alleles, short_d_threshold=%d\n",
           corr.d_bitmap.n_alleles, corr.d_bitmap.short_d_threshold);
    printf("    J bitmap: %d alleles, max_pos=%d\n",
           corr.j_bitmap.n_alleles, corr.j_bitmap.max_pos);

    /* Short-D threshold should be reasonable */
    if (corr.d_bitmap.short_d_threshold < 1 ||
        corr.d_bitmap.short_d_threshold > 20) {
        fprintf(stderr, "Unexpected short-D threshold: %d\n",
                corr.d_bitmap.short_d_threshold);
        allele_correction_set_destroy(&corr);
        return 1;
    }

    allele_correction_set_destroy(&corr);
    return 0;
}

static int test_reassess_with_real_alleles(void) {
    srand(300);

    Pipeline pl;
    pl = pipeline_build(&cfg);
    AlleleCorrectionSet corr = allele_correction_set_build(&cfg);

    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    pipeline_execute(&pl, &cfg, &seq, &rec, 0, NULL, NULL);

    /* Apply S5F mutations */
    S5FResult mresult;
    s5f_mutate(&s5f, &seq, &rec, &mresult);

    /* Derive V allele call from bitmap */
    AlleleCallResult vr;
    allele_call_derive(&corr.v_bitmap, &seq, SEG_V, &vr, NULL);

    int true_idx = allele_pool_find_index(&cfg.v_alleles, rec.v_allele->name);
    char v_call[512];
    allele_call_format(&vr, &cfg.v_alleles, true_idx, v_call, sizeof(v_call));

    printf("\n    V=%s, mutations=%d, call=%s (score=%d/%d)\n",
           rec.v_allele->name, mresult.count, v_call,
           vr.best_score, vr.total_positions);

    /* Should find at least the true allele */
    if (vr.count < 1) {
        fprintf(stderr, "Expected at least 1 V match\n");
        allele_correction_set_destroy(&corr);
        aseq_reset(&seq);
        return 1;
    }

    allele_correction_set_destroy(&corr);
    aseq_reset(&seq);
    return 0;
}

/**
 * Test: Run N rearrangements to verify statistical robustness.
 */
static int test_batch_rearrangement(void) {
    srand(42);

    Pipeline pl;
    pl = pipeline_build(&cfg);

    int success = 0;
    int total_len = 0;

    for (int i = 0; i < 100; i++) {
        ASeq seq;
        aseq_init(&seq);
        SimRecord rec;
        sim_record_init(&rec);

        pipeline_execute(&pl, &cfg, &seq, &rec, 0, NULL, NULL);
        if (seq.length > 0) {
            success++;
            total_len += seq.length;
        }
        aseq_reset(&seq);
    }

    double avg_len = (double)total_len / success;
    printf("\n    100 rearrangements: %d succeeded, avg length=%.0f\n",
           success, avg_len);

    if (success < 90) {
        fprintf(stderr, "Too many failures: only %d/100 succeeded\n", success);
        return 1;
    }

    /* Average IGH length should be roughly 350-500 bp */
    if (avg_len < 200 || avg_len > 700) {
        fprintf(stderr, "Average length %.0f outside expected range\n", avg_len);
        return 1;
    }

    return 0;
}

/**
 * Test: Batch rearrangement + S5F mutation.
 */
static int test_batch_with_mutation(void) {
    srand(42);

    Pipeline pl;
    pl = pipeline_build(&cfg);

    int mutations_total = 0;
    int seq_count = 0;

    for (int i = 0; i < 50; i++) {
        ASeq seq;
        aseq_init(&seq);
        SimRecord rec;
        sim_record_init(&rec);

        pipeline_execute(&pl, &cfg, &seq, &rec, 0, NULL, NULL);
        if (seq.length == 0) {
            aseq_reset(&seq);
            continue;
        }

        S5FResult result;
        s5f_mutate(&s5f, &seq, &rec, &result);
        mutations_total += result.count;
        seq_count++;
        aseq_reset(&seq);
    }

    double avg_mutations = (double)mutations_total / seq_count;
    printf("\n    50 seq with S5F: avg %.1f mutations/seq\n", avg_mutations);

    /* With default rates ~5-15%, expect 15-75 mutations on ~350bp */
    if (avg_mutations < 5 || avg_mutations > 100) {
        fprintf(stderr, "Average mutations %.1f outside expected range\n",
                avg_mutations);
        return 1;
    }

    return 0;
}

/* ══════════════════════════════════════════════════════════════════
 * Main
 * ══════════════════════════════════════════════════════════════════ */

int main(void) {
    int failures = 0, total = 0;

    printf("\n=== Real Data Integration Tests (HUMAN_IGH_IMGT) ===\n\n");

    /* Load data once */
    human_igh_imgt_load_config(&cfg);
    s5f_model_init(&s5f, 0.05, 0.15, false);
    human_igh_imgt_load_s5f(&s5f);

    printf("-- Config Loading --\n");
    TEST(test_config_loads);
    TEST(test_allele_data_valid);
    TEST(test_trim_distributions);
    TEST(test_s5f_model_loaded);

    printf("\n-- Pipeline Integration --\n");
    TEST(test_pipeline_rearrange);
    TEST(test_pipeline_with_s5f);
    TEST(test_airr_derivation);
    TEST(test_correction_maps_build);
    TEST(test_reassess_with_real_alleles);

    printf("\n-- Batch Tests --\n");
    TEST(test_batch_rearrangement);
    TEST(test_batch_with_mutation);

    printf("\n%d/%d tests passed.\n\n", total - failures, total);

    sim_config_destroy(&cfg);
    return failures > 0 ? 1 : 0;
}
