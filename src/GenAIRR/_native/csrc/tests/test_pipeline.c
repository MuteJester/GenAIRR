/**
 * test_pipeline.c — Tests for pipeline build/execute with mock alleles.
 *
 * Creates a minimal SimConfig with fake V/D/J alleles to verify
 * that the pipeline builds correctly and produces valid ASeq output.
 */

#include "genairr/genairr.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define TEST(name) \
    do { \
        printf("  %-50s", #name); \
        int _r = name(); \
        printf("%s\n", _r == 0 ? "PASS" : "FAIL"); \
        if (_r) failures++; \
        total++; \
    } while (0)

/* ── Helper: create a mock allele ─────────────────────────────── */

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

/* ── Helper: create a minimal IGH config ──────────────────────── */

static SimConfig make_igh_config(void) {
    SimConfig cfg;
    sim_config_init(&cfg, CHAIN_IGH);

    /* A simple V allele: 30bp, anchor at position 24 (conserved Cys) */
    Allele v = make_allele(
        "IGHV1-2*01",
        "ATCGATCGATCGATCGATCGATCGATCGAT",  /* 30bp */
        24, SEG_V
    );

    /* A simple D allele: 12bp */
    Allele d = make_allele(
        "IGHD3-3*01",
        "GTATTACTGTGC",  /* 12bp */
        0, SEG_D
    );

    /* A simple J allele: 20bp, anchor at position 5 (conserved Trp) */
    Allele j = make_allele(
        "IGHJ4*01",
        "ACTACTTTGACTACTGGGGC",  /* 20bp */
        5, SEG_J
    );

    cfg.v_alleles = allele_pool_create(4);
    cfg.d_alleles = allele_pool_create(4);
    cfg.j_alleles = allele_pool_create(4);

    allele_pool_add(&cfg.v_alleles, &v);
    allele_pool_add(&cfg.d_alleles, &d);
    allele_pool_add(&cfg.j_alleles, &j);

    return cfg;
}

/* ── Test: pipeline build for IGH ─────────────────────────────── */

static int test_pipeline_build_igh(void) {
    SimConfig cfg = make_igh_config();
    Pipeline p = pipeline_build(&cfg);

    /* IGH (VDJ) should have: sample_v, sample_d, sample_j,
     * trim_v, trim_d, trim_j, assemble, assess = 8 steps */
    if (p.n_steps != 8) {
        fprintf(stderr, "Expected 8 steps, got %d\n", p.n_steps);
        sim_config_destroy(&cfg);
        return 1;
    }
    if (p.retry_boundary != -1) {
        sim_config_destroy(&cfg);
        return 1;  /* productive not enabled */
    }

    sim_config_destroy(&cfg);
    return 0;
}

/* ── Test: pipeline build for IGK (VJ, no D) ──────────────────── */

static int test_pipeline_build_igk(void) {
    SimConfig cfg;
    sim_config_init(&cfg, CHAIN_IGK);

    Allele v = make_allele("IGKV1-5*01", "ATCGATCGATCGATCGATCG", 15, SEG_V);
    Allele j = make_allele("IGKJ1*01", "TTTGGCCAAGGG", 3, SEG_J);

    cfg.v_alleles = allele_pool_create(4);
    cfg.j_alleles = allele_pool_create(4);
    allele_pool_add(&cfg.v_alleles, &v);
    allele_pool_add(&cfg.j_alleles, &j);

    Pipeline p = pipeline_build(&cfg);

    /* VJ chain: sample_v, sample_j, trim_v, trim_j, assemble, assess = 6 */
    if (p.n_steps != 6) {
        fprintf(stderr, "Expected 6 steps for VJ, got %d\n", p.n_steps);
        sim_config_destroy(&cfg);
        return 1;
    }

    sim_config_destroy(&cfg);
    return 0;
}

/* ── Test: pipeline execute produces a sequence ───────────────── */

static int test_pipeline_execute(void) {
    srand(42);  /* deterministic */

    SimConfig cfg = make_igh_config();
    Pipeline p = pipeline_build(&cfg);

    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;

    pipeline_execute(&p, &cfg, &seq, &rec, 0, NULL, NULL);

    /* Should have produced a non-empty sequence */
    if (seq.length == 0) {
        fprintf(stderr, "Empty sequence after execution\n");
        sim_config_destroy(&cfg);
        return 1;
    }

    /* Should have V, NP1, D, NP2, J segments */
    if (!aseq_has_segment(&seq, SEG_V)) { sim_config_destroy(&cfg); return 1; }
    if (!aseq_has_segment(&seq, SEG_J)) { sim_config_destroy(&cfg); return 1; }
    /* NP1 might be empty (length 0), so don't require it */

    /* Alleles should be set */
    if (rec.v_allele == NULL) { sim_config_destroy(&cfg); return 1; }
    if (rec.d_allele == NULL) { sim_config_destroy(&cfg); return 1; }
    if (rec.j_allele == NULL) { sim_config_destroy(&cfg); return 1; }

    /* Print for visual inspection */
    char buf[GENAIRR_MAX_SEQ_LEN];
    aseq_to_string(&seq, buf, sizeof(buf));
    printf("\n    seq[%d]: %.40s%s", seq.length, buf,
           seq.length > 40 ? "..." : "");
    printf("\n    V=%d NP1=%d D=%d NP2=%d J=%d",
           aseq_segment_length(&seq, SEG_V),
           aseq_segment_length(&seq, SEG_NP1),
           aseq_segment_length(&seq, SEG_D),
           aseq_segment_length(&seq, SEG_NP2),
           aseq_segment_length(&seq, SEG_J));
    printf("\n    productive=%d stop=%d frame=%d\n    ",
           rec.productive, rec.stop_codon, rec.vj_in_frame);

    sim_config_destroy(&cfg);
    return 0;
}

/* ── Test: productive mode retries ────────────────────────────── */

static int test_productive_mode(void) {
    srand(123);

    SimConfig cfg = make_igh_config();
    cfg.features.productive = true;
    cfg.max_productive_attempts = 25;

    Pipeline p = pipeline_build(&cfg);

    if (p.retry_boundary < 0) {
        fprintf(stderr, "retry_boundary not set\n");
        sim_config_destroy(&cfg);
        return 1;
    }

    /* Execute — should retry until productive (or exhaust attempts) */
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;

    pipeline_execute(&p, &cfg, &seq, &rec, 0, NULL, NULL);

    printf("\n    productive=%d after pipeline (retry_boundary=%d)\n    ",
           rec.productive, p.retry_boundary);

    sim_config_destroy(&cfg);
    return 0;  /* Pass regardless — just verify no crash */
}

/* ── Test: debug print doesn't crash ──────────────────────────── */

static int test_debug_print(void) {
    srand(99);

    SimConfig cfg = make_igh_config();
    Pipeline p = pipeline_build(&cfg);

    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    pipeline_execute(&p, &cfg, &seq, &rec, 0, NULL, NULL);

    /* Should not crash */
    aseq_debug_print(&seq);

    sim_config_destroy(&cfg);
    return 0;
}

/* ── Main ─────────────────────────────────────────────────────── */

int main(void) {
    int failures = 0;
    int total = 0;

    printf("=== Pipeline Tests ===\n");

    TEST(test_pipeline_build_igh);
    TEST(test_pipeline_build_igk);
    TEST(test_pipeline_execute);
    TEST(test_productive_mode);
    TEST(test_debug_print);

    printf("\n%d/%d tests passed.\n", total - failures, total);
    return failures > 0 ? 1 : 0;
}
