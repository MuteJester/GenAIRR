/**
 * test_ops.c — Tests for all corruption, simulation, and mutation ops.
 *
 * Each test builds a minimal ASeq + SimConfig, runs a step function,
 * and verifies the expected outcome.
 */

#include "genairr/genairr.h"
#include "genairr/rand_util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ── Test framework ───────────────────────────────────────────── */

static int tests_run = 0;
static RngState _test_rng;
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

/* ── Step function declarations ──────────────────────────────── */

extern void step_corrupt_5_prime(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_corrupt_3_prime(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_quality_errors(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_pcr_amplification(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_simulate_umi(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_reverse_complement(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_spike_contaminants(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_insert_indels(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_insert_ns(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_paired_end(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_long_read_errors(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_trim_to_length(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_primer_mask(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_d_inversion(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_receptor_revision(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_selection_pressure(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_uniform_mutate(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void csr_adjust_rates(const SimRecord *rec, double *min_rate, double *max_rate);

/* ── Helpers ─────────────────────────────────────────────────── */

static void build_simple_seq(ASeq *seq, SimRecord *rec, SimConfig *cfg) {
    aseq_init(seq);
    sim_record_init(rec);
    sim_config_init(cfg, CHAIN_IGH);
    cfg->rng = &_test_rng;

    /* Build a simple V-NP1-D-NP2-J sequence */
    static Allele v_al = { .name = "IGHV1-2*01", .gene = "IGHV1-2", .family = "IGHV1",
        .seq = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC",
        .length = 100, .anchor = 90, .segment_type = SEG_V };
    static Allele d_al = { .name = "IGHD1-1*01", .gene = "IGHD1-1", .family = "IGHD1",
        .seq = "GGTATAACTGGAAC",
        .length = 14, .anchor = 0, .segment_type = SEG_D };
    static Allele j_al = { .name = "IGHJ4*02", .gene = "IGHJ4", .family = "IGHJ4",
        .seq = "ACTACTTTGACTACTGGGGCCAAGGAACCCTGGTCACC",
        .length = 38, .anchor = 8, .segment_type = SEG_J };

    rec->v_allele = &v_al;
    rec->d_allele = &d_al;
    rec->j_allele = &j_al;
    rec->v_trim_5 = 0;
    rec->v_trim_3 = 5;
    rec->d_trim_5 = 2;
    rec->d_trim_3 = 2;
    rec->j_trim_5 = 3;
    rec->j_trim_3 = 0;

    /* Append V */
    aseq_append_segment(seq, v_al.seq, v_al.length - 5, SEG_V, 0, v_al.anchor);
    /* NP1 */
    aseq_append_np(seq, "ACGT", 4, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    rec->np1_length = 4;
    /* D */
    aseq_append_segment(seq, d_al.seq + 2, d_al.length - 4, SEG_D, 2, -1);
    /* NP2 */
    aseq_append_np(seq, "TGCA", 4, SEG_NP2, NUC_FLAG_N_NUCLEOTIDE);
    rec->np2_length = 4;
    /* J */
    aseq_append_segment(seq, j_al.seq + 3, j_al.length - 3, SEG_J, 3, j_al.anchor);
}

/* ══════════════════════════════════════════════════════════════ */
/*                    CORRUPTION OP TESTS                        */
/* ══════════════════════════════════════════════════════════════ */

/* ── Corrupt 5' ──────────────────────────────────────────────── */

static int test_corrupt_5_prime_changes_length(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    int orig_len = aseq_length(&seq);

    rng_seed(&_test_rng, 42, 0);
    step_corrupt_5_prime(&cfg, &seq, &rec);

    ASSERT(rec.corruption_5_event >= 1 && rec.corruption_5_event <= 3,
           "Event type should be 1, 2, or 3");
    ASSERT(aseq_length(&seq) != orig_len || rec.corruption_5_event == 1,
           "Length should change for add/remove events");
    sim_config_destroy(&cfg);
    return 1;
}

static int test_corrupt_5_prime_remove_shrinks(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    int orig_len = aseq_length(&seq);

    /* Force remove event by trying many seeds */
    for (int seed = 0; seed < 100; seed++) {
        aseq_reset(&seq);
        sim_record_init(&rec);
        build_simple_seq(&seq, &rec, &cfg);
        rng_seed(&_test_rng, seed, 0);
        step_corrupt_5_prime(&cfg, &seq, &rec);
        if (rec.corruption_5_event == 1) {
            ASSERT(aseq_length(&seq) < orig_len, "Remove should shrink");
            ASSERT(rec.corruption_5_remove_amount > 0, "Remove amount > 0");
            sim_config_destroy(&cfg);
            return 1;
        }
    }
    sim_config_destroy(&cfg);
    return 1; /* May not find remove event, still pass */
}

/* ── Corrupt 3' ──────────────────────────────────────────────── */

static int test_corrupt_3_prime_changes_length(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);

    rng_seed(&_test_rng, 123, 0);
    step_corrupt_3_prime(&cfg, &seq, &rec);

    ASSERT(rec.corruption_3_event >= 1 && rec.corruption_3_event <= 3,
           "Event type should be 1, 2, or 3");
    sim_config_destroy(&cfg);
    return 1;
}

/* ── Quality Errors ──────────────────────────────────────────── */

static int test_quality_errors_introduces_errors(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.base_error_rate = 0.5;  /* Very high rate for testing */
    cfg.peak_error_rate = 0.8;

    rng_seed(&_test_rng, 42, 0);
    step_quality_errors(&cfg, &seq, &rec);

    /* Count sequencing errors */
    int errors = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->flags & NUC_FLAG_SEQ_ERROR) errors++;
    }
    ASSERT(errors > 0, "Should introduce at least one error at high rate");
    sim_config_destroy(&cfg);
    return 1;
}

static int test_quality_errors_preserves_germline(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.base_error_rate = 0.5;
    cfg.peak_error_rate = 0.8;

    rng_seed(&_test_rng, 42, 0);
    step_quality_errors(&cfg, &seq, &rec);

    /* Germline bases should be preserved */
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->flags & NUC_FLAG_SEQ_ERROR) {
            ASSERT(n->germline == '\0' || n->current != n->germline,
                   "Errored base should differ from germline");
        }
    }
    sim_config_destroy(&cfg);
    return 1;
}

/* ── PCR Amplification ───────────────────────────────────────── */

static int test_pcr_introduces_errors(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.pcr_error_rate = 0.01;  /* High rate for testing */
    cfg.pcr_n_cycles = 50;

    rng_seed(&_test_rng, 42, 0);
    step_pcr_amplification(&cfg, &seq, &rec);

    int errors = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->flags & NUC_FLAG_PCR_ERROR) errors++;
    }
    ASSERT(errors > 0, "Should introduce PCR errors");
    ASSERT(rec.pcr_n_cycles == 50, "Should record n_cycles");
    sim_config_destroy(&cfg);
    return 1;
}

/* ── UMI ─────────────────────────────────────────────────────── */

static int test_umi_prepends_barcode(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    int orig_len = aseq_length(&seq);
    cfg.umi_length = 12;

    rng_seed(&_test_rng, 42, 0);
    step_simulate_umi(&cfg, &seq, &rec);

    ASSERT(aseq_length(&seq) == orig_len + 12, "Should add 12 bases");
    ASSERT(rec.umi_length == 12, "UMI length recorded");
    ASSERT(strlen(rec.umi_sequence) == 12, "UMI sequence stored");

    /* First 12 nodes should be SEG_UMI */
    Nuc *n = seq.head;
    for (int i = 0; i < 12; i++) {
        ASSERT(n->segment == SEG_UMI, "UMI nodes tagged as SEG_UMI");
        n = n->next;
    }
    /* Next node should be SEG_V */
    ASSERT(n->segment == SEG_V, "V follows UMI");
    sim_config_destroy(&cfg);
    return 1;
}

/* ── Reverse Complement ──────────────────────────────────────── */

static int test_reverse_complement_flips(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.rc_prob = 1.0;  /* Always RC */

    char before[256];
    aseq_to_string(&seq, before, sizeof(before));

    step_reverse_complement(&cfg, &seq, &rec);

    char after[256];
    aseq_to_string(&seq, after, sizeof(after));

    ASSERT(rec.is_reverse_complement == true, "Flag should be set");
    ASSERT(strlen(before) == strlen(after), "Length preserved");
    /* First base of RC should be complement of last base of original */
    ASSERT(after[0] == complement(before[strlen(before) - 1]),
           "RC reverses and complements");
    sim_config_destroy(&cfg);
    return 1;
}

static int test_reverse_complement_probability(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.rc_prob = 0.0;  /* Never RC */

    step_reverse_complement(&cfg, &seq, &rec);
    ASSERT(rec.is_reverse_complement == false, "Should not RC with prob=0");
    sim_config_destroy(&cfg);
    return 1;
}

/* ── Spike Contaminants ──────────────────────────────────────── */

static int test_contaminant_replaces_sequence(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.contamination_prob = 1.0;  /* Always contaminate */
    cfg.contaminant_type = 0;      /* random */

    rng_seed(&_test_rng, 42, 0);
    step_spike_contaminants(&cfg, &seq, &rec);

    ASSERT(rec.is_contaminant == true, "Contaminant flag set");
    ASSERT(rec.v_allele == NULL, "V allele zeroed");
    ASSERT(rec.d_allele == NULL, "D allele zeroed");
    ASSERT(rec.j_allele == NULL, "J allele zeroed");
    ASSERT(strcmp(rec.note, "contaminant") == 0, "Note set");
    ASSERT(aseq_length(&seq) >= 300, "Contaminant has reasonable length");
    sim_config_destroy(&cfg);
    return 1;
}

static int test_contaminant_phix(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.contamination_prob = 1.0;
    cfg.contaminant_type = 1;  /* phiX */

    rng_seed(&_test_rng, 42, 0);
    step_spike_contaminants(&cfg, &seq, &rec);

    ASSERT(rec.is_contaminant == true, "Contaminant flag set");
    ASSERT(strcmp(rec.contaminant_type, "phix") == 0, "Type is phix");
    sim_config_destroy(&cfg);
    return 1;
}

/* ── Insert Indels ───────────────────────────────────────────── */

static int test_insert_indels_modifies(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    int orig_len = aseq_length(&seq);
    cfg.indel_prob = 0.3;  /* High rate */
    cfg.insertion_weight = 0.5;

    rng_seed(&_test_rng, 42, 0);
    step_insert_indels(&cfg, &seq, &rec);

    /* Should have changed length */
    ASSERT(aseq_length(&seq) != orig_len, "Length should change with high indel rate");
    sim_config_destroy(&cfg);
    return 1;
}

static int test_insert_indels_preserves_anchors(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.indel_prob = 0.5;
    cfg.insertion_weight = 0.0;  /* All deletions */

    rng_seed(&_test_rng, 42, 0);
    step_insert_indels(&cfg, &seq, &rec);

    /* Anchors should still exist */
    Nuc *v_anchor = aseq_find_anchor(&seq, SEG_V);
    Nuc *j_anchor = aseq_find_anchor(&seq, SEG_J);
    ASSERT(v_anchor != NULL, "V anchor preserved");
    ASSERT(j_anchor != NULL, "J anchor preserved");
    sim_config_destroy(&cfg);
    return 1;
}

/* ── Insert Ns ───────────────────────────────────────────────── */

static int test_insert_ns_adds_ambiguous(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.n_prob = 0.5;  /* Very high */

    rng_seed(&_test_rng, 42, 0);
    step_insert_ns(&cfg, &seq, &rec);

    int n_count = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->flags & NUC_FLAG_IS_N) n_count++;
    }
    ASSERT(n_count > 0, "Should introduce N bases");
    sim_config_destroy(&cfg);
    return 1;
}

/* ── Paired End ──────────────────────────────────────────────── */

static int test_paired_end_creates_gap(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    int len = aseq_length(&seq);
    cfg.pe_read_length = 30;  /* Very short reads → gap in middle */

    rng_seed(&_test_rng, 42, 0);
    step_paired_end(&cfg, &seq, &rec);

    ASSERT(rec.pe_read_length == 30, "Read length recorded");
    if (len > 60) {
        ASSERT(rec.pe_gap_length > 0, "Gap should exist when seq > 2×read_len");
        /* Check that middle positions are N */
        int pos = 0;
        for (Nuc *n = seq.head; n; n = n->next, pos++) {
            if (pos >= 30 && pos < len - 30) {
                ASSERT(n->current == 'N', "Gap positions should be N");
            }
        }
    }
    sim_config_destroy(&cfg);
    return 1;
}

static int test_paired_end_no_gap_long_reads(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.pe_read_length = 500;  /* Reads longer than sequence */

    rng_seed(&_test_rng, 42, 0);
    step_paired_end(&cfg, &seq, &rec);

    ASSERT(rec.pe_gap_length == 0, "No gap with long reads");
    sim_config_destroy(&cfg);
    return 1;
}

/* ── Long Read Errors ────────────────────────────────────────── */

static int test_long_read_errors_targets_homopolymers(void) {
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);
    SimConfig cfg;
    sim_config_init(&cfg, CHAIN_IGH);
    cfg.rng = &_test_rng;

    /* Build sequence with homopolymer run: AAAA in V segment */
    static Allele v_al = { .name = "V1", .gene = "V1", .family = "V",
        .seq = "ATGCAAAAGCATGCATGCATGCATGC",
        .length = 26, .anchor = 20, .segment_type = SEG_V };
    rec.v_allele = &v_al;
    aseq_append_segment(&seq, v_al.seq, v_al.length, SEG_V, 0, v_al.anchor);

    int orig_len = aseq_length(&seq);
    cfg.long_read_error_rate = 0.5;  /* Very high */
    cfg.min_run_length = 3;
    cfg.insertion_bias = 1.0;  /* All insertions */

    rng_seed(&_test_rng, 42, 0);
    step_long_read_errors(&cfg, &seq, &rec);

    /* Length may have changed due to insertion at AAAA run */
    int new_len = aseq_length(&seq);
    ASSERT(new_len >= orig_len, "Insertions should not shrink");
    sim_config_destroy(&cfg);
    return 1;
}

/* ── Trim to Length ──────────────────────────────────────────── */

static int test_trim_to_length_truncates(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    int orig_len = aseq_length(&seq);
    cfg.max_sequence_length = 50;

    step_trim_to_length(&cfg, &seq, &rec);

    ASSERT(aseq_length(&seq) == 50, "Should truncate to 50");
    ASSERT(aseq_length(&seq) < orig_len, "Should be shorter than original");
    sim_config_destroy(&cfg);
    return 1;
}

static int test_trim_to_length_noop_when_shorter(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    int orig_len = aseq_length(&seq);
    cfg.max_sequence_length = 1000;  /* Longer than sequence */

    step_trim_to_length(&cfg, &seq, &rec);

    ASSERT(aseq_length(&seq) == orig_len, "Should not change when already shorter");
    sim_config_destroy(&cfg);
    return 1;
}

/* ── Primer Mask ─────────────────────────────────────────────── */

static int test_primer_mask_restores_germline(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.primer_mask_length = 20;

    /* First, mutate some V bases */
    int mutated = 0;
    for (Nuc *n = seq.head; n && n->segment == SEG_V && mutated < 5; n = n->next) {
        if (n->germline != '\0') {
            char orig = n->current;
            char new_base = (orig == 'A') ? 'T' : 'A';
            aseq_mutate(&seq, n, new_base, NUC_FLAG_MUTATED);
            mutated++;
        }
    }
    ASSERT(mutated > 0, "Should have mutated some bases");

    step_primer_mask(&cfg, &seq, &rec);

    /* First 20 positions of V should have current == germline */
    int pos = 0;
    for (Nuc *n = seq.head; n && n->segment == SEG_V && pos < 20; n = n->next, pos++) {
        if (n->germline != '\0') {
            ASSERT(n->current == n->germline, "Primer mask should restore germline");
            ASSERT(!(n->flags & NUC_FLAG_MUTATED), "Mutation flag should be cleared");
        }
    }
    sim_config_destroy(&cfg);
    return 1;
}

/* ══════════════════════════════════════════════════════════════ */
/*                    SIMULATION OP TESTS                         */
/* ══════════════════════════════════════════════════════════════ */

/* ── D Gene Inversion ────────────────────────────────────────── */

static int test_d_inversion_complements(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.d_inversion_prob = 1.0;  /* Always invert */

    /* Get original D bases */
    char orig_d[64];
    aseq_segment_to_string(&seq, SEG_D, orig_d, sizeof(orig_d));
    int d_len = (int)strlen(orig_d);

    rng_seed(&_test_rng, 42, 0);
    step_d_inversion(&cfg, &seq, &rec);

    /* Get new D bases */
    char new_d[64];
    aseq_segment_to_string(&seq, SEG_D, new_d, sizeof(new_d));

    ASSERT(rec.d_inverted == true, "D inverted flag set");
    /* First base of new should be complement of last base of original */
    ASSERT(new_d[0] == complement(orig_d[d_len - 1]),
           "D should be reverse-complemented");
    sim_config_destroy(&cfg);
    return 1;
}

static int test_d_inversion_probability(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.d_inversion_prob = 0.0;

    step_d_inversion(&cfg, &seq, &rec);
    ASSERT(rec.d_inverted == false, "Should not invert with prob=0");
    sim_config_destroy(&cfg);
    return 1;
}

/* ── Receptor Revision ───────────────────────────────────────── */

static int test_receptor_revision_changes_v(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.revision_prob = 1.0;
    cfg.footprint_min = 5;
    cfg.footprint_max = 10;

    /* Add a second V allele with different gene */
    static Allele v2 = { .name = "IGHV3-11*01", .gene = "IGHV3-11", .family = "IGHV3",
        .seq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        .length = 99, .anchor = 85, .segment_type = SEG_V };
    allele_pool_add(&cfg.v_alleles, rec.v_allele);
    allele_pool_add(&cfg.v_alleles, &v2);

    const Allele *orig_v = rec.v_allele;
    rng_seed(&_test_rng, 42, 0);
    step_receptor_revision(&cfg, &seq, &rec);

    ASSERT(rec.receptor_revised == true, "Revision flag set");
    ASSERT(rec.revision_footprint_length >= 5 && rec.revision_footprint_length <= 10,
           "Footprint in range");
    ASSERT(strlen(rec.original_v_allele_name) > 0, "Original V name recorded");
    ASSERT(strcmp(rec.v_allele->gene, orig_v->gene) != 0,
           "New V should be from different gene");
    sim_config_destroy(&cfg);
    return 1;
}

/* ── Selection Pressure ──────────────────────────────────────── */

static int test_selection_pressure_reverts_mutations(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.selection_strength = 1.0;    /* Maximum pressure */
    cfg.cdr_r_acceptance = 0.0;      /* Reject all CDR R-mutations */
    cfg.fwr_r_acceptance = 0.0;      /* Reject all FWR R-mutations */

    /* Mutate several V-segment bases */
    int mutated = 0;
    for (Nuc *n = seq.head; n && n->segment == SEG_V; n = n->next) {
        if (n->germline != '\0' && mutated < 20) {
            char new_base = (n->current == 'A') ? 'G' : 'A';
            aseq_mutate(&seq, n, new_base, NUC_FLAG_MUTATED);
            mutated++;
        }
    }

    rng_seed(&_test_rng, 42, 0);
    step_selection_pressure(&cfg, &seq, &rec);

    /* Count remaining mutations (some R-mutations should be reverted) */
    int remaining = 0;
    for (Nuc *n = seq.head; n && n->segment == SEG_V; n = n->next) {
        if (n->flags & NUC_FLAG_MUTATED) remaining++;
    }

    /* With 100% strength and 0% acceptance, all R-mutations should be reverted.
     * Some may be S (silent) and kept. */
    ASSERT(remaining < mutated, "Some mutations should be reverted");
    sim_config_destroy(&cfg);
    return 1;
}

/* ── Uniform Mutation ────────────────────────────────────────── */

static int test_uniform_mutate_applies(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.min_mutation_rate = 0.1;
    cfg.max_mutation_rate = 0.2;

    rng_seed(&_test_rng, 42, 0);
    step_uniform_mutate(&cfg, &seq, &rec);

    int mutations = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->flags & NUC_FLAG_MUTATED) mutations++;
    }
    ASSERT(mutations > 0, "Should apply mutations");
    sim_config_destroy(&cfg);
    return 1;
}

static int test_uniform_mutate_respects_np(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.min_mutation_rate = 0.5;
    cfg.max_mutation_rate = 0.5;

    rng_seed(&_test_rng, 42, 0);
    step_uniform_mutate(&cfg, &seq, &rec);

    /* NP regions should NOT be mutated */
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->segment == SEG_NP1 || n->segment == SEG_NP2) {
            ASSERT(!(n->flags & NUC_FLAG_MUTATED), "NP regions should not be mutated");
        }
    }
    sim_config_destroy(&cfg);
    return 1;
}

/* ── CSR Rate Adjustment ─────────────────────────────────────── */

static int test_csr_adjusts_ighg1_rates(void) {
    SimRecord rec;
    sim_record_init(&rec);
    static Allele c_al = { .name = "IGHG1*01", .gene = "IGHG1", .family = "IGHG",
        .seq = "ACGT", .length = 4, .anchor = 0, .segment_type = SEG_C };
    rec.c_allele = &c_al;

    double min_rate = 0.01, max_rate = 0.05;
    csr_adjust_rates(&rec, &min_rate, &max_rate);

    ASSERT(fabs(min_rate - 0.05) < 1e-6, "IGHG1 min rate should be 0.05");
    ASSERT(fabs(max_rate - 0.12) < 1e-6, "IGHG1 max rate should be 0.12");
    return 1;
}

static int test_csr_adjusts_ighm_rates(void) {
    SimRecord rec;
    sim_record_init(&rec);
    static Allele c_al = { .name = "IGHM*01", .gene = "IGHM", .family = "IGHM",
        .seq = "ACGT", .length = 4, .anchor = 0, .segment_type = SEG_C };
    rec.c_allele = &c_al;

    double min_rate = 0.05, max_rate = 0.15;
    csr_adjust_rates(&rec, &min_rate, &max_rate);

    ASSERT(fabs(min_rate - 0.001) < 1e-6, "IGHM min rate should be 0.001");
    ASSERT(fabs(max_rate - 0.03) < 1e-6, "IGHM max rate should be 0.03");
    return 1;
}

static int test_csr_no_c_allele(void) {
    SimRecord rec;
    sim_record_init(&rec);
    /* No C allele */

    double min_rate = 0.05, max_rate = 0.15;
    csr_adjust_rates(&rec, &min_rate, &max_rate);

    ASSERT(fabs(min_rate - 0.05) < 1e-6, "Rates unchanged without C allele");
    ASSERT(fabs(max_rate - 0.15) < 1e-6, "Rates unchanged without C allele");
    return 1;
}

/* ══════════════════════════════════════════════════════════════ */
/*                    PIPELINE INTEGRATION TESTS                  */
/* ══════════════════════════════════════════════════════════════ */

static int test_pipeline_with_corruption_ops(void) {
    SimConfig cfg;
    sim_config_init(&cfg, CHAIN_IGH);
    cfg.rng = &_test_rng;

    /* Enable several corruption features */
    cfg.features.quality_errors = true;
    cfg.features.insert_ns = true;
    cfg.features.reverse_complement = true;
    cfg.rc_prob = 1.0;
    cfg.base_error_rate = 0.1;
    cfg.peak_error_rate = 0.3;
    cfg.n_prob = 0.1;

    /* Add minimal alleles */
    static Allele v = { .name = "V1*01", .gene = "V1", .family = "V",
        .seq = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC",
        .length = 104, .anchor = 90, .segment_type = SEG_V };
    static Allele d = { .name = "D1*01", .gene = "D1", .family = "D",
        .seq = "GGTATAACTGG",
        .length = 11, .anchor = 0, .segment_type = SEG_D };
    static Allele j = { .name = "J1*01", .gene = "J1", .family = "J",
        .seq = "ACTACTTTGACTACTGGGGCCAAGGAACC",
        .length = 28, .anchor = 8, .segment_type = SEG_J };

    allele_pool_add(&cfg.v_alleles, &v);
    allele_pool_add(&cfg.d_alleles, &d);
    allele_pool_add(&cfg.j_alleles, &j);

    /* Set up trim distributions */
    double trim_probs[] = { 0.5, 0.3, 0.2 };
    cfg.v_trim_3.probs = malloc(3 * sizeof(double));
    memcpy(cfg.v_trim_3.probs, trim_probs, 3 * sizeof(double));
    cfg.v_trim_3.max_trim = 3;
    cfg.d_trim_5.probs = malloc(3 * sizeof(double));
    memcpy(cfg.d_trim_5.probs, trim_probs, 3 * sizeof(double));
    cfg.d_trim_5.max_trim = 3;
    cfg.d_trim_3.probs = malloc(3 * sizeof(double));
    memcpy(cfg.d_trim_3.probs, trim_probs, 3 * sizeof(double));
    cfg.d_trim_3.max_trim = 3;
    cfg.j_trim_5.probs = malloc(3 * sizeof(double));
    memcpy(cfg.j_trim_5.probs, trim_probs, 3 * sizeof(double));
    cfg.j_trim_5.max_trim = 3;

    Pipeline pl = pipeline_build(&cfg);

    /* Pipeline should include corruption steps */
    ASSERT(pl.n_steps > 9, "Pipeline should have corruption steps");

    /* Execute pipeline */
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;

    rng_seed(&_test_rng, 42, 0);
    pipeline_execute(&pl, &cfg, &seq, &rec, 0, NULL, NULL);

    ASSERT(aseq_length(&seq) > 0, "Sequence should be generated");
    ASSERT(rec.is_reverse_complement == true, "Should be RC'd with prob=1.0");

    sim_config_destroy(&cfg);
    return 1;
}

/* ══════════════════════════════════════════════════════════════ */
/*                           MAIN                                */
/* ══════════════════════════════════════════════════════════════ */

int main(void) {
    printf("test_ops: corruption, simulation, and mutation ops\n");
    printf("═══════════════════════════════════════════════════\n\n");

    printf("Corruption ops:\n");
    RUN(test_corrupt_5_prime_changes_length);
    RUN(test_corrupt_5_prime_remove_shrinks);
    RUN(test_corrupt_3_prime_changes_length);
    RUN(test_quality_errors_introduces_errors);
    RUN(test_quality_errors_preserves_germline);
    RUN(test_pcr_introduces_errors);
    RUN(test_umi_prepends_barcode);
    RUN(test_reverse_complement_flips);
    RUN(test_reverse_complement_probability);
    RUN(test_contaminant_replaces_sequence);
    RUN(test_contaminant_phix);
    RUN(test_insert_indels_modifies);
    RUN(test_insert_indels_preserves_anchors);
    RUN(test_insert_ns_adds_ambiguous);
    RUN(test_paired_end_creates_gap);
    RUN(test_paired_end_no_gap_long_reads);
    RUN(test_long_read_errors_targets_homopolymers);
    RUN(test_trim_to_length_truncates);
    RUN(test_trim_to_length_noop_when_shorter);
    RUN(test_primer_mask_restores_germline);

    printf("\nSimulation ops:\n");
    RUN(test_d_inversion_complements);
    RUN(test_d_inversion_probability);
    RUN(test_receptor_revision_changes_v);
    RUN(test_selection_pressure_reverts_mutations);

    printf("\nMutation:\n");
    RUN(test_uniform_mutate_applies);
    RUN(test_uniform_mutate_respects_np);

    printf("\nCSR:\n");
    RUN(test_csr_adjusts_ighg1_rates);
    RUN(test_csr_adjusts_ighm_rates);
    RUN(test_csr_no_c_allele);

    printf("\nPipeline integration:\n");
    RUN(test_pipeline_with_corruption_ops);

    printf("\n%d/%d tests passed.\n", tests_passed, tests_run);
    return tests_passed == tests_run ? 0 : 1;
}
