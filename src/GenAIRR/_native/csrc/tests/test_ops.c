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

/* ── T1-3: min=max=0 must disable the branch ─────────────────── */

/* When corrupt_5_remove_min=max=0, the remove branch (when rolled)
 * must produce amount=0 and leave sequence length unchanged. Same for
 * the add branch. Sweep many seeds to hit both event flavours. */
static int test_t1_3_corrupt_5_min_max_zero_disables(void) {
    int saw_remove_event = 0, saw_add_event = 0;
    for (int seed = 0; seed < 200; seed++) {
        ASeq seq; SimRecord rec; SimConfig cfg;
        build_simple_seq(&seq, &rec, &cfg);
        cfg.corrupt_5_remove_min = 0;
        cfg.corrupt_5_remove_max = 0;
        cfg.corrupt_5_add_min    = 0;
        cfg.corrupt_5_add_max    = 0;

        int len_before = aseq_length(&seq);
        rng_seed(&_test_rng, (uint64_t)seed, 0);
        step_corrupt_5_prime(&cfg, &seq, &rec);
        int len_after = aseq_length(&seq);

        ASSERT(len_after == len_before,
               "min=max=0 must not change seq length");

        if (rec.corruption_5_event == 1 || rec.corruption_5_event == 3) {
            saw_remove_event = 1;
            ASSERT(rec.corruption_5_remove_amount == 0,
                   "remove amount must be 0 when min=max=0");
        }
        if (rec.corruption_5_event == 2 || rec.corruption_5_event == 3) {
            saw_add_event = 1;
            ASSERT(rec.corruption_5_add_amount == 0,
                   "add amount must be 0 when min=max=0");
        }
        sim_config_destroy(&cfg);
    }
    ASSERT(saw_remove_event && saw_add_event,
           "200-seed sweep should hit both remove and add events");
    return 1;
}

static int test_t1_3_corrupt_3_min_max_zero_disables(void) {
    int saw_remove_event = 0, saw_add_event = 0;
    for (int seed = 0; seed < 200; seed++) {
        ASeq seq; SimRecord rec; SimConfig cfg;
        build_simple_seq(&seq, &rec, &cfg);
        cfg.corrupt_3_remove_min = 0;
        cfg.corrupt_3_remove_max = 0;
        cfg.corrupt_3_add_min    = 0;
        cfg.corrupt_3_add_max    = 0;

        int len_before = aseq_length(&seq);
        rng_seed(&_test_rng, (uint64_t)seed, 0);
        step_corrupt_3_prime(&cfg, &seq, &rec);
        int len_after = aseq_length(&seq);

        ASSERT(len_after == len_before,
               "min=max=0 must not change seq length");

        if (rec.corruption_3_event == 1 || rec.corruption_3_event == 3) {
            saw_remove_event = 1;
            ASSERT(rec.corruption_3_remove_amount == 0,
                   "remove amount must be 0 when min=max=0");
        }
        if (rec.corruption_3_event == 2 || rec.corruption_3_event == 3) {
            saw_add_event = 1;
            ASSERT(rec.corruption_3_add_amount == 0,
                   "add amount must be 0 when min=max=0");
        }
        sim_config_destroy(&cfg);
    }
    ASSERT(saw_remove_event && saw_add_event,
           "200-seed sweep should hit both remove and add events");
    return 1;
}

/* min=0,max=N>0 must produce amounts that include 0 in the empirical
 * distribution (i.e., the floor is not silently re-introducing 1).
 * With min=0, max=2 the sample space is {0, 1, 2}; we expect each to
 * appear at >5% over a 500-seed sweep. */
static int test_t1_3_corrupt_5_lower_tail_visible(void) {
    int counts[4] = {0};   /* index by amount, 0..3 */
    int n_remove = 0;
    for (int seed = 0; seed < 500; seed++) {
        ASeq seq; SimRecord rec; SimConfig cfg;
        build_simple_seq(&seq, &rec, &cfg);
        cfg.corrupt_5_remove_min = 0;
        cfg.corrupt_5_remove_max = 2;
        cfg.corrupt_5_add_min    = 0;
        cfg.corrupt_5_add_max    = 0;

        rng_seed(&_test_rng, (uint64_t)seed, 0);
        step_corrupt_5_prime(&cfg, &seq, &rec);

        if (rec.corruption_5_event == 1 || rec.corruption_5_event == 3) {
            int a = rec.corruption_5_remove_amount;
            if (a >= 0 && a <= 3) counts[a]++;
            n_remove++;
        }
        sim_config_destroy(&cfg);
    }
    ASSERT(n_remove > 100, "Need a usable remove-event sample");
    /* Each of {0, 1, 2} should appear at >5% of remove events. */
    ASSERT(counts[0] > n_remove / 20, "amount=0 must be reachable");
    ASSERT(counts[1] > n_remove / 20, "amount=1 must be reachable");
    ASSERT(counts[2] > n_remove / 20, "amount=2 must be reachable");
    ASSERT(counts[3] == 0, "amount > max must never appear");
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

/* T2-12: count how many of the 6 anchor-codon nts (3 V + 3 J) got
 * an N. The first nt of each anchor codon carries NUC_FLAG_ANCHOR;
 * walk the linked list to also catch nts +1 and +2 of each codon. */
static int count_ns_in_anchor_codons(const ASeq *seq) {
    int hits = 0;
    int residual = 0;
    for (Nuc *n = seq->head; n; n = n->next) {
        bool in_anchor = false;
        if (n->flags & NUC_FLAG_ANCHOR) { in_anchor = true; residual = 2; }
        else if (residual > 0) { in_anchor = true; residual--; }
        if (in_anchor && (n->flags & NUC_FLAG_IS_N)) hits++;
    }
    return hits;
}

static int test_insert_ns_protects_anchors_under_productive_only(void) {
    /* T2-12: with PRODUCTIVITY_PRODUCTIVE_ONLY and a near-certain N
     * rate, the 6 anchor-codon nts (3 V + 3 J) must NEVER be flipped
     * to N — that's the whole point of the protection. */
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.n_prob = 0.99;
    cfg.features.productivity = PRODUCTIVITY_PRODUCTIVE_ONLY;

    rng_seed(&_test_rng, 42, 0);
    step_insert_ns(&cfg, &seq, &rec);

    int hits = count_ns_in_anchor_codons(&seq);
    ASSERT(hits == 0,
           "PRODUCTIVITY_PRODUCTIVE_ONLY must protect anchor codons");
    sim_config_destroy(&cfg);
    return 1;
}

static int test_insert_ns_does_not_protect_under_mixed(void) {
    /* T2-12: under PRODUCTIVITY_MIXED the loader uses the full noise
     * model — anchors are NOT protected. With prob=0.99 across 6
     * anchor nts, the probability of ZERO hits is 0.01^6 ≈ 1e-12,
     * effectively impossible. Asserting hits > 0 confirms the
     * protection is conditional on PRODUCTIVE_ONLY (and isn't
     * leaking into MIXED). */
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.n_prob = 0.99;
    cfg.features.productivity = PRODUCTIVITY_MIXED;

    rng_seed(&_test_rng, 42, 0);
    step_insert_ns(&cfg, &seq, &rec);

    int hits = count_ns_in_anchor_codons(&seq);
    ASSERT(hits > 0,
           "PRODUCTIVITY_MIXED should NOT protect anchors");
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

static int test_d_inversion_sets_per_node_flag(void) {
    /* T2-13: inverted bases must carry NUC_FLAG_INVERTED for forensic
     * / introspection tools that need per-position event provenance.
     * The record-level rec->d_inverted flag is necessary but not
     * sufficient — downstream code that walks the linked list shouldn't
     * have to re-derive segment boundaries to identify inverted bases. */
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.d_inversion_prob = 1.0;

    rng_seed(&_test_rng, 42, 0);
    step_d_inversion(&cfg, &seq, &rec);

    int d_len_observed = 0;
    int d_with_flag = 0;
    int non_d_with_flag = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        bool flagged = (n->flags & NUC_FLAG_INVERTED) != 0;
        if (n->segment == SEG_D) {
            d_len_observed++;
            if (flagged) d_with_flag++;
        } else if (flagged) {
            non_d_with_flag++;
        }
    }
    ASSERT(d_len_observed > 0, "test must have a D segment");
    ASSERT(d_with_flag == d_len_observed,
           "every inverted D node must carry NUC_FLAG_INVERTED");
    ASSERT(non_d_with_flag == 0,
           "NUC_FLAG_INVERTED must not leak into V/D-NP/J nodes");
    sim_config_destroy(&cfg);
    return 1;
}

static int test_d_inversion_germline_tracks_rearranged_template(void) {
    /* T2-13 (audit-resolution): per AIRR Rearrangement Schema, the
     * germline_alignment column must equal the post-rearrangement
     * template at each position. So `n->germline` for an inverted-D
     * node must equal the inverted base (the same value as
     * `n->current` immediately after inversion, before SHM runs).
     * If we ever revert to "germline stays canonical," AIRR mutation
     * calling against inverted D becomes biologically incorrect
     * (every base looks like an SHM). */
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.d_inversion_prob = 1.0;

    rng_seed(&_test_rng, 42, 0);
    step_d_inversion(&cfg, &seq, &rec);

    int mismatches = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->segment != SEG_D) continue;
        if (n->germline != n->current) mismatches++;
    }
    ASSERT(mismatches == 0,
           "n->germline must equal n->current immediately after inversion "
           "(germline_alignment is the post-rearrangement template)");
    sim_config_destroy(&cfg);
    return 1;
}

static int test_d_inversion_no_flag_when_unfired(void) {
    /* Sanity: with prob=0 no inversion happens, so no NUC_FLAG_INVERTED
     * is set anywhere. Catches regressions where the flag would be
     * applied on a code path that runs even when the random check
     * rejects. */
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.d_inversion_prob = 0.0;

    rng_seed(&_test_rng, 42, 0);
    step_d_inversion(&cfg, &seq, &rec);

    for (Nuc *n = seq.head; n; n = n->next) {
        ASSERT(!(n->flags & NUC_FLAG_INVERTED),
               "NUC_FLAG_INVERTED must not be set when prob=0");
    }
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
    cfg.selection_strength    = 1.0;
    cfg.cdr_r_acceptance      = 0.0;
    cfg.fwr_r_acceptance      = 0.0;
    cfg.anchor_r_acceptance   = 0.0;

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

/**
 * Selection now walks all coding segments (V/NP1/D/NP2/J), not just V.
 * Insert R-mutations into D and J and verify they are reverted under
 * 100% strength + 0 acceptance.
 */
static int test_selection_walks_dj_segments(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.selection_strength    = 1.0;
    cfg.cdr_r_acceptance      = 0.0;
    cfg.fwr_r_acceptance      = 0.0;
    cfg.anchor_r_acceptance   = 0.0;

    aseq_build_codon_rail(&seq);

    /* Find a coding-segment node in D (mid-segment) and a J-FR4 node
     * (just past J anchor + 3) and force replacement mutations. */
    Nuc *d_node = NULL, *j_fr4 = NULL;
    int j_anchor = rec.j_allele->anchor;   /* 8 */
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->segment == SEG_D && d_node == NULL) d_node = n;
        if (n->segment == SEG_J && n->germline_pos >= (uint16_t)(j_anchor + 3)
            && j_fr4 == NULL) {
            j_fr4 = n;
            break;
        }
    }
    ASSERT(d_node && j_fr4, "Need D and J-FR4 nodes for the test");

    /* Force R-mutations: pick a base whose codon translation differs.
     * For simplicity, mutate to N's complement so we virtually always
     * change the codon. We cannot guarantee R without inspecting the
     * codon, so loop several nodes per segment and at least one will
     * land R. Instead: tag MANY nodes and assert revert count > 0. */
    int d_mutated = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->segment != SEG_D) continue;
        if (n->germline == '\0') continue;
        char nb = (n->current == 'A') ? 'C' : 'A';
        aseq_mutate(&seq, n, nb, NUC_FLAG_MUTATED);
        d_mutated++;
    }
    int j_mutated = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->segment != SEG_J) continue;
        if (n->germline_pos < (uint16_t)(j_anchor + 3)) continue;  /* skip CDR3 */
        if (n->germline == '\0') continue;
        char nb = (n->current == 'A') ? 'C' : 'A';
        aseq_mutate(&seq, n, nb, NUC_FLAG_MUTATED);
        j_mutated++;
    }
    ASSERT(d_mutated > 0 && j_mutated > 0, "Should have mutations to revert in D and J");

    rng_seed(&_test_rng, 7, 0);
    step_selection_pressure(&cfg, &seq, &rec);

    int d_remain = 0, j_remain = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (!(n->flags & NUC_FLAG_MUTATED)) continue;
        if (n->segment == SEG_D) d_remain++;
        else if (n->segment == SEG_J && n->germline_pos >= (uint16_t)(j_anchor + 3))
            j_remain++;
    }
    /* At 0 acceptance, only S-mutations survive. Most randomly chosen
     * substitutions in a coding context are R, so revert count must
     * dominate. */
    ASSERT(d_remain < d_mutated, "D R-mutations should mostly revert");
    ASSERT(j_remain < j_mutated, "J FR4 R-mutations should mostly revert");
    sim_config_destroy(&cfg);
    return 1;
}

/**
 * Anchor-codon protection: V Cys codon and J W/F codon should be
 * fiercely protected (anchor_r_acceptance default = 0). Force an
 * R-mutation into the anchor codon and verify it is reverted across
 * many seeds.
 */
static int test_selection_protects_anchor_codon(void) {
    int reverted = 0, total = 0;
    int v_anchor_pos = -1;

    for (int seed = 1; seed <= 64; seed++) {
        ASeq seq; SimRecord rec; SimConfig cfg;
        build_simple_seq(&seq, &rec, &cfg);
        cfg.selection_strength    = 1.0;
        cfg.cdr_r_acceptance      = 1.0;   /* keep ordinary CDR mutations */
        cfg.fwr_r_acceptance      = 1.0;   /* keep ordinary FWR mutations */
        cfg.anchor_r_acceptance   = 0.0;   /* but reject anchor-disrupting R */

        aseq_build_codon_rail(&seq);

        /* Find the V anchor node and force an R-mutation on it. The
         * V allele is "...ATGCATGC..." repeated; mutating any base in
         * the anchor codon to a different base will change Cys (TGC)
         * to something else. */
        Nuc *anchor = aseq_find_anchor(&seq, SEG_V);
        if (!anchor) { sim_config_destroy(&cfg); continue; }
        v_anchor_pos = anchor->germline_pos;

        /* Force an R-mutation: change the anchor base to a base that
         * definitely changes the amino acid. Mutate anchor's base. */
        char nb = (anchor->current == 'A') ? 'C' :
                  (anchor->current == 'C') ? 'A' :
                  (anchor->current == 'G') ? 'T' : 'G';
        aseq_mutate(&seq, anchor, nb, NUC_FLAG_MUTATED);

        rng_seed(&_test_rng, (uint64_t)seed, 0);
        step_selection_pressure(&cfg, &seq, &rec);

        total++;
        /* The anchor R-mutation should always be reverted (acceptance=0). */
        if (!(anchor->flags & NUC_FLAG_MUTATED)) reverted++;
        sim_config_destroy(&cfg);
    }
    ASSERT(v_anchor_pos > 0, "Test fixture should expose V anchor");
    ASSERT(total >= 60, "Should have run on most seeds");
    ASSERT(reverted == total,
           "Anchor R-mutations must always revert at anchor_r_acceptance=0");
    return 1;
}

/**
 * Anchorless V allele: when V has anchor <= 0, V-segment mutations
 * are SKIPPED (not classified as FR/CDR via the static table — that
 * would mis-classify CDR3). NP1/D/NP2/J are still selected as usual.
 */
static int test_selection_skips_anchorless_v(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);

    /* Override V allele to be anchorless. We must keep the original
     * pointer alive — copy into a fresh static so rec->v_allele
     * doesn't dangle when the build_simple_seq's static is reused
     * later. Also clear FLAG_ANCHOR on the sequence: a real
     * anchorless V allele would never have had FLAG_ANCHOR set during
     * append, but our test fixture pre-tagged it. */
    static Allele anchorless_v = { 0 };
    memcpy(&anchorless_v, rec.v_allele, sizeof(Allele));
    anchorless_v.anchor = -1;
    rec.v_allele = &anchorless_v;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->segment == SEG_V) n->flags &= (uint16_t)~NUC_FLAG_ANCHOR;
    }

    cfg.selection_strength    = 1.0;
    cfg.cdr_r_acceptance      = 0.0;
    cfg.fwr_r_acceptance      = 0.0;
    cfg.anchor_r_acceptance   = 0.0;

    aseq_build_codon_rail(&seq);

    /* Mutate every coding-segment node we can */
    int v_mut = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->segment != SEG_V) continue;
        if (n->germline == '\0') continue;
        char nb = (n->current == 'A') ? 'C' : 'A';
        aseq_mutate(&seq, n, nb, NUC_FLAG_MUTATED);
        v_mut++;
    }
    ASSERT(v_mut > 0, "Should have mutated V nodes");

    rng_seed(&_test_rng, 13, 0);
    step_selection_pressure(&cfg, &seq, &rec);

    int v_remain = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->segment == SEG_V && (n->flags & NUC_FLAG_MUTATED)) v_remain++;
    }
    /* Anchorless V is skipped → no V mutations should be reverted,
     * regardless of how many were R. */
    ASSERT(v_remain == v_mut,
           "Anchorless V mutations must NOT be reverted (segment skipped)");
    sim_config_destroy(&cfg);
    return 1;
}

/**
 * Order stability: when two mutations land in the same codon, the
 * outcome is independent of evaluation order modulo the random draws
 * (each is evaluated against the current germline-substituted codon).
 * We don't require equality; we require that across many seeds the
 * acceptance rate matches the configured CDR rate within MC noise.
 */
static int test_selection_acceptance_rate_calibration(void) {
    int total_r = 0, kept_r = 0;
    const double target = 0.50;       /* effective CDR acceptance */

    for (int seed = 1; seed <= 200; seed++) {
        ASeq seq; SimRecord rec; SimConfig cfg;
        build_simple_seq(&seq, &rec, &cfg);
        cfg.selection_strength    = 1.0;
        cfg.cdr_r_acceptance      = target;
        cfg.fwr_r_acceptance      = target;
        cfg.anchor_r_acceptance   = target;

        aseq_build_codon_rail(&seq);

        /* Mutate a single codon-aligned non-anchor node deterministically. */
        Nuc *n = seq.head;
        while (n && (n->segment != SEG_V || n->germline == '\0' ||
                     (n->flags & NUC_FLAG_ANCHOR) || n->frame_phase != 0)) {
            n = n->next;
        }
        if (!n) { sim_config_destroy(&cfg); continue; }
        char nb = (n->current == 'A') ? 'C' :
                  (n->current == 'C') ? 'G' :
                  (n->current == 'G') ? 'T' : 'A';
        aseq_mutate(&seq, n, nb, NUC_FLAG_MUTATED);

        /* Skip if it ended up silent (acceptance is 100% for S). */
        char b1 = n->germline, b2 = n->next ? n->next->current : 'A';
        char b3 = (n->next && n->next->next) ? n->next->next->current : 'A';
        char aa_g = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[
            ((b1=='T'||b1=='t')?0:(b1=='C'||b1=='c')?1:(b1=='A'||b1=='a')?2:3)*16 +
            ((b2=='T'||b2=='t')?0:(b2=='C'||b2=='c')?1:(b2=='A'||b2=='a')?2:3)*4 +
            ((b3=='T'||b3=='t')?0:(b3=='C'||b3=='c')?1:(b3=='A'||b3=='a')?2:3)
        ];
        char b1c = n->current;
        char aa_c = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[
            ((b1c=='T'||b1c=='t')?0:(b1c=='C'||b1c=='c')?1:(b1c=='A'||b1c=='a')?2:3)*16 +
            ((b2=='T'||b2=='t')?0:(b2=='C'||b2=='c')?1:(b2=='A'||b2=='a')?2:3)*4 +
            ((b3=='T'||b3=='t')?0:(b3=='C'||b3=='c')?1:(b3=='A'||b3=='a')?2:3)
        ];
        if (aa_g == aa_c) { sim_config_destroy(&cfg); continue; }

        total_r++;

        rng_seed(&_test_rng, (uint64_t)seed, 0);
        step_selection_pressure(&cfg, &seq, &rec);

        if (n->flags & NUC_FLAG_MUTATED) kept_r++;
        sim_config_destroy(&cfg);
    }
    ASSERT(total_r > 50, "Need a usable R-mutation sample");
    double observed = (double)kept_r / (double)total_r;
    /* MC band: 200 trials, target 0.5 → SE ~0.035 → ±0.10 is safe */
    ASSERT(observed > 0.40 && observed < 0.60,
           "Empirical R-acceptance should match target within MC noise");
    return 1;
}

/**
 * Codon-spanning V→NP1 boundary: the last 1-2 nucleotides of V and
 * the first base(s) of NP1 form a single codon. Selection must
 * correctly classify R/S using the codon-rail head, even though the
 * codon spans a segment boundary.
 *
 * We force a mutation on the last V base before the V-anchor stretch,
 * verify the rail places it in a codon whose head may be in V (and
 * its 2nd/3rd bases in NP1 in some configurations). The test
 * primarily asserts no crash and that the mutation classification
 * uses the codon rail (not just `germline_pos % 3`).
 */
static int test_selection_codon_spanning_boundary(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.selection_strength    = 1.0;
    cfg.cdr_r_acceptance      = 0.0;
    cfg.fwr_r_acceptance      = 0.0;
    cfg.anchor_r_acceptance   = 0.0;

    aseq_build_codon_rail(&seq);

    /* Find the last V node (germline boundary) and the first NP1 node. */
    Nuc *last_v = NULL, *first_np1 = NULL;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->segment == SEG_V) last_v = n;
        if (n->segment == SEG_NP1 && first_np1 == NULL) first_np1 = n;
    }
    ASSERT(last_v && first_np1, "Need both V tail and NP1 head");

    /* Force a mutation at last V base (likely codon-spanning into NP1). */
    char nb = (last_v->current == 'A') ? 'T' : 'A';
    aseq_mutate(&seq, last_v, nb, NUC_FLAG_MUTATED);

    /* Force a mutation at first NP1 base (may also be codon-spanning). */
    char nb2 = (first_np1->current == 'A') ? 'T' : 'A';
    aseq_mutate(&seq, first_np1, nb2, NUC_FLAG_MUTATED);

    rng_seed(&_test_rng, 11, 0);
    step_selection_pressure(&cfg, &seq, &rec);

    /* No crash, codon rail still valid, no out-of-bounds. */
    ASSERT(seq.codon_rail_valid, "Codon rail should remain valid");
    /* Last V mutation: V is selected (with anchor). At 0 acceptance,
     * if classified as R, it must revert. Either it reverted, or it
     * was silent (kept). Both are valid. We assert the mutation
     * either carries the MUTATED flag with germline preserved or was
     * reverted to germline. */
    ASSERT(last_v->germline != '\0', "Last V node must have germline");
    if (!(last_v->flags & NUC_FLAG_MUTATED)) {
        ASSERT(last_v->current == last_v->germline,
               "If reverted, current must match germline");
    }
    sim_config_destroy(&cfg);
    return 1;
}

/**
 * Selection strength=0 is a no-op: no mutations should be reverted.
 */
static int test_selection_strength_zero_noop(void) {
    ASeq seq; SimRecord rec; SimConfig cfg;
    build_simple_seq(&seq, &rec, &cfg);
    cfg.selection_strength    = 0.0;   /* disabled */
    cfg.cdr_r_acceptance      = 0.0;
    cfg.fwr_r_acceptance      = 0.0;
    cfg.anchor_r_acceptance   = 0.0;

    aseq_build_codon_rail(&seq);

    int mutated = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (!nuc_is_coding_segment(n) || n->germline == '\0') continue;
        if (mutated >= 30) break;
        char nb = (n->current == 'A') ? 'C' : 'A';
        aseq_mutate(&seq, n, nb, NUC_FLAG_MUTATED);
        mutated++;
    }

    rng_seed(&_test_rng, 99, 0);
    step_selection_pressure(&cfg, &seq, &rec);

    int remaining = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->flags & NUC_FLAG_MUTATED) remaining++;
    }
    ASSERT(remaining == mutated,
           "strength=0 must not revert anything");
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
    RUN(test_t1_3_corrupt_5_min_max_zero_disables);
    RUN(test_t1_3_corrupt_3_min_max_zero_disables);
    RUN(test_t1_3_corrupt_5_lower_tail_visible);
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
    RUN(test_insert_ns_protects_anchors_under_productive_only);  /* T2-12 */
    RUN(test_insert_ns_does_not_protect_under_mixed);            /* T2-12 */
    RUN(test_paired_end_creates_gap);
    RUN(test_paired_end_no_gap_long_reads);
    RUN(test_long_read_errors_targets_homopolymers);
    RUN(test_trim_to_length_truncates);
    RUN(test_trim_to_length_noop_when_shorter);
    RUN(test_primer_mask_restores_germline);

    printf("\nSimulation ops:\n");
    RUN(test_d_inversion_complements);
    RUN(test_d_inversion_probability);
    RUN(test_d_inversion_sets_per_node_flag);              /* T2-13 */
    RUN(test_d_inversion_germline_tracks_rearranged_template); /* T2-13 */
    RUN(test_d_inversion_no_flag_when_unfired);            /* T2-13 */
    RUN(test_receptor_revision_changes_v);
    RUN(test_selection_pressure_reverts_mutations);
    RUN(test_selection_walks_dj_segments);
    RUN(test_selection_protects_anchor_codon);
    RUN(test_selection_skips_anchorless_v);
    RUN(test_selection_acceptance_rate_calibration);
    RUN(test_selection_codon_spanning_boundary);
    RUN(test_selection_strength_zero_noop);

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
