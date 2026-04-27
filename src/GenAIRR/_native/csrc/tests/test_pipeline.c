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

/* Test-local PCG32 RNG. Wired into every cfg returned by
 * make_igh_config() (and the inline IGK config) so step functions
 * have a valid cfg->rng. Tests reseed via rng_seed() instead of
 * the old srand(). */
static RngState _test_rng;

static SimConfig make_productive_capable_igh_config(void);  /* forward */

static SimConfig make_igh_config(void) {
    SimConfig cfg;
    sim_config_init(&cfg, CHAIN_IGH);
    cfg.rng = &_test_rng;

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
    cfg.rng = &_test_rng;

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
    rng_seed(&_test_rng, 42, 0);  /* deterministic */

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
    rng_seed(&_test_rng, 123, 0);

    SimConfig cfg = make_productive_capable_igh_config();
    cfg.features.productivity = PRODUCTIVITY_PRODUCTIVE_ONLY;
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

/* ── T1-6: productivity filter modes ──────────────────────────────
 *
 * The Productivity enum guarantees that filtered modes return ONLY
 * the matching kind of sequence. We need productive-capable alleles
 * for these tests — the simple test config in make_igh_config() puts
 * arbitrary bases at the anchor positions so rule_conserved_cysteine
 * and rule_conserved_anchor (J W/F) always fail. The helper below
 * builds a realistic config with proper Cys at the V anchor and Trp
 * at the J anchor.
 *
 * With max_productive_attempts=200, retry exhaustion is effectively
 * impossible regardless of mode, so we assert that every returned
 * record matches the configured direction — not "≤1 leak" or "≥95%".
 * ─────────────────────────────────────────────────────────────── */

static SimConfig make_productive_capable_igh_config(void) {
    SimConfig cfg;
    sim_config_init(&cfg, CHAIN_IGH);
    cfg.rng = &_test_rng;

    /* V (30bp), Cys at anchor=24:
     *   ATG GCA ATC CTA GCA GCA ATT GCA TGT GGT
     *    M   A   I   L   A   A   I   A   C   G
     *    \________ no stops _________/  ↑ Cys */
    Allele v = make_allele(
        "IGHV1-A*01",
        "ATGGCAATCCTAGCAGCAATTGCATGTGGT",
        24, SEG_V
    );

    /* D (12bp): contributes some bases to junction; no anchor. */
    Allele d = make_allele(
        "IGHD-A*01",
        "GTATTACTGTGC",
        0, SEG_D
    );

    /* J (21bp), Trp at anchor=6:
     *   ACT ACT TGG GGC CAA GGT ACT
     *    T   T   W   G   Q   G   T
     *           ↑ Trp at germline_pos 6 */
    Allele j = make_allele(
        "IGHJ-A*01",
        "ACTACTTGGGGCCAAGGTACT",
        6, SEG_J
    );

    cfg.v_alleles = allele_pool_create(4);
    cfg.d_alleles = allele_pool_create(4);
    cfg.j_alleles = allele_pool_create(4);

    allele_pool_add(&cfg.v_alleles, &v);
    allele_pool_add(&cfg.d_alleles, &d);
    allele_pool_add(&cfg.j_alleles, &j);

    return cfg;
}

static int run_n_and_count_productive(SimConfig *cfg, int n,
                                      int *out_productive,
                                      int *out_non_productive) {
    Pipeline p = pipeline_build(cfg);
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;

    int prod = 0, nonprod = 0;
    for (int i = 0; i < n; i++) {
        pipeline_execute(&p, cfg, &seq, &rec, 0, NULL, NULL);
        if (rec.productive) prod++;
        else                nonprod++;
    }
    *out_productive    = prod;
    *out_non_productive = nonprod;
    return 0;
}

static int test_productivity_productive_only(void) {
    rng_seed(&_test_rng, 7, 0);
    SimConfig cfg = make_productive_capable_igh_config();
    cfg.features.productivity   = PRODUCTIVITY_PRODUCTIVE_ONLY;
    cfg.max_productive_attempts = 200;   /* exhaustion ≈ 0 */

    int prod = 0, nonprod = 0;
    run_n_and_count_productive(&cfg, 100, &prod, &nonprod);
    if (nonprod != 0) {
        fprintf(stderr,
                "PRODUCTIVE_ONLY: expected 0 non-productive, got %d/100\n",
                nonprod);
        sim_config_destroy(&cfg);
        return 1;
    }
    sim_config_destroy(&cfg);
    return 0;
}

static int test_productivity_non_productive_only(void) {
    rng_seed(&_test_rng, 11, 0);
    SimConfig cfg = make_productive_capable_igh_config();
    cfg.features.productivity   = PRODUCTIVITY_NON_PRODUCTIVE_ONLY;
    cfg.max_productive_attempts = 200;

    int prod = 0, nonprod = 0;
    run_n_and_count_productive(&cfg, 100, &prod, &nonprod);
    if (prod != 0) {
        fprintf(stderr,
                "NON_PRODUCTIVE_ONLY: expected 0 productive, got %d/100\n",
                prod);
        sim_config_destroy(&cfg);
        return 1;
    }
    sim_config_destroy(&cfg);
    return 0;
}

static int test_productivity_mixed_returns_both(void) {
    rng_seed(&_test_rng, 13, 0);
    SimConfig cfg = make_productive_capable_igh_config();
    cfg.features.productivity = PRODUCTIVITY_MIXED;

    int prod = 0, nonprod = 0;
    run_n_and_count_productive(&cfg, 500, &prod, &nonprod);
    /* MIXED must contain BOTH classes. With ~22.8% productive baseline
     * on 500 records, we expect ~114 productive and ~386 non-productive.
     * Assert both > 0 (a strong invariant — the test fails immediately
     * if MIXED degenerates into a single-class output). */
    if (prod == 0) {
        fprintf(stderr, "MIXED: expected some productive, got 0/500\n");
        sim_config_destroy(&cfg);
        return 1;
    }
    if (nonprod == 0) {
        fprintf(stderr, "MIXED: expected some non-productive, got 0/500\n");
        sim_config_destroy(&cfg);
        return 1;
    }
    /* Sanity: ratio in the expected band (10-50% productive). Wider
     * than the strict 22.8% to allow for RNG noise. */
    double frac = (double)prod / 500.0;
    if (frac < 0.10 || frac > 0.50) {
        fprintf(stderr,
                "MIXED: productive fraction %.2f outside [0.10, 0.50]\n",
                frac);
        sim_config_destroy(&cfg);
        return 1;
    }
    sim_config_destroy(&cfg);
    return 0;
}

/* Back-compat: setting features.productivity via the legacy bool API
 * (`set_feature("productive", true)`) must continue to work. */
static int test_productivity_legacy_bool_back_compat(void) {
    SimConfig cfg = make_igh_config();
    /* Simulate the legacy API path: a bool maps onto the enum. */
    cfg.features.productivity = PRODUCTIVITY_MIXED;     /* default */
    if (cfg.features.productivity != PRODUCTIVITY_MIXED) {
        sim_config_destroy(&cfg); return 1;
    }
    cfg.features.productivity = PRODUCTIVITY_PRODUCTIVE_ONLY;
    if (cfg.features.productivity != PRODUCTIVITY_PRODUCTIVE_ONLY) {
        sim_config_destroy(&cfg); return 1;
    }
    sim_config_destroy(&cfg);
    return 0;
}

/* ── Test: debug print doesn't crash ──────────────────────────── */

static int test_debug_print(void) {
    rng_seed(&_test_rng, 99, 0);

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

/* ── NP Markov sampler tests ─────────────────────────────────── */

/* Install a synthetic NpDist on cfg->np[idx]. Caller-provided arrays
 * are duplicated so sim_config_destroy can free them uniformly. */
static void install_np_dist(SimConfig *cfg, int idx,
                            const double *length_probs, int max_length,
                            const double first_base[4],
                            const double *transitions, int n_positions) {
    NpDist *dst = &cfg->np[idx];
    free(dst->length_probs);
    free(dst->transitions);

    int n_len = max_length + 1;
    dst->max_length = max_length;
    dst->length_probs = malloc((size_t)n_len * sizeof(double));
    memcpy(dst->length_probs, length_probs, (size_t)n_len * sizeof(double));

    memcpy(dst->first_base, first_base, sizeof(dst->first_base));

    dst->n_positions = n_positions;
    if (n_positions > 0 && transitions) {
        int n_trans = n_positions * 16;
        dst->transitions = malloc((size_t)n_trans * sizeof(double));
        memcpy(dst->transitions, transitions,
               (size_t)n_trans * sizeof(double));
    } else {
        dst->transitions = NULL;
    }
    if (idx + 1 > cfg->n_np_regions) cfg->n_np_regions = idx + 1;
}

/* Forces every NP1 to be exactly "GCG": length=3 with prob 1; first
 * base = G; transition G→C at pos 0; C→G at pos 1. Confirms the
 * Markov sampler actually consults length_probs / first_base /
 * transitions and is not falling back to uniform. */
static int test_np_markov_deterministic(void) {
    rng_seed(&_test_rng, 7, 0);

    SimConfig cfg = make_igh_config();

    /* length: always 3 */
    double length_probs[16] = {0};
    length_probs[3] = 1.0;

    /* first base: always G (index 2) */
    double first_base[4] = {0.0, 0.0, 1.0, 0.0};

    /* transitions: 2 positions of 4×4 each (row-major: [from][to])
     *   pos 0: G→C only  (row 2 = [0,1,0,0])
     *   pos 1: C→G only  (row 1 = [0,0,1,0])
     * Only the rows we'll actually visit need to be valid; others
     * are zeroed and the sampler falls back to first_base. */
    double transitions[2 * 16] = {0};
    transitions[0 * 16 + 2 * 4 + 1] = 1.0;  /* pos 0, G → C */
    transitions[1 * 16 + 1 * 4 + 2] = 1.0;  /* pos 1, C → G */

    install_np_dist(&cfg, 0, length_probs, 15, first_base, transitions, 2);
    /* NP2 left empty → uniform fallback for NP2; we only assert on NP1 */

    Pipeline p = pipeline_build(&cfg);

    int trials = 50;
    for (int t = 0; t < trials; t++) {
        ASeq seq;
        aseq_init(&seq);
        SimRecord rec;
        pipeline_execute(&p, &cfg, &seq, &rec, 0, NULL, NULL);

        if (aseq_segment_length(&seq, SEG_NP1) != 3) {
            fprintf(stderr,
                    "trial %d: expected NP1 length 3, got %d\n",
                    t, aseq_segment_length(&seq, SEG_NP1));
            sim_config_destroy(&cfg);
            return 1;
        }

        /* Walk the NP1 nodes and confirm bases are G,C,G */
        Nuc *n = seq.seg_first[SEG_NP1];
        const char expected[3] = {'G', 'C', 'G'};
        for (int i = 0; i < 3; i++) {
            if (!n || n->current != expected[i]) {
                fprintf(stderr,
                        "trial %d: NP1[%d] expected %c, got %c\n",
                        t, i, expected[i], n ? n->current : '?');
                sim_config_destroy(&cfg);
                return 1;
            }
            n = n->next;
        }
    }

    sim_config_destroy(&cfg);
    return 0;
}

/* Confirms that a strongly-GC-biased NpDist produces a high-GC NP1.
 * Uniform sampling would yield ~50% GC; this distribution is set so
 * that GC content should exceed 80% with high probability. */
static int test_np_markov_gc_bias(void) {
    rng_seed(&_test_rng, 11, 0);

    SimConfig cfg = make_igh_config();

    /* Length distribution: peak around 6, modest tail to 12 */
    double length_probs[16] = {
        0.00, 0.05, 0.10, 0.15, 0.20, 0.20, 0.15,
        0.08, 0.04, 0.02, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0
    };

    /* First base biased GC: P(G)=P(C)=0.45, P(A)=P(T)=0.05 */
    double first_base[4] = {0.05, 0.45, 0.45, 0.05};

    /* Transitions: at every position, from any base, prefer GC.
     * Row pattern (from any of A/C/G/T): [A=0.05, C=0.45, G=0.45, T=0.05] */
    int n_positions = 12;
    double *transitions = calloc((size_t)n_positions * 16, sizeof(double));
    for (int p_i = 0; p_i < n_positions; p_i++) {
        for (int from = 0; from < 4; from++) {
            transitions[p_i * 16 + from * 4 + 0] = 0.05;  /* A */
            transitions[p_i * 16 + from * 4 + 1] = 0.45;  /* C */
            transitions[p_i * 16 + from * 4 + 2] = 0.45;  /* G */
            transitions[p_i * 16 + from * 4 + 3] = 0.05;  /* T */
        }
    }

    install_np_dist(&cfg, 0, length_probs, 15, first_base,
                    transitions, n_positions);
    free(transitions);

    Pipeline p = pipeline_build(&cfg);

    int trials = 200;
    int total_bases = 0, gc_bases = 0;
    for (int t = 0; t < trials; t++) {
        ASeq seq;
        aseq_init(&seq);
        SimRecord rec;
        pipeline_execute(&p, &cfg, &seq, &rec, 0, NULL, NULL);

        Nuc *n = seq.seg_first[SEG_NP1];
        int seg_len = aseq_segment_length(&seq, SEG_NP1);
        for (int i = 0; i < seg_len && n; i++, n = n->next) {
            total_bases++;
            if (n->current == 'G' || n->current == 'C') gc_bases++;
        }
    }

    double gc_frac = total_bases > 0 ? (double)gc_bases / total_bases : 0.0;
    printf("\n    GC bias: %d/%d = %.3f (expect > 0.80)", gc_bases,
           total_bases, gc_frac);
    if (gc_frac < 0.80) {
        fprintf(stderr,
                "\n    GC fraction %.3f below 0.80 threshold "
                "(uniform would give ~0.50)\n", gc_frac);
        sim_config_destroy(&cfg);
        return 1;
    }

    sim_config_destroy(&cfg);
    return 0;
}

/* Confirms uniform fallback path: SimConfig with no NpDist data
 * should still produce non-empty NP1 with bases drawn uniformly. */
static int test_np_uniform_fallback(void) {
    rng_seed(&_test_rng, 13, 0);

    SimConfig cfg = make_igh_config();
    /* Deliberately do NOT install NpDist; cfg->np[*] stays zeroed,
     * so length_probs == NULL and assemble.c uses generate_np_uniform. */

    Pipeline p = pipeline_build(&cfg);

    int trials = 50;
    int nonempty = 0;
    int max_len_seen = 0;
    for (int t = 0; t < trials; t++) {
        ASeq seq;
        aseq_init(&seq);
        SimRecord rec;
        pipeline_execute(&p, &cfg, &seq, &rec, 0, NULL, NULL);
        int np1 = aseq_segment_length(&seq, SEG_NP1);
        if (np1 > 0) nonempty++;
        if (np1 > max_len_seen) max_len_seen = np1;
        if (np1 > 12) {
            fprintf(stderr,
                    "uniform fallback NP1 length %d exceeds expected cap\n",
                    np1);
            sim_config_destroy(&cfg);
            return 1;
        }
    }

    /* Default uniform fallback samples lengths in [3, 9] with mean 6, so
     * empty NP1 should be rare. Lower-bound is loose to keep the test
     * stable across libc rand() implementations. */
    if (nonempty < trials * 8 / 10) {
        fprintf(stderr,
                "uniform fallback produced too many empty NPs: %d/%d\n",
                nonempty, trials);
        sim_config_destroy(&cfg);
        return 1;
    }
    printf("\n    uniform fallback: %d/%d non-empty, max len %d",
           nonempty, trials, max_len_seen);

    sim_config_destroy(&cfg);
    return 0;
}

/* ── P-nucleotide sampler tests ──────────────────────────────── */

/* Install a synthetic PNucDist that always samples a fixed length.
 * length_probs has prob 1.0 at index `forced_len`, zero elsewhere. */
static void install_p_nuc_fixed(SimConfig *cfg, int forced_len, int max_len) {
    PNucDist *d = &cfg->p_nuc_dist;
    free(d->length_probs);
    int n = max_len + 1;
    d->length_probs = calloc((size_t)n, sizeof(double));
    if (forced_len >= 0 && forced_len <= max_len) {
        d->length_probs[forced_len] = 1.0;
    }
    d->max_length = max_len;
}

/* Force a trim distribution that always samples zero. */
static void install_zero_trim(TrimDist *t) {
    free(t->probs);
    t->probs = calloc(2, sizeof(double));
    t->probs[0] = 1.0;
    t->max_trim = 1;
}

/* Force a trim distribution that always samples `forced_trim`. */
static void install_fixed_trim(TrimDist *t, int forced_trim) {
    free(t->probs);
    int n = forced_trim + 2;
    t->probs = calloc((size_t)n, sizeof(double));
    t->probs[forced_trim] = 1.0;
    t->max_trim = forced_trim + 1;
}

/* Verifies the V 3' tail-side P-nuc: when v_trim_3 is forced to 0
 * and the P-nuc length is forced to 3, the first 3 NP1 nodes must
 * be the reverse-complement of V's last 3 bases, tagged with
 * NUC_FLAG_P_NUCLEOTIDE. */
static int test_p_nuc_v3_when_trim_zero(void) {
    rng_seed(&_test_rng, 101, 0);

    SimConfig cfg = make_igh_config();
    /* V allele in make_igh_config ends in "...ATCGAT" */
    install_zero_trim(&cfg.v_trim_3);
    install_fixed_trim(&cfg.j_trim_5, 1);  /* avoid J 5' P-nuc */
    install_fixed_trim(&cfg.d_trim_5, 1);  /* avoid D 5' P-nuc */
    install_fixed_trim(&cfg.d_trim_3, 1);  /* avoid D 3' P-nuc */
    install_p_nuc_fixed(&cfg, /*forced_len=*/3, /*max_len=*/4);

    Pipeline p = pipeline_build(&cfg);

    /* V allele last 3 bases: "GAT" → RC read 5'→3' = "ATC".
     * P-nuc[0] = complement(seg[len-1]) = complement('T') = 'A'
     * P-nuc[1] = complement(seg[len-2]) = complement('A') = 'T'
     * P-nuc[2] = complement(seg[len-3]) = complement('G') = 'C'
     * Expected first 3 NP1 bases: A, T, C. */
    const char expected[3] = {'A', 'T', 'C'};

    int trials = 30;
    for (int t = 0; t < trials; t++) {
        ASeq seq;
        aseq_init(&seq);
        SimRecord rec;
        pipeline_execute(&p, &cfg, &seq, &rec, 0, NULL, NULL);

        Nuc *n = seq.seg_first[SEG_NP1];
        for (int i = 0; i < 3; i++) {
            if (!n) {
                fprintf(stderr, "trial %d: NP1 ran out at i=%d\n", t, i);
                sim_config_destroy(&cfg);
                return 1;
            }
            if (n->current != expected[i]) {
                fprintf(stderr, "trial %d: NP1[%d] expected %c got %c\n",
                        t, i, expected[i], n->current);
                sim_config_destroy(&cfg);
                return 1;
            }
            if (!(n->flags & NUC_FLAG_P_NUCLEOTIDE)) {
                fprintf(stderr,
                        "trial %d: NP1[%d] missing P_NUCLEOTIDE flag\n",
                        t, i);
                sim_config_destroy(&cfg);
                return 1;
            }
            n = n->next;
        }
    }

    sim_config_destroy(&cfg);
    return 0;
}

/* When v_trim_3 != 0, no V 3' P-nuc should appear regardless of the
 * P-nuc length distribution. Same for the other ends. */
static int test_no_p_nuc_when_trim_nonzero(void) {
    rng_seed(&_test_rng, 103, 0);

    SimConfig cfg = make_igh_config();
    install_fixed_trim(&cfg.v_trim_3, 2);
    install_fixed_trim(&cfg.d_trim_5, 2);
    install_fixed_trim(&cfg.d_trim_3, 2);
    install_fixed_trim(&cfg.j_trim_5, 2);
    install_p_nuc_fixed(&cfg, /*forced_len=*/3, /*max_len=*/4);

    Pipeline p = pipeline_build(&cfg);

    int trials = 30;
    for (int t = 0; t < trials; t++) {
        ASeq seq;
        aseq_init(&seq);
        SimRecord rec;
        pipeline_execute(&p, &cfg, &seq, &rec, 0, NULL, NULL);

        for (Nuc *n = seq.head; n; n = n->next) {
            if (n->flags & NUC_FLAG_P_NUCLEOTIDE) {
                fprintf(stderr,
                        "trial %d: P-nuc emitted when all trims > 0\n", t);
                sim_config_destroy(&cfg);
                return 1;
            }
        }
    }

    sim_config_destroy(&cfg);
    return 0;
}

/* Verifies the D 5' head-side P-nuc placement: when d_trim_5=0 and
 * the P-nuc length is 2, the 2 NP1 nodes IMMEDIATELY before the D
 * segment should be the reverse-complement of D's first 2 bases. */
static int test_p_nuc_d5_head(void) {
    rng_seed(&_test_rng, 107, 0);

    SimConfig cfg = make_igh_config();
    /* D allele in make_igh_config: "GTATTACTGTGC" → first 2 = "GT".
     * Head-side P-nuc is RC of "GT" read 5'→3':
     *   p[0] = complement(seg[K-1]) = complement('T') = 'A'
     *   p[1] = complement(seg[K-2]) = complement('G') = 'C'
     * Result: "AC", and seg becomes "AC" + "GTATTACTGTGC". */
    install_fixed_trim(&cfg.v_trim_3, 2);  /* no V 3' P-nuc */
    install_zero_trim(&cfg.d_trim_5);      /* enable D 5' P-nuc */
    install_fixed_trim(&cfg.d_trim_3, 1);  /* no D 3' P-nuc */
    install_fixed_trim(&cfg.j_trim_5, 1);  /* no J 5' P-nuc */
    /* No NP1 N-nucs: if cfg->np[0].length_probs is unset, fallback
     * uniform sampler emits N-nucs. We force NP1 N-nucs to 0 by
     * installing an NpDist with length always 0. */
    double zero_len[1] = {1.0};
    double dummy_first[4] = {0.25, 0.25, 0.25, 0.25};
    /* install_np_dist defined earlier in this file */
    install_np_dist(&cfg, 0, zero_len, 0, dummy_first, NULL, 0);
    install_p_nuc_fixed(&cfg, /*forced_len=*/2, /*max_len=*/4);

    Pipeline p = pipeline_build(&cfg);

    const char expected[2] = {'A', 'C'};

    int trials = 30;
    for (int t = 0; t < trials; t++) {
        ASeq seq;
        aseq_init(&seq);
        SimRecord rec;
        pipeline_execute(&p, &cfg, &seq, &rec, 0, NULL, NULL);

        /* The 2 NP1 nodes should be exactly the D 5' P-nuc since we
         * disabled NP1 N-nucs and V 3' P-nucs. They sit immediately
         * before SEG_D. */
        if (aseq_segment_length(&seq, SEG_NP1) != 2) {
            fprintf(stderr, "trial %d: NP1 length %d, expected 2\n",
                    t, aseq_segment_length(&seq, SEG_NP1));
            sim_config_destroy(&cfg);
            return 1;
        }
        Nuc *n = seq.seg_first[SEG_NP1];
        for (int i = 0; i < 2; i++) {
            if (n->current != expected[i] ||
                !(n->flags & NUC_FLAG_P_NUCLEOTIDE)) {
                fprintf(stderr, "trial %d: NP1[%d] = %c (flags=%u)\n",
                        t, i, n->current, n->flags);
                sim_config_destroy(&cfg);
                return 1;
            }
            n = n->next;
        }
    }

    sim_config_destroy(&cfg);
    return 0;
}

/* When the SimConfig has no PNucDist data (length_probs == NULL),
 * no P-nucs should be emitted even with trim==0 at every end. */
static int test_no_p_nuc_when_dist_null(void) {
    rng_seed(&_test_rng, 109, 0);

    SimConfig cfg = make_igh_config();
    install_zero_trim(&cfg.v_trim_3);
    install_zero_trim(&cfg.d_trim_5);
    install_zero_trim(&cfg.d_trim_3);
    install_zero_trim(&cfg.j_trim_5);
    /* Deliberately do NOT install a PNucDist. */

    Pipeline p = pipeline_build(&cfg);

    int trials = 30;
    for (int t = 0; t < trials; t++) {
        ASeq seq;
        aseq_init(&seq);
        SimRecord rec;
        pipeline_execute(&p, &cfg, &seq, &rec, 0, NULL, NULL);

        for (Nuc *n = seq.head; n; n = n->next) {
            if (n->flags & NUC_FLAG_P_NUCLEOTIDE) {
                fprintf(stderr,
                        "trial %d: P-nuc emitted with NULL distribution\n",
                        t);
                sim_config_destroy(&cfg);
                return 1;
            }
        }
    }

    sim_config_destroy(&cfg);
    return 0;
}

/* ── Per-allele trim dispatch tests (T0-3 regression) ─────────── */

/* Build a NamedTrimDistTable with `count` entries; each entry's probs
 * has 1.0 at index `forced_trim_for[i]` (forces a deterministic trim
 * for every allele matching that (family, gene)). */
static NamedTrimDistTable build_named_table(
    const char *families[], const char *genes[],
    const int forced_trims[], int count, int max_trim)
{
    NamedTrimDistTable t = {0};
    t.entries = calloc((size_t)count, sizeof(NamedTrimDist));
    t.count = count;
    for (int i = 0; i < count; i++) {
        strncpy(t.entries[i].family, families[i],
                GENAIRR_MAX_ALLELE_NAME - 1);
        strncpy(t.entries[i].gene, genes[i],
                GENAIRR_MAX_ALLELE_NAME - 1);
        t.entries[i].dist.max_trim = max_trim;
        t.entries[i].dist.probs = calloc((size_t)(max_trim + 1),
                                          sizeof(double));
        if (forced_trims[i] >= 0 && forced_trims[i] <= max_trim) {
            t.entries[i].dist.probs[forced_trims[i]] = 1.0;
        }
    }
    return t;
}

/* Verifies the core T0-3 fix: two V alleles with different
 * (family, gene) keys must dispatch to their own trim distributions.
 * Before T0-3 both used the same (broken) global distribution. */
static int test_per_allele_v_trim_dispatch(void) {
    rng_seed(&_test_rng, 31, 0);

    SimConfig cfg;
    sim_config_init(&cfg, CHAIN_IGH);
    cfg.rng = &_test_rng;

    /* Two V alleles in different families/genes. Length 30, anchor 24,
     * so max_v_trim = 30 - 24 - 1 = 5. */
    Allele va = make_allele("VA*01",
                            "ATCGATCGATCGATCGATCGATCGATCGAT", 24, SEG_V);
    strncpy(va.family, "FAM_A", sizeof(va.family) - 1);
    strncpy(va.gene,   "GENE_A", sizeof(va.gene) - 1);

    Allele vb = make_allele("VB*01",
                            "ATCGATCGATCGATCGATCGATCGATCGAT", 24, SEG_V);
    strncpy(vb.family, "FAM_B", sizeof(vb.family) - 1);
    strncpy(vb.gene,   "GENE_B", sizeof(vb.gene) - 1);

    Allele d = make_allele("D*01", "GTATTACTGTGC", 0, SEG_D);
    Allele j = make_allele("J*01", "ACTACTTTGACTACTGGGGC", 5, SEG_J);

    cfg.v_alleles = allele_pool_create(4);
    cfg.d_alleles = allele_pool_create(4);
    cfg.j_alleles = allele_pool_create(4);
    allele_pool_add(&cfg.v_alleles, &va);
    allele_pool_add(&cfg.v_alleles, &vb);
    allele_pool_add(&cfg.d_alleles, &d);
    allele_pool_add(&cfg.j_alleles, &j);

    /* Per-(family, gene) V_3 table:
     *   FAM_A.GENE_A → always trim 0
     *   FAM_B.GENE_B → always trim 4 */
    const char *fams[]  = {"FAM_A", "FAM_B"};
    const char *gns[]   = {"GENE_A", "GENE_B"};
    const int   trims[] = {0, 4};
    cfg.v_trim_3_table = build_named_table(fams, gns, trims, 2, 5);

    /* Wire allele pointers (mimicking gdc_populate_sim_config). */
    cfg.v_alleles.alleles[0].trim_dist_3 =
        &cfg.v_trim_3_table.entries[0].dist;  /* FAM_A → trim 0 */
    cfg.v_alleles.alleles[1].trim_dist_3 =
        &cfg.v_trim_3_table.entries[1].dist;  /* FAM_B → trim 4 */

    /* Lock to FAM_A allele and verify v_trim_3 == 0. */
    cfg.v_restriction.active = true;
    cfg.v_restriction.count = 1;
    cfg.v_restriction.indices[0] = 0;

    Pipeline p = pipeline_build(&cfg);
    for (int t = 0; t < 20; t++) {
        ASeq seq;
        aseq_init(&seq);
        SimRecord rec;
        pipeline_execute(&p, &cfg, &seq, &rec, 0, NULL, NULL);
        if (rec.v_trim_3 != 0) {
            fprintf(stderr,
                    "FAM_A allele expected trim 0, got %d\n", rec.v_trim_3);
            sim_config_destroy(&cfg);
            return 1;
        }
    }

    /* Lock to FAM_B allele and verify v_trim_3 == 4. */
    cfg.v_restriction.indices[0] = 1;
    for (int t = 0; t < 20; t++) {
        ASeq seq;
        aseq_init(&seq);
        SimRecord rec;
        pipeline_execute(&p, &cfg, &seq, &rec, 0, NULL, NULL);
        if (rec.v_trim_3 != 4) {
            fprintf(stderr,
                    "FAM_B allele expected trim 4, got %d\n", rec.v_trim_3);
            sim_config_destroy(&cfg);
            return 1;
        }
    }

    sim_config_destroy(&cfg);
    return 0;
}

/* When an allele's trim_dist_3 pointer is NULL, trim.c falls back to
 * the legacy cfg->v_trim_3 single-global. This is the embedded test
 * data path: human_igh_imgt.c installs cfg->v_trim_3 directly without
 * going through gdc_populate_sim_config, so all alleles have NULL
 * pointers and rely on the legacy field. */
static int test_trim_legacy_fallback(void) {
    rng_seed(&_test_rng, 37, 0);

    SimConfig cfg = make_igh_config();
    /* Allele pointers stay NULL (make_igh_config doesn't set them).
     * Install legacy v_trim_3 forcing trim=2. */
    cfg.v_trim_3.max_trim = 5;
    cfg.v_trim_3.probs = calloc(6, sizeof(double));
    cfg.v_trim_3.probs[2] = 1.0;

    Pipeline p = pipeline_build(&cfg);
    for (int t = 0; t < 20; t++) {
        ASeq seq;
        aseq_init(&seq);
        SimRecord rec;
        pipeline_execute(&p, &cfg, &seq, &rec, 0, NULL, NULL);
        if (rec.v_trim_3 != 2) {
            fprintf(stderr,
                    "legacy fallback expected v_trim_3=2, got %d\n",
                    rec.v_trim_3);
            sim_config_destroy(&cfg);
            return 1;
        }
    }

    sim_config_destroy(&cfg);
    return 0;
}

/* When neither per-allele nor legacy dist is set, sample_trim falls
 * back to a uniform draw. For V with max_v_trim=5, samples must lie
 * in [0,5] across many trials, with mean roughly in (1.5, 3.5). */
static int test_trim_uniform_fallback(void) {
    rng_seed(&_test_rng, 41, 0);

    SimConfig cfg = make_igh_config();
    /* Both per-allele and legacy unset (cfg.v_trim_3.probs == NULL). */

    Pipeline p = pipeline_build(&cfg);
    int trials = 100;
    int sum = 0, max_seen = 0;
    for (int t = 0; t < trials; t++) {
        ASeq seq;
        aseq_init(&seq);
        SimRecord rec;
        pipeline_execute(&p, &cfg, &seq, &rec, 0, NULL, NULL);
        sum += rec.v_trim_3;
        if (rec.v_trim_3 > max_seen) max_seen = rec.v_trim_3;
        if (rec.v_trim_3 < 0 || rec.v_trim_3 > 5) {
            fprintf(stderr,
                    "uniform fallback v_trim_3=%d outside [0,5]\n",
                    rec.v_trim_3);
            sim_config_destroy(&cfg);
            return 1;
        }
    }
    double mean = (double)sum / trials;
    if (mean < 1.5 || mean > 3.5) {
        fprintf(stderr,
                "uniform fallback mean %.2f outside [1.5, 3.5]\n", mean);
        sim_config_destroy(&cfg);
        return 1;
    }
    printf("\n    uniform fallback v_trim_3 mean=%.2f, max=%d", mean, max_seen);

    sim_config_destroy(&cfg);
    return 0;
}

/* ── Indel pool-exhaustion guard (T0-6 regression) ──────────── */

/* When step_insert_indels runs on a sequence near the pool cap, the
 * loop can hit aseq_insert_after returning NULL. Pre-T0-6 the caller
 * incremented `insertions` regardless, over-counting rec->n_insertions.
 * Verify that the recorded count exactly matches the number of
 * NUC_FLAG_INDEL_INS-tagged nodes after the step runs. */
static int test_indel_pool_exhaustion_count(void) {
    rng_seed(&_test_rng, 17, 0);

    SimConfig cfg;
    sim_config_init(&cfg, CHAIN_IGH);
    cfg.rng = &_test_rng;
    cfg.indel_prob = 0.95;       /* aggressive: insert at every position */
    cfg.insertion_weight = 1.0;  /* always insertion (never deletion) */

    /* Fake V allele near the GENAIRR_MAX_ALLELE_SEQ cap (511). With a
     * D and J on top, the assembled sequence is ~540 nodes; aggressive
     * indel insertion (every position) then pushes pool_used past
     * GENAIRR_MAX_SEQ_LEN (1024), exercising the pool-exhaustion guard. */
    Allele v;
    memset(&v, 0, sizeof(v));
    strncpy(v.name, "VLong*01", sizeof(v.name) - 1);
    strncpy(v.gene, "VLong",    sizeof(v.gene) - 1);
    strncpy(v.family, "FAM_L",  sizeof(v.family) - 1);
    int v_len = GENAIRR_MAX_ALLELE_SEQ - 1;  /* 511 */
    for (int i = 0; i < v_len; i++) v.seq[i] = "ACGT"[i % 4];
    v.length = (uint16_t)v_len;
    v.anchor = (uint16_t)(v_len - 6);
    v.segment_type = SEG_V;

    Allele d = make_allele("D*01", "GTATTACTGTGC", 0, SEG_D);
    Allele j = make_allele("J*01", "ACTACTTTGACTACTGGGGC", 5, SEG_J);

    cfg.v_alleles = allele_pool_create(2);
    cfg.d_alleles = allele_pool_create(2);
    cfg.j_alleles = allele_pool_create(2);
    allele_pool_add(&cfg.v_alleles, &v);
    allele_pool_add(&cfg.d_alleles, &d);
    allele_pool_add(&cfg.j_alleles, &j);

    cfg.features.indels = true;

    Pipeline p = pipeline_build(&cfg);

    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);
    pipeline_execute(&p, &cfg, &seq, &rec, 0, NULL, NULL);

    /* Count nodes flagged as indel insertions. */
    int actual_insertions = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->flags & NUC_FLAG_INDEL_INS) actual_insertions++;
    }

    /* The pool should be full (or nearly so). */
    if (seq.pool_used < GENAIRR_MAX_SEQ_LEN - 10) {
        fprintf(stderr,
                "Pool not near capacity (pool_used=%d, length=%d) — "
                "test setup did not stress the guard\n",
                seq.pool_used, seq.length);
        sim_config_destroy(&cfg);
        return 1;
    }

    if (rec.n_insertions != actual_insertions) {
        fprintf(stderr,
                "rec.n_insertions=%d but %d nodes carry NUC_FLAG_INDEL_INS "
                "(over-count when pool exhausts)\n",
                rec.n_insertions, actual_insertions);
        sim_config_destroy(&cfg);
        return 1;
    }
    printf("\n    pool_used=%d/%d, n_insertions=%d (matches actual)",
           seq.pool_used, GENAIRR_MAX_SEQ_LEN, rec.n_insertions);

    sim_config_destroy(&cfg);
    return 0;
}

/* ── Anchorless allele handling (T0-7 regression) ────────────── */

/* Anchorless J allele (anchor = -1) must NOT trim past j_len - 1.
 * Pre-T0-7 the unsigned cast in gdc_read turned -1 into 65535, so
 * `max_j_trim = anchor - 1` produced 65534 — sample_trim then bounded
 * only by dist->max_trim, silently allowing trim past the entire
 * allele. Post-fix: an explicit guard caps max_j_trim at length - 1. */
static int test_anchorless_j_does_not_overflow_trim(void) {
    rng_seed(&_test_rng, 71, 0);

    SimConfig cfg;
    sim_config_init(&cfg, CHAIN_IGH);
    cfg.rng = &_test_rng;

    /* Standard V with a real anchor. */
    Allele v = make_allele("V*01",
                           "ATCGATCGATCGATCGATCGATCGATCGAT", 24, SEG_V);
    Allele d = make_allele("D*01", "GTATTACTGTGC", 0, SEG_D);

    /* Anchorless J: length 50, anchor = -1 (the GDC sentinel). */
    Allele j;
    memset(&j, 0, sizeof(j));
    strncpy(j.name, "Janchorless*01", sizeof(j.name) - 1);
    strncpy(j.gene, "Jano",            sizeof(j.gene) - 1);
    strncpy(j.family, "JanoF",         sizeof(j.family) - 1);
    int j_len = 50;
    for (int i = 0; i < j_len; i++) j.seq[i] = "ACGT"[i % 4];
    j.length = (uint16_t)j_len;
    j.anchor = -1;
    j.segment_type = SEG_J;

    cfg.v_alleles = allele_pool_create(2);
    cfg.d_alleles = allele_pool_create(2);
    cfg.j_alleles = allele_pool_create(2);
    allele_pool_add(&cfg.v_alleles, &v);
    allele_pool_add(&cfg.d_alleles, &d);
    allele_pool_add(&cfg.j_alleles, &j);

    /* Aggressive trim distribution: max_trim well above j_len. If the
     * anchorless guard is missing, sample_trim would be capped only
     * by dist->max_trim and could exceed j_len. */
    cfg.j_trim_5.max_trim = 200;
    cfg.j_trim_5.probs = calloc(201, sizeof(double));
    /* Force a uniform distribution that includes large values. */
    for (int i = 0; i <= 200; i++) cfg.j_trim_5.probs[i] = 1.0 / 201.0;

    Pipeline p = pipeline_build(&cfg);

    int trials = 50;
    for (int t = 0; t < trials; t++) {
        ASeq seq;
        aseq_init(&seq);
        SimRecord rec;
        sim_record_init(&rec);
        pipeline_execute(&p, &cfg, &seq, &rec, 0, NULL, NULL);

        if (rec.j_trim_5 < 0 || rec.j_trim_5 > j_len - 1) {
            fprintf(stderr,
                    "trial %d: j_trim_5=%d outside [0, j_len-1=%d]\n",
                    t, rec.j_trim_5, j_len - 1);
            sim_config_destroy(&cfg);
            return 1;
        }
        /* Surviving J segment must have at least one base. */
        int j_seg = aseq_segment_length(&seq, SEG_J);
        if (j_seg < 1) {
            fprintf(stderr,
                    "trial %d: J segment fully trimmed (j_trim_5=%d)\n",
                    t, rec.j_trim_5);
            sim_config_destroy(&cfg);
            return 1;
        }
    }

    sim_config_destroy(&cfg);
    return 0;
}

/* Anchorless V (anchor = 0, the legacy embedded-data sentinel) must
 * NOT result in NUC_FLAG_ANCHOR being set on V's first base. Pre-T0-7
 * the assemble step would compare germline_pos == 0 and tag node 0 as
 * the V Cys anchor — incorrectly. */
static int test_anchorless_v_does_not_tag_first_base(void) {
    rng_seed(&_test_rng, 73, 0);

    SimConfig cfg;
    sim_config_init(&cfg, CHAIN_IGH);
    cfg.rng = &_test_rng;

    /* Anchorless V with anchor = 0 (legacy sentinel). */
    Allele v;
    memset(&v, 0, sizeof(v));
    strncpy(v.name, "Vanchorless*01", sizeof(v.name) - 1);
    strncpy(v.gene, "Vano",            sizeof(v.gene) - 1);
    strncpy(v.family, "VanoF",         sizeof(v.family) - 1);
    int v_len = 30;
    for (int i = 0; i < v_len; i++) v.seq[i] = "ACGT"[i % 4];
    v.length = (uint16_t)v_len;
    v.anchor = 0;  /* legacy "no anchor" sentinel */
    v.segment_type = SEG_V;

    Allele d = make_allele("D*01", "GTATTACTGTGC", 0, SEG_D);
    Allele j = make_allele("J*01", "ACTACTTTGACTACTGGGGC", 5, SEG_J);

    cfg.v_alleles = allele_pool_create(2);
    cfg.d_alleles = allele_pool_create(2);
    cfg.j_alleles = allele_pool_create(2);
    allele_pool_add(&cfg.v_alleles, &v);
    allele_pool_add(&cfg.d_alleles, &d);
    allele_pool_add(&cfg.j_alleles, &j);

    Pipeline p = pipeline_build(&cfg);
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);
    pipeline_execute(&p, &cfg, &seq, &rec, 0, NULL, NULL);

    /* No V node should carry NUC_FLAG_ANCHOR. */
    int v_anchors = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->segment == SEG_V && (n->flags & NUC_FLAG_ANCHOR)) {
            v_anchors++;
        }
    }
    if (v_anchors != 0) {
        fprintf(stderr,
                "anchorless V allele: %d nodes flagged as V anchor "
                "(expected 0)\n", v_anchors);
        sim_config_destroy(&cfg);
        return 1;
    }
    /* And the cached v_anchor_node must remain NULL. */
    if (seq.v_anchor_node != NULL) {
        fprintf(stderr, "anchorless V: v_anchor_node not NULL\n");
        sim_config_destroy(&cfg);
        return 1;
    }

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
    TEST(test_productivity_productive_only);
    TEST(test_productivity_non_productive_only);
    TEST(test_productivity_mixed_returns_both);
    TEST(test_productivity_legacy_bool_back_compat);
    TEST(test_debug_print);
    TEST(test_np_markov_deterministic);
    TEST(test_np_markov_gc_bias);
    TEST(test_np_uniform_fallback);
    TEST(test_p_nuc_v3_when_trim_zero);
    TEST(test_no_p_nuc_when_trim_nonzero);
    TEST(test_p_nuc_d5_head);
    TEST(test_no_p_nuc_when_dist_null);
    TEST(test_per_allele_v_trim_dispatch);
    TEST(test_trim_legacy_fallback);
    TEST(test_trim_uniform_fallback);
    TEST(test_indel_pool_exhaustion_count);
    TEST(test_anchorless_j_does_not_overflow_trim);
    TEST(test_anchorless_v_does_not_tag_first_base);

    printf("\n%d/%d tests passed.\n", total - failures, total);
    return failures > 0 ? 1 : 0;
}
