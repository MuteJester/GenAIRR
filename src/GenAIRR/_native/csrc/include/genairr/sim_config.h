/**
 * sim_config.h — Simulation configuration (built once, used many times).
 *
 * SimConfig holds all reference data (allele pools, trim distributions)
 * and feature flags. It is the C equivalent of Python's DataConfig +
 * convenience function parameters combined.
 */

#ifndef GENAIRR_SIM_CONFIG_H
#define GENAIRR_SIM_CONFIG_H

#include "types.h"
#include "allele.h"
#include "rand_util.h"

/* ── Trim distribution (per gene family) ──────────────────────── */

typedef struct TrimDist {
    double  *probs;     /* probability for each trim amount */
    int      max_trim;  /* length of probs array            */
} TrimDist;

/* ── Per-(family, gene) trim distribution table ───────────────── */
/* Each Allele's trim_dist_5/3 pointer is set at gdc_populate_sim_config
 * time to point at the matching entry's `dist`. Pointers stay valid
 * for the SimConfig's lifetime (the entries[] array is owned). */
typedef struct {
    char     family[GENAIRR_MAX_ALLELE_NAME];
    char     gene[GENAIRR_MAX_ALLELE_NAME];
    TrimDist dist;     /* probs[] is owned */
} NamedTrimDist;

typedef struct {
    NamedTrimDist *entries;   /* owned, length = count */
    int            count;
} NamedTrimDistTable;

/* ── NP region distribution (TdT Markov model) ────────────────── */
/* Position-inhomogeneous order-1 Markov chain over A/C/G/T.
 * Length is sampled from the empirical length distribution; the
 * first base from `first_base`; subsequent bases from
 * `transitions[pos][prev]`. When length_probs is NULL the assemble
 * step falls back to a uniform sampler (legacy / manual SimConfig). */
typedef struct {
    double  *length_probs;     /* dense [0..max_length], owned, may be NULL */
    int      max_length;
    double   first_base[4];    /* P(A), P(C), P(G), P(T) — normalized */
    int      n_positions;      /* rows in `transitions` */
    double  *transitions;      /* flat [n_positions × 16], row-major
                                  [pos × 16 + prev × 4 + next], owned */
} NpDist;

/* ── P-nucleotide length distribution ─────────────────────────── */
/* Palindromic nucleotides (P-nucs) arise from RAG/Artemis hairpin
 * opening before exonuclease trim. They appear at a segment end
 * only when that end has zero exonuclease trim. The base sequence
 * of a K-base P-nuc is the reverse-complement of the K bases at
 * that segment edge.
 *
 * One shared length distribution is sampled independently at each
 * eligible end (V 3', D 5', D 3', J 5'). Per-end distributions are
 * a possible future extension when training data is available. When
 * length_probs is NULL no P-nucs are emitted. */
typedef struct {
    double  *length_probs;     /* dense [0..max_length], owned, may be NULL */
    int      max_length;
} PNucDist;

/* ── Feature flags ────────────────────────────────────────────── */

typedef struct {
    /* Productivity filter mode. Zero-initialized (memset) defaults to
     * PRODUCTIVITY_MIXED (no filtering). The legacy bool field name is
     * gone — use `productivity` everywhere; back-compat is provided at
     * the API boundary (genairr_set_feature("productive", val)). */
    ProductivityMode productivity;
    bool  mutate;
    bool  selection_pressure;
    bool  csr;
    bool  receptor_revision;
    bool  d_inversion;

    bool  corrupt_5_prime;
    bool  corrupt_3_prime;
    bool  quality_errors;
    bool  paired_end;
    bool  pcr;
    bool  umi;
    bool  primer_mask;
    bool  reverse_complement;
    bool  contaminants;
    bool  long_read_errors;
    bool  indels;
    bool  insert_ns;
} SimFeatures;

/* ── Main config ──────────────────────────────────────────────── */

typedef struct {
    ChainType     chain_type;
    SimFeatures   features;

    /* Allele pools (one per segment type that exists) */
    AllelePool    v_alleles;
    AllelePool    d_alleles;
    AllelePool    j_alleles;
    AllelePool    c_alleles;

    /* Trim distributions — legacy single-global per segment, kept as
     * fallback for the embedded test data path (which doesn't go
     * through gdc_populate_sim_config) and for manual SimConfigs in
     * C tests. Production GDC-loaded SimConfigs use the per-allele
     * pointers (see Allele.trim_dist_5/3) which override these. */
    TrimDist      v_trim_3;
    TrimDist      d_trim_5;
    TrimDist      d_trim_3;
    TrimDist      j_trim_5;

    /* Per-(family, gene) trim tables — populated by GDC load.
     * Allele.trim_dist_5/3 pointers index into these tables. */
    NamedTrimDistTable  v_trim_3_table;
    NamedTrimDistTable  d_trim_5_table;
    NamedTrimDistTable  d_trim_3_table;
    NamedTrimDistTable  j_trim_5_table;

    /* NP region distributions ([0]=NP1, [1]=NP2). When n_np_regions
     * is 1 (VJ chains: kappa, lambda, TCRA, TCRG), only np[0] is used. */
    NpDist        np[2];
    int           n_np_regions;

    /* P-nucleotide length distribution (shared across all four
     * eligible segment ends: V 3', D 5', D 3', J 5'). */
    PNucDist      p_nuc_dist;

    /* Mutation parameters */
    double        min_mutation_rate;
    double        max_mutation_rate;

    /* 5'/3' corruption parameters (base counts) */
    int           corrupt_5_remove_min;
    int           corrupt_5_remove_max;
    int           corrupt_5_add_min;
    int           corrupt_5_add_max;
    int           corrupt_3_remove_min;
    int           corrupt_3_remove_max;
    int           corrupt_3_add_min;
    int           corrupt_3_add_max;

    /* Quality error parameters */
    double        base_error_rate;
    double        peak_error_rate;
    double        transition_weight;

    /* PCR parameters */
    double        pcr_error_rate;
    int           pcr_n_cycles;

    /* UMI parameters */
    int           umi_length;

    /* Paired-end parameters */
    int           pe_read_length;

    /* Contaminant parameters */
    double        contamination_prob;
    int           contaminant_type;      /* 0=random, 1=phix */

    /* Indel parameters */
    double        indel_prob;
    double        insertion_weight;      /* fraction that are insertions */

    /* Insert N parameters */
    double        n_prob;

    /* Trim-to-length */
    int           max_sequence_length;

    /* Selection pressure */
    double        selection_strength;
    double        cdr_r_acceptance;
    double        fwr_r_acceptance;
    double        anchor_r_acceptance;   /* V-Cys / J-W or F anchor codon */

    /* D inversion probability */
    double        d_inversion_prob;

    /* Receptor revision */
    double        revision_prob;
    int           footprint_min;
    int           footprint_max;

    /* Long-read / base composition skew */
    double        long_read_error_rate;
    int           min_run_length;
    double        insertion_bias;

    /* Reverse complement probability */
    double        rc_prob;

    /* Primer mask length (0 = full FR1) */
    int           primer_mask_length;

    /* Productive mode */
    int           max_productive_attempts;

    /* Allele restrictions (locking) */
    AlleleRestriction  v_restriction;
    AlleleRestriction  d_restriction;
    AlleleRestriction  j_restriction;
    AlleleRestriction  c_restriction;

    /* Per-simulator RNG. Owned externally (by GenAIRRSimulator); the
     * SimConfig only borrows the pointer. step functions take
     * `const SimConfig *cfg` but `cfg->rng->state` is mutable because
     * what is const is the SimConfig fields, not the pointed-to state.
     * Manual SimConfigs (C tests) must set this before pipeline_execute;
     * sim_config_init initializes it to NULL. */
    RngState     *rng;
} SimConfig;

/* ── Lifecycle ────────────────────────────────────────────────── */

/** Initialize a SimConfig with sensible defaults. */
void  sim_config_init(SimConfig *cfg, ChainType chain_type);

/** Free all dynamically allocated data within a SimConfig. */
void  sim_config_destroy(SimConfig *cfg);

#endif /* GENAIRR_SIM_CONFIG_H */
