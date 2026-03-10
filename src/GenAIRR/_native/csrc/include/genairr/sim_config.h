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

/* ── Trim distribution (per gene family) ──────────────────────── */

typedef struct {
    double  *probs;     /* probability for each trim amount */
    int      max_trim;  /* length of probs array            */
} TrimDist;

/* ── Feature flags ────────────────────────────────────────────── */

typedef struct {
    bool  productive;
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

    /* Trim distributions — indexed by family or gene.
     * For simplicity, we start with a single distribution per
     * segment; family-specific distributions can be added later. */
    TrimDist      v_trim_3;
    TrimDist      d_trim_5;
    TrimDist      d_trim_3;
    TrimDist      j_trim_5;

    /* NP region parameters */
    int           np1_length_mean;
    int           np1_length_max;
    int           np2_length_mean;
    int           np2_length_max;

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
} SimConfig;

/* ── Lifecycle ────────────────────────────────────────────────── */

/** Initialize a SimConfig with sensible defaults. */
void  sim_config_init(SimConfig *cfg, ChainType chain_type);

/** Free all dynamically allocated data within a SimConfig. */
void  sim_config_destroy(SimConfig *cfg);

#endif /* GENAIRR_SIM_CONFIG_H */
