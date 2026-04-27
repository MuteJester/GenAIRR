/**
 * sim_config.c — SimConfig initialization and teardown.
 */

#include "genairr/sim_config.h"
#include <stdlib.h>
#include <string.h>

void sim_config_init(SimConfig *cfg, ChainType chain_type) {
    memset(cfg, 0, sizeof(*cfg));
    cfg->chain_type = chain_type;

    /* Sensible defaults */
    cfg->max_productive_attempts = 25;
    cfg->min_mutation_rate = 0.01;
    cfg->max_mutation_rate = 0.10;
    /* NP regions: leave length_probs/transitions as NULL by default.
     * The assemble step falls back to a uniform sampler in that case;
     * GDC-loaded SimConfigs overwrite these with the empirical Markov
     * tables in gdc_populate_sim_config. */
    cfg->n_np_regions = 0;

    /* 5'/3' corruption: uniform range (base counts) */
    cfg->corrupt_5_remove_min = 1;
    cfg->corrupt_5_remove_max = 20;
    cfg->corrupt_5_add_min    = 1;
    cfg->corrupt_5_add_max    = 10;
    cfg->corrupt_3_remove_min = 1;
    cfg->corrupt_3_remove_max = 20;
    cfg->corrupt_3_add_min    = 1;
    cfg->corrupt_3_add_max    = 10;

    /* Quality errors (Illumina linear profile) */
    cfg->base_error_rate   = 0.001;   /* Q30 */
    cfg->peak_error_rate   = 0.02;    /* Q17 */
    cfg->transition_weight = 0.7;

    /* PCR amplification */
    cfg->pcr_error_rate = 1e-4;       /* Taq default */
    cfg->pcr_n_cycles   = 30;

    /* UMI */
    cfg->umi_length = 12;

    /* Paired-end */
    cfg->pe_read_length = 300;

    /* Contaminants */
    cfg->contamination_prob = 0.01;
    cfg->contaminant_type   = 0;      /* random */

    /* Indels */
    cfg->indel_prob       = 0.01;
    cfg->insertion_weight = 0.5;

    /* Insert Ns */
    cfg->n_prob = 0.02;

    /* Trim to length */
    cfg->max_sequence_length = 0;     /* 0 = no trimming */

    /* Selection pressure. The anchor (V Cys, J W/F) codon is a conserved
     * structural residue — disrupting it breaks the V/J domain fold, so
     * R-mutations in the anchor codon are nearly always reverted. */
    cfg->selection_strength   = 0.5;
    cfg->cdr_r_acceptance     = 0.85;
    cfg->fwr_r_acceptance     = 0.40;
    cfg->anchor_r_acceptance  = 0.0;

    /* D inversion */
    cfg->d_inversion_prob = 0.02;

    /* Receptor revision */
    cfg->revision_prob = 0.05;
    cfg->footprint_min = 5;
    cfg->footprint_max = 20;

    /* Long-read errors */
    cfg->long_read_error_rate = 0.03;
    cfg->min_run_length       = 3;
    cfg->insertion_bias       = 0.6;

    /* Reverse complement */
    cfg->rc_prob = 0.5;

    /* Primer mask */
    cfg->primer_mask_length = 0;      /* 0 = full FR1 */

    /* All features off by default */
    memset(&cfg->features, 0, sizeof(cfg->features));
}

void sim_config_destroy(SimConfig *cfg) {
    allele_pool_destroy(&cfg->v_alleles);
    allele_pool_destroy(&cfg->d_alleles);
    allele_pool_destroy(&cfg->j_alleles);
    allele_pool_destroy(&cfg->c_alleles);

    /* Free trim distributions (legacy single-global) */
    if (cfg->v_trim_3.probs) { free(cfg->v_trim_3.probs); cfg->v_trim_3.probs = NULL; }
    if (cfg->d_trim_5.probs) { free(cfg->d_trim_5.probs); cfg->d_trim_5.probs = NULL; }
    if (cfg->d_trim_3.probs) { free(cfg->d_trim_3.probs); cfg->d_trim_3.probs = NULL; }
    if (cfg->j_trim_5.probs) { free(cfg->j_trim_5.probs); cfg->j_trim_5.probs = NULL; }

    /* Free per-(family, gene) trim tables. Each entry's TrimDist owns
     * its probs[] array. */
    NamedTrimDistTable *tables[4] = {
        &cfg->v_trim_3_table, &cfg->d_trim_5_table,
        &cfg->d_trim_3_table, &cfg->j_trim_5_table,
    };
    for (int t = 0; t < 4; t++) {
        for (int i = 0; i < tables[t]->count; i++) {
            free(tables[t]->entries[i].dist.probs);
        }
        free(tables[t]->entries);
        tables[t]->entries = NULL;
        tables[t]->count = 0;
    }

    /* Free NP distributions */
    for (int i = 0; i < 2; i++) {
        if (cfg->np[i].length_probs) {
            free(cfg->np[i].length_probs);
            cfg->np[i].length_probs = NULL;
        }
        if (cfg->np[i].transitions) {
            free(cfg->np[i].transitions);
            cfg->np[i].transitions = NULL;
        }
    }

    /* Free P-nucleotide length distribution */
    if (cfg->p_nuc_dist.length_probs) {
        free(cfg->p_nuc_dist.length_probs);
        cfg->p_nuc_dist.length_probs = NULL;
    }
}
