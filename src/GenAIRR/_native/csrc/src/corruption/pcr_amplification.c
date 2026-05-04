/**
 * pcr_amplification.c — Uniform PCR polymerase errors.
 *
 * Uniform error rate across ALL positions (including NP regions).
 * Effective rate: 1 - (1 - error_rate)^n_cycles.
 * Uniform substitution: equal probability for any non-self base.
 */

#include "genairr/pipeline.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <ctype.h>
#include <math.h>
#include <stdlib.h>

static char random_other_base(RngState *rng, char base) {
    /* Lowercase to keep case consistent with V/D/J reference and NP
     * regions; the loop below uses a case-insensitive compare so it
     * still rejects "same letter" regardless of caller's case. */
    static const char bases[] = "acgt";
    char base_lo = (char)tolower((unsigned char)base);
    char alt;
    do {
        alt = bases[rng_range(rng, 4)];
    } while (alt == base_lo);
    return alt;
}

void step_pcr_amplification(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    int len = aseq_length(seq);
    if (len == 0) return;

    /* Effective error rate after n_cycles of PCR */
    double eff_rate = 1.0 - pow(1.0 - cfg->pcr_error_rate, cfg->pcr_n_cycles);
    rec->pcr_n_cycles = cfg->pcr_n_cycles;

    TRACE("[pcr] seq_len=%d, per_cycle_rate=%.6f, n_cycles=%d, effective_rate=%.6f",
          len, cfg->pcr_error_rate, cfg->pcr_n_cycles, eff_rate);

    int error_count = 0;
    for (Nuc *n = seq->head; n; n = n->next) {
        if (n->current == 'N') continue;
        if (rng_uniform(cfg->rng) >= eff_rate) continue;

        char new_base = random_other_base(cfg->rng, n->current);
        aseq_mutate(seq, n, new_base, NUC_FLAG_PCR_ERROR);
        error_count++;
    }

    TRACE("[pcr] introduced %d PCR polymerase errors", error_count);
}
