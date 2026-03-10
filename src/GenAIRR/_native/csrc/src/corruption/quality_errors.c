/**
 * quality_errors.c — Illumina position-dependent sequencing errors.
 *
 * Linear error profile: rate(pos) = base_rate + (pos/len) * (peak - base).
 * Transition bias: ~70% A↔G / C↔T transitions (configurable).
 */

#include "genairr/pipeline.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>

/* Transition partners: A↔G, C↔T */
static char transition_of(char base) {
    switch (base) {
        case 'A': return 'G';
        case 'G': return 'A';
        case 'C': return 'T';
        case 'T': return 'C';
        default:  return base;
    }
}

/* Random transversion (not self, not transition) */
static char transversion_of(char base) {
    switch (base) {
        case 'A': return (rand() % 2) ? 'C' : 'T';
        case 'G': return (rand() % 2) ? 'C' : 'T';
        case 'C': return (rand() % 2) ? 'A' : 'G';
        case 'T': return (rand() % 2) ? 'A' : 'G';
        default:  return base;
    }
}

void step_quality_errors(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)rec;
    int len = aseq_length(seq);
    if (len == 0) return;

    double base_rate = cfg->base_error_rate;
    double peak_rate = cfg->peak_error_rate;
    double tw = cfg->transition_weight;

    TRACE("[quality_errors] seq_len=%d, base_rate=%.4f, peak_rate=%.4f, transition_weight=%.2f",
          len, base_rate, peak_rate, tw);

    int error_count = 0;
    int pos = 0;
    for (Nuc *n = seq->head; n; n = n->next, pos++) {
        if (n->current == 'N') continue;

        double rate = base_rate + ((double)pos / len) * (peak_rate - base_rate);
        if (rand_uniform() >= rate) continue;

        /* Apply error: transition or transversion */
        char new_base;
        if (rand_uniform() < tw) {
            new_base = transition_of(n->current);
        } else {
            new_base = transversion_of(n->current);
        }

        aseq_mutate(seq, n, new_base, NUC_FLAG_SEQ_ERROR);
        error_count++;
    }

    TRACE("[quality_errors] introduced %d sequencing errors", error_count);
}
