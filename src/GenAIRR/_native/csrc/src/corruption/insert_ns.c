/**
 * insert_ns.c — Replace random bases with ambiguous N.
 *
 * Each position has a small probability of being replaced with 'N'.
 * Only operates on germline segment positions.
 */

#include "genairr/pipeline.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>

void step_insert_ns(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)rec;

    double prob = cfg->n_prob;
    int count = 0;

    for (Nuc *n = seq->head; n; n = n->next) {
        if (n->current == 'N') continue;
        if (rand_uniform() >= prob) continue;

        aseq_mutate(seq, n, 'N', NUC_FLAG_IS_N);
        count++;
    }

    TRACE("[insert_ns] replaced %d bases with N (prob=%.4f)", count, prob);
}
