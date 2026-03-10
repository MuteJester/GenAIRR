/**
 * uniform.c — Uniform mutation model.
 *
 * Each mutable position (V/D/J segment) has equal probability of
 * mutation. Mutation rate is uniformly sampled between min and max.
 * Substitution is uniform: equal probability for any non-self base.
 *
 * NP-aware: only V/D/J positions are mutable (same as S5F).
 */

#include "genairr/pipeline.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>

void step_uniform_mutate(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)rec;
    int len = aseq_length(seq);
    if (len == 0) return;

    /* Count mutable positions (V/D/J only) */
    int n_mutable = 0;
    for (Nuc *n = seq->head; n; n = n->next) {
        if (nuc_is_germline_segment(n) && n->segment != SEG_C) {
            n_mutable++;
        }
    }
    if (n_mutable == 0) return;

    /* Sample mutation rate */
    double rate = cfg->min_mutation_rate +
                  rand_uniform() * (cfg->max_mutation_rate - cfg->min_mutation_rate);

    int target_mutations = (int)(rate * n_mutable + 0.5);

    TRACE("[uniform] rate_range=[%.4f, %.4f], sampled_rate=%.4f, mutable=%d/%d, target=%d",
          cfg->min_mutation_rate, cfg->max_mutation_rate, rate, n_mutable, len, target_mutations);

    if (target_mutations <= 0) return;

    /* Randomly select positions to mutate */
    static const char bases[] = "ACGT";
    int applied = 0;

    for (Nuc *n = seq->head; n && applied < target_mutations; n = n->next) {
        if (!nuc_is_germline_segment(n) || n->segment == SEG_C) continue;

        /* Probability = remaining_mutations / remaining_positions */
        int remaining = n_mutable - applied;
        if (remaining <= 0) break;

        double p = (double)(target_mutations - applied) / remaining;
        if (rand_uniform() >= p) {
            /* Count this as a passed position even if not mutated */
            continue;
        }

        /* Mutate to a different base */
        char new_base;
        do {
            new_base = bases[rand() % 4];
        } while (new_base == n->current);

        aseq_mutate(seq, n, new_base, NUC_FLAG_MUTATED);
        applied++;
    }

    TRACE("[uniform] done: %d mutations applied (%.4f effective rate)",
          applied, (double)applied / n_mutable);
}
