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
#include "genairr/productivity_guard.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>

void step_uniform_mutate(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    int len = aseq_length(seq);
    if (len == 0) return;

    const bool productive_guard = simcfg_productive_only(cfg);

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
                  rng_uniform(cfg->rng) * (cfg->max_mutation_rate - cfg->min_mutation_rate);

    int target_mutations = (int)(rate * n_mutable + 0.5);

    TRACE("[uniform] rate_range=[%.4f, %.4f], sampled_rate=%.4f, mutable=%d/%d, target=%d, productive_guard=%s",
          cfg->min_mutation_rate, cfg->max_mutation_rate, rate, n_mutable, len,
          target_mutations, productive_guard ? "yes" : "no");

    if (target_mutations <= 0) return;

    /* Randomly select positions to mutate */
    static const char bases[] = "acgt";
    int applied = 0;

    for (Nuc *n = seq->head; n && applied < target_mutations; n = n->next) {
        if (!nuc_is_germline_segment(n) || n->segment == SEG_C) continue;

        /* Probability = remaining_mutations / remaining_positions */
        int remaining = n_mutable - applied;
        if (remaining <= 0) break;

        double p = (double)(target_mutations - applied) / remaining;
        if (rng_uniform(cfg->rng) >= p) {
            /* Count this as a passed position even if not mutated */
            continue;
        }

        /* Mutate to a different base */
        char new_base;
        if (productive_guard) {
            char allowed[3];
            int allowed_count = 0;

            for (int i = 0; i < 4; i++) {
                char candidate = bases[i];
                if (candidate == n->current) continue;

                ProductivityDecision decision = productivity_guard_substitution(
                    cfg, seq, rec, PROD_STAGE_MOLECULE, n, candidate);
                if (decision == PROD_DECISION_ALLOW) {
                    allowed[allowed_count++] = candidate;
                }
            }

            if (allowed_count == 0) {
                continue;
            }

            new_base = allowed[rng_range(cfg->rng, (uint32_t)allowed_count)];
        } else {
            do {
                new_base = bases[rng_range(cfg->rng, 4)];
            } while (new_base == n->current);
        }

        aseq_mutate(seq, n, new_base, NUC_FLAG_MUTATED);
        applied++;
    }

    TRACE("[uniform] done: %d mutations applied (%.4f effective rate)",
          applied, (double)applied / n_mutable);
}
