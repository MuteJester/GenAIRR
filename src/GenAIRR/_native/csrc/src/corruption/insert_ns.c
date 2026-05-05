/**
 * insert_ns.c — Replace random bases with ambiguous N.
 *
 * Each position has a small probability of being replaced with 'N'.
 * Only operates on germline segment positions.
 *
 * T2-12 / V5 step 5: when the user explicitly requested
 * PRODUCTIVITY_PRODUCTIVE_ONLY, candidate N substitutions are routed
 * through the shared productivity guard at OBSERVED scope. In
 * practice this still protects the V/J anchor codons, but the
 * decision now comes from a cloned productivity reevaluation rather
 * than a step-local anchor special case.
 *
 * Background: AIRR `productive` is recomputed AFTER corruption
 * (airr.c::derive_final_productivity). If an N lands on the conserved
 * Cys/Trp/Phe codon, the codon retranslates to '?', the
 * conserved-residue rule fails, and `productive` silently flips to
 * False — even though the underlying rearrangement was productive.
 *
 * Other productivity modes (MIXED, NON_PRODUCTIVE_ONLY) keep the
 * full noise model — N can land anywhere, including anchors. Real
 * AIRR-seq pipelines do see Ns at conserved residues; users who
 * don't have a productive contract get realistic noise behavior.
 */

#include "genairr/pipeline.h"
#include "genairr/productivity_guard.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>

void step_insert_ns(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    double prob = cfg->n_prob;
    int count = 0;

    for (Nuc *n = seq->head; n; n = n->next) {
        if (n->current == 'N') continue;
        if (rng_uniform(cfg->rng) >= prob) continue;

        ProductivityDecision decision = productivity_guard_substitution(
            cfg, seq, rec, PROD_STAGE_OBSERVED, n, 'N');
        if (decision != PROD_DECISION_ALLOW) continue;

        aseq_mutate(seq, n, 'N', NUC_FLAG_IS_N);
        count++;
    }

    TRACE("[insert_ns] replaced %d bases with N (prob=%.4f)", count, prob);
}
