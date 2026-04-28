/**
 * insert_ns.c — Replace random bases with ambiguous N.
 *
 * Each position has a small probability of being replaced with 'N'.
 * Only operates on germline segment positions.
 *
 * T2-12: when the user explicitly requested PRODUCTIVITY_PRODUCTIVE_ONLY,
 * the V and J anchor codons are protected from N corruption. This
 * preserves the user's contract that every output is productive.
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
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdbool.h>
#include <stdlib.h>

void step_insert_ns(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)rec;

    double prob = cfg->n_prob;
    int count = 0;

    /* Only the FIRST nt of each anchor codon carries NUC_FLAG_ANCHOR
     * (see aseq.c:107). To protect the full 3-nt codon we use a small
     * residual counter that ticks down for the next 2 nodes after a
     * flagged position is hit. */
    const bool protect_anchors =
        (cfg->features.productivity == PRODUCTIVITY_PRODUCTIVE_ONLY);
    int anchor_codon_remaining = 0;

    for (Nuc *n = seq->head; n; n = n->next) {
        /* Compute anchor-protection state FIRST so the residual
         * counter advances on every node, regardless of whether
         * we end up drawing rng for this node. Keeps the codon
         * boundary tracking correct even when interspersed with
         * already-N nodes. */
        bool in_anchor_codon = false;
        if (protect_anchors) {
            if (n->flags & NUC_FLAG_ANCHOR) {
                in_anchor_codon = true;
                anchor_codon_remaining = 2;       /* next 2 nts also part of this codon */
            } else if (anchor_codon_remaining > 0) {
                in_anchor_codon = true;
                anchor_codon_remaining--;
            }
        }

        if (n->current == 'N') continue;
        if (rng_uniform(cfg->rng) >= prob) continue;

        /* T2-12: under PRODUCTIVITY_PRODUCTIVE_ONLY, suppress the
         * mutation when the position falls inside an anchor codon.
         * We've already drawn rng_uniform above, so this preserves
         * RNG stream stability — only the application of the
         * mutation is skipped. */
        if (in_anchor_codon) continue;

        aseq_mutate(seq, n, 'N', NUC_FLAG_IS_N);
        count++;
    }

    TRACE("[insert_ns] replaced %d bases with N (prob=%.4f)", count, prob);
}
