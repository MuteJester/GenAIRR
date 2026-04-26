/**
 * insert_indels.c — Random insertions and deletions.
 *
 * Applies random indels at valid positions in the sequence.
 * Insertions add a random base after the target node.
 * Deletions remove the target node.
 * Only operates on germline segment positions (V/D/J/C).
 */

#include "genairr/pipeline.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>

void step_insert_indels(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    int len = aseq_length(seq);
    if (len < 20) return;

    double prob = cfg->indel_prob;
    double ins_weight = cfg->insertion_weight;

    TRACE("[indels] seq_len=%d, indel_prob=%.4f, insertion_weight=%.2f", len, prob, ins_weight);

    int insertions = 0, deletions = 0;
    for (Nuc *n = seq->head; n; ) {
        Nuc *next = n->next;  /* save before potential deletion */

        if (!nuc_is_germline_segment(n)) {
            n = next;
            continue;
        }
        if (nuc_is_anchor(n)) {
            n = next;
            continue;
        }

        if (rng_uniform(cfg->rng) >= prob) {
            n = next;
            continue;
        }

        if (rng_uniform(cfg->rng) < ins_weight) {
            /* Insertion: add a random base after this node. Pool may
             * be exhausted on very long sequences (>= GENAIRR_MAX_SEQ_LEN
             * nodes already allocated); aseq_insert_after returns NULL
             * in that case. Break out — further insertions will also
             * fail and we'd over-count rec->n_insertions. */
            static const char bases[] = "ACGT";
            char base = bases[rng_range(cfg->rng, 4)];
            if (aseq_insert_after(seq, n, base, n->segment,
                                  NUC_FLAG_INDEL_INS) == NULL) {
                TRACE("[indels] pool exhausted at insertion; stopping "
                      "(pool_used=%d/%d)",
                      seq->pool_used, GENAIRR_MAX_SEQ_LEN);
                break;
            }
            insertions++;
        } else {
            /* Deletion: remove this node */
            aseq_delete(seq, n);
            deletions++;
        }

        n = next;
    }

    rec->n_insertions = insertions;
    rec->n_deletions  = deletions;

    TRACE("[indels] applied %d insertions, %d deletions", insertions, deletions);
}
