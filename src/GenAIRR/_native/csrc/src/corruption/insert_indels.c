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
    (void)rec;
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

        if (rand_uniform() >= prob) {
            n = next;
            continue;
        }

        if (rand_uniform() < ins_weight) {
            /* Insertion: add a random base after this node */
            static const char bases[] = "ACGT";
            char base = bases[rand() % 4];
            aseq_insert_after(seq, n, base, n->segment, NUC_FLAG_INDEL_INS);
            insertions++;
        } else {
            /* Deletion: remove this node */
            aseq_delete(seq, n);
            deletions++;
        }

        n = next;
    }

    TRACE("[indels] applied %d insertions, %d deletions", insertions, deletions);
}
