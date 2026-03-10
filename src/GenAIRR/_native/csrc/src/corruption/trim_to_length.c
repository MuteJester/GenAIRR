/**
 * trim_to_length.c — Enforce maximum sequence length from 5' end.
 *
 * If the sequence exceeds max_sequence_length, remove nodes from
 * the 3' end until the length constraint is satisfied.
 */

#include "genairr/pipeline.h"
#include "genairr/trace.h"

void step_trim_to_length(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)rec;
    int max_len = cfg->max_sequence_length;
    if (max_len <= 0) return;

    int before = aseq_length(seq);
    while (aseq_length(seq) > max_len && seq->tail) {
        aseq_delete(seq, seq->tail);
    }
    int removed = before - aseq_length(seq);
    if (removed > 0) {
        TRACE("[trim_to_length] trimmed %d bases from 3' end (max=%d, was=%d)", removed, max_len, before);
    }
}
