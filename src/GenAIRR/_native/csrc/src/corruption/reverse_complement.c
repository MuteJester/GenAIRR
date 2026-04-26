/**
 * reverse_complement.c — Antisense read orientation.
 *
 * Simply delegates to aseq_reverse_complement() which handles the
 * linked-list reversal and base complementation. Flags and germline
 * annotations are preserved by aseq_reverse_complement().
 */

#include "genairr/pipeline.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"

void step_reverse_complement(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    if (rng_uniform(cfg->rng) >= cfg->rc_prob) return;

    TRACE("[reverse_complement] applying reverse-complement (prob=%.2f, seq_len=%d)",
          cfg->rc_prob, aseq_length(seq));
    aseq_reverse_complement(seq);
    rec->is_reverse_complement = true;
}
