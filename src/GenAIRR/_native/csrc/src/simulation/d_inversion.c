/**
 * d_inversion.c — D gene inversion (reverse-complement D segment in-place).
 *
 * With probability d_inversion_prob, reverse-complements only the D
 * segment of the assembled sequence. This is done by walking the D
 * segment, collecting bases, then writing them back in reverse-complement.
 *
 * This must run BEFORE mutation (biological event during rearrangement).
 */

#include "genairr/pipeline.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <string.h>

void step_d_inversion(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    if (!rec->d_allele) return;
    if (rng_uniform(cfg->rng) >= cfg->d_inversion_prob) return;

    Nuc *first = seq->seg_first[SEG_D];
    Nuc *last  = seq->seg_last[SEG_D];
    if (!first || !last) return;

    /* Collect D bases and germline positions */
    int d_len = aseq_segment_length(seq, SEG_D);
    if (d_len <= 1) return;

    char bases[GENAIRR_MAX_ALLELE_SEQ];
    uint16_t gpos[GENAIRR_MAX_ALLELE_SEQ];
    int i = 0;
    for (Nuc *n = first; n; n = n->next) {
        if (n->segment != SEG_D) break;
        bases[i] = n->current;
        gpos[i]  = n->germline_pos;
        i++;
    }

    /* Write back in reverse-complement, reversing germline_pos
     * so that AlleleBitmap lookups remain correct.
     * D inversion is a rearrangement event: the germline reference
     * for this D is now the reverse complement of the allele. */
    int j = d_len - 1;
    for (Nuc *n = first; n; n = n->next) {
        if (n->segment != SEG_D) break;
        n->current     = complement(bases[j]);
        n->germline    = complement(bases[j]);
        n->germline_pos = gpos[j];
        j--;
    }

    rec->d_inverted = true;
    TRACE("[d_inversion] reverse-complemented D segment (%dbp, allele=%s)",
          d_len, rec->d_allele->name);
}
