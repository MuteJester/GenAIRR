/**
 * primer_mask.c — FR1 overwritten with germline primer sequence.
 *
 * Replaces the FR1 region (first ~78 germline positions of V) with
 * the original germline V allele sequence. This removes any mutations
 * in the masked region.
 *
 * Uses IMGT numbering: FR1 = germline positions 0..77 (first 78 bases).
 */

#include "genairr/pipeline.h"
#include "genairr/trace.h"
#include <string.h>

#define IMGT_FR1_END  78  /* IMGT gapped FR1 boundary */

void step_primer_mask(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)cfg;
    if (!rec->v_allele) return;

    /* Determine mask length: either configured or full FR1 */
    int mask_len = cfg->primer_mask_length;
    if (mask_len <= 0) mask_len = IMGT_FR1_END;

    /* Walk V segment and restore germline bases for positions < mask_len */
    Nuc *n = seq->seg_first[SEG_V];
    int masked = 0;

    while (n && n->segment == SEG_V && n->germline_pos < (uint16_t)mask_len) {
        if (n->germline != '\0' && n->current != n->germline) {
            aseq_revert(seq, n);
            n->flags &= ~(NUC_FLAG_SEQ_ERROR | NUC_FLAG_PCR_ERROR);
            masked++;
        }
        n = n->next;
    }

    rec->primer_masked_length = masked;
    TRACE("[primer_mask] reverted %d mutations in FR1 (mask_len=%d)", masked, mask_len);
}
