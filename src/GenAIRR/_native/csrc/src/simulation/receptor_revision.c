/**
 * receptor_revision.c — V-replacement with old-V footprint.
 *
 * Simulates receptor editing: replaces V segment content with a
 * different V allele while preserving the last 5-20nt of the
 * original V at the V-D/V-J junction (the "footprint").
 *
 * Length-preserving: replaces V bases in-place, no position shifts.
 */

#include "genairr/pipeline.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>
#include <string.h>

void step_receptor_revision(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    if (!rec->v_allele) return;
    if (rng_uniform(cfg->rng) >= cfg->revision_prob) return;
    if (cfg->v_alleles.count < 2) return;

    /* Save original V allele name */
    strncpy(rec->original_v_allele_name, rec->v_allele->name,
            sizeof(rec->original_v_allele_name) - 1);
    rec->original_v_allele_name[sizeof(rec->original_v_allele_name) - 1] = '\0';

    /* Pick a different V allele (different gene) */
    const Allele *new_v = NULL;
    int attempts = 0;
    while (attempts < 50) {
        new_v = allele_pool_random(&cfg->v_alleles, cfg->rng);
        if (strcmp(new_v->gene, rec->v_allele->gene) != 0) break;
        new_v = NULL;
        attempts++;
    }
    if (!new_v) return;

    /* Determine footprint length */
    int fp_min = cfg->footprint_min;
    int fp_max = cfg->footprint_max;
    int fp_len = fp_min + (int)rng_range(cfg->rng,
                                          (uint32_t)(fp_max - fp_min + 1));

    /* Count V segment length */
    int v_len = aseq_segment_length(seq, SEG_V);
    if (fp_len >= v_len) fp_len = v_len - 1;
    if (fp_len < 1) return;

    int replace_len = v_len - fp_len;  /* bases to replace with new V */

    /* Replace V bases (skipping the footprint at the 3' end of V) */
    Nuc *n = seq->seg_first[SEG_V];
    int pos = 0;
    for (; n && n->segment == SEG_V && pos < replace_len; n = n->next, pos++) {
        /* Map to new V allele sequence */
        int new_v_pos = rec->v_trim_5 + pos;
        if (new_v_pos < new_v->length) {
            n->current  = new_v->seq[new_v_pos];
            n->germline = new_v->seq[new_v_pos];
        }
    }

    rec->receptor_revised = true;
    rec->revision_footprint_length = fp_len;
    rec->v_allele = new_v;

    TRACE("[receptor_revision] replaced V: %s → %s (footprint=%dbp, replaced=%dbp)",
          rec->original_v_allele_name, new_v->name, fp_len, replace_len);
}
