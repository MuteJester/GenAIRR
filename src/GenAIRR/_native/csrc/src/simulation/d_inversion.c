/**
 * d_inversion.c — D gene inversion (reverse-complement D segment in-place).
 *
 * BIOLOGY
 * -------
 * V(D)J recombination by RAG1/RAG2 cleaves at the two RSSs flanking
 * the chosen segments. When the RSS orientations match (both 12-bp
 * spacers on D, so bidirectional D usage IS allowed), joining
 * proceeds by INVERSION rather than excision: the intervening DNA is
 * physically flipped end-to-end and retained in the chromosome on
 * the opposite strand. (Lewis 1994 PMC231929; Chi et al. 2020 review
 * PMC7341547.)
 *
 * The lymphocyte's genomic DNA at that locus then reads the REVERSE
 * COMPLEMENT of the canonical D allele on the sense strand. There
 * is no version of the cell that retains the canonical-orientation
 * D bases — that is the defining feature of inversional vs.
 * deletional V(D)J. Recent human BCR repertoire data confirm the
 * mechanism is biologically real (Hu et al. 2025, Comm Biol —
 * 25 unique inverted D usages observed).
 *
 * METADATA INVARIANT
 * ------------------
 * Per AIRR Rearrangement Schema (and matching what every real
 * AIRR-seq tool — IgBLAST, partis, MiXCR — emits): `germline_alignment`
 * must be column-wise consistent with `sequence_alignment`. Since the
 * read carries inverted bases at inverted-D positions (that is what
 * is actually templated in the cell's DNA after rearrangement), the
 * germline column must also be the inverted base. AID/SHM acts on
 * the rearranged template, so the "wild-type" reference for any
 * downstream mutation event on an inverted D base IS the inverted
 * base — NOT the canonical D allele base.
 *
 * Therefore in this step we set BOTH `n->current` and `n->germline`
 * to the reverse-complement of the canonical base. The canonical-
 * orientation provenance is preserved separately via:
 *   - `rec->d_inverted = true`         (record-level flag)
 *   - `n->flags |= NUC_FLAG_INVERTED` (per-position flag, T2-13)
 *   - `rec->d_allele->name`            (the canonical allele used)
 *
 * Audit note (T2-13 resolution): the original audit suggested keeping
 * `n->germline` canonical and only flipping `n->current`. That would
 * produce AIRR records where `germline_alignment != sequence_alignment`
 * over the entire inverted D region even with no SHM — every inverted
 * base would look like an SHM mutation against the canonical reference,
 * and column-wise mutation calling would be biologically incorrect.
 * The fix kept here matches biology + AIRR spec; the audit's per-
 * position flag is added for forensic / introspection use only.
 *
 * Pipeline ordering: this step must run BEFORE mutation (S5F /
 * uniform), because SHM acts on the post-rearrangement template.
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

    /* Collect D bases and germline positions. At this point in the
     * pipeline (pre-SHM) `n->current == n->germline` for every D
     * node, so collecting `n->current` captures the canonical bases
     * for use as the new (inverted) germline reference below. */
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

    /* Write back in reverse-complement, reversing germline_pos so
     * AlleleBitmap lookups remain correct. Both `current` and
     * `germline` get the inverted base — the new germline reference
     * for this region IS the inverted form (see file header). */
    int j = d_len - 1;
    for (Nuc *n = first; n; n = n->next) {
        if (n->segment != SEG_D) break;
        char inverted = complement(bases[j]);
        n->current     = inverted;
        n->germline    = inverted;
        n->germline_pos = gpos[j];
        n->flags      |= NUC_FLAG_INVERTED;   /* T2-13 provenance flag */
        j--;
    }

    rec->d_inverted = true;
    TRACE("[d_inversion] reverse-complemented D segment (%dbp, allele=%s)",
          d_len, rec->d_allele->name);
}
