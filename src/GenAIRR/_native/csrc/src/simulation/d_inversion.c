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
#include "genairr/productivity_guard.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <string.h>

/**
 * Compute the inverted-D base array and (optionally) the reversed
 * germline_pos array without mutating `seq`. Returns the D segment
 * length, or 0 if the segment is missing or trivially short.
 */
static int compute_d_inversion(const ASeq *seq, char *out_inverted,
                               uint16_t *out_gpos_reversed) {
    const Nuc *first = seq->seg_first[SEG_D];
    const Nuc *last  = seq->seg_last[SEG_D];
    if (!first || !last) return 0;

    int d_len = aseq_segment_length(seq, SEG_D);
    if (d_len <= 1) return 0;

    char bases[GENAIRR_MAX_ALLELE_SEQ];
    uint16_t gpos[GENAIRR_MAX_ALLELE_SEQ];
    int i = 0;
    for (const Nuc *n = first; n; n = n->next) {
        if (n->segment != SEG_D) break;
        bases[i] = n->current;
        gpos[i]  = n->germline_pos;
        i++;
    }

    for (int k = 0; k < d_len; k++) {
        out_inverted[k] = complement(bases[d_len - 1 - k]);
        if (out_gpos_reversed) {
            out_gpos_reversed[k] = gpos[d_len - 1 - k];
        }
    }
    return d_len;
}

static int apply_d_inversion(ASeq *seq, SimRecord *rec) {
    Nuc *first = seq->seg_first[SEG_D];
    Nuc *last  = seq->seg_last[SEG_D];
    if (!first || !last) return 0;

    /* Collect D bases and germline positions. At this point in the
     * pipeline (pre-SHM) `n->current == n->germline` for every D
     * node, so collecting `n->current` captures the canonical bases
     * for use as the new (inverted) germline reference below. */
    int d_len = aseq_segment_length(seq, SEG_D);
    if (d_len <= 1) return 0;

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
        n->current      = inverted;
        n->germline     = inverted;
        n->germline_pos = gpos[j];
        n->flags       |= NUC_FLAG_INVERTED;   /* T2-13 provenance flag */
        j--;
    }

    /* D inversion is a bulk rewrite over multiple coding nodes.
     * Rebuild the codon rail so cached amino acids, stop counts, and
     * per-node productive flags reflect the inverted template. */
    if (seq->codon_rail_valid) {
        aseq_build_codon_rail(seq);
    }

    rec->d_inverted = true;
    return d_len;
}

void step_d_inversion(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    if (!rec->d_allele) return;
    if (rng_uniform(cfg->rng) >= cfg->d_inversion_prob) return;

    if (simcfg_productive_only(cfg)) {
        char inverted[GENAIRR_MAX_ALLELE_SEQ];
        int d_len = compute_d_inversion(seq, inverted, NULL);
        if (d_len <= 1) return;

        ProductivityDecision decision = productivity_guard_span_rewrite(
            cfg, seq, rec, PROD_STAGE_MOLECULE, SEG_D,
            0, d_len, inverted);
        if (decision != PROD_DECISION_ALLOW) {
            TRACE("[d_inversion] skipped productive-unsafe inversion (allele=%s)",
                  rec->d_allele->name);
            return;
        }
    }

    int d_len = apply_d_inversion(seq, rec);
    if (d_len <= 1) return;

    TRACE("[d_inversion] reverse-complemented D segment (%dbp, allele=%s)",
          d_len, rec->d_allele->name);
}
