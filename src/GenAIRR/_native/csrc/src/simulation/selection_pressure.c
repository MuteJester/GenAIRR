/**
 * selection_pressure.c — Antigen-driven selection after SHM.
 *
 * Classifies each mutation as R (replacement) or S (silent) via codon
 * translation. Silent mutations are always kept. R mutations are
 * accepted with region-dependent probability:
 *   - CDR regions: higher acceptance (default 0.85)
 *   - FWR regions: lower acceptance (default 0.40)
 *
 * Reverted R mutations restore the original germline base.
 * Effective acceptance = 1.0 - strength * (1.0 - base_acceptance).
 *
 * IMGT boundaries (ungapped): FR1=0:78, CDR1=78:114, FR2=114:165,
 * CDR2=165:195, FR3=195:312. Positions ≥312 are CDR3/FR4.
 */

#include "genairr/pipeline.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>

/* ── IMGT boundaries ─────────────────────────────────────────────── */

typedef enum { REGION_FWR, REGION_CDR } RegionType;

static RegionType classify_region(uint16_t germ_pos) {
    if (germ_pos < 78)  return REGION_FWR;  /* FR1 */
    if (germ_pos < 114) return REGION_CDR;  /* CDR1 */
    if (germ_pos < 165) return REGION_FWR;  /* FR2 */
    if (germ_pos < 195) return REGION_CDR;  /* CDR2 */
    if (germ_pos < 312) return REGION_FWR;  /* FR3 */
    return REGION_CDR;                       /* CDR3/FR4 */
}

/* ── Codon translation (compact) ────────────────────────────────── */

static int base_idx(char c) {
    switch (c) {
        case 'T': case 't': return 0;
        case 'C': case 'c': return 1;
        case 'A': case 'a': return 2;
        case 'G': case 'g': return 3;
        default: return -1;
    }
}

static const char CODON_TABLE[64] = {
    'F','F','L','L', 'S','S','S','S', 'Y','Y','*','*', 'C','C','*','W',
    'L','L','L','L', 'P','P','P','P', 'H','H','Q','Q', 'R','R','R','R',
    'I','I','I','M', 'T','T','T','T', 'N','N','K','K', 'S','S','R','R',
    'V','V','V','V', 'A','A','A','A', 'D','D','E','E', 'G','G','G','G',
};

static char translate3(char b1, char b2, char b3) {
    int i1 = base_idx(b1), i2 = base_idx(b2), i3 = base_idx(b3);
    if (i1 < 0 || i2 < 0 || i3 < 0) return '?';
    return CODON_TABLE[i1 * 16 + i2 * 4 + i3];
}

/* ── Check if mutation is R or S ────────────────────────────────── */

/**
 * A mutation is "silent" (S) if the amino acid doesn't change.
 * We need the codon context — find the reading-frame codon that
 * contains this position. Since we need the V-anchor for the
 * reading frame, we check codon in the V-anchor frame.
 */
static bool is_replacement(Nuc *node) {
    /* Find codon boundaries by walking back to a position divisible by 3
     * relative to the V reading frame. We approximate by using the
     * germline position if available. */
    if (node->germline == '\0') return true;  /* NP mutation = replacement */

    /* Collect 3-base codon around this node (approximate: current + neighbors) */
    Nuc *n0 = node;
    Nuc *n1 = node->next;
    Nuc *n2 = n1 ? n1->next : NULL;

    /* Try to align to codon boundary using germline_pos % 3 */
    int offset = node->germline_pos % 3;
    if (offset == 1 && node->prev) {
        n0 = node->prev;
        n1 = node;
        n2 = node->next;
    } else if (offset == 2 && node->prev && node->prev->prev) {
        n0 = node->prev->prev;
        n1 = node->prev;
        n2 = node;
    }

    if (!n0 || !n1 || !n2) return true;

    /* Translate with current bases */
    char aa_current = translate3(n0->current, n1->current, n2->current);

    /* Translate with germline bases (only the mutated one reverted) */
    char c0 = n0->current, c1 = n1->current, c2 = n2->current;
    if (n0 == node) c0 = node->germline;
    if (n1 == node) c1 = node->germline;
    if (n2 == node) c2 = node->germline;
    char aa_germline = translate3(c0, c1, c2);

    return aa_current != aa_germline;
}

/* ── Main step ──────────────────────────────────────────────────── */

void step_selection_pressure(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)rec;
    double strength = cfg->selection_strength;
    if (strength <= 0.0) return;

    double cdr_accept = 1.0 - strength * (1.0 - cfg->cdr_r_acceptance);
    double fwr_accept = 1.0 - strength * (1.0 - cfg->fwr_r_acceptance);

    TRACE("[selection] strength=%.2f, CDR_R_accept=%.2f, FWR_R_accept=%.2f",
          strength, cdr_accept, fwr_accept);

    int reverted_cdr = 0, reverted_fwr = 0, kept_r = 0, kept_s = 0;

    /* Walk all V-segment nodes with mutations */
    for (Nuc *n = seq->seg_first[SEG_V]; n && n->segment == SEG_V; n = n->next) {
        if (!(n->flags & NUC_FLAG_MUTATED)) continue;

        /* Silent mutations always kept */
        if (!is_replacement(n)) { kept_s++; continue; }

        /* Determine acceptance probability based on region */
        RegionType region = classify_region(n->germline_pos);
        double accept = (region == REGION_CDR) ? cdr_accept : fwr_accept;

        if (rng_uniform(cfg->rng) < accept) { kept_r++; continue; }  /* accepted */

        /* Revert: restore germline base */
        aseq_revert(seq, n);
        if (region == REGION_CDR) reverted_cdr++;
        else reverted_fwr++;
    }

    TRACE("[selection] kept %d silent + %d replacement, reverted %d FWR + %d CDR replacements",
          kept_s, kept_r, reverted_fwr, reverted_cdr);
}
