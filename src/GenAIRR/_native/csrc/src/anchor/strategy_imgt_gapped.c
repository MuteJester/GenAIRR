/**
 * strategy_imgt_gapped.c — IMGT-gapped position derivation.
 *
 * The IMGT V-QUEST alignment puts the conserved 2nd-Cys (V anchor) at
 * IMGT amino-acid position 104, which corresponds to gapped nucleotide
 * positions 309–311 (0-based). To get the equivalent UNGAPPED position
 * we count non-gap characters in gapped[0:309].
 *
 * This is deterministic — no `rfind` heuristic, no coincidental-match
 * risk. Unlike the legacy code, we never accept an anchor we couldn't
 * actually validate as a Cys/Trp codon.
 *
 * For J anchors we don't have an analogous fixed-position derivation
 * here — the J motif is inherently 5'-anchored to the conserved Trp/Phe
 * residue, but its absolute IMGT position varies by length of upstream
 * J framework. Use motif search (strategy_motif_j.c) instead.
 */

#include <stddef.h>   /* NULL */
#include "strategies.h"
#include "codon.h"

/* IMGT V-anchor: the 2nd-Cys codon's first base is at 0-based gapped
 * position 309 (IMGT amino-acid 104, codon position 1). */
#define IMGT_V_ANCHOR_GAPPED_START  309
#define IMGT_V_ANCHOR_GAPPED_END    312    /* exclusive */

/* Count non-gap chars in gapped[0:end). Gap character is '.' in IMGT
 * convention. Returns the equivalent ungapped position for whatever
 * is at gapped[end]. */
static int ungapped_position_of(const char *gapped, int end) {
    int n = 0;
    for (int i = 0; i < end; i++) {
        if (gapped[i] != '.') n++;
    }
    return n;
}

AnchorResult try_v_imgt_gapped(const AnchorResolverConfig *cfg,
                               const LoadedAlleleRecord *rec) {
    if (!rec || !rec->gapped_sequence || !rec->gap_convention_imgt) {
        return anchor_make_rejected(NULL, METHOD_IMGT_GAPPED);
    }
    if (rec->gapped_length < IMGT_V_ANCHOR_GAPPED_END) {
        return anchor_make_rejected(
            "gapped sequence shorter than IMGT V anchor position 312",
            METHOD_IMGT_GAPPED);
    }

    /* Verify the gapped codon at the anchor position is itself
     * non-gap — a fully-gapped codon at position 104 means the
     * allele was never mapped to that IMGT cell, so the entire
     * derivation is meaningless. */
    const char *g = rec->gapped_sequence;
    if (g[IMGT_V_ANCHOR_GAPPED_START + 0] == '.' ||
        g[IMGT_V_ANCHOR_GAPPED_START + 1] == '.' ||
        g[IMGT_V_ANCHOR_GAPPED_START + 2] == '.') {
        return anchor_make_rejected(
            "IMGT-gapped position 309-311 contains gap chars",
            METHOD_IMGT_GAPPED);
    }

    int pos = ungapped_position_of(g, IMGT_V_ANCHOR_GAPPED_START);
    if (rec->sequence == NULL || pos + 3 > rec->sequence_length) {
        return anchor_make_rejected(
            "IMGT-derived position is out of bounds for ungapped seq",
            METHOD_IMGT_GAPPED);
    }

    const char *codon = rec->sequence + pos;
    if (codon_is_stop(codon)) {
        return anchor_make_rejected(
            "IMGT-derived V anchor lands on a stop codon",
            METHOD_IMGT_GAPPED);
    }

    AnchorResidue residue = codon_classify(codon);
    if (residue == ANCHOR_CYS) {
        return anchor_make_result(pos, codon, residue,
                                  CONF_CANONICAL, METHOD_IMGT_GAPPED);
    }
    if (residue == ANCHOR_TRP) {
        if (cfg && !cfg->allow_trp_v) {
            return anchor_make_rejected(
                "Trp at IMGT V anchor disallowed by resolver config",
                METHOD_IMGT_GAPPED);
        }
        return anchor_make_result(pos, codon, residue,
                                  CONF_ALTERNATIVE, METHOD_IMGT_GAPPED);
    }
    if (codon_has_iupac(codon)) {
        return anchor_make_result(pos, codon, ANCHOR_RES_UNKNOWN,
                                  CONF_ALTERNATIVE, METHOD_IMGT_GAPPED);
    }
    return anchor_make_rejected(
        "IMGT-derived V anchor codon is neither Cys nor Trp",
        METHOD_IMGT_GAPPED);
}
