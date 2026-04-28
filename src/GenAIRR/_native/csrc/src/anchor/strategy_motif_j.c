/**
 * strategy_motif_j.c — locus-aware J-anchor motif scan.
 *
 * The J-region anchor at IMGT amino-acid position 118 is followed by
 * a Gly-X-Gly motif. This gives us the "WGXG" or "FGXG" motif:
 *
 *   J anchor codon  +  Gly codon  +  any codon  +  Gly codon
 *
 * where the J anchor is Trp (TGG) for IGH/TRB/TRD loci and Phe
 * (TTT/TTC) for IGK/IGL/TRA/TRG. We scan all 3 frames; the FIRST
 * (5'-most) in-frame match wins because CDR3 ends at the FIRST
 * conserved residue from the 5' end of J.
 *
 * Locus-awareness controls confidence: matching the expected residue
 * for the locus → CONF_CANONICAL; the opposite residue → CONF_ALTERNATIVE
 * (still biologically meaningful, often happens on edge alleles).
 *
 * If `cfg->locus == LOCUS_UNKNOWN`, both Trp and Phe are accepted with
 * CONF_BEST_GUESS — the result's `residue` field tells the caller which
 * one was found (per design Q1).
 */

#include <stddef.h>   /* NULL */
#include "strategies.h"
#include "codon.h"

static bool is_glycine_codon(const char *c) {
    char buf[4];
    codon_normalize(c, buf);
    /* Gly = GGN → ggt, ggc, gga, ggg, ggn, ggr, ggy, etc.
     * We accept any codon starting with "gg" — that includes the
     * IUPAC-ambiguous third position (ggn, ggw, etc.). */
    return buf[0] == 'g' && buf[1] == 'g';
}

AnchorResult try_j_motif(const AnchorResolverConfig *cfg,
                         const LoadedAlleleRecord *rec) {
    if (!rec || !rec->sequence) {
        return anchor_make_rejected(NULL, METHOD_MOTIF_SEARCH);
    }
    /* Need at least 4 codons (12 nt) for the WGXG / FGXG motif. */
    if (rec->sequence_length < 12) {
        return anchor_make_rejected(
            "J sequence too short for WGXG/FGXG motif",
            METHOD_MOTIF_SEARCH);
    }

    const Locus locus = cfg ? cfg->locus : LOCUS_UNKNOWN;
    /* Conserved residue at IMGT J position 118, per Lefranc et al.
     * (and verified empirically against the shipped IMGT pickles):
     *   IGH                                    → Trp (TGG)
     *   IGK / IGL / TRA / TRB / TRD / TRG     → Phe (TT[CT])
     * TRA does have ~2 rare Trp-anchor alleles in human; those
     * resolve to CONF_ALTERNATIVE under this rule, which is
     * biologically correct labelling for non-canonical anchors. */
    const bool expects_trp = (locus == LOCUS_IGH);
    const bool expects_phe =
        (locus == LOCUS_IGK || locus == LOCUS_IGL ||
         locus == LOCUS_TRA || locus == LOCUS_TRB ||
         locus == LOCUS_TRD || locus == LOCUS_TRG);

    /* Scan all 3 frames; pick the earliest in-frame match. */
    for (int frame = 0; frame < 3; frame++) {
        for (int i = frame; i + 12 <= rec->sequence_length; i += 3) {
            const char *anchor = rec->sequence + i;
            const char *gly1   = rec->sequence + i + 3;
            const char *xxx    = rec->sequence + i + 6;
            const char *gly2   = rec->sequence + i + 9;

            (void)xxx;   /* X is any codon — no constraint */

            AnchorResidue residue = codon_classify(anchor);
            if (residue != ANCHOR_TRP && residue != ANCHOR_PHE) continue;
            if (!is_glycine_codon(gly1) || !is_glycine_codon(gly2)) continue;

            /* Locus-aware confidence: */
            AnchorConfidence conf;
            if (locus == LOCUS_UNKNOWN) {
                /* Q1 decision: ambiguous locus → BEST_GUESS, residue
                 * field reveals which residue was matched. */
                conf = CONF_BEST_GUESS;
            } else if ((residue == ANCHOR_TRP && expects_trp) ||
                       (residue == ANCHOR_PHE && expects_phe)) {
                conf = CONF_CANONICAL;
            } else {
                conf = CONF_ALTERNATIVE;
            }

            return anchor_make_result(i, anchor, residue, conf,
                                      METHOD_MOTIF_SEARCH);
        }
    }

    return anchor_make_rejected(
        "no in-frame WGXG/FGXG motif found in J sequence",
        METHOD_MOTIF_SEARCH);
}
