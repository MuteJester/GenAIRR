/**
 * strategy_explicit.c — use the explicit anchor coordinate from the
 * loader if it carried one (AIRR-C JSON sidecar, IgBLAST .aux file).
 *
 * If no explicit anchor was provided, the strategy returns
 * CONF_REJECTED so the next strategy gets a chance.
 *
 * If an explicit anchor IS provided, this strategy validates the
 * codon at that position rather than blindly trusting it. Reasons:
 *   - protects against off-by-one errors when adapting external
 *     coordinate systems (e.g. IgBLAST is 0-based; some custom
 *     formats are 1-based)
 *   - flags pseudogenes whose recorded "anchor" lands on a stop codon
 *     so downstream code can route them to the pseudogene path
 */

#include <stddef.h>   /* NULL */
#include "strategies.h"
#include "codon.h"

static AnchorResult validate_v_codon(int position,
                                     const char *codon,
                                     const AnchorResolverConfig *cfg,
                                     AnchorMethod method) {
    if (codon_is_stop(codon)) {
        return anchor_make_rejected(
            "explicit anchor lands on a stop codon", method);
    }
    AnchorResidue residue = codon_classify(codon);
    if (residue == ANCHOR_CYS) {
        return anchor_make_result(position, codon, residue,
                                  CONF_CANONICAL, method);
    }
    if (residue == ANCHOR_TRP) {
        if (cfg && !cfg->allow_trp_v) {
            return anchor_make_rejected(
                "Trp anchor disallowed by resolver config", method);
        }
        return anchor_make_result(position, codon, residue,
                                  CONF_ALTERNATIVE, method);
    }
    if (codon_has_iupac(codon)) {
        /* IUPAC ambiguity at the anchor — accept with reduced
         * confidence since the residue could still be Cys/Trp. */
        return anchor_make_result(position, codon,
                                  ANCHOR_RES_UNKNOWN,
                                  CONF_ALTERNATIVE, method);
    }
    return anchor_make_rejected(
        "explicit V anchor codon is neither Cys nor Trp", method);
}

AnchorResult try_v_explicit(const AnchorResolverConfig *cfg,
                            const LoadedAlleleRecord *rec) {
    if (!rec || rec->explicit_anchor < 0) {
        return anchor_make_rejected(NULL, METHOD_EXPLICIT);
    }
    const int pos = rec->explicit_anchor;
    if (rec->sequence == NULL || pos + 3 > rec->sequence_length) {
        return anchor_make_rejected(
            "explicit anchor is out of bounds for the ungapped sequence",
            METHOD_EXPLICIT);
    }
    return validate_v_codon(pos, rec->sequence + pos, cfg, METHOD_EXPLICIT);
}

AnchorResult try_j_explicit(const AnchorResolverConfig *cfg,
                            const LoadedAlleleRecord *rec) {
    if (!rec || rec->explicit_anchor < 0) {
        return anchor_make_rejected(NULL, METHOD_EXPLICIT);
    }
    const int pos = rec->explicit_anchor;
    if (rec->sequence == NULL || pos + 3 > rec->sequence_length) {
        return anchor_make_rejected(
            "explicit J anchor is out of bounds", METHOD_EXPLICIT);
    }
    const char *codon = rec->sequence + pos;
    if (codon_is_stop(codon)) {
        return anchor_make_rejected(
            "explicit J anchor lands on a stop codon", METHOD_EXPLICIT);
    }
    AnchorResidue residue = codon_classify(codon);

    /* Locus-aware: TGG (Trp) is canonical only for IGH J;
     * TT[CT] (Phe) is canonical for every other locus
     * (IGK / IGL / TRA / TRB / TRD / TRG). See
     * strategy_motif_j.c for the full citation. */
    bool expects_trp = (cfg && cfg->locus == LOCUS_IGH);
    if (residue == ANCHOR_TRP) {
        return anchor_make_result(pos, codon, residue,
                                  expects_trp ? CONF_CANONICAL : CONF_ALTERNATIVE,
                                  METHOD_EXPLICIT);
    }
    if (residue == ANCHOR_PHE) {
        return anchor_make_result(pos, codon, residue,
                                  expects_trp ? CONF_ALTERNATIVE : CONF_CANONICAL,
                                  METHOD_EXPLICIT);
    }
    if (codon_has_iupac(codon)) {
        return anchor_make_result(pos, codon, ANCHOR_RES_UNKNOWN,
                                  CONF_ALTERNATIVE, METHOD_EXPLICIT);
    }
    return anchor_make_rejected(
        "explicit J anchor codon is neither Trp nor Phe", METHOD_EXPLICIT);
}
