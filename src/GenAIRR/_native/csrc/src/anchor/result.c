/**
 * result.c — AnchorResult constructors and human-readable enum names.
 *
 * Tiny module — just value constructors and string lookup tables.
 * Lives next to the strategy code because it shares the `codon` helper
 * for normalizing the codon field.
 */

#include <string.h>
#include "genairr/anchor.h"
#include "genairr/ref_record.h"
#include "codon.h"

/* ── AnchorResolverConfig defaults ─────────────────────────────── */

void anchor_resolver_config_init(AnchorResolverConfig *cfg) {
    if (!cfg) return;
    cfg->locus = LOCUS_UNKNOWN;
    cfg->allow_trp_v = true;     /* Trp at V anchor is biology, not error */
    cfg->strict = false;         /* allow CONF_BEST_GUESS results */
    cfg->custom_finder = NULL;
    cfg->custom_user_data = NULL;
}

/* ── LoadedAlleleRecord defaults ───────────────────────────────── */

void loaded_allele_record_init(LoadedAlleleRecord *r) {
    if (!r) return;
    memset(r, 0, sizeof(*r));
    r->segment = SEG_UNKNOWN;
    r->locus = LOCUS_UNKNOWN;
    r->functional_status = FUNC_UNKNOWN;
    r->explicit_anchor = -1;     /* sentinel for "no explicit anchor" */
}

/* ── AnchorResult constructors ─────────────────────────────────── */

AnchorResult anchor_make_rejected(const char *reason, AnchorMethod method) {
    AnchorResult r;
    r.position   = -1;
    r.codon[0]   = '?';
    r.codon[1]   = '?';
    r.codon[2]   = '?';
    r.codon[3]   = '\0';
    r.residue    = ANCHOR_RES_UNKNOWN;
    r.confidence = CONF_REJECTED;
    r.method     = method;
    r.reason     = reason;       /* static literal or NULL — never owned */
    return r;
}

AnchorResult anchor_make_result(int position,
                                const char *codon,
                                AnchorResidue residue,
                                AnchorConfidence confidence,
                                AnchorMethod method) {
    AnchorResult r;
    r.position   = position;
    if (codon) {
        codon_normalize(codon, r.codon);   /* lowercases first 3 chars + NUL */
    } else {
        r.codon[0] = r.codon[1] = r.codon[2] = '?';
        r.codon[3] = '\0';
    }
    r.residue    = residue;
    r.confidence = confidence;
    r.method     = method;
    r.reason     = NULL;
    return r;
}

/* ── Enum → name lookups ───────────────────────────────────────── */

const char *anchor_residue_name(AnchorResidue r) {
    switch (r) {
        case ANCHOR_RES_UNKNOWN: return "unknown";
        case ANCHOR_CYS:         return "Cys";
        case ANCHOR_TRP:         return "Trp";
        case ANCHOR_PHE:         return "Phe";
    }
    return "unknown";
}

const char *anchor_confidence_name(AnchorConfidence c) {
    switch (c) {
        case CONF_REJECTED:    return "rejected";
        case CONF_BEST_GUESS:  return "best_guess";
        case CONF_ALTERNATIVE: return "alternative";
        case CONF_CANONICAL:   return "canonical";
    }
    return "unknown";
}

const char *anchor_method_name(AnchorMethod m) {
    switch (m) {
        case METHOD_NONE:         return "none";
        case METHOD_OVERRIDE:     return "override";
        case METHOD_CUSTOM:       return "custom";
        case METHOD_EXPLICIT:     return "explicit";
        case METHOD_IMGT_GAPPED:  return "imgt_gapped";
        case METHOD_MOTIF_SEARCH: return "motif_search";
    }
    return "unknown";
}

const char *locus_name(Locus l) {
    switch (l) {
        case LOCUS_UNKNOWN: return "unknown";
        case LOCUS_IGH:     return "IGH";
        case LOCUS_IGK:     return "IGK";
        case LOCUS_IGL:     return "IGL";
        case LOCUS_TRA:     return "TRA";
        case LOCUS_TRB:     return "TRB";
        case LOCUS_TRD:     return "TRD";
        case LOCUS_TRG:     return "TRG";
    }
    return "unknown";
}

const char *segment_name(Segment s) {
    /* Only V/D/J/C are meaningful for reference alleles; the simulator-
     * internal node tags (SEG_NP1, SEG_NP2, SEG_UMI, SEG_ADAPTER)
     * never appear on records and collapse to "unknown". SEG_COUNT
     * serves as the SEG_UNKNOWN sentinel (see anchor.h). */
    switch (s) {
        case SEG_V: return "V";
        case SEG_D: return "D";
        case SEG_J: return "J";
        case SEG_C: return "C";
        default:    return "unknown";
    }
}

const char *functional_status_name(FunctionalStatus s) {
    switch (s) {
        case FUNC_UNKNOWN: return "unknown";
        case FUNC_F:       return "F";
        case FUNC_ORF:     return "ORF";
        case FUNC_PSEUDO:  return "P";
        case FUNC_PARTIAL: return "partial";
    }
    return "unknown";
}
