/**
 * codon.c — codon classification for V/J anchor resolution.
 *
 * Tiny, self-contained helpers used by every anchor strategy. The
 * surface is small (5 functions) and intentionally case-insensitive so
 * upstream loaders can hand us either case without an extra normalize
 * pass.
 *
 * IUPAC codes are handled conservatively: Y (C-or-T) at codon position
 * 3 is the most common ambiguity in real IMGT data and resolves
 * unambiguously to Cys (TG[CT]) or Phe (TT[CT]); anything else makes
 * the classifier return ANCHOR_RES_UNKNOWN — the caller's strategy can
 * still emit CONF_ALTERNATIVE if the position is otherwise canonical.
 */

#include <ctype.h>
#include "codon.h"

void codon_normalize(const char *src, char out[4]) {
    for (int i = 0; i < 3; i++) {
        out[i] = (char)tolower((unsigned char)src[i]);
    }
    out[3] = '\0';
}

static bool is_canonical_base(char c) {
    return c == 'a' || c == 'c' || c == 'g' || c == 't';
}

bool codon_is_canonical(const char *codon) {
    char c[4];
    codon_normalize(codon, c);
    return is_canonical_base(c[0])
        && is_canonical_base(c[1])
        && is_canonical_base(c[2]);
}

static bool is_iupac_ambiguity(char c) {
    /* Subset of IUPAC ambiguity codes that we accept as "ambiguous but
     * known". 'z' / numerics / punctuation fall through and the
     * caller treats them as garbage via codon_classify. */
    switch (c) {
        case 'r': case 'y': case 's': case 'w': case 'k': case 'm':
        case 'b': case 'd': case 'h': case 'v': case 'n':
            return true;
        default:
            return false;
    }
}

bool codon_has_iupac(const char *codon) {
    char c[4];
    codon_normalize(codon, c);
    for (int i = 0; i < 3; i++) {
        if (is_iupac_ambiguity(c[i])) return true;
    }
    return false;
}

AnchorResidue codon_classify(const char *codon) {
    char c[4];
    codon_normalize(codon, c);

    /* All anchor residues we care about start with 't'. */
    if (c[0] != 't') return ANCHOR_RES_UNKNOWN;

    if (c[1] == 'g') {
        /* TG[TC] = Cys, TG[Y] = Cys with Y resolving to T-or-C. */
        if (c[2] == 't' || c[2] == 'c' || c[2] == 'y') return ANCHOR_CYS;
        if (c[2] == 'g') return ANCHOR_TRP;
        return ANCHOR_RES_UNKNOWN;
    }
    if (c[1] == 't') {
        /* TT[TC] = Phe, TT[Y] = Phe with Y resolving to T-or-C. */
        if (c[2] == 't' || c[2] == 'c' || c[2] == 'y') return ANCHOR_PHE;
        return ANCHOR_RES_UNKNOWN;
    }
    return ANCHOR_RES_UNKNOWN;
}

bool codon_is_stop(const char *codon) {
    char c[4];
    codon_normalize(codon, c);
    if (c[0] != 't') return false;

    /* TAA, TAG → "ta[ag]"; TAR (R=A/G) collapses both. */
    if (c[1] == 'a' && (c[2] == 'a' || c[2] == 'g' || c[2] == 'r')) {
        return true;
    }
    /* TGA */
    if (c[1] == 'g' && c[2] == 'a') {
        return true;
    }
    return false;
}
