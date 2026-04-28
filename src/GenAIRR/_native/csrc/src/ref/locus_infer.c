/**
 * locus_infer.c — heuristic locus + segment inference from names.
 *
 * IMGT names follow a strict pattern:
 *
 *   <locus 2-3 letters><segment letter><number>[-<sub>]*[ALLELE]
 *
 * Examples:
 *   IGHV1-2*01      → locus=IGH, segment=V
 *   IGKJ5*01        → locus=IGK, segment=J
 *   TRBV20-1*01     → locus=TRB, segment=V
 *   TRBC1*01        → locus=TRB, segment=C
 *   IGHC*01         → locus=IGH, segment=C
 *   IGHE*01         → locus=IGH, segment=C  (special-case: C-region constants)
 *
 * Filename inference is heuristic only — used as a hint when name-based
 * inference fails. We look at the basename and return the first locus
 * marker we find.
 */

#include <ctype.h>
#include <string.h>
#include "genairr/ref_loader.h"

/* Match a 3-character locus prefix at the start of `name`, then return
 * (a) the matching Locus and (b) advance *cursor past the prefix. */
static Locus match_locus_prefix(const char *name, int *cursor) {
    /* IG loci */
    if (!strncmp(name, "IGH", 3)) { *cursor = 3; return LOCUS_IGH; }
    if (!strncmp(name, "IGK", 3)) { *cursor = 3; return LOCUS_IGK; }
    if (!strncmp(name, "IGL", 3)) { *cursor = 3; return LOCUS_IGL; }
    /* TR loci */
    if (!strncmp(name, "TRA", 3)) { *cursor = 3; return LOCUS_TRA; }
    if (!strncmp(name, "TRB", 3)) { *cursor = 3; return LOCUS_TRB; }
    if (!strncmp(name, "TRD", 3)) { *cursor = 3; return LOCUS_TRD; }
    if (!strncmp(name, "TRG", 3)) { *cursor = 3; return LOCUS_TRG; }
    return LOCUS_UNKNOWN;
}

Locus locus_from_gene_name(const char *name) {
    if (!name) return LOCUS_UNKNOWN;
    int cursor = 0;
    return match_locus_prefix(name, &cursor);
}

Segment segment_from_gene_name(const char *name) {
    if (!name) return SEG_UNKNOWN;
    int cursor = 0;
    if (match_locus_prefix(name, &cursor) == LOCUS_UNKNOWN) {
        return SEG_UNKNOWN;
    }
    if (name[cursor] == '\0') return SEG_UNKNOWN;
    switch (name[cursor]) {
        case 'V': return SEG_V;
        case 'D': return SEG_D;
        case 'J': return SEG_J;
        case 'C':
        case 'M': case 'G': case 'A': case 'E':   /* IGHM/IGHG/IGHA/IGHE */
            return SEG_C;
        default:  return SEG_UNKNOWN;
    }
}

/* Find the basename (last path component) of `path`. */
static const char *path_basename(const char *path) {
    if (!path) return NULL;
    const char *p = path;
    const char *last = path;
    for (; *p; p++) {
        if (*p == '/' || *p == '\\') last = p + 1;
    }
    return last;
}

/* Case-insensitive substring search. */
static bool ci_contains(const char *haystack, const char *needle) {
    int hlen = (int)strlen(haystack);
    int nlen = (int)strlen(needle);
    if (nlen == 0) return true;
    for (int i = 0; i + nlen <= hlen; i++) {
        bool match = true;
        for (int j = 0; j < nlen; j++) {
            char h = (char)tolower((unsigned char)haystack[i + j]);
            char n = (char)tolower((unsigned char)needle[j]);
            if (h != n) { match = false; break; }
        }
        if (match) return true;
    }
    return false;
}

Locus locus_from_filename(const char *path) {
    if (!path) return LOCUS_UNKNOWN;
    const char *base = path_basename(path);
    if (!base) return LOCUS_UNKNOWN;

    /* Order matters: longer / more specific first to avoid e.g.
     * "IG" matching inside "IGHV" before "IGH" gets a chance. */
    if (ci_contains(base, "IGH")) return LOCUS_IGH;
    if (ci_contains(base, "IGK")) return LOCUS_IGK;
    if (ci_contains(base, "IGL")) return LOCUS_IGL;
    if (ci_contains(base, "TRA")) return LOCUS_TRA;
    if (ci_contains(base, "TRB")) return LOCUS_TRB;
    if (ci_contains(base, "TRD")) return LOCUS_TRD;
    if (ci_contains(base, "TRG")) return LOCUS_TRG;
    /* Common TCR aliases without TR prefix */
    if (ci_contains(base, "TCRA")) return LOCUS_TRA;
    if (ci_contains(base, "TCRB")) return LOCUS_TRB;
    if (ci_contains(base, "TCRD")) return LOCUS_TRD;
    if (ci_contains(base, "TCRG")) return LOCUS_TRG;
    return LOCUS_UNKNOWN;
}
