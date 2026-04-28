/**
 * loader_auto.c — format-detection dispatcher.
 *
 * Inspects the input path (and a small sniff of file contents) to
 * choose the right concrete ReferenceLoader. Falls back to plain
 * FASTA when nothing else matches.
 *
 * Decision rules (in order):
 *   1. Extension `.json` / `.JSON`              → AIRR-C
 *   2. Sniff first non-whitespace char:
 *        '{' or '['                             → AIRR-C (JSON)
 *        '>'                                    → FASTA family
 *   3. FASTA: look for a sibling `<base>.json`  → OGRDB combined
 *   4. FASTA: header has '|' delimiters         → IMGT V-QUEST
 *   5. Otherwise                                 → plain FASTA
 *
 * Caller passes a `LoaderHints` struct (locus_hint, segment_hint)
 * that's only used by the plain-FASTA path; ignored otherwise.
 */

#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "genairr/ref_loader.h"

/* ── Helpers ─────────────────────────────────────────────────── */

static bool ends_with_ci(const char *s, const char *suffix) {
    size_t n = strlen(s);
    size_t m = strlen(suffix);
    if (m > n) return false;
    for (size_t i = 0; i < m; i++) {
        char a = (char)tolower((unsigned char)s[n - m + i]);
        char b = (char)tolower((unsigned char)suffix[i]);
        if (a != b) return false;
    }
    return true;
}

/* Read up to `max` non-whitespace chars from `path` into `out`.
 * Returns the number of chars written (NUL-terminated). On error,
 * returns -1. Doesn't open the file if max==0. */
static int sniff_first_non_whitespace(const char *path, char *out, int max) {
    FILE *fp = fopen(path, "rb");
    if (!fp) return -1;
    int n = 0;
    int c;
    /* Skip leading whitespace. */
    while ((c = fgetc(fp)) != EOF && isspace((unsigned char)c)) {}
    if (c == EOF) { fclose(fp); out[0] = '\0'; return 0; }
    out[n++] = (char)c;
    while (n < max - 1 && (c = fgetc(fp)) != EOF) {
        out[n++] = (char)c;
    }
    out[n] = '\0';
    fclose(fp);
    return n;
}

/* Build sidecar JSON path by replacing the extension of `fasta_path`
 * with ".json". Caller frees. Returns NULL on alloc failure or paths
 * with no extension (we don't try to add one). */
static char *sidecar_json_path(const char *fasta_path) {
    if (!fasta_path) return NULL;
    size_t n = strlen(fasta_path);
    /* Find the last '.' that's after the last '/' or '\\'. */
    size_t dot = (size_t)-1;
    for (size_t i = n; i > 0; i--) {
        char c = fasta_path[i - 1];
        if (c == '/' || c == '\\') break;
        if (c == '.' && dot == (size_t)-1) dot = i - 1;
    }
    if (dot == (size_t)-1) return NULL;
    /* base = fasta_path[0..dot] + ".json" */
    char *out = (char *)malloc(dot + 6);
    if (!out) return NULL;
    memcpy(out, fasta_path, dot);
    memcpy(out + dot, ".json", 6);
    return out;
}

static bool file_exists(const char *path) {
    if (!path) return false;
    FILE *fp = fopen(path, "r");
    if (!fp) return false;
    fclose(fp);
    return true;
}

/* ── Public dispatcher ─────────────────────────────────────────── */

ReferenceLoader *reference_loader_open_auto(const char *path,
                                             const LoaderHints *hints,
                                             const char **err_msg) {
    if (!path) {
        if (err_msg) *err_msg = "path is NULL";
        return NULL;
    }

    /* 1. JSON extension → AIRR-C. */
    if (ends_with_ci(path, ".json")) {
        return airrc_germline_loader_open(path, err_msg);
    }

    /* 2. Sniff first non-whitespace char. */
    char sniff[8] = {0};
    sniff_first_non_whitespace(path, sniff, sizeof(sniff));
    if (sniff[0] == '{' || sniff[0] == '[') {
        return airrc_germline_loader_open(path, err_msg);
    }

    /* 3-5. FASTA family. We need to peek at the FIRST '>' header line
     * to decide between IMGT (pipe-delimited) and plain. Also probe
     * for an OGRDB sibling JSON. */
    if (sniff[0] != '>') {
        if (err_msg) *err_msg = "could not detect reference format "
                                "(file does not start with '>' or '{'/'[' "
                                "for FASTA / JSON)";
        return NULL;
    }

    /* Check for OGRDB sibling JSON. */
    char *sidecar = sidecar_json_path(path);
    bool has_sidecar = sidecar && file_exists(sidecar);
    if (has_sidecar) {
        Locus locus = (hints && hints->locus_hint != LOCUS_UNKNOWN)
                      ? hints->locus_hint
                      : locus_from_filename(path);
        Segment seg = hints ? hints->segment_hint : SEG_UNKNOWN;
        ReferenceLoader *L = ogrdb_loader_open(path, sidecar,
                                                locus, seg, err_msg);
        free(sidecar);
        return L;
    }
    free(sidecar);

    /* Read the first header to disambiguate IMGT (pipe-delimited) vs
     * plain (no pipes). Reuse the sniff buffer with more room. */
    {
        FILE *fp = fopen(path, "rb");
        if (!fp) {
            if (err_msg) *err_msg = "could not open FASTA file for sniff";
            return NULL;
        }
        char header[512];
        int n = 0;
        int c;
        /* Skip leading whitespace. */
        while ((c = fgetc(fp)) != EOF && isspace((unsigned char)c)) {}
        if (c == EOF || c != '>') {
            fclose(fp);
            if (err_msg) *err_msg = "FASTA missing header line";
            return NULL;
        }
        while (n < (int)sizeof(header) - 1
               && (c = fgetc(fp)) != EOF
               && c != '\n' && c != '\r') {
            header[n++] = (char)c;
        }
        header[n] = '\0';
        fclose(fp);

        if (strchr(header, '|')) {
            return imgt_vquest_loader_open(path,
                hints ? hints->segment_hint : SEG_UNKNOWN, err_msg);
        }
    }

    /* 5. Plain FASTA fallback. */
    {
        Locus locus = (hints && hints->locus_hint != LOCUS_UNKNOWN)
                      ? hints->locus_hint
                      : locus_from_filename(path);
        Segment seg = hints ? hints->segment_hint : SEG_UNKNOWN;
        return plain_fasta_loader_open(path, locus, seg, err_msg);
    }
}
