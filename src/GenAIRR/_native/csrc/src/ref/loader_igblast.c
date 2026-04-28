/**
 * loader_igblast.c — IgBLAST custom-DB bundle reader.
 *
 * IgBLAST (NCBI) ships its references as a 3-component bundle, not a
 * single file:
 *
 *   - <bundle>/<species>_V.fasta    (and ..._D.fasta, ..._J.fasta)
 *     Plain ungapped FASTA, headers are just allele names. Functional
 *     status is NOT carried in the FASTA.
 *
 *   - <bundle>/internal_data/<species>/<species>.ndm.imgt
 *     Tab-delimited V-region delineation: per-allele FWR/CDR
 *     coordinates and the three conserved residue positions
 *     (Cys23, Trp41, Cys104). One row per V allele.
 *
 *   - <bundle>/optional_file/<species>_gl.aux
 *     Tab-delimited J-anchor aux file. Columns:
 *       1: J allele name
 *       2: 0-based coding-frame start
 *       3: chain type (e.g. "VH", "VK", "VL", "VB", "VA")
 *       4: 0-based CDR3 end position (= position of last nt of the
 *          conserved Trp/Phe codon)
 *
 * NB: IgBLAST coordinates are 0-based; AIRR-C and IMGT are 1-based.
 * Our internal `explicit_anchor` is 0-based ungapped, so the .aux
 * conversion is just (cdr3_end_0based - 2) — the conserved codon
 * starts 2 nt before its last nt.
 *
 * Phase 1: support a SINGLE allele type per loader (V or J or D).
 * The caller passes the FASTA path AND optionally the matching aux
 * / ndm sidecar paths. Auto-derivation of sidecar paths from the
 * FASTA filename is left for the auto-dispatcher.
 *
 * For now, this loader is a thin wrapper around the plain FASTA
 * loader plus a per-allele anchor lookup from the aux/ndm sidecar.
 * The plain loader handles the FASTA streaming; we patch the
 * `explicit_anchor` field in each yielded record.
 */

#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "genairr/ref_loader.h"

#define IGBLAST_LINE_INITIAL_CAP   512

/* ── Sidecar map (allele name → 0-based anchor) ──────────────── */

typedef struct {
    char *name;
    int   anchor_0based;
} IgblastSidecarEntry;

typedef struct {
    IgblastSidecarEntry *items;
    int count;
    int cap;
} IgblastSidecarIndex;

static int igblast_sidecar_grow(IgblastSidecarIndex *idx) {
    int new_cap = idx->cap ? idx->cap * 2 : 16;
    IgblastSidecarEntry *grown = (IgblastSidecarEntry *)realloc(
        idx->items, sizeof(IgblastSidecarEntry) * (size_t)new_cap);
    if (!grown) return -1;
    idx->items = grown;
    idx->cap = new_cap;
    return 0;
}

static const IgblastSidecarEntry *igblast_sidecar_find(
    const IgblastSidecarIndex *idx, const char *name) {
    for (int i = 0; i < idx->count; i++) {
        if (strcmp(idx->items[i].name, name) == 0) return &idx->items[i];
    }
    return NULL;
}

static void igblast_sidecar_destroy(IgblastSidecarIndex *idx) {
    if (!idx->items) return;
    for (int i = 0; i < idx->count; i++) free(idx->items[i].name);
    free(idx->items);
    idx->items = NULL;
    idx->count = idx->cap = 0;
}

/* Parse a `_gl.aux` J-anchor file. Format (tab-delimited):
 *   J_allele_name <TAB> coding_frame_start <TAB> chain_type <TAB> cdr3_end
 * cdr3_end is 0-based; conserved codon starts at cdr3_end - 2. */
static int load_aux_j(IgblastSidecarIndex *idx,
                       const char *aux_path,
                       const char **err_msg) {
    FILE *fp = fopen(aux_path, "r");
    if (!fp) {
        if (err_msg) *err_msg = "could not open IgBLAST .aux J file";
        return -1;
    }
    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        /* Skip comments / blanks. */
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0') continue;

        /* tokens[0] = name, [3] = cdr3_end. Tab-split. */
        char *name = strtok(p, "\t");
        char *frame_s = strtok(NULL, "\t");   (void)frame_s;
        char *chain_s = strtok(NULL, "\t");   (void)chain_s;
        char *end_s   = strtok(NULL, "\t\r\n");
        if (!name || !end_s) continue;

        int cdr3_end_0based = atoi(end_s);
        if (cdr3_end_0based < 2) continue;

        /* Conserved Trp/Phe codon STARTS 2 nt before cdr3_end. */
        int anchor_0based = cdr3_end_0based - 2;

        if (idx->count >= idx->cap && igblast_sidecar_grow(idx) != 0) {
            fclose(fp);
            if (err_msg) *err_msg = "OOM building J aux index";
            return -1;
        }
        IgblastSidecarEntry *e = &idx->items[idx->count++];
        e->name = (char *)malloc(strlen(name) + 1);
        if (!e->name) { fclose(fp); if (err_msg) *err_msg = "OOM"; return -1; }
        memcpy(e->name, name, strlen(name) + 1);
        e->anchor_0based = anchor_0based;
    }
    fclose(fp);
    return 0;
}

/* Parse a `.ndm.imgt` V-delineation file. Format (tab-delimited, 1-based):
 *   allele_name <TAB> chain_type <TAB> fwr1_start <TAB> fwr1_end
 *   <TAB> cdr1_start ... <TAB> fwr3_end <TAB> cys_start <TAB> cys_end ... etc.
 *
 * The exact column layout differs across IgBLAST releases; we parse
 * forgivingly: take the LAST integer column on the line as the
 * conserved Cys position (1-based start) — this is correct for the
 * common IgBLAST format where the trailing column is `cdr3_start`.
 * If the line has fewer than 4 fields, skip. */
static int load_ndm_v(IgblastSidecarIndex *idx,
                       const char *ndm_path,
                       const char **err_msg) {
    FILE *fp = fopen(ndm_path, "r");
    if (!fp) {
        if (err_msg) *err_msg = "could not open IgBLAST .ndm.imgt V file";
        return -1;
    }
    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0') continue;

        /* Strip trailing newline/CR. */
        size_t L = strlen(p);
        while (L > 0 && (p[L-1] == '\n' || p[L-1] == '\r')) p[--L] = '\0';

        /* Split by tabs; remember the first token (allele name) and
         * the last token (the conserved-Cys start, 1-based). */
        char *saveptr = NULL;
        (void)saveptr;
        char *name = strtok(p, "\t");
        if (!name) continue;
        char *last = NULL, *tok;
        while ((tok = strtok(NULL, "\t")) != NULL) last = tok;
        if (!last) continue;
        int cys_1based = atoi(last);
        if (cys_1based <= 0) continue;
        int anchor_0based = cys_1based - 1;

        if (idx->count >= idx->cap && igblast_sidecar_grow(idx) != 0) {
            fclose(fp);
            if (err_msg) *err_msg = "OOM building V ndm index";
            return -1;
        }
        IgblastSidecarEntry *e = &idx->items[idx->count++];
        e->name = (char *)malloc(strlen(name) + 1);
        if (!e->name) { fclose(fp); if (err_msg) *err_msg = "OOM"; return -1; }
        memcpy(e->name, name, strlen(name) + 1);
        e->anchor_0based = anchor_0based;
    }
    fclose(fp);
    return 0;
}

/* ── State ─────────────────────────────────────────────────────── */

typedef struct {
    /* Wrapped plain-FASTA loader handles the streaming. */
    ReferenceLoader *inner;
    IgblastSidecarIndex sidecar;
    Locus locus_hint;
    Segment segment_hint;
} IgblastState;

/* ── Vtable ops ────────────────────────────────────────────────── */

static int igblast_next(ReferenceLoader *self,
                        LoadedAlleleRecord *out,
                        const char **err_msg) {
    IgblastState *st = (IgblastState *)self->state;
    if (!st || !st->inner) {
        if (err_msg) *err_msg = "loader has been closed";
        return -1;
    }
    int rc = reference_loader_next(st->inner, out, err_msg);
    if (rc != 0) return rc;

    /* Decorate: locus / segment from hint, anchor from sidecar. */
    if (out->locus == LOCUS_UNKNOWN) out->locus = st->locus_hint;
    if (out->segment == SEG_UNKNOWN) out->segment = st->segment_hint;
    out->source = "igblast-bundle";

    const IgblastSidecarEntry *e = igblast_sidecar_find(&st->sidecar,
                                                         out->name);
    if (e) out->explicit_anchor = e->anchor_0based;
    return 0;
}

static void igblast_close(ReferenceLoader *self) {
    if (!self) return;
    IgblastState *st = (IgblastState *)self->state;
    if (st) {
        if (st->inner) reference_loader_close(st->inner);
        igblast_sidecar_destroy(&st->sidecar);
        free(st);
    }
    free(self);
}

static const ReferenceLoaderVTable IGBLAST_VTABLE = {
    .next  = igblast_next,
    .close = igblast_close,
};

/* ── Factory ───────────────────────────────────────────────────── */

/* Phase 1 IgBLAST loader: caller passes a single FASTA + optional
 * .ndm.imgt OR _gl.aux sidecar. Segment-segregated databases
 * (separate V / J / D FASTAs) are the IgBLAST norm, so the bundle
 * boundary is resolved at the dispatcher level (auto / Python). */
GENAIRR_EXPORT ReferenceLoader *igblast_loader_open(
    const char *fasta_path,
    const char *aux_or_ndm_path,    /* may be NULL */
    Locus locus_hint,
    Segment segment_hint,
    const char **err_msg) {

    if (!fasta_path) {
        if (err_msg) *err_msg = "fasta_path is NULL";
        return NULL;
    }

    /* Sidecar dispatch: ".aux" → J, ".ndm" → V, NULL → no sidecar. */
    IgblastSidecarIndex sidecar = {0};
    if (aux_or_ndm_path) {
        size_t L = strlen(aux_or_ndm_path);
        if (L >= 4 && (strcmp(aux_or_ndm_path + L - 4, ".aux") == 0
                       || strstr(aux_or_ndm_path, "_gl.aux"))) {
            if (load_aux_j(&sidecar, aux_or_ndm_path, err_msg) != 0) {
                igblast_sidecar_destroy(&sidecar);
                return NULL;
            }
        } else {
            if (load_ndm_v(&sidecar, aux_or_ndm_path, err_msg) != 0) {
                igblast_sidecar_destroy(&sidecar);
                return NULL;
            }
        }
    }

    /* Inner = plain FASTA loader. */
    ReferenceLoader *inner = plain_fasta_loader_open(fasta_path,
        locus_hint, segment_hint, err_msg);
    if (!inner) {
        igblast_sidecar_destroy(&sidecar);
        return NULL;
    }

    IgblastState *st = (IgblastState *)calloc(1, sizeof(IgblastState));
    ReferenceLoader *loader = (ReferenceLoader *)calloc(1, sizeof(ReferenceLoader));
    if (!st || !loader) {
        reference_loader_close(inner);
        igblast_sidecar_destroy(&sidecar);
        free(st); free(loader);
        if (err_msg) *err_msg = "out of memory";
        return NULL;
    }
    st->inner = inner;
    st->sidecar = sidecar;
    st->locus_hint = locus_hint;
    st->segment_hint = segment_hint;

    loader->vt = &IGBLAST_VTABLE;
    loader->state = st;
    return loader;
}
