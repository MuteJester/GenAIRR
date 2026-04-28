/**
 * loader_ogrdb.c — OGRDB combined FASTA + JSON loader.
 *
 * OGRDB ships per-locus references as a pair:
 *   - <locus>.fasta — minimal headers (just allele names), gapped
 *     sequences using IMGT positional convention.
 *   - <locus>.json — AIRR-C-format GermlineSet sidecar with
 *     SequenceDelineationV / j_cdr3_end coordinates and the
 *     `functional` boolean.
 *
 * This loader reads the FASTA in streaming fashion (sequences are
 * the bulk of the data) and decorates each record with anchor +
 * functional-status metadata pulled from the JSON. If the sidecar
 * is missing or doesn't contain a given allele, we fall through to
 * the FASTA-only data and let the AnchorResolver handle anchor
 * resolution via its motif fallback.
 *
 * Aliases / paralogs from the JSON are surfaced on
 * LoadedAlleleRecord.aliases — OGRDB encodes the same biological
 * allele under multiple names sometimes (per agent's research).
 */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "genairr/ref_loader.h"
#include "genairr/ref_record.h"
#include "cJSON.h"

#define OGRDB_LINE_INITIAL_CAP   512
#define OGRDB_NAME_BUF_LEN       128

/* ── JSON sidecar — small in-memory map (name → metadata) ───── */

typedef struct {
    char  *name;
    int    explicit_anchor;   /* -1 if not in JSON */
    Segment segment_hint;     /* from delineation type, may stay UNKNOWN */
    FunctionalStatus functional_status;
    /* Aliases — pointers to JSON-owned strings; owned via the same
     * cJSON tree as the rest of the metadata. */
    char  *aliases_buf;       /* '\0' separated list */
    int    n_aliases;
    /* Pointers into aliases_buf, populated lazily. */
} OgrdbMeta;

typedef struct {
    OgrdbMeta *items;
    int        count;
    int        cap;
} OgrdbMetaIndex;

static int ogrdb_meta_grow(OgrdbMetaIndex *idx) {
    int new_cap = idx->cap ? idx->cap * 2 : 16;
    OgrdbMeta *grown = (OgrdbMeta *)realloc(
        idx->items, sizeof(OgrdbMeta) * (size_t)new_cap);
    if (!grown) return -1;
    idx->items = grown;
    idx->cap = new_cap;
    return 0;
}

static void ogrdb_meta_destroy(OgrdbMetaIndex *idx) {
    if (!idx->items) return;
    for (int i = 0; i < idx->count; i++) {
        free(idx->items[i].name);
        free(idx->items[i].aliases_buf);
    }
    free(idx->items);
    idx->items = NULL;
    idx->count = idx->cap = 0;
}

/* Linear-search lookup. Allele counts are small (~100s per locus) so
 * O(N) is fine; if it ever bottlenecks, swap for a hash table. */
static const OgrdbMeta *ogrdb_meta_find(const OgrdbMetaIndex *idx,
                                         const char *name) {
    for (int i = 0; i < idx->count; i++) {
        if (strcmp(idx->items[i].name, name) == 0) {
            return &idx->items[i];
        }
    }
    return NULL;
}

static char *plain_dup_local(const char *src) {
    if (!src) return NULL;
    size_t n = strlen(src);
    char *out = (char *)malloc(n + 1);
    if (!out) return NULL;
    memcpy(out, src, n + 1);
    return out;
}

/* Extract V anchor (cdr3_start, 1-based) from v_gene_delineations
 * (preferring IMGT scheme), or 0 if not present. */
static int extract_v_anchor(cJSON *allele) {
    cJSON *delins = cJSON_GetObjectItem(allele, "v_gene_delineations");
    if (!cJSON_IsArray(delins)) return -1;
    cJSON *d;
    cJSON_ArrayForEach(d, delins) {
        const char *scheme = cJSON_GetStringValue(
            cJSON_GetObjectItem(d, "delineation_scheme"));
        if (scheme && strcasecmp(scheme, "IMGT") == 0) {
            cJSON *cs = cJSON_GetObjectItem(d, "cdr3_start");
            if (cJSON_IsNumber(cs)) {
                int v = (int)cs->valuedouble;
                if (v > 0) return v - 1;
            }
        }
    }
    return -1;
}

/* Extract J anchor (j_cdr3_end - 3, 0-based codon start) or -1. */
static int extract_j_anchor(cJSON *allele) {
    cJSON *je = cJSON_GetObjectItem(allele, "j_cdr3_end");
    if (cJSON_IsNumber(je)) {
        int v = (int)je->valuedouble;
        if (v >= 3) return v - 3;
    }
    return -1;
}

/* Scan a JSON file for `allele_descriptions[]` and build the lookup
 * map. Tolerates the same top-level shapes as loader_airrc.c
 * (GermlineSet wrapper, single object, bare array). */
static int load_sidecar(OgrdbMetaIndex *idx,
                       const char *json_path,
                       const char **err_msg) {
    FILE *fp = fopen(json_path, "rb");
    if (!fp) {
        /* Sidecar missing is OK — caller decides the policy. */
        if (err_msg) *err_msg = "could not open OGRDB JSON sidecar";
        return -1;
    }
    fseek(fp, 0, SEEK_END);
    long sz = ftell(fp);
    rewind(fp);
    char *raw = (char *)malloc((size_t)sz + 1);
    if (!raw) { fclose(fp); if (err_msg) *err_msg = "OOM"; return -1; }
    fread(raw, 1, (size_t)sz, fp);
    raw[sz] = '\0';
    fclose(fp);

    cJSON *root = cJSON_Parse(raw);
    free(raw);
    if (!root) {
        if (err_msg) *err_msg = "could not parse OGRDB JSON sidecar";
        return -1;
    }

    /* Find allele_descriptions array — try (a) GermlineSet[].allele_descriptions,
     * (b) root.allele_descriptions, (c) root as array. */
    cJSON *allele_iters[8] = {0};
    int n_iters = 0;

    cJSON *gset = cJSON_GetObjectItem(root, "GermlineSet");
    if (cJSON_IsArray(gset)) {
        cJSON *gs;
        cJSON_ArrayForEach(gs, gset) {
            cJSON *ad = cJSON_GetObjectItem(gs, "allele_descriptions");
            if (cJSON_IsArray(ad) && n_iters < 8) {
                allele_iters[n_iters++] = ad;
            }
        }
    } else {
        cJSON *ad = cJSON_GetObjectItem(root, "allele_descriptions");
        if (cJSON_IsArray(ad)) {
            allele_iters[n_iters++] = ad;
        } else if (cJSON_IsArray(root)) {
            allele_iters[n_iters++] = root;
        }
    }

    /* Walk each allele descriptions array. */
    for (int it = 0; it < n_iters; it++) {
        cJSON *allele;
        cJSON_ArrayForEach(allele, allele_iters[it]) {
            const char *label = cJSON_GetStringValue(
                cJSON_GetObjectItem(allele, "label"));
            if (!label) continue;

            if (idx->count >= idx->cap) {
                if (ogrdb_meta_grow(idx) != 0) {
                    cJSON_Delete(root);
                    if (err_msg) *err_msg = "OOM";
                    return -1;
                }
            }
            OgrdbMeta *m = &idx->items[idx->count];
            memset(m, 0, sizeof(*m));
            m->name = plain_dup_local(label);
            m->segment_hint = SEG_UNKNOWN;
            /* segment hint from gene_designation if present */
            const char *gd = cJSON_GetStringValue(
                cJSON_GetObjectItem(allele, "gene_designation"));
            if (gd) {
                if (strstr(gd, "V")) m->segment_hint = SEG_V;
                else if (strstr(gd, "D")) m->segment_hint = SEG_D;
                else if (strstr(gd, "J")) m->segment_hint = SEG_J;
                else if (strstr(gd, "C")) m->segment_hint = SEG_C;
            }
            /* anchor: try V, then J. */
            int anchor = extract_v_anchor(allele);
            if (anchor < 0) anchor = extract_j_anchor(allele);
            m->explicit_anchor = anchor;

            /* functional */
            cJSON *func_j = cJSON_GetObjectItem(allele, "functional");
            if (cJSON_IsBool(func_j)) {
                m->functional_status = cJSON_IsTrue(func_j)
                    ? FUNC_F : FUNC_PSEUDO;
            } else {
                m->functional_status = FUNC_UNKNOWN;
            }
            const char *itype = cJSON_GetStringValue(
                cJSON_GetObjectItem(allele, "inference_type"));
            if (m->functional_status == FUNC_F && itype && strstr(itype, "ORF")) {
                m->functional_status = FUNC_ORF;
            }

            /* aliases */
            cJSON *aliases_j = cJSON_GetObjectItem(allele, "sequence_alias");
            if (cJSON_IsArray(aliases_j)) {
                int total_len = 0;
                cJSON *a;
                cJSON_ArrayForEach(a, aliases_j) {
                    const char *s = cJSON_GetStringValue(a);
                    if (s) total_len += (int)strlen(s) + 1;
                }
                if (total_len > 0) {
                    m->aliases_buf = (char *)malloc((size_t)total_len);
                    if (m->aliases_buf) {
                        char *p = m->aliases_buf;
                        cJSON_ArrayForEach(a, aliases_j) {
                            const char *s = cJSON_GetStringValue(a);
                            if (s) {
                                size_t len = strlen(s);
                                memcpy(p, s, len + 1);
                                p += len + 1;
                                m->n_aliases++;
                            }
                        }
                    }
                }
            }

            if (!m->name) {
                cJSON_Delete(root);
                if (err_msg) *err_msg = "OOM";
                return -1;
            }
            idx->count++;
        }
    }

    cJSON_Delete(root);
    return 0;
}

/* ── State ─────────────────────────────────────────────────────── */

typedef struct {
    char  name[OGRDB_NAME_BUF_LEN];
    bool  populated;
} OgrdbHeaderSlot;

typedef struct {
    FILE   *fp;
    char   *line;
    size_t  line_cap;

    char   *gapped;
    int     gapped_len;
    int     gapped_cap;
    char   *ungapped;
    int     ungapped_len;
    int     ungapped_cap;

    OgrdbHeaderSlot cur, next_header;

    Locus   locus_hint;
    Segment segment_hint;
    bool    eof;

    /* Sidecar metadata — empty if no sidecar provided. */
    OgrdbMetaIndex meta;

    /* Per-record alias array used for the LoadedAlleleRecord output —
     * pointers into the meta's aliases_buf (which is owned across
     * loader lifetime, so safe to expose). */
    const char *aliases_out[REF_MAX_ALIASES];
} OgrdbState;

/* ── FASTA streaming (mirror loader_imgt minus the pipe parsing) ─ */

static int ogrdb_ensure_cap(char **buf, int *cap, int need) {
    if (need <= *cap) return 0;
    int new_cap = *cap ? *cap : 64;
    while (new_cap < need) new_cap *= 2;
    char *grown = (char *)realloc(*buf, (size_t)new_cap);
    if (!grown) return -1;
    *buf = grown;
    *cap = new_cap;
    return 0;
}

static int ogrdb_read_line(FILE *fp, char **buf, size_t *cap) {
    if (*cap < OGRDB_LINE_INITIAL_CAP) {
        char *grown = (char *)realloc(*buf, OGRDB_LINE_INITIAL_CAP);
        if (!grown) return -1;
        *buf = grown;
        *cap = OGRDB_LINE_INITIAL_CAP;
    }
    size_t len = 0;
    int c;
    while ((c = fgetc(fp)) != EOF) {
        if (len + 2 > *cap) {
            size_t new_cap = *cap * 2;
            char *grown = (char *)realloc(*buf, new_cap);
            if (!grown) return -1;
            *buf = grown;
            *cap = new_cap;
        }
        if (c == '\n') {
            (*buf)[len] = '\0';
            if (len > 0 && (*buf)[len - 1] == '\r') (*buf)[len - 1] = '\0';
            return 0;
        }
        (*buf)[len++] = (char)c;
    }
    if (len == 0) return 1;
    (*buf)[len] = '\0';
    if (len > 0 && (*buf)[len - 1] == '\r') (*buf)[len - 1] = '\0';
    return 0;
}

static void ogrdb_parse_header(OgrdbHeaderSlot *slot, const char *line) {
    const char *p = line[0] == '>' ? line + 1 : line;
    while (*p && isspace((unsigned char)*p)) p++;
    int n = 0;
    while (*p && !isspace((unsigned char)*p) && n < OGRDB_NAME_BUF_LEN - 1) {
        slot->name[n++] = *p++;
    }
    slot->name[n] = '\0';
    slot->populated = true;
}

static int ogrdb_append_seq(OgrdbState *st, const char *src) {
    int src_len = (int)strlen(src);
    if (ogrdb_ensure_cap(&st->gapped, &st->gapped_cap,
                         st->gapped_len + src_len + 1) != 0) return -1;
    if (ogrdb_ensure_cap(&st->ungapped, &st->ungapped_cap,
                         st->ungapped_len + src_len + 1) != 0) return -1;
    for (int i = 0; src[i]; i++) {
        char c = (char)tolower((unsigned char)src[i]);
        if (isspace((unsigned char)c)) continue;
        st->gapped[st->gapped_len++] = c;
        if (c != '.' && c != '-') st->ungapped[st->ungapped_len++] = c;
    }
    st->gapped[st->gapped_len] = '\0';
    st->ungapped[st->ungapped_len] = '\0';
    return 0;
}

static void ogrdb_reset_seq(OgrdbState *st) {
    st->gapped_len = 0;
    st->ungapped_len = 0;
    if (st->gapped) st->gapped[0] = '\0';
    if (st->ungapped) st->ungapped[0] = '\0';
}

/* Copy alias pointers from meta into out. */
static int populate_aliases(OgrdbState *st,
                            const OgrdbMeta *meta) {
    int n = 0;
    if (meta && meta->aliases_buf && meta->n_aliases > 0) {
        const char *p = meta->aliases_buf;
        while (n < meta->n_aliases && n < REF_MAX_ALIASES) {
            st->aliases_out[n++] = p;
            p += strlen(p) + 1;
        }
    }
    return n;
}

static void ogrdb_emit_record(OgrdbState *st, LoadedAlleleRecord *out) {
    loaded_allele_record_init(out);
    out->name = st->cur.name;

    /* Look up sidecar metadata. */
    const OgrdbMeta *meta = ogrdb_meta_find(&st->meta, st->cur.name);
    int na = populate_aliases(st, meta);
    for (int i = 0; i < na; i++) out->aliases[i] = st->aliases_out[i];
    out->n_aliases = na;

    Segment seg = segment_from_gene_name(st->cur.name);
    if (seg == SEG_UNKNOWN && meta) seg = meta->segment_hint;
    if (seg == SEG_UNKNOWN) seg = st->segment_hint;
    Locus loc = locus_from_gene_name(st->cur.name);
    if (loc == LOCUS_UNKNOWN) loc = st->locus_hint;

    out->segment = seg;
    out->locus = loc;
    out->species = NULL;
    out->sequence = st->ungapped;
    out->sequence_length = st->ungapped_len;
    out->gapped_sequence = st->gapped;
    out->gapped_length = st->gapped_len;
    /* OGRDB FASTAs use IMGT positional convention when gapped. */
    out->gap_convention_imgt = true;

    out->functional_status = meta ? meta->functional_status : FUNC_UNKNOWN;
    out->explicit_anchor = meta ? meta->explicit_anchor : -1;
    out->source = "ogrdb";
}

/* ── Vtable ops ────────────────────────────────────────────────── */

/* T2-9: filter logic moved to the Python policy layer. Loader now
 * passes every record through with whatever functional_status the
 * sidecar reported (or FUNC_UNKNOWN if no sidecar). */

static int ogrdb_next(ReferenceLoader *self,
                      LoadedAlleleRecord *out,
                      const char **err_msg) {
    OgrdbState *st = (OgrdbState *)self->state;
    if (!st || !st->fp) {
        if (err_msg) *err_msg = "loader has been closed";
        return -1;
    }

    while (1) {
        if (!st->cur.populated) {
            if (st->next_header.populated) {
                st->cur = st->next_header;
                st->next_header.populated = false;
            } else {
                if (st->eof) return 1;
                int rc;
                do {
                    rc = ogrdb_read_line(st->fp, &st->line, &st->line_cap);
                    if (rc < 0) {
                        if (err_msg) *err_msg = "I/O error reading FASTA";
                        return -1;
                    }
                    if (rc == 1) { st->eof = true; return 1; }
                } while (st->line[0] != '>');
                ogrdb_parse_header(&st->cur, st->line);
            }
        }

        ogrdb_reset_seq(st);
        bool hit_new_header = false;
        while (!hit_new_header) {
            int rc = ogrdb_read_line(st->fp, &st->line, &st->line_cap);
            if (rc < 0) {
                if (err_msg) *err_msg = "I/O error reading FASTA";
                return -1;
            }
            if (rc == 1) { st->eof = true; break; }
            if (st->line[0] == '>') {
                ogrdb_parse_header(&st->next_header, st->line);
                hit_new_header = true;
                break;
            }
            if (ogrdb_append_seq(st, st->line) != 0) {
                if (err_msg) *err_msg = "out of memory";
                return -1;
            }
        }

        if (st->ungapped_len == 0) {
            st->cur.populated = false;
            continue;
        }
        /* T2-9: emit every record regardless of functional_status —
         * Python policy layer decides. */
        ogrdb_emit_record(st, out);
        st->cur.populated = false;
        return 0;
    }
}

static void ogrdb_close(ReferenceLoader *self) {
    if (!self) return;
    OgrdbState *st = (OgrdbState *)self->state;
    if (st) {
        if (st->fp) fclose(st->fp);
        free(st->line);
        free(st->gapped);
        free(st->ungapped);
        ogrdb_meta_destroy(&st->meta);
        free(st);
    }
    free(self);
}

static const ReferenceLoaderVTable OGRDB_VTABLE = {
    .next  = ogrdb_next,
    .close = ogrdb_close,
};

/* ── Factory ───────────────────────────────────────────────────── */

ReferenceLoader *ogrdb_loader_open(const char *fasta_path,
                                    const char *json_sidecar_path,
                                    Locus locus_hint,
                                    Segment segment_hint,
                                    const char **err_msg) {
    if (!fasta_path) {
        if (err_msg) *err_msg = "fasta_path is NULL";
        return NULL;
    }
    FILE *fp = fopen(fasta_path, "r");
    if (!fp) {
        if (err_msg) *err_msg = "could not open OGRDB FASTA file";
        return NULL;
    }
    OgrdbState *st = (OgrdbState *)calloc(1, sizeof(OgrdbState));
    ReferenceLoader *loader = (ReferenceLoader *)calloc(1, sizeof(ReferenceLoader));
    if (!st || !loader) {
        fclose(fp);
        free(st);
        free(loader);
        if (err_msg) *err_msg = "out of memory";
        return NULL;
    }
    st->fp = fp;
    st->locus_hint = locus_hint;
    st->segment_hint = segment_hint;

    /* Load sidecar if provided. Failure to parse the sidecar is fatal
     * (caller asked for it explicitly); a NULL sidecar path is fine. */
    if (json_sidecar_path) {
        if (load_sidecar(&st->meta, json_sidecar_path, err_msg) != 0) {
            fclose(fp);
            free(st->line);
            ogrdb_meta_destroy(&st->meta);
            free(st);
            free(loader);
            return NULL;
        }
    }

    loader->vt = &OGRDB_VTABLE;
    loader->state = st;
    return loader;
}
