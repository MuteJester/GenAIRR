/**
 * loader_airrc.c — AIRR-Community GermlineSet JSON reader.
 *
 * AIRR-C JSON is the modern source-of-truth for AIRR-format references
 * (https://docs.airr-community.org/en/stable/datarep/germline.html).
 * Unlike IMGT FASTA, it carries explicit anchor coordinates inline.
 * Top-level structure (excerpted):
 *
 *   {
 *     "GermlineSet": [
 *       {
 *         "germline_set_name": "Homo sapiens IGH 1.0.0",
 *         "species": { "label": "Homo sapiens", ... },
 *         "allele_descriptions": [
 *           {
 *             "label": "IGHV1-2*02",
 *             "sequence_alias": [...],
 *             "locus": "IGH",
 *             "v_gene_delineations": [
 *               {
 *                 "delineation_scheme": "IMGT",
 *                 "fwr1_start": 1, "fwr1_end": 78,
 *                 ...
 *                 "cdr3_start": 295         // V anchor (1-based) — first nt of CDR3
 *               }, ...
 *             ],
 *             "j_codon_frame": 1,
 *             "j_cdr3_end": 12,             // J anchor (1-based ungapped)
 *             "functional": true,
 *             "inference_type": "Genomic and rearranged",
 *             ...
 *           }, ...
 *         ]
 *       }
 *     ]
 *   }
 *
 * NB: AIRR-C coordinates are 1-based and refer to ``unaligned_sequence``
 * (the ungapped form). We convert to 0-based on emit. The V "anchor"
 * we use is `cdr3_start - 1` (0-based first nt of CDR3 = first nt of
 * the conserved Cys/Trp codon). The J "anchor" we use is the codon
 * start derived from `j_cdr3_end` minus 3 (the conserved Trp/Phe
 * codon ENDS at j_cdr3_end), then converted to 0-based.
 *
 * Allele records are read all at once into memory at open time
 * (JSON is small enough — even large germline sets are < 5 MB).
 * The loader streams them out one at a time on next() calls.
 */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "genairr/ref_loader.h"
#include "genairr/ref_record.h"
#include "cJSON.h"

/* ── Per-allele staged record (loader-owned strings) ──────────── */

typedef struct {
    char  *name;
    char  *species;
    char  *sequence;          /* lowercased ungapped */
    int    sequence_length;
    Locus  locus;
    Segment segment;
    FunctionalStatus functional_status;
    int    explicit_anchor;   /* -1 if not in JSON */
} AirrcStaged;

typedef struct {
    AirrcStaged *records;
    int          count;
    int          cap;

    /* Latch buffers for the currently-emitted record so
     * LoadedAlleleRecord pointers stay valid through the next() call. */
    int          cursor;       /* index of the NEXT record to emit */
    /* The currently-emitted record's strings are pointed at directly
     * from records[cursor-1]; nothing extra to copy. */
} AirrcState;

/* ── Helpers ─────────────────────────────────────────────────── */

/* Read the whole file into a malloc'd NUL-terminated buffer. Caller
 * frees. Returns NULL on error. */
static char *slurp_file(const char *path, const char **err_msg) {
    FILE *fp = fopen(path, "rb");
    if (!fp) {
        if (err_msg) *err_msg = "could not open JSON file";
        return NULL;
    }
    if (fseek(fp, 0, SEEK_END) != 0) {
        fclose(fp);
        if (err_msg) *err_msg = "could not seek JSON file";
        return NULL;
    }
    long size = ftell(fp);
    if (size < 0) {
        fclose(fp);
        if (err_msg) *err_msg = "could not size JSON file";
        return NULL;
    }
    rewind(fp);
    char *buf = (char *)malloc((size_t)size + 1);
    if (!buf) {
        fclose(fp);
        if (err_msg) *err_msg = "out of memory reading JSON";
        return NULL;
    }
    size_t got = fread(buf, 1, (size_t)size, fp);
    fclose(fp);
    if (got != (size_t)size) {
        free(buf);
        if (err_msg) *err_msg = "short read on JSON file";
        return NULL;
    }
    buf[size] = '\0';
    return buf;
}

/* Lowercase + copy. Returns malloc'd NUL-terminated string, or NULL
 * on alloc failure. Returns NULL for NULL input too. */
static char *lower_dup(const char *src) {
    if (!src) return NULL;
    size_t n = strlen(src);
    char *out = (char *)malloc(n + 1);
    if (!out) return NULL;
    for (size_t i = 0; i < n; i++) out[i] = (char)tolower((unsigned char)src[i]);
    out[n] = '\0';
    return out;
}

static char *plain_dup(const char *src) {
    if (!src) return NULL;
    size_t n = strlen(src);
    char *out = (char *)malloc(n + 1);
    if (!out) return NULL;
    memcpy(out, src, n + 1);
    return out;
}

/* AIRR-C strings the locus as a plain "IGH"/"IGK"/.../"TRG" string.
 * Map to the Locus enum. Tolerate "IGHV"/"TRBV"-prefixed variants if
 * they ever appear (they shouldn't per spec, but safety net). */
static Locus parse_locus_string(const char *s) {
    if (!s) return LOCUS_UNKNOWN;
    if (!strncasecmp(s, "IGH", 3)) return LOCUS_IGH;
    if (!strncasecmp(s, "IGK", 3)) return LOCUS_IGK;
    if (!strncasecmp(s, "IGL", 3)) return LOCUS_IGL;
    if (!strncasecmp(s, "TRA", 3)) return LOCUS_TRA;
    if (!strncasecmp(s, "TRB", 3)) return LOCUS_TRB;
    if (!strncasecmp(s, "TRD", 3)) return LOCUS_TRD;
    if (!strncasecmp(s, "TRG", 3)) return LOCUS_TRG;
    return LOCUS_UNKNOWN;
}

/* Best-effort segment from an allele label. Re-uses
 * `segment_from_gene_name`. */

/* Parse one allele description from JSON, append to state. Returns 0
 * on success, -1 on error (state is not rolled back — caller should
 * close the loader on error). Skipped (returns 0 with no append) when
 * the allele lacks the minimum fields. */
static int stage_allele(AirrcState *st,
                        cJSON *allele,
                        const char *species,
                        const char **err_msg) {
    /* Required: label + sequence (or coding_sequence). */
    cJSON *label_j = cJSON_GetObjectItem(allele, "label");
    const char *label = cJSON_GetStringValue(label_j);
    if (!label || !*label) return 0;   /* skip nameless */

    /* Sequence: prefer "sequence" (full), fallback to "coding_sequence".
     * AIRR-C "sequence" is the unaligned (ungapped) form. */
    const char *seq = cJSON_GetStringValue(
        cJSON_GetObjectItem(allele, "sequence"));
    if (!seq) {
        seq = cJSON_GetStringValue(
            cJSON_GetObjectItem(allele, "coding_sequence"));
    }
    if (!seq || !*seq) return 0;        /* skip empty */

    /* Locus: "locus" is a plain string in modern AIRR-C. */
    Locus locus = parse_locus_string(
        cJSON_GetStringValue(cJSON_GetObjectItem(allele, "locus")));
    /* Fallback: infer from label prefix. */
    if (locus == LOCUS_UNKNOWN) locus = locus_from_gene_name(label);

    /* Segment: AIRR-C carries it as part of the gene_designation field
     * for V/D/J typing. Just sniff the label letter — every modern
     * AIRR-C name embeds it. */
    Segment segment = segment_from_gene_name(label);

    /* Functional status: AIRR-C is boolean, but ORF distinction is in
     * `inference_type`. We project to our tri-state enum:
     *   functional==true               → FUNC_F (or FUNC_ORF if hinted)
     *   functional==false              → FUNC_PSEUDO
     *   missing                        → FUNC_UNKNOWN
     * Only kept records (F/ORF) actually emit later — the policy in
     * Phase 1A skips P; we keep that.
     */
    FunctionalStatus status = FUNC_UNKNOWN;
    cJSON *func_j = cJSON_GetObjectItem(allele, "functional");
    if (cJSON_IsBool(func_j)) {
        status = cJSON_IsTrue(func_j) ? FUNC_F : FUNC_PSEUDO;
    }
    /* "inference_type" hint for ORF (rare but real). */
    const char *itype = cJSON_GetStringValue(
        cJSON_GetObjectItem(allele, "inference_type"));
    if (status == FUNC_F && itype && strstr(itype, "ORF")) {
        status = FUNC_ORF;
    }

    /* Explicit anchor extraction (AIRR-C, 1-based on unaligned seq):
     *   - V: scan v_gene_delineations[] for the IMGT scheme; take
     *     `cdr3_start` (1-based first nt of CDR3 = first nt of conserved Cys).
     *     Convert to 0-based by subtracting 1.
     *   - J: top-level `j_cdr3_end` is the 1-based last nt of CDR3.
     *     The conserved Trp/Phe codon ends at j_cdr3_end, so the codon
     *     starts at j_cdr3_end - 2 (1-based) = j_cdr3_end - 3 (0-based).
     *     But we want the position OF the anchor codon, so 0-based start.
     */
    int explicit_anchor = -1;
    if (segment == SEG_V) {
        cJSON *delins = cJSON_GetObjectItem(allele, "v_gene_delineations");
        if (cJSON_IsArray(delins)) {
            cJSON *d;
            cJSON_ArrayForEach(d, delins) {
                const char *scheme = cJSON_GetStringValue(
                    cJSON_GetObjectItem(d, "delineation_scheme"));
                if (scheme && strcasecmp(scheme, "IMGT") == 0) {
                    cJSON *cs = cJSON_GetObjectItem(d, "cdr3_start");
                    if (cJSON_IsNumber(cs)) {
                        int v = (int)cs->valuedouble;
                        if (v > 0) explicit_anchor = v - 1;  /* 1-based → 0-based */
                    }
                    break;
                }
            }
        }
    } else if (segment == SEG_J) {
        cJSON *je = cJSON_GetObjectItem(allele, "j_cdr3_end");
        if (cJSON_IsNumber(je)) {
            int v = (int)je->valuedouble;
            /* j_cdr3_end is the 1-based last nt of CDR3 (i.e. last nt
             * of the conserved Trp/Phe codon). Anchor is the START of
             * that codon: (j_cdr3_end - 2) 1-based = (j_cdr3_end - 3)
             * 0-based. */
            if (v >= 3) explicit_anchor = v - 3;
        }
    }

    /* Stage the record. Grow the array if needed. */
    if (st->count >= st->cap) {
        int new_cap = st->cap ? st->cap * 2 : 16;
        AirrcStaged *grown = (AirrcStaged *)realloc(
            st->records, sizeof(AirrcStaged) * (size_t)new_cap);
        if (!grown) {
            if (err_msg) *err_msg = "out of memory growing AIRR-C record array";
            return -1;
        }
        st->records = grown;
        st->cap = new_cap;
    }

    AirrcStaged *r = &st->records[st->count];
    memset(r, 0, sizeof(*r));
    r->name = plain_dup(label);
    r->species = species ? plain_dup(species) : NULL;
    r->sequence = lower_dup(seq);
    r->sequence_length = r->sequence ? (int)strlen(r->sequence) : 0;
    r->locus = locus;
    r->segment = segment;
    r->functional_status = status;
    r->explicit_anchor = explicit_anchor;

    if (!r->name || !r->sequence) {
        if (err_msg) *err_msg = "out of memory copying AIRR-C strings";
        return -1;
    }
    st->count++;
    return 0;
}

/* ── Vtable ops ────────────────────────────────────────────────── */

static int airrc_next(ReferenceLoader *self,
                      LoadedAlleleRecord *out,
                      const char **err_msg) {
    AirrcState *st = (AirrcState *)self->state;
    if (!st) {
        if (err_msg) *err_msg = "loader has been closed";
        return -1;
    }
    /* T2-9: emit every record, regardless of status. The Python
     * policy layer (load_segment_alleles) decides whether to include
     * F/ORF/P/partial based on user flags. Loader is dumb pass-through. */
    while (st->cursor < st->count) {
        AirrcStaged *r = &st->records[st->cursor++];
        loaded_allele_record_init(out);
        out->name = r->name;
        out->n_aliases = 0;
        out->species = r->species;
        out->segment = r->segment;
        out->locus = r->locus;
        out->sequence = r->sequence;
        out->sequence_length = r->sequence_length;
        out->gapped_sequence = NULL;
        out->gapped_length = 0;
        out->gap_convention_imgt = false;
        out->functional_status = r->functional_status;
        out->explicit_anchor = r->explicit_anchor;
        out->source = "airrc-germline-set";
        return 0;
    }
    return 1;   /* EOF */
}

static void airrc_close(ReferenceLoader *self) {
    if (!self) return;
    AirrcState *st = (AirrcState *)self->state;
    if (st) {
        for (int i = 0; i < st->count; i++) {
            free(st->records[i].name);
            free(st->records[i].species);
            free(st->records[i].sequence);
        }
        free(st->records);
        free(st);
    }
    free(self);
}

static const ReferenceLoaderVTable AIRRC_VTABLE = {
    .next  = airrc_next,
    .close = airrc_close,
};

/* ── Factory ───────────────────────────────────────────────────── */

ReferenceLoader *airrc_germline_loader_open(const char *json_path,
                                            const char **err_msg) {
    if (!json_path) {
        if (err_msg) *err_msg = "json_path is NULL";
        return NULL;
    }
    char *raw = slurp_file(json_path, err_msg);
    if (!raw) return NULL;

    cJSON *root = cJSON_Parse(raw);
    free(raw);   /* cJSON copies what it needs */
    if (!root) {
        if (err_msg) *err_msg = "failed to parse AIRR-C JSON";
        return NULL;
    }

    AirrcState *st = (AirrcState *)calloc(1, sizeof(AirrcState));
    ReferenceLoader *loader = (ReferenceLoader *)calloc(1, sizeof(ReferenceLoader));
    if (!st || !loader) {
        cJSON_Delete(root);
        free(st);
        free(loader);
        if (err_msg) *err_msg = "out of memory";
        return NULL;
    }

    /* AIRR-C structure: top-level may be either:
     *   (a) { "GermlineSet": [ { allele_descriptions: [...] }, ... ] }
     *   (b) a bare GermlineSet object
     *   (c) a list of allele descriptions (older / minimal exports)
     *
     * We handle (a) and (b) directly; (c) is detected as "the root
     * IS the array". */
    cJSON *germline_array = cJSON_GetObjectItem(root, "GermlineSet");
    if (cJSON_IsArray(germline_array)) {
        cJSON *germline_set;
        cJSON_ArrayForEach(germline_set, germline_array) {
            const char *species = cJSON_GetStringValue(
                cJSON_GetObjectItem(
                    cJSON_GetObjectItem(germline_set, "species"),
                    "label"));
            cJSON *alleles = cJSON_GetObjectItem(germline_set,
                                                 "allele_descriptions");
            if (cJSON_IsArray(alleles)) {
                cJSON *allele;
                cJSON_ArrayForEach(allele, alleles) {
                    if (stage_allele(st, allele, species, err_msg) != 0) {
                        cJSON_Delete(root);
                        airrc_close(loader);
                        return NULL;
                    }
                }
            }
        }
    } else if (cJSON_IsObject(root)) {
        /* (b): root is a single GermlineSet. */
        const char *species = cJSON_GetStringValue(
            cJSON_GetObjectItem(
                cJSON_GetObjectItem(root, "species"), "label"));
        cJSON *alleles = cJSON_GetObjectItem(root, "allele_descriptions");
        if (cJSON_IsArray(alleles)) {
            cJSON *allele;
            cJSON_ArrayForEach(allele, alleles) {
                if (stage_allele(st, allele, species, err_msg) != 0) {
                    cJSON_Delete(root);
                    airrc_close(loader);
                    return NULL;
                }
            }
        }
    } else if (cJSON_IsArray(root)) {
        /* (c): bare list of alleles. */
        cJSON *allele;
        cJSON_ArrayForEach(allele, root) {
            if (stage_allele(st, allele, NULL, err_msg) != 0) {
                cJSON_Delete(root);
                airrc_close(loader);
                return NULL;
            }
        }
    }

    cJSON_Delete(root);

    loader->vt = &AIRRC_VTABLE;
    loader->state = st;
    return loader;
}
