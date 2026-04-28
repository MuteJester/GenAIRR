/**
 * loader_imgt.c — IMGT V-QUEST FASTA reader as a ReferenceLoader.
 *
 * Parses IMGT V-QUEST reference FASTAs (gapped sequences, pipe-delimited
 * headers) and yields LoadedAlleleRecord values one at a time.
 *
 * Header format (real example):
 *   >X62106|IGHV1-2*02|Homo sapiens|F|V-REGION|...|296 nt|1| | | | |296+0=296| |
 *
 * Pipe fields (1-based):
 *   1   accession (e.g. X62106)
 *   2   allele name (e.g. IGHV1-2*02)
 *   3   species (e.g. Homo sapiens)
 *   4   functional status (F / ORF / P), may also contain "partial"
 *   5   feature label (V-REGION, J-REGION, etc.)
 *   6+  various coordinates / metadata we mostly ignore
 *
 * Sequence: lowercase, '.' = gap (IMGT-aligned), continues across
 * subsequent lines until the next '>' or EOF.
 *
 * Filtering policy (Phase 1A):
 *   - "partial" anywhere in the header → skipped silently
 *   - functional status not in {F, ORF} → skipped silently (pseudogenes
 *     deferred to Phase 3 with explicit keep_pseudogenes flag)
 *
 * State machine — the inner loop runs over three states:
 *
 *   AWAIT_HEADER       — looking for the next '>' line
 *   AWAIT_SEQUENCE     — header parsed, accumulating sequence
 *   READY_TO_EMIT      — virtual state; next() returns here
 *
 * When sequence accumulation hits another '>' line, the new header is
 * stashed in `next_header_*` slots so the following next() call picks
 * it up without rewinding the FILE*.
 */

#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "genairr/ref_loader.h"
#include "genairr/ref_record.h"

#define LINE_INITIAL_CAP   512
#define NAME_BUF_LEN       128
#define SPECIES_BUF_LEN    128

typedef struct {
    char             name[NAME_BUF_LEN];
    char             species[SPECIES_BUF_LEN];
    FunctionalStatus status;
    bool             partial;
    bool             populated;
} HeaderSlot;

typedef struct {
    FILE       *fp;

    char       *line;
    size_t      line_cap;

    char       *gapped;
    int         gapped_len;
    int         gapped_cap;
    char       *ungapped;
    int         ungapped_len;
    int         ungapped_cap;

    /* Header for the record currently being built / about to emit.
     * Stays valid across the next() return so the LoadedAlleleRecord's
     * pointers don't dangle. */
    HeaderSlot  cur;

    /* Header read AHEAD when we encountered a '>' while accumulating
     * the previous record's sequence. The next call to next() picks
     * this up directly without re-reading from the FILE*. */
    HeaderSlot  next_header;

    Segment     segment_hint;
    bool        eof;
} ImgtState;

/* ── Buffer helpers ────────────────────────────────────────────── */

static int ensure_cap(char **buf, int *cap, int need) {
    if (need <= *cap) return 0;
    int new_cap = *cap ? *cap : 64;
    while (new_cap < need) new_cap *= 2;
    char *grown = (char *)realloc(*buf, (size_t)new_cap);
    if (!grown) return -1;
    *buf = grown;
    *cap = new_cap;
    return 0;
}

/* getline replacement (POSIX getline isn't on every toolchain we
 * target — older Windows in particular). Returns:
 *   0  success — *buf NUL-terminated, trailing \r\n stripped
 *   1  EOF before any char read
 *  -1  error
 */
static int read_line(FILE *fp, char **buf, size_t *cap) {
    if (*cap < LINE_INITIAL_CAP) {
        char *grown = (char *)realloc(*buf, LINE_INITIAL_CAP);
        if (!grown) return -1;
        *buf = grown;
        *cap = LINE_INITIAL_CAP;
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

/* ── Header parsing ────────────────────────────────────────────── */

static int extract_pipe_field(const char *header, int idx,
                              char *dst, int dst_size) {
    int field = 1;
    const char *start = header;
    if (*start == '>') start++;
    while (*start && field < idx) {
        if (*start == '|') field++;
        start++;
    }
    if (field != idx) {
        if (dst && dst_size > 0) dst[0] = '\0';
        return 0;
    }
    int n = 0;
    while (start[n] && start[n] != '|' && n < dst_size - 1) {
        dst[n] = start[n];
        n++;
    }
    dst[n] = '\0';
    return n;
}

static FunctionalStatus parse_status(const char *s) {
    if (!s) return FUNC_UNKNOWN;
    char buf[16];
    int n = 0;
    for (int i = 0; s[i] && n + 1 < (int)sizeof(buf); i++) {
        if (s[i] == '(' || s[i] == ')' ||
            s[i] == '[' || s[i] == ']' ||
            isspace((unsigned char)s[i])) continue;
        buf[n++] = s[i];
    }
    buf[n] = '\0';
    if (!strcmp(buf, "F"))   return FUNC_F;
    if (!strcmp(buf, "ORF")) return FUNC_ORF;
    if (!strcmp(buf, "P"))   return FUNC_PSEUDO;
    return FUNC_UNKNOWN;
}

static void parse_header_into(HeaderSlot *slot, const char *header) {
    /* Field 2: allele name (mandatory; fall back to whole-line on parse failure). */
    if (extract_pipe_field(header, 2, slot->name, NAME_BUF_LEN) == 0) {
        const char *p = header[0] == '>' ? header + 1 : header;
        strncpy(slot->name, p, NAME_BUF_LEN - 1);
        slot->name[NAME_BUF_LEN - 1] = '\0';
    }
    extract_pipe_field(header, 3, slot->species, SPECIES_BUF_LEN);

    char status_buf[16];
    extract_pipe_field(header, 4, status_buf, (int)sizeof(status_buf));
    slot->status  = parse_status(status_buf);
    slot->partial = (strstr(header, "partial") != NULL);
    /* T2-9: a record marked "partial" overrides whatever the IMGT
     * status field said — partial sequences cannot represent a
     * complete F or ORF allele. Tag them as FUNC_PARTIAL so the
     * Python policy layer can filter them out (default) or
     * include them (opt-in). Pseudogenes (P) flow through with
     * FUNC_PSEUDO so the policy layer can decide. */
    if (slot->partial) {
        slot->status = FUNC_PARTIAL;
    }
    slot->populated = true;
}

/* ── Sequence accumulation ─────────────────────────────────────── */

static int append_seq_chars(ImgtState *st, const char *src) {
    int src_len = (int)strlen(src);
    if (ensure_cap(&st->gapped, &st->gapped_cap,
                   st->gapped_len + src_len + 1) != 0) return -1;
    if (ensure_cap(&st->ungapped, &st->ungapped_cap,
                   st->ungapped_len + src_len + 1) != 0) return -1;
    for (int i = 0; src[i]; i++) {
        char c = (char)tolower((unsigned char)src[i]);
        if (isspace((unsigned char)c)) continue;
        st->gapped[st->gapped_len++] = c;
        if (c != '.') st->ungapped[st->ungapped_len++] = c;
    }
    st->gapped[st->gapped_len] = '\0';
    st->ungapped[st->ungapped_len] = '\0';
    return 0;
}

static void reset_sequence(ImgtState *st) {
    st->gapped_len = 0;
    st->ungapped_len = 0;
    if (st->gapped) st->gapped[0] = '\0';
    if (st->ungapped) st->ungapped[0] = '\0';
}

/* ── Record emission ───────────────────────────────────────────── */

static void emit_record(ImgtState *st, LoadedAlleleRecord *out) {
    loaded_allele_record_init(out);
    out->name    = st->cur.name;
    out->species = st->cur.species[0] ? st->cur.species : NULL;
    out->n_aliases = 0;

    Segment seg = segment_from_gene_name(st->cur.name);
    if (seg == SEG_UNKNOWN) seg = st->segment_hint;
    out->segment = seg;
    out->locus = locus_from_gene_name(st->cur.name);

    out->sequence        = st->ungapped;
    out->sequence_length = st->ungapped_len;
    out->gapped_sequence = st->gapped;
    out->gapped_length   = st->gapped_len;
    out->gap_convention_imgt = true;

    out->functional_status = st->cur.status;
    out->explicit_anchor   = -1;
    out->source            = "imgt-vquest";
}

/* ── Vtable ops ────────────────────────────────────────────────── */

/* Read until we have either a kept record ready to emit, or EOF.
 * On entry we either have:
 *   (a) cur.populated=false, next_header.populated=false → must read header
 *   (b) cur.populated=false, next_header.populated=true  → promote
 *   (c) (cur populated mid-stream cannot occur here — only inside the
 *        inner loop, and we always reset at function entry)
 */
static int imgt_next(ReferenceLoader *self,
                     LoadedAlleleRecord *out,
                     const char **err_msg) {
    ImgtState *st = (ImgtState *)self->state;
    if (!st || !st->fp) {
        if (err_msg) *err_msg = "loader has been closed";
        return -1;
    }

    while (1) {
        /* Step 1: ensure cur is populated. */
        if (!st->cur.populated) {
            if (st->next_header.populated) {
                /* Promote pre-read header. */
                st->cur = st->next_header;
                st->next_header.populated = false;
            } else {
                /* Scan forward for the next '>' line. */
                if (st->eof) return 1;
                int rc;
                do {
                    rc = read_line(st->fp, &st->line, &st->line_cap);
                    if (rc < 0) {
                        if (err_msg) *err_msg = "I/O error reading IMGT FASTA";
                        return -1;
                    }
                    if (rc == 1) { st->eof = true; return 1; }
                } while (st->line[0] != '>');
                parse_header_into(&st->cur, st->line);
            }
        }

        /* Step 2: accumulate sequence lines until a new '>' or EOF. */
        reset_sequence(st);
        bool hit_new_header = false;
        while (!hit_new_header) {
            int rc = read_line(st->fp, &st->line, &st->line_cap);
            if (rc < 0) {
                if (err_msg) *err_msg = "I/O error reading IMGT FASTA";
                return -1;
            }
            if (rc == 1) { st->eof = true; break; }
            if (st->line[0] == '>') {
                /* End of current record; stash header for next iteration. */
                parse_header_into(&st->next_header, st->line);
                hit_new_header = true;
                break;
            }
            if (append_seq_chars(st, st->line) != 0) {
                if (err_msg) *err_msg = "out of memory accumulating sequence";
                return -1;
            }
        }

        /* Step 3: emit. T2-9: we no longer drop records by status —
         * every parsed record flows through, with `functional_status`
         * tagged so the Python policy layer (load_segment_alleles)
         * can decide whether to include F/ORF/P/partial. The only
         * loader-side filter that survives is the empty-sequence
         * check (a header with no associated sequence is a parse
         * artifact, not a biological signal). */
        if (st->ungapped_len == 0) {
            st->cur.populated = false;
            continue;
        }
        emit_record(st, out);
        /* Mark cur as consumed; next call will promote next_header (if
         * any) or advance. */
        st->cur.populated = false;
        return 0;
    }
}

static void imgt_close(ReferenceLoader *self) {
    if (!self) return;
    ImgtState *st = (ImgtState *)self->state;
    if (st) {
        if (st->fp) fclose(st->fp);
        free(st->line);
        free(st->gapped);
        free(st->ungapped);
        free(st);
    }
    free(self);
}

static const ReferenceLoaderVTable IMGT_VTABLE = {
    .next  = imgt_next,
    .close = imgt_close,
};

/* ── Factory ───────────────────────────────────────────────────── */

ReferenceLoader *imgt_vquest_loader_open(const char *fasta_path,
                                         Segment segment_hint,
                                         const char **err_msg) {
    if (!fasta_path) {
        if (err_msg) *err_msg = "fasta_path is NULL";
        return NULL;
    }
    FILE *fp = fopen(fasta_path, "r");
    if (!fp) {
        if (err_msg) *err_msg = "could not open FASTA file";
        return NULL;
    }
    ImgtState *st = (ImgtState *)calloc(1, sizeof(ImgtState));
    ReferenceLoader *loader = (ReferenceLoader *)calloc(1, sizeof(ReferenceLoader));
    if (!st || !loader) {
        fclose(fp);
        free(st);
        free(loader);
        if (err_msg) *err_msg = "out of memory";
        return NULL;
    }
    st->fp = fp;
    st->segment_hint = segment_hint;
    st->cur.status = FUNC_UNKNOWN;
    st->next_header.status = FUNC_UNKNOWN;

    loader->vt = &IMGT_VTABLE;
    loader->state = st;
    return loader;
}
