/**
 * loader_plain.c — bare FASTA reader as a ReferenceLoader.
 *
 * The minimum-input format for custom lab references and exotic
 * species: a FASTA file with one allele per record, header is just
 * the allele name (no pipe-delimited metadata), sequence is ungapped.
 * No functional-status filter (everything is kept), no anchor
 * coordinate (must come from the AnchorResolver's motif fallback).
 *
 * The caller MUST provide a locus hint at open time — without it the
 * J anchor scan can't decide between Trp (IGH/TRB/TRD) and Phe
 * (IGK/IGL/TRA/TRG) confidence. A `locus_hint = LOCUS_UNKNOWN` is
 * accepted (callers may not always know the locus); resolution
 * downgrades to CONF_BEST_GUESS in that case.
 *
 * Design contract is identical to loader_imgt.c:
 *   - all string fields owned by the loader, valid until next
 *     next() call
 *   - close() is NULL-safe and idempotent
 *   - ungapped sequence is lowercased
 *
 * The state machine uses the same `cur` + `next_header` slot pair
 * as the IMGT loader so the inner-loop logic stays simple.
 */

#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "genairr/ref_loader.h"
#include "genairr/ref_record.h"

#define PLAIN_LINE_INITIAL_CAP   512
#define PLAIN_NAME_BUF_LEN       128

typedef struct {
    char  name[PLAIN_NAME_BUF_LEN];
    bool  populated;
} PlainHeaderSlot;

typedef struct {
    FILE   *fp;

    char   *line;
    size_t  line_cap;

    char   *seq;
    int     seq_len;
    int     seq_cap;

    PlainHeaderSlot cur;
    PlainHeaderSlot next_header;

    Locus   locus_hint;
    Segment segment_hint;

    bool    eof;
} PlainState;

/* ── Buffer helpers (mirror loader_imgt.c) ─────────────────────── */

static int plain_ensure_cap(char **buf, int *cap, int need) {
    if (need <= *cap) return 0;
    int new_cap = *cap ? *cap : 64;
    while (new_cap < need) new_cap *= 2;
    char *grown = (char *)realloc(*buf, (size_t)new_cap);
    if (!grown) return -1;
    *buf = grown;
    *cap = new_cap;
    return 0;
}

static int plain_read_line(FILE *fp, char **buf, size_t *cap) {
    if (*cap < PLAIN_LINE_INITIAL_CAP) {
        char *grown = (char *)realloc(*buf, PLAIN_LINE_INITIAL_CAP);
        if (!grown) return -1;
        *buf = grown;
        *cap = PLAIN_LINE_INITIAL_CAP;
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

/* ── Header parsing — name only ────────────────────────────────── */

static void plain_parse_header(PlainHeaderSlot *slot, const char *line) {
    /* Strip the leading '>' and any leading whitespace. The "name" is
     * everything up to the first whitespace character (so a description
     * suffix like ">IGHV1-2*01 Homo sapiens" still picks the name).
     * If the line is just '>', name becomes empty string. */
    const char *p = line[0] == '>' ? line + 1 : line;
    while (*p && isspace((unsigned char)*p)) p++;

    int n = 0;
    while (*p && !isspace((unsigned char)*p) && n < PLAIN_NAME_BUF_LEN - 1) {
        slot->name[n++] = *p++;
    }
    slot->name[n] = '\0';
    slot->populated = true;
}

/* ── Sequence accumulation ─────────────────────────────────────── */

static int plain_append_seq(PlainState *st, const char *src) {
    int src_len = (int)strlen(src);
    if (plain_ensure_cap(&st->seq, &st->seq_cap,
                         st->seq_len + src_len + 1) != 0) return -1;
    for (int i = 0; src[i]; i++) {
        char c = (char)tolower((unsigned char)src[i]);
        if (isspace((unsigned char)c)) continue;
        /* Plain FASTA may contain '.' or '-' (unlikely but tolerate) —
         * strip them silently. */
        if (c == '.' || c == '-') continue;
        st->seq[st->seq_len++] = c;
    }
    st->seq[st->seq_len] = '\0';
    return 0;
}

static void plain_reset_seq(PlainState *st) {
    st->seq_len = 0;
    if (st->seq) st->seq[0] = '\0';
}

/* ── Record emission ───────────────────────────────────────────── */

static void plain_emit_record(PlainState *st, LoadedAlleleRecord *out) {
    loaded_allele_record_init(out);
    out->name = st->cur.name;
    out->n_aliases = 0;
    out->species = NULL;

    /* Try to infer segment from name; fall back to caller's hint.
     * Locus from name first, fall back to caller's hint. */
    Segment seg = segment_from_gene_name(st->cur.name);
    if (seg == SEG_UNKNOWN) seg = st->segment_hint;
    Locus loc = locus_from_gene_name(st->cur.name);
    if (loc == LOCUS_UNKNOWN) loc = st->locus_hint;

    out->segment = seg;
    out->locus = loc;
    out->sequence = st->seq;
    out->sequence_length = st->seq_len;
    out->gapped_sequence = NULL;
    out->gapped_length = 0;
    out->gap_convention_imgt = false;
    out->functional_status = FUNC_UNKNOWN;
    out->explicit_anchor = -1;
    out->source = "plain-fasta";
}

/* ── Vtable ops ────────────────────────────────────────────────── */

static int plain_next(ReferenceLoader *self,
                      LoadedAlleleRecord *out,
                      const char **err_msg) {
    PlainState *st = (PlainState *)self->state;
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
                    rc = plain_read_line(st->fp, &st->line, &st->line_cap);
                    if (rc < 0) {
                        if (err_msg) *err_msg = "I/O error reading FASTA";
                        return -1;
                    }
                    if (rc == 1) { st->eof = true; return 1; }
                } while (st->line[0] != '>');
                plain_parse_header(&st->cur, st->line);
            }
        }

        plain_reset_seq(st);
        bool hit_new_header = false;
        while (!hit_new_header) {
            int rc = plain_read_line(st->fp, &st->line, &st->line_cap);
            if (rc < 0) {
                if (err_msg) *err_msg = "I/O error reading FASTA";
                return -1;
            }
            if (rc == 1) { st->eof = true; break; }
            if (st->line[0] == '>') {
                plain_parse_header(&st->next_header, st->line);
                hit_new_header = true;
                break;
            }
            if (plain_append_seq(st, st->line) != 0) {
                if (err_msg) *err_msg = "out of memory";
                return -1;
            }
        }

        /* Skip empty records (header without sequence) — but don't
         * filter on any other criterion. Plain FASTA may contain
         * pseudogenes, ORFs, novel alleles — caller decides. */
        if (st->seq_len == 0) {
            st->cur.populated = false;
            continue;
        }

        plain_emit_record(st, out);
        st->cur.populated = false;
        return 0;
    }
}

static void plain_close(ReferenceLoader *self) {
    if (!self) return;
    PlainState *st = (PlainState *)self->state;
    if (st) {
        if (st->fp) fclose(st->fp);
        free(st->line);
        free(st->seq);
        free(st);
    }
    free(self);
}

static const ReferenceLoaderVTable PLAIN_VTABLE = {
    .next  = plain_next,
    .close = plain_close,
};

/* ── Factory ───────────────────────────────────────────────────── */

ReferenceLoader *plain_fasta_loader_open(const char *fasta_path,
                                         Locus locus_hint,
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
    PlainState *st = (PlainState *)calloc(1, sizeof(PlainState));
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

    loader->vt = &PLAIN_VTABLE;
    loader->state = st;
    return loader;
}
