/**
 * gdc_write.c — Write GdcData to a .gdc binary file.
 *
 * Layout: FileHeader → SectionTable → Section payloads (in ID order).
 * Each section is written to a contiguous block; the section table
 * records each section's offset and size for random access.
 */

#include "genairr/gdc_io.h"
#include "genairr/gdc_format.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ── Buffered write helpers ───────────────────────────────────── */

typedef struct {
    FILE   *fp;
    int     error;     /* sticky error flag */
} Writer;

static void w_bytes(Writer *w, const void *buf, size_t n) {
    if (w->error) return;
    if (fwrite(buf, 1, n, w->fp) != n) w->error = 1;
}

static void w_u8(Writer *w, uint8_t v)   { w_bytes(w, &v, 1); }
static void w_u16(Writer *w, uint16_t v) { w_bytes(w, &v, 2); }
static void w_u32(Writer *w, uint32_t v) { w_bytes(w, &v, 4); }
/* w_u64 reserved for future use with 64-bit section offsets. */
static void w_i16(Writer *w, int16_t v)  { w_bytes(w, &v, 2); }
static void w_f64(Writer *w, double v)   { w_bytes(w, &v, 8); }

static void w_str(Writer *w, const char *s) {
    uint16_t len = (uint16_t)strlen(s);
    w_u16(w, len);
    w_bytes(w, s, len);
}

static void w_f64_array(Writer *w, const double *arr, int n) {
    w_bytes(w, arr, (size_t)n * sizeof(double));
}

/* ── Section writers ──────────────────────────────────────────── */

static void write_metadata(Writer *w, const GdcData *data) {
    w_str(w, data->name);
    w_u8(w, data->chain_type);
    w_u8(w, data->has_d ? 1 : 0);
    w_str(w, data->species);
    w_str(w, data->reference_set);
    w_u32(w, data->build_date);
}

static void write_allele_pool(Writer *w, const AllelePool *pool) {
    w_u16(w, (uint16_t)pool->count);
    for (int i = 0; i < pool->count; i++) {
        const Allele *a = &pool->alleles[i];
        w_str(w, a->name);
        w_str(w, a->gene);
        w_str(w, a->family);
        w_u16(w, a->length);
        w_bytes(w, a->seq, a->length);
        w_i16(w, (int16_t)a->anchor);
    }
}

static void write_alleles(Writer *w, const GdcData *data) {
    write_allele_pool(w, &data->v_alleles);
    write_allele_pool(w, &data->d_alleles);
    write_allele_pool(w, &data->j_alleles);
    write_allele_pool(w, &data->c_alleles);
}

static void write_gene_use(Writer *w, const GdcData *data) {
    /* Count how many segment types have entries */
    uint8_t n_types = 0;
    for (int i = 0; i < 4; i++) {
        if (data->gene_use[i].count > 0) n_types++;
    }
    w_u8(w, n_types);

    for (int seg = 0; seg < 4; seg++) {
        const GeneUseTable *t = &data->gene_use[seg];
        if (t->count == 0) continue;

        w_u8(w, (uint8_t)seg);
        w_u32(w, (uint32_t)t->count);
        for (int i = 0; i < t->count; i++) {
            w_str(w, t->entries[i].gene_name);
            w_f64(w, t->entries[i].probability);
        }
    }
}

static void write_trim_dists(Writer *w, const GdcData *data) {
    uint8_t n_types = 0;
    for (int i = 0; i < 4; i++) {
        if (data->trim_dists[i].count > 0) n_types++;
    }
    w_u8(w, n_types);

    for (int tt = 0; tt < 4; tt++) {
        const TrimDistTable *t = &data->trim_dists[tt];
        if (t->count == 0) continue;

        w_u8(w, (uint8_t)tt);
        w_u16(w, (uint16_t)t->count);
        for (int i = 0; i < t->count; i++) {
            const GeneTrimDist *d = &t->dists[i];
            w_str(w, d->family_name);
            w_str(w, d->gene_name);
            w_u16(w, (uint16_t)d->max_trim);
            w_f64_array(w, d->probs, d->max_trim + 1);
        }
    }
}

static void write_np_params(Writer *w, const GdcData *data) {
    w_u8(w, (uint8_t)data->n_np_regions);

    for (int r = 0; r < data->n_np_regions; r++) {
        const NpParams *np = &data->np[r];

        /* Length distribution */
        w_u16(w, (uint16_t)np->max_length);
        w_f64_array(w, np->length_probs, np->max_length + 1);

        /* First base probabilities */
        w_f64_array(w, np->first_base, 4);

        /* Markov transitions: n_positions × 16 doubles */
        w_u16(w, (uint16_t)np->n_positions);
        w_f64_array(w, np->transitions, np->n_positions * 16);
    }
}

static void write_p_nuc_probs(Writer *w, const GdcData *data) {
    w_u8(w, (uint8_t)data->p_nuc_max);
    w_f64_array(w, data->p_nuc_probs, data->p_nuc_max + 1);
}

static void write_mutation_model(Writer *w, const GdcData *data) {
    w_u8(w, data->mutation_model_type);

    if (data->mutation_model_type == GDC_MUTMODEL_S5F) {
        /* Mutability: 3125 doubles */
        w_f64_array(w, data->mutability, GDC_S5F_KMER_SPACE);

        /* Substitution: per 5-mer, count + bases + weights */
        for (int k = 0; k < GDC_S5F_KMER_SPACE; k++) {
            uint8_t cnt = data->sub_counts[k];
            w_u8(w, cnt);
            w_bytes(w, &data->sub_bases[k * 4], cnt);
            w_f64_array(w, &data->substitution[k * 4], cnt);
        }
    }
}

/* ── Main write function ──────────────────────────────────────── */

int gdc_save(const char *path, const GdcData *data) {
    FILE *fp = fopen(path, "wb");
    if (!fp) return -1;

    Writer w = { .fp = fp, .error = 0 };

    /* We need to know section offsets before writing the header,
     * so we write sections to a temporary buffer approach:
     * 1. Write header + section table placeholders
     * 2. Write each section, recording offset/size
     * 3. Seek back and patch the section table */

    int n_sections = GDC_SECTION_COUNT;

    /* Write file header */
    GdcFileHeader hdr;
    memset(&hdr, 0, sizeof(hdr));
    hdr.magic[0] = GDC_MAGIC_0;
    hdr.magic[1] = GDC_MAGIC_1;
    hdr.magic[2] = GDC_MAGIC_2;
    hdr.magic[3] = GDC_MAGIC_3;
    hdr.format_version = GDC_FORMAT_VERSION;
    hdr.n_sections = (uint16_t)n_sections;
    w_bytes(&w, &hdr, sizeof(hdr));

    /* Write section table placeholders (will patch later) */
    long toc_pos = ftell(fp);
    GdcSectionEntry entries[GDC_SECTION_COUNT];
    memset(entries, 0, sizeof(entries));
    w_bytes(&w, entries, sizeof(entries));

    /* Write each section, record offset and size */
    typedef void (*SectionWriter)(Writer *, const GdcData *);
    static const SectionWriter writers[] = {
        write_metadata,
        write_alleles,
        write_gene_use,
        write_trim_dists,
        write_np_params,
        write_p_nuc_probs,
        write_mutation_model,
    };

    for (int i = 0; i < n_sections; i++) {
        long start = ftell(fp);
        entries[i].section_id = (uint32_t)i;
        entries[i].offset = (uint64_t)start;
        writers[i](&w, data);
        long end = ftell(fp);
        entries[i].size = (uint32_t)(end - start);
    }

    /* Patch section table */
    fseek(fp, toc_pos, SEEK_SET);
    w_bytes(&w, entries, sizeof(entries));

    fclose(fp);
    return w.error ? -1 : 0;
}
