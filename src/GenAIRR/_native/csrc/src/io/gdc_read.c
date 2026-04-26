/**
 * gdc_read.c — Read a .gdc binary file into GdcData.
 *
 * Validates magic/version, reads section table, then processes
 * each known section. Unknown sections are silently skipped
 * (forward compatibility).
 */

#include "genairr/gdc_io.h"
#include "genairr/gdc_format.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ── Read helpers ─────────────────────────────────────────────── */

typedef struct {
    FILE *fp;
    int   error;    /* sticky: 0 = ok, -1 = I/O, -2 = format */
} Reader;

static void r_bytes(Reader *r, void *buf, size_t n) {
    if (r->error) return;
    if (fread(buf, 1, n, r->fp) != n) r->error = -1;
}

static uint8_t r_u8(Reader *r) {
    uint8_t v = 0;
    r_bytes(r, &v, 1);
    return v;
}

static uint16_t r_u16(Reader *r) {
    uint16_t v = 0;
    r_bytes(r, &v, 2);
    return v;
}

static uint32_t r_u32(Reader *r) {
    uint32_t v = 0;
    r_bytes(r, &v, 4);
    return v;
}

static int16_t r_i16(Reader *r) {
    int16_t v = 0;
    r_bytes(r, &v, 2);
    return v;
}

static double r_f64(Reader *r) {
    double v = 0.0;
    r_bytes(r, &v, 8);
    return v;
}

static void r_f64_array(Reader *r, double *arr, int n) {
    r_bytes(r, arr, (size_t)n * sizeof(double));
}

/** Read a length-prefixed string into a fixed buffer. */
static void r_str(Reader *r, char *buf, int buf_size) {
    uint16_t len = r_u16(r);
    if (r->error) return;

    if (len >= buf_size) {
        /* String too long for buffer — read and truncate */
        r_bytes(r, buf, (size_t)(buf_size - 1));
        buf[buf_size - 1] = '\0';
        /* Skip remaining bytes */
        int skip = len - (buf_size - 1);
        if (skip > 0) fseek(r->fp, skip, SEEK_CUR);
    } else {
        r_bytes(r, buf, len);
        buf[len] = '\0';
    }
}

/* ── Section readers ──────────────────────────────────────────── */

static void read_metadata(Reader *r, GdcData *data) {
    r_str(r, data->name, (int)sizeof(data->name));
    data->chain_type = r_u8(r);
    data->has_d = r_u8(r) != 0;
    r_str(r, data->species, (int)sizeof(data->species));
    r_str(r, data->reference_set, (int)sizeof(data->reference_set));
    data->build_date = r_u32(r);
}

static void read_allele_pool(Reader *r, AllelePool *pool) {
    uint16_t count = r_u16(r);
    if (r->error) return;

    *pool = allele_pool_create((int)count > 0 ? (int)count : 1);

    for (int i = 0; i < (int)count && !r->error; i++) {
        Allele a;
        memset(&a, 0, sizeof(a));

        r_str(r, a.name, GENAIRR_MAX_ALLELE_NAME);
        r_str(r, a.gene, GENAIRR_MAX_ALLELE_NAME);
        r_str(r, a.family, GENAIRR_MAX_ALLELE_NAME);

        a.length = r_u16(r);
        if (a.length > GENAIRR_MAX_ALLELE_SEQ - 1) {
            r->error = -2;
            return;
        }
        r_bytes(r, a.seq, a.length);
        a.seq[a.length] = '\0';

        /* Preserve the signed sentinel: GDC writes -1 for anchorless
         * V/J alleles. Pre-T0-7 this was cast to uint16_t, turning -1
         * into 65535 and triggering a 65534-base J trim overflow. */
        a.anchor = r_i16(r);

        allele_pool_add(pool, &a);
    }
}

static void read_alleles(Reader *r, GdcData *data) {
    read_allele_pool(r, &data->v_alleles);
    if (r->error) return;

    /* Set segment types */
    for (int i = 0; i < data->v_alleles.count; i++)
        data->v_alleles.alleles[i].segment_type = SEG_V;

    read_allele_pool(r, &data->d_alleles);
    if (r->error) return;
    for (int i = 0; i < data->d_alleles.count; i++)
        data->d_alleles.alleles[i].segment_type = SEG_D;

    read_allele_pool(r, &data->j_alleles);
    if (r->error) return;
    for (int i = 0; i < data->j_alleles.count; i++)
        data->j_alleles.alleles[i].segment_type = SEG_J;

    read_allele_pool(r, &data->c_alleles);
    if (r->error) return;
    for (int i = 0; i < data->c_alleles.count; i++)
        data->c_alleles.alleles[i].segment_type = SEG_C;
}

static void gene_use_table_add(GeneUseTable *t, const char *name, double prob) {
    if (t->count >= t->capacity) {
        int new_cap = t->capacity ? t->capacity * 2 : 16;
        GeneUseEntry *new_entries = realloc(t->entries,
                                            (size_t)new_cap * sizeof(GeneUseEntry));
        if (!new_entries) return;
        t->entries = new_entries;
        t->capacity = new_cap;
    }
    strncpy(t->entries[t->count].gene_name, name, GENAIRR_MAX_ALLELE_NAME - 1);
    t->entries[t->count].gene_name[GENAIRR_MAX_ALLELE_NAME - 1] = '\0';
    t->entries[t->count].probability = prob;
    t->count++;
}

static void read_gene_use(Reader *r, GdcData *data) {
    uint8_t n_types = r_u8(r);

    for (int t = 0; t < (int)n_types && !r->error; t++) {
        uint8_t seg = r_u8(r);
        uint32_t count = r_u32(r);

        if (seg >= 4) { r->error = -2; return; }

        for (uint32_t i = 0; i < count && !r->error; i++) {
            char name[GENAIRR_MAX_ALLELE_NAME];
            r_str(r, name, GENAIRR_MAX_ALLELE_NAME);
            double prob = r_f64(r);
            gene_use_table_add(&data->gene_use[seg], name, prob);
        }
    }
}

static void trim_dist_table_add(TrimDistTable *t, const char *family,
                                 const char *gene, double *probs, int max_trim) {
    if (t->count >= t->capacity) {
        int new_cap = t->capacity ? t->capacity * 2 : 16;
        GeneTrimDist *new_dists = realloc(t->dists,
                                          (size_t)new_cap * sizeof(GeneTrimDist));
        if (!new_dists) { free(probs); return; }
        t->dists = new_dists;
        t->capacity = new_cap;
    }
    GeneTrimDist *d = &t->dists[t->count];
    strncpy(d->family_name, family, GENAIRR_MAX_ALLELE_NAME - 1);
    d->family_name[GENAIRR_MAX_ALLELE_NAME - 1] = '\0';
    strncpy(d->gene_name, gene, GENAIRR_MAX_ALLELE_NAME - 1);
    d->gene_name[GENAIRR_MAX_ALLELE_NAME - 1] = '\0';
    d->probs = probs;
    d->max_trim = max_trim;
    t->count++;
}

static void read_trim_dists(Reader *r, GdcData *data) {
    uint8_t n_types = r_u8(r);

    for (int t = 0; t < (int)n_types && !r->error; t++) {
        uint8_t tt = r_u8(r);
        uint16_t count = r_u16(r);

        if (tt >= 4) { r->error = -2; return; }

        for (int i = 0; i < (int)count && !r->error; i++) {
            char family[GENAIRR_MAX_ALLELE_NAME];
            char gene[GENAIRR_MAX_ALLELE_NAME];
            r_str(r, family, GENAIRR_MAX_ALLELE_NAME);
            r_str(r, gene, GENAIRR_MAX_ALLELE_NAME);

            uint16_t max_trim = r_u16(r);
            int n_probs = (int)max_trim + 1;
            double *probs = calloc((size_t)n_probs, sizeof(double));
            if (!probs) { r->error = -1; return; }
            r_f64_array(r, probs, n_probs);

            trim_dist_table_add(&data->trim_dists[tt], family, gene,
                                probs, (int)max_trim);
        }
    }
}

static void read_np_params(Reader *r, GdcData *data) {
    data->n_np_regions = r_u8(r);
    if (data->n_np_regions > 2) { r->error = -2; return; }

    for (int i = 0; i < data->n_np_regions && !r->error; i++) {
        NpParams *np = &data->np[i];

        /* Length distribution */
        np->max_length = r_u16(r);
        int n_probs = np->max_length + 1;
        np->length_probs = calloc((size_t)n_probs, sizeof(double));
        if (!np->length_probs) { r->error = -1; return; }
        r_f64_array(r, np->length_probs, n_probs);

        /* First base probabilities */
        r_f64_array(r, np->first_base, 4);

        /* Markov transitions */
        np->n_positions = r_u16(r);
        int n_trans = np->n_positions * 16;
        np->transitions = calloc((size_t)n_trans, sizeof(double));
        if (!np->transitions) { r->error = -1; return; }
        r_f64_array(r, np->transitions, n_trans);
    }
}

static void read_p_nuc_probs(Reader *r, GdcData *data) {
    data->p_nuc_max = r_u8(r);
    int n_probs = data->p_nuc_max + 1;
    data->p_nuc_probs = calloc((size_t)n_probs, sizeof(double));
    if (!data->p_nuc_probs) { r->error = -1; return; }
    r_f64_array(r, data->p_nuc_probs, n_probs);
}

static void read_mutation_model(Reader *r, GdcData *data) {
    data->mutation_model_type = r_u8(r);

    if (data->mutation_model_type == GDC_MUTMODEL_S5F) {
        /* Mutability: 3125 doubles */
        r_f64_array(r, data->mutability, GDC_S5F_KMER_SPACE);

        /* Substitution: per 5-mer, count + bases + weights */
        for (int k = 0; k < GDC_S5F_KMER_SPACE && !r->error; k++) {
            uint8_t cnt = r_u8(r);
            data->sub_counts[k] = cnt;
            if (cnt > 4) { r->error = -2; return; }
            r_bytes(r, &data->sub_bases[k * 4], cnt);
            r_f64_array(r, &data->substitution[k * 4], cnt);
        }
    }
}

/* ── Lifecycle ────────────────────────────────────────────────── */

void gdc_data_init(GdcData *data) {
    memset(data, 0, sizeof(*data));
}

void gdc_data_destroy(GdcData *data) {
    allele_pool_destroy(&data->v_alleles);
    allele_pool_destroy(&data->d_alleles);
    allele_pool_destroy(&data->j_alleles);
    allele_pool_destroy(&data->c_alleles);

    for (int i = 0; i < 4; i++) {
        free(data->gene_use[i].entries);
        for (int j = 0; j < data->trim_dists[i].count; j++) {
            free(data->trim_dists[i].dists[j].probs);
        }
        free(data->trim_dists[i].dists);
    }

    for (int i = 0; i < 2; i++) {
        free(data->np[i].length_probs);
        free(data->np[i].transitions);
    }

    free(data->p_nuc_probs);
    memset(data, 0, sizeof(*data));
}

/* ── Main load function ───────────────────────────────────────── */

int gdc_load(const char *path, GdcData *data) {
    FILE *fp = fopen(path, "rb");
    if (!fp) return -1;

    gdc_data_init(data);
    Reader r = { .fp = fp, .error = 0 };

    /* Read and validate header */
    GdcFileHeader hdr;
    r_bytes(&r, &hdr, sizeof(hdr));
    if (r.error) goto fail;

    if (hdr.magic[0] != GDC_MAGIC_0 || hdr.magic[1] != GDC_MAGIC_1 ||
        hdr.magic[2] != GDC_MAGIC_2 || hdr.magic[3] != GDC_MAGIC_3) {
        r.error = -2;
        goto fail;
    }
    if (hdr.format_version > GDC_FORMAT_VERSION) {
        r.error = -2;
        goto fail;
    }

    /* Read section table */
    int n_sections = (int)hdr.n_sections;
    GdcSectionEntry *entries = calloc((size_t)n_sections, sizeof(GdcSectionEntry));
    if (!entries) { r.error = -1; goto fail; }
    r_bytes(&r, entries, (size_t)n_sections * sizeof(GdcSectionEntry));
    if (r.error) { free(entries); goto fail; }

    /* Process each section by seeking to its offset */
    typedef void (*SectionReader)(Reader *, GdcData *);
    static const SectionReader readers[] = {
        [GDC_SECTION_METADATA]       = read_metadata,
        [GDC_SECTION_ALLELES]        = read_alleles,
        [GDC_SECTION_GENE_USE]       = read_gene_use,
        [GDC_SECTION_TRIM_DISTS]     = read_trim_dists,
        [GDC_SECTION_NP_PARAMS]      = read_np_params,
        [GDC_SECTION_P_NUC_PROBS]    = read_p_nuc_probs,
        [GDC_SECTION_MUTATION_MODEL] = read_mutation_model,
    };

    for (int i = 0; i < n_sections && !r.error; i++) {
        uint32_t id = entries[i].section_id;
        if (id >= GDC_SECTION_COUNT) continue;  /* skip unknown sections */

        fseek(fp, (long)entries[i].offset, SEEK_SET);
        readers[id](&r, data);
    }

    free(entries);

fail:
    fclose(fp);
    if (r.error) {
        gdc_data_destroy(data);
        return r.error;
    }
    return 0;
}

/* ── Populate SimConfig from GdcData ──────────────────────────── */

void gdc_populate_sim_config(SimConfig *cfg, const GdcData *data) {
    cfg->chain_type = (ChainType)data->chain_type;

    /* Copy allele pools (deep copy: allele_pool_add copies each allele) */
    allele_pool_destroy(&cfg->v_alleles);
    allele_pool_destroy(&cfg->d_alleles);
    allele_pool_destroy(&cfg->j_alleles);
    allele_pool_destroy(&cfg->c_alleles);

    cfg->v_alleles = allele_pool_create(data->v_alleles.count);
    cfg->d_alleles = allele_pool_create(data->d_alleles.count);
    cfg->j_alleles = allele_pool_create(data->j_alleles.count);
    cfg->c_alleles = allele_pool_create(data->c_alleles.count);

    for (int i = 0; i < data->v_alleles.count; i++)
        allele_pool_add(&cfg->v_alleles, &data->v_alleles.alleles[i]);
    for (int i = 0; i < data->d_alleles.count; i++)
        allele_pool_add(&cfg->d_alleles, &data->d_alleles.alleles[i]);
    for (int i = 0; i < data->j_alleles.count; i++)
        allele_pool_add(&cfg->j_alleles, &data->j_alleles.alleles[i]);
    for (int i = 0; i < data->c_alleles.count; i++)
        allele_pool_add(&cfg->c_alleles, &data->c_alleles.alleles[i]);

    /* Build per-(family, gene) trim tables AND legacy single-global
     * fallback. The legacy field stays populated from the first GDC
     * entry so any code path that reads cfg->v_trim_3 directly still
     * sees a valid distribution. The per-allele pointers (set below)
     * override the legacy field for production simulations. */
    NamedTrimDistTable *cfg_tables[4] = {
        &cfg->v_trim_3_table, &cfg->d_trim_5_table,
        &cfg->d_trim_3_table, &cfg->j_trim_5_table,
    };
    TrimDist *cfg_legacy[4] = {
        &cfg->v_trim_3, &cfg->d_trim_5,
        &cfg->d_trim_3, &cfg->j_trim_5,
    };
    static const int tt_idx[4] = {
        GDC_TRIM_V3, GDC_TRIM_D5, GDC_TRIM_D3, GDC_TRIM_J5,
    };

    for (int slot = 0; slot < 4; slot++) {
        int tt = tt_idx[slot];
        const TrimDistTable *src_table = &data->trim_dists[tt];
        NamedTrimDistTable *dst_table = cfg_tables[slot];
        TrimDist *legacy = cfg_legacy[slot];

        /* Free any previously built table (idempotent populate). */
        for (int i = 0; i < dst_table->count; i++) {
            free(dst_table->entries[i].dist.probs);
        }
        free(dst_table->entries);
        dst_table->entries = NULL;
        dst_table->count = 0;

        if (src_table->count == 0) continue;

        /* Deep-copy each (family, gene, probs) entry into the table. */
        dst_table->entries = calloc((size_t)src_table->count,
                                    sizeof(NamedTrimDist));
        if (!dst_table->entries) continue;
        dst_table->count = src_table->count;

        for (int i = 0; i < src_table->count; i++) {
            const GeneTrimDist *src = &src_table->dists[i];
            NamedTrimDist *dst = &dst_table->entries[i];
            strncpy(dst->family, src->family_name,
                    GENAIRR_MAX_ALLELE_NAME - 1);
            dst->family[GENAIRR_MAX_ALLELE_NAME - 1] = '\0';
            strncpy(dst->gene, src->gene_name,
                    GENAIRR_MAX_ALLELE_NAME - 1);
            dst->gene[GENAIRR_MAX_ALLELE_NAME - 1] = '\0';

            int n_probs = src->max_trim + 1;
            dst->dist.max_trim = src->max_trim;
            dst->dist.probs = calloc((size_t)n_probs, sizeof(double));
            if (dst->dist.probs) {
                memcpy(dst->dist.probs, src->probs,
                       (size_t)n_probs * sizeof(double));
            }
        }

        /* Legacy fallback: keep populated from the first GDC entry, so
         * code that reads cfg->v_trim_3 (etc.) directly still works. */
        const GeneTrimDist *first = &src_table->dists[0];
        free(legacy->probs);
        legacy->max_trim = first->max_trim;
        int n_legacy = first->max_trim + 1;
        legacy->probs = calloc((size_t)n_legacy, sizeof(double));
        if (legacy->probs) {
            memcpy(legacy->probs, first->probs,
                   (size_t)n_legacy * sizeof(double));
        }
    }

    /* Wire each allele's trim_dist pointers to the matching (family,
     * gene) entry. On miss, fall back to the first table entry — never
     * NULL when the table is non-empty. This guarantees per-segment
     * sampling consistency even for alleles unrepresented in training
     * data. Empty tables leave pointers NULL → trim.c falls back to
     * cfg->v_trim_3 etc. */
    AllelePool *pools[3] = {
        &cfg->v_alleles, &cfg->d_alleles, &cfg->j_alleles,
    };
    for (int p = 0; p < 3; p++) {
        AllelePool *pool = pools[p];
        for (int a = 0; a < pool->count; a++) {
            Allele *al = &pool->alleles[a];
            al->trim_dist_5 = NULL;
            al->trim_dist_3 = NULL;

            const NamedTrimDistTable *t5 = NULL, *t3 = NULL;
            switch (al->segment_type) {
                case SEG_V: t3 = &cfg->v_trim_3_table; break;
                case SEG_D: t5 = &cfg->d_trim_5_table;
                            t3 = &cfg->d_trim_3_table; break;
                case SEG_J: t5 = &cfg->j_trim_5_table; break;
                default:    break;
            }

            const NamedTrimDistTable *needs[2]   = { t5, t3 };
            const TrimDist        **slot[2] = {
                &al->trim_dist_5, &al->trim_dist_3,
            };
            for (int s = 0; s < 2; s++) {
                if (!needs[s] || needs[s]->count == 0) continue;
                const NamedTrimDist *match = NULL;
                for (int i = 0; i < needs[s]->count; i++) {
                    const NamedTrimDist *e = &needs[s]->entries[i];
                    if (strcmp(e->family, al->family) == 0 &&
                        strcmp(e->gene,   al->gene)   == 0) {
                        match = e;
                        break;
                    }
                }
                if (!match) match = &needs[s]->entries[0];
                *slot[s] = &match->dist;
            }
        }
    }

    /* NP region distributions (TdT Markov model). Copy from GdcData
     * into SimConfig so the SimConfig owns its own buffers and
     * outlives the loader's GdcData if needed. Normalize first_base
     * once at load time so the assemble sampler can avoid drift. */
    cfg->n_np_regions = data->n_np_regions;
    for (int i = 0; i < data->n_np_regions && i < 2; i++) {
        const NpParams *src = &data->np[i];
        NpDist *dst = &cfg->np[i];

        /* Length distribution */
        dst->max_length = src->max_length;
        int n_len = src->max_length + 1;
        free(dst->length_probs);
        dst->length_probs = calloc((size_t)n_len, sizeof(double));
        if (!dst->length_probs) continue;
        memcpy(dst->length_probs, src->length_probs,
               (size_t)n_len * sizeof(double));

        /* First-base distribution: copy and normalize */
        double fb_sum = 0.0;
        for (int b = 0; b < 4; b++) fb_sum += src->first_base[b];
        if (fb_sum > 1e-12) {
            for (int b = 0; b < 4; b++)
                dst->first_base[b] = src->first_base[b] / fb_sum;
        } else {
            for (int b = 0; b < 4; b++) dst->first_base[b] = 0.25;
        }

        /* Markov transitions */
        dst->n_positions = src->n_positions;
        int n_trans = src->n_positions * 16;
        free(dst->transitions);
        dst->transitions = NULL;
        if (n_trans > 0) {
            dst->transitions = calloc((size_t)n_trans, sizeof(double));
            if (dst->transitions) {
                memcpy(dst->transitions, src->transitions,
                       (size_t)n_trans * sizeof(double));
            }
        }
    }

    /* P-nucleotide length distribution. Same shared distribution
     * is sampled independently at each eligible segment end. When
     * the GDC carries no p_nuc data, length_probs stays NULL and
     * the assemble step emits no P-nucs. */
    free(cfg->p_nuc_dist.length_probs);
    cfg->p_nuc_dist.length_probs = NULL;
    cfg->p_nuc_dist.max_length = 0;
    if (data->p_nuc_probs && data->p_nuc_max >= 0) {
        int n_p = data->p_nuc_max + 1;
        cfg->p_nuc_dist.length_probs = calloc((size_t)n_p, sizeof(double));
        if (cfg->p_nuc_dist.length_probs) {
            memcpy(cfg->p_nuc_dist.length_probs, data->p_nuc_probs,
                   (size_t)n_p * sizeof(double));
            cfg->p_nuc_dist.max_length = data->p_nuc_max;
        }
    }
}

/* ── Populate S5FModel from GdcData ───────────────────────────── */

void gdc_populate_s5f_model(S5FModel *model, const GdcData *data) {
    if (data->mutation_model_type != GDC_MUTMODEL_S5F) return;

    for (int k = 0; k < GDC_S5F_KMER_SPACE; k++) {
        model->mutability[k] = data->mutability[k];

        int cnt = data->sub_counts[k];
        if (cnt > 0) {
            s5f_set_substitution(model, k,
                                 &data->sub_bases[k * 4],
                                 &data->substitution[k * 4],
                                 cnt);
        }
    }
}
