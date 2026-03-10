/**
 * gdc_io.h — Read and write GenAIRR DataConfig binary files (.gdc).
 *
 * gdc_load() reads a .gdc file and populates a GdcData struct with all
 * reference data needed for simulation: alleles, gene usage, trim
 * distributions, NP parameters, and mutation model tables.
 *
 * gdc_save() writes a GdcData struct to a .gdc file.
 *
 * Both functions return 0 on success, negative on error.
 */

#ifndef GENAIRR_GDC_IO_H
#define GENAIRR_GDC_IO_H

#include "gdc_format.h"
#include "allele.h"
#include "sim_config.h"
#include "s5f.h"
#include <stdbool.h>

/* ── Gene usage entry ─────────────────────────────────────────── */

typedef struct {
    char   gene_name[GENAIRR_MAX_ALLELE_NAME];
    double probability;
} GeneUseEntry;

typedef struct {
    GeneUseEntry *entries;
    int           count;
    int           capacity;
} GeneUseTable;

/* ── Per-gene trim distribution ───────────────────────────────── */

typedef struct {
    char    gene_name[GENAIRR_MAX_ALLELE_NAME];
    char    family_name[GENAIRR_MAX_ALLELE_NAME];
    double *probs;      /* dense array [0..max_trim] */
    int     max_trim;
} GeneTrimDist;

typedef struct {
    GeneTrimDist *dists;
    int           count;
    int           capacity;
} TrimDistTable;

/* ── NP region parameters ─────────────────────────────────────── */

typedef struct {
    double *length_probs;       /* dense [0..max_length] */
    int     max_length;
    double  first_base[4];      /* P(A), P(C), P(G), P(T) */
    int     n_positions;        /* number of Markov transition positions */
    double *transitions;        /* flat [n_positions × 16] row-major 4×4 */
} NpParams;

/* ── Loaded config ────────────────────────────────────────────── */

typedef struct {
    /* Metadata */
    char      name[256];
    uint8_t   chain_type;       /* GDC_CHAIN_* */
    bool      has_d;
    char      species[64];
    char      reference_set[64];
    uint32_t  build_date;       /* days since 2000-01-01 */

    /* Allele pools (one per segment type) */
    AllelePool v_alleles;
    AllelePool d_alleles;
    AllelePool j_alleles;
    AllelePool c_alleles;

    /* Gene usage (one table per segment type) */
    GeneUseTable gene_use[4];   /* indexed by GDC_SEG_V/D/J/C */

    /* Trim distributions (one table per trim type) */
    TrimDistTable trim_dists[4]; /* indexed by GDC_TRIM_V3/D5/D3/J5 */

    /* NP parameters */
    NpParams np[2];             /* [0]=NP1, [1]=NP2 */
    int      n_np_regions;      /* 1 for VJ, 2 for VDJ */

    /* P-nucleotide length probs */
    double  *p_nuc_probs;       /* dense [0..p_nuc_max] */
    int      p_nuc_max;

    /* Mutation model */
    uint8_t  mutation_model_type;  /* GDC_MUTMODEL_* */
    double   mutability[GDC_S5F_KMER_SPACE];
    double   substitution[GDC_S5F_KMER_SPACE * 4]; /* flat [kmer*4 + base] */
    uint8_t  sub_counts[GDC_S5F_KMER_SPACE];       /* valid entries per kmer */
    char     sub_bases[GDC_S5F_KMER_SPACE * 4];    /* target bases per kmer */
} GdcData;

/* ── Public API ───────────────────────────────────────────────── */

/** Initialize a GdcData struct to zero/empty state. */
void gdc_data_init(GdcData *data);

/** Free all dynamically allocated memory in a GdcData struct. */
void gdc_data_destroy(GdcData *data);

/**
 * Load a .gdc file into a GdcData struct.
 * @return 0 on success, -1 on I/O error, -2 on format error.
 */
int gdc_load(const char *path, GdcData *data);

/**
 * Write a GdcData struct to a .gdc file.
 * @return 0 on success, -1 on I/O error.
 */
int gdc_save(const char *path, const GdcData *data);

/**
 * Populate a SimConfig's reference data fields from a loaded GdcData.
 * Copies allele pools and trim distributions into the SimConfig.
 * Runtime parameters (feature flags, error rates) are not touched.
 */
void gdc_populate_sim_config(SimConfig *cfg, const GdcData *data);

/**
 * Populate an S5FModel's mutability and substitution tables from GdcData.
 */
void gdc_populate_s5f_model(S5FModel *model, const GdcData *data);

#endif /* GENAIRR_GDC_IO_H */
