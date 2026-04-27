/**
 * test_gdc.c — Tests for GDC binary format read/write round-trip.
 */

#include "genairr/gdc_io.h"
#include "genairr/gdc_format.h"
#include "genairr/genairr_api.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PASS(name) printf("  PASS  %s\n", name)

static const char *TEST_FILE = "/tmp/test_genairr.gdc";

/* ── Helper: build a minimal GdcData for testing ──────────────── */

static void add_gene_use(GeneUseTable *t, const char *name, double prob) {
    if (t->count >= t->capacity) {
        int new_cap = t->capacity ? t->capacity * 2 : 8;
        t->entries = realloc(t->entries, (size_t)new_cap * sizeof(GeneUseEntry));
        t->capacity = new_cap;
    }
    strncpy(t->entries[t->count].gene_name, name, GENAIRR_MAX_ALLELE_NAME - 1);
    t->entries[t->count].gene_name[GENAIRR_MAX_ALLELE_NAME - 1] = '\0';
    t->entries[t->count].probability = prob;
    t->count++;
}

static void add_trim_dist(TrimDistTable *t, const char *family,
                           const char *gene, int max_trim) {
    if (t->count >= t->capacity) {
        int new_cap = t->capacity ? t->capacity * 2 : 8;
        t->dists = realloc(t->dists, (size_t)new_cap * sizeof(GeneTrimDist));
        t->capacity = new_cap;
    }
    GeneTrimDist *d = &t->dists[t->count];
    strncpy(d->family_name, family, GENAIRR_MAX_ALLELE_NAME - 1);
    d->family_name[GENAIRR_MAX_ALLELE_NAME - 1] = '\0';
    strncpy(d->gene_name, gene, GENAIRR_MAX_ALLELE_NAME - 1);
    d->gene_name[GENAIRR_MAX_ALLELE_NAME - 1] = '\0';
    d->max_trim = max_trim;
    d->probs = calloc((size_t)(max_trim + 1), sizeof(double));
    double total = 0.0;
    for (int i = 0; i <= max_trim; i++) {
        d->probs[i] = 1.0 / (1.0 + i);
        total += d->probs[i];
    }
    for (int i = 0; i <= max_trim; i++) d->probs[i] /= total;
    t->count++;
}

static GdcData build_test_data(void) {
    GdcData data;
    gdc_data_init(&data);

    /* Metadata */
    strcpy(data.name, "TEST_IGH");
    data.chain_type = GDC_CHAIN_IGH;
    data.has_d = true;
    strcpy(data.species, "Human");
    strcpy(data.reference_set, "IMGT");
    data.build_date = 9557; /* ~2026-03-08 */

    /* Alleles: 2 V, 1 D, 1 J */
    data.v_alleles = allele_pool_create(2);
    Allele v1 = { "IGHV1-2*01", "IGHV1-2", "IGHV1",
                  "ATGATGATGATGATGATG", 18, 15, SEG_V };
    Allele v2 = { "IGHV1-2*02", "IGHV1-2", "IGHV1",
                  "ATGATGATGATGATGATC", 18, 15, SEG_V };
    allele_pool_add(&data.v_alleles, &v1);
    allele_pool_add(&data.v_alleles, &v2);

    data.d_alleles = allele_pool_create(1);
    Allele d1 = { "IGHD1-1*01", "IGHD1-1", "IGHD1",
                  "GATTACA", 7, 0, SEG_D };
    allele_pool_add(&data.d_alleles, &d1);

    data.j_alleles = allele_pool_create(1);
    Allele j1 = { "IGHJ4*02", "IGHJ4", "IGHJ4",
                  "ACTACTTTGACTACTGG", 17, 14, SEG_J };
    allele_pool_add(&data.j_alleles, &j1);

    data.c_alleles = allele_pool_create(0);

    /* Gene use */
    add_gene_use(&data.gene_use[GDC_SEG_V], "IGHV1-2", 0.7);
    add_gene_use(&data.gene_use[GDC_SEG_V], "IGHV3-11", 0.3);
    add_gene_use(&data.gene_use[GDC_SEG_J], "IGHJ4", 1.0);
    add_gene_use(&data.gene_use[GDC_SEG_D], "IGHD1-1", 1.0);

    /* Trim dists */
    add_trim_dist(&data.trim_dists[GDC_TRIM_V3], "IGHV1", "IGHV1-2", 10);
    add_trim_dist(&data.trim_dists[GDC_TRIM_D5], "IGHD1", "IGHD1-1", 5);
    add_trim_dist(&data.trim_dists[GDC_TRIM_D3], "IGHD1", "IGHD1-1", 5);
    add_trim_dist(&data.trim_dists[GDC_TRIM_J5], "IGHJ4", "IGHJ4", 8);

    /* NP params */
    data.n_np_regions = 2;
    for (int r = 0; r < 2; r++) {
        NpParams *np = &data.np[r];
        np->max_length = 20;
        np->length_probs = calloc(21, sizeof(double));
        double total = 0.0;
        for (int i = 0; i <= 20; i++) {
            np->length_probs[i] = (i < 15) ? 1.0 / (1.0 + i) : 0.0;
            total += np->length_probs[i];
        }
        for (int i = 0; i <= 20; i++) np->length_probs[i] /= total;

        np->first_base[0] = 0.25; /* A */
        np->first_base[1] = 0.25; /* C */
        np->first_base[2] = 0.25; /* G */
        np->first_base[3] = 0.25; /* T */

        np->n_positions = 3; /* small for test */
        np->transitions = calloc((size_t)np->n_positions * 16, sizeof(double));
        for (int p = 0; p < np->n_positions; p++) {
            for (int f = 0; f < 4; f++) {
                for (int t = 0; t < 4; t++) {
                    np->transitions[p * 16 + f * 4 + t] = 0.25;
                }
            }
        }
    }

    /* P-nucleotide probs */
    data.p_nuc_max = 4;
    data.p_nuc_probs = calloc(5, sizeof(double));
    data.p_nuc_probs[0] = 0.50;
    data.p_nuc_probs[1] = 0.25;
    data.p_nuc_probs[2] = 0.15;
    data.p_nuc_probs[3] = 0.07;
    data.p_nuc_probs[4] = 0.03;

    /* Mutation model: minimal S5F with a few entries */
    data.mutation_model_type = GDC_MUTMODEL_S5F;
    memset(data.mutability, 0, sizeof(data.mutability));
    memset(data.substitution, 0, sizeof(data.substitution));
    memset(data.sub_counts, 0, sizeof(data.sub_counts));
    memset(data.sub_bases, 0, sizeof(data.sub_bases));

    /* Set a few test entries: AAAAA (key=0), ACGTA (key=140) */
    data.mutability[0] = 0.005;
    data.sub_counts[0] = 3;
    data.sub_bases[0] = 'C'; data.sub_bases[1] = 'G'; data.sub_bases[2] = 'T';
    data.substitution[0] = 0.354; data.substitution[1] = 0.415; data.substitution[2] = 0.231;

    data.mutability[140] = 0.0023;
    data.sub_counts[140] = 2;
    data.sub_bases[140 * 4] = 'A'; data.sub_bases[140 * 4 + 1] = 'T';
    data.substitution[140 * 4] = 0.6; data.substitution[140 * 4 + 1] = 0.4;

    return data;
}

/* ── Tests ─────────────────────────────────────────────────────── */

static void test_write_and_read_roundtrip(void) {
    GdcData orig = build_test_data();

    /* Write */
    int rc = gdc_save(TEST_FILE, &orig);
    assert(rc == 0);

    /* Read back */
    GdcData loaded;
    rc = gdc_load(TEST_FILE, &loaded);
    assert(rc == 0);

    /* Verify metadata */
    assert(strcmp(loaded.name, "TEST_IGH") == 0);
    assert(loaded.chain_type == GDC_CHAIN_IGH);
    assert(loaded.has_d == true);
    assert(strcmp(loaded.species, "Human") == 0);
    assert(strcmp(loaded.reference_set, "IMGT") == 0);
    assert(loaded.build_date == 9557);

    PASS("metadata round-trip");

    /* Verify alleles */
    assert(loaded.v_alleles.count == 2);
    assert(loaded.d_alleles.count == 1);
    assert(loaded.j_alleles.count == 1);
    assert(loaded.c_alleles.count == 0);

    assert(strcmp(loaded.v_alleles.alleles[0].name, "IGHV1-2*01") == 0);
    assert(loaded.v_alleles.alleles[0].length == 18);
    assert(loaded.v_alleles.alleles[0].anchor == 15);
    assert(strcmp(loaded.v_alleles.alleles[0].seq, "ATGATGATGATGATGATG") == 0);
    assert(strcmp(loaded.v_alleles.alleles[0].gene, "IGHV1-2") == 0);
    assert(strcmp(loaded.v_alleles.alleles[0].family, "IGHV1") == 0);
    assert(loaded.v_alleles.alleles[0].segment_type == SEG_V);

    assert(strcmp(loaded.j_alleles.alleles[0].name, "IGHJ4*02") == 0);
    assert(loaded.j_alleles.alleles[0].anchor == 14);
    assert(loaded.j_alleles.alleles[0].segment_type == SEG_J);

    PASS("alleles round-trip");

    /* Verify gene use */
    assert(loaded.gene_use[GDC_SEG_V].count == 2);
    assert(strcmp(loaded.gene_use[GDC_SEG_V].entries[0].gene_name, "IGHV1-2") == 0);
    assert(fabs(loaded.gene_use[GDC_SEG_V].entries[0].probability - 0.7) < 1e-10);
    assert(loaded.gene_use[GDC_SEG_J].count == 1);
    assert(loaded.gene_use[GDC_SEG_D].count == 1);

    PASS("gene_use round-trip");

    /* Verify trim dists */
    assert(loaded.trim_dists[GDC_TRIM_V3].count == 1);
    assert(loaded.trim_dists[GDC_TRIM_V3].dists[0].max_trim == 10);
    assert(strcmp(loaded.trim_dists[GDC_TRIM_V3].dists[0].family_name, "IGHV1") == 0);
    assert(strcmp(loaded.trim_dists[GDC_TRIM_V3].dists[0].gene_name, "IGHV1-2") == 0);
    /* Verify first prob is nonzero and sum ≈ 1.0 */
    double sum = 0.0;
    for (int i = 0; i <= 10; i++) {
        sum += loaded.trim_dists[GDC_TRIM_V3].dists[0].probs[i];
    }
    assert(fabs(sum - 1.0) < 1e-10);

    assert(loaded.trim_dists[GDC_TRIM_D5].count == 1);
    assert(loaded.trim_dists[GDC_TRIM_D3].count == 1);
    assert(loaded.trim_dists[GDC_TRIM_J5].count == 1);

    PASS("trim_dists round-trip");

    /* Verify NP params */
    assert(loaded.n_np_regions == 2);
    assert(loaded.np[0].max_length == 20);
    assert(loaded.np[0].n_positions == 3);
    double np_sum = 0.0;
    for (int i = 0; i <= 20; i++) np_sum += loaded.np[0].length_probs[i];
    assert(fabs(np_sum - 1.0) < 1e-10);
    assert(fabs(loaded.np[0].first_base[0] - 0.25) < 1e-10);
    /* Check uniform transitions */
    assert(fabs(loaded.np[0].transitions[0] - 0.25) < 1e-10);

    PASS("np_params round-trip");

    /* Verify P-nucleotide probs */
    assert(loaded.p_nuc_max == 4);
    assert(fabs(loaded.p_nuc_probs[0] - 0.50) < 1e-10);
    assert(fabs(loaded.p_nuc_probs[4] - 0.03) < 1e-10);

    PASS("p_nuc_probs round-trip");

    /* Verify mutation model */
    assert(loaded.mutation_model_type == GDC_MUTMODEL_S5F);
    assert(fabs(loaded.mutability[0] - 0.005) < 1e-10);
    assert(loaded.sub_counts[0] == 3);
    assert(loaded.sub_bases[0] == 'C');
    assert(loaded.sub_bases[1] == 'G');
    assert(loaded.sub_bases[2] == 'T');
    assert(fabs(loaded.substitution[0] - 0.354) < 1e-10);

    assert(fabs(loaded.mutability[140] - 0.0023) < 1e-10);
    assert(loaded.sub_counts[140] == 2);
    assert(loaded.sub_bases[140 * 4] == 'A');
    assert(fabs(loaded.substitution[140 * 4] - 0.6) < 1e-10);

    /* Zero entries should remain zero */
    assert(loaded.mutability[1] == 0.0);
    assert(loaded.sub_counts[1] == 0);

    PASS("mutation_model round-trip");

    gdc_data_destroy(&loaded);
    gdc_data_destroy(&orig);
    remove(TEST_FILE);
}

static void test_invalid_magic(void) {
    /* Write a full 32-byte header with wrong magic */
    FILE *fp = fopen(TEST_FILE, "wb");
    assert(fp);
    char buf[32];
    memset(buf, 0, sizeof(buf));
    buf[0] = 'N'; buf[1] = 'O'; buf[2] = 'P'; buf[3] = 'E';
    fwrite(buf, 1, 32, fp);
    fclose(fp);

    GdcData data;
    int rc = gdc_load(TEST_FILE, &data);
    assert(rc == -2);  /* format error */

    remove(TEST_FILE);
    PASS("invalid magic rejected");
}

static void test_populate_sim_config(void) {
    GdcData orig = build_test_data();
    int rc = gdc_save(TEST_FILE, &orig);
    assert(rc == 0);

    GdcData loaded;
    rc = gdc_load(TEST_FILE, &loaded);
    assert(rc == 0);

    SimConfig cfg;
    sim_config_init(&cfg, CHAIN_IGH);
    gdc_populate_sim_config(&cfg, &loaded);

    assert(cfg.chain_type == CHAIN_IGH);
    assert(cfg.v_alleles.count == 2);
    assert(cfg.d_alleles.count == 1);
    assert(cfg.j_alleles.count == 1);
    assert(strcmp(cfg.v_alleles.alleles[0].name, "IGHV1-2*01") == 0);

    /* Trim dist should be populated */
    assert(cfg.v_trim_3.probs != NULL);
    assert(cfg.v_trim_3.max_trim == 10);

    /* NP region distributions should be copied from GdcData. The
     * synthetic build_test_data() sets max_length=20 and uniform
     * first_base = 0.25 each. */
    assert(cfg.n_np_regions == 2);
    assert(cfg.np[0].max_length == 20);
    assert(cfg.np[0].length_probs != NULL);
    assert(fabs(cfg.np[0].first_base[0] - 0.25) < 1e-9);
    assert(fabs(cfg.np[0].first_base[3] - 0.25) < 1e-9);

    PASS("populate_sim_config");

    sim_config_destroy(&cfg);
    gdc_data_destroy(&loaded);
    gdc_data_destroy(&orig);
    remove(TEST_FILE);
}

static void test_populate_s5f_model(void) {
    GdcData orig = build_test_data();
    int rc = gdc_save(TEST_FILE, &orig);
    assert(rc == 0);

    GdcData loaded;
    rc = gdc_load(TEST_FILE, &loaded);
    assert(rc == 0);

    S5FModel model;
    s5f_model_init(&model, 0.01, 0.15, true);
    gdc_populate_s5f_model(&model, &loaded);

    assert(fabs(model.mutability[0] - 0.005) < 1e-10);
    assert(model.substitution[0].count == 3);
    assert(model.substitution[0].bases[0] == 'C');

    /* Key 140 */
    assert(fabs(model.mutability[140] - 0.0023) < 1e-10);
    assert(model.substitution[140].count == 2);

    PASS("populate_s5f_model");

    gdc_data_destroy(&loaded);
    gdc_data_destroy(&orig);
    remove(TEST_FILE);
}

static void test_empty_sections(void) {
    GdcData data;
    gdc_data_init(&data);

    strcpy(data.name, "EMPTY");
    data.chain_type = GDC_CHAIN_IGK;
    data.has_d = false;
    strcpy(data.species, "Mouse");
    strcpy(data.reference_set, "Test");

    /* Create empty allele pools */
    data.v_alleles = allele_pool_create(0);
    data.d_alleles = allele_pool_create(0);
    data.j_alleles = allele_pool_create(0);
    data.c_alleles = allele_pool_create(0);

    /* Empty NP */
    data.n_np_regions = 0;

    /* Empty P-nuc */
    data.p_nuc_max = 0;
    data.p_nuc_probs = calloc(1, sizeof(double));
    data.p_nuc_probs[0] = 1.0;

    /* No mutation model */
    data.mutation_model_type = GDC_MUTMODEL_NONE;

    int rc = gdc_save(TEST_FILE, &data);
    assert(rc == 0);

    GdcData loaded;
    rc = gdc_load(TEST_FILE, &loaded);
    assert(rc == 0);

    assert(strcmp(loaded.name, "EMPTY") == 0);
    assert(loaded.chain_type == GDC_CHAIN_IGK);
    assert(loaded.has_d == false);
    assert(loaded.v_alleles.count == 0);
    assert(loaded.n_np_regions == 0);
    assert(loaded.mutation_model_type == GDC_MUTMODEL_NONE);

    PASS("empty sections round-trip");

    gdc_data_destroy(&loaded);
    gdc_data_destroy(&data);
    remove(TEST_FILE);
}

static void test_file_size_reasonable(void) {
    GdcData data = build_test_data();
    int rc = gdc_save(TEST_FILE, &data);
    assert(rc == 0);

    FILE *fp = fopen(TEST_FILE, "rb");
    assert(fp);
    fseek(fp, 0, SEEK_END);
    long size = ftell(fp);
    fclose(fp);

    /* Header (32) + TOC (7×16=112) + sections.
     * Mutation model alone: 3125×8 (mutability) + substitution = ~25-30 KB.
     * Total should be reasonable, under 100 KB for this test data. */
    assert(size > 100);      /* not trivially small */
    assert(size < 200000);   /* not bloated */

    printf("  INFO  GDC file size: %ld bytes\n", size);
    PASS("file size reasonable");

    gdc_data_destroy(&data);
    remove(TEST_FILE);
}

/* ═══════════════════════════════════════════════════════════════
 * T1-5: genairr_create_from_memory must reject bad inputs and
 * succeed on valid GDC bytes. Pre-fix it ignored fwrite/fclose
 * returns; on disk-full the temp file would be silently truncated
 * and genairr_create would read corrupted data.
 *
 * We can't portably inject an fwrite failure from a unit test, so
 * these tests cover (a) the happy path (regression: the fix didn't
 * break valid input) and (b) explicit invalid-input rejection paths.
 * ═══════════════════════════════════════════════════════════════ */

static void test_create_from_memory_happy_path(void) {
    /* Build a minimal valid GDC, save it, slurp the bytes back, then
     * call genairr_create_from_memory. */
    GdcData orig = build_test_data();
    int rc = gdc_save(TEST_FILE, &orig);
    assert(rc == 0);

    FILE *fp = fopen(TEST_FILE, "rb");
    assert(fp != NULL);
    fseek(fp, 0, SEEK_END);
    long size = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    assert(size > 0);

    void *bytes = malloc((size_t)size);
    assert(bytes);
    size_t got = fread(bytes, 1, (size_t)size, fp);
    assert(got == (size_t)size);
    fclose(fp);

    GenAIRRSimulator *sim = genairr_create_from_memory(bytes, (int)size);
    assert(sim != NULL);
    genairr_destroy(sim);

    free(bytes);
    gdc_data_destroy(&orig);
    remove(TEST_FILE);
    PASS("genairr_create_from_memory: happy path returns non-NULL");
}

static void test_create_from_memory_invalid_inputs(void) {
    /* NULL gdc_bytes → NULL */
    GenAIRRSimulator *sim = genairr_create_from_memory(NULL, 100);
    assert(sim == NULL);

    /* gdc_len <= 0 → NULL */
    char buf[16] = {0};
    sim = genairr_create_from_memory(buf, 0);
    assert(sim == NULL);
    sim = genairr_create_from_memory(buf, -1);
    assert(sim == NULL);

    /* Garbage bytes (valid pointer, non-zero length, but not GDC) → NULL.
     * gdc_load's magic check rejects this; the fix is not what makes
     * this case work, but we keep it as a regression. */
    char garbage[64];
    memset(garbage, 'x', sizeof(garbage));
    sim = genairr_create_from_memory(garbage, (int)sizeof(garbage));
    assert(sim == NULL);

    PASS("genairr_create_from_memory: invalid inputs return NULL");
}

/* ── Main ──────────────────────────────────────────────────────── */

int main(void) {
    printf("test_gdc: GDC binary format\n");

    test_write_and_read_roundtrip();
    test_invalid_magic();
    test_populate_sim_config();
    test_populate_s5f_model();
    test_empty_sections();
    test_file_size_reasonable();
    test_create_from_memory_happy_path();
    test_create_from_memory_invalid_inputs();

    printf("\nAll GDC tests passed.\n");
    return 0;
}
