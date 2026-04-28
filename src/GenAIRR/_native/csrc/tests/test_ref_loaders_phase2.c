/**
 * test_ref_loaders_phase2.c — tests for the Phase 2 + Phase 3a loaders.
 *
 * Covers:
 *   - plain FASTA loader (no JSON, no metadata)
 *   - AIRR-C GermlineSet JSON loader (modern source-of-truth format)
 *   - OGRDB combined FASTA + JSON sidecar
 *   - IgBLAST bundle (FASTA + .aux for J / .ndm.imgt for V)
 *   - reference_loader_open_auto dispatcher
 */

#include "genairr/ref_loader.h"
#include "genairr/anchor.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define TEST(name) \
    do { \
        printf("  %-55s", #name); \
        int _r = name(); \
        printf("%s\n", _r == 0 ? "PASS" : "FAIL"); \
        if (_r) failures++; \
        total++; \
    } while (0)

/* ── Test helpers ──────────────────────────────────────────────── */

static char *write_temp_file(const char *body, const char *suffix) {
    char templ[512];
    snprintf(templ, sizeof(templ),
             "/tmp/genairr_loader_test_XXXXXX%s", suffix ? suffix : "");
    /* mkstemps takes the suffix length to preserve. */
    int suffix_len = suffix ? (int)strlen(suffix) : 0;
    int fd = mkstemps(templ, suffix_len);
    if (fd < 0) return NULL;
    FILE *fp = fdopen(fd, "w");
    if (!fp) { close(fd); unlink(templ); return NULL; }
    fputs(body, fp);
    fclose(fp);
    char *path = (char *)malloc(strlen(templ) + 1);
    if (!path) { unlink(templ); return NULL; }
    memcpy(path, templ, strlen(templ) + 1);
    return path;
}

/* ── Plain FASTA loader ─────────────────────────────────────────── */

static int test_plain_fasta_basic(void) {
    const char *body =
        ">IGHV1-2*01\n"
        "acgtacgtacgtacgt\n"
        ">IGHJ4*02\n"
        "tggggcaaaggc\n";
    char *path = write_temp_file(body, ".fasta");
    if (!path) return 1;

    const char *err = NULL;
    ReferenceLoader *L = plain_fasta_loader_open(path, LOCUS_IGH,
                                                  SEG_UNKNOWN, &err);
    int rc = (L != NULL) ? 0 : 1;
    if (rc) { unlink(path); free(path); return 1; }

    LoadedAlleleRecord rec;
    if (reference_loader_next(L, &rec, &err) != 0) goto fail;
    if (strcmp(rec.name, "IGHV1-2*01") != 0) goto fail;
    if (rec.segment != SEG_V) goto fail;
    if (rec.locus != LOCUS_IGH) goto fail;
    if (rec.gap_convention_imgt) goto fail;     /* plain: never IMGT */
    if (rec.gapped_sequence != NULL) goto fail;  /* plain: no gapped */
    if (rec.sequence_length != 16) goto fail;

    if (reference_loader_next(L, &rec, &err) != 0) goto fail;
    if (strcmp(rec.name, "IGHJ4*02") != 0) goto fail;
    if (rec.segment != SEG_J) goto fail;

    if (reference_loader_next(L, &rec, &err) != 1) goto fail;   /* EOF */

    reference_loader_close(L);
    unlink(path); free(path);
    return 0;
fail:
    reference_loader_close(L);
    unlink(path); free(path);
    return 1;
}

static int test_plain_fasta_locus_hint_used_when_name_unknown(void) {
    /* Allele name with no recognizable locus → hint should be used. */
    const char *body =
        ">novel_allele_X\n"
        "acgtacgtacgtacgt\n";
    char *path = write_temp_file(body, ".fasta");
    if (!path) return 1;

    const char *err = NULL;
    ReferenceLoader *L = plain_fasta_loader_open(path, LOCUS_TRB,
                                                  SEG_V, &err);
    if (!L) { unlink(path); free(path); return 1; }

    LoadedAlleleRecord rec;
    int n = reference_loader_next(L, &rec, &err);
    int rc = 0;
    if (n != 0) rc = 1;
    else if (rec.locus != LOCUS_TRB) rc = 1;
    else if (rec.segment != SEG_V) rc = 1;
    reference_loader_close(L);
    unlink(path); free(path);
    return rc;
}

/* ── AIRR-C JSON loader ─────────────────────────────────────────── */

static int test_airrc_basic(void) {
    /* Minimal AIRR-C-shaped JSON with one V allele and one J allele,
     * including explicit anchor coordinates. */
    const char *body =
        "{\n"
        "  \"GermlineSet\": [\n"
        "    {\n"
        "      \"germline_set_name\": \"test\",\n"
        "      \"species\": {\"label\": \"Homo sapiens\"},\n"
        "      \"allele_descriptions\": [\n"
        "        {\n"
        "          \"label\": \"IGHV1-2*02\",\n"
        "          \"locus\": \"IGH\",\n"
        "          \"sequence\": \"acgtacgtacgttgtaaa\",\n"
        "          \"functional\": true,\n"
        "          \"v_gene_delineations\": [\n"
        "            {\"delineation_scheme\": \"IMGT\","
        "             \"cdr3_start\": 13}\n"
        "          ]\n"
        "        },\n"
        "        {\n"
        "          \"label\": \"IGHJ4*01\",\n"
        "          \"locus\": \"IGH\",\n"
        "          \"sequence\": \"tggggcaaaggc\",\n"
        "          \"functional\": true,\n"
        "          \"j_cdr3_end\": 3\n"
        "        }\n"
        "      ]\n"
        "    }\n"
        "  ]\n"
        "}\n";
    char *path = write_temp_file(body, ".json");
    if (!path) return 1;

    const char *err = NULL;
    ReferenceLoader *L = airrc_germline_loader_open(path, &err);
    if (!L) { unlink(path); free(path); return 1; }

    LoadedAlleleRecord rec;
    int rc = 0;
    if (reference_loader_next(L, &rec, &err) != 0) rc = 1;
    else if (strcmp(rec.name, "IGHV1-2*02") != 0) rc = 1;
    else if (rec.locus != LOCUS_IGH) rc = 1;
    else if (rec.segment != SEG_V) rc = 1;
    /* cdr3_start = 13 (1-based) → explicit_anchor = 12 (0-based). */
    else if (rec.explicit_anchor != 12) rc = 1;
    else if (rec.functional_status != FUNC_F) rc = 1;
    else if (strcmp(rec.species, "Homo sapiens") != 0) rc = 1;
    else if (strcmp(rec.source, "airrc-germline-set") != 0) rc = 1;

    if (!rc) {
        if (reference_loader_next(L, &rec, &err) != 0) rc = 1;
        else if (strcmp(rec.name, "IGHJ4*01") != 0) rc = 1;
        else if (rec.segment != SEG_J) rc = 1;
        /* j_cdr3_end = 3 (1-based) → anchor = 0 (0-based codon start). */
        else if (rec.explicit_anchor != 0) rc = 1;
    }

    if (!rc && reference_loader_next(L, &rec, &err) != 1) rc = 1;

    reference_loader_close(L);
    unlink(path); free(path);
    return rc;
}

/* ── OGRDB combined loader ──────────────────────────────────────── */

static int test_ogrdb_with_sidecar(void) {
    const char *fasta_body =
        ">IGHV1-2*02\n"
        "acgtacgtacgttgtaaa\n";
    const char *json_body =
        "{\"GermlineSet\":[{\"allele_descriptions\":["
        "{\"label\":\"IGHV1-2*02\",\"functional\":true,"
        " \"v_gene_delineations\":["
        "  {\"delineation_scheme\":\"IMGT\",\"cdr3_start\":13}]}]}]}";
    char *fasta_path = write_temp_file(fasta_body, ".fasta");
    char *json_path  = write_temp_file(json_body, ".json");
    if (!fasta_path || !json_path) {
        if (fasta_path) { unlink(fasta_path); free(fasta_path); }
        if (json_path)  { unlink(json_path);  free(json_path);  }
        return 1;
    }
    const char *err = NULL;
    ReferenceLoader *L = ogrdb_loader_open(fasta_path, json_path,
                                            LOCUS_IGH, SEG_V, &err);
    int rc = 0;
    if (!L) rc = 1;
    else {
        LoadedAlleleRecord rec;
        if (reference_loader_next(L, &rec, &err) != 0) rc = 1;
        else if (strcmp(rec.name, "IGHV1-2*02") != 0) rc = 1;
        else if (rec.explicit_anchor != 12) rc = 1;       /* from sidecar */
        else if (rec.functional_status != FUNC_F) rc = 1;  /* from sidecar */
        else if (strcmp(rec.source, "ogrdb") != 0) rc = 1;
        if (!rc && reference_loader_next(L, &rec, &err) != 1) rc = 1;
        reference_loader_close(L);
    }
    unlink(fasta_path); free(fasta_path);
    unlink(json_path);  free(json_path);
    return rc;
}

static int test_ogrdb_without_sidecar(void) {
    const char *fasta_body =
        ">IGHV1-2*02\n"
        "acgtacgt\n";
    char *fasta_path = write_temp_file(fasta_body, ".fasta");
    if (!fasta_path) return 1;

    const char *err = NULL;
    ReferenceLoader *L = ogrdb_loader_open(fasta_path, NULL,
                                            LOCUS_IGH, SEG_V, &err);
    int rc = 0;
    if (!L) rc = 1;
    else {
        LoadedAlleleRecord rec;
        if (reference_loader_next(L, &rec, &err) != 0) rc = 1;
        else if (rec.explicit_anchor != -1) rc = 1;       /* no sidecar */
        else if (rec.functional_status != FUNC_UNKNOWN) rc = 1;
        reference_loader_close(L);
    }
    unlink(fasta_path); free(fasta_path);
    return rc;
}

/* ── IgBLAST bundle loader ──────────────────────────────────────── */

static int test_igblast_with_aux_j(void) {
    const char *fasta_body =
        ">IGHJ4*01\n"
        "tggggcaaaggc\n";
    /* .aux columns: name <TAB> coding_frame <TAB> chain <TAB> cdr3_end (0-based)
     * For "tggggcaaaggc", the conserved Trp ends at position 2 (TGG = pos 0-2);
     * cdr3_end = 2 → anchor = cdr3_end - 2 = 0. */
    const char *aux_body =
        "IGHJ4*01\t0\tVH\t2\n";
    char *fasta = write_temp_file(fasta_body, ".fasta");
    char *aux   = write_temp_file(aux_body, ".aux");
    if (!fasta || !aux) {
        if (fasta) { unlink(fasta); free(fasta); }
        if (aux)   { unlink(aux);   free(aux);   }
        return 1;
    }
    const char *err = NULL;
    ReferenceLoader *L = igblast_loader_open(fasta, aux,
                                              LOCUS_IGH, SEG_J, &err);
    int rc = 0;
    if (!L) rc = 1;
    else {
        LoadedAlleleRecord rec;
        if (reference_loader_next(L, &rec, &err) != 0) rc = 1;
        else if (rec.explicit_anchor != 0) rc = 1;
        else if (strcmp(rec.source, "igblast-bundle") != 0) rc = 1;
        reference_loader_close(L);
    }
    unlink(fasta); free(fasta);
    unlink(aux);   free(aux);
    return rc;
}

/* ── Auto-dispatcher ────────────────────────────────────────────── */

static int test_auto_dispatch_json_to_airrc(void) {
    const char *body = "[]";   /* empty array — valid AIRR-C path */
    char *path = write_temp_file(body, ".json");
    if (!path) return 1;
    const char *err = NULL;
    LoaderHints hints = {.locus_hint = LOCUS_IGH, .segment_hint = SEG_V};
    ReferenceLoader *L = reference_loader_open_auto(path, &hints, &err);
    int rc = (L != NULL) ? 0 : 1;
    if (L) reference_loader_close(L);
    unlink(path); free(path);
    return rc;
}

static int test_auto_dispatch_imgt_fasta(void) {
    const char *body =
        ">X1|IGHV1-2*02|Homo sapiens|F|V-REGION|...\n"
        "acgtacgt\n";
    char *path = write_temp_file(body, ".fasta");
    if (!path) return 1;
    const char *err = NULL;
    ReferenceLoader *L = reference_loader_open_auto(path, NULL, &err);
    int rc = 0;
    if (!L) rc = 1;
    else {
        LoadedAlleleRecord rec;
        if (reference_loader_next(L, &rec, &err) != 0) rc = 1;
        else if (strcmp(rec.source, "imgt-vquest") != 0) rc = 1;
        reference_loader_close(L);
    }
    unlink(path); free(path);
    return rc;
}

static int test_auto_dispatch_plain_fasta(void) {
    const char *body =
        ">IGHV1-2*02\n"
        "acgtacgt\n";
    char *path = write_temp_file(body, ".fasta");
    if (!path) return 1;
    const char *err = NULL;
    LoaderHints hints = {.locus_hint = LOCUS_IGH, .segment_hint = SEG_V};
    ReferenceLoader *L = reference_loader_open_auto(path, &hints, &err);
    int rc = 0;
    if (!L) rc = 1;
    else {
        LoadedAlleleRecord rec;
        if (reference_loader_next(L, &rec, &err) != 0) rc = 1;
        else if (strcmp(rec.source, "plain-fasta") != 0) rc = 1;
        reference_loader_close(L);
    }
    unlink(path); free(path);
    return rc;
}

int main(void) {
    int total = 0, failures = 0;
    printf("Phase 2/3 reference loader tests\n");
    printf("─────────────────────────────────\n");

    TEST(test_plain_fasta_basic);
    TEST(test_plain_fasta_locus_hint_used_when_name_unknown);
    TEST(test_airrc_basic);
    TEST(test_ogrdb_with_sidecar);
    TEST(test_ogrdb_without_sidecar);
    TEST(test_igblast_with_aux_j);
    TEST(test_auto_dispatch_json_to_airrc);
    TEST(test_auto_dispatch_imgt_fasta);
    TEST(test_auto_dispatch_plain_fasta);

    printf("\n%d/%d tests passed\n", total - failures, total);
    return failures == 0 ? 0 : 1;
}
