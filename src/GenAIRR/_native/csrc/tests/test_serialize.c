/**
 * test_serialize.c — Tests for AIRR serialization (airr_serialize).
 *
 * Validates that airr_serialize() produces a complete AirrRecord
 * from ASeq + SimRecord, including:
 *   1. Sequence and germline extraction
 *   2. Position derivation
 *   3. Junction extraction
 *   4. NP region extraction
 *   5. Mutation collection
 *   6. Allele call construction (with corrections)
 *   7. TSV output
 *   8. AlleleCorrectionSet lifecycle
 */

#include "genairr/genairr.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

#define TEST(name) \
    do { \
        printf("  %-55s", #name); \
        int _r = name(); \
        printf("%s\n", _r == 0 ? "PASS" : "FAIL"); \
        if (_r) failures++; \
        total++; \
    } while (0)

/* ── Helper: build a test allele ──────────────────────────────── */

static Allele make_allele(const char *name, const char *seq,
                          int anchor, Segment seg) {
    Allele a = {0};
    strncpy(a.name, name, sizeof(a.name) - 1);
    strncpy(a.seq, seq, sizeof(a.seq) - 1);
    a.length = strlen(seq);
    a.anchor = anchor;
    a.segment_type = seg;
    return a;
}

/* ═══════════════════════════════════════════════════════════════
 * Test 1: basic serialization without corrections
 * ═══════════════════════════════════════════════════════════════ */

static int test_basic_serialize(void) {
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    /* V(8bp) + NP1(3bp) + D(4bp) + NP2(3bp) + J(6bp) = 24bp */
    Allele v = make_allele("IGHV1-2*01", "AAAAACCCCC", 7, SEG_V);
    Allele d = make_allele("IGHD3-3*01", "GGTTTTGG",   0, SEG_D);
    Allele j = make_allele("IGHJ4*02",   "GGGTTTTTT",  2, SEG_J);

    rec.v_allele = &v;
    rec.d_allele = &d;
    rec.j_allele = &j;
    rec.v_trim_3 = 2;
    rec.d_trim_5 = 2;
    rec.d_trim_3 = 2;
    rec.j_trim_5 = 3;
    rec.productive = true;
    rec.vj_in_frame = true;

    aseq_append_segment(&seq, "AAAAACCC", 8, SEG_V, 0, 7);
    aseq_append_np(&seq, "ATA", 3, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "TTTT", 4, SEG_D, 2, -1);
    aseq_append_np(&seq, "GCG", 3, SEG_NP2, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "TTTTTT", 6, SEG_J, 3, 2);

    AirrRecord airr;
    SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    /* Sequence */
    if (strcmp(airr.sequence, "AAAAACCCATATTTTGCGTTTTTT") != 0) {
        fprintf(stderr, "Expected AAAAACCCATATTTTGCGTTTTTT, got %s\n",
                airr.sequence);
        return 1;
    }
    if (airr.sequence_length != 24) { fprintf(stderr, "FAIL: seq_len=%d\n", airr.sequence_length); return 1; }

    /* Germline alignment: V+D+J use germline, NP uses N (but boundary extension
     * may patch some NP positions with real germline bases) */
    if (strncmp(airr.germline_alignment, "AAAAACCC", 8) != 0) {
        fprintf(stderr, "FAIL: germ V='%.8s'\n", airr.germline_alignment); return 1;
    }
    if (strncmp(airr.germline_alignment + 8, "NNN", 3) != 0) {
        fprintf(stderr, "FAIL: germ NP1='%.3s'\n", airr.germline_alignment + 8); return 1;
    }
    if (strncmp(airr.germline_alignment + 11, "TTTT", 4) != 0) {
        fprintf(stderr, "FAIL: germ D='%.4s'\n", airr.germline_alignment + 11); return 1;
    }
    /* NP2: D 3' extends 1bp (pos 15→G), J 5' extends 1bp (pos 17→G) */
    if (strncmp(airr.germline_alignment + 15, "GNG", 3) != 0) {
        fprintf(stderr, "FAIL: germ NP2='%.3s'\n", airr.germline_alignment + 15); return 1;
    }

    /* Allele calls (no corrections → just true names) */
    if (strcmp(airr.v_call, "IGHV1-2*01") != 0) { fprintf(stderr, "FAIL: v_call=%s\n", airr.v_call); return 1; }
    if (strcmp(airr.d_call, "IGHD3-3*01") != 0) { fprintf(stderr, "FAIL: d_call=%s\n", airr.d_call); return 1; }
    if (strcmp(airr.j_call, "IGHJ4*02") != 0) { fprintf(stderr, "FAIL: j_call=%s\n", airr.j_call); return 1; }

    /* NP regions */
    if (strcmp(airr.np1_region, "ATA") != 0) { fprintf(stderr, "FAIL: np1=%s\n", airr.np1_region); return 1; }
    if (airr.np1_length != 3) { fprintf(stderr, "FAIL: np1_len=%d\n", airr.np1_length); return 1; }
    if (strcmp(airr.np2_region, "GCG") != 0) { fprintf(stderr, "FAIL: np2=%s\n", airr.np2_region); return 1; }
    if (airr.np2_length != 3) { fprintf(stderr, "FAIL: np2_len=%d\n", airr.np2_length); return 1; }

    /* Productivity */
    if (!airr.productive) { fprintf(stderr, "FAIL: productive=%d\n", airr.productive); return 1; }
    if (!airr.vj_in_frame) { fprintf(stderr, "FAIL: vj_in_frame=%d\n", airr.vj_in_frame); return 1; }
    if (airr.stop_codon) { fprintf(stderr, "FAIL: stop_codon=%d\n", airr.stop_codon); return 1; }

    /* Mutations (boundary-extended regions that match don't count as mutations) */
    if (airr.n_mutations != 0) { fprintf(stderr, "FAIL: n_mut=%d rate=%f\n", airr.n_mutations, airr.mutation_rate); return 1; }
    if (airr.mutation_rate != 0.0) { fprintf(stderr, "FAIL: mut_rate=%f\n", airr.mutation_rate); return 1; }

    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * Test 2: mutation annotations collected correctly
 * ═══════════════════════════════════════════════════════════════ */

static int test_mutation_collection(void) {
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    Allele v = make_allele("V1", "AAAAAAAAAA", 7, SEG_V);
    Allele j = make_allele("J1", "TTTTTT", 2, SEG_J);
    rec.v_allele = &v;
    rec.j_allele = &j;

    aseq_append_segment(&seq, "AAAAAAAAAA", 10, SEG_V, 0, 7);
    aseq_append_np(&seq, "CC", 2, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "TTTTTT", 6, SEG_J, 0, 2);

    /* Mutate positions 2 and 5 in V */
    Nuc *n = seq.head->next->next;  /* pos 2 */
    aseq_mutate(&seq, n, 'G', NUC_FLAG_MUTATED);
    n = n->next->next->next;         /* pos 5 */
    aseq_mutate(&seq, n, 'C', NUC_FLAG_MUTATED);

    AirrRecord airr;
    SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    if (airr.n_mutations != 2) {
        fprintf(stderr, "Expected 2 mutations, got %d\n", airr.n_mutations);
        return 1;
    }

    /* Mutation string should be "2:A>G,5:A>C" */
    if (strcmp(airr.mutations, "2:A>G,5:A>C") != 0) {
        fprintf(stderr, "Expected '2:A>G,5:A>C', got '%s'\n", airr.mutations);
        return 1;
    }

    /* Mutation rate: 2 mutated / 16 germline (10 V + 6 J) */
    if (fabs(airr.mutation_rate - 2.0/16.0) > 0.001) {
        fprintf(stderr, "Expected rate ~%.4f, got %.4f\n",
                2.0/16.0, airr.mutation_rate);
        return 1;
    }

    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * Test 3: sequencing and PCR errors
 * ═══════════════════════════════════════════════════════════════ */

static int test_error_collection(void) {
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    Allele v = make_allele("V1", "AAAAAA", 3, SEG_V);
    Allele j = make_allele("J1", "TTTTTT", 2, SEG_J);
    rec.v_allele = &v;
    rec.j_allele = &j;

    aseq_append_segment(&seq, "AAAAAA", 6, SEG_V, 0, 3);
    aseq_append_segment(&seq, "TTTTTT", 6, SEG_J, 0, 2);

    /* Add a sequencing error at pos 1 and PCR error at pos 8 */
    Nuc *n1 = seq.head->next;
    aseq_mutate(&seq, n1, 'G', NUC_FLAG_SEQ_ERROR);
    Nuc *n8 = seq.head;
    for (int i = 0; i < 8; i++) n8 = n8->next;
    aseq_mutate(&seq, n8, 'C', NUC_FLAG_PCR_ERROR);

    AirrRecord airr;
    SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    if (airr.n_sequencing_errors != 1) return 1;
    if (airr.n_pcr_errors != 1) return 1;
    if (strcmp(airr.sequencing_errors, "1:A>G") != 0) {
        fprintf(stderr, "Expected '1:A>G', got '%s'\n", airr.sequencing_errors);
        return 1;
    }
    if (strcmp(airr.pcr_errors, "8:T>C") != 0) {
        fprintf(stderr, "Expected '8:T>C', got '%s'\n", airr.pcr_errors);
        return 1;
    }

    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * Test 4: junction extraction
 * ═══════════════════════════════════════════════════════════════ */

static int test_junction_extraction(void) {
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    /* V anchor at germline pos 4, J anchor at germline pos 2 */
    Allele v = make_allele("V1", "AAAAACCCCC", 4, SEG_V);
    Allele j = make_allele("J1", "GGGTTTTTT",  2, SEG_J);
    rec.v_allele = &v;
    rec.j_allele = &j;

    aseq_append_segment(&seq, "AAAAACCCCC", 10, SEG_V, 0, 4);
    aseq_append_np(&seq, "AT", 2, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "GGGTTTTTT", 9, SEG_J, 0, 2);

    AirrRecord airr;
    SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    /* Junction: V-anchor (pos 4) to J-anchor+3 (pos 14+3=17, AIRR convention) */
    if (airr.junction_length != 13) {
        fprintf(stderr, "Expected junction_length=13, got %d\n",
                airr.junction_length);
        return 1;
    }

    /* Junction: V[4..10) + NP1 "AT" + J[0..5) = "ACCCCCATGGGTT" */
    if (strncmp(airr.junction_nt, "ACCCCCATGGGTT", 13) != 0) {
        fprintf(stderr, "Expected junction 'ACCCCCATGGGTT', got '%.*s'\n",
                airr.junction_length, airr.junction_nt);
        return 1;
    }

    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * Test 5: VJ chain (no D segment)
 * ═══════════════════════════════════════════════════════════════ */

static int test_vj_chain(void) {
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    Allele v = make_allele("IGKV1-5*01", "AAAAACCCCC", 7, SEG_V);
    Allele j = make_allele("IGKJ2*01",   "GGGTTTTTT",  2, SEG_J);
    rec.v_allele = &v;
    rec.j_allele = &j;
    rec.v_trim_3 = 1;
    rec.j_trim_5 = 2;
    /* No D allele */

    aseq_append_segment(&seq, "AAAAACCCC", 9, SEG_V, 0, 7);
    aseq_append_np(&seq, "TAT", 3, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "GTTTTTT", 7, SEG_J, 2, 2);

    AirrRecord airr;
    SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    /* D fields should be zero */
    if (airr.d_sequence_start != 0 || airr.d_sequence_end != 0) return 1;
    if (strlen(airr.d_call) != 0) return 1;

    /* NP2 should be empty */
    if (airr.np2_length != 0) return 1;

    /* V and J should be populated */
    if (strcmp(airr.v_call, "IGKV1-5*01") != 0) return 1;
    if (strcmp(airr.j_call, "IGKJ2*01") != 0) return 1;

    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * Test 6: position derivation with V 3' ambiguity
 * ═══════════════════════════════════════════════════════════════ */

static int test_positions_with_ambiguity(void) {
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    /* V trimmed by 3, but NP1 starts with 2 matching bases */
    Allele v = make_allele("V1", "AAAAAAACCC", 7, SEG_V);
    rec.v_allele = &v;
    rec.v_trim_3 = 3;

    Allele j = make_allele("J1", "TTTTTT", 2, SEG_J);
    rec.j_allele = &j;
    rec.j_trim_5 = 0;

    aseq_append_segment(&seq, "AAAAAAA", 7, SEG_V, 0, 7);
    aseq_append_np(&seq, "CCT", 3, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "TTTTTT", 6, SEG_J, 0, 2);

    AirrRecord airr;
    SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    /* V should extend by 2 due to ambiguity */
    if (airr.v_sequence_end != 9) {
        fprintf(stderr, "Expected v_sequence_end=9, got %d\n",
                airr.v_sequence_end);
        return 1;
    }
    /* Adjusted trim: 3 - 2 = 1 */
    if (airr.v_trim_3 != 1) {
        fprintf(stderr, "Expected v_trim_3=1, got %d\n", airr.v_trim_3);
        return 1;
    }

    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * Test 7: flags propagation
 * ═══════════════════════════════════════════════════════════════ */

static int test_flags_propagation(void) {
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    Allele v = make_allele("V1", "AAAAAA", 3, SEG_V);
    Allele j = make_allele("J1", "TTTTTT", 2, SEG_J);
    rec.v_allele = &v;
    rec.j_allele = &j;
    rec.productive = false;
    rec.stop_codon = true;
    rec.vj_in_frame = false;
    rec.is_reverse_complement = true;
    rec.is_contaminant = false;
    strncpy(rec.note, "Stop codon at pos 42", sizeof(rec.note) - 1);

    aseq_append_segment(&seq, "AAAAAA", 6, SEG_V, 0, 3);
    aseq_append_segment(&seq, "TTTTTT", 6, SEG_J, 0, 2);

    AirrRecord airr;
    SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    if (airr.productive) return 1;
    if (!airr.stop_codon) return 1;
    if (airr.vj_in_frame) return 1;
    if (!airr.is_reverse_complement) return 1;
    if (airr.is_contaminant) return 1;
    if (strcmp(airr.note, "Stop codon at pos 42") != 0) return 1;

    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * Test 8: AlleleCorrectionSet build and destroy
 * ═══════════════════════════════════════════════════════════════ */

static int test_correction_context_lifecycle(void) {
    SimConfig cfg;
    sim_config_init(&cfg, CHAIN_IGH);

    /* Add a few alleles */
    Allele v1 = make_allele("V1", "AAAAACCCCC", 7, SEG_V);
    Allele v2 = make_allele("V2", "AAAAAGGGGG", 7, SEG_V);
    allele_pool_add(&cfg.v_alleles, &v1);
    allele_pool_add(&cfg.v_alleles, &v2);

    Allele d1 = make_allele("D1", "TTTTGG", 0, SEG_D);
    Allele d2 = make_allele("D2", "TTTTCC", 0, SEG_D);
    allele_pool_add(&cfg.d_alleles, &d1);
    allele_pool_add(&cfg.d_alleles, &d2);

    Allele j1 = make_allele("J1", "GGGTTTTTT", 2, SEG_J);
    allele_pool_add(&cfg.j_alleles, &j1);

    AlleleCorrectionSet ctx = allele_correction_set_build(&cfg);

    /* V bitmap should have 2 alleles */
    if (ctx.v_bitmap.n_alleles != 2) {
        fprintf(stderr, "Expected 2 V alleles in bitmap, got %d\n",
                ctx.v_bitmap.n_alleles);
        allele_correction_set_destroy(&ctx);
        sim_config_destroy(&cfg);
        return 1;
    }

    /* D bitmap should have 2 alleles */
    if (ctx.d_bitmap.n_alleles != 2) {
        allele_correction_set_destroy(&ctx);
        sim_config_destroy(&cfg);
        return 1;
    }

    /* J bitmap should have 1 allele */
    if (ctx.j_bitmap.n_alleles != 1) {
        allele_correction_set_destroy(&ctx);
        sim_config_destroy(&cfg);
        return 1;
    }

    /* Short-D threshold should be set */
    if (ctx.d_bitmap.short_d_threshold <= 0) {
        allele_correction_set_destroy(&ctx);
        sim_config_destroy(&cfg);
        return 1;
    }

    allele_correction_set_destroy(&ctx);
    sim_config_destroy(&cfg);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * Test 9: serialize with corrections applied
 * ═══════════════════════════════════════════════════════════════ */

static int test_serialize_with_corrections(void) {
    SimConfig cfg;
    sim_config_init(&cfg, CHAIN_IGH);

    /* Two V alleles that share the same first 5 bases */
    Allele v1 = make_allele("V1", "AAAAACCCCC", 7, SEG_V);
    Allele v2 = make_allele("V2", "AAAAAGGGGG", 7, SEG_V);
    allele_pool_add(&cfg.v_alleles, &v1);
    allele_pool_add(&cfg.v_alleles, &v2);

    Allele j1 = make_allele("J1", "GGGTTTTTT", 2, SEG_J);
    allele_pool_add(&cfg.j_alleles, &j1);

    Allele d1 = make_allele("D1", "TTTTGG", 0, SEG_D);
    allele_pool_add(&cfg.d_alleles, &d1);

    AlleleCorrectionSet ctx = allele_correction_set_build(&cfg);

    /* Build sequence using V1 trimmed by 5 (→ only "AAAAA" remains,
     * which is also a substring of V2's sequence) */
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    rec.v_allele = &v1;
    rec.j_allele = &j1;
    rec.d_allele = &d1;
    rec.v_trim_3 = 5;
    rec.d_trim_5 = 2;
    rec.d_trim_3 = 2;
    rec.j_trim_5 = 0;

    aseq_append_segment(&seq, "AAAAA", 5, SEG_V, 0, 7);
    aseq_append_np(&seq, "AT", 2, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "TT", 2, SEG_D, 2, -1);
    aseq_append_np(&seq, "GA", 2, SEG_NP2, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "GGGTTTTTT", 9, SEG_J, 0, 2);

    AirrRecord airr;
    airr_serialize(&seq, &rec, &cfg, &ctx, &airr);

    /* V call should contain both V1 and V2 (trim-based correction:
     * "AAAAA" is substring of both alleles) */
    if (strstr(airr.v_call, "V1") == NULL) {
        fprintf(stderr, "v_call missing V1: '%s'\n", airr.v_call);
        allele_correction_set_destroy(&ctx);
        sim_config_destroy(&cfg);
        return 1;
    }
    if (strstr(airr.v_call, "V2") == NULL) {
        fprintf(stderr, "v_call missing V2: '%s'\n", airr.v_call);
        allele_correction_set_destroy(&ctx);
        sim_config_destroy(&cfg);
        return 1;
    }

    /* V1 should be first (true allele) */
    if (strncmp(airr.v_call, "V1", 2) != 0) {
        fprintf(stderr, "V1 should be first in v_call: '%s'\n", airr.v_call);
        allele_correction_set_destroy(&ctx);
        sim_config_destroy(&cfg);
        return 1;
    }

    allele_correction_set_destroy(&ctx);
    sim_config_destroy(&cfg);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * Test 10: TSV output
 * ═══════════════════════════════════════════════════════════════ */

static int test_tsv_output(void) {
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    Allele v = make_allele("V1", "AAAAAA", 3, SEG_V);
    Allele j = make_allele("J1", "TTTTTT", 2, SEG_J);
    rec.v_allele = &v;
    rec.j_allele = &j;
    rec.productive = true;

    aseq_append_segment(&seq, "AAAAAA", 6, SEG_V, 0, 3);
    aseq_append_np(&seq, "CC", 2, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "TTTTTT", 6, SEG_J, 0, 2);

    AirrRecord airr;
    SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    /* Write to temp file and verify it's parseable */
    FILE *fp = tmpfile();
    if (!fp) return 1;

    int ncols = airr_write_tsv_header(fp);
    if (ncols <= 0) { fclose(fp); return 1; }

    airr_write_tsv_row(fp, &airr);

    /* Read back and verify header has tabs */
    fseek(fp, 0, SEEK_SET);
    char header[4096];
    if (!fgets(header, sizeof(header), fp)) { fclose(fp); return 1; }

    /* Count tabs in header — should be ncols-1 */
    int tabs = 0;
    for (char *p = header; *p; p++) if (*p == '\t') tabs++;
    if (tabs != ncols - 1) {
        fprintf(stderr, "Expected %d tabs in header, got %d\n",
                ncols - 1, tabs);
        fclose(fp);
        return 1;
    }

    /* Read data row */
    char row[8192];
    if (!fgets(row, sizeof(row), fp)) { fclose(fp); return 1; }

    /* Data row should have same number of tabs */
    int row_tabs = 0;
    for (char *p = row; *p; p++) if (*p == '\t') row_tabs++;
    if (row_tabs != tabs) {
        fprintf(stderr, "Expected %d tabs in row, got %d\n", tabs, row_tabs);
        fclose(fp);
        return 1;
    }

    /* First field should be the sequence */
    if (strncmp(row, "AAAAAACCTTTTTT", 14) != 0) {
        fprintf(stderr, "Row doesn't start with expected sequence\n");
        fclose(fp);
        return 1;
    }

    fclose(fp);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * Test 11: germline alignment has correct NP handling
 * ═══════════════════════════════════════════════════════════════ */

static int test_germline_alignment(void) {
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    Allele v = make_allele("V1", "AACCC", 3, SEG_V);
    Allele j = make_allele("J1", "GGTTT", 2, SEG_J);
    rec.v_allele = &v;
    rec.j_allele = &j;

    aseq_append_segment(&seq, "AACCC", 5, SEG_V, 0, 3);
    aseq_append_np(&seq, "TT", 2, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "GGTTT", 5, SEG_J, 0, 2);

    /* Mutate V pos 1 from A to G */
    aseq_mutate(&seq, seq.head->next, 'G', NUC_FLAG_MUTATED);

    AirrRecord airr;
    SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    /* Current sequence should reflect mutation */
    if (airr.sequence[1] != 'G') {
        fprintf(stderr, "Expected 'G' at pos 1, got '%c'\n", airr.sequence[1]);
        return 1;
    }

    /* Germline alignment should have original base */
    if (airr.germline_alignment[1] != 'A') {
        fprintf(stderr, "Expected germline 'A' at pos 1, got '%c'\n",
                airr.germline_alignment[1]);
        return 1;
    }

    /* NP positions should be N in germline */
    if (airr.germline_alignment[5] != 'N' ||
        airr.germline_alignment[6] != 'N') {
        fprintf(stderr, "Expected 'N' at NP positions in germline\n");
        return 1;
    }

    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * Test 12: short-D triggers "Short-D" call
 * ═══════════════════════════════════════════════════════════════ */

static int test_short_d_in_serialize(void) {
    SimConfig cfg;
    sim_config_init(&cfg, CHAIN_IGH);

    /* Create two D alleles with shared substrings (threshold > 1) */
    Allele d1 = make_allele("D1", "AATTGGCCAATT", 0, SEG_D);
    Allele d2 = make_allele("D2", "AATTGGCCTTAA", 0, SEG_D);
    allele_pool_add(&cfg.d_alleles, &d1);
    allele_pool_add(&cfg.d_alleles, &d2);

    Allele v1 = make_allele("V1", "AAAAAA", 3, SEG_V);
    Allele j1 = make_allele("J1", "TTTTTT", 2, SEG_J);
    allele_pool_add(&cfg.v_alleles, &v1);
    allele_pool_add(&cfg.j_alleles, &j1);

    AlleleCorrectionSet ctx = allele_correction_set_build(&cfg);

    /* Build sequence with D trimmed to just 1 base (< threshold) */
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    rec.v_allele = &v1;
    rec.d_allele = &d1;
    rec.j_allele = &j1;
    rec.d_trim_5 = 5;
    rec.d_trim_3 = 6;

    aseq_append_segment(&seq, "AAAAAA", 6, SEG_V, 0, 3);
    aseq_append_np(&seq, "TT", 2, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "A", 1, SEG_D, 5, -1);
    aseq_append_np(&seq, "AA", 2, SEG_NP2, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "TTTTTT", 6, SEG_J, 0, 2);

    AirrRecord airr;
    airr_serialize(&seq, &rec, &cfg, &ctx, &airr);

    /* D call should be "Short-D" since 1 bp < threshold */
    if (strcmp(airr.d_call, "Short-D") != 0) {
        fprintf(stderr, "Expected 'Short-D', got '%s' (threshold=%d)\n",
                airr.d_call, ctx.d_bitmap.short_d_threshold);
        allele_correction_set_destroy(&ctx);
        sim_config_destroy(&cfg);
        return 1;
    }

    allele_correction_set_destroy(&ctx);
    sim_config_destroy(&cfg);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * Test 13: C allele propagation
 * ═══════════════════════════════════════════════════════════════ */

static int test_c_allele_call(void) {
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    Allele v = make_allele("V1", "AAAAAA", 3, SEG_V);
    Allele j = make_allele("J1", "TTTTTT", 2, SEG_J);
    Allele c = make_allele("IGHM*01", "CCCCCCCC", 0, SEG_C);
    rec.v_allele = &v;
    rec.j_allele = &j;
    rec.c_allele = &c;

    aseq_append_segment(&seq, "AAAAAA", 6, SEG_V, 0, 3);
    aseq_append_segment(&seq, "TTTTTT", 6, SEG_J, 0, 2);
    aseq_append_segment(&seq, "CCCCCCCC", 8, SEG_C, 0, -1);

    AirrRecord airr;
    SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    if (strcmp(airr.c_call, "IGHM*01") != 0) {
        fprintf(stderr, "Expected c_call='IGHM*01', got '%s'\n", airr.c_call);
        return 1;
    }

    return 0;
}

/* ── Main ─────────────────────────────────────────────────────── */

int main(void) {
    int failures = 0;
    int total = 0;

    printf("=== AIRR Serialization Tests ===\n");

    TEST(test_basic_serialize);
    TEST(test_mutation_collection);
    TEST(test_error_collection);
    TEST(test_junction_extraction);
    TEST(test_vj_chain);
    TEST(test_positions_with_ambiguity);
    TEST(test_flags_propagation);
    TEST(test_correction_context_lifecycle);
    TEST(test_serialize_with_corrections);
    TEST(test_tsv_output);
    TEST(test_germline_alignment);
    TEST(test_short_d_in_serialize);
    TEST(test_c_allele_call);

    printf("\n%d/%d tests passed.\n", total - failures, total);
    return failures > 0 ? 1 : 0;
}
