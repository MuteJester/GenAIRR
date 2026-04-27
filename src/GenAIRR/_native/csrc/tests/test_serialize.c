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

    /* Productivity: airr_serialize derives productive/stop/in-frame from
     * the codon rail (post-T0-4, V-anchored frame). The synthetic test
     * setup is biologically inconsistent (V Cys anchor at germline pos
     * 7 in "AAAAACCCCC" → not a Cys codon; J anchor at germline pos 2
     * but j_trim_5=3 → anchor trimmed off). The serializer correctly
     * reports productive=false for this. The test focuses on
     * SERIALIZATION mechanics (sequence, germline, calls, NP, mutations);
     * productivity-derivation correctness is covered by test_codon_rail
     * (V-anchored frame tests) and test_simulation_integrity. */

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
    aseq_mark_junction(&seq);
    rec.original_junction_length = 13;

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
    rec.is_reverse_complement = true;
    rec.is_contaminant = false;

    aseq_append_segment(&seq, "AAAAAA", 6, SEG_V, 0, 3);
    aseq_append_segment(&seq, "TTTTTT", 6, SEG_J, 0, 2);

    AirrRecord airr;
    SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    /* Test the flags that are still pass-through from rec → airr.
     * `productive`, `stop_codon`, `vj_in_frame`, and `note` are now
     * derived from the codon rail (post-T0-4), not copied from rec. */
    if (!airr.is_reverse_complement) return 1;
    if (airr.is_contaminant) return 1;

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
 * T1-11: TSV output sanitizes embedded tab/newline/CR characters
 *
 * Pre-fix, both writers passed string fields verbatim to %s\t.
 * A tab in `note` (or any other string field) would silently shift
 * the column count and corrupt downstream csv.DictReader / pandas
 * parsing. The fix sanitizes all string fields inline (\t/\n/\r
 * → space) so the column count is always exactly N_COLUMNS.
 * ═══════════════════════════════════════════════════════════════ */

/* Helper: build a minimal valid AirrRecord and snprintf a row, then
 * verify the column-count and absence-of-stray-newlines invariants.
 * Returns the number of tab separators in the produced row. */
static int t1_11_count_tabs_in_row(const AirrRecord *r, char *out_row,
                                    int out_size) {
    int n = airr_snprintf_tsv_row(out_row, out_size, r);
    if (n <= 0) return -1;
    int tabs = 0;
    for (int i = 0; i < n; i++) {
        if (out_row[i] == '\t') tabs++;
        if (out_row[i] == '\n') return -2;   /* never embed \n */
        if (out_row[i] == '\r') return -3;
    }
    return tabs;
}

static int test_t1_11_tab_in_note_does_not_split_columns(void) {
    AirrRecord r = {0};
    /* Hostile note: tabs and newlines that would shift column count. */
    snprintf(r.note, sizeof(r.note),
             "rule failed.\thidden\tcolumns\nand a fake row\n42:G>A");
    /* Minimal required string fields. */
    strcpy(r.sequence, "ACGT");
    strcpy(r.germline_alignment, "ACGT");
    strcpy(r.v_call,        "IGHV1*01");
    strcpy(r.d_call,        "IGHD1*01");
    strcpy(r.j_call,        "IGHJ1*01");
    strcpy(r.c_call,        "");
    strcpy(r.v_call_true,   "IGHV1*01");
    strcpy(r.d_call_true,   "IGHD1*01");
    strcpy(r.j_call_true,   "IGHJ1*01");
    strcpy(r.original_v_allele_name, "IGHV1*01");
    strcpy(r.junction_nt,   "ACGT");
    strcpy(r.junction_aa,   "X");
    strcpy(r.np1_region,    "AC");
    strcpy(r.np2_region,    "GT");

    char row[8192];
    int tabs = t1_11_count_tabs_in_row(&r, row, sizeof(row));
    if (tabs == -2) {
        fprintf(stderr, "Embedded \\n leaked through into TSV row\n");
        return 1;
    }
    if (tabs == -3) {
        fprintf(stderr, "Embedded \\r leaked through into TSV row\n");
        return 1;
    }
    /* 56 columns → exactly 55 tab separators. */
    if (tabs != 55) {
        fprintf(stderr,
                "Expected 55 tabs (56 columns), got %d. Note-field "
                "tabs leaked into the column structure.\n", tabs);
        return 1;
    }
    /* Sanitization placeholders should be visible. */
    if (strstr(row, "rule failed. hidden columns and a fake row 42:G>A") == NULL) {
        fprintf(stderr,
                "Expected sanitized note (tabs/newlines → spaces); "
                "got row: %.200s\n", row);
        return 1;
    }
    return 0;
}

static int test_t1_11_tab_in_allele_name_handled(void) {
    AirrRecord r = {0};
    /* A maliciously crafted GDC could embed a tab in an allele name.
     * The TSV writer must absorb it. */
    snprintf(r.v_call, sizeof(r.v_call), "IGHV1\tBOGUS*01");
    /* Other required strings */
    strcpy(r.sequence, "ACGT");
    strcpy(r.germline_alignment, "ACGT");
    strcpy(r.d_call, "");
    strcpy(r.j_call, "IGHJ1*01");
    strcpy(r.c_call, "");
    strcpy(r.v_call_true, "IGHV1*01");
    strcpy(r.d_call_true, "");
    strcpy(r.j_call_true, "IGHJ1*01");

    char row[8192];
    int tabs = t1_11_count_tabs_in_row(&r, row, sizeof(row));
    if (tabs != 55) {
        fprintf(stderr,
                "Expected 55 tabs after v_call sanitization, got %d\n", tabs);
        return 1;
    }
    return 0;
}

static int test_t1_11_normal_record_unchanged(void) {
    /* Regression: a record with no hostile chars produces a row
     * indistinguishable (modulo whitespace) from the legacy fprintf
     * output. We just assert column count and presence of expected
     * sequence prefix — byte-for-byte format equality is checked
     * implicitly by every other test that round-trips through TSV. */
    ASeq seq; SimRecord rec;
    aseq_init(&seq); sim_record_init(&rec);
    Allele v = make_allele("V1", "AAAAAA", 3, SEG_V);
    Allele j = make_allele("J1", "TTTTTT", 2, SEG_J);
    rec.v_allele = &v; rec.j_allele = &j; rec.productive = true;

    aseq_append_segment(&seq, "AAAAAA", 6, SEG_V, 0, 3);
    aseq_append_np(&seq, "CC", 2, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "TTTTTT", 6, SEG_J, 0, 2);

    AirrRecord airr;
    SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    char row[8192];
    int tabs = t1_11_count_tabs_in_row(&airr, row, sizeof(row));
    if (tabs != 55) {
        fprintf(stderr,
                "Normal record: expected 55 tabs, got %d (regression)\n", tabs);
        aseq_reset(&seq);
        return 1;
    }
    if (strncmp(row, "AAAAAACCTTTTTT", 14) != 0) {
        fprintf(stderr, "Normal record row doesn't start with sequence\n");
        aseq_reset(&seq);
        return 1;
    }
    aseq_reset(&seq);
    return 0;
}

static int test_t1_11_buffer_truncation_returns_minus1(void) {
    /* If the snprintf buffer is too small, the writer must signal
     * failure (-1) rather than emit a partial row. */
    AirrRecord r = {0};
    strcpy(r.sequence, "ACGTACGT");
    strcpy(r.germline_alignment, "ACGTACGT");
    char tiny[16];   /* far too small for 56 columns */
    int n = airr_snprintf_tsv_row(tiny, sizeof(tiny), &r);
    if (n != -1) {
        fprintf(stderr,
                "Tiny buffer: expected -1 (truncated), got %d\n", n);
        return 1;
    }
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * Test 11 (legacy): germline alignment has correct NP handling
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

/* ═══════════════════════════════════════════════════════════════
 * T1-2: junction tracked as a node-flagged region.
 *
 * Build a small VDJ-like sequence, mark the junction at "assembly",
 * then exercise: anchor mutation (region intact), 3' tail trim past
 * W/F (partial), 5' head trim past V Cys (partial), indel inside,
 * contaminant wipe (fully removed).
 * ═══════════════════════════════════════════════════════════════ */

/* Helper: build a V-NP1-J sequence with V anchor at germline_pos=4,
 * J anchor at germline_pos=2. Mark junction. */
static void build_jn_seq(ASeq *seq, SimRecord *rec) {
    aseq_init(seq);
    sim_record_init(rec);
    static Allele v, j;
    v = make_allele("V1", "AAAAACCCCC", 4, SEG_V);
    j = make_allele("J1", "GGGTTTTTT",  2, SEG_J);
    rec->v_allele = &v;
    rec->j_allele = &j;
    aseq_append_segment(seq, "AAAAACCCCC", 10, SEG_V, 0, 4);
    aseq_append_np(seq, "AT", 2, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(seq, "GGGTTTTTT", 9, SEG_J, 0, 2);
    rec->original_junction_length = aseq_mark_junction(seq);
}

static int junction_node_count(const ASeq *seq) {
    int c = 0;
    for (Nuc *n = seq->head; n; n = n->next)
        if (n->flags & NUC_FLAG_JUNCTION) c++;
    return c;
}

/* T1-2 case 1: normal — full junction emitted, original == live. */
static int test_junction_full_emitted(void) {
    ASeq seq; SimRecord rec;
    build_jn_seq(&seq, &rec);
    if (rec.original_junction_length != 13) {
        fprintf(stderr, "Expected original=13, got %d\n",
                rec.original_junction_length);
        return 1;
    }
    if (junction_node_count(&seq) != 13) return 1;

    AirrRecord airr; SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);
    if (airr.junction_length != 13) return 1;
    if (strncmp(airr.junction_nt, "ACCCCCATGGGTT", 13) != 0) return 1;
    /* No truncation note */
    if (strstr(airr.note, "junction") != NULL) {
        fprintf(stderr, "Unexpected note: %s\n", airr.note);
        return 1;
    }
    return 0;
}

/* T1-2 case 2: anchor base mutated (Cys->Ser at V anchor). Region
 * fully intact — coordinates and content emitted normally. */
static int test_junction_anchor_mutated_kept(void) {
    ASeq seq; SimRecord rec;
    build_jn_seq(&seq, &rec);

    /* Change V-anchor base to a different base. Node + flag intact. */
    Nuc *va = aseq_find_anchor(&seq, SEG_V);
    if (!va) return 1;
    char before = va->current;
    aseq_mutate(&seq, va, (before == 'C') ? 'A' : 'C', NUC_FLAG_MUTATED);

    AirrRecord airr; SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    if (airr.junction_length != 13) {
        fprintf(stderr, "Anchor-mutated junction should still be full, got %d\n",
                airr.junction_length);
        return 1;
    }
    /* No truncation marker */
    if (strstr(airr.note, "junction partial") != NULL ||
        strstr(airr.note, "junction fully") != NULL) {
        fprintf(stderr, "Unexpected truncation note: %s\n", airr.note);
        return 1;
    }
    return 0;
}

/* T1-2 case 3: 3' trim removes 1 base of W/F codon. Junction shrinks
 * by 1 nt; truncation note appended. */
static int test_junction_3prime_partial_wf(void) {
    ASeq seq; SimRecord rec;
    build_jn_seq(&seq, &rec);

    /* Delete the LAST junction-flagged node (third base of W/F). */
    Nuc *last_jn = NULL;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->flags & NUC_FLAG_JUNCTION) last_jn = n;
    }
    if (!last_jn) return 1;
    aseq_delete(&seq, last_jn);

    AirrRecord airr; SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    if (airr.junction_length != 12) {
        fprintf(stderr, "Expected junction_length=12 after 1-base trim, got %d\n",
                airr.junction_length);
        return 1;
    }
    if (strstr(airr.note, "junction partial: 12/13 nt") == NULL) {
        fprintf(stderr, "Expected truncation note, got: %s\n", airr.note);
        return 1;
    }
    return 0;
}

/* T1-2 case 4: 5' trim removes the V anchor (junction shrinks at start). */
static int test_junction_5prime_past_v_anchor(void) {
    ASeq seq; SimRecord rec;
    build_jn_seq(&seq, &rec);
    int orig = rec.original_junction_length;

    /* Delete the V anchor node. The first surviving junction-flagged
     * node is the next one. */
    Nuc *va = aseq_find_anchor(&seq, SEG_V);
    if (!va) return 1;
    aseq_delete(&seq, va);

    AirrRecord airr; SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    if (airr.junction_length != orig - 1) {
        fprintf(stderr, "Expected junction_length=%d, got %d\n",
                orig - 1, airr.junction_length);
        return 1;
    }
    if (strstr(airr.note, "junction partial") == NULL) {
        fprintf(stderr, "Expected partial note, got: %s\n", airr.note);
        return 1;
    }
    return 0;
}

/* T1-2 case 5: indel inserts inside junction → flag inherited, span
 * grows by 1. */
static int test_junction_indel_inside_grows(void) {
    ASeq seq; SimRecord rec;
    build_jn_seq(&seq, &rec);
    int orig = rec.original_junction_length;

    /* Find a junction-flagged node mid-span (skip first to ensure both
     * neighbors have JUNCTION flag). */
    Nuc *target = NULL;
    int seen = 0;
    for (Nuc *n = seq.head; n; n = n->next) {
        if (n->flags & NUC_FLAG_JUNCTION) {
            seen++;
            if (seen == 5) { target = n; break; }
        }
    }
    if (!target) return 1;

    /* Insert after target — both target and target->next are junction. */
    Nuc *new_node = aseq_insert_after(&seq, target, 'X', SEG_NP1,
                                       NUC_FLAG_INDEL_INS);
    if (!new_node) return 1;
    if (!(new_node->flags & NUC_FLAG_JUNCTION)) {
        fprintf(stderr, "Inserted indel did not inherit JUNCTION flag\n");
        return 1;
    }

    AirrRecord airr; SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    if (airr.junction_length != orig + 1) {
        fprintf(stderr, "Expected junction_length=%d after indel, got %d\n",
                orig + 1, airr.junction_length);
        return 1;
    }
    /* Note signals divergence from original */
    if (strstr(airr.note, "junction partial") == NULL) {
        fprintf(stderr, "Expected partial note (live != original), got: %s\n",
                airr.note);
        return 1;
    }
    return 0;
}

/* T1-2 case 6: insert OUTSIDE junction does NOT inherit the flag. */
static int test_junction_insert_outside_no_inherit(void) {
    ASeq seq; SimRecord rec;
    build_jn_seq(&seq, &rec);

    /* Insert before head (first V base, before V anchor → outside junction). */
    Nuc *first = seq.head;
    if (first->flags & NUC_FLAG_JUNCTION) return 1;  /* sanity */
    Nuc *new_node = aseq_insert_before(&seq, first, 'N', SEG_V,
                                        NUC_FLAG_INDEL_INS);
    if (!new_node) return 1;
    if (new_node->flags & NUC_FLAG_JUNCTION) {
        fprintf(stderr, "Outside-junction insert wrongly inherited JUNCTION\n");
        return 1;
    }
    return 0;
}

/* T1-2 case 7: every junction-flagged node deleted → fully removed. */
static int test_junction_all_removed(void) {
    ASeq seq; SimRecord rec;
    build_jn_seq(&seq, &rec);

    /* Delete every junction-flagged node. */
    Nuc *n = seq.head;
    while (n) {
        Nuc *next = n->next;
        if (n->flags & NUC_FLAG_JUNCTION) aseq_delete(&seq, n);
        n = next;
    }

    AirrRecord airr; SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    if (airr.junction_length != 0) return 1;
    if (airr.junction_start != 0 || airr.junction_end != 0) return 1;
    if (airr.junction_nt[0] != '\0') return 1;
    if (airr.junction_aa[0] != '\0') return 1;
    if (strstr(airr.note, "junction fully removed") == NULL) {
        fprintf(stderr, "Expected 'junction fully removed', got: %s\n", airr.note);
        return 1;
    }
    return 0;
}

/* T1-2 case 8: no anchors → no junction marked, original=0, no
 * spurious truncation note. */
static int test_junction_no_anchors(void) {
    ASeq seq; SimRecord rec;
    aseq_init(&seq); sim_record_init(&rec);
    static Allele v, j;
    v = make_allele("V1", "AAAAA", -1, SEG_V);     /* anchorless */
    j = make_allele("J1", "TTTTT", -1, SEG_J);
    rec.v_allele = &v;
    rec.j_allele = &j;
    aseq_append_segment(&seq, "AAAAA", 5, SEG_V, 0, -1);
    aseq_append_segment(&seq, "TTTTT", 5, SEG_J, 0, -1);
    rec.original_junction_length = aseq_mark_junction(&seq);

    if (rec.original_junction_length != 0) return 1;

    AirrRecord airr; SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    if (airr.junction_length != 0) return 1;
    /* No truncation note (we never had a junction to truncate). */
    if (strstr(airr.note, "junction partial") != NULL ||
        strstr(airr.note, "junction fully") != NULL) {
        return 1;
    }
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * T1-4: when D is missing from the live sequence (e.g., fully
 * deleted by indels) but rec->d_allele is non-NULL with non-zero
 * d_trim_5/3, the D-side adjustments in airr_derive_positions used
 * to be skipped, leaving d_trim_*_adjusted at memset 0. That made
 * d{5,3}_ambig in patch_boundary_extensions / collect_boundary_
 * mutations equal rec->d_trim_X (full trim), and the patch/scan
 * fired at d_sequence_start = 0 (memset default), corrupting V's
 * region with D's germline and emitting spurious mutation entries.
 *
 * The fix: initialize d_trim_*_adjusted to rec->d_trim_X
 * unconditionally, before the has_d block.
 * ═══════════════════════════════════════════════════════════════ */

static int test_t1_4_no_spurious_mut_when_d_missing(void) {
    ASeq seq; SimRecord rec;
    aseq_init(&seq); sim_record_init(&rec);

    static Allele v, d, j;
    v = make_allele("V1", "ACGTACGTACGTACGTACGTACGTACGTAC", 4, SEG_V);  /* 30bp */
    d = make_allele("D1", "GGGCCCAAATTTGGG",                  0, SEG_D);  /* 15bp */
    j = make_allele("J1", "TTTAAATTT",                        2, SEG_J);  /* 9bp */
    rec.v_allele = &v;
    rec.d_allele = &d;
    rec.j_allele = &j;
    /* Critical: simulator-side trims are non-zero, but D will be
     * absent from the live sequence (we don't append it). This is
     * the post-indel-D-wipeout scenario. */
    rec.v_trim_3 = 0;
    rec.d_trim_5 = 5;     /* would normally extend into NP1 */
    rec.d_trim_3 = 3;
    rec.j_trim_5 = 0;

    /* Build sequence with V, NP1, NP2, J — NO D. */
    aseq_append_segment(&seq, v.seq, v.length, SEG_V, 0, 4);
    aseq_append_np(&seq, "ATGC", 4, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_np(&seq, "TGCA", 4, SEG_NP2, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, j.seq, j.length, SEG_J, 0, 2);
    aseq_mark_junction(&seq);

    AirrRecord airr; SimConfig cfg = {0};
    airr_serialize(&seq, &rec, &cfg, NULL, &airr);

    /* No SHM, no quality errors, no PCR errors → 0 mutations. */
    if (airr.n_mutations != 0) {
        fprintf(stderr,
                "Expected n_mutations=0 with D missing, got %d. "
                "Mutations field: %s\n",
                airr.n_mutations, airr.mutations);
        return 1;
    }
    if (airr.mutations[0] != '\0') {
        fprintf(stderr, "Expected empty mutations string, got: %s\n",
                airr.mutations);
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
    TEST(test_t1_11_tab_in_note_does_not_split_columns);
    TEST(test_t1_11_tab_in_allele_name_handled);
    TEST(test_t1_11_normal_record_unchanged);
    TEST(test_t1_11_buffer_truncation_returns_minus1);
    TEST(test_germline_alignment);
    TEST(test_short_d_in_serialize);
    TEST(test_c_allele_call);

    printf("\n--- T1-2 junction-as-tracked-region ---\n");
    TEST(test_junction_full_emitted);
    TEST(test_junction_anchor_mutated_kept);
    TEST(test_junction_3prime_partial_wf);
    TEST(test_junction_5prime_past_v_anchor);
    TEST(test_junction_indel_inside_grows);
    TEST(test_junction_insert_outside_no_inherit);
    TEST(test_junction_all_removed);
    TEST(test_junction_no_anchors);

    printf("\n--- T1-4 boundary mutations with missing D ---\n");
    TEST(test_t1_4_no_spurious_mut_when_d_missing);

    printf("\n%d/%d tests passed.\n", total - failures, total);
    return failures > 0 ? 1 : 0;
}
