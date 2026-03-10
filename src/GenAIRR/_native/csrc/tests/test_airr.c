/**
 * test_airr.c — Tests for AIRR position derivation.
 *
 * Validates that airr_derive_positions() produces the same results
 * as the Python FixVPosition + FixDPosition + FixJPosition corrections.
 *
 * Key scenarios:
 *   1. No ambiguity: positions match ground truth exactly
 *   2. V 3' ambiguity: V trimmed bases match NP1 start → V extends
 *   3. J 5' ambiguity: J trimmed bases match NP2/NP1 end → J extends
 *   4. D ambiguity: both 5' and 3' overlap with NP regions
 *   5. VJ chain (no D): only V and J ambiguity
 */

#include "genairr/genairr.h"
#include <stdio.h>
#include <string.h>

#define TEST(name) \
    do { \
        printf("  %-55s", #name); \
        int _r = name(); \
        printf("%s\n", _r == 0 ? "PASS" : "FAIL"); \
        if (_r) failures++; \
        total++; \
    } while (0)

/* ── Helper: build a known sequence manually ──────────────────── */

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

/* ── Test: basic positions, no ambiguity ──────────────────────── */

static int test_basic_positions(void) {
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    /* V allele: "AAAAACCCCC" (10bp), anchor at 7, trim 3' by 2 → "AAAAACCC" */
    Allele v = make_allele("V1", "AAAAACCCCC", 7, SEG_V);
    rec.v_allele = &v;
    rec.v_trim_3 = 2;  /* trimmed "CC" */

    /* J allele: "GGGGGTTTTTT" (11bp), anchor at 4, trim 5' by 3 → "GGTTTTTT" */
    Allele j = make_allele("J1", "GGGGGTTTTTT", 4, SEG_J);
    rec.j_allele = &j;
    rec.j_trim_5 = 3;  /* trimmed "GGG" */

    /* Build: V(8bp) + NP1(3bp "ATA") + J(8bp)
     * NP1 "ATA": first 'A' != V trimmed 'C' → no V ambiguity
     *            last 'A' != J trimmed 'G' → no J ambiguity */
    aseq_append_segment(&seq, "AAAAACCC", 8, SEG_V, 0, 7);
    aseq_append_np(&seq, "ATA", 3, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "GGTTTTTT", 8, SEG_J, 3, 4);

    AirrPositions pos;
    airr_derive_positions(&seq, &rec, &pos);

    /* No ambiguity: NP1 = "ATA", V trimmed = "CC", J trimmed = "GGG" — no match */
    if (pos.v_sequence_start != 0) return 1;
    if (pos.v_sequence_end != 8) return 1;    /* ground truth, no extension */
    if (pos.j_sequence_start != 11) return 1; /* ground truth, no extension */
    if (pos.j_sequence_end != 19) return 1;
    if (pos.v_trim_3_adjusted != 2) return 1;
    if (pos.j_trim_5_adjusted != 3) return 1;

    return 0;
}

/* ── Test: V 3' boundary ambiguity ────────────────────────────── */

static int test_v_3prime_ambiguity(void) {
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    /* V allele: "AAAAACCCCC" (10bp), anchor at 7, trim 3' by 3 → "AAAAACCC"
     * Trimmed bases = "CCC" (removed from 3' end) */
    Allele v = make_allele("V1", "AAAAAAACCC", 7, SEG_V);
    rec.v_allele = &v;
    rec.v_trim_3 = 3;

    /* NP1 starts with "CC" — matches first 2 of trimmed "CCC" */
    Allele j = make_allele("J1", "GGGGGTTTTTT", 4, SEG_J);
    rec.j_allele = &j;
    rec.j_trim_5 = 0;

    /* Build: V(7bp) + NP1("CCT" - first 2 match V trim) + J(11bp) */
    aseq_append_segment(&seq, "AAAAAAA", 7, SEG_V, 0, 7);
    aseq_append_np(&seq, "CCT", 3, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "GGGGGTTTTTT", 11, SEG_J, 0, 4);

    AirrPositions pos;
    airr_derive_positions(&seq, &rec, &pos);

    /* V 3' ambiguity: trimmed "CCC", NP1 starts "CC" → 2 bases overlap */
    /* v_sequence_end should extend by 2: 7 → 9 */
    if (pos.v_sequence_end != 9) {
        fprintf(stderr, "Expected v_sequence_end=9, got %d\n", pos.v_sequence_end);
        return 1;
    }
    /* v_germline_end should also extend by 2 */
    if (pos.v_germline_end != 9) {
        fprintf(stderr, "Expected v_germline_end=9, got %d\n", pos.v_germline_end);
        return 1;
    }
    /* Adjusted trim should decrease by 2: 3 → 1 */
    if (pos.v_trim_3_adjusted != 1) {
        fprintf(stderr, "Expected v_trim_3_adjusted=1, got %d\n", pos.v_trim_3_adjusted);
        return 1;
    }

    return 0;
}

/* ── Test: J 5' boundary ambiguity ────────────────────────────── */

static int test_j_5prime_ambiguity(void) {
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    Allele v = make_allele("V1", "AAAAAA", 4, SEG_V);
    rec.v_allele = &v;
    rec.v_trim_3 = 0;

    /* J allele: "GGGTTTTTT" (9bp), trim 5' by 3 → "TTTTTT"
     * Trimmed bases = "GGG" */
    Allele j = make_allele("J1", "GGGTTTTTT", 4, SEG_J);
    rec.j_allele = &j;
    rec.j_trim_5 = 3;

    /* NP1 ends with "GG" — matches last 2 of trimmed "GGG" (compared backwards) */
    aseq_append_segment(&seq, "AAAAAA", 6, SEG_V, 0, 4);
    aseq_append_np(&seq, "AAGG", 4, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "TTTTTT", 6, SEG_J, 3, 4);

    AirrPositions pos;
    airr_derive_positions(&seq, &rec, &pos);

    /* J 5' ambiguity: trimmed "GGG", NP1 ends "GG" → 2 bases overlap */
    /* j_sequence_start should retract by 2: 10 → 8 */
    if (pos.j_sequence_start != 8) {
        fprintf(stderr, "Expected j_sequence_start=8, got %d\n", pos.j_sequence_start);
        return 1;
    }
    if (pos.j_trim_5_adjusted != 1) {
        fprintf(stderr, "Expected j_trim_5_adjusted=1, got %d\n", pos.j_trim_5_adjusted);
        return 1;
    }

    return 0;
}

/* ── Test: D ambiguity on both ends ───────────────────────────── */

static int test_d_ambiguity(void) {
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    Allele v = make_allele("V1", "AAAAAA", 4, SEG_V);
    rec.v_allele = &v;
    rec.v_trim_3 = 0;

    /* D allele: "CCTTTTCC" (8bp), trim 5' by 2, trim 3' by 2 → "TTTT"
     * Trimmed 5' = "CC", trimmed 3' = "CC" */
    Allele d = make_allele("D1", "CCTTTTCC", 0, SEG_D);
    rec.d_allele = &d;
    rec.d_trim_5 = 2;
    rec.d_trim_3 = 2;

    Allele j = make_allele("J1", "GGGGGG", 2, SEG_J);
    rec.j_allele = &j;
    rec.j_trim_5 = 0;

    /* NP1 ends with "C" (matches 1 of trimmed D 5' "CC")
     * NP2 starts with "CC" (matches 2 of trimmed D 3' "CC") */
    aseq_append_segment(&seq, "AAAAAA", 6, SEG_V, 0, 4);
    aseq_append_np(&seq, "AAC", 3, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "TTTT", 4, SEG_D, 2, -1);
    aseq_append_np(&seq, "CCA", 3, SEG_NP2, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "GGGGGG", 6, SEG_J, 0, 2);

    AirrPositions pos;
    airr_derive_positions(&seq, &rec, &pos);

    /* D 5' ambiguity: trimmed "CC", NP1 ends "C" → 1 base overlap */
    /* d_sequence_start should retract by 1 */
    int d_ground_truth_start = 6 + 3;  /* V=6 + NP1=3 = 9 */
    if (pos.d_sequence_start != d_ground_truth_start - 1) {
        fprintf(stderr, "Expected d_start=%d, got %d\n", d_ground_truth_start - 1, pos.d_sequence_start);
        return 1;
    }
    if (pos.d_trim_5_adjusted != 1) {
        fprintf(stderr, "Expected d_trim_5_adjusted=1, got %d\n", pos.d_trim_5_adjusted);
        return 1;
    }

    /* D 3' ambiguity: trimmed "CC", NP2 starts "CC" → 2 bases overlap */
    int d_ground_truth_end = 9 + 4;  /* 13 */
    if (pos.d_sequence_end != d_ground_truth_end + 2) {
        fprintf(stderr, "Expected d_end=%d, got %d\n", d_ground_truth_end + 2, pos.d_sequence_end);
        return 1;
    }
    if (pos.d_trim_3_adjusted != 0) {
        fprintf(stderr, "Expected d_trim_3_adjusted=0, got %d\n", pos.d_trim_3_adjusted);
        return 1;
    }

    return 0;
}

/* ── Test: junction spans V-anchor to J-anchor ────────────────── */

static int test_junction_from_anchors(void) {
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    /* V with anchor at germline pos 4, J with anchor at germline pos 2 */
    Allele v = make_allele("V1", "AAAAACCCCC", 4, SEG_V);
    rec.v_allele = &v;
    rec.v_trim_3 = 0;

    Allele j = make_allele("J1", "GGGTTTTTT", 2, SEG_J);
    rec.j_allele = &j;
    rec.j_trim_5 = 0;

    aseq_append_segment(&seq, "AAAAACCCCC", 10, SEG_V, 0, 4);
    aseq_append_np(&seq, "AT", 2, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "GGGTTTTTT", 9, SEG_J, 0, 2);

    AirrPositions pos;
    airr_derive_positions(&seq, &rec, &pos);

    /* V anchor is at node with germline_pos=4, which is position 4 in the sequence */
    if (pos.junction_start != 4) {
        fprintf(stderr, "Expected junction_start=4, got %d\n", pos.junction_start);
        return 1;
    }
    /* J anchor is at germline_pos=2, which is at seq pos 10+2+2 = 14.
     * AIRR convention: junction includes full W/F codon → end = 14+3 = 17. */
    if (pos.junction_end != 17) {
        fprintf(stderr, "Expected junction_end=17, got %d\n", pos.junction_end);
        return 1;
    }
    if (pos.junction_length != 13) {
        fprintf(stderr, "Expected junction_length=13, got %d\n", pos.junction_length);
        return 1;
    }

    return 0;
}

/* ── Test: indel doesn't break positions ──────────────────────── */

static int test_positions_after_indel(void) {
    ASeq seq;
    aseq_init(&seq);
    SimRecord rec;
    sim_record_init(&rec);

    Allele v = make_allele("V1", "AAAAAA", 3, SEG_V);
    rec.v_allele = &v;
    rec.v_trim_3 = 0;

    Allele j = make_allele("J1", "TTTTTT", 2, SEG_J);
    rec.j_allele = &j;
    rec.j_trim_5 = 0;

    aseq_append_segment(&seq, "AAAAAA", 6, SEG_V, 0, 3);
    aseq_append_np(&seq, "CC", 2, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "TTTTTT", 6, SEG_J, 0, 2);

    /* Insert an indel in V segment — this would break Python's position tracking */
    Nuc *v_third = seq.head->next->next;
    aseq_insert_after(&seq, v_third, 'G', SEG_V, NUC_FLAG_INDEL_INS);
    /* Sequence is now: AAAGAAA CC TTTTTT (15 nodes) */

    AirrPositions pos;
    airr_derive_positions(&seq, &rec, &pos);

    /* V is now 7 nodes (6 original + 1 inserted) */
    if (pos.v_sequence_start != 0) return 1;
    if (pos.v_sequence_end != 7) {
        fprintf(stderr, "Expected v_end=7 after indel, got %d\n", pos.v_sequence_end);
        return 1;
    }

    /* J should start at position 9 (V=7 + NP1=2) */
    if (pos.j_sequence_start != 9) {
        fprintf(stderr, "Expected j_start=9 after indel, got %d\n", pos.j_sequence_start);
        return 1;
    }

    /* Junction should still be found via anchors */
    if (pos.junction_start < 0) {
        fprintf(stderr, "Junction not found after indel\n");
        return 1;
    }

    return 0;
}

/* ── Test: Short-D check ──────────────────────────────────────── */

static int test_short_d(void) {
    ASeq seq;
    aseq_init(&seq);

    /* D segment with only 2 bases */
    aseq_append_segment(&seq, "AAAAAA", 6, SEG_V, 0, -1);
    aseq_append_np(&seq, "CC", 2, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "GT", 2, SEG_D, 0, -1);
    aseq_append_np(&seq, "AA", 2, SEG_NP2, NUC_FLAG_N_NUCLEOTIDE);
    aseq_append_segment(&seq, "TTTTTT", 6, SEG_J, 0, -1);

    /* Threshold 5: D length 2 < 5 → Short-D */
    if (!(aseq_segment_length(&seq, SEG_D) < 5)) {
        fprintf(stderr, "Expected Short-D for 2bp D with threshold 5\n");
        return 1;
    }

    /* Threshold 1: D length 2 >= 1 → not Short-D */
    if ((aseq_segment_length(&seq, SEG_D) < 1)) {
        fprintf(stderr, "Expected non-Short-D for 2bp D with threshold 1\n");
        return 1;
    }

    return 0;
}

/* ── Main ─────────────────────────────────────────────────────── */

int main(void) {
    int failures = 0;
    int total = 0;

    printf("=== AIRR Derivation Tests ===\n");

    TEST(test_basic_positions);
    TEST(test_v_3prime_ambiguity);
    TEST(test_j_5prime_ambiguity);
    TEST(test_d_ambiguity);
    TEST(test_junction_from_anchors);
    TEST(test_positions_after_indel);
    TEST(test_short_d);

    printf("\n%d/%d tests passed.\n", total - failures, total);
    return failures > 0 ? 1 : 0;
}
