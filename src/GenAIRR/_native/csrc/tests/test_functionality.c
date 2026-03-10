/**
 * test_functionality.c — Tests for the rule-based FunctionalityValidator.
 *
 * Validates all 5 rules individually and the full assessment flow,
 * matching the Python FunctionalityValidator behavior exactly.
 *
 * Also tests codon rail integration with the functionality rules.
 */

#include "genairr/functionality.h"
#include "genairr/pipeline.h"
#include "genairr/aseq.h"
#include "genairr/codon.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

static int tests_passed = 0;
static int tests_total  = 0;

#define RUN_TEST(fn)  do {                                          \
    tests_total++;                                                  \
    printf("  %-50s ", #fn);                                        \
    fn();                                                           \
    printf("PASS\n");                                               \
    tests_passed++;                                                 \
} while (0)

/* ── Helpers ─────────────────────────────────────────────────── */

/**
 * Build a minimal ASeq with V + NP1 + J segments for testing.
 * V anchor is placed at v_anchor_pos within V (germline_pos).
 * J anchor is placed at j_anchor_pos within J (germline_pos).
 *
 * IMPORTANT: j_anchor_pos must leave room for 2 more bases after it
 * in the J segment, since the junction includes the full 3-base W/F
 * codon starting at the anchor.
 */
static void build_test_seq(ASeq *seq,
                           const char *v_seq, int v_len, int v_anchor_pos,
                           const char *np1, int np1_len,
                           const char *j_seq, int j_len, int j_anchor_pos) {
    aseq_init(seq);
    if (v_len > 0)
        aseq_append_segment(seq, v_seq, v_len, SEG_V, 0, v_anchor_pos);
    if (np1_len > 0)
        aseq_append_np(seq, np1, np1_len, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    if (j_len > 0)
        aseq_append_segment(seq, j_seq, j_len, SEG_J, 0, j_anchor_pos);
}

/**
 * Build a VDJ sequence (with D and NP2).
 */
static void build_test_seq_vdj(ASeq *seq,
                               const char *v_seq, int v_len, int v_anchor_pos,
                               const char *np1, int np1_len,
                               const char *d_seq, int d_len,
                               const char *np2, int np2_len,
                               const char *j_seq, int j_len, int j_anchor_pos) {
    aseq_init(seq);
    if (v_len > 0)
        aseq_append_segment(seq, v_seq, v_len, SEG_V, 0, v_anchor_pos);
    if (np1_len > 0)
        aseq_append_np(seq, np1, np1_len, SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
    if (d_len > 0)
        aseq_append_segment(seq, d_seq, d_len, SEG_D, 0, -1);
    if (np2_len > 0)
        aseq_append_np(seq, np2, np2_len, SEG_NP2, NUC_FLAG_N_NUCLEOTIDE);
    if (j_len > 0)
        aseq_append_segment(seq, j_seq, j_len, SEG_J, 0, j_anchor_pos);
}

/**
 * Initialize a FuncContext from an ASeq, exactly matching func_context_init
 * in functionality.c. This is a test-visible mirror.
 *
 * Junction boundary: j_anchor_pos + 3 (inclusive of full W/F codon).
 */
static void init_test_context(FuncContext *ctx, const ASeq *seq) {
    memset(ctx, 0, sizeof(*ctx));
    ctx->seq = seq;

    int pos = 0;
    for (Nuc *n = seq->head; n; n = n->next) {
        if ((n->flags & NUC_FLAG_ANCHOR) && n->segment == SEG_V && !ctx->v_anchor) {
            ctx->v_anchor = n;
            ctx->junction_start = pos;
        }
        if ((n->flags & NUC_FLAG_ANCHOR) && n->segment == SEG_J && !ctx->j_anchor) {
            ctx->j_anchor = n;
            ctx->junction_end = pos + 3;
        }
        pos++;
    }

    ctx->anchors_present = (ctx->v_anchor != NULL && ctx->j_anchor != NULL);
    if (ctx->anchors_present) {
        ctx->junction_len = ctx->junction_end - ctx->junction_start;
    }
}

/*
 * Standard test sequences used throughout:
 *
 * Productive (junction_len=12, 4 codons: C G P W):
 *   V = "ATGATGTGT" (9bp), anchor@6 → assembled pos 6, codon TGT=C
 *   NP1 = "GGG" (3bp)
 *   J = "CCCTGG" (6bp), anchor@3 → assembled pos 15, codon TGG=W
 *   junction_start=6, junction_end=18, junction_len=12
 *   Full sequence: ATG ATG TGT | GGG | CCC TGG → reading frame codons
 *   Junction: TGT GGG CCC TGG → C G P W
 *
 * M-start (junction starts with M instead of C):
 *   V = "ATGATGATG" (9bp), anchor@6 → codon ATG=M
 *
 * Stop codon in V:
 *   V = "TAATGTTGT" (9bp) → TAA at pos 0 = stop
 *
 * F-anchor:
 *   J = "CCCTTT" (6bp), anchor@3 → codon TTT=F
 *
 * G-anchor (non-W/F):
 *   J = "CCCGGG" (6bp), anchor@3 → codon GGG=G
 */

/* ═══════════════════════════════════════════════════════════════
 * Rule 1: Stop Codon
 * ═══════════════════════════════════════════════════════════════ */

static void test_stop_codon_absent(void) {
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGTGT", 9, 6,
                   "GGG", 3,
                   "CCCTGG", 6, 3);

    FuncContext ctx;
    init_test_context(&ctx, &seq);

    FuncRuleResult r = rule_stop_codon(&ctx);
    assert(r.passed);
    assert(!ctx.stop_codon_found);

    aseq_reset(&seq);
}

static void test_stop_codon_in_v_region(void) {
    /* TAA at position 0-2 in V = stop codon */
    ASeq seq;
    build_test_seq(&seq,
                   "TAATGTTGT", 9, 6,
                   "GGG", 3,
                   "CCCTGG", 6, 3);

    FuncContext ctx;
    init_test_context(&ctx, &seq);

    FuncRuleResult r = rule_stop_codon(&ctx);
    assert(!r.passed);
    assert(ctx.stop_codon_found);

    aseq_reset(&seq);
}

static void test_stop_codon_in_junction(void) {
    /* Stop codon TGA within the junction */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGTGA", 9, 6,   /* TGA at pos 6 = stop in junction */
                   "GGG", 3,
                   "CCCTGG", 6, 3);

    FuncContext ctx;
    init_test_context(&ctx, &seq);

    FuncRuleResult r = rule_stop_codon(&ctx);
    assert(!r.passed);
    assert(ctx.stop_codon_found);

    aseq_reset(&seq);
}

static void test_stop_codon_scans_entire_sequence(void) {
    /* Stop codon only in J region (after junction) */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGTGT", 9, 6,
                   "GGG", 3,
                   "CCCTGA", 6, 0);    /* TGA at J pos 3 (assembled pos 15) = stop */

    FuncContext ctx;
    init_test_context(&ctx, &seq);

    FuncRuleResult r = rule_stop_codon(&ctx);
    assert(!r.passed);   /* Must catch stop codon anywhere */
    assert(ctx.stop_codon_found);

    aseq_reset(&seq);
}

/* ═══════════════════════════════════════════════════════════════
 * Rule 2: Frame Alignment
 * ═══════════════════════════════════════════════════════════════ */

static void test_frame_aligned(void) {
    /* V anchor@6 (6%3=0), J anchor@3 → assembled pos 15 → junction_end=18 (18%3=0)
     * junction_len=12 (12%3=0). All divisible by 3. */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGATG", 9, 6,
                   "GGG", 3,
                   "CCCTGG", 6, 3);

    FuncContext ctx;
    init_test_context(&ctx, &seq);

    assert(ctx.junction_start == 6);
    assert(ctx.junction_start % 3 == 0);
    assert(ctx.junction_end % 3 == 0);
    assert(ctx.junction_len % 3 == 0);

    FuncRuleResult r = rule_frame_alignment(&ctx);
    assert(r.passed);

    aseq_reset(&seq);
}

static void test_frame_misaligned(void) {
    /* V anchor at position 7 (not divisible by 3) */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGATGA", 10, 7,   /* V: 10bp, anchor at pos 7 */
                   "GGG", 3,
                   "CCCTGG", 6, 3);

    FuncContext ctx;
    init_test_context(&ctx, &seq);

    assert(ctx.junction_start == 7);
    assert(ctx.junction_start % 3 != 0);

    FuncRuleResult r = rule_frame_alignment(&ctx);
    assert(!r.passed);

    aseq_reset(&seq);
}

/* ═══════════════════════════════════════════════════════════════
 * Rule 3: Junction Translatable
 * ═══════════════════════════════════════════════════════════════ */

static void test_junction_translatable_ok(void) {
    /* Junction length = 12 (divisible by 3), should translate */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGTGT", 9, 6,   /* anchor at pos 6, TGT=C */
                   "GGG", 3,
                   "CCCTGG", 6, 3);     /* anchor at pos 3, TGG=W */

    FuncContext ctx;
    init_test_context(&ctx, &seq);

    assert(ctx.junction_len == 12);
    assert(ctx.junction_len % 3 == 0);

    FuncRuleResult r = rule_junction_translatable(&ctx);
    assert(r.passed);
    assert(ctx.junction_aa_len > 0);

    aseq_reset(&seq);
}

static void test_junction_translatable_bad_length(void) {
    /* NP1=4bp breaks frame → junction_len not divisible by 3 */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGTGT", 9, 6,
                   "GGGG", 4,           /* 4bp NP1 breaks frame */
                   "CCCTGG", 6, 3);

    FuncContext ctx;
    init_test_context(&ctx, &seq);

    /* junction_start=6, J starts at 13, J anchor at 16, junction_end=19 */
    /* junction_len = 19-6 = 13 (not divisible by 3) */
    assert(ctx.junction_len % 3 != 0);

    FuncRuleResult r = rule_junction_translatable(&ctx);
    assert(!r.passed);
    assert(ctx.junction_aa_len == 0);  /* no translation performed */

    aseq_reset(&seq);
}

/* ═══════════════════════════════════════════════════════════════
 * Rule 4: Conserved Cysteine
 * ═══════════════════════════════════════════════════════════════ */

static void test_conserved_cysteine_present(void) {
    /* Junction starts with TGT → translates to C (Cysteine) */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGTGT", 9, 6,   /* V: anchor at 6, codon TGT=C */
                   "GGG", 3,
                   "CCCTGG", 6, 3);     /* J: anchor at 3, codon TGG=W */

    FuncContext ctx;
    init_test_context(&ctx, &seq);

    /* First translate */
    rule_junction_translatable(&ctx);
    assert(ctx.junction_aa_len > 0);
    assert(ctx.junction_aa[0] == 'C');

    FuncRuleResult r = rule_conserved_cysteine(&ctx);
    assert(r.passed);

    aseq_reset(&seq);
}

static void test_conserved_cysteine_absent(void) {
    /* Junction starts with ATG → translates to M (Methionine), not C */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGATG", 9, 6,   /* V: anchor at 6, codon ATG=M */
                   "GGG", 3,
                   "CCCTGG", 6, 3);

    FuncContext ctx;
    init_test_context(&ctx, &seq);

    rule_junction_translatable(&ctx);
    assert(ctx.junction_aa[0] == 'M');

    FuncRuleResult r = rule_conserved_cysteine(&ctx);
    assert(!r.passed);
    assert(strstr(r.note, "C not present") != NULL);

    aseq_reset(&seq);
}

/* ═══════════════════════════════════════════════════════════════
 * Rule 5: Conserved Anchor (W/F)
 * ═══════════════════════════════════════════════════════════════ */

static void test_conserved_anchor_w(void) {
    /* Last junction codon translates to W (TGG) */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGTGT", 9, 6,
                   "GGG", 3,
                   "CCCTGG", 6, 3);     /* J anchor → TGG = W */

    FuncContext ctx;
    init_test_context(&ctx, &seq);

    rule_junction_translatable(&ctx);

    FuncRuleResult r = rule_conserved_anchor(&ctx);
    assert(r.passed);

    aseq_reset(&seq);
}

static void test_conserved_anchor_f(void) {
    /* Last junction codon translates to F (TTT) */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGTGT", 9, 6,
                   "GGG", 3,
                   "CCCTTT", 6, 3);     /* J anchor → TTT = F */

    FuncContext ctx;
    init_test_context(&ctx, &seq);

    rule_junction_translatable(&ctx);

    FuncRuleResult r = rule_conserved_anchor(&ctx);
    assert(r.passed);

    aseq_reset(&seq);
}

static void test_conserved_anchor_absent(void) {
    /* Last junction codon = GGG → G (Glycine), not W/F */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGTGT", 9, 6,
                   "GGG", 3,
                   "CCCGGG", 6, 3);     /* J anchor → GGG = G */

    FuncContext ctx;
    init_test_context(&ctx, &seq);

    rule_junction_translatable(&ctx);

    FuncRuleResult r = rule_conserved_anchor(&ctx);
    assert(!r.passed);
    assert(strstr(r.note, "W/F") != NULL);

    aseq_reset(&seq);
}

/* ═══════════════════════════════════════════════════════════════
 * Full assessment integration
 * ═══════════════════════════════════════════════════════════════ */

static void test_assess_productive(void) {
    /* Fully productive: no stop codon, in-frame, C at start, W at end */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGTGT", 9, 6,
                   "GGG", 3,
                   "CCCTGG", 6, 3);

    SimRecord rec;
    sim_record_init(&rec);

    FunctionalityValidator v = functionality_validator_default();
    functionality_assess(&v, &seq, &rec);

    assert(rec.productive);
    assert(!rec.stop_codon);
    assert(rec.vj_in_frame);
    assert(rec.note[0] == '\0');

    aseq_reset(&seq);
}

static void test_assess_stop_codon_is_fatal(void) {
    /* Stop codon → productive=false, stop_codon=true, vj_in_frame=false */
    ASeq seq;
    /* TAA at pos 0-2 = stop codon */
    build_test_seq(&seq,
                   "TAATGTTGT", 9, 6,
                   "GGG", 3,
                   "CCCTGG", 6, 3);

    SimRecord rec;
    sim_record_init(&rec);

    FunctionalityValidator v = functionality_validator_default();
    functionality_assess(&v, &seq, &rec);

    assert(!rec.productive);
    assert(rec.stop_codon);
    assert(!rec.vj_in_frame);  /* stop codon contributes to vj_in_frame */

    aseq_reset(&seq);
}

static void test_assess_out_of_frame(void) {
    /* No stop codon, but junction length not divisible by 3 */
    ASeq seq;
    /* V=9bp anchor@6, NP1=4bp (breaks frame), J=6bp anchor@3 */
    build_test_seq(&seq,
                   "ATGATGTGT", 9, 6,
                   "GGGG", 4,           /* 4bp NP1 → out of frame */
                   "CCCTGG", 6, 3);

    SimRecord rec;
    sim_record_init(&rec);

    FunctionalityValidator v = functionality_validator_default();
    functionality_assess(&v, &seq, &rec);

    assert(!rec.productive);
    assert(!rec.stop_codon);
    assert(!rec.vj_in_frame);

    aseq_reset(&seq);
}

static void test_assess_missing_cysteine(void) {
    /* In frame, no stop codon, but first junction AA is not C */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGATG", 9, 6,   /* ATG=M, not C */
                   "GGG", 3,
                   "CCCTGG", 6, 3);

    SimRecord rec;
    sim_record_init(&rec);

    FunctionalityValidator v = functionality_validator_default();
    functionality_assess(&v, &seq, &rec);

    assert(!rec.productive);       /* NOT productive — missing C */
    assert(!rec.stop_codon);
    assert(rec.vj_in_frame);       /* frame is OK */
    assert(strstr(rec.note, "C not present") != NULL);

    aseq_reset(&seq);
}

static void test_assess_missing_anchor_wf(void) {
    /* In frame, no stop, C present, but last AA is not W/F */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGTGT", 9, 6,
                   "GGG", 3,
                   "CCCGGG", 6, 3);     /* GGG = G, not W/F */

    SimRecord rec;
    sim_record_init(&rec);

    FunctionalityValidator v = functionality_validator_default();
    functionality_assess(&v, &seq, &rec);

    assert(!rec.productive);
    assert(!rec.stop_codon);
    assert(rec.vj_in_frame);
    assert(strstr(rec.note, "W/F") != NULL);

    aseq_reset(&seq);
}

static void test_assess_missing_anchors(void) {
    /* No V or J anchor → immediate failure */
    ASeq seq;
    aseq_init(&seq);
    aseq_append_segment(&seq, "ATGATG", 6, SEG_V, 0, -1);

    SimRecord rec;
    sim_record_init(&rec);

    FunctionalityValidator v = functionality_validator_default();
    functionality_assess(&v, &seq, &rec);

    assert(!rec.productive);
    assert(!rec.stop_codon);
    assert(!rec.vj_in_frame);
    assert(strstr(rec.note, "Missing anchor") != NULL);

    aseq_reset(&seq);
}

static void test_assess_vdj_productive(void) {
    /* Full VDJ sequence that is productive */
    ASeq seq;
    /* V=9bp anchor@6, NP1=3bp, D=6bp, NP2=3bp, J=6bp anchor@3 */
    /* Total=27, junction_start=6, junction_end=27, junction_len=21 (7 codons) */
    build_test_seq_vdj(&seq,
                       "ATGATGTGT", 9, 6,    /* V anchor@6 → TGT=C */
                       "GGG", 3,             /* NP1 */
                       "ATGATG", 6,          /* D */
                       "GGG", 3,             /* NP2 */
                       "CCCTGG", 6, 3);      /* J anchor@3 → TGG=W */

    SimRecord rec;
    sim_record_init(&rec);

    FunctionalityValidator v = functionality_validator_default();
    functionality_assess(&v, &seq, &rec);

    assert(rec.productive);
    assert(!rec.stop_codon);
    assert(rec.vj_in_frame);

    aseq_reset(&seq);
}

/* ═══════════════════════════════════════════════════════════════
 * Codon rail integration
 * ═══════════════════════════════════════════════════════════════ */

static void test_codon_rail_stop_codon_fast_path(void) {
    /* Build codon rail, then use rule_stop_codon which should use O(1) path */
    ASeq seq;
    build_test_seq(&seq,
                   "TAATGTTGT", 9, 6,   /* TAA at pos 0 = stop */
                   "GGG", 3,
                   "CCCTGG", 6, 3);

    aseq_build_codon_rail(&seq);
    assert(seq.codon_rail_valid);
    assert(seq.n_stop_codons > 0);

    FuncContext ctx;
    init_test_context(&ctx, &seq);

    FuncRuleResult r = rule_stop_codon(&ctx);
    assert(!r.passed);
    assert(ctx.stop_codon_found);

    aseq_reset(&seq);
}

static void test_codon_rail_no_stop_fast_path(void) {
    /* Codon rail with no stop codons */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGTGT", 9, 6,
                   "GGG", 3,
                   "CCCTGG", 6, 3);

    aseq_build_codon_rail(&seq);
    assert(seq.codon_rail_valid);
    assert(seq.n_stop_codons == 0);

    FuncContext ctx;
    init_test_context(&ctx, &seq);

    FuncRuleResult r = rule_stop_codon(&ctx);
    assert(r.passed);

    aseq_reset(&seq);
}

static void test_codon_rail_junction_translation(void) {
    /* Codon rail provides amino acids directly for junction translation */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGTGT", 9, 6,
                   "GGG", 3,
                   "CCCTGG", 6, 3);

    aseq_build_codon_rail(&seq);

    FuncContext ctx;
    init_test_context(&ctx, &seq);

    /* Junction translation should use codon rail fast path */
    FuncRuleResult r = rule_junction_translatable(&ctx);
    assert(r.passed);
    assert(ctx.junction_aa_len == 4);  /* CGPW = 4 codons */
    assert(ctx.junction_aa[0] == 'C');
    assert(ctx.junction_aa[ctx.junction_aa_len - 1] == 'W');

    aseq_reset(&seq);
}

static void test_codon_rail_mutation_updates_stop_count(void) {
    /* Mutating a base should update the stop codon count via codon rail */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGTGT", 9, 6,
                   "GGG", 3,
                   "CCCTGG", 6, 3);

    aseq_build_codon_rail(&seq);
    assert(seq.n_stop_codons == 0);

    /* Mutate first base A→T to create TAA stop codon at pos 0 */
    /* ATG → TTG (not stop), but let's mutate to create TAA */
    /* Pos 0='A', pos 1='T', pos 2='G'. Change pos 2 to 'A' → ATA=I (not stop) */
    /* Actually, to create TAA: change pos 0 to T → T,T,G = TTG (not stop).
     * Need pos 0=T, pos 1=A, pos 2=A. Let's just mutate the bases directly. */

    /* Original: ATG ATG TGT GGG CCC TGG */
    /* Mutate pos 0: A→T, pos 2: G→A → TAA at pos 0 = stop */
    Nuc *n0 = seq.head;
    Nuc *n2 = n0->next->next;
    aseq_mutate(&seq, n0, 'T', NUC_FLAG_MUTATED);
    aseq_mutate(&seq, n2, 'A', NUC_FLAG_MUTATED);

    /* Now codon at pos 0 = T,T,A = not stop (TTG was original...) */
    /* Actually: n0='T', n1='T' (original), n2='A' → TTA = Leu, not stop */
    /* Let me fix: need T,A,A. pos 0→T already, pos 1→A */
    Nuc *n1 = n0->next;
    aseq_mutate(&seq, n1, 'A', NUC_FLAG_MUTATED);
    /* Now pos 0='T', pos 1='A', pos 2='A' → TAA = stop */
    assert(seq.n_stop_codons == 1);

    /* Revert pos 1 back to T → TTA = Leu, not stop */
    aseq_mutate(&seq, n1, 'T', NUC_FLAG_MUTATED);
    assert(seq.n_stop_codons == 0);

    aseq_reset(&seq);
}

static void test_codon_rail_assess_with_rail(void) {
    /* Full assess with codon rail pre-built should give same results */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGTGT", 9, 6,
                   "GGG", 3,
                   "CCCTGG", 6, 3);

    /* Build codon rail before assess */
    aseq_build_codon_rail(&seq);

    SimRecord rec;
    sim_record_init(&rec);

    FunctionalityValidator v = functionality_validator_default();
    functionality_assess(&v, &seq, &rec);

    assert(rec.productive);
    assert(!rec.stop_codon);
    assert(rec.vj_in_frame);

    aseq_reset(&seq);
}

/* ═══════════════════════════════════════════════════════════════
 * Utilities
 * ═══════════════════════════════════════════════════════════════ */

static void test_translate_codon_utility(void) {
    /* Test the shared translate_codon function */
    assert(translate_codon('A', 'T', 'G') == 'M');
    assert(translate_codon('T', 'G', 'T') == 'C');
    assert(translate_codon('T', 'G', 'G') == 'W');
    assert(translate_codon('T', 'T', 'T') == 'F');
    assert(translate_codon('T', 'T', 'C') == 'F');
    assert(translate_codon('T', 'A', 'A') == '*');
    assert(translate_codon('T', 'A', 'G') == '*');
    assert(translate_codon('T', 'G', 'A') == '*');
    assert(translate_codon('N', 'A', 'A') == '?');
}

static void test_validator_customization(void) {
    /* Build a custom validator with only 2 rules */
    FunctionalityValidator v;
    v.n_rules = 2;
    v.rules[0] = (FuncRule){
        .name = "Stop Codon Check",
        .check = rule_stop_codon,
        .is_fatal = true,
        .contributes_to_vj_in_frame = true,
        .requires_translation = false,
    };
    v.rules[1] = (FuncRule){
        .name = "Frame Alignment Check",
        .check = rule_frame_alignment,
        .is_fatal = false,
        .contributes_to_vj_in_frame = true,
        .requires_translation = false,
    };

    /* A sequence that would fail conserved cysteine but passes the 2 rules */
    ASeq seq;
    build_test_seq(&seq,
                   "ATGATGATG", 9, 6,   /* ATG = M, not C */
                   "GGG", 3,
                   "CCCTGG", 6, 3);

    SimRecord rec;
    sim_record_init(&rec);
    functionality_assess(&v, &seq, &rec);

    /* With only stop_codon + frame rules, this should be "productive" */
    assert(rec.productive);

    aseq_reset(&seq);
}

/* ═══════════════════════════════════════════════════════════════
 * Main
 * ═══════════════════════════════════════════════════════════════ */

int main(void) {
    printf("test_functionality: Rule-based FunctionalityValidator\n");
    printf("═══════════════════════════════════════════════════════════\n\n");

    printf("Rule 1 — Stop Codon:\n");
    RUN_TEST(test_stop_codon_absent);
    RUN_TEST(test_stop_codon_in_v_region);
    RUN_TEST(test_stop_codon_in_junction);
    RUN_TEST(test_stop_codon_scans_entire_sequence);

    printf("\nRule 2 — Frame Alignment:\n");
    RUN_TEST(test_frame_aligned);
    RUN_TEST(test_frame_misaligned);

    printf("\nRule 3 — Junction Translatable:\n");
    RUN_TEST(test_junction_translatable_ok);
    RUN_TEST(test_junction_translatable_bad_length);

    printf("\nRule 4 — Conserved Cysteine:\n");
    RUN_TEST(test_conserved_cysteine_present);
    RUN_TEST(test_conserved_cysteine_absent);

    printf("\nRule 5 — Conserved Anchor (W/F):\n");
    RUN_TEST(test_conserved_anchor_w);
    RUN_TEST(test_conserved_anchor_f);
    RUN_TEST(test_conserved_anchor_absent);

    printf("\nFull Assessment:\n");
    RUN_TEST(test_assess_productive);
    RUN_TEST(test_assess_stop_codon_is_fatal);
    RUN_TEST(test_assess_out_of_frame);
    RUN_TEST(test_assess_missing_cysteine);
    RUN_TEST(test_assess_missing_anchor_wf);
    RUN_TEST(test_assess_missing_anchors);
    RUN_TEST(test_assess_vdj_productive);

    printf("\nCodon Rail Integration:\n");
    RUN_TEST(test_codon_rail_stop_codon_fast_path);
    RUN_TEST(test_codon_rail_no_stop_fast_path);
    RUN_TEST(test_codon_rail_junction_translation);
    RUN_TEST(test_codon_rail_mutation_updates_stop_count);
    RUN_TEST(test_codon_rail_assess_with_rail);

    printf("\nUtilities:\n");
    RUN_TEST(test_translate_codon_utility);
    RUN_TEST(test_validator_customization);

    printf("\n%d/%d tests passed.\n", tests_passed, tests_total);
    return tests_passed == tests_total ? 0 : 1;
}
