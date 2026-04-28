/**
 * test_anchor.c — Unit tests for the anchor resolution subsystem.
 *
 * Covers:
 *   - codon classification (canonical, IUPAC, stop codons)
 *   - explicit-anchor strategy (Cys, Trp, IUPAC, stop, out-of-bounds)
 *   - IMGT-gapped strategy (deterministic gap counting; short-gapped
 *     rejection; gap chars at anchor position)
 *   - V motif fallback (single candidate, ambiguous, none)
 *   - J motif (locus-aware Trp vs Phe; first-match-wins; frame scan)
 *   - resolver orchestration (custom finder priority; strict mode)
 */

#include "genairr/anchor.h"
#include "genairr/ref_record.h"
#include "genairr/ref_loader.h"
#include <stdio.h>
#include <string.h>

#define TEST(name) \
    do { \
        printf("  %-50s", #name); \
        int _r = name(); \
        printf("%s\n", _r == 0 ? "PASS" : "FAIL"); \
        if (_r) failures++; \
        total++; \
    } while (0)

/* Helper: build a minimal LoadedAlleleRecord for V tests. */
static LoadedAlleleRecord make_v_record(const char *seq,
                                        const char *gapped,
                                        Locus locus,
                                        int explicit_anchor) {
    LoadedAlleleRecord r;
    loaded_allele_record_init(&r);
    r.name = "TEST*01";
    r.segment = SEG_V;
    r.locus = locus;
    r.sequence = seq;
    r.sequence_length = (int)strlen(seq);
    if (gapped) {
        r.gapped_sequence = gapped;
        r.gapped_length = (int)strlen(gapped);
        r.gap_convention_imgt = true;
    }
    r.functional_status = FUNC_F;
    r.explicit_anchor = explicit_anchor;
    r.source = "test";
    return r;
}

/* ── Explicit strategy ─────────────────────────────────────────── */

static int test_explicit_v_canonical_cys(void) {
    /* Cys (TGT) at position 0 of a 3-char sequence. */
    LoadedAlleleRecord rec = make_v_record("tgt", NULL, LOCUS_IGH, 0);
    AnchorResolverConfig cfg;
    anchor_resolver_config_init(&cfg);
    AnchorResult r = anchor_resolve_v(&cfg, &rec);
    if (r.confidence != CONF_CANONICAL) return 1;
    if (r.position != 0) return 1;
    if (r.residue != ANCHOR_CYS) return 1;
    if (r.method != METHOD_EXPLICIT) return 1;
    return 0;
}

static int test_explicit_v_canonical_tgc(void) {
    LoadedAlleleRecord rec = make_v_record("tgc", NULL, LOCUS_IGH, 0);
    AnchorResolverConfig cfg;
    anchor_resolver_config_init(&cfg);
    AnchorResult r = anchor_resolve_v(&cfg, &rec);
    if (r.confidence != CONF_CANONICAL) return 1;
    if (r.residue != ANCHOR_CYS) return 1;
    return 0;
}

static int test_explicit_v_alternative_trp(void) {
    /* Trp (TGG) is biologically valid for some V genes — should
     * resolve as ALTERNATIVE confidence. */
    LoadedAlleleRecord rec = make_v_record("tgg", NULL, LOCUS_IGH, 0);
    AnchorResolverConfig cfg;
    anchor_resolver_config_init(&cfg);
    AnchorResult r = anchor_resolve_v(&cfg, &rec);
    if (r.confidence != CONF_ALTERNATIVE) return 1;
    if (r.residue != ANCHOR_TRP) return 1;
    return 0;
}

static int test_explicit_v_trp_disallowed(void) {
    /* When allow_trp_v=false, Trp at V is rejected. */
    LoadedAlleleRecord rec = make_v_record("tgg", NULL, LOCUS_IGH, 0);
    AnchorResolverConfig cfg;
    anchor_resolver_config_init(&cfg);
    cfg.allow_trp_v = false;
    /* Disable IMGT-gapped and motif fallbacks by clearing gapped + short seq.
     * Already done — gapped is NULL and seq is only 3 chars.
     * But motif fallback runs on the 3-char seq too — let's lengthen it
     * to ensure motif also rejects, isolating the explicit strategy.
     */
    /* Use a longer sequence with no Cys anywhere to ensure motif fallback fails. */
    rec.sequence = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"; /* 66 'a's */
    rec.sequence_length = (int)strlen(rec.sequence);
    rec.explicit_anchor = -1; /* deliberately unset to test motif path */
    /* Now set explicit anchor on a Trp codon by overlaying it. */
    /* Simpler: keep tgg as the seq + set explicit_anchor=0. */
    rec.sequence = "tgg";
    rec.sequence_length = 3;
    rec.explicit_anchor = 0;

    AnchorResult r = anchor_resolve_v(&cfg, &rec);
    if (r.confidence != CONF_REJECTED) return 1;
    return 0;
}

static int test_explicit_v_iupac_alternative(void) {
    /* TGY (Y = C/T) is functionally Cys but the classifier returns
     * an explicit IUPAC path? Actually our classifier handles TGY as
     * ANCHOR_CYS — let's test a more ambiguous case: TGN (N at pos 3). */
    LoadedAlleleRecord rec = make_v_record("tgn", NULL, LOCUS_IGH, 0);
    AnchorResolverConfig cfg;
    anchor_resolver_config_init(&cfg);
    AnchorResult r = anchor_resolve_v(&cfg, &rec);
    /* "tgn" — codon_classify returns UNKNOWN (n is not in t/c/g/y).
     * codon_has_iupac returns true. So the strategy returns
     * ALTERNATIVE with residue=UNKNOWN. */
    if (r.confidence != CONF_ALTERNATIVE) return 1;
    if (r.residue != ANCHOR_RES_UNKNOWN) return 1;
    return 0;
}

static int test_explicit_v_stop_codon_rejected(void) {
    /* Stop codon at the explicit anchor — pseudogene territory. */
    LoadedAlleleRecord rec = make_v_record("tagaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                                           NULL, LOCUS_IGH, 0);
    AnchorResolverConfig cfg;
    anchor_resolver_config_init(&cfg);
    AnchorResult r = anchor_resolve_v(&cfg, &rec);
    if (r.confidence != CONF_REJECTED) return 1;
    return 0;
}

static int test_explicit_v_out_of_bounds(void) {
    /* explicit_anchor + 3 > sequence_length. */
    LoadedAlleleRecord rec = make_v_record("tgt", NULL, LOCUS_IGH, 5);
    AnchorResolverConfig cfg;
    anchor_resolver_config_init(&cfg);
    AnchorResult r = anchor_resolve_v(&cfg, &rec);
    /* Out of bounds for explicit; long seq for motif? No — seq is only 3 chars,
     * also too short for motif (requires >= 60). So all strategies reject. */
    if (r.confidence != CONF_REJECTED) return 1;
    return 0;
}

/* ── IMGT-gapped strategy ──────────────────────────────────────── */

/* Build a synthetic IMGT-style sequence: 309 chars of '.' filler
 * interspersed with real bases such that the ungapped position 0
 * corresponds to gapped position 0, and there's room for the Cys
 * codon at gapped 309-311. We construct the gapped string with
 * 309 'a's followed by the test codon followed by a few padding bases. */
static int test_imgt_gapped_canonical_cys(void) {
    char gapped[400];
    char seq[400];
    /* 309 a's, then "tgt", then "aaa". */
    for (int i = 0; i < 309; i++) gapped[i] = 'a';
    gapped[309] = 't'; gapped[310] = 'g'; gapped[311] = 't';
    gapped[312] = 'a'; gapped[313] = 'a'; gapped[314] = 'a';
    gapped[315] = '\0';
    /* Ungapped is identical (no '.') */
    memcpy(seq, gapped, 316);

    LoadedAlleleRecord rec;
    loaded_allele_record_init(&rec);
    rec.name = "TEST*01";
    rec.segment = SEG_V;
    rec.locus = LOCUS_IGH;
    rec.sequence = seq;
    rec.sequence_length = (int)strlen(seq);
    rec.gapped_sequence = gapped;
    rec.gapped_length = (int)strlen(gapped);
    rec.gap_convention_imgt = true;
    rec.functional_status = FUNC_F;
    rec.explicit_anchor = -1;
    rec.source = "test";

    AnchorResolverConfig cfg;
    anchor_resolver_config_init(&cfg);
    AnchorResult r = anchor_resolve_v(&cfg, &rec);
    if (r.confidence != CONF_CANONICAL) return 1;
    if (r.position != 309) return 1;
    if (r.residue != ANCHOR_CYS) return 1;
    if (r.method != METHOD_IMGT_GAPPED) return 1;
    return 0;
}

static int test_imgt_gapped_with_dots_counts_correctly(void) {
    /* Insert 9 gap chars early; ungapped position should be 309 - 9 = 300. */
    char gapped[400];
    int gi = 0;
    for (int i = 0; i < 50; i++) gapped[gi++] = 'a';
    for (int i = 0; i < 9;  i++) gapped[gi++] = '.';
    while (gi < 309) gapped[gi++] = 'a';
    gapped[309] = 't'; gapped[310] = 'g'; gapped[311] = 'c';
    gapped[312] = '\0';

    /* Ungapped = gapped without the dots. */
    char seq[400];
    int si = 0;
    for (int i = 0; gapped[i]; i++) {
        if (gapped[i] != '.') seq[si++] = gapped[i];
    }
    seq[si] = '\0';

    LoadedAlleleRecord rec;
    loaded_allele_record_init(&rec);
    rec.name = "TEST*01"; rec.segment = SEG_V; rec.locus = LOCUS_IGH;
    rec.sequence = seq; rec.sequence_length = si;
    rec.gapped_sequence = gapped; rec.gapped_length = (int)strlen(gapped);
    rec.gap_convention_imgt = true;
    rec.functional_status = FUNC_F; rec.explicit_anchor = -1;
    rec.source = "test";

    AnchorResolverConfig cfg;
    anchor_resolver_config_init(&cfg);
    AnchorResult r = anchor_resolve_v(&cfg, &rec);
    if (r.confidence != CONF_CANONICAL) return 1;
    if (r.position != 300) return 1;     /* 309 minus the 9 dots */
    if (r.method != METHOD_IMGT_GAPPED) return 1;
    return 0;
}

static int test_imgt_gapped_short_rejected(void) {
    /* gapped_length < 312 → strategy rejects. */
    LoadedAlleleRecord rec;
    loaded_allele_record_init(&rec);
    rec.name = "SHORT*01"; rec.segment = SEG_V; rec.locus = LOCUS_IGH;
    char short_seq[100]; for (int i = 0; i < 99; i++) short_seq[i] = 'a';
    short_seq[99] = '\0';
    rec.sequence = short_seq; rec.sequence_length = 99;
    rec.gapped_sequence = short_seq; rec.gapped_length = 99;
    rec.gap_convention_imgt = true;
    rec.functional_status = FUNC_F; rec.explicit_anchor = -1;
    rec.source = "test";

    AnchorResolverConfig cfg; anchor_resolver_config_init(&cfg);
    AnchorResult r = anchor_resolve_v(&cfg, &rec);
    /* Motif fallback also rejects (no Cys in 99 a's). Result: REJECTED. */
    if (r.confidence != CONF_REJECTED) return 1;
    return 0;
}

static int test_imgt_gapped_dot_at_anchor_rejected(void) {
    /* Anchor codon position contains dots → strategy rejects. */
    char gapped[400];
    for (int i = 0; i < 309; i++) gapped[i] = 'a';
    gapped[309] = '.'; gapped[310] = '.'; gapped[311] = '.';
    gapped[312] = 'a'; gapped[313] = '\0';
    /* Build seq with no dots */
    char seq[400];
    int si = 0;
    for (int i = 0; gapped[i]; i++) if (gapped[i] != '.') seq[si++] = gapped[i];
    seq[si] = '\0';

    LoadedAlleleRecord rec; loaded_allele_record_init(&rec);
    rec.name = "T*01"; rec.segment = SEG_V; rec.locus = LOCUS_IGH;
    rec.sequence = seq; rec.sequence_length = si;
    rec.gapped_sequence = gapped; rec.gapped_length = (int)strlen(gapped);
    rec.gap_convention_imgt = true;
    rec.functional_status = FUNC_F; rec.explicit_anchor = -1;
    rec.source = "test";

    AnchorResolverConfig cfg; anchor_resolver_config_init(&cfg);
    AnchorResult r = anchor_resolve_v(&cfg, &rec);
    if (r.confidence != CONF_REJECTED) return 1;
    return 0;
}

/* ── V motif fallback ──────────────────────────────────────────── */

static int test_motif_v_single_cys(void) {
    /* 60-char sequence with a single Cys codon in the last 30 codons,
     * frame 0. No gapped, no explicit. Strategy must pick it up as
     * BEST_GUESS. */
    /* Build: 60 chars, "tgt" at position 30 (codon 10), rest 'a'. */
    char seq[64];
    for (int i = 0; i < 60; i++) seq[i] = 'a';
    seq[30] = 't'; seq[31] = 'g'; seq[32] = 't';
    seq[60] = '\0';

    LoadedAlleleRecord rec; loaded_allele_record_init(&rec);
    rec.name = "MOTIFV*01"; rec.segment = SEG_V; rec.locus = LOCUS_IGH;
    rec.sequence = seq; rec.sequence_length = 60;
    rec.functional_status = FUNC_F; rec.explicit_anchor = -1;
    rec.source = "test";

    AnchorResolverConfig cfg; anchor_resolver_config_init(&cfg);
    AnchorResult r = anchor_resolve_v(&cfg, &rec);
    if (r.confidence != CONF_BEST_GUESS) return 1;
    if (r.position != 30) return 1;
    if (r.residue != ANCHOR_CYS) return 1;
    if (r.method != METHOD_MOTIF_SEARCH) return 1;
    return 0;
}

static int test_motif_v_multiple_rejected(void) {
    /* Two Cys codons in the search window — ambiguous. */
    char seq[100];
    for (int i = 0; i < 90; i++) seq[i] = 'a';
    seq[30] = 't'; seq[31] = 'g'; seq[32] = 't';     /* codon 10 */
    seq[60] = 't'; seq[61] = 'g'; seq[62] = 'c';     /* codon 20 */
    seq[90] = '\0';

    LoadedAlleleRecord rec; loaded_allele_record_init(&rec);
    rec.name = "AMBIG*01"; rec.segment = SEG_V; rec.locus = LOCUS_IGH;
    rec.sequence = seq; rec.sequence_length = 90;
    rec.functional_status = FUNC_F; rec.explicit_anchor = -1;
    rec.source = "test";

    AnchorResolverConfig cfg; anchor_resolver_config_init(&cfg);
    AnchorResult r = anchor_resolve_v(&cfg, &rec);
    if (r.confidence != CONF_REJECTED) return 1;
    return 0;
}

/* ── J motif (locus-aware) ─────────────────────────────────────── */

static int test_motif_j_canonical_trp_for_igh(void) {
    /* IGH J: TGG-GG?-???-GG? motif → CANONICAL (Trp expected). */
    /* Build: "tgg ggt aaa ggc" then padding. 12 nt total. */
    char seq[24];
    memcpy(seq, "tggggtaaaggc", 12);
    memcpy(seq + 12, "aaaaaaaaaaaa", 12);
    seq[24] = '\0'; seq[23] = '\0';

    LoadedAlleleRecord rec; loaded_allele_record_init(&rec);
    rec.name = "IGHJ4*01"; rec.segment = SEG_J; rec.locus = LOCUS_IGH;
    rec.sequence = seq; rec.sequence_length = 23;
    rec.functional_status = FUNC_F; rec.explicit_anchor = -1;
    rec.source = "test";

    AnchorResolverConfig cfg; anchor_resolver_config_init(&cfg);
    cfg.locus = LOCUS_IGH;
    AnchorResult r = anchor_resolve_j(&cfg, &rec);
    if (r.confidence != CONF_CANONICAL) return 1;
    if (r.residue != ANCHOR_TRP) return 1;
    if (r.position != 0) return 1;
    return 0;
}

static int test_motif_j_phe_for_igk(void) {
    /* IGK J expects Phe (TT[CT]); a Phe match should be CANONICAL. */
    char seq[24];
    memcpy(seq, "ttcggcaaagga", 12);
    seq[12] = '\0';

    LoadedAlleleRecord rec; loaded_allele_record_init(&rec);
    rec.name = "IGKJ1*01"; rec.segment = SEG_J; rec.locus = LOCUS_IGK;
    rec.sequence = seq; rec.sequence_length = 12;
    rec.functional_status = FUNC_F; rec.explicit_anchor = -1;
    rec.source = "test";

    AnchorResolverConfig cfg; anchor_resolver_config_init(&cfg);
    cfg.locus = LOCUS_IGK;
    AnchorResult r = anchor_resolve_j(&cfg, &rec);
    if (r.confidence != CONF_CANONICAL) return 1;
    if (r.residue != ANCHOR_PHE) return 1;
    return 0;
}

static int test_motif_j_canonical_phe_for_trb(void) {
    /* T2-9 real-data finding: TRB J uses Phe at IMGT position 118
     * (NOT Trp — that's IGH only, per Lefranc et al.). Verifies
     * the locus-aware confidence labelling is correct for TRB. */
    char seq[24];
    memcpy(seq, "ttcggcaaagga", 12);
    seq[12] = '\0';

    LoadedAlleleRecord rec; loaded_allele_record_init(&rec);
    rec.name = "TRBJ1-1*01"; rec.segment = SEG_J; rec.locus = LOCUS_TRB;
    rec.sequence = seq; rec.sequence_length = 12;
    rec.functional_status = FUNC_F; rec.explicit_anchor = -1;
    rec.source = "test";

    AnchorResolverConfig cfg; anchor_resolver_config_init(&cfg);
    cfg.locus = LOCUS_TRB;
    AnchorResult r = anchor_resolve_j(&cfg, &rec);
    if (r.confidence != CONF_CANONICAL) return 1;
    if (r.residue != ANCHOR_PHE) return 1;
    return 0;
}

static int test_motif_j_canonical_phe_for_trd(void) {
    /* T2-9 real-data finding: TRD J uses Phe (verified empirically
     * on the shipped HUMAN_TRD_IMGT pickle: 4/4 J alleles are Phe). */
    char seq[24];
    memcpy(seq, "ttcggcaaagga", 12);
    seq[12] = '\0';

    LoadedAlleleRecord rec; loaded_allele_record_init(&rec);
    rec.name = "TRDJ1*01"; rec.segment = SEG_J; rec.locus = LOCUS_TRD;
    rec.sequence = seq; rec.sequence_length = 12;
    rec.functional_status = FUNC_F; rec.explicit_anchor = -1;
    rec.source = "test";

    AnchorResolverConfig cfg; anchor_resolver_config_init(&cfg);
    cfg.locus = LOCUS_TRD;
    AnchorResult r = anchor_resolve_j(&cfg, &rec);
    if (r.confidence != CONF_CANONICAL) return 1;
    if (r.residue != ANCHOR_PHE) return 1;
    return 0;
}

static int test_motif_j_alternative_trp_for_trb(void) {
    /* And the inverse: a Trp anchor in TRB is non-canonical. */
    char seq[24];
    memcpy(seq, "tggggcaaaggc", 12);
    seq[12] = '\0';

    LoadedAlleleRecord rec; loaded_allele_record_init(&rec);
    rec.name = "TRBJ?*01"; rec.segment = SEG_J; rec.locus = LOCUS_TRB;
    rec.sequence = seq; rec.sequence_length = 12;
    rec.functional_status = FUNC_F; rec.explicit_anchor = -1;
    rec.source = "test";

    AnchorResolverConfig cfg; anchor_resolver_config_init(&cfg);
    cfg.locus = LOCUS_TRB;
    AnchorResult r = anchor_resolve_j(&cfg, &rec);
    if (r.confidence != CONF_ALTERNATIVE) return 1;
    if (r.residue != ANCHOR_TRP) return 1;
    return 0;
}

static int test_motif_j_alternative_phe_for_igh(void) {
    /* IGH expects Trp; a Phe match in IGH is ALTERNATIVE. */
    char seq[24];
    memcpy(seq, "ttcggcaaagga", 12);
    seq[12] = '\0';

    LoadedAlleleRecord rec; loaded_allele_record_init(&rec);
    rec.name = "IGHJ?*01"; rec.segment = SEG_J; rec.locus = LOCUS_IGH;
    rec.sequence = seq; rec.sequence_length = 12;
    rec.functional_status = FUNC_F; rec.explicit_anchor = -1;
    rec.source = "test";

    AnchorResolverConfig cfg; anchor_resolver_config_init(&cfg);
    cfg.locus = LOCUS_IGH;
    AnchorResult r = anchor_resolve_j(&cfg, &rec);
    if (r.confidence != CONF_ALTERNATIVE) return 1;
    if (r.residue != ANCHOR_PHE) return 1;
    return 0;
}

static int test_motif_j_unknown_locus_best_guess(void) {
    /* No locus → BEST_GUESS, residue still reflects what was found. */
    char seq[24];
    memcpy(seq, "tggggcaaaggt", 12);
    seq[12] = '\0';

    LoadedAlleleRecord rec; loaded_allele_record_init(&rec);
    rec.name = "?*01"; rec.segment = SEG_J; rec.locus = LOCUS_UNKNOWN;
    rec.sequence = seq; rec.sequence_length = 12;
    rec.functional_status = FUNC_F; rec.explicit_anchor = -1;
    rec.source = "test";

    AnchorResolverConfig cfg; anchor_resolver_config_init(&cfg);
    /* cfg.locus stays UNKNOWN */
    AnchorResult r = anchor_resolve_j(&cfg, &rec);
    if (r.confidence != CONF_BEST_GUESS) return 1;
    if (r.residue != ANCHOR_TRP) return 1;
    return 0;
}

/* ── Resolver orchestration ────────────────────────────────────── */

static AnchorResult my_custom_finder(const AnchorResolverConfig *cfg,
                                     const LoadedAlleleRecord *rec,
                                     void *user_data) {
    (void)cfg; (void)rec;
    int *called = (int *)user_data;
    *called = 1;
    /* Always succeed at position 42 with Trp. */
    return anchor_make_result(42, "tgg", ANCHOR_TRP,
                              CONF_ALTERNATIVE, METHOD_NONE);
    /* Resolver overrides method to METHOD_CUSTOM. */
}

static int test_custom_finder_runs_first(void) {
    LoadedAlleleRecord rec = make_v_record("tgt", NULL, LOCUS_IGH, 0);
    AnchorResolverConfig cfg;
    anchor_resolver_config_init(&cfg);
    int called = 0;
    cfg.custom_finder = my_custom_finder;
    cfg.custom_user_data = &called;

    AnchorResult r = anchor_resolve_v(&cfg, &rec);
    if (!called) return 1;
    if (r.position != 42) return 1;       /* custom's value, not explicit's 0 */
    if (r.method != METHOD_CUSTOM) return 1;  /* resolver forces method */
    return 0;
}

static int test_strict_mode_rejects_best_guess(void) {
    char seq[64];
    for (int i = 0; i < 60; i++) seq[i] = 'a';
    seq[30] = 't'; seq[31] = 'g'; seq[32] = 't';
    seq[60] = '\0';

    LoadedAlleleRecord rec; loaded_allele_record_init(&rec);
    rec.name = "STRICT*01"; rec.segment = SEG_V; rec.locus = LOCUS_IGH;
    rec.sequence = seq; rec.sequence_length = 60;
    rec.functional_status = FUNC_F; rec.explicit_anchor = -1;
    rec.source = "test";

    AnchorResolverConfig cfg; anchor_resolver_config_init(&cfg);
    cfg.strict = true;
    AnchorResult r = anchor_resolve_v(&cfg, &rec);
    /* In non-strict, this would be CONF_BEST_GUESS. Strict downgrades to REJECTED. */
    if (r.confidence != CONF_REJECTED) return 1;
    return 0;
}

/* ── Locus / segment inference ─────────────────────────────────── */

static int test_locus_from_gene_name(void) {
    if (locus_from_gene_name("IGHV1-2*01") != LOCUS_IGH) return 1;
    if (locus_from_gene_name("IGKJ1*01")   != LOCUS_IGK) return 1;
    if (locus_from_gene_name("IGLV3-1*01") != LOCUS_IGL) return 1;
    if (locus_from_gene_name("TRBV20-1*01")!= LOCUS_TRB) return 1;
    if (locus_from_gene_name("TRDD1*01")   != LOCUS_TRD) return 1;
    if (locus_from_gene_name("TRGV5*01")   != LOCUS_TRG) return 1;
    if (locus_from_gene_name("TRAJ23*01")  != LOCUS_TRA) return 1;
    if (locus_from_gene_name("UNKNOWN")    != LOCUS_UNKNOWN) return 1;
    return 0;
}

static int test_segment_from_gene_name(void) {
    if (segment_from_gene_name("IGHV1-2*01") != SEG_V) return 1;
    if (segment_from_gene_name("IGHD3-3*01") != SEG_D) return 1;
    if (segment_from_gene_name("IGHJ4*02")   != SEG_J) return 1;
    if (segment_from_gene_name("IGHC*01")    != SEG_C) return 1;
    if (segment_from_gene_name("IGHM*01")    != SEG_C) return 1;  /* C-region */
    if (segment_from_gene_name("UNK")        != SEG_UNKNOWN) return 1;
    return 0;
}

/* ── Codon classification (direct) ─────────────────────────────── */

#include "../src/anchor/codon.h"

static int test_codon_classify_basics(void) {
    if (codon_classify("tgt") != ANCHOR_CYS) return 1;
    if (codon_classify("tgc") != ANCHOR_CYS) return 1;
    if (codon_classify("tgg") != ANCHOR_TRP) return 1;
    if (codon_classify("ttt") != ANCHOR_PHE) return 1;
    if (codon_classify("ttc") != ANCHOR_PHE) return 1;
    if (codon_classify("aaa") != ANCHOR_RES_UNKNOWN) return 1;
    return 0;
}

static int test_codon_classify_iupac_y(void) {
    /* TGY ambiguity → Cys (Y is C/T, both make Cys). */
    if (codon_classify("tgy") != ANCHOR_CYS) return 1;
    /* TTY ambiguity → Phe. */
    if (codon_classify("tty") != ANCHOR_PHE) return 1;
    return 0;
}

static int test_codon_is_stop(void) {
    if (!codon_is_stop("taa")) return 1;
    if (!codon_is_stop("tag")) return 1;
    if (!codon_is_stop("tga")) return 1;
    if (!codon_is_stop("tar")) return 1;   /* IUPAC R = A/G covers TAA+TAG */
    if (codon_is_stop("tgt")) return 1;    /* not stop */
    if (codon_is_stop("tgg")) return 1;
    return 0;
}

/* ── Driver ────────────────────────────────────────────────────── */

int main(void) {
    int total = 0, failures = 0;
    printf("Anchor resolution subsystem tests\n");
    printf("─────────────────────────────────\n");

    /* Codon classifier */
    TEST(test_codon_classify_basics);
    TEST(test_codon_classify_iupac_y);
    TEST(test_codon_is_stop);

    /* Explicit V */
    TEST(test_explicit_v_canonical_cys);
    TEST(test_explicit_v_canonical_tgc);
    TEST(test_explicit_v_alternative_trp);
    TEST(test_explicit_v_trp_disallowed);
    TEST(test_explicit_v_iupac_alternative);
    TEST(test_explicit_v_stop_codon_rejected);
    TEST(test_explicit_v_out_of_bounds);

    /* IMGT-gapped V */
    TEST(test_imgt_gapped_canonical_cys);
    TEST(test_imgt_gapped_with_dots_counts_correctly);
    TEST(test_imgt_gapped_short_rejected);
    TEST(test_imgt_gapped_dot_at_anchor_rejected);

    /* V motif */
    TEST(test_motif_v_single_cys);
    TEST(test_motif_v_multiple_rejected);

    /* J motif */
    TEST(test_motif_j_canonical_trp_for_igh);
    TEST(test_motif_j_phe_for_igk);
    TEST(test_motif_j_canonical_phe_for_trb);    /* T2-9 biology fix */
    TEST(test_motif_j_canonical_phe_for_trd);    /* T2-9 biology fix */
    TEST(test_motif_j_alternative_trp_for_trb);  /* T2-9 biology fix */
    TEST(test_motif_j_alternative_phe_for_igh);
    TEST(test_motif_j_unknown_locus_best_guess);

    /* Resolver */
    TEST(test_custom_finder_runs_first);
    TEST(test_strict_mode_rejects_best_guess);

    /* Locus / segment inference */
    TEST(test_locus_from_gene_name);
    TEST(test_segment_from_gene_name);

    printf("\n%d/%d tests passed\n", total - failures, total);
    return failures == 0 ? 0 : 1;
}
