/**
 * functionality.h — Rule-based productivity/functionality validator.
 *
 * Matches the Python FunctionalityValidator exactly: 5 sequential rules,
 * each with is_fatal and contributes_to_vj_in_frame metadata.
 *
 * Rules (in order):
 *   1. StopCodon          — scans ENTIRE sequence for stop codons (fatal, vj)
 *   2. FrameAlignment     — junction boundaries codon-aligned (vj)
 *   3. JunctionTranslatable — junction length % 3 == 0, translates (fatal)
 *   4. ConservedCysteine  — first junction AA must be C
 *   5. ConservedAnchor    — last junction AA must be F or W
 *
 * productive = all 5 rules pass.
 */

#ifndef GENAIRR_FUNCTIONALITY_H
#define GENAIRR_FUNCTIONALITY_H

#include "aseq.h"
#include "codon.h"
#include <stdbool.h>

/* ── Validation context (shared across rules) ────────────────── */

#define FUNC_MAX_JUNCTION_AA  128

typedef struct FuncContext {
    /* Sequence (read-only) */
    const ASeq *seq;

    /* Anchors (found during context init) */
    Nuc  *v_anchor;
    Nuc  *j_anchor;
    bool  anchors_present;

    /* Junction positions (0-based in assembled sequence) */
    int   junction_start;      /* position of V anchor */
    int   junction_end;        /* position after J anchor (exclusive) */
    int   junction_len;        /* nucleotide count (end - start) */

    /* Translated junction (set by JunctionTranslatable rule) */
    char  junction_aa[FUNC_MAX_JUNCTION_AA];
    int   junction_aa_len;

    /* Stop codon flag (set by StopCodon rule) */
    bool  stop_codon_found;
} FuncContext;

/* ── Rule result ─────────────────────────────────────────────── */

typedef struct {
    bool        passed;
    const char *note;          /* static string, NULL if passed */
} FuncRuleResult;

/* ── Rule function pointer ───────────────────────────────────── */

typedef FuncRuleResult (*FuncRuleFn)(FuncContext *ctx);

/* ── Rule descriptor ─────────────────────────────────────────── */

typedef struct {
    const char *name;
    FuncRuleFn  check;
    bool        is_fatal;                  /* stop validation on failure */
    bool        contributes_to_vj_in_frame;
    bool        requires_translation;      /* skip if junction_aa empty */
} FuncRule;

/* ── Validator (array of rules) ──────────────────────────────── */

#define FUNC_MAX_RULES  16

typedef struct {
    FuncRule  rules[FUNC_MAX_RULES];
    int       n_rules;
} FunctionalityValidator;

/**
 * Build the default validator with all 5 rules.
 */
FunctionalityValidator  functionality_validator_default(void);

/**
 * Run all rules against the sequence, fill SimRecord productivity fields.
 *
 * Sets: rec->productive, rec->stop_codon, rec->vj_in_frame, rec->note.
 *
 * @param v    The validator (rule set).
 * @param seq  The annotated sequence (also sets seq->junction_in_frame).
 * @param rec  Output: productivity fields filled.
 */
struct SimRecord;  /* forward declaration */
void  functionality_assess(const FunctionalityValidator *v,
                           ASeq *seq, struct SimRecord *rec);

/* ── Individual rules (exposed for testing and customization) ── */

FuncRuleResult  rule_stop_codon(FuncContext *ctx);
FuncRuleResult  rule_frame_alignment(FuncContext *ctx);
FuncRuleResult  rule_junction_translatable(FuncContext *ctx);
FuncRuleResult  rule_conserved_cysteine(FuncContext *ctx);
FuncRuleResult  rule_conserved_anchor(FuncContext *ctx);

#endif /* GENAIRR_FUNCTIONALITY_H */
