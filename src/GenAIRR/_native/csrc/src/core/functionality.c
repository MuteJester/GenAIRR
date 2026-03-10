/**
 * functionality.c — Rule-based productivity/functionality validator.
 *
 * Matches the Python FunctionalityValidator exactly:
 *   1. StopCodon          — full-sequence stop codon scan (fatal, vj)
 *   2. FrameAlignment     — junction start/end/length % 3 (vj)
 *   3. JunctionTranslatable — translate junction (fatal)
 *   4. ConservedCysteine  — junction_aa[0] == 'C'
 *   5. ConservedAnchor    — junction_aa[-1] in {F, W}
 *
 * productive = all 5 rules pass.
 *
 * When the codon rail is valid, rules 1 and 3 use cached data
 * for faster evaluation. When invalid, they fall back to walking
 * the linked list directly.
 */

#include "genairr/functionality.h"
#include "genairr/pipeline.h"
#include "genairr/codon.h"
#include <string.h>
#include <stdio.h>

/* ═══════════════════════════════════════════════════════════════
 * Context initialization
 *
 * Single walk over the linked list to find anchors and compute
 * their integer positions. Everything else is derived.
 * ═══════════════════════════════════════════════════════════════ */

static void func_context_init(FuncContext *ctx, const ASeq *seq) {
    memset(ctx, 0, sizeof(*ctx));
    ctx->seq = seq;

    /* Walk once to find V and J anchors and their positions */
    int pos = 0;
    for (Nuc *n = seq->head; n; n = n->next) {
        if ((n->flags & NUC_FLAG_ANCHOR) && n->segment == SEG_V && !ctx->v_anchor) {
            ctx->v_anchor = n;
            ctx->junction_start = pos;
        }
        if ((n->flags & NUC_FLAG_ANCHOR) && n->segment == SEG_J && !ctx->j_anchor) {
            ctx->j_anchor = n;
            /* Junction includes the full W/F codon (3 bases starting at
             * J anchor), matching AIRR convention and Python behavior.
             * junction_end = j_anchor_pos + 3 (exclusive). */
            ctx->junction_end = pos + 3;
        }
        pos++;
    }

    ctx->anchors_present = (ctx->v_anchor != NULL && ctx->j_anchor != NULL);

    if (ctx->anchors_present) {
        ctx->junction_len = ctx->junction_end - ctx->junction_start;
    }
}

/* ═══════════════════════════════════════════════════════════════
 * Rule 1: Stop Codon
 *
 * Scans the ENTIRE sequence from position 0 in reading frame
 * (every 3 bases). Matches Python StopCodonRule exactly.
 *
 * When codon rail is valid, uses the cached n_stop_codons for O(1).
 *
 * Fatal: yes.  Contributes to vj_in_frame: yes.
 * ═══════════════════════════════════════════════════════════════ */

FuncRuleResult rule_stop_codon(FuncContext *ctx) {
    const ASeq *seq = ctx->seq;

    /* Fast path: use codon rail aggregate */
    if (seq->codon_rail_valid) {
        if (seq->n_stop_codons > 0) {
            ctx->stop_codon_found = true;
            return (FuncRuleResult){ false, "Stop codon present." };
        }
        return (FuncRuleResult){ true, NULL };
    }

    /* Slow path: walk linked list in triplets */
    Nuc *n = seq->head;

    while (n) {
        Nuc *n2 = n->next;
        Nuc *n3 = n2 ? n2->next : NULL;
        if (!n2 || !n3) break;

        char aa = translate_codon(n->current, n2->current, n3->current);
        if (aa == '*') {
            ctx->stop_codon_found = true;
            return (FuncRuleResult){ false, "Stop codon present." };
        }

        n = n3->next;
    }

    return (FuncRuleResult){ true, NULL };
}

/* ═══════════════════════════════════════════════════════════════
 * Rule 2: Frame Alignment
 *
 * Checks junction_start % 3, junction_end % 3, and
 * junction_length % 3. Matches Python FrameAlignmentRule.
 *
 * Fatal: no.  Contributes to vj_in_frame: yes.
 * ═══════════════════════════════════════════════════════════════ */

FuncRuleResult rule_frame_alignment(FuncContext *ctx) {
    if (!ctx->anchors_present)
        return (FuncRuleResult){ false, "Missing anchor(s)." };

    bool start_ok = (ctx->junction_start % 3 == 0);
    bool end_ok   = (ctx->junction_end   % 3 == 0);
    bool len_ok   = (ctx->junction_len   % 3 == 0);

    if (start_ok && end_ok && len_ok)
        return (FuncRuleResult){ true, NULL };

    return (FuncRuleResult){ false, "VJ out of frame." };
}

/* ═══════════════════════════════════════════════════════════════
 * Rule 3: Junction Translatable
 *
 * Verifies junction length is divisible by 3, then translates
 * the junction (V-anchor to J-anchor inclusive) and stores the
 * amino acid sequence in ctx->junction_aa for subsequent rules.
 *
 * When codon rail is valid, reads amino acids directly from
 * codon head nodes (no retranslation needed).
 *
 * Fatal: yes.  Contributes to vj_in_frame: no.
 * ═══════════════════════════════════════════════════════════════ */

FuncRuleResult rule_junction_translatable(FuncContext *ctx) {
    if (!ctx->anchors_present)
        return (FuncRuleResult){ false, "Missing anchor(s)." };

    if (ctx->junction_len % 3 != 0)
        return (FuncRuleResult){ false, "Junction length not divisible by 3." };

    const ASeq *seq = ctx->seq;

    /* Fast path: codon rail is valid — read amino acids from phase-0 nodes */
    if (seq->codon_rail_valid) {
        ctx->junction_aa_len = 0;
        Nuc *n = ctx->v_anchor;

        /* Walk in steps of 3 using codon rail */
        while (n && ctx->junction_aa_len < FUNC_MAX_JUNCTION_AA - 1) {
            if (n->frame_phase == 0 && n->amino_acid != '\0') {
                ctx->junction_aa[ctx->junction_aa_len++] = n->amino_acid;

                /* Check if we've passed the J anchor */
                Nuc *n2 = n->next;
                Nuc *n3 = n2 ? n2->next : NULL;
                if (n == ctx->j_anchor ||
                    (n2 && n2 == ctx->j_anchor) ||
                    (n3 && n3 == ctx->j_anchor))
                    break;

                n = n->codon_next;
            } else {
                n = n->next;
            }
        }
        ctx->junction_aa[ctx->junction_aa_len] = '\0';
        return (FuncRuleResult){ true, NULL };
    }

    /* Slow path: translate manually */
    ctx->junction_aa_len = 0;
    Nuc *n = ctx->v_anchor;

    while (n && ctx->junction_aa_len < FUNC_MAX_JUNCTION_AA - 1) {
        Nuc *n2 = n->next;
        Nuc *n3 = n2 ? n2->next : NULL;
        if (!n2 || !n3) break;

        ctx->junction_aa[ctx->junction_aa_len++] =
            translate_codon(n->current, n2->current, n3->current);

        /* Stop after the codon that contains J anchor */
        if (n == ctx->j_anchor || n2 == ctx->j_anchor || n3 == ctx->j_anchor)
            break;

        n = n3->next;
    }
    ctx->junction_aa[ctx->junction_aa_len] = '\0';

    return (FuncRuleResult){ true, NULL };
}

/* ═══════════════════════════════════════════════════════════════
 * Rule 4: Conserved Cysteine
 *
 * First amino acid of junction must be Cysteine (C).
 * This is the V region's 2nd conserved cysteine at CDR3 start.
 *
 * Fatal: no.  Requires translation: yes.
 * ═══════════════════════════════════════════════════════════════ */

FuncRuleResult rule_conserved_cysteine(FuncContext *ctx) {
    if (ctx->junction_aa_len == 0)
        return (FuncRuleResult){ false, "Cannot check: junction not translated." };

    if (ctx->junction_aa[0] != 'C')
        return (FuncRuleResult){ false, "V second C not present." };

    return (FuncRuleResult){ true, NULL };
}

/* ═══════════════════════════════════════════════════════════════
 * Rule 5: Conserved Anchor
 *
 * Last amino acid of junction must be Phenylalanine (F) or
 * Tryptophan (W). This is the J region anchor residue.
 *
 * Fatal: no.  Requires translation: yes.
 * ═══════════════════════════════════════════════════════════════ */

FuncRuleResult rule_conserved_anchor(FuncContext *ctx) {
    if (ctx->junction_aa_len == 0)
        return (FuncRuleResult){ false, "Cannot check: junction not translated." };

    char last_aa = ctx->junction_aa[ctx->junction_aa_len - 1];
    if (last_aa != 'F' && last_aa != 'W')
        return (FuncRuleResult){ false, "J anchor (W/F) not present." };

    return (FuncRuleResult){ true, NULL };
}

/* ═══════════════════════════════════════════════════════════════
 * Default validator
 * ═══════════════════════════════════════════════════════════════ */

FunctionalityValidator functionality_validator_default(void) {
    FunctionalityValidator v;
    v.n_rules = 5;

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
    v.rules[2] = (FuncRule){
        .name = "Junction Translation Check",
        .check = rule_junction_translatable,
        .is_fatal = true,
        .contributes_to_vj_in_frame = false,
        .requires_translation = false,
    };
    v.rules[3] = (FuncRule){
        .name = "Conserved Cysteine (C) Check",
        .check = rule_conserved_cysteine,
        .is_fatal = false,
        .contributes_to_vj_in_frame = false,
        .requires_translation = true,
    };
    v.rules[4] = (FuncRule){
        .name = "Conserved J Anchor (F/W) Check",
        .check = rule_conserved_anchor,
        .is_fatal = false,
        .contributes_to_vj_in_frame = false,
        .requires_translation = true,
    };

    return v;
}

/* ═══════════════════════════════════════════════════════════════
 * Assess functionality
 *
 * Runs all rules in order, respecting metadata:
 *   - requires_translation: skip if junction_aa_len == 0
 *   - is_fatal: stop on failure
 *   - contributes_to_vj_in_frame: failure sets vj_in_frame = false
 *
 * Fills: rec->productive, rec->stop_codon, rec->vj_in_frame,
 *        rec->note.
 * ═══════════════════════════════════════════════════════════════ */

void functionality_assess(const FunctionalityValidator *v,
                          ASeq *seq, SimRecord *rec) {
    FuncContext ctx;
    func_context_init(&ctx, seq);

    /* If anchors missing, fail immediately */
    if (!ctx.anchors_present) {
        rec->productive  = false;
        rec->stop_codon  = false;
        rec->vj_in_frame = false;
        seq->junction_in_frame = false;
        snprintf(rec->note, sizeof(rec->note), "Missing anchor(s)");
        return;
    }

    bool all_passed       = true;
    bool vj_in_frame      = true;
    const char *first_note = NULL;

    for (int i = 0; i < v->n_rules; i++) {
        const FuncRule *rule = &v->rules[i];

        /* Skip translation-dependent rules if no translation yet */
        if (rule->requires_translation && ctx.junction_aa_len == 0)
            continue;

        FuncRuleResult result = rule->check(&ctx);

        if (!result.passed) {
            all_passed = false;

            if (!first_note)
                first_note = result.note;

            if (rule->contributes_to_vj_in_frame)
                vj_in_frame = false;

            if (rule->is_fatal)
                break;
        }
    }

    rec->productive  = all_passed;
    rec->stop_codon  = ctx.stop_codon_found;
    rec->vj_in_frame = vj_in_frame;

    /* Set ASeq-level reactive productivity state */
    seq->junction_in_frame = vj_in_frame;

    if (first_note)
        snprintf(rec->note, sizeof(rec->note), "%s", first_note);
    else
        rec->note[0] = '\0';
}
