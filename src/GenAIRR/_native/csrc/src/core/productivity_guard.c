#include "genairr/productivity_guard.h"
#include "genairr/trace.h"
#include <stddef.h>

static const char *guard_stage_name(ProductivityStage stage) {
    switch (stage) {
        case PROD_STAGE_REARRANGEMENT: return "rearrangement";
        case PROD_STAGE_MOLECULE:      return "molecule";
        case PROD_STAGE_OBSERVED:      return "observed";
        default:                       return "unknown";
    }
}

static bool productive_guard_enabled(const SimConfig *cfg) {
    return simcfg_productive_only(cfg);
}

static void evaluate_stage_status(const ASeq *seq, const SimRecord *rec,
                                  ProductivityStage stage,
                                  ProductivityStatus *out) {
    if (!seq || !rec || !out) return;

    if (stage == PROD_STAGE_OBSERVED && rec->is_reverse_complement) {
        ASeq analysis_seq = *seq;
        aseq_rebase(&analysis_seq, seq);
        aseq_reverse_complement(&analysis_seq);
        productivity_evaluate_status(&analysis_seq, rec, out);
        return;
    }

    productivity_evaluate_status(seq, rec, out);
}

ProductivityDecision productivity_guard_state(
        const SimConfig *cfg, const ASeq *seq, const SimRecord *rec,
        ProductivityStage stage) {
    if (!cfg || !seq || !rec) {
        return PROD_DECISION_CONFIG_ERROR;
    }

    if (!productive_guard_enabled(cfg)) {
        return PROD_DECISION_ALLOW;
    }

    ProductivityStatus status = {0};
    evaluate_stage_status(seq, rec, stage, &status);
    if (!status.valid) {
        return PROD_DECISION_CONFIG_ERROR;
    }

    if (!status.productive) {
        TRACE("[guard] reject state stage=%s note=%s",
              guard_stage_name(stage),
              status.note[0] ? status.note : "non-productive");
        return PROD_DECISION_REJECT_EDIT;
    }

    return PROD_DECISION_ALLOW;
}

ProductivityDecision productivity_guard_substitution(
        const SimConfig *cfg, const ASeq *seq, const SimRecord *rec,
        ProductivityStage stage, const Nuc *node, char new_base) {
    if (!cfg || !seq || !rec || !node) {
        return PROD_DECISION_CONFIG_ERROR;
    }

    if (!productive_guard_enabled(cfg)) {
        return PROD_DECISION_ALLOW;
    }

    ptrdiff_t idx = node - seq->pool;
    if (idx < 0 || idx >= seq->pool_used) {
        return PROD_DECISION_CONFIG_ERROR;
    }

    if (node->current == new_base) {
        return PROD_DECISION_ALLOW;
    }

    ASeq seq_copy = *seq;
    aseq_rebase(&seq_copy, seq);

    SimRecord rec_copy = *rec;
    Nuc *node_copy = &seq_copy.pool[idx];

    aseq_mutate(&seq_copy, node_copy, new_base, 0);

    ProductivityDecision decision =
        productivity_guard_state(cfg, &seq_copy, &rec_copy, stage);
    if (decision != PROD_DECISION_ALLOW) {
        TRACE("[guard] reject substitution stage=%s pos=%td %c->%c note=%s",
              guard_stage_name(stage), idx, node->current, new_base,
              "candidate would violate productive contract");
        return decision;
    }

    return PROD_DECISION_ALLOW;
}

ProductivityDecision productivity_guard_insertion(
        const SimConfig *cfg, const ASeq *seq, const SimRecord *rec,
        ProductivityStage stage, const Nuc *after, char base,
        Segment segment, uint16_t flags) {
    if (!cfg || !seq || !rec || !after) {
        return PROD_DECISION_CONFIG_ERROR;
    }

    if (!productive_guard_enabled(cfg)) {
        return PROD_DECISION_ALLOW;
    }

    ptrdiff_t idx = after - seq->pool;
    if (idx < 0 || idx >= seq->pool_used) {
        return PROD_DECISION_CONFIG_ERROR;
    }

    ASeq seq_copy = *seq;
    aseq_rebase(&seq_copy, seq);

    SimRecord rec_copy = *rec;
    Nuc *after_copy = &seq_copy.pool[idx];

    if (aseq_insert_after(&seq_copy, after_copy, base, segment, flags) == NULL) {
        return PROD_DECISION_CONFIG_ERROR;
    }

    ProductivityDecision decision =
        productivity_guard_state(cfg, &seq_copy, &rec_copy, stage);
    if (decision != PROD_DECISION_ALLOW) {
        TRACE("[guard] reject insertion stage=%s pos=%td +%c note=%s",
              guard_stage_name(stage), idx, base,
              "candidate would violate productive contract");
        return decision;
    }

    return PROD_DECISION_ALLOW;
}

ProductivityDecision productivity_guard_span_rewrite(
        const SimConfig *cfg, const ASeq *seq, const SimRecord *rec,
        ProductivityStage stage, Segment seg,
        int local_start, int local_len, const char *replacement) {
    if (!cfg || !seq || !rec || !replacement) {
        return PROD_DECISION_CONFIG_ERROR;
    }

    if (!productive_guard_enabled(cfg)) {
        return PROD_DECISION_ALLOW;
    }

    if (seg < 0 || seg >= SEG_COUNT) {
        return PROD_DECISION_CONFIG_ERROR;
    }

    if (local_start < 0 || local_len <= 0) {
        return PROD_DECISION_CONFIG_ERROR;
    }

    int seg_len = aseq_segment_length(seq, seg);
    if (local_start + local_len > seg_len) {
        return PROD_DECISION_CONFIG_ERROR;
    }

    Nuc *first = seq->seg_first[seg];
    if (!first) {
        return PROD_DECISION_CONFIG_ERROR;
    }
    ptrdiff_t first_idx = first - seq->pool;
    if (first_idx < 0 || first_idx >= seq->pool_used) {
        return PROD_DECISION_CONFIG_ERROR;
    }

    ASeq seq_copy = *seq;
    aseq_rebase(&seq_copy, seq);

    SimRecord rec_copy = *rec;

    Nuc *cursor = &seq_copy.pool[first_idx];
    int pos = 0;
    while (cursor && cursor->segment == seg && pos < local_start) {
        cursor = cursor->next;
        pos++;
    }
    if (pos != local_start || !cursor) {
        return PROD_DECISION_CONFIG_ERROR;
    }

    int written = 0;
    while (cursor && cursor->segment == seg && written < local_len) {
        cursor->current  = replacement[written];
        cursor->germline = replacement[written];
        cursor = cursor->next;
        written++;
    }
    if (written != local_len) {
        return PROD_DECISION_CONFIG_ERROR;
    }

    if (seq_copy.codon_rail_valid) {
        aseq_build_codon_rail(&seq_copy);
    }

    ProductivityDecision decision =
        productivity_guard_state(cfg, &seq_copy, &rec_copy, stage);
    if (decision != PROD_DECISION_ALLOW) {
        TRACE("[guard] reject span_rewrite stage=%s seg=%d start=%d len=%d note=%s",
              guard_stage_name(stage), (int)seg, local_start, local_len,
              "candidate would violate productive contract");
        return decision;
    }

    return PROD_DECISION_ALLOW;
}

ProductivityDecision productivity_guard_deletion(
        const SimConfig *cfg, const ASeq *seq, const SimRecord *rec,
        ProductivityStage stage, const Nuc *node) {
    if (!cfg || !seq || !rec || !node) {
        return PROD_DECISION_CONFIG_ERROR;
    }

    if (!productive_guard_enabled(cfg)) {
        return PROD_DECISION_ALLOW;
    }

    ptrdiff_t idx = node - seq->pool;
    if (idx < 0 || idx >= seq->pool_used) {
        return PROD_DECISION_CONFIG_ERROR;
    }

    ASeq seq_copy = *seq;
    aseq_rebase(&seq_copy, seq);

    SimRecord rec_copy = *rec;
    Nuc *node_copy = &seq_copy.pool[idx];

    aseq_delete(&seq_copy, node_copy);

    ProductivityDecision decision =
        productivity_guard_state(cfg, &seq_copy, &rec_copy, stage);
    if (decision != PROD_DECISION_ALLOW) {
        TRACE("[guard] reject deletion stage=%s pos=%td note=%s",
              guard_stage_name(stage), idx,
              "candidate would violate productive contract");
        return decision;
    }

    return PROD_DECISION_ALLOW;
}
