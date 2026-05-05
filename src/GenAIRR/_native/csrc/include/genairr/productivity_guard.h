#ifndef GENAIRR_PRODUCTIVITY_GUARD_H
#define GENAIRR_PRODUCTIVITY_GUARD_H

#include "pipeline.h"

typedef enum {
    PROD_DECISION_ALLOW = 0,
    PROD_DECISION_REJECT_EDIT = 1,
    PROD_DECISION_RETRY_LOCAL = 2,
    PROD_DECISION_RESAMPLE_RECORD = 3,
    PROD_DECISION_CONFIG_ERROR = 4,
} ProductivityDecision;

ProductivityDecision productivity_guard_substitution(
    const SimConfig *cfg, const ASeq *seq, const SimRecord *rec,
    ProductivityStage stage, const Nuc *node, char new_base);

ProductivityDecision productivity_guard_insertion(
    const SimConfig *cfg, const ASeq *seq, const SimRecord *rec,
    ProductivityStage stage, const Nuc *after, char base,
    Segment segment, uint16_t flags);

ProductivityDecision productivity_guard_deletion(
    const SimConfig *cfg, const ASeq *seq, const SimRecord *rec,
    ProductivityStage stage, const Nuc *node);

ProductivityDecision productivity_guard_state(
    const SimConfig *cfg, const ASeq *seq, const SimRecord *rec,
    ProductivityStage stage);

/**
 * Preflight a bulk segment rewrite that replaces `local_len` bases
 * starting at `local_start` (offset in nodes from the segment's first
 * node) within segment `seg`. The simulated rewrite sets both
 * `current` and `germline` to `replacement[i]` for each node in range
 * and rebuilds the codon rail before evaluating the contract — this
 * matches a template-level change (D inversion, V replacement) where
 * the post-rewrite bases become the new germline reference.
 *
 * Productivity-only metadata (germline_pos remap, per-node provenance
 * flags) is intentionally not modeled here because it does not affect
 * codon translation. Callers still apply the real rewrite themselves
 * if the decision is ALLOW.
 */
ProductivityDecision productivity_guard_span_rewrite(
    const SimConfig *cfg, const ASeq *seq, const SimRecord *rec,
    ProductivityStage stage, Segment seg,
    int local_start, int local_len, const char *replacement);

#endif /* GENAIRR_PRODUCTIVITY_GUARD_H */
