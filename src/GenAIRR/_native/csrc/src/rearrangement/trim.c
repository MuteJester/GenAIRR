/**
 * trim.c — Segment trimming step functions.
 */

#include "genairr/pipeline.h"
#include "genairr/trace.h"
#include <stdlib.h>

/* ── Internal: sample from a trim distribution ────────────────── */

static int sample_trim(const TrimDist *dist, int max_allowed) {
    if (dist->probs == NULL || dist->max_trim == 0) {
        int upper = max_allowed < 10 ? max_allowed : 10;
        if (upper <= 0) return 0;
        return rand() % (upper + 1);
    }

    double r = (double)rand() / RAND_MAX;
    double cumulative = 0.0;
    int limit = dist->max_trim < max_allowed ? dist->max_trim : max_allowed;

    for (int i = 0; i <= limit; i++) {
        cumulative += dist->probs[i];
        if (r <= cumulative) return i;
    }
    return limit;
}

/* ── Step functions ───────────────────────────────────────────── */

void step_trim_v(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)seq;
    if (!rec->v_allele) return;

    int max_v_trim = rec->v_allele->length - rec->v_allele->anchor - 1;
    if (max_v_trim < 0) max_v_trim = 0;

    rec->v_trim_5 = 0;
    rec->v_trim_3 = sample_trim(&cfg->v_trim_3, max_v_trim);
    TRACE("[trim_v] 5'=%d, 3'=%d (max_3'=%d, allele_len=%d, anchor=%d)",
          rec->v_trim_5, rec->v_trim_3, max_v_trim,
          rec->v_allele->length, rec->v_allele->anchor);
}

void step_trim_d(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)seq;
    if (!rec->d_allele) return;

    int d_len = rec->d_allele->length;

    rec->d_trim_5 = sample_trim(&cfg->d_trim_5, d_len - 1);

    int remaining = d_len - rec->d_trim_5;
    rec->d_trim_3 = sample_trim(&cfg->d_trim_3, remaining - 1);

    if (rec->d_trim_5 + rec->d_trim_3 > d_len) {
        rec->d_trim_3 = d_len - rec->d_trim_5;
    }

    int d_remaining = d_len - rec->d_trim_5 - rec->d_trim_3;
    TRACE("[trim_d] 5'=%d, 3'=%d (d_len=%d, remaining=%d)",
          rec->d_trim_5, rec->d_trim_3, d_len, d_remaining);
}

void step_trim_j(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)seq;
    if (!rec->j_allele) return;

    int max_j_trim = rec->j_allele->anchor - 1;
    if (max_j_trim < 0) max_j_trim = 0;

    rec->j_trim_5 = sample_trim(&cfg->j_trim_5, max_j_trim);
    rec->j_trim_3 = 0;
    TRACE("[trim_j] 5'=%d, 3'=%d (max_5'=%d, allele_len=%d, anchor=%d)",
          rec->j_trim_5, rec->j_trim_3, max_j_trim,
          rec->j_allele->length, rec->j_allele->anchor);
}
