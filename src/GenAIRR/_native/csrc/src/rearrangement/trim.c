/**
 * trim.c — Segment trimming step functions.
 */

#include "genairr/pipeline.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>

/* ── Internal: sample from a trim distribution ────────────────── */

static int sample_trim(RngState *rng, const TrimDist *dist, int max_allowed) {
    if (!dist || dist->probs == NULL || dist->max_trim == 0) {
        int upper = max_allowed < 10 ? max_allowed : 10;
        if (upper <= 0) return 0;
        return (int)rng_range(rng, (uint32_t)(upper + 1));
    }

    double r = rng_uniform(rng);
    double cumulative = 0.0;
    int limit = dist->max_trim < max_allowed ? dist->max_trim : max_allowed;

    for (int i = 0; i <= limit; i++) {
        cumulative += dist->probs[i];
        if (r <= cumulative) return i;
    }
    return limit;
}

/* Choose the trim distribution for an allele: prefer the per-allele
 * pointer (set by gdc_populate_sim_config from the per-(family, gene)
 * table), fall back to the legacy single-global cfg field for the
 * embedded test data path or manual SimConfigs that never went through
 * GDC populate. sample_trim handles a NULL dist with a uniform draw. */
static const TrimDist *resolve_dist(const TrimDist *allele_dist,
                                    const TrimDist *cfg_legacy) {
    if (allele_dist && allele_dist->probs) return allele_dist;
    if (cfg_legacy && cfg_legacy->probs)   return cfg_legacy;
    return NULL;
}

/* ── Step functions ───────────────────────────────────────────── */

void step_trim_v(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)seq;
    if (!rec->v_allele) return;

    /* Anchored V: cap trim at length - anchor - 1 so the conserved
     * Cys is preserved. Anchorless V (pseudogene/partial) has no
     * anchor to protect; bound only by physical length. sample_trim
     * further caps to dist->max_trim. */
    int max_v_trim;
    if (allele_has_anchor(rec->v_allele)) {
        max_v_trim = rec->v_allele->length - rec->v_allele->anchor - 1;
    } else {
        max_v_trim = rec->v_allele->length - 1;
    }
    if (max_v_trim < 0) max_v_trim = 0;

    const TrimDist *d = resolve_dist(rec->v_allele->trim_dist_3,
                                     &cfg->v_trim_3);
    rec->v_trim_5 = 0;
    rec->v_trim_3 = sample_trim(cfg->rng, d, max_v_trim);
    TRACE("[trim_v] 5'=%d, 3'=%d (max_3'=%d, allele_len=%d, anchor=%d, "
          "per_allele=%s)",
          rec->v_trim_5, rec->v_trim_3, max_v_trim,
          rec->v_allele->length, rec->v_allele->anchor,
          rec->v_allele->trim_dist_3 ? "yes" : "no");
}

void step_trim_d(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)seq;
    if (!rec->d_allele) return;

    int d_len = rec->d_allele->length;

    const TrimDist *d5 = resolve_dist(rec->d_allele->trim_dist_5,
                                      &cfg->d_trim_5);
    const TrimDist *d3 = resolve_dist(rec->d_allele->trim_dist_3,
                                      &cfg->d_trim_3);

    rec->d_trim_5 = sample_trim(cfg->rng, d5, d_len - 1);

    int remaining = d_len - rec->d_trim_5;
    rec->d_trim_3 = sample_trim(cfg->rng, d3, remaining - 1);

    if (rec->d_trim_5 + rec->d_trim_3 > d_len) {
        rec->d_trim_3 = d_len - rec->d_trim_5;
    }

    int d_remaining = d_len - rec->d_trim_5 - rec->d_trim_3;
    TRACE("[trim_d] 5'=%d, 3'=%d (d_len=%d, remaining=%d, "
          "per_allele_5=%s, per_allele_3=%s)",
          rec->d_trim_5, rec->d_trim_3, d_len, d_remaining,
          rec->d_allele->trim_dist_5 ? "yes" : "no",
          rec->d_allele->trim_dist_3 ? "yes" : "no");
}

void step_trim_j(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)seq;
    if (!rec->j_allele) return;

    /* Anchored J: cap trim at anchor - 1 so the conserved W/F codon
     * is preserved. Anchorless J (pseudogene/partial) has no anchor
     * to protect; bound only by physical length - 1 so at least one
     * base remains. This is the T0-7 fix: pre-fix, an anchorless J
     * with anchor=0xFFFF (the unsigned-cast bug) gave max_j_trim =
     * 65534, which sample_trim then bounded only by dist->max_trim
     * — silently allowing trim past the entire J allele. */
    int max_j_trim;
    if (allele_has_anchor(rec->j_allele)) {
        max_j_trim = rec->j_allele->anchor - 1;
    } else {
        max_j_trim = rec->j_allele->length - 1;
    }
    if (max_j_trim < 0) max_j_trim = 0;

    const TrimDist *d = resolve_dist(rec->j_allele->trim_dist_5,
                                     &cfg->j_trim_5);
    rec->j_trim_5 = sample_trim(cfg->rng, d, max_j_trim);
    rec->j_trim_3 = 0;
    TRACE("[trim_j] 5'=%d, 3'=%d (max_5'=%d, allele_len=%d, anchor=%d, "
          "per_allele=%s)",
          rec->j_trim_5, rec->j_trim_3, max_j_trim,
          rec->j_allele->length, rec->j_allele->anchor,
          rec->j_allele->trim_dist_5 ? "yes" : "no");
}
