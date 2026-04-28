/**
 * resolver.c — anchor resolution orchestrator (Chain of Responsibility).
 *
 * The resolver runs strategies in priority order and returns the first
 * non-rejected result. Custom finder hooks run BEFORE the built-in
 * stack so callers can override behavior for non-IMGT references or
 * species-specific edge cases (per design Q2).
 *
 * Strategy order — V:
 *   1. custom_finder (if configured)
 *   2. try_v_explicit       — sidecar coordinate from loader
 *   3. try_v_imgt_gapped    — deterministic IMGT-104 derivation
 *   4. try_v_motif          — last-30-codon scan (best-guess fallback)
 *
 * Strategy order — J:
 *   1. custom_finder (if configured)
 *   2. try_j_explicit       — sidecar coordinate
 *   3. try_j_motif          — locus-aware WGXG/FGXG scan
 *
 * (J has no IMGT-gapped equivalent because the J anchor's absolute
 * gapped position varies with framework length.)
 *
 * In strict mode, CONF_BEST_GUESS results from any strategy are
 * downgraded to REJECTED so callers can refuse heuristic anchors.
 */

#include <stddef.h>   /* NULL */
#include "genairr/anchor.h"
#include "strategies.h"

/* Apply strict-mode filter: a CONF_BEST_GUESS result is treated as
 * REJECTED when the resolver is configured for strict resolution.
 * Other confidences pass through unchanged. */
static AnchorResult apply_strict(const AnchorResolverConfig *cfg,
                                 AnchorResult r) {
    if (cfg && cfg->strict && r.confidence == CONF_BEST_GUESS) {
        return anchor_make_rejected(
            "strict mode: CONF_BEST_GUESS result downgraded",
            r.method);
    }
    return r;
}

/* Run the custom finder hook first if one is configured. Returns the
 * hook's result if non-rejected; otherwise REJECTED so the built-in
 * stack runs. The hook's result method is forced to METHOD_CUSTOM
 * (regardless of what the hook claims) so result provenance is honest. */
static AnchorResult run_custom(const AnchorResolverConfig *cfg,
                               const struct LoadedAlleleRecord *rec) {
    if (!cfg || !cfg->custom_finder) {
        return anchor_make_rejected(NULL, METHOD_CUSTOM);
    }
    AnchorResult r = cfg->custom_finder(cfg, rec, cfg->custom_user_data);
    if (r.confidence == CONF_REJECTED) return r;
    r.method = METHOD_CUSTOM;
    return r;
}

AnchorResult anchor_resolve_v(const AnchorResolverConfig *cfg,
                              const struct LoadedAlleleRecord *rec) {
    AnchorResult r;

    r = run_custom(cfg, rec);
    if (r.confidence != CONF_REJECTED) return apply_strict(cfg, r);

    r = try_v_explicit(cfg, rec);
    if (r.confidence != CONF_REJECTED) return apply_strict(cfg, r);

    r = try_v_imgt_gapped(cfg, rec);
    if (r.confidence != CONF_REJECTED) return apply_strict(cfg, r);

    r = try_v_motif(cfg, rec);
    if (r.confidence != CONF_REJECTED) return apply_strict(cfg, r);

    return anchor_make_rejected(
        "V anchor: no strategy yielded a result", METHOD_NONE);
}

AnchorResult anchor_resolve_j(const AnchorResolverConfig *cfg,
                              const struct LoadedAlleleRecord *rec) {
    AnchorResult r;

    r = run_custom(cfg, rec);
    if (r.confidence != CONF_REJECTED) return apply_strict(cfg, r);

    r = try_j_explicit(cfg, rec);
    if (r.confidence != CONF_REJECTED) return apply_strict(cfg, r);

    r = try_j_motif(cfg, rec);
    if (r.confidence != CONF_REJECTED) return apply_strict(cfg, r);

    return anchor_make_rejected(
        "J anchor: no strategy yielded a result", METHOD_NONE);
}
