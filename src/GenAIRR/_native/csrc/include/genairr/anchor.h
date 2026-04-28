/**
 * anchor.h — biology-informed V/J anchor resolution for AIRR references.
 *
 * Replaces the legacy `_find_anchor` heuristics with a single, structured
 * resolver that produces an `AnchorResult` carrying position, codon,
 * residue identity, confidence, method, and reason. See
 * `.private/anchor_subsystem_plan_2026-04-28.md` for the full design.
 *
 * The resolver applies a strategy stack:
 *   1. custom_finder (caller-supplied hook, runs first)
 *   2. explicit anchor from the LoadedAlleleRecord
 *   3. IMGT-gapped derivation (deterministic gap counting)
 *   4. motif search on ungapped sequence (locus-aware for J)
 *
 * The first strategy that returns a non-rejected result wins.
 */

#ifndef GENAIRR_ANCHOR_H
#define GENAIRR_ANCHOR_H

#include <stdbool.h>
#include "genairr/types.h"   /* Segment enum (SEG_V, SEG_D, SEG_J, SEG_C, …) */
#include "genairr/export.h"  /* GENAIRR_EXPORT */

#ifdef __cplusplus
extern "C" {
#endif

/* Forward declaration — full type lives in ref_record.h. */
struct LoadedAlleleRecord;

/* ── Locus ─────────────────────────────────────────────────────── */

typedef enum {
    LOCUS_UNKNOWN = 0,
    LOCUS_IGH,
    LOCUS_IGK,
    LOCUS_IGL,
    LOCUS_TRA,
    LOCUS_TRB,
    LOCUS_TRD,
    LOCUS_TRG,
} Locus;

/* For the anchor subsystem we only care about V/D/J/C segments;
 * SEG_NP1 / SEG_NP2 / SEG_UMI / SEG_ADAPTER from types.h are runtime
 * node tags and never appear on a reference allele. We treat
 * SEG_COUNT as the "segment unknown / not yet classified" sentinel
 * because it's already the size marker for the enum. */
#define SEG_UNKNOWN SEG_COUNT

/* ── Anchor result components ──────────────────────────────────── */

/* Conserved residue at the anchor position. UNKNOWN covers IUPAC
 * ambiguity codes and any codon outside the canonical/alternative
 * residue set. */
typedef enum {
    ANCHOR_RES_UNKNOWN = 0,
    ANCHOR_CYS,    /* TGT/TGC — canonical V */
    ANCHOR_TRP,    /* TGG — alternative V; canonical J for IGH/TRB/TRD */
    ANCHOR_PHE,    /* TT[CT] — canonical J for IGK/IGL/TRA/TRG */
} AnchorResidue;

/* How confident we are in the resolved anchor.
 *
 *   CANONICAL    — expected codon at expected position
 *   ALTERNATIVE  — non-canonical residue but biologically valid (e.g.
 *                  Trp at V anchor, Phe at IGH J anchor)
 *   BEST_GUESS   — motif fallback succeeded with a unique candidate
 *                  but the position derivation was heuristic
 *   REJECTED     — no anchor could be assigned
 */
typedef enum {
    CONF_REJECTED = 0,
    CONF_BEST_GUESS,
    CONF_ALTERNATIVE,
    CONF_CANONICAL,
} AnchorConfidence;

/* Which strategy produced the result. */
typedef enum {
    METHOD_NONE = 0,
    METHOD_OVERRIDE,        /* caller-supplied anchor_override (highest precedence) */
    METHOD_CUSTOM,          /* user-provided custom_finder hook */
    METHOD_EXPLICIT,        /* sidecar coordinate (AIRR-C, IgBLAST .aux) */
    METHOD_IMGT_GAPPED,     /* gap-counting from gapped position 309/118 */
    METHOD_MOTIF_SEARCH,    /* ungapped motif scan */
} AnchorMethod;

/* Result of a single anchor resolution attempt.
 *
 * Value type — freely copyable, no destructor needed. The `reason`
 * field is either a pointer to a static string literal or NULL; the
 * resolver never allocates heap strings for it. */
typedef struct {
    int               position;       /* -1 if rejected; else 0-based ungapped */
    char              codon[4];       /* NUL-terminated lowercase, or "???" */
    AnchorResidue     residue;
    AnchorConfidence  confidence;
    AnchorMethod      method;
    const char       *reason;         /* NULL unless rejected */
} AnchorResult;

/* ── Resolver configuration ────────────────────────────────────── */

typedef struct AnchorResolverConfig AnchorResolverConfig;

/* Custom finder hook signature. Returning a result with confidence
 * != CONF_REJECTED short-circuits the strategy stack. The hook runs
 * BEFORE the built-in strategies so callers can fully override
 * resolution behavior for species-specific or non-IMGT references. */
typedef AnchorResult (*AnchorCustomFinder)(
    const AnchorResolverConfig *cfg,
    const struct LoadedAlleleRecord *rec,
    void *user_data);

struct AnchorResolverConfig {
    Locus locus;             /* required for J resolution; optional for V */
    bool  allow_trp_v;       /* default true — Trp at V anchor is biology */
    bool  strict;            /* default false — when true, REJECT best_guess results */
    AnchorCustomFinder custom_finder;  /* may be NULL */
    void  *custom_user_data;           /* opaque pointer passed to custom_finder */
};

/* Initialize a config to safe defaults: locus=UNKNOWN, allow_trp_v=true,
 * strict=false, no custom finder. */
GENAIRR_EXPORT void anchor_resolver_config_init(AnchorResolverConfig *cfg);

/* ── Resolution entry points ───────────────────────────────────── */

/* Resolve the V anchor for an allele record. Returns an AnchorResult
 * by value; check `confidence` to see whether resolution succeeded. */
GENAIRR_EXPORT AnchorResult anchor_resolve_v(
    const AnchorResolverConfig *cfg,
    const struct LoadedAlleleRecord *rec);

/* Resolve the J anchor (locus-aware Phe/Trp residue). */
GENAIRR_EXPORT AnchorResult anchor_resolve_j(
    const AnchorResolverConfig *cfg,
    const struct LoadedAlleleRecord *rec);

/* ── Result constructors ──────────────────────────────────────── */

/* Construct a rejected result. `reason` must be a static string
 * literal (or NULL); the resolver never copies it. */
GENAIRR_EXPORT AnchorResult anchor_make_rejected(const char *reason,
                                                 AnchorMethod method);

/* Construct a successful result. `codon` is copied into the result's
 * inline buffer (truncated to 3 chars + NUL). */
GENAIRR_EXPORT AnchorResult anchor_make_result(int position,
                                               const char *codon,
                                               AnchorResidue residue,
                                               AnchorConfidence confidence,
                                               AnchorMethod method);

/* ── Human-readable enum names (for diagnostics, logs, JSON) ──── */

GENAIRR_EXPORT const char *anchor_residue_name(AnchorResidue r);
GENAIRR_EXPORT const char *anchor_confidence_name(AnchorConfidence c);
GENAIRR_EXPORT const char *anchor_method_name(AnchorMethod m);
GENAIRR_EXPORT const char *locus_name(Locus l);
GENAIRR_EXPORT const char *segment_name(Segment s);

#ifdef __cplusplus
}
#endif

#endif /* GENAIRR_ANCHOR_H */
