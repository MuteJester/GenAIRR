/**
 * types.h — Fundamental type definitions for GenAIRR-C.
 *
 * Every header in the project includes this file. It defines the
 * enums, flags, and small value types used throughout.
 */

#ifndef GENAIRR_TYPES_H
#define GENAIRR_TYPES_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

/* ── Segment identity ─────────────────────────────────────────── */

typedef enum {
    SEG_V     = 0,
    SEG_NP1   = 1,
    SEG_D     = 2,
    SEG_NP2   = 3,
    SEG_J     = 4,
    SEG_C     = 5,
    SEG_UMI   = 6,
    SEG_ADAPTER = 7,
    SEG_COUNT           /* sentinel: number of segment types */
} Segment;

/* ── Per-nucleotide annotation flags ──────────────────────────── */

typedef enum {
    NUC_FLAG_NONE         = 0,
    NUC_FLAG_MUTATED      = 1 << 0,   /* somatic hypermutation          */
    NUC_FLAG_SEQ_ERROR    = 1 << 1,   /* Illumina sequencing error      */
    NUC_FLAG_PCR_ERROR    = 1 << 2,   /* polymerase error               */
    NUC_FLAG_P_NUCLEOTIDE = 1 << 3,   /* palindromic nucleotide         */
    NUC_FLAG_N_NUCLEOTIDE = 1 << 4,   /* TdT-generated N-nucleotide    */
    NUC_FLAG_IS_N         = 1 << 5,   /* corrupted to ambiguous N       */
    NUC_FLAG_ANCHOR       = 1 << 6,   /* V-cys or J-trp/gly anchor     */
    NUC_FLAG_INDEL_INS    = 1 << 7,   /* inserted by indel operation    */
    NUC_FLAG_JUNCTION     = 1 << 8,   /* part of the junction region    */
    NUC_FLAG_INVERTED     = 1 << 9,   /* T2-13: D-segment inversion     */
} NucFlags;

/* ── Productivity filter mode ─────────────────────────────────── */

/**
 * Filter mode for which V(D)J rearrangements appear in the output.
 *
 *   PRODUCTIVITY_MIXED               (= 0, default after memset)
 *       No filtering. Pure V(D)J recombination output —
 *       biologically faithful. ~22.8% productive on human IGH
 *       (T1-15: 33% in-frame × ~68% stop-free).
 *   PRODUCTIVITY_PRODUCTIVE_ONLY     (= 1)
 *       Retry rearrangement (up to max_productive_attempts) until
 *       rec.productive == true. S5F mutation, uniform mutation,
 *       D inversion, receptor revision, selection pressure, PCR
 *       amplification, quality-profile sequencing errors, paired-end
 *       read modeling, long-read homopolymer indels, 5' end
 *       corruption, 3' end corruption, observed indels, contaminants,
 *       and insert_ns now honor that contract downstream, but other
 *       later edit paths can still flip the final observed record
 *       until the guard rollout is wider.
 *   PRODUCTIVITY_NON_PRODUCTIVE_ONLY (= 2)
 *       Symmetric — retry until rec.productive == false. Useful
 *       for negative training data and out-of-frame statistics.
 *
 * Mirrors GenAIRR.Productivity in Python; integer values stay
 * stable across the ABI for safe ctypes passthrough.
 */
typedef enum {
    PRODUCTIVITY_MIXED               = 0,
    PRODUCTIVITY_PRODUCTIVE_ONLY     = 1,
    PRODUCTIVITY_NON_PRODUCTIVE_ONLY = 2,
} ProductivityMode;

/* ── Productivity contract (V5 native-first model) ───────────── */

typedef enum {
    PROD_TARGET_MIXED               = 0,
    PROD_TARGET_PRODUCTIVE_ONLY     = 1,
    PROD_TARGET_NON_PRODUCTIVE_ONLY = 2,
} ProductivityTarget;

typedef enum {
    PROD_SCOPE_REARRANGEMENT = 0,
    PROD_SCOPE_MOLECULE      = 1,
    PROD_SCOPE_OBSERVED      = 2,
} ProductivityScope;

typedef enum {
    PROD_STRICTNESS_BEST_EFFORT = 0,
    PROD_STRICTNESS_STRICT      = 1,
} ProductivityStrictness;

typedef struct {
    ProductivityTarget      target;
    ProductivityScope       scope;
    ProductivityStrictness  strictness;
} ProductivityContract;

/* Bridge helper for the public/API boundary, which still accepts the
 * legacy ProductivityMode and converts it into the native contract. */
static inline ProductivityContract
productivity_contract_from_mode(ProductivityMode mode) {
    ProductivityContract c = {
        .target = PROD_TARGET_MIXED,
        .scope = PROD_SCOPE_REARRANGEMENT,
        .strictness = PROD_STRICTNESS_BEST_EFFORT,
    };
    if (mode == PRODUCTIVITY_PRODUCTIVE_ONLY) {
        c.target = PROD_TARGET_PRODUCTIVE_ONLY;
    } else if (mode == PRODUCTIVITY_NON_PRODUCTIVE_ONLY) {
        c.target = PROD_TARGET_NON_PRODUCTIVE_ONLY;
    }
    return c;
}

/* ── Chain type ───────────────────────────────────────────────── */

typedef enum {
    CHAIN_IGH,
    CHAIN_IGK,
    CHAIN_IGL,
    CHAIN_TCRA,
    CHAIN_TCRB,
    CHAIN_TCRD,
    CHAIN_TCRG,
    CHAIN_COUNT
} ChainType;

static inline bool chain_has_d(ChainType ct) {
    return ct == CHAIN_IGH || ct == CHAIN_TCRB || ct == CHAIN_TCRD;
}

/* ── Nucleotide helpers ───────────────────────────────────────── */

static inline char complement(char base) {
    switch (base) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'a': return 't';
        case 't': return 'a';
        case 'c': return 'g';
        case 'g': return 'c';
        default:  return 'N';
    }
}

/* ── Pipeline hook points ────────────────────────────────────── */

typedef enum {
    HOOK_POST_ASSEMBLY = 0,
    HOOK_POST_FUNCTIONALITY,
    HOOK_POST_D_INVERSION,
    HOOK_POST_RECEPTOR_REV,
    HOOK_POST_MUTATION,
    HOOK_POST_SELECTION,
    HOOK_POST_CORRUPT_5,
    HOOK_POST_CORRUPT_3,
    HOOK_POST_INDELS,
    HOOK_POST_NS,
    HOOK_POST_PCR,
    HOOK_POST_QUALITY,
    HOOK_FINAL,
    HOOK_COUNT
} HookPoint;

/* ── Common limits ────────────────────────────────────────────── */

#define GENAIRR_MAX_SEQ_LEN      1024
#define GENAIRR_MAX_ALLELE_NAME  128
#define GENAIRR_MAX_ALLELE_SEQ   512
#define GENAIRR_MAX_PIPELINE     64
#define GENAIRR_MAX_ALLELES      512

#endif /* GENAIRR_TYPES_H */
