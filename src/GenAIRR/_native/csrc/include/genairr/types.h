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
} NucFlags;

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
