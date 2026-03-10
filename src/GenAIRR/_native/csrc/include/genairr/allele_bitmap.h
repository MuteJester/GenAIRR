/**
 * allele_bitmap.h — Per-position allele membership bitmaps.
 *
 * Replaces the O(N²·L) pairwise diff tables and trim correction maps
 * with a reactive O(L) allele call derivation.
 *
 * CORE IDEA: For each germline position and each possible base, a
 * bitmask records which alleles have that base at that position.
 * The allele call is derived by walking the linked list, looking up
 * each node's (germline_pos, current_base) in the bitmap, and
 * counting hits per allele via popcount.
 *
 * The allele(s) with the most matching positions = the call.
 *
 * REACTIVITY: Because the call is derived from node state, any
 * mutation, deletion, or insertion automatically changes the call.
 * No bookkeeping, no stale state.
 *
 * STORAGE: For N alleles, we need ceil(N/64) uint64_t words per
 * position per base. With GENAIRR_MAX_ALLELES=512, that's 8 words.
 * At ~300 positions × 4 bases = 9600 words = ~75KB per segment type.
 *
 * SHORT-D: If the D segment has too few matching positions to
 * distinguish any allele uniquely, the call is "Short-D".
 */

#ifndef GENAIRR_ALLELE_BITMAP_H
#define GENAIRR_ALLELE_BITMAP_H

#include "types.h"
#include "allele.h"
#include "aseq.h"
#include "sim_config.h"

/* Forward declaration — full definition in pipeline.h */
struct SimRecord;

/* ── Bitmap dimensions ───────────────────────────────────────── */

/** Number of uint64_t words needed for N alleles. */
#define BITMAP_WORDS(n)  (((n) + 63) / 64)

/** Maximum words per bitmap entry (supports up to 512 alleles). */
#define BITMAP_MAX_WORDS  BITMAP_WORDS(GENAIRR_MAX_ALLELES)  /* 8 */

/** Base index: A=0, C=1, G=2, T=3 */
static inline int bitmap_base_idx(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return -1;  /* N or unknown: no match */
    }
}

/* ── Bitmap for one segment type ─────────────────────────────── */

/**
 * AlleleBitmap: per-position, per-base membership bitmaps for one
 * segment type (V, D, or J).
 *
 * Access: bitmap[position][base] → array of `n_words` uint64_t.
 *
 * If bit `a` is set in bitmap[pos][base], then allele `a` has
 * base `base` at germline position `pos`.
 *
 * The allele-to-index mapping is the same as AllelePool order.
 */
typedef struct {
    uint64_t   *data;          /* flat array: [max_pos][4][n_words] */
    int         max_pos;       /* maximum germline position stored  */
    int         n_alleles;     /* number of alleles in the pool     */
    int         n_words;       /* ceil(n_alleles / 64)              */
    int         short_d_threshold;  /* D-only: min length for valid call */
} AlleleBitmap;

/**
 * Build a bitmap from an allele pool.
 *
 * Scans every allele's ungapped sequence and sets the appropriate
 * bits. O(N × L) construction — replaces O(N² × L) pairwise diffs.
 */
AlleleBitmap  allele_bitmap_build(const AllelePool *pool);

/**
 * Free bitmap data.
 */
void  allele_bitmap_destroy(AlleleBitmap *bm);

/* ── Bitmap lookup ───────────────────────────────────────────── */

/**
 * Get the membership bitmask for a given (position, base) pair.
 *
 * Returns a pointer to n_words uint64_t values. Each set bit
 * indicates an allele that has `base` at `germline_pos`.
 *
 * Returns NULL if position or base is out of range.
 */
static inline const uint64_t *allele_bitmap_lookup(
        const AlleleBitmap *bm, int germline_pos, char base) {
    int bi = bitmap_base_idx(base);
    if (bi < 0 || germline_pos < 0 || germline_pos >= bm->max_pos)
        return NULL;
    /* Layout: data[(pos * 4 + base_idx) * n_words ... +n_words] */
    return &bm->data[(germline_pos * 4 + bi) * bm->n_words];
}

/* ── Allele call derivation ──────────────────────────────────── */

/** Maximum number of tied alleles to report. */
#define ALLELE_CALL_MAX  64

/**
 * Result of deriving an allele call from the bitmap.
 */
typedef struct {
    int   indices[ALLELE_CALL_MAX];  /* allele pool indices (best matches) */
    int   count;                      /* number of best-matching alleles   */
    int   best_score;                 /* number of matching positions      */
    int   total_positions;            /* total positions evaluated         */
} AlleleCallResult;

/**
 * Derive the allele call for a segment by walking the ASeq.
 *
 * GOD-ALIGNER SEMANTICS: Uses node->current (the actual base in the
 * sequence) for the bitmap lookup. The question answered is:
 * "given what we can OBSERVE in the sequence, which germline alleles
 * best match at each retained position?"
 *
 * This means:
 *   - Mutations CHANGE the call. If SHM mutates a distinguishing
 *     position to match another allele, the call reflects that.
 *   - Sequencing/PCR errors also affect the call (as they would
 *     affect any real alignment tool).
 *   - N bases match nothing (bitmap_base_idx returns -1 for 'N'),
 *     providing neutral evidence — exactly what an aligner would do.
 *   - Trimmed/deleted nodes don't exist → votes vanish automatically.
 *
 * Only nodes with a valid germline origin (germline != '\0') have
 * a germline_pos for the bitmap lookup. NP nodes, indel-inserted
 * nodes, and adapter nodes are skipped — UNLESS boundary extension
 * claims them for the segment (when `rec` is provided).
 *
 * O(L × N/64) where L = segment length, N = number of alleles.
 *
 * @param bm    The bitmap for this segment type.
 * @param seq   The annotated sequence.
 * @param seg   Which segment to evaluate (SEG_V, SEG_D, SEG_J).
 * @param out   Output: best-matching allele indices and scores.
 * @param rec   SimRecord (for boundary extension). NULL to skip.
 */
void  allele_call_derive(const AlleleBitmap *bm, const ASeq *seq,
                          Segment seg, AlleleCallResult *out,
                          const struct SimRecord *rec);

/**
 * Format the allele call result as a comma-separated string.
 *
 * Uses the AllelePool to look up names. The true (simulated) allele
 * is listed first, followed by any ties.
 *
 * @param result    The derived call result.
 * @param pool      The allele pool (for name lookup).
 * @param true_idx  Index of the true (simulated) allele (-1 if unknown).
 * @param buf       Output buffer.
 * @param bufsize   Size of output buffer.
 */
void  allele_call_format(const AlleleCallResult *result,
                          const AllelePool *pool, int true_idx,
                          char *buf, int bufsize);

/**
 * Check if a D segment call should be "Short-D".
 *
 * Returns true if too few positions remain to distinguish alleles.
 */
static inline bool allele_call_is_short_d(const AlleleCallResult *result,
                                           const AlleleBitmap *bm) {
    return result->total_positions < bm->short_d_threshold;
}

/* ── Full correction context (replaces CorrectionContext) ────── */

/**
 * AlleleCorrectionSet: bitmaps for all segment types, built once
 * from SimConfig. Replaces the old CorrectionContext with its
 * O(N²·L) pairwise diff tables.
 *
 * Construction: O(N·L) per segment type.
 * Per-sequence query: O(L·N/64) per segment type.
 */
typedef struct {
    AlleleBitmap  v_bitmap;
    AlleleBitmap  d_bitmap;
    AlleleBitmap  j_bitmap;
} AlleleCorrectionSet;

/**
 * Build correction bitmaps from SimConfig allele pools.
 */
AlleleCorrectionSet  allele_correction_set_build(const SimConfig *cfg);

/**
 * Free all bitmap data.
 */
void  allele_correction_set_destroy(AlleleCorrectionSet *set);

#endif /* GENAIRR_ALLELE_BITMAP_H */
