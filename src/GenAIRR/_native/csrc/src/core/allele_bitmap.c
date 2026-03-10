/**
 * allele_bitmap.c — Per-position allele membership bitmaps.
 *
 * Replaces the O(N²·L) pairwise diff tables and trim correction maps
 * with O(N·L) bitmap construction and O(L·N/64) per-query derivation.
 *
 * GOD-ALIGNER SEMANTICS:
 *   The allele call uses node->current (the actual observable base),
 *   NOT node->germline. This ensures the call answers "which reference
 *   alleles best match the OBSERVABLE sequence?" — exactly what a
 *   perfect aligner would derive from the sequence alone.
 *
 *   Mutations, sequencing errors, and PCR errors all affect the call,
 *   just as they would affect any real alignment tool. N bases provide
 *   no evidence (neutral). Trimmed/deleted nodes don't exist in the
 *   list, so their votes vanish automatically. Inserted nodes (NP,
 *   indels) have no germline_pos and don't contribute.
 */

#include "genairr/allele_bitmap.h"
#include "genairr/pipeline.h"
#include "genairr/sim_config.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ── Portable popcount ───────────────────────────────────────── */

static inline int popcount64(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(x);
#else
    /* Fallback: Hamming weight */
    x = x - ((x >> 1) & 0x5555555555555555ULL);
    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
    return (int)((x * 0x0101010101010101ULL) >> 56);
#endif
}

/* ── Build bitmap from AllelePool ────────────────────────────── */

AlleleBitmap allele_bitmap_build(const AllelePool *pool) {
    AlleleBitmap bm;
    memset(&bm, 0, sizeof(bm));

    if (pool->count == 0) return bm;

    bm.n_alleles = pool->count;
    bm.n_words = BITMAP_WORDS(pool->count);

    /* Find maximum allele length */
    int max_len = 0;
    for (int i = 0; i < pool->count; i++) {
        if (pool->alleles[i].length > max_len)
            max_len = pool->alleles[i].length;
    }
    bm.max_pos = max_len;

    /* Allocate: max_pos × 4 bases × n_words uint64_t */
    size_t total = (size_t)max_len * 4 * bm.n_words;
    bm.data = calloc(total, sizeof(uint64_t));
    if (!bm.data) return bm;

    /* Populate: for each allele, for each position, set the bit */
    for (int a = 0; a < pool->count; a++) {
        const Allele *al = &pool->alleles[a];
        int word_idx = a / 64;
        uint64_t bit = (uint64_t)1 << (a % 64);

        for (int p = 0; p < al->length; p++) {
            int bi = bitmap_base_idx(al->seq[p]);
            if (bi < 0) continue;  /* skip N's in reference */

            /* data[(p * 4 + bi) * n_words + word_idx] |= bit */
            bm.data[(p * 4 + bi) * bm.n_words + word_idx] |= bit;
        }
    }

    /* Compute short-D threshold (for D bitmaps only):
     * Find the minimum substring length at which at least one
     * allele can be uniquely identified. We approximate this by
     * finding positions where alleles diverge. */
    if (pool->count > 1) {
        /* Simple heuristic: minimum allele length / 3, clamped to [3, 10] */
        int min_len = pool->alleles[0].length;
        for (int i = 1; i < pool->count; i++) {
            if (pool->alleles[i].length < min_len)
                min_len = pool->alleles[i].length;
        }
        bm.short_d_threshold = min_len / 3;
        if (bm.short_d_threshold < 3) bm.short_d_threshold = 3;
        if (bm.short_d_threshold > 10) bm.short_d_threshold = 10;
    } else {
        bm.short_d_threshold = 1;
    }

    return bm;
}

void allele_bitmap_destroy(AlleleBitmap *bm) {
    if (bm->data) {
        free(bm->data);
        bm->data = NULL;
    }
    bm->n_alleles = 0;
    bm->max_pos = 0;
}

/* ── Allele call derivation (carry-save / positional popcount) ─ */

/*
 * MAX_DIGIT_PLANES: number of binary digit planes needed for the
 * carry-save accumulator. With max ~500 positions per segment,
 * we need ceil(log2(500)) = 9 planes. We use 10 for safety
 * (supports counts up to 1023).
 */
#define MAX_DIGIT_PLANES 10

void allele_call_derive(const AlleleBitmap *bm, const ASeq *seq,
                         Segment seg, AlleleCallResult *out,
                         const SimRecord *rec) {
    memset(out, 0, sizeof(*out));

    if (!bm->data || bm->n_alleles == 0) return;
    if (!seq->seg_first[seg]) return;

    int n_words = bm->n_words;

    /*
     * Carry-save accumulator: digit_planes[d][w] holds the d-th
     * binary digit of the per-allele match count for alleles in
     * word w. After processing all positions, the count for allele
     * (w*64 + bit) is the binary number formed by reading bit `bit`
     * from digit_planes[0][w], digit_planes[1][w], ..., with
     * plane 0 being the LSB.
     */
    uint64_t planes[MAX_DIGIT_PLANES][BITMAP_MAX_WORDS];
    memset(planes, 0, sizeof(planes));

    int total_positions = 0;

    /* ── Helper macro: vote one (germline_pos, base) pair ────────── */
#define VOTE(germ_pos, base) do { \
        const uint64_t *_bits = allele_bitmap_lookup(bm, germ_pos, base); \
        if (_bits) { \
            total_positions++; \
            for (int _w = 0; _w < n_words; _w++) { \
                uint64_t _carry = _bits[_w]; \
                for (int _d = 0; _d < MAX_DIGIT_PLANES && _carry; _d++) { \
                    uint64_t _nv = planes[_d][_w] ^ _carry; \
                    _carry = planes[_d][_w] & _carry; \
                    planes[_d][_w] = _nv; \
                } \
            } \
        } \
    } while (0)

    /* ── Walk the segment nodes (core voting) ────────────────────── */
    for (Nuc *n = seq->seg_first[seg]; n; n = n->next) {
        if (n->segment != seg) break;
        if (n->germline == '\0') continue;

        VOTE(n->germline_pos, n->current);
    }

    /* ── Boundary extension voting ───────────────────────────────── *
     *
     * When NP bases at a segment boundary match the trimmed germline,
     * a god-aligner would attribute them to the segment. Include those
     * NP bases in the allele vote with their synthetic germline_pos.
     *
     * This replicates the same boundary extension logic used in
     * airr_derive_positions() (check_v_3prime_ambiguity, etc.).
     */
    if (rec) {
        /* ── V 3' extension into NP1 ─────────────────────────────── */
        if (seg == SEG_V && rec->v_allele && rec->v_trim_3 > 0
            && seq->seg_first[SEG_NP1])
        {
            int v_len = rec->v_allele->length;
            int trim_3 = rec->v_trim_3;
            const char *trimmed = rec->v_allele->seq + (v_len - trim_3);

            Nuc *np = seq->seg_first[SEG_NP1];
            for (int i = 0; i < trim_3 && np && np->segment == SEG_NP1;
                 i++, np = np->next) {
                if (trimmed[i] != np->current) break;
                VOTE(v_len - trim_3 + i, np->current);
            }
        }

        /* ── J 5' extension into NP (NP2 or NP1 for VJ chains) ──── */
        if (seg == SEG_J && rec->j_allele && rec->j_trim_5 > 0)
        {
            Segment np_seg = seq->seg_first[SEG_NP2] ? SEG_NP2 : SEG_NP1;
            if (seq->seg_last[np_seg]) {
                int trim_5 = rec->j_trim_5;
                const char *j_ref = rec->j_allele->seq;

                Nuc *np = seq->seg_last[np_seg];
                for (int i = trim_5 - 1;
                     i >= 0 && np && np->segment == np_seg;
                     i--, np = np->prev) {
                    if (j_ref[i] != np->current) break;
                    VOTE(i, np->current);
                }
            }
        }

        /* ── D 5' extension (NP1 end → D start) ─────────────────── */
        if (seg == SEG_D && rec->d_allele && rec->d_trim_5 > 0
            && seq->seg_last[SEG_NP1])
        {
            int trim_5 = rec->d_trim_5;
            const char *d_ref = rec->d_allele->seq;

            Nuc *np = seq->seg_last[SEG_NP1];
            for (int i = trim_5 - 1;
                 i >= 0 && np && np->segment == SEG_NP1;
                 i--, np = np->prev) {
                if (d_ref[i] != np->current) break;
                VOTE(i, np->current);
            }
        }

        /* ── D 3' extension (NP2 start → D end) ─────────────────── */
        if (seg == SEG_D && rec->d_allele && rec->d_trim_3 > 0
            && seq->seg_first[SEG_NP2])
        {
            int d_len = rec->d_allele->length;
            int trim_3 = rec->d_trim_3;
            const char *trimmed = rec->d_allele->seq + (d_len - trim_3);

            Nuc *np = seq->seg_first[SEG_NP2];
            for (int i = 0; i < trim_3 && np && np->segment == SEG_NP2;
                 i++, np = np->next) {
                if (trimmed[i] != np->current) break;
                VOTE(d_len - trim_3 + i, np->current);
            }
        }
    }

#undef VOTE

    out->total_positions = total_positions;

    /*
     * Extract per-allele counts from digit planes and find best.
     *
     * For each allele a (word w, bit b), reconstruct its count
     * from the digit planes: count = sum(bit_b(planes[d][w]) << d).
     *
     * We process one word at a time for cache friendliness.
     */
    int best = 0;
    int counts[GENAIRR_MAX_ALLELES];

    for (int w = 0; w < n_words; w++) {
        int base = w * 64;
        int limit = bm->n_alleles - base;
        if (limit > 64) limit = 64;

        for (int b = 0; b < limit; b++) {
            int count = 0;
            for (int d = 0; d < MAX_DIGIT_PLANES; d++) {
                if (planes[d][w] & ((uint64_t)1 << b))
                    count |= (1 << d);
            }
            counts[base + b] = count;
            if (count > best) best = count;
        }
    }

    out->best_score = best;

    /* Collect all alleles that tie with the best */
    if (best > 0) {
        for (int a = 0; a < bm->n_alleles && out->count < ALLELE_CALL_MAX; a++) {
            if (counts[a] == best) {
                out->indices[out->count++] = a;
            }
        }
    }
}

/* ── Format allele call as comma-separated string ────────────── */

void allele_call_format(const AlleleCallResult *result,
                         const AllelePool *pool, int true_idx,
                         char *buf, int bufsize) {
    if (result->count == 0 || bufsize <= 0) {
        buf[0] = '\0';
        return;
    }

    int pos = 0;

    /* True allele first, but ONLY if it's actually among the best matches.
     * With god-aligner semantics, mutations can shift the best match
     * away from the true allele — in that case, don't force-list it. */
    bool true_in_results = false;
    if (true_idx >= 0 && true_idx < pool->count) {
        for (int i = 0; i < result->count; i++) {
            if (result->indices[i] == true_idx) {
                true_in_results = true;
                break;
            }
        }
        if (true_in_results) {
            pos += snprintf(buf + pos, bufsize - pos, "%s",
                            pool->alleles[true_idx].name);
        }
    }

    /* Then all others that tied */
    for (int i = 0; i < result->count && pos < bufsize - 1; i++) {
        int idx = result->indices[i];
        if (idx == true_idx && true_in_results) continue;  /* already listed */
        if (idx < 0 || idx >= pool->count) continue;
        const char *sep = (pos > 0) ? "," : "";
        pos += snprintf(buf + pos, bufsize - pos, "%s%s",
                        sep, pool->alleles[idx].name);
    }
}

/* ── AlleleCorrectionSet ─────────────────────────────────────── */

AlleleCorrectionSet allele_correction_set_build(const SimConfig *cfg) {
    AlleleCorrectionSet set;
    set.v_bitmap = allele_bitmap_build(&cfg->v_alleles);
    set.d_bitmap = allele_bitmap_build(&cfg->d_alleles);
    set.j_bitmap = allele_bitmap_build(&cfg->j_alleles);
    return set;
}

void allele_correction_set_destroy(AlleleCorrectionSet *set) {
    allele_bitmap_destroy(&set->v_bitmap);
    allele_bitmap_destroy(&set->d_bitmap);
    allele_bitmap_destroy(&set->j_bitmap);
}
