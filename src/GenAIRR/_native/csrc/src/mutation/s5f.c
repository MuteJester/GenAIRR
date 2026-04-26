/**
 * s5f.c — S5F mutation model implementation.
 *
 * Core algorithm (mirrors Python S5F._apply_standard_flat):
 *
 *   1. Extract flat base array from ASeq, padded with NN on each end.
 *   2. Compute integer 5-mer keys for each position.
 *   3. Build mutable mask from segment tags (V/D/J only).
 *   4. Set weights[i] = mutability[kmer_key[i]] for mutable positions.
 *   5. Build Fenwick tree over weights.
 *   6. Loop:
 *      a. Weighted-select a position via Fenwick tree.
 *      b. Look up substitution table for the position's 5-mer key.
 *      c. Choose a target base from the cumulative distribution.
 *      d. Apply mutation: update base array, recompute affected 5-mer keys.
 *      e. Update Fenwick tree for the ≤5 affected positions.
 *   7. Write mutations back to ASeq nodes.
 */

#include "genairr/s5f.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ═══════════════════════════════════════════════════════════════
 * Fenwick Tree
 * ═══════════════════════════════════════════════════════════════ */

void fenwick_init(FenwickTree *ft, const double *weights, int n) {
    ft->n = n;
    ft->total = 0.0;
    ft->update_count = 0;
    ft->tree = calloc(n + 1, sizeof(double));

    for (int i = 0; i < n; i++) {
        ft->tree[i + 1] = weights[i];
        ft->total += weights[i];
    }
    /* O(N) bottom-up propagation */
    for (int i = 1; i <= n; i++) {
        int j = i + (i & -i);
        if (j <= n) {
            ft->tree[j] += ft->tree[i];
        }
    }
}

void fenwick_destroy(FenwickTree *ft) {
    free(ft->tree);
    ft->tree = NULL;
    ft->n = 0;
    ft->total = 0.0;
}

void fenwick_update(FenwickTree *ft, int i, double new_w, double old_w) {
    double delta = new_w - old_w;
    if (delta == 0.0) return;
    ft->total += delta;
    ft->update_count++;
    i++;  /* 1-indexed */
    int n = ft->n;
    double *tree = ft->tree;
    while (i <= n) {
        tree[i] += delta;
        i += i & (-i);
    }
}

void fenwick_rebuild(FenwickTree *ft, const double *weights) {
    int n = ft->n;
    double *tree = ft->tree;
    double total = 0.0;
    for (int i = 0; i < n; i++) {
        tree[i + 1] = weights[i];
        total += weights[i];
    }
    for (int i = 1; i <= n; i++) {
        int j = i + (i & -i);
        if (j <= n) {
            tree[j] += tree[i];
        }
    }
    ft->total = total;
    ft->update_count = 0;
}

int fenwick_select(const FenwickTree *ft, double r) {
    int pos = 0;
    int n = ft->n;
    const double *tree = ft->tree;
    /* Find highest power of 2 <= n */
    int bit = 1;
    while (bit <= n) bit <<= 1;
    bit >>= 1;
    while (bit > 0) {
        int nxt = pos + bit;
        if (nxt <= n && tree[nxt] <= r) {
            r -= tree[nxt];
            pos = nxt;
        }
        bit >>= 1;
    }
    return pos;  /* 0-indexed */
}

/* ═══════════════════════════════════════════════════════════════
 * S5F Model init and data loading
 * ═══════════════════════════════════════════════════════════════ */

void s5f_model_init(S5FModel *model,
                     double min_rate, double max_rate,
                     bool productive) {
    memset(model, 0, sizeof(*model));
    model->min_mutation_rate = min_rate;
    model->max_mutation_rate = max_rate;
    model->productive = productive;
}

void s5f_set_substitution(S5FModel *model, int key,
                           const char *bases, const double *weights, int count) {
    if (key < 0 || key >= S5F_KMER_SPACE || count <= 0) return;
    if (count > 4) count = 4;

    S5FSubstitution *sub = &model->substitution[key];
    sub->count = count;

    /* Normalize to cumulative distribution */
    double total = 0.0;
    for (int i = 0; i < count; i++) total += weights[i];
    if (total <= 0.0) { sub->count = 0; return; }

    double running = 0.0;
    for (int i = 0; i < count; i++) {
        sub->bases[i] = bases[i];
        running += weights[i] / total;
        sub->cumulative[i] = running;
    }
    sub->cumulative[count - 1] = 1.0;  /* avoid bisect edge case */
}

void s5f_set_mutability_str(S5FModel *model, const char *kmer, double value) {
    if (!kmer || strlen(kmer) != 5) return;
    int key = s5f_encode_kmer(
        s5f_base_to_int(kmer[0]), s5f_base_to_int(kmer[1]),
        s5f_base_to_int(kmer[2]), s5f_base_to_int(kmer[3]),
        s5f_base_to_int(kmer[4])
    );
    s5f_set_mutability(model, key, value);
}

void s5f_set_substitution_str(S5FModel *model, const char *kmer,
                                const char *bases, const double *weights, int count) {
    if (!kmer || strlen(kmer) != 5) return;
    int key = s5f_encode_kmer(
        s5f_base_to_int(kmer[0]), s5f_base_to_int(kmer[1]),
        s5f_base_to_int(kmer[2]), s5f_base_to_int(kmer[3]),
        s5f_base_to_int(kmer[4])
    );
    s5f_set_substitution(model, key, bases, weights, count);
}

/* ═══════════════════════════════════════════════════════════════
 * Internal: update affected positions after a mutation
 * ═══════════════════════════════════════════════════════════════ */

static void update_affected(
    const int *ibases,    /* integer-encoded padded bases (length n+4) */
    int n,                /* sequence length (unpadded) */
    int pos,              /* position that was mutated */
    double *weights,      /* per-position weights */
    int *kmer_keys,       /* per-position 5-mer keys */
    const bool *mutable_mask,
    FenwickTree *fenwick,
    const double *mutability
) {
    int lo = pos >= 2 ? pos - 2 : 0;
    int hi = pos + 3 <= n ? pos + 3 : n;

    for (int p = lo; p < hi; p++) {
        int new_key = s5f_encode_kmer(
            ibases[p], ibases[p + 1], ibases[p + 2],
            ibases[p + 3], ibases[p + 4]
        );
        kmer_keys[p] = new_key;
        if (mutable_mask[p]) {
            double old_w = weights[p];
            double new_w = mutability[new_key];
            weights[p] = new_w;
            fenwick_update(fenwick, p, new_w, old_w);
        }
    }
}

/* ── Internal: bisect-left on cumulative array ───────────────── */

static int bisect_left(const double *arr, int n, double value) {
    int lo = 0, hi = n;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (arr[mid] < value) lo = mid + 1;
        else hi = mid;
    }
    return lo < n ? lo : n - 1;
}

/* The per-simulator PCG32 RNG is now plumbed through s5f_mutate's
 * `rng` parameter; the local rand_uniform/rand_01 shadow defns that
 * wrapped libc rand() were removed in T0-5. */

/* ── Codon stop check for productive mode ────────────────────── */

static bool is_stop_codon(char b1, char b2, char b3) {
    /* TAG, TAA, TGA */
    if (b1 == 'T') {
        if (b2 == 'A' && (b3 == 'G' || b3 == 'A')) return true;
        if (b2 == 'G' && b3 == 'A') return true;
    }
    return false;
}

/* ═══════════════════════════════════════════════════════════════
 * s5f_mutate — main entry point
 * ═══════════════════════════════════════════════════════════════ */

void s5f_mutate(const S5FModel *model, ASeq *seq, const SimRecord *rec,
                 struct RngState *rng, S5FResult *result) {
    (void)rec;  /* used in productive mode for anchor positions */

    memset(result, 0, sizeof(*result));

    /* ── Step 1: extract flat arrays from ASeq ───────────────── */
    int n = seq->length;
    if (n == 0) return;

    /* Sample mutation rate */
    double mutation_rate = rng_uniform_lh(rng, model->min_mutation_rate,
                                          model->max_mutation_rate);
    int target_mutations = (int)(mutation_rate * n);

    TRACE("[s5f] rate_range=[%.4f, %.4f], sampled_rate=%.4f, seq_len=%d, target=%d, productive_guard=%s",
          model->min_mutation_rate, model->max_mutation_rate,
          mutation_rate, n, target_mutations,
          model->productive ? "yes" : "no");

    if (target_mutations == 0) {
        result->mutation_rate = mutation_rate;
        TRACE("[s5f] target=0, skipping mutation");
        return;
    }

    /* Padded integer bases: ibases[i+2] = position i in sequence.
     * Pad with N (index 4) on each side for 5-mer context. */
    int ibases[GENAIRR_MAX_SEQ_LEN + 4];
    char bases[GENAIRR_MAX_SEQ_LEN + 4];
    char naive[GENAIRR_MAX_SEQ_LEN];
    bool mutable_mask[GENAIRR_MAX_SEQ_LEN];
    double weights[GENAIRR_MAX_SEQ_LEN];
    int kmer_keys[GENAIRR_MAX_SEQ_LEN];

    /* Node pointer map: to write mutations back to ASeq */
    Nuc *node_map[GENAIRR_MAX_SEQ_LEN];

    ibases[0] = 4; ibases[1] = 4;  /* NN prefix */
    bases[0] = 'N'; bases[1] = 'N';

    int pos = 0;
    for (Nuc *node = seq->head; node && pos < GENAIRR_MAX_SEQ_LEN; node = node->next) {
        char c = node->current;
        ibases[pos + 2] = s5f_base_to_int(c);
        bases[pos + 2] = c;
        naive[pos] = node->germline;
        node_map[pos] = node;

        /* Mutable if V, D, or J segment (not NP1, NP2, etc.)
         * AND not an indel-inserted node (artifact, not biological). */
        Segment seg = node->segment;
        mutable_mask[pos] = (seg == SEG_V || seg == SEG_D || seg == SEG_J)
                            && !(node->flags & NUC_FLAG_INDEL_INS);

        pos++;
    }
    n = pos;  /* actual length */

    ibases[n + 2] = 4; ibases[n + 3] = 4;  /* NN suffix */
    bases[n + 2] = 'N'; bases[n + 3] = 'N';

    /* ── Step 2: compute initial 5-mer keys and weights ──────── */
    const double *mutability = model->mutability;

    for (int i = 0; i < n; i++) {
        int key = s5f_encode_kmer(
            ibases[i], ibases[i + 1], ibases[i + 2],
            ibases[i + 3], ibases[i + 4]
        );
        kmer_keys[i] = key;
        weights[i] = mutable_mask[i] ? mutability[key] : 0.0;
    }

    /* ── Step 3: build Fenwick tree ──────────────────────────── */
    FenwickTree fenwick;
    fenwick_init(&fenwick, weights, n);

    /* ── Step 4: mutation loop ───────────────────────────────── */
    const int MAX_CONSECUTIVE_FAILURES = 500;
    const int max_total_attempts = target_mutations * 100 > 1000
                                   ? target_mutations * 100 : 1000;
    int consecutive_failures = 0;
    int total_attempts = 0;
    int mutation_count = 0;

    while (mutation_count < target_mutations) {
        total_attempts++;
        if (total_attempts >= max_total_attempts) break;

        /* Weighted selection via Fenwick tree */
        double total_w = fenwick.total;
        if (total_w <= 0.0) {
            consecutive_failures++;
            if (consecutive_failures >= MAX_CONSECUTIVE_FAILURES) break;
            continue;
        }

        int sel = fenwick_select(&fenwick, rng_uniform_lh(rng, 0.0, total_w));
        if (sel < 0 || sel >= n) {
            consecutive_failures++;
            if (consecutive_failures >= MAX_CONSECUTIVE_FAILURES) break;
            continue;
        }

        /* Look up substitution for this 5-mer */
        int key = kmer_keys[sel];
        const S5FSubstitution *sub = &model->substitution[key];
        if (sub->count == 0) {
            consecutive_failures++;
            if (consecutive_failures >= MAX_CONSECUTIVE_FAILURES) break;
            continue;
        }

        /* Choose target base from cumulative distribution */
        int chosen = bisect_left(sub->cumulative, sub->count, rng_uniform(rng));
        char mutation_to = sub->bases[chosen];

        /* Skip if substitution matrix drew the same base (silent non-mutation) */
        if (mutation_to == bases[sel + 2]) {
            continue;
        }

        /* ── Productive mode checks ──────────────────────────── */
        if (model->productive) {
            /* Check if mutation would create a stop codon.
             * Check the 3 codons that overlap this position. */
            bool creates_stop = false;
            for (int offset = 0; offset < 3 && !creates_stop; offset++) {
                int codon_start = sel - offset;
                if (codon_start < 0 || codon_start + 2 >= n) continue;
                char c1 = (codon_start == sel)     ? mutation_to : bases[codon_start + 2];
                char c2 = (codon_start + 1 == sel) ? mutation_to : bases[codon_start + 3];
                char c3 = (codon_start + 2 == sel) ? mutation_to : bases[codon_start + 4];
                if (is_stop_codon(c1, c2, c3)) creates_stop = true;
            }
            if (creates_stop) {
                TRACE("[s5f] REJECTED pos=%d %c→%c (would create stop codon)", sel, bases[sel + 2], mutation_to);
                consecutive_failures++;
                if (consecutive_failures >= MAX_CONSECUTIVE_FAILURES) break;
                continue;
            }
        }

        /* ── Apply mutation ──────────────────────────────────── */
        TRACE("[s5f] mutate pos=%d %c→%c (kmer_key=%d, weight=%.4f)",
              sel, bases[sel + 2], mutation_to, key, weights[sel]);
        bases[sel + 2] = mutation_to;
        ibases[sel + 2] = s5f_base_to_int(mutation_to);
        update_affected(ibases, n, sel, weights, kmer_keys,
                        mutable_mask, &fenwick, mutability);

        /* Record if not a reversion to naive */
        if (mutation_to != naive[sel]) {
            /* Check if we already mutated this position */
            bool found = false;
            for (int m = 0; m < mutation_count; m++) {
                if (result->mutations[m].position == sel) {
                    /* Update existing entry */
                    result->mutations[m].to_base = mutation_to;
                    found = true;
                    break;
                }
            }
            if (!found && mutation_count < GENAIRR_MAX_SEQ_LEN) {
                result->mutations[mutation_count].position = sel;
                result->mutations[mutation_count].from_base = naive[sel];
                result->mutations[mutation_count].to_base = mutation_to;
                mutation_count++;
            }
        } else {
            /* Reverted to naive — remove from results */
            for (int m = 0; m < mutation_count; m++) {
                if (result->mutations[m].position == sel) {
                    result->mutations[m] = result->mutations[mutation_count - 1];
                    mutation_count--;
                    break;
                }
            }
        }

        /* Periodic rebuild */
        if (fenwick_needs_rebuild(&fenwick)) {
            fenwick_rebuild(&fenwick, weights);
        }

        consecutive_failures = 0;
    }

    /* ── Step 5: write mutations back to ASeq ────────────────── */
    for (int m = 0; m < mutation_count; m++) {
        int p = result->mutations[m].position;
        if (p >= 0 && p < n) {
            
            aseq_mutate(seq, node_map[p], result->mutations[m].to_base, NUC_FLAG_MUTATED);
        }
    }

    result->count = mutation_count;
    result->mutation_rate = (double)mutation_count / (double)n;

    TRACE("[s5f] done: %d mutations applied (%.4f effective rate), %d total attempts",
          mutation_count, result->mutation_rate, total_attempts);

    fenwick_destroy(&fenwick);
}
