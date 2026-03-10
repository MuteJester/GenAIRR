/**
 * s5f.h — S5F (Somatic 5-mer Frequency) mutation model.
 *
 * Context-dependent 5-mer mutation model: each position's mutation
 * probability depends on its surrounding 2+1+2 nucleotide context.
 *
 * Performance internals:
 *   - Integer-encoded 5-mers (base-5, key space 3125)
 *   - Flat arrays for mutability and substitution lookup
 *   - Fenwick tree for O(log N) weighted position selection
 *   - NP-aware: only V/D/J positions are mutable (from ASeq segment tags)
 *
 * This mirrors the Python S5F class in mutation/s5f.py.
 */

#ifndef GENAIRR_S5F_H
#define GENAIRR_S5F_H

#include "types.h"
#include "aseq.h"
#include "pipeline.h"
#include <stdbool.h>

/* ── 5-mer integer encoding ──────────────────────────────────── */

#define S5F_KMER_SPACE  3125   /* 5^5 */

/* Base encoding: A=0, C=1, G=2, T=3, N=4 */
static inline int s5f_base_to_int(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return 4;  /* N or unknown */
    }
}

static inline char s5f_int_to_base(int i) {
    static const char bases[] = "ACGTN";
    return (i >= 0 && i < 5) ? bases[i] : 'N';
}

/* Encode 5 integer bases into a single key (base-5 positional). */
static inline int s5f_encode_kmer(int b0, int b1, int b2, int b3, int b4) {
    return b0 * 625 + b1 * 125 + b2 * 25 + b3 * 5 + b4;
}

/* ── Fenwick Tree ────────────────────────────────────────────── */

/**
 * Binary Indexed Tree for O(log N) weighted random selection.
 *
 * O(N) construction, O(log N) point update, O(log N) select.
 * Includes periodic rebuild to prevent floating-point drift.
 */
typedef struct {
    double *tree;          /* 1-indexed tree array (size n+1) */
    int     n;             /* number of elements */
    double  total;         /* sum of all weights */
    int     update_count;  /* number of updates since last rebuild */
} FenwickTree;

#define FENWICK_REBUILD_INTERVAL  200

void   fenwick_init(FenwickTree *ft, const double *weights, int n);
void   fenwick_destroy(FenwickTree *ft);
void   fenwick_update(FenwickTree *ft, int i, double new_w, double old_w);
void   fenwick_rebuild(FenwickTree *ft, const double *weights);
int    fenwick_select(const FenwickTree *ft, double r);

static inline bool fenwick_needs_rebuild(const FenwickTree *ft) {
    return ft->update_count >= FENWICK_REBUILD_INTERVAL;
}

/* ── Substitution entry ──────────────────────────────────────── */

/**
 * Pre-cached substitution data for one 5-mer context.
 * Contains sorted cumulative probabilities for bisect-based selection.
 */
typedef struct {
    char   bases[4];       /* possible target bases */
    double cumulative[4];  /* cumulative probability [0..1] */
    int    count;          /* number of valid entries (0 = no data) */
} S5FSubstitution;

/* ── S5F Model ───────────────────────────────────────────────── */

typedef struct {
    double          mutability[S5F_KMER_SPACE];    /* per-5-mer mutation likelihood */
    S5FSubstitution substitution[S5F_KMER_SPACE];  /* per-5-mer substitution targets */

    double min_mutation_rate;
    double max_mutation_rate;
    bool   productive;
} S5FModel;

/**
 * Initialize an S5F model with default parameters.
 * Mutability and substitution tables are zeroed — must be populated
 * via s5f_set_mutability() and s5f_set_substitution().
 */
void s5f_model_init(S5FModel *model,
                     double min_rate, double max_rate,
                     bool productive);

/**
 * Set mutability for a single 5-mer key.
 */
static inline void s5f_set_mutability(S5FModel *model, int key, double value) {
    if (key >= 0 && key < S5F_KMER_SPACE) {
        model->mutability[key] = value;
    }
}

/**
 * Set substitution data for a single 5-mer key.
 *
 * @param bases   Array of target bases (e.g., "ACG" for 3 alternatives).
 * @param weights Array of raw weights (will be normalized to cumulative).
 * @param count   Number of alternatives.
 */
void s5f_set_substitution(S5FModel *model, int key,
                           const char *bases, const double *weights, int count);

/**
 * Set mutability from a string 5-mer key.
 * Convenience: encodes the string and calls s5f_set_mutability().
 */
void s5f_set_mutability_str(S5FModel *model, const char *kmer, double value);

/**
 * Set substitution from a string 5-mer key.
 * Convenience: encodes the string and calls s5f_set_substitution().
 */
void s5f_set_substitution_str(S5FModel *model, const char *kmer,
                                const char *bases, const double *weights, int count);

/* ── Mutation result ─────────────────────────────────────────── */

/** Per-position mutation record. */
typedef struct {
    int  position;     /* position in the assembled sequence */
    char from_base;    /* original (germline) base */
    char to_base;      /* mutated base */
} S5FMutation;

/** Result of applying S5F mutations. */
typedef struct {
    S5FMutation mutations[GENAIRR_MAX_SEQ_LEN];
    int         count;          /* number of mutations applied */
    double      mutation_rate;  /* actual = count / seq_length */
} S5FResult;

/* ── Main mutation function ──────────────────────────────────── */

/**
 * Apply S5F somatic hypermutation to an ASeq.
 *
 * Mutates nucleotides IN-PLACE on the ASeq (sets n->current and
 * NUC_FLAG_MUTATED). Only positions with segment V/D/J are eligible.
 *
 * @param model  The S5F model with loaded mutability/substitution data.
 * @param seq    The annotated sequence to mutate.
 * @param rec    Simulation record (junction anchors for productive mode).
 * @param result Output: mutation details.
 */
void s5f_mutate(const S5FModel *model, ASeq *seq, const SimRecord *rec,
                 S5FResult *result);

#endif /* GENAIRR_S5F_H */
