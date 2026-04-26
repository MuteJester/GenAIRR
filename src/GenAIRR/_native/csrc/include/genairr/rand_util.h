/**
 * rand_util.h — Per-instance PCG32 RNG and sampling utilities.
 *
 * Replaces the previous process-global libc rand()/srand() RNG.
 * Each GenAIRRSimulator owns an `RngState` and exposes it through
 * SimConfig.rng so step functions can draw without touching shared
 * state. Two simulators in the same process can run independently;
 * with the same seed they produce byte-identical output.
 *
 * Algorithm: PCG32 (Permuted Congruential Generator, M.E. O'Neill,
 * 2014). 64-bit state, 64-bit increment, 32-bit output. Public
 * domain. ~30 lines of impl.
 */

#ifndef GENAIRR_RAND_UTIL_H
#define GENAIRR_RAND_UTIL_H

#include <stdint.h>

typedef struct RngState {
    uint64_t state;   /* RNG state; advances every draw */
    uint64_t inc;     /* odd increment; selects the stream */
} RngState;

/* Seed an RngState. `seed` is the user-supplied seed; `stream` is a
 * stream selector (0 is fine for single-stream use). After seeding,
 * the same (seed, stream) always reproduces the same sequence. */
void     rng_seed(RngState *r, uint64_t seed, uint64_t stream);

/* Raw 32-bit draw. */
uint32_t rng_next_u32(RngState *r);

/* Uniform in [0, 1). 32-bit precision. */
double   rng_uniform(RngState *r);

/* Uniform in [lo, hi). */
double   rng_uniform_lh(RngState *r, double lo, double hi);

/* Uniform integer in [0, bound). Unbiased via Lemire's method. */
uint32_t rng_range(RngState *r, uint32_t bound);

/* Standard normal (Box-Muller, single sample per call). */
double   rng_normal(RngState *r);

/* Gamma(alpha, 1) via Marsaglia-Tsang. */
double   rng_gamma(RngState *r, double alpha);

/* Beta(alpha, beta) in [0, 1]. */
double   rng_beta(RngState *r, double alpha, double beta);

/* Fill `buf` with `len` uniform random A/C/G/T bases, NUL-terminated. */
void     rng_nucleotides(RngState *r, char *buf, int len);

#endif /* GENAIRR_RAND_UTIL_H */
