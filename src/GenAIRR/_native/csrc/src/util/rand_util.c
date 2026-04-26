/**
 * rand_util.c — Per-instance PCG32 RNG and sampling utilities.
 *
 * See rand_util.h for the API contract. The PCG32 reference impl is
 * Melissa O'Neill's public-domain version (pcg-random.org); the
 * Marsaglia-Tsang gamma sampler and Box-Muller normal are textbook.
 */

#include "genairr/rand_util.h"

/* M_PI is not part of the C standard; MSVC requires _USE_MATH_DEFINES */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <math.h>

/* ── PCG32 core ───────────────────────────────────────────────── */

#define PCG_MULT 6364136223846793005ULL

void rng_seed(RngState *r, uint64_t seed, uint64_t stream) {
    r->state = 0;
    r->inc   = (stream << 1u) | 1u;     /* must be odd */
    rng_next_u32(r);
    r->state += seed;
    rng_next_u32(r);
}

uint32_t rng_next_u32(RngState *r) {
    uint64_t old = r->state;
    r->state = old * PCG_MULT + (r->inc | 1u);
    uint32_t xorshifted = (uint32_t)(((old >> 18u) ^ old) >> 27u);
    uint32_t rot = (uint32_t)(old >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31u));
}

/* ── Uniform [0, 1) ───────────────────────────────────────────── */

double rng_uniform(RngState *r) {
    /* 2^32 in the denominator; gives ~9.3e-10 spacing. */
    return (double)rng_next_u32(r) * (1.0 / 4294967296.0);
}

double rng_uniform_lh(RngState *r, double lo, double hi) {
    return lo + rng_uniform(r) * (hi - lo);
}

/* ── Bounded integer (Lemire's debiased modulo) ───────────────── */

uint32_t rng_range(RngState *r, uint32_t bound) {
    if (bound == 0) return 0;
    uint64_t m = (uint64_t)rng_next_u32(r) * (uint64_t)bound;
    uint32_t l = (uint32_t)m;
    if (l < bound) {
        uint32_t t = (uint32_t)(-bound) % bound;
        while (l < t) {
            m = (uint64_t)rng_next_u32(r) * (uint64_t)bound;
            l = (uint32_t)m;
        }
    }
    return (uint32_t)(m >> 32);
}

/* ── Standard normal via Box-Muller ──────────────────────────── */

double rng_normal(RngState *r) {
    double u1 = rng_uniform(r);
    double u2 = rng_uniform(r);
    if (u1 < 1e-15) u1 = 1e-15;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* ── Gamma(alpha, 1) via Marsaglia-Tsang ─────────────────────── */

static double rng_gamma_ge1(RngState *r, double alpha) {
    double d = alpha - 1.0 / 3.0;
    double c = 1.0 / sqrt(9.0 * d);

    for (;;) {
        double x, v;
        do {
            x = rng_normal(r);
            v = 1.0 + c * x;
        } while (v <= 0.0);

        v = v * v * v;
        double u = rng_uniform(r);

        if (u < 1.0 - 0.0331 * (x * x) * (x * x))
            return d * v;

        if (log(u) < 0.5 * x * x + d * (1.0 - v + log(v)))
            return d * v;
    }
}

double rng_gamma(RngState *r, double alpha) {
    if (alpha < 1.0) {
        double u = rng_uniform(r);
        if (u < 1e-15) u = 1e-15;
        return rng_gamma_ge1(r, alpha + 1.0) * pow(u, 1.0 / alpha);
    }
    return rng_gamma_ge1(r, alpha);
}

/* ── Beta(alpha, beta) ───────────────────────────────────────── */

double rng_beta(RngState *r, double alpha, double beta) {
    double x = rng_gamma(r, alpha);
    double y = rng_gamma(r, beta);
    double sum = x + y;
    if (sum < 1e-15) return 0.5;
    return x / sum;
}

/* ── Random nucleotide string ────────────────────────────────── */

void rng_nucleotides(RngState *r, char *buf, int len) {
    static const char bases[] = "ACGT";
    for (int i = 0; i < len; i++) {
        buf[i] = bases[rng_range(r, 4)];
    }
    buf[len] = '\0';
}
