/**
 * rand_util.c — Random number utilities.
 *
 * Provides beta distribution sampling, uniform doubles, and random
 * nucleotide generation used by corruption/simulation ops.
 */

#include "genairr/rand_util.h"

/* M_PI is not part of the C standard; MSVC requires _USE_MATH_DEFINES */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <math.h>
#include <stdlib.h>

/* ── Seed ─────────────────────────────────────────────────────────── */

void rand_seed(unsigned int seed) {
    srand(seed);
}

/* ── Uniform [0,1) ───────────────────────────────────────────────── */

double rand_uniform(void) {
    return (double)rand() / ((double)RAND_MAX + 1.0);
}

/* ── Gamma via Marsaglia–Tsang ──────────────────────────────────── */

/**
 * Generate a standard normal using Box-Muller.
 */
static double rand_normal(void) {
    double u1 = rand_uniform();
    double u2 = rand_uniform();
    if (u1 < 1e-15) u1 = 1e-15;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/**
 * Gamma(alpha, 1) for alpha >= 1 using Marsaglia–Tsang method.
 */
static double rand_gamma_ge1(double alpha) {
    double d = alpha - 1.0 / 3.0;
    double c = 1.0 / sqrt(9.0 * d);

    for (;;) {
        double x, v;
        do {
            x = rand_normal();
            v = 1.0 + c * x;
        } while (v <= 0.0);

        v = v * v * v;
        double u = rand_uniform();

        if (u < 1.0 - 0.0331 * (x * x) * (x * x))
            return d * v;

        if (log(u) < 0.5 * x * x + d * (1.0 - v + log(v)))
            return d * v;
    }
}

double rand_gamma(double alpha) {
    if (alpha < 1.0) {
        /* Use Gamma(alpha+1) * U^(1/alpha) */
        double u = rand_uniform();
        if (u < 1e-15) u = 1e-15;
        return rand_gamma_ge1(alpha + 1.0) * pow(u, 1.0 / alpha);
    }
    return rand_gamma_ge1(alpha);
}

/* ── Beta distribution ──────────────────────────────────────────── */

double rand_beta(double alpha, double beta) {
    double x = rand_gamma(alpha);
    double y = rand_gamma(beta);
    double sum = x + y;
    if (sum < 1e-15) return 0.5;
    return x / sum;
}

/* ── Random nucleotide string ───────────────────────────────────── */

void rand_nucleotides(char *buf, int len) {
    static const char bases[] = "ACGT";
    for (int i = 0; i < len; i++) {
        buf[i] = bases[rand() % 4];
    }
    buf[len] = '\0';
}
