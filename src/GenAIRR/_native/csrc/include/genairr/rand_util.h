/**
 * rand_util.h — Random number utilities.
 */

#ifndef GENAIRR_RAND_UTIL_H
#define GENAIRR_RAND_UTIL_H

/** Seed the RNG. */
void rand_seed(unsigned int seed);

/** Uniform random double in [0, 1). */
double rand_uniform(void);

/** Gamma(alpha, 1) random variate. */
double rand_gamma(double alpha);

/** Beta(alpha, beta) random variate in [0, 1]. */
double rand_beta(double alpha, double beta);

/** Fill buf with len random ACGT nucleotides, NUL-terminated. */
void rand_nucleotides(char *buf, int len);

#endif /* GENAIRR_RAND_UTIL_H */
