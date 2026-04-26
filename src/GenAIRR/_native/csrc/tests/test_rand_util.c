/**
 * test_rand_util.c — PCG32 RNG unit tests (T0-5).
 *
 * Verifies the per-instance PCG32 implementation:
 *   - Determinism: same seed → identical sequence
 *   - Independence: two states never collide on the same seed
 *   - Distribution: uniform mean ≈ 0.5; rng_range bounded
 *   - Independence under interleaving: two RngStates with different
 *     seeds, drawn from in alternating order, each produce the same
 *     sequence as if drawn alone.
 */

#include "genairr/rand_util.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define RUN_TEST(fn) do {                                  \
    printf("  %-55s ", #fn);                               \
    fn();                                                  \
    printf("PASS\n");                                      \
} while (0)

/* Determinism: two RngStates seeded with the same (seed, stream)
 * must produce a bit-identical 1000-draw sequence. */
static void test_pcg32_determinism(void) {
    RngState a, b;
    rng_seed(&a, 42, 0);
    rng_seed(&b, 42, 0);
    for (int i = 0; i < 1000; i++) {
        uint32_t va = rng_next_u32(&a);
        uint32_t vb = rng_next_u32(&b);
        assert(va == vb);
    }
}

/* Different seeds must produce different sequences (probabilistically
 * certain across 100 draws). */
static void test_pcg32_different_seeds(void) {
    RngState a, b;
    rng_seed(&a, 1, 0);
    rng_seed(&b, 2, 0);
    int diffs = 0;
    for (int i = 0; i < 100; i++) {
        if (rng_next_u32(&a) != rng_next_u32(&b)) diffs++;
    }
    assert(diffs > 90);  /* nearly all should differ */
}

/* Different streams (same seed, different stream id) must also
 * produce different sequences. */
static void test_pcg32_different_streams(void) {
    RngState a, b;
    rng_seed(&a, 42, 0);
    rng_seed(&b, 42, 1);
    int diffs = 0;
    for (int i = 0; i < 100; i++) {
        if (rng_next_u32(&a) != rng_next_u32(&b)) diffs++;
    }
    assert(diffs > 90);
}

/* Independence under interleaving: drawing alternately from two
 * states must NOT change either state's sequence. Two states with
 * the same seed both stand-alone and interleaved must produce the
 * same stream — this is the property that makes the per-simulator
 * RNG correct under multi-simulator workflows. */
static void test_pcg32_interleaving_independence(void) {
    RngState a_solo, b_solo, a_inter, b_inter;
    rng_seed(&a_solo,  111, 0);
    rng_seed(&b_solo,  222, 0);
    rng_seed(&a_inter, 111, 0);
    rng_seed(&b_inter, 222, 0);

    for (int i = 0; i < 500; i++) {
        uint32_t s_a = rng_next_u32(&a_solo);
        uint32_t i_a = rng_next_u32(&a_inter);
        assert(s_a == i_a);

        uint32_t s_b = rng_next_u32(&b_solo);
        uint32_t i_b = rng_next_u32(&b_inter);
        assert(s_b == i_b);
    }
}

/* rng_uniform must produce values in [0, 1) with mean ≈ 0.5. */
static void test_uniform_distribution(void) {
    RngState r;
    rng_seed(&r, 7, 0);
    double sum = 0.0;
    int n = 50000;
    for (int i = 0; i < n; i++) {
        double u = rng_uniform(&r);
        assert(u >= 0.0 && u < 1.0);
        sum += u;
    }
    double mean = sum / n;
    assert(fabs(mean - 0.5) < 0.02);
}

/* rng_range(r, bound) must always be in [0, bound). With 100k draws
 * across bound=10, every bucket should be hit and the spread should
 * not collapse to 0 or 9 only. */
static void test_range_bounded(void) {
    RngState r;
    rng_seed(&r, 13, 0);
    int counts[10] = {0};
    int n = 100000;
    for (int i = 0; i < n; i++) {
        uint32_t v = rng_range(&r, 10);
        assert(v < 10);
        counts[v]++;
    }
    for (int i = 0; i < 10; i++) {
        /* Each bucket should be ~10000; allow ±15%. */
        assert(counts[i] >= 8500 && counts[i] <= 11500);
    }
}

/* rng_nucleotides must fill exactly `len` ACGT bases plus NUL. */
static void test_nucleotides(void) {
    RngState r;
    rng_seed(&r, 23, 0);
    char buf[64];
    rng_nucleotides(&r, buf, 50);
    assert(buf[50] == '\0');
    for (int i = 0; i < 50; i++) {
        assert(buf[i] == 'A' || buf[i] == 'C' ||
               buf[i] == 'G' || buf[i] == 'T');
    }
}

/* rng_uniform_lh must respect [lo, hi) bounds. */
static void test_uniform_lh(void) {
    RngState r;
    rng_seed(&r, 29, 0);
    for (int i = 0; i < 1000; i++) {
        double v = rng_uniform_lh(&r, -1.5, 2.5);
        assert(v >= -1.5 && v < 2.5);
    }
}

/* rng_normal: mean ≈ 0, std ≈ 1 over many draws. */
static void test_normal_distribution(void) {
    RngState r;
    rng_seed(&r, 31, 0);
    int n = 50000;
    double sum = 0.0, sq = 0.0;
    for (int i = 0; i < n; i++) {
        double x = rng_normal(&r);
        sum += x;
        sq  += x * x;
    }
    double mean = sum / n;
    double var  = sq / n - mean * mean;
    assert(fabs(mean) < 0.05);
    assert(fabs(sqrt(var) - 1.0) < 0.05);
}

int main(void) {
    printf("test_rand_util: PCG32 RNG and sampler tests\n");
    printf("═══════════════════════════════════════════════════════════\n\n");

    printf("Core PCG32:\n");
    RUN_TEST(test_pcg32_determinism);
    RUN_TEST(test_pcg32_different_seeds);
    RUN_TEST(test_pcg32_different_streams);
    RUN_TEST(test_pcg32_interleaving_independence);

    printf("\nSamplers:\n");
    RUN_TEST(test_uniform_distribution);
    RUN_TEST(test_range_bounded);
    RUN_TEST(test_nucleotides);
    RUN_TEST(test_uniform_lh);
    RUN_TEST(test_normal_distribution);

    printf("\nAll PCG32 tests passed.\n");
    return 0;
}
