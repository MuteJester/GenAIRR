/**
 * allele_pool.c — Allele pool management.
 *
 * Simple growable array of Allele records. Pools are built once
 * during config loading and remain immutable during simulation.
 */

#include "genairr/allele.h"
#include <stdlib.h>
#include <string.h>

AllelePool allele_pool_create(int initial_capacity) {
    AllelePool pool = {0};
    pool.capacity = initial_capacity > 0 ? initial_capacity : 16;
    pool.alleles = (Allele *)calloc(pool.capacity, sizeof(Allele));
    pool.count = 0;
    return pool;
}

void allele_pool_destroy(AllelePool *pool) {
    if (pool->alleles) {
        free(pool->alleles);
        pool->alleles = NULL;
    }
    pool->count = 0;
    pool->capacity = 0;
}

void allele_pool_add(AllelePool *pool, const Allele *allele) {
    if (pool->count >= pool->capacity) {
        int new_cap = pool->capacity > 0 ? pool->capacity * 2 : 16;
        pool->alleles = (Allele *)realloc(
            pool->alleles, new_cap * sizeof(Allele)
        );
        pool->capacity = new_cap;
    }
    memcpy(&pool->alleles[pool->count], allele, sizeof(Allele));
    pool->count++;
}

const Allele *allele_pool_random(const AllelePool *pool) {
    if (pool->count == 0) return NULL;
    int idx = rand() % pool->count;
    return &pool->alleles[idx];
}

const Allele *allele_pool_pick(const AllelePool *pool,
                                const AlleleRestriction *restriction) {
    if (pool->count == 0) return NULL;

    if (restriction && restriction->active && restriction->count > 0) {
        int pick = rand() % restriction->count;
        int idx = restriction->indices[pick];
        if (idx >= 0 && idx < pool->count) {
            return &pool->alleles[idx];
        }
    }

    return allele_pool_random(pool);
}
