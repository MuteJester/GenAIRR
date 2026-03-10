/**
 * allele.h — Allele representation.
 *
 * An Allele is an immutable reference record loaded from the DataConfig.
 * It stores the ungapped germline sequence, its length, the conserved
 * anchor position, gene/family names, and segment type.
 */

#ifndef GENAIRR_ALLELE_H
#define GENAIRR_ALLELE_H

#include "types.h"
#include <string.h>

typedef struct {
    char     name[GENAIRR_MAX_ALLELE_NAME];     /* e.g. "IGHV1-2*01"       */
    char     gene[GENAIRR_MAX_ALLELE_NAME];     /* e.g. "IGHV1-2"          */
    char     family[GENAIRR_MAX_ALLELE_NAME];   /* e.g. "IGHV1"            */
    char     seq[GENAIRR_MAX_ALLELE_SEQ];       /* ungapped nucleotide seq */
    uint16_t length;                             /* strlen(seq)             */
    uint16_t anchor;                             /* conserved anchor pos    */
    Segment  segment_type;                       /* SEG_V, SEG_D, SEG_J, SEG_C */
} Allele;

/* ── Allele pool: a flat array of alleles for one segment type ── */

typedef struct {
    Allele  *alleles;
    int      count;
    int      capacity;
} AllelePool;

AllelePool  allele_pool_create(int initial_capacity);
void        allele_pool_destroy(AllelePool *pool);
void        allele_pool_add(AllelePool *pool, const Allele *allele);
const Allele *allele_pool_random(const AllelePool *pool);

/** Find the index of an allele by name. Returns -1 if not found. */
static inline int allele_pool_find_index(const AllelePool *pool,
                                          const char *name) {
    for (int i = 0; i < pool->count; i++) {
        if (strcmp(pool->alleles[i].name, name) == 0) return i;
    }
    return -1;
}

/* ── Allele restriction (for locking specific alleles) ─────────── */

#define GENAIRR_MAX_LOCKED 64

typedef struct {
    int    indices[GENAIRR_MAX_LOCKED];   /* indices into the AllelePool */
    int    count;                          /* number of locked alleles   */
    bool   active;                        /* false = no restriction     */
} AlleleRestriction;

/**
 * Pick a random allele, respecting restrictions if active.
 * If restricted, picks uniformly from the allowed indices.
 * If not restricted, picks uniformly from the full pool.
 */
const Allele *allele_pool_pick(const AllelePool *pool,
                                const AlleleRestriction *restriction);

#endif /* GENAIRR_ALLELE_H */
