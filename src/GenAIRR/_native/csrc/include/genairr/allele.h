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

/* Forward declaration: full TrimDist definition lives in sim_config.h.
 * Allele only stores non-owning pointers, so we don't need the full
 * type here. This breaks the include cycle (sim_config.h → allele.h). */
struct TrimDist;

typedef struct {
    char     name[GENAIRR_MAX_ALLELE_NAME];     /* e.g. "IGHV1-2*01"       */
    char     gene[GENAIRR_MAX_ALLELE_NAME];     /* e.g. "IGHV1-2"          */
    char     family[GENAIRR_MAX_ALLELE_NAME];   /* e.g. "IGHV1"            */
    char     seq[GENAIRR_MAX_ALLELE_SEQ];       /* ungapped nucleotide seq */
    uint16_t length;                             /* strlen(seq)             */
    int16_t  anchor;                             /* conserved anchor pos.
                                                  * Signed: -1 (or 0) means
                                                  * "no anchor" (anchorless
                                                  * V/J pseudogene). Real
                                                  * anchor positions are
                                                  * always > 0. D alleles
                                                  * leave this 0 (anchor is
                                                  * meaningless for D).    */
    Segment  segment_type;                       /* SEG_V, SEG_D, SEG_J, SEG_C */

    /* Per-allele trim distributions, looked up by (family, gene) and
     * cached at gdc_populate_sim_config time. NULL means the allele
     * has no specific match — trim.c then falls back to the legacy
     * cfg->v_trim_3 / d_trim_5 / d_trim_3 / j_trim_5 globals.
     *   - V uses trim_dist_3 only
     *   - D uses both
     *   - J uses trim_dist_5 only
     *   - C uses neither */
    const struct TrimDist *trim_dist_5;
    const struct TrimDist *trim_dist_3;
} Allele;

/* True if the allele has a usable anchor (V Cys / J W or F).
 * Treats both `<= 0` sentinels as anchorless: GDC writes -1 for
 * anchorless V/J alleles, while the legacy embedded test data path
 * leaves anchor=0 for the same case. Real V/J anchors are always at
 * a positive germline position (V Cys is at IMGT pos 104; J anchor
 * is mid-allele). D alleles always have anchor=0 and never call
 * this — assemble.c passes a literal -1 for D. */
static inline bool allele_has_anchor(const Allele *a) {
    return a && a->anchor > 0;
}

/* ── Allele pool: a flat array of alleles for one segment type ── */

typedef struct {
    Allele  *alleles;
    int      count;
    int      capacity;
} AllelePool;

/* Forward decl for the per-simulator RNG (defined in rand_util.h). */
struct RngState;

AllelePool  allele_pool_create(int initial_capacity);
void        allele_pool_destroy(AllelePool *pool);
void        allele_pool_add(AllelePool *pool, const Allele *allele);
const Allele *allele_pool_random(const AllelePool *pool, struct RngState *rng);

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
                                const AlleleRestriction *restriction,
                                struct RngState *rng);

#endif /* GENAIRR_ALLELE_H */
