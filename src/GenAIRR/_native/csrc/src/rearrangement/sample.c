/**
 * sample.c — Allele sampling step functions.
 *
 * Each step picks a random allele from the appropriate pool
 * in the SimConfig and stores it in the SimRecord.
 */

#include "genairr/pipeline.h"
#include "genairr/trace.h"
#include <stdlib.h>

void step_sample_v(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)seq;
    rec->v_allele = allele_pool_pick(&cfg->v_alleles, &cfg->v_restriction);
    TRACE("[sample_v] %s (len=%d, anchor=%d, pool=%d, locked=%s)",
          rec->v_allele ? rec->v_allele->name : "NULL",
          rec->v_allele ? rec->v_allele->length : 0,
          rec->v_allele ? rec->v_allele->anchor : 0,
          cfg->v_alleles.count,
          cfg->v_restriction.active ? "yes" : "no");
}

void step_sample_d(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)seq;
    if (cfg->d_alleles.count > 0) {
        rec->d_allele = allele_pool_pick(&cfg->d_alleles, &cfg->d_restriction);
        TRACE("[sample_d] %s (len=%d, pool=%d, locked=%s)",
              rec->d_allele ? rec->d_allele->name : "NULL",
              rec->d_allele ? rec->d_allele->length : 0,
              cfg->d_alleles.count,
              cfg->d_restriction.active ? "yes" : "no");
    } else {
        rec->d_allele = NULL;
        TRACE("[sample_d] no D alleles in pool (VJ chain)");
    }
}

void step_sample_j(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)seq;
    rec->j_allele = allele_pool_pick(&cfg->j_alleles, &cfg->j_restriction);
    TRACE("[sample_j] %s (len=%d, anchor=%d, pool=%d, locked=%s)",
          rec->j_allele ? rec->j_allele->name : "NULL",
          rec->j_allele ? rec->j_allele->length : 0,
          rec->j_allele ? rec->j_allele->anchor : 0,
          cfg->j_alleles.count,
          cfg->j_restriction.active ? "yes" : "no");
}

void step_sample_c(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)seq;
    if (cfg->c_alleles.count > 0) {
        rec->c_allele = allele_pool_pick(&cfg->c_alleles, &cfg->c_restriction);
        TRACE("[sample_c] %s (pool=%d)",
              rec->c_allele ? rec->c_allele->name : "NULL",
              cfg->c_alleles.count);
    } else {
        rec->c_allele = NULL;
    }
}
