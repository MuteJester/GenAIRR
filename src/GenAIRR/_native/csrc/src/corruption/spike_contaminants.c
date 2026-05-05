/**
 * spike_contaminants.c — Cross-sample / phiX contamination.
 *
 * Replaces the entire sequence with a contaminant and zeros out all
 * AIRR annotations. Two types: "random" (random ACGT) and "phix"
 * (fragment of phiX174 genome).
 */

#include "genairr/pipeline.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ~600bp of phiX174 genome */
static const char PHIX_FRAGMENT[] =
    "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAA"
    "AATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTG"
    "CTGGCGGAAAATGAGAAAATTCGACCTATCCTTGCGCAGCTCGAGAAGCTCTTACTTTGCG"
    "ACCTTTCGCCATCAACTAACGATTCTGTCAAAAACTGACGCGTTGGATGAGGAGAAGTGGC"
    "TTAATATGCTTGGCACGTTCGTCAAGGACTGGTTTAGATATGAGTCACATTTTGTTCATGG"
    "TAGAGATTCTCTTGTTGACATTTTAAAAGAGCGTGGATTACTATCTGAGTCCGATGCTGTT"
    "CAACCACTAATAGGTAAGAAATCATGAGTCAAGTTACTGAACAATCCGTACGTTTCCAGAC"
    "CGCTTTGGCCTCTATTAAGCTCATTCAGGCTTCTGCCGTTTTGGATTTAACCGAAGATGAT"
    "TTCGATTTTCTGACGAGTAACAAAGTTTGGATTGCTACTGACCGCTCTCGTGCTCGTCGCT"
    "GCGTTGAGGCTTGCGTTTATGGTACGCTGGACTTTGTGGGATACCCTCGCTTTCCTGCTCC";

static int sample_contaminant_payload(const SimConfig *cfg,
                                      char *seq_buf,
                                      size_t seq_buf_cap,
                                      char *type_buf,
                                      size_t type_buf_cap,
                                      int *fragment_offset) {
    int cont_len = 300 + (int)rng_range(cfg->rng, 200);  /* 300-500 bp */

    if (cfg->contaminant_type == 1) {
        int phix_len = (int)strlen(PHIX_FRAGMENT);
        int start = (int)rng_range(cfg->rng,
            (uint32_t)(phix_len > cont_len ? phix_len - cont_len : 1));
        if (start + cont_len > phix_len) cont_len = phix_len - start;
        if ((size_t)cont_len >= seq_buf_cap) cont_len = (int)seq_buf_cap - 1;
        memcpy(seq_buf, PHIX_FRAGMENT + start, (size_t)cont_len);
        seq_buf[cont_len] = '\0';
        strncpy(type_buf, "phix", type_buf_cap);
        if (type_buf_cap > 0) type_buf[type_buf_cap - 1] = '\0';
        *fragment_offset = start;
        return cont_len;
    }

    if ((size_t)cont_len >= seq_buf_cap) cont_len = (int)seq_buf_cap - 1;
    rng_nucleotides(cfg->rng, seq_buf, cont_len);
    seq_buf[cont_len] = '\0';
    strncpy(type_buf, "random", type_buf_cap);
    if (type_buf_cap > 0) type_buf[type_buf_cap - 1] = '\0';
    *fragment_offset = -1;
    return cont_len;
}

void step_spike_contaminants(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    if (rng_uniform(cfg->rng) >= cfg->contamination_prob) return;

    TRACE("[contaminant] TRIGGERED — replacing sequence with contaminant (prob=%.4f)", cfg->contamination_prob);

    char cont_buf[512];
    char cont_type[32];
    int fragment_offset = -1;
    int cont_len = sample_contaminant_payload(
        cfg, cont_buf, sizeof(cont_buf),
        cont_type, sizeof(cont_type),
        &fragment_offset);

    if (simcfg_productive_only(cfg)) {
        TRACE("[contaminant] skipped productive-unsafe contamination type=%s, length=%d",
              cont_type, cont_len);
        return;
    }

    /* Clear the entire sequence */
    aseq_reset(seq);

    /* Mark as contaminant */
    rec->is_contaminant = true;
    rec->productive = false;
    rec->stop_codon = false;
    rec->vj_in_frame = false;
    rec->v_allele = NULL;
    rec->d_allele = NULL;
    rec->j_allele = NULL;
    rec->c_allele = NULL;
    snprintf(rec->note, sizeof(rec->note), "contaminant");
    strncpy(rec->contaminant_type, cont_type, sizeof(rec->contaminant_type));
    rec->contaminant_type[sizeof(rec->contaminant_type) - 1] = '\0';
    aseq_append_segment(seq, cont_buf, cont_len, SEG_ADAPTER, 0, -1);
    if (fragment_offset >= 0) {
        TRACE("[contaminant] type=phiX, length=%d, fragment_offset=%d",
              cont_len, fragment_offset);
    } else {
        TRACE("[contaminant] type=random, length=%d", cont_len);
    }
}
