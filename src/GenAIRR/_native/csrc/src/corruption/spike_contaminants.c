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

void step_spike_contaminants(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    if (rand_uniform() >= cfg->contamination_prob) return;

    TRACE("[contaminant] TRIGGERED — replacing sequence with contaminant (prob=%.4f)", cfg->contamination_prob);

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

    /* Generate contaminant sequence */
    int cont_len = 300 + rand() % 200;  /* 300-500 bp */

    if (cfg->contaminant_type == 1) {
        /* phiX: random substring of the fragment */
        int phix_len = (int)strlen(PHIX_FRAGMENT);
        int start = rand() % (phix_len > cont_len ? phix_len - cont_len : 1);
        if (start + cont_len > phix_len) cont_len = phix_len - start;
        strncpy(rec->contaminant_type, "phix", sizeof(rec->contaminant_type));
        aseq_append_segment(seq, PHIX_FRAGMENT + start, cont_len,
                            SEG_ADAPTER, 0, -1);
        TRACE("[contaminant] type=phiX, length=%d, fragment_offset=%d", cont_len, start);
    } else {
        /* Random nucleotides */
        char buf[512];
        if (cont_len > 511) cont_len = 511;
        rand_nucleotides(buf, cont_len);
        strncpy(rec->contaminant_type, "random", sizeof(rec->contaminant_type));
        aseq_append_segment(seq, buf, cont_len, SEG_ADAPTER, 0, -1);
        TRACE("[contaminant] type=random, length=%d", cont_len);
    }
}
