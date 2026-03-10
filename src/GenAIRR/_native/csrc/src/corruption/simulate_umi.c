/**
 * simulate_umi.c — Prepend a random UMI barcode.
 *
 * Generates a random ACGT barcode and prepends it to the sequence.
 * Uses batch prepend for efficiency (single codon rail rebuild).
 */

#include "genairr/pipeline.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>
#include <string.h>

void step_simulate_umi(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    int umi_len = cfg->umi_length;
    if (umi_len <= 0 || umi_len > 30) return;

    /* Generate random UMI and store in record */
    char umi[32];
    rand_nucleotides(umi, umi_len);
    strncpy(rec->umi_sequence, umi, sizeof(rec->umi_sequence) - 1);
    rec->umi_length = umi_len;

    /* Batch prepend UMI nodes before head, tagged as SEG_UMI */
    aseq_prepend_bases(seq, umi, umi_len, SEG_UMI, 0);

    TRACE("[umi] prepended %dbp UMI: %s", umi_len, umi);
}
