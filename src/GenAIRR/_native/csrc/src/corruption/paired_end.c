/**
 * paired_end.c — Paired-end merge quality profile.
 *
 * V-shaped (bathtub) error profile simulating merged paired-end reads:
 *   - R1 error increases from 5' toward center
 *   - R2 error increases from 3' toward center
 *   - Overlap zone: merged_err = R1_err × R2_err (very low)
 *   - Gap zone (seq > 2×read_length): positions become 'N'
 *
 * Uses same transition/transversion logic as quality_errors.
 */

#include "genairr/pipeline.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>

/* Transition partners */
static char pe_transition_of(char base) {
    switch (base) {
        case 'A': return 'G'; case 'G': return 'A';
        case 'C': return 'T'; case 'T': return 'C';
        default:  return base;
    }
}

static char pe_transversion_of(char base) {
    switch (base) {
        case 'A': return (rand() % 2) ? 'C' : 'T';
        case 'G': return (rand() % 2) ? 'C' : 'T';
        case 'C': return (rand() % 2) ? 'A' : 'G';
        case 'T': return (rand() % 2) ? 'A' : 'G';
        default:  return base;
    }
}

void step_paired_end(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    int len = aseq_length(seq);
    if (len == 0) return;

    int read_len = cfg->pe_read_length;
    rec->pe_read_length = read_len;

    TRACE("[paired_end] seq_len=%d, read_length=%d, overlap=%d, gap=%d",
          len, read_len,
          (len < 2 * read_len) ? 2 * read_len - len : 0,
          (len > 2 * read_len) ? len - 2 * read_len : 0);

    double base_err = 0.001;
    double peak_err = 0.02;

    int pos = 0;
    for (Nuc *n = seq->head; n; n = n->next, pos++) {
        /* Check if in gap zone */
        if (pos >= read_len && pos < len - read_len) {
            /* Gap zone: replace with N */
            aseq_mutate(seq, n, 'N', NUC_FLAG_IS_N);
            continue;
        }

        if (n->current == 'N') continue;

        /* R1 error rate: increases from 5' */
        double r1_err = (pos < read_len)
            ? base_err + ((double)pos / read_len) * (peak_err - base_err)
            : 1.0;  /* beyond R1, no R1 correction */

        /* R2 error rate: increases from 3' */
        int pos_from_3 = len - 1 - pos;
        double r2_err = (pos_from_3 < read_len)
            ? base_err + ((double)pos_from_3 / read_len) * (peak_err - base_err)
            : 1.0;

        /* Overlap zone: both reads cover this position */
        double eff_err;
        if (pos < read_len && pos_from_3 < read_len) {
            eff_err = r1_err * r2_err;  /* both must be wrong */
        } else if (pos < read_len) {
            eff_err = r1_err;
        } else {
            eff_err = r2_err;
        }

        if (rand_uniform() >= eff_err) continue;

        char new_base;
        if (rand_uniform() < 0.7) {
            new_base = pe_transition_of(n->current);
        } else {
            new_base = pe_transversion_of(n->current);
        }
        aseq_mutate(seq, n, new_base, NUC_FLAG_SEQ_ERROR);
    }

    rec->pe_gap_length = (len > 2 * read_len) ? len - 2 * read_len : 0;
    TRACE("[paired_end] gap_length=%d", rec->pe_gap_length);
}
