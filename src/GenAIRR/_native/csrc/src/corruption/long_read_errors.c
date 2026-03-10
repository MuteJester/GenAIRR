/**
 * long_read_errors.c — Nanopore/PacBio homopolymer-targeted indels.
 *
 * Scans for homopolymer runs (≥min_run_length). Error probability
 * scales with run length: min(0.8, error_rate × run_length).
 * Insertions add one copy of the run base; deletions remove the last.
 *
 * Processes right-to-left so position changes don't affect earlier runs.
 */

#include "genairr/pipeline.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>

void step_long_read_errors(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)rec;
    int len = aseq_length(seq);
    if (len < 5) return;

    double error_rate = cfg->long_read_error_rate;
    int min_run = cfg->min_run_length;
    double ins_bias = cfg->insertion_bias;

    /* Build array of node pointers for right-to-left traversal */
    /* We identify homopolymer runs by walking forward, then
     * process them in reverse order. */

    /* Walk forward, find homopolymer runs */
    typedef struct { Nuc *start; Nuc *end; int length; char base; } Run;
    Run runs[256];
    int n_runs = 0;

    Nuc *run_start = seq->head;
    int run_len = 1;
    char run_base = run_start ? run_start->current : 0;

    for (Nuc *n = seq->head ? seq->head->next : NULL; n; n = n->next) {
        if (n->current == run_base) {
            run_len++;
        } else {
            if (run_len >= min_run && n_runs < 256) {
                runs[n_runs].start = run_start;
                runs[n_runs].end = n->prev;
                runs[n_runs].length = run_len;
                runs[n_runs].base = run_base;
                n_runs++;
            }
            run_start = n;
            run_len = 1;
            run_base = n->current;
        }
    }
    /* Final run */
    if (run_len >= min_run && n_runs < 256) {
        runs[n_runs].start = run_start;
        runs[n_runs].end = seq->tail;
        runs[n_runs].length = run_len;
        runs[n_runs].base = run_base;
        n_runs++;
    }

    TRACE("[long_read] found %d homopolymer runs (min_run=%d, error_rate=%.4f, ins_bias=%.2f)",
          n_runs, min_run, error_rate, ins_bias);

    /* Process runs right-to-left */
    int insertions = 0, deletions_applied = 0;
    for (int i = n_runs - 1; i >= 0; i--) {
        double p = error_rate * runs[i].length;
        if (p > 0.8) p = 0.8;
        if (rand_uniform() >= p) continue;

        if (rand_uniform() < ins_bias) {
            /* Insertion: add one copy at end of run */
            aseq_insert_after(seq, runs[i].end, runs[i].base,
                              runs[i].end->segment, NUC_FLAG_INDEL_INS);
            insertions++;
        } else {
            /* Deletion: remove last base of run */
            if (runs[i].length > 1 && !nuc_is_anchor(runs[i].end)) {
                aseq_delete(seq, runs[i].end);
                deletions_applied++;
            }
        }
    }

    TRACE("[long_read] applied %d homopolymer insertions, %d deletions", insertions, deletions_applied);
}
