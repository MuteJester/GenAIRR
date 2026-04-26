/**
 * corrupt_3_prime.c — 3' end corruption (add/remove/remove-before-add).
 *
 * Mirror of corrupt_5_prime but operating on the 3' (tail) end.
 * No position shifting needed since 3' removal preserves indices.
 */

#include "genairr/pipeline.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>

void step_corrupt_3_prime(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    int len = aseq_length(seq);
    if (len < 10) return;

    double r = rng_uniform(cfg->rng);
    int event;
    if (r < 0.4)      event = 1;
    else if (r < 0.7) event = 2;
    else               event = 3;

    rec->corruption_3_event = event;
    TRACE("[corrupt_3'] event=%s, seq_len=%d",
          event == 1 ? "remove" : event == 2 ? "add" : "remove+add", len);

    /* Remove from 3' end */
    if (event == 1 || event == 3) {
        int lo = cfg->corrupt_3_remove_min;
        int hi = cfg->corrupt_3_remove_max;
        int amount = lo + (int)(rng_uniform(cfg->rng) * (hi - lo + 1));
        if (amount > hi) amount = hi;
        if (amount >= len - 5) amount = len - 5;
        if (amount < 1) amount = 1;

        aseq_delete_tail_n(seq, amount);
        rec->corruption_3_remove_amount = amount;
        TRACE("[corrupt_3'] removed %d bases from 3' end", amount);
    }

    /* Add random bases at 3' end */
    if (event == 2 || event == 3) {
        len = aseq_length(seq);
        int lo = cfg->corrupt_3_add_min;
        int hi = cfg->corrupt_3_add_max;
        int amount = lo + (int)(rng_uniform(cfg->rng) * (hi - lo + 1));
        if (amount > hi) amount = hi;
        if (amount > 50) amount = 50;
        if (amount < 1) amount = 1;

        char buf[50];
        for (int i = 0; i < amount; i++) {
            static const char bases[] = "ACGT";
            buf[i] = bases[rng_range(cfg->rng, 4)];
        }
        aseq_append_bases(seq, buf, amount, SEG_ADAPTER, 0);
        rec->corruption_3_add_amount = amount;
        TRACE("[corrupt_3'] added %d random bases at 3' end", amount);
    }
}
