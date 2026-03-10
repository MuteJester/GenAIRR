/**
 * corrupt_5_prime.c — 5' end corruption (add/remove/remove-before-add).
 *
 * Simulates 5' end artifacts: degradation (remove), primer overshoot
 * (add), or adapter read-through (remove + add). Amounts are drawn
 * uniformly from configurable [min, max] ranges.
 */

#include "genairr/pipeline.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>
#include <string.h>

void step_corrupt_5_prime(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    int len = aseq_length(seq);
    if (len < 10) return;

    /* Pick event type: 40% remove, 30% add, 30% remove+add */
    double r = rand_uniform();
    int event;
    if (r < 0.4)      event = 1;  /* remove */
    else if (r < 0.7) event = 2;  /* add */
    else               event = 3;  /* remove + add */

    rec->corruption_5_event = event;
    TRACE("[corrupt_5'] event=%s, seq_len=%d",
          event == 1 ? "remove" : event == 2 ? "add" : "remove+add", len);

    /* Remove from 5' end */
    if (event == 1 || event == 3) {
        int lo = cfg->corrupt_5_remove_min;
        int hi = cfg->corrupt_5_remove_max;
        int amount = lo + (int)(rand_uniform() * (hi - lo + 1));
        if (amount > hi) amount = hi;
        if (amount >= len - 5) amount = len - 5;
        if (amount < 1) amount = 1;

        aseq_delete_head_n(seq, amount);
        rec->corruption_5_remove_amount = amount;
        TRACE("[corrupt_5'] removed %d bases from 5' end", amount);
    }

    /* Add random bases at 5' end */
    if (event == 2 || event == 3) {
        len = aseq_length(seq);
        int lo = cfg->corrupt_5_add_min;
        int hi = cfg->corrupt_5_add_max;
        int amount = lo + (int)(rand_uniform() * (hi - lo + 1));
        if (amount > hi) amount = hi;
        if (amount > 50) amount = 50;
        if (amount < 1) amount = 1;

        char buf[50];
        for (int i = 0; i < amount; i++) {
            static const char bases[] = "ACGT";
            buf[i] = bases[rand() % 4];
        }
        aseq_prepend_bases(seq, buf, amount, SEG_ADAPTER, 0);
        rec->corruption_5_add_amount = amount;
        TRACE("[corrupt_5'] added %d random bases at 5' end", amount);
    }
}
