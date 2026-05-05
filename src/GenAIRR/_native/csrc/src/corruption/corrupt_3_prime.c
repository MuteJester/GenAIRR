/**
 * corrupt_3_prime.c — 3' end corruption (add/remove/remove-before-add).
 *
 * Mirror of corrupt_5_prime but operating on the 3' (tail) end.
 * No position shifting needed since 3' removal preserves indices.
 */

#include "genairr/pipeline.h"
#include "genairr/productivity_guard.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>

static void apply_corrupt_3_candidate(ASeq *seq, SimRecord *rec,
                                      int event,
                                      int remove_amount,
                                      int add_amount,
                                      const char *add_bases) {
    rec->corruption_3_event = event;

    if (event == 1 || event == 3) {
        if (remove_amount > 0) {
            aseq_delete_tail_n(seq, remove_amount);
        }
        rec->corruption_3_remove_amount = remove_amount;
    }

    if (event == 2 || event == 3) {
        if (add_amount > 0) {
            aseq_append_bases(seq, add_bases, add_amount, SEG_ADAPTER, 0);
        }
        rec->corruption_3_add_amount = add_amount;
    }
}

void step_corrupt_3_prime(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    int len = aseq_length(seq);
    if (len < 10) return;

    double r = rng_uniform(cfg->rng);
    int event;
    if (r < 0.4)      event = 1;
    else if (r < 0.7) event = 2;
    else               event = 3;

    TRACE("[corrupt_3'] event=%s, seq_len=%d",
          event == 1 ? "remove" : event == 2 ? "add" : "remove+add", len);

    int remove_amount = 0;
    int add_amount = 0;
    char add_buf[50];

    /* Remove from 3' end. Respects min=max=0 as "no-op" — see
     * corrupt_5_prime.c for the same fix rationale (T1-3). */
    if (event == 1 || event == 3) {
        int lo = cfg->corrupt_3_remove_min;
        int hi = cfg->corrupt_3_remove_max;
        remove_amount = lo + (int)(rng_uniform(cfg->rng) * (hi - lo + 1));
        if (remove_amount > hi) remove_amount = hi;
        if (remove_amount >= len - 5) remove_amount = len - 5;
        if (remove_amount < 0) remove_amount = 0;
    }

    /* Add random bases at 3' end. */
    if (event == 2 || event == 3) {
        len = aseq_length(seq);
        int lo = cfg->corrupt_3_add_min;
        int hi = cfg->corrupt_3_add_max;
        add_amount = lo + (int)(rng_uniform(cfg->rng) * (hi - lo + 1));
        if (add_amount > hi) add_amount = hi;
        if (add_amount > 50) add_amount = 50;
        if (add_amount < 0) add_amount = 0;

        if (add_amount > 0) {
            for (int i = 0; i < add_amount; i++) {
                static const char bases[] = "acgt";
                add_buf[i] = bases[rng_range(cfg->rng, 4)];
            }
        }
    }

    if (simcfg_productive_only(cfg)) {
        ASeq seq_copy = *seq;
        aseq_rebase(&seq_copy, seq);
        SimRecord rec_copy = *rec;

        apply_corrupt_3_candidate(
            &seq_copy, &rec_copy, event, remove_amount, add_amount, add_buf);

        ProductivityDecision decision = productivity_guard_state(
            cfg, &seq_copy, &rec_copy, PROD_STAGE_OBSERVED);
        if (decision != PROD_DECISION_ALLOW) {
            TRACE("[corrupt_3'] skipped productive-unsafe event (type=%d, remove=%d, add=%d)",
                  event, remove_amount, add_amount);
            return;
        }
    }

    apply_corrupt_3_candidate(seq, rec, event, remove_amount, add_amount, add_buf);

    if (event == 1 || event == 3) {
        TRACE("[corrupt_3'] removed %d bases from 3' end", remove_amount);
    }
    if (event == 2 || event == 3) {
        TRACE("[corrupt_3'] added %d random bases at 3' end", add_amount);
    }
}
