/**
 * pipeline.c — Pipeline builder and executor.
 *
 * Builds the simulation pipeline by conditionally registering steps
 * based on enabled features. The correction logic is encoded in
 * the step order.
 *
 * Step registration order:
 *   1. Core rearrangement (sample → trim → assemble → assess)
 *   2. Simulation ops (d_inversion, receptor_revision)
 *   3. Mutation (S5F or Uniform — handled externally via s5f_mutate)
 *   4. Selection pressure (post-mutation)
 *   5. Corruption ops (in wet-lab order):
 *      primer_mask → umi → pcr → quality/paired_end/long_read →
 *      reverse_complement → contaminants → corrupt_5'/3' →
 *      indels → insert_ns → trim_to_length
 */

#include "genairr/pipeline.h"
#include "genairr/trace.h"
#include <string.h>
#include <stdio.h>

/* ── Forward declarations of all step functions ───────────────── */

/* Rearrangement primitives (src/rearrangement/) */
extern void step_sample_v(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_sample_d(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_sample_j(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_sample_c(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_trim_v(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_trim_d(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_trim_j(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_assemble(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_assess_functionality(const SimConfig *cfg, ASeq *seq, SimRecord *rec);

/* Simulation ops (src/simulation/) */
extern void step_d_inversion(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_receptor_revision(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_selection_pressure(const SimConfig *cfg, ASeq *seq, SimRecord *rec);

/* Mutation (src/mutation/) */
extern void step_uniform_mutate(const SimConfig *cfg, ASeq *seq, SimRecord *rec);

/* Corruption ops (src/corruption/) */
extern void step_primer_mask(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_simulate_umi(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_pcr_amplification(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_quality_errors(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_paired_end(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_long_read_errors(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_reverse_complement(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_spike_contaminants(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_corrupt_5_prime(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_corrupt_3_prime(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_insert_indels(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_insert_ns(const SimConfig *cfg, ASeq *seq, SimRecord *rec);
extern void step_trim_to_length(const SimConfig *cfg, ASeq *seq, SimRecord *rec);


/* ── Hook point names ─────────────────────────────────────────── */

static const char *hook_names[HOOK_COUNT] = {
    "post_assembly",       /* HOOK_POST_ASSEMBLY     */
    "post_functionality",  /* HOOK_POST_FUNCTIONALITY */
    "post_d_inversion",    /* HOOK_POST_D_INVERSION  */
    "post_receptor_rev",   /* HOOK_POST_RECEPTOR_REV */
    "post_mutation",       /* HOOK_POST_MUTATION     */
    "post_selection",      /* HOOK_POST_SELECTION    */
    "post_corrupt_5",      /* HOOK_POST_CORRUPT_5    */
    "post_corrupt_3",      /* HOOK_POST_CORRUPT_3    */
    "post_indels",         /* HOOK_POST_INDELS       */
    "post_ns",             /* HOOK_POST_NS           */
    "post_pcr",            /* HOOK_POST_PCR          */
    "post_quality",        /* HOOK_POST_QUALITY      */
    "final",               /* HOOK_FINAL             */
};

const char *hook_point_name(int hook_id) {
    if (hook_id < 0 || hook_id >= HOOK_COUNT) return "unknown";
    return hook_names[hook_id];
}

/* ── Helper: register a step with optional hook ──────────────── */

static inline void add_step(Pipeline *p, int *i, StepFn fn, int hook_id) {
    p->steps[*i] = fn;
    p->hook_ids[*i] = hook_id;
    (*i)++;
}

/* ── Pipeline builder ─────────────────────────────────────────── */

Pipeline pipeline_build(const SimConfig *cfg) {
    Pipeline p;
    memset(&p, 0, sizeof(p));
    p.retry_boundary = -1;
    p.post_mutation_start = 0;
    /* Initialize all hook_ids to -1 (no hook) */
    for (int j = 0; j < GENAIRR_MAX_PIPELINE; j++)
        p.hook_ids[j] = -1;

    int i = 0;
    bool has_d = chain_has_d(cfg->chain_type);

    /* ── Phase 1: Core rearrangement ─────────────────────────── */

    add_step(&p, &i, step_sample_v, -1);
    if (has_d) add_step(&p, &i, step_sample_d, -1);
    add_step(&p, &i, step_sample_j, -1);
    if (cfg->c_alleles.count > 0) add_step(&p, &i, step_sample_c, -1);

    add_step(&p, &i, step_trim_v, -1);
    if (has_d) add_step(&p, &i, step_trim_d, -1);
    add_step(&p, &i, step_trim_j, -1);

    add_step(&p, &i, step_assemble, HOOK_POST_ASSEMBLY);
    add_step(&p, &i, step_assess_functionality, HOOK_POST_FUNCTIONALITY);

    /* PRODUCTIVITY_MIXED disables retry (single rearrangement attempt).
     * PRODUCTIVE_ONLY / NON_PRODUCTIVE_ONLY both retry, with the
     * direction selected at execution time in pipeline_execute_pre_mutation. */
    if (cfg->features.productivity != PRODUCTIVITY_MIXED) {
        p.retry_boundary = i - 1;
    }

    /* ── Phase 2: Pre-mutation simulation ops ────────────────── */

    if (cfg->features.d_inversion && has_d) {
        add_step(&p, &i, step_d_inversion, HOOK_POST_D_INVERSION);
    }

    if (cfg->features.receptor_revision) {
        add_step(&p, &i, step_receptor_revision, HOOK_POST_RECEPTOR_REV);
    }

    /* ── Phase 3: Mutation ───────────────────────────────────── */
    /* S5F mutation is handled externally by the API, which executes
     * pre-mutation steps, applies S5F, then executes post-mutation
     * steps. For uniform mutation, we register the step directly.
     * Mutate feature flag is checked externally — uniform_mutate
     * uses the same min/max rates. */

    /* ── Phase 4: Post-mutation simulation ───────────────────── */
    p.post_mutation_start = i;

    if (cfg->features.selection_pressure) {
        add_step(&p, &i, step_selection_pressure, HOOK_POST_SELECTION);
    }

    /* ── Phase 5: Corruption ops (wet-lab order) ─────────────── */

    if (cfg->features.primer_mask) {
        add_step(&p, &i, step_primer_mask, -1);
    }

    if (cfg->features.umi) {
        add_step(&p, &i, step_simulate_umi, -1);
    }

    if (cfg->features.pcr) {
        add_step(&p, &i, step_pcr_amplification, HOOK_POST_PCR);
    }

    /* Quality/paired-end/long-read are mutually exclusive */
    if (cfg->features.quality_errors) {
        add_step(&p, &i, step_quality_errors, HOOK_POST_QUALITY);
    } else if (cfg->features.paired_end) {
        add_step(&p, &i, step_paired_end, -1);
    } else if (cfg->features.long_read_errors) {
        add_step(&p, &i, step_long_read_errors, -1);
    }

    if (cfg->features.reverse_complement) {
        add_step(&p, &i, step_reverse_complement, -1);
    }

    if (cfg->features.contaminants) {
        add_step(&p, &i, step_spike_contaminants, -1);
    }

    if (cfg->features.corrupt_5_prime) {
        add_step(&p, &i, step_corrupt_5_prime, HOOK_POST_CORRUPT_5);
    }

    if (cfg->features.corrupt_3_prime) {
        add_step(&p, &i, step_corrupt_3_prime, HOOK_POST_CORRUPT_3);
    }

    if (cfg->features.indels) {
        add_step(&p, &i, step_insert_indels, HOOK_POST_INDELS);
    }

    if (cfg->features.insert_ns) {
        add_step(&p, &i, step_insert_ns, HOOK_POST_NS);
    }

    /* Trim to length is always last */
    if (cfg->max_sequence_length > 0) {
        add_step(&p, &i, step_trim_to_length, -1);
    }

    p.n_steps = i;
    return p;
}

/* ── Snapshot helpers ─────────────────────────────────────────── */

void aseq_rebase(ASeq *dst, const ASeq *src) {
    /* After memcpy, all pointers in dst still reference src's pool.
     * Compute the offset and adjust every pointer. */
    ptrdiff_t offset = (char *)dst->pool - (char *)src->pool;

    #define REBASE(ptr) do { if (ptr) ptr = (Nuc *)((char *)(ptr) + offset); } while(0)

    REBASE(dst->head);
    REBASE(dst->tail);
    REBASE(dst->codon_head);
    REBASE(dst->v_anchor_node);
    REBASE(dst->j_anchor_node);

    for (int s = 0; s < SEG_COUNT; s++) {
        REBASE(dst->seg_first[s]);
        REBASE(dst->seg_last[s]);
    }

    /* Rebase per-node pointers */
    for (int n = 0; n < dst->pool_used; n++) {
        REBASE(dst->pool[n].next);
        REBASE(dst->pool[n].prev);
        REBASE(dst->pool[n].codon_next);
    }

    #undef REBASE
}

void take_snapshot(Snapshot *snaps, int *n, int hook_id,
                   const ASeq *seq, const SimRecord *rec) {
    if (*n >= MAX_SNAPSHOTS) return;
    Snapshot *s = &snaps[*n];
    memcpy(&s->seq, seq, sizeof(ASeq));
    aseq_rebase(&s->seq, seq);
    memcpy(&s->rec, rec, sizeof(SimRecord));
    s->hook_id = hook_id;
    s->used = true;
    (*n)++;
}

/* ── Pipeline executor ────────────────────────────────────────── */

static inline void check_hook(const Pipeline *p, int step_idx,
                               const ASeq *seq, const SimRecord *rec,
                               uint32_t hook_mask, Snapshot *snaps, int *n_snap) {
    if (!snaps) return;
    int hid = p->hook_ids[step_idx];
    if (hid >= 0 && (hook_mask & (1u << hid))) {
        take_snapshot(snaps, n_snap, hid, seq, rec);
    }
}

static inline void run_step_range(const Pipeline *p, const SimConfig *cfg,
                                  ASeq *seq, SimRecord *rec,
                                  int start, int end,
                                  uint32_t hook_mask,
                                  Snapshot *snaps, int *n_snap) {
    for (int i = start; i < end; i++) {
        p->steps[i](cfg, seq, rec);
        check_hook(p, i, seq, rec, hook_mask, snaps, n_snap);
    }
}

void pipeline_execute_pre_mutation(const Pipeline *p, const SimConfig *cfg,
                                   ASeq *seq, SimRecord *rec,
                                   uint32_t hook_mask, Snapshot *snaps, int *n_snap) {
    int pre_end = p->post_mutation_start;
    if (pre_end < 0) pre_end = 0;
    if (pre_end > p->n_steps) pre_end = p->n_steps;

    aseq_reset(seq);
    sim_record_init(rec);

    if (p->retry_boundary < 0 || pre_end == 0) {
        /* Non-productive: straight run */
        TRACE("[pipeline] mode=non-productive, pre-mutation steps=%d", pre_end);
        run_step_range(p, cfg, seq, rec, 0, pre_end, hook_mask, snaps, n_snap);
        return;
    }

    /* Filtered mode: retry rearrangement phase until the configured
     * direction matches. PRODUCTIVE_ONLY succeeds on rec->productive,
     * NON_PRODUCTIVE_ONLY succeeds on !rec->productive. */
    int max_attempts = cfg->max_productive_attempts;
    if (max_attempts <= 0) max_attempts = 25;

    bool want_productive =
        (cfg->features.productivity == PRODUCTIVITY_PRODUCTIVE_ONLY);

    TRACE("[pipeline] mode=%s, max_attempts=%d, pre-mutation steps=%d",
          want_productive ? "productive_only" : "non_productive_only",
          max_attempts, pre_end);

    int attempt;
    for (attempt = 0; attempt < max_attempts; attempt++) {
        aseq_reset(seq);
        sim_record_init(rec);

        TRACE("[pipeline] rearrangement attempt %d/%d", attempt + 1, max_attempts);

        for (int i = 0; i <= p->retry_boundary; i++) {
            p->steps[i](cfg, seq, rec);
            /* Only check hooks on final successful attempt */
        }

        bool target_hit = want_productive ? rec->productive : !rec->productive;
        if (target_hit) {
            TRACE("[pipeline] target hit on attempt %d (productive=%s)",
                  attempt + 1, rec->productive ? "yes" : "no");
            break;
        } else {
            TRACE("[pipeline] target miss (productive=%s, reason: %s), retrying...",
                  rec->productive ? "yes" : "no",
                  rec->note[0] ? rec->note : "ok");
        }
    }

    /* Check hooks for rearrangement steps on the final attempt */
    if (snaps && hook_mask) {
        for (int i = 0; i <= p->retry_boundary; i++) {
            check_hook(p, i, seq, rec, hook_mask, snaps, n_snap);
        }
    }

    /* Run remaining pre-mutation steps once */
    run_step_range(p, cfg, seq, rec,
                   p->retry_boundary + 1, pre_end,
                   hook_mask, snaps, n_snap);
}

void pipeline_execute_post_mutation(const Pipeline *p, const SimConfig *cfg,
                                    ASeq *seq, SimRecord *rec,
                                    uint32_t hook_mask, Snapshot *snaps, int *n_snap) {
    int start = p->post_mutation_start;
    if (start < 0) start = 0;
    if (start > p->n_steps) start = p->n_steps;

    TRACE("[pipeline] post-mutation steps=%d", p->n_steps - start);
    run_step_range(p, cfg, seq, rec, start, p->n_steps, hook_mask, snaps, n_snap);
}

void pipeline_execute(const Pipeline *p, const SimConfig *cfg,
                      ASeq *seq, SimRecord *rec,
                      uint32_t hook_mask, Snapshot *snaps, int *n_snap) {
    pipeline_execute_pre_mutation(p, cfg, seq, rec, hook_mask, snaps, n_snap);
    pipeline_execute_post_mutation(p, cfg, seq, rec, hook_mask, snaps, n_snap);
}

/* ── SimRecord ────────────────────────────────────────────────── */

void sim_record_init(SimRecord *rec) {
    memset(rec, 0, sizeof(*rec));
}
