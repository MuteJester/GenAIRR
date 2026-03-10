/**
 * pipeline.h — Function-pointer pipeline builder and executor.
 *
 * The Pipeline is an array of step functions with a uniform signature.
 * Steps are conditionally registered at build time based on enabled
 * features. The pipeline is built once and executed N times.
 */

#ifndef GENAIRR_PIPELINE_H
#define GENAIRR_PIPELINE_H

#include "types.h"
#include "aseq.h"
#include "sim_config.h"

/* Forward-declare the simulation record */
typedef struct SimRecord SimRecord;

/* ── Step function signature ──────────────────────────────────── */

/**
 * Every pipeline step has the same signature:
 *   - cfg:    read-only simulation config (alleles, distributions, flags)
 *   - seq:    the annotated sequence being built/modified
 *   - rec:    the simulation record (scalar metadata alongside the sequence)
 */
typedef void (*StepFn)(const SimConfig *cfg, ASeq *seq, SimRecord *rec);

/* ── Pipeline ─────────────────────────────────────────────────── */

typedef struct {
    StepFn    steps[GENAIRR_MAX_PIPELINE];
    int       hook_ids[GENAIRR_MAX_PIPELINE];  /* HookPoint or -1 */
    int       n_steps;
    int       retry_boundary;    /* index of assess_functionality, or -1 */
} Pipeline;

/**
 * Build a pipeline from the SimConfig feature flags.
 *
 * This is the C equivalent of the Python 5-pass compiler.
 * The correction/injection logic is encoded directly in the
 * conditional construction order.
 */
Pipeline  pipeline_build(const SimConfig *cfg);

/* ── Simulation record ────────────────────────────────────────── */

/**
 * SimRecord holds scalar metadata that lives alongside the ASeq.
 * Things like the selected allele pointers, trim amounts, and
 * boolean flags that don't map to individual nucleotides.
 */
struct SimRecord {
    /* Selected alleles (pointers into SimConfig pools, not owned) */
    const Allele *v_allele;
    const Allele *d_allele;
    const Allele *j_allele;
    const Allele *c_allele;

    /* Trim amounts */
    int  v_trim_5;
    int  v_trim_3;
    int  d_trim_5;
    int  d_trim_3;
    int  j_trim_5;
    int  j_trim_3;
    int  c_trim_3;

    /* NP region lengths */
    int  np1_length;
    int  np2_length;

    /* Functionality assessment */
    bool productive;
    bool stop_codon;
    bool vj_in_frame;

    /* Sequence-level flags */
    bool is_reverse_complement;
    bool is_contaminant;
    bool receptor_revised;
    bool d_inverted;

    /* 5' corruption */
    int  corruption_5_event;      /* 0=none, 1=remove, 2=add, 3=remove+add */
    int  corruption_5_add_amount;
    int  corruption_5_remove_amount;

    /* 3' corruption */
    int  corruption_3_event;
    int  corruption_3_add_amount;
    int  corruption_3_remove_amount;

    /* UMI */
    char umi_sequence[32];
    int  umi_length;

    /* Primer mask */
    int  primer_masked_length;

    /* Paired-end */
    int  pe_read_length;
    int  pe_gap_length;

    /* PCR */
    int  pcr_n_cycles;

    /* Contaminant */
    char contaminant_type[16];    /* "random" or "phix" */

    /* Receptor revision */
    int  revision_footprint_length;
    char original_v_allele_name[GENAIRR_MAX_ALLELE_NAME];

    /* Note/status (small fixed buffer) */
    char note[256];
};

/** Reset a SimRecord to zero state. */
void  sim_record_init(SimRecord *rec);

/* ── Snapshot (for introspection hooks) ──────────────────────── */

#define MAX_SNAPSHOTS 16

typedef struct {
    ASeq      seq;           /* copied + rebased ASeq */
    SimRecord rec;           /* SimRecord at this point */
    int       hook_id;       /* which HookPoint */
    bool      used;
} Snapshot;

/**
 * Execute the pipeline once, producing one simulated sequence.
 *
 * If retry_boundary >= 0 (productive mode), retries the rearrangement
 * phase until productive or max attempts exhausted.
 *
 * Hook parameters (all optional — pass 0/NULL/NULL to disable):
 *   hook_mask: bitmask of HookPoint values to capture
 *   snaps:     array of MAX_SNAPSHOTS Snapshot slots
 *   n_snap:    pointer to snapshot count (incremented per capture)
 */
void  pipeline_execute(const Pipeline *p, const SimConfig *cfg,
                       ASeq *seq, SimRecord *rec,
                       uint32_t hook_mask, Snapshot *snaps, int *n_snap);

/* ── Snapshot helpers ────────────────────────────────────────── */

/** Take a snapshot of the current ASeq + SimRecord state. */
void take_snapshot(Snapshot *snaps, int *n, int hook_id,
                   const ASeq *seq, const SimRecord *rec);

/** Rebase all internal pointers in a copied ASeq to the new pool. */
void aseq_rebase(ASeq *dst, const ASeq *src);

/** Hook point name lookup. */
const char *hook_point_name(int hook_id);

/* ── CSR rate adjustment (simulate_csr.c) ──────────────────── */

/**
 * Adjust mutation rates based on the C allele isotype.
 * If rec->c_allele is NULL or unknown, rates are left unchanged.
 */
void  csr_adjust_rates(const SimRecord *rec,
                       double *min_rate, double *max_rate);

#endif /* GENAIRR_PIPELINE_H */
