/**
 * genairr_api.c — Public C API implementation.
 *
 * Wraps the internal SimConfig/Pipeline/ASeq/S5F/AIRR machinery
 * behind a clean opaque-handle interface for Python bindings.
 */

#include "genairr/genairr_api.h"
#include "genairr/genairr.h"
#include "genairr/gdc_io.h"
#include "genairr/sim_config.h"
#include "genairr/pipeline.h"
#include "genairr/airr.h"
#include "genairr/s5f.h"
#include "genairr/allele_bitmap.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
  #include <io.h>
  #include <process.h>
  #define GENAIRR_UNLINK  _unlink
  #define GENAIRR_CLOSE   _close
  #define GENAIRR_FDOPEN  _fdopen
#else
  #include <unistd.h>
  #define GENAIRR_UNLINK  unlink
  #define GENAIRR_CLOSE   close
  #define GENAIRR_FDOPEN  fdopen
#endif

/* ── Opaque handle definition ────────────────────────────────── */

struct GenAIRRSimulator {
    GdcData       gdc;
    SimConfig     cfg;
    S5FModel      s5f;
    bool          s5f_loaded;     /* true if GDC had S5F data */
    bool          compiled;       /* true after first simulate() */
    bool          seeded;         /* true after RNG seeded (for streaming) */
    Pipeline      pipeline;
    AlleleCorrectionSet corr;
    uint64_t      seed;
    RngState      rng_state;      /* per-simulator PCG32; pointed to by cfg.rng */
    /* Trace */
    bool          trace_enabled;
    TraceLog      trace;
    bool          trace_initialized;
    /* Introspection hooks */
    uint32_t      hook_mask;
    Snapshot      snapshots[MAX_SNAPSHOTS];
    int           n_snapshots;
};

/* ── Helpers ─────────────────────────────────────────────────── */

static int init_from_gdc(GenAIRRSimulator *sim) {
    /* Map GDC chain_type to internal ChainType enum */
    ChainType ct;
    switch (sim->gdc.chain_type) {
        case 0: ct = CHAIN_IGH;  break;
        case 1: ct = CHAIN_IGK;  break;
        case 2: ct = CHAIN_IGL;  break;
        case 3: ct = CHAIN_TCRA; break;
        case 4: ct = CHAIN_TCRB; break;
        case 5: ct = CHAIN_TCRD; break;
        case 6: ct = CHAIN_TCRG; break;
        default: return -1;
    }

    sim_config_init(&sim->cfg, ct);
    /* Wire the simulator's owned RngState into the SimConfig so step
     * functions can draw via cfg->rng. ensure_seeded() seeds it on
     * the first simulate() call. */
    sim->cfg.rng = &sim->rng_state;
    gdc_populate_sim_config(&sim->cfg, &sim->gdc);

    /* Load S5F if present */
    if (sim->gdc.mutation_model_type == 1) { /* GDC_MUTMODEL_S5F */
        s5f_model_init(&sim->s5f, 0.01, 0.05, false);
        gdc_populate_s5f_model(&sim->s5f, &sim->gdc);
        sim->s5f_loaded = true;
    }

    return 0;
}

static void compile(GenAIRRSimulator *sim) {
    if (sim->compiled) return;

    /* Sync S5F params with SimConfig mutation rates */
    if (sim->s5f_loaded && sim->cfg.features.mutate) {
        sim->s5f.min_mutation_rate = sim->cfg.min_mutation_rate;
        sim->s5f.max_mutation_rate = sim->cfg.max_mutation_rate;
        sim->s5f.productive = sim->cfg.features.productive;
    }

    sim->pipeline = pipeline_build(&sim->cfg);
    sim->corr = allele_correction_set_build(&sim->cfg);
    sim->compiled = true;
}

/* ── Lifecycle ───────────────────────────────────────────────── */

GenAIRRSimulator *genairr_create(const char *gdc_path) {
    if (!gdc_path) return NULL;

    GenAIRRSimulator *sim = calloc(1, sizeof(*sim));
    if (!sim) return NULL;

    gdc_data_init(&sim->gdc);

    int rc = gdc_load(gdc_path, &sim->gdc);
    if (rc != 0) {
        free(sim);
        return NULL;
    }

    if (init_from_gdc(sim) != 0) {
        gdc_data_destroy(&sim->gdc);
        free(sim);
        return NULL;
    }

    return sim;
}

GenAIRRSimulator *genairr_create_from_memory(const void *gdc_bytes,
                                              int gdc_len) {
    if (!gdc_bytes || gdc_len <= 0) return NULL;

    /* Write to temp file, load via gdc_load, remove temp. */
#ifdef _WIN32
    /* Windows: use tmpnam_s for a temp path */
    char tmppath[L_tmpnam_s + 1];
    if (tmpnam_s(tmppath, sizeof(tmppath)) != 0) return NULL;
    FILE *fp = fopen(tmppath, "wb");
#else
    char tmppath[] = "/tmp/genairr_gdc_XXXXXX";
    int fd = mkstemp(tmppath);
    if (fd < 0) return NULL;
    FILE *fp = GENAIRR_FDOPEN(fd, "wb");
#endif
    if (!fp) {
#ifndef _WIN32
        GENAIRR_CLOSE(fd);
#endif
        GENAIRR_UNLINK(tmppath);
        return NULL;
    }
    fwrite(gdc_bytes, 1, gdc_len, fp);
    fclose(fp);

    GenAIRRSimulator *sim = genairr_create(tmppath);
    GENAIRR_UNLINK(tmppath);
    return sim;
}

void genairr_destroy(GenAIRRSimulator *sim) {
    if (!sim) return;
    if (sim->compiled) {
        allele_correction_set_destroy(&sim->corr);
    }
    if (sim->trace_initialized) {
        trace_log_destroy(&sim->trace);
    }
    sim_config_destroy(&sim->cfg);
    gdc_data_destroy(&sim->gdc);
    free(sim);
}

/* ── Feature flags ───────────────────────────────────────────── */

int genairr_set_feature(GenAIRRSimulator *sim,
                         const char *name, int enabled) {
    if (!sim || !name) return -1;
    sim->compiled = false;  /* invalidate compiled state */

    SimFeatures *f = &sim->cfg.features;
    bool val = (enabled != 0);

    if      (strcmp(name, "productive") == 0)         f->productive = val;
    else if (strcmp(name, "mutate") == 0)             f->mutate = val;
    else if (strcmp(name, "selection_pressure") == 0) f->selection_pressure = val;
    else if (strcmp(name, "csr") == 0)                f->csr = val;
    else if (strcmp(name, "d_inversion") == 0)        f->d_inversion = val;
    else if (strcmp(name, "receptor_revision") == 0)  f->receptor_revision = val;
    else if (strcmp(name, "corrupt_5_prime") == 0)    f->corrupt_5_prime = val;
    else if (strcmp(name, "corrupt_3_prime") == 0)    f->corrupt_3_prime = val;
    else if (strcmp(name, "quality_errors") == 0)     f->quality_errors = val;
    else if (strcmp(name, "paired_end") == 0)         f->paired_end = val;
    else if (strcmp(name, "pcr") == 0)                f->pcr = val;
    else if (strcmp(name, "umi") == 0)                f->umi = val;
    else if (strcmp(name, "primer_mask") == 0)        f->primer_mask = val;
    else if (strcmp(name, "reverse_complement") == 0) f->reverse_complement = val;
    else if (strcmp(name, "contaminants") == 0)       f->contaminants = val;
    else if (strcmp(name, "long_read_errors") == 0)   f->long_read_errors = val;
    else if (strcmp(name, "indels") == 0)             f->indels = val;
    else if (strcmp(name, "insert_ns") == 0)          f->insert_ns = val;
    else return -1;

    return 0;
}

/* ── Parameters ──────────────────────────────────────────────── */

int genairr_set_param_f64(GenAIRRSimulator *sim,
                           const char *name, double value) {
    if (!sim || !name) return -1;
    sim->compiled = false;

    SimConfig *c = &sim->cfg;

    if      (strcmp(name, "min_mutation_rate") == 0)    c->min_mutation_rate = value;
    else if (strcmp(name, "max_mutation_rate") == 0)    c->max_mutation_rate = value;
    else if (strcmp(name, "selection_strength") == 0)   c->selection_strength = value;
    else if (strcmp(name, "cdr_r_acceptance") == 0)     c->cdr_r_acceptance = value;
    else if (strcmp(name, "fwr_r_acceptance") == 0)     c->fwr_r_acceptance = value;
    else if (strcmp(name, "d_inversion_prob") == 0)     c->d_inversion_prob = value;
    else if (strcmp(name, "revision_prob") == 0)        c->revision_prob = value;
    else if (strcmp(name, "base_error_rate") == 0)      c->base_error_rate = value;
    else if (strcmp(name, "peak_error_rate") == 0)      c->peak_error_rate = value;
    else if (strcmp(name, "transition_weight") == 0)    c->transition_weight = value;
    else if (strcmp(name, "pcr_error_rate") == 0)       c->pcr_error_rate = value;
    else if (strcmp(name, "contamination_prob") == 0)   c->contamination_prob = value;
    else if (strcmp(name, "indel_prob") == 0)           c->indel_prob = value;
    else if (strcmp(name, "n_prob") == 0)               c->n_prob = value;
    else if (strcmp(name, "rc_prob") == 0)              c->rc_prob = value;
    else if (strcmp(name, "long_read_error_rate") == 0) c->long_read_error_rate = value;
    else if (strcmp(name, "insertion_bias") == 0)       c->insertion_bias = value;
    else return -1;

    return 0;
}

int genairr_set_param_i32(GenAIRRSimulator *sim,
                           const char *name, int value) {
    if (!sim || !name) return -1;
    sim->compiled = false;

    SimConfig *c = &sim->cfg;

    if      (strcmp(name, "pcr_cycles") == 0)              c->pcr_n_cycles = value;
    else if (strcmp(name, "umi_length") == 0)              c->umi_length = value;
    else if (strcmp(name, "pe_read_length") == 0)          c->pe_read_length = value;
    else if (strcmp(name, "max_sequence_length") == 0)     c->max_sequence_length = value;
    else if (strcmp(name, "footprint_min") == 0)           c->footprint_min = value;
    else if (strcmp(name, "footprint_max") == 0)           c->footprint_max = value;
    else if (strcmp(name, "min_run_length") == 0)          c->min_run_length = value;
    else if (strcmp(name, "primer_mask_length") == 0)      c->primer_mask_length = value;
    else if (strcmp(name, "max_productive_attempts") == 0) c->max_productive_attempts = value;
    else if (strcmp(name, "contaminant_type") == 0)        c->contaminant_type = value;
    else if (strcmp(name, "corrupt_5_remove_min") == 0)   c->corrupt_5_remove_min = value;
    else if (strcmp(name, "corrupt_5_remove_max") == 0)   c->corrupt_5_remove_max = value;
    else if (strcmp(name, "corrupt_5_add_min") == 0)      c->corrupt_5_add_min = value;
    else if (strcmp(name, "corrupt_5_add_max") == 0)      c->corrupt_5_add_max = value;
    else if (strcmp(name, "corrupt_3_remove_min") == 0)   c->corrupt_3_remove_min = value;
    else if (strcmp(name, "corrupt_3_remove_max") == 0)   c->corrupt_3_remove_max = value;
    else if (strcmp(name, "corrupt_3_add_min") == 0)      c->corrupt_3_add_min = value;
    else if (strcmp(name, "corrupt_3_add_max") == 0)      c->corrupt_3_add_max = value;
    else return -1;

    return 0;
}

/* ── Allele locking ─────────────────────────────────────────── */

static int resolve_segment(const char *segment,
                            AllelePool **pool, AlleleRestriction **restr,
                            GenAIRRSimulator *sim) {
    if      (strcmp(segment, "v") == 0 || strcmp(segment, "V") == 0)
        { *pool = &sim->cfg.v_alleles; *restr = &sim->cfg.v_restriction; }
    else if (strcmp(segment, "d") == 0 || strcmp(segment, "D") == 0)
        { *pool = &sim->cfg.d_alleles; *restr = &sim->cfg.d_restriction; }
    else if (strcmp(segment, "j") == 0 || strcmp(segment, "J") == 0)
        { *pool = &sim->cfg.j_alleles; *restr = &sim->cfg.j_restriction; }
    else if (strcmp(segment, "c") == 0 || strcmp(segment, "C") == 0)
        { *pool = &sim->cfg.c_alleles; *restr = &sim->cfg.c_restriction; }
    else
        return -1;
    return 0;
}

int genairr_lock_allele(GenAIRRSimulator *sim,
                         const char *segment, const char *name) {
    if (!sim || !segment || !name) return -1;

    AllelePool *pool = NULL;
    AlleleRestriction *restr = NULL;
    if (resolve_segment(segment, &pool, &restr, sim) != 0) return -1;

    int idx = allele_pool_find_index(pool, name);
    if (idx < 0) return -1;

    if (restr->count >= GENAIRR_MAX_LOCKED) return -2;

    /* Avoid duplicates */
    for (int i = 0; i < restr->count; i++) {
        if (restr->indices[i] == idx) return 0;
    }

    restr->indices[restr->count++] = idx;
    restr->active = true;
    sim->compiled = false;
    return 0;
}

void genairr_clear_locks(GenAIRRSimulator *sim, const char *segment) {
    if (!sim) return;

    if (!segment) {
        /* Clear all */
        sim->cfg.v_restriction = (AlleleRestriction){0};
        sim->cfg.d_restriction = (AlleleRestriction){0};
        sim->cfg.j_restriction = (AlleleRestriction){0};
        sim->cfg.c_restriction = (AlleleRestriction){0};
    } else {
        AllelePool *pool = NULL;
        AlleleRestriction *restr = NULL;
        if (resolve_segment(segment, &pool, &restr, sim) == 0) {
            *restr = (AlleleRestriction){0};
        }
    }
    sim->compiled = false;
}

/* ── Seed ────────────────────────────────────────────────────── */

void genairr_set_seed(GenAIRRSimulator *sim, uint64_t seed) {
    if (!sim) return;
    sim->seed = seed;
    sim->seeded = false;  /* re-seed on next simulate_one() */
}

/* ── Simulate ────────────────────────────────────────────────── */

static void ensure_seeded(GenAIRRSimulator *sim);

int genairr_simulate(GenAIRRSimulator *sim, int n,
                      const char *output_path) {
    if (!sim || n <= 0 || !output_path) return -1;

    FILE *fp = fopen(output_path, "w");
    if (!fp) return -1;

    /* Compile if needed */
    compile(sim);

    /* Seed RNG (once only — avoid re-seeding on batched calls) */
    ensure_seeded(sim);

    /* Write header */
    airr_write_tsv_header(fp);

    /* Simulate */
    ASeq seq;
    aseq_init(&seq);

    SimRecord rec;
    AirrRecord airr;
    S5FResult s5f_result;

    bool use_s5f = sim->s5f_loaded && sim->cfg.features.mutate;
    bool use_csr = use_s5f && sim->cfg.features.csr;

    /* Save base mutation rates (CSR overrides per-sequence) */
    double base_min = sim->s5f.min_mutation_rate;
    double base_max = sim->s5f.max_mutation_rate;

    for (int i = 0; i < n; i++) {
        pipeline_execute(&sim->pipeline, &sim->cfg, &seq, &rec,
                         0, NULL, NULL);

        /* S5F mutation (runs after pipeline, before AIRR serialization) */
        if (use_s5f) {
            if (use_csr) {
                /* Adjust rates per isotype, then restore */
                double min_r = base_min, max_r = base_max;
                csr_adjust_rates(&rec, &min_r, &max_r);
                sim->s5f.min_mutation_rate = min_r;
                sim->s5f.max_mutation_rate = max_r;
            }
            s5f_mutate(&sim->s5f, &seq, &rec, &sim->rng_state, &s5f_result);
            if (use_csr) {
                sim->s5f.min_mutation_rate = base_min;
                sim->s5f.max_mutation_rate = base_max;
            }
        }

        airr_serialize(&seq, &rec, &sim->cfg, &sim->corr, &airr);
        airr_write_tsv_row(fp, &airr);
    }

    fclose(fp);
    return n;
}

int genairr_simulate_to_fd(GenAIRRSimulator *sim, int n, int fd) {
    if (!sim || n <= 0 || fd < 0) return -1;

    FILE *fp = GENAIRR_FDOPEN(fd, "w");
    if (!fp) return -1;

    compile(sim);

    /* Seed RNG (once only — avoid re-seeding on batched calls) */
    ensure_seeded(sim);

    airr_write_tsv_header(fp);

    ASeq seq;
    aseq_init(&seq);

    SimRecord rec;
    AirrRecord airr;
    S5FResult s5f_result;

    bool use_s5f = sim->s5f_loaded && sim->cfg.features.mutate;
    bool use_csr = use_s5f && sim->cfg.features.csr;

    double base_min = sim->s5f.min_mutation_rate;
    double base_max = sim->s5f.max_mutation_rate;

    for (int i = 0; i < n; i++) {
        pipeline_execute(&sim->pipeline, &sim->cfg, &seq, &rec,
                         0, NULL, NULL);

        if (use_s5f) {
            if (use_csr) {
                double min_r = base_min, max_r = base_max;
                csr_adjust_rates(&rec, &min_r, &max_r);
                sim->s5f.min_mutation_rate = min_r;
                sim->s5f.max_mutation_rate = max_r;
            }
            s5f_mutate(&sim->s5f, &seq, &rec, &sim->rng_state, &s5f_result);
            if (use_csr) {
                sim->s5f.min_mutation_rate = base_min;
                sim->s5f.max_mutation_rate = base_max;
            }
        }

        airr_serialize(&seq, &rec, &sim->cfg, &sim->corr, &airr);
        airr_write_tsv_row(fp, &airr);
    }

    /* Don't fclose — caller owns the fd */
    fflush(fp);
    return n;
}

/* ── Streaming (single-sequence) ─────────────────────────────── */

static void ensure_seeded(GenAIRRSimulator *sim) {
    if (!sim->seeded) {
        uint64_t seed;
        if (sim->seed != 0) {
            /* Explicit seed: deterministic, independent of simulator
             * identity. Two simulators with the same SimConfig and
             * same explicit seed produce byte-identical output. */
            seed = sim->seed;
        } else {
            /* Auto-seed: mix wall time with simulator address so two
             * concurrent simulators in the same process get distinct
             * streams even when both default. */
            seed = (uint64_t)time(NULL) ^ (uint64_t)(uintptr_t)sim;
        }
        rng_seed(&sim->rng_state, seed, 0);
        sim->seeded = true;
    }
}

int genairr_simulate_one(GenAIRRSimulator *sim,
                          char *buffer, int buffer_size) {
    if (!sim || !buffer || buffer_size <= 0) return -1;

    compile(sim);
    ensure_seeded(sim);

    /* Set up tracing if enabled */
    if (sim->trace_enabled) {
        if (!sim->trace_initialized) {
            trace_log_init(&sim->trace);
            sim->trace_initialized = true;
        }
        trace_log_clear(&sim->trace);
        trace_set_current(&sim->trace);
    }

    ASeq seq;
    aseq_init(&seq);

    SimRecord rec;
    AirrRecord airr;
    S5FResult s5f_result;

    bool use_s5f = sim->s5f_loaded && sim->cfg.features.mutate;
    bool use_csr = use_s5f && sim->cfg.features.csr;

    /* Reset snapshot state */
    sim->n_snapshots = 0;

    TRACE("[begin] simulate_one: chain=%d, s5f=%s, csr=%s, productive=%s",
          sim->cfg.chain_type,
          use_s5f ? "yes" : "no",
          use_csr ? "yes" : "no",
          sim->cfg.features.productive ? "yes" : "no");

    pipeline_execute(&sim->pipeline, &sim->cfg, &seq, &rec,
                     sim->hook_mask, sim->snapshots, &sim->n_snapshots);

    if (use_s5f) {
        double base_min = sim->s5f.min_mutation_rate;
        double base_max = sim->s5f.max_mutation_rate;
        if (use_csr) {
            double min_r = base_min, max_r = base_max;
            TRACE("[csr] isotype=%s, base_rates=[%.4f, %.4f]",
                  rec.c_allele ? rec.c_allele->name : "none", base_min, base_max);
            csr_adjust_rates(&rec, &min_r, &max_r);
            sim->s5f.min_mutation_rate = min_r;
            sim->s5f.max_mutation_rate = max_r;
            TRACE("[csr] adjusted_rates=[%.4f, %.4f]", min_r, max_r);
        }
        s5f_mutate(&sim->s5f, &seq, &rec, &sim->rng_state, &s5f_result);
        if (use_csr) {
            sim->s5f.min_mutation_rate = base_min;
            sim->s5f.max_mutation_rate = base_max;
        }
    }

    /* Post-mutation hook (S5F runs outside pipeline) */
    if (sim->hook_mask & (1u << HOOK_POST_MUTATION)) {
        take_snapshot(sim->snapshots, &sim->n_snapshots,
                      HOOK_POST_MUTATION, &seq, &rec);
    }

    TRACE("[serialize] begin AIRR serialization");
    airr_serialize(&seq, &rec, &sim->cfg, &sim->corr, &airr);
    TRACE("[serialize] complete: seq_len=%d, v_call=%s, d_call=%s, j_call=%s",
          airr.sequence_length, airr.v_call, airr.d_call, airr.j_call);

    /* Final hook (after serialization) */
    if (sim->hook_mask & (1u << HOOK_FINAL)) {
        take_snapshot(sim->snapshots, &sim->n_snapshots,
                      HOOK_FINAL, &seq, &rec);
    }

    TRACE("[end] simulate_one complete");

    /* Clear thread-local trace pointer */
    if (sim->trace_enabled) {
        trace_set_current(NULL);
    }

    return airr_snprintf_tsv_row(buffer, buffer_size, &airr);
}

int genairr_get_header(char *buffer, int buffer_size) {
    return airr_snprintf_tsv_header(buffer, buffer_size);
}

/* ── Tracing ──────────────────────────────────────────────────── */

void genairr_set_trace(GenAIRRSimulator *sim, int enabled) {
    if (!sim) return;
    sim->trace_enabled = (enabled != 0);
    if (sim->trace_enabled && !sim->trace_initialized) {
        trace_log_init(&sim->trace);
        sim->trace_initialized = true;
    }
}

int genairr_get_trace(GenAIRRSimulator *sim,
                       char *buffer, int buffer_size) {
    if (!sim || !buffer || buffer_size <= 0) return 0;
    if (!sim->trace_initialized || !sim->trace.buf) return 0;

    int len = sim->trace.len;
    int copy = len < buffer_size - 1 ? len : buffer_size - 1;
    memcpy(buffer, sim->trace.buf, copy);
    buffer[copy] = '\0';
    return len;
}

void genairr_clear_trace(GenAIRRSimulator *sim) {
    if (!sim || !sim->trace_initialized) return;
    trace_log_clear(&sim->trace);
}

/* ── Introspection hooks ─────────────────────────────────────── */

void genairr_set_hooks(GenAIRRSimulator *sim, uint32_t hook_mask) {
    if (!sim) return;
    sim->hook_mask = hook_mask;
}

int genairr_get_snapshot_count(GenAIRRSimulator *sim) {
    if (!sim) return 0;
    return sim->n_snapshots;
}

const char *genairr_get_snapshot_name(GenAIRRSimulator *sim, int index) {
    if (!sim || index < 0 || index >= sim->n_snapshots) return "";
    return hook_point_name(sim->snapshots[index].hook_id);
}

/* ── ASeq JSON serialization ─────────────────────────────────── */

static const char *seg_names[] = {
    "V", "NP1", "D", "NP2", "J", "C", "UMI", "ADAPTER"
};

static int aseq_to_json(const ASeq *seq, char *buf, int bufsize) {
    int pos = 0;
    int seq_pos = 0;

    #define EMIT(...) do { \
        int _w = snprintf(buf + pos, bufsize - pos, __VA_ARGS__); \
        if (_w < 0 || pos + _w >= bufsize) goto truncated; \
        pos += _w; \
    } while(0)

    EMIT("[");

    for (const Nuc *n = seq->head; n; n = n->next) {
        if (seq_pos > 0) EMIT(",");

        /* Build flags array */
        char flags_buf[128];
        int fpos = 0;
        flags_buf[fpos++] = '[';
        int nflags = 0;

        #define FLAG(mask, name) do { \
            if (n->flags & (mask)) { \
                if (nflags > 0) flags_buf[fpos++] = ','; \
                int fn = snprintf(flags_buf + fpos, sizeof(flags_buf) - fpos, \
                                  "\"%s\"", name); \
                fpos += fn; nflags++; \
            } \
        } while(0)

        FLAG(NUC_FLAG_MUTATED,      "mutated");
        FLAG(NUC_FLAG_SEQ_ERROR,    "seq_error");
        FLAG(NUC_FLAG_PCR_ERROR,    "pcr_error");
        FLAG(NUC_FLAG_P_NUCLEOTIDE, "p_nuc");
        FLAG(NUC_FLAG_N_NUCLEOTIDE, "n_nuc");
        FLAG(NUC_FLAG_IS_N,         "is_n");
        FLAG(NUC_FLAG_ANCHOR,       "anchor");
        FLAG(NUC_FLAG_INDEL_INS,    "indel_ins");

        #undef FLAG

        flags_buf[fpos++] = ']';
        flags_buf[fpos] = '\0';

        EMIT("{\"pos\":%d,\"cur\":\"%c\",\"germ\":\"%c\","
             "\"seg\":\"%s\",\"flags\":%s,"
             "\"gp\":%d,\"ph\":%d,\"aa\":\"%c\",\"prod\":%s}",
             seq_pos,
             n->current ? n->current : '?',
             n->germline ? n->germline : '.',
             (n->segment < SEG_COUNT) ? seg_names[n->segment] : "?",
             flags_buf,
             (int)n->germline_pos,
             (int)n->frame_phase,
             n->amino_acid ? n->amino_acid : '.',
             n->productive ? "true" : "false");

        seq_pos++;
    }

    EMIT("]");
    #undef EMIT
    return pos;

truncated:
    /* Ensure null-terminated even on truncation */
    if (bufsize > 0) buf[bufsize - 1] = '\0';
    return -1;
}

int genairr_dump_snapshot(GenAIRRSimulator *sim, int index,
                          char *buffer, int buffer_size) {
    if (!sim || !buffer || buffer_size <= 0) return -1;
    if (index < 0 || index >= sim->n_snapshots) return -1;

    return aseq_to_json(&sim->snapshots[index].seq, buffer, buffer_size);
}

int genairr_dump_snapshot_codon_rail(GenAIRRSimulator *sim, int index,
                                     char *buffer, int buffer_size) {
    if (!sim || !buffer || buffer_size <= 0) return -1;
    if (index < 0 || index >= sim->n_snapshots) return -1;

    const ASeq *seq = &sim->snapshots[index].seq;
    int pos = 0;

    #define EMIT(...) do { \
        int _w = snprintf(buffer + pos, buffer_size - pos, __VA_ARGS__); \
        if (_w < 0 || pos + _w >= buffer_size) goto truncated; \
        pos += _w; \
    } while(0)

    EMIT("[");

    int codon_idx = 0;
    for (const Nuc *n = seq->codon_head; n; n = n->codon_next) {
        if (codon_idx > 0) EMIT(",");

        /* Extract the 3 bases of this codon */
        char b0 = n->current ? n->current : '?';
        char b1 = n->next ? n->next->current : '?';
        char b2 = (n->next && n->next->next) ? n->next->next->current : '?';

        EMIT("{\"idx\":%d,\"bases\":\"%c%c%c\",\"aa\":\"%c\","
             "\"is_stop\":%s,\"seg\":\"%s\"}",
             codon_idx,
             b0, b1, b2,
             n->amino_acid ? n->amino_acid : '.',
             (n->amino_acid == '*') ? "true" : "false",
             (n->segment < SEG_COUNT) ? seg_names[n->segment] : "?");

        codon_idx++;
    }

    EMIT("]");
    #undef EMIT
    return pos;

truncated:
    if (buffer_size > 0) buffer[buffer_size - 1] = '\0';
    return -1;
}

/* ── Version ─────────────────────────────────────────────────── */

const char *genairr_version(void) {
    return GENAIRR_C_VERSION_STRING;
}
