/**
 * genairr_api.h — Public C API for GenAIRR Python bindings.
 *
 * Opaque-handle API: no struct layout knowledge needed by callers.
 * Python calls these via ctypes; all complexity is internal.
 *
 * Lifecycle:
 *   1. genairr_create(gdc_path)       → opaque handle
 *   2. genairr_set_*()                → configure features + params
 *   3. genairr_simulate(handle, n, …) → run, write TSV
 *   4. genairr_destroy(handle)        → cleanup
 */

#ifndef GENAIRR_API_H
#define GENAIRR_API_H

#include <stdint.h>

/* ── DLL export macro (Windows needs __declspec(dllexport)) ──── */
#ifdef _WIN32
  #define GENAIRR_EXPORT __declspec(dllexport)
#else
  #define GENAIRR_EXPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ── Opaque handle ───────────────────────────────────────────── */

typedef struct GenAIRRSimulator GenAIRRSimulator;

/* ── Lifecycle ───────────────────────────────────────────────── */

/**
 * Create a simulator from a .gdc config file.
 * Returns NULL on error (bad path or corrupt file).
 */
GENAIRR_EXPORT GenAIRRSimulator *genairr_create(const char *gdc_path);

/**
 * Create a simulator from in-memory .gdc bytes.
 * Returns NULL on error.
 */
GENAIRR_EXPORT GenAIRRSimulator *genairr_create_from_memory(const void *gdc_bytes,
                                              int gdc_len);

/**
 * Destroy a simulator and free all resources.
 */
GENAIRR_EXPORT void genairr_destroy(GenAIRRSimulator *sim);

/* ── Feature flags ───────────────────────────────────────────── */

/**
 * Enable/disable a named feature. Feature names match Python op names:
 *
 *   "productive", "mutate", "selection_pressure", "csr",
 *   "d_inversion", "receptor_revision",
 *   "corrupt_5_prime", "corrupt_3_prime", "quality_errors",
 *   "paired_end", "pcr", "umi", "primer_mask",
 *   "reverse_complement", "contaminants", "long_read_errors",
 *   "indels", "insert_ns"
 *
 * Returns 0 on success, -1 if name is unrecognized.
 */
GENAIRR_EXPORT int genairr_set_feature(GenAIRRSimulator *sim,
                         const char *name, int enabled);

/* ── Parameters ──────────────────────────────────────────────── */

/**
 * Set a named float64 parameter. Parameter names:
 *
 *   "min_mutation_rate", "max_mutation_rate",
 *   "selection_strength", "cdr_r_acceptance", "fwr_r_acceptance",
 *   "d_inversion_prob", "revision_prob",
 *   "base_error_rate", "peak_error_rate", "transition_weight",
 *   "pcr_error_rate", "contamination_prob",
 *   "indel_prob", "n_prob", "rc_prob",
 *   "long_read_error_rate", "insertion_bias",
 *   "long_read_error_rate", "insertion_bias"
 *
 * Returns 0 on success, -1 if name is unrecognized.
 */
GENAIRR_EXPORT int genairr_set_param_f64(GenAIRRSimulator *sim,
                           const char *name, double value);

/**
 * Set a named int32 parameter. Parameter names:
 *
 *   "pcr_cycles", "umi_length", "pe_read_length",
 *   "max_sequence_length", "footprint_min", "footprint_max",
 *   "min_run_length", "primer_mask_length",
 *   "max_productive_attempts", "contaminant_type",
 *   "corrupt_5_remove_min/max", "corrupt_5_add_min/max",
 *   "corrupt_3_remove_min/max", "corrupt_3_add_min/max"
 *
 * Returns 0 on success, -1 if name is unrecognized.
 */
GENAIRR_EXPORT int genairr_set_param_i32(GenAIRRSimulator *sim,
                           const char *name, int value);

/* ── Allele locking ─────────────────────────────────────────── */

/**
 * Lock sampling of a segment to specific allele(s) by name.
 *
 * @param sim      The simulator handle.
 * @param segment  Segment name: "v", "d", "j", or "c".
 * @param name     Allele name (e.g. "IGHV1-2*01"). Must exist in the pool.
 * @return 0 on success, -1 if segment or allele name is invalid,
 *         -2 if too many alleles locked (max 64).
 *
 * Call multiple times to allow multiple alleles for one segment.
 * The simulator will pick uniformly among locked alleles.
 */
GENAIRR_EXPORT int genairr_lock_allele(GenAIRRSimulator *sim,
                         const char *segment, const char *name);

/**
 * Clear all locked alleles for a segment, restoring random sampling.
 * Pass NULL for segment to clear all segments.
 */
GENAIRR_EXPORT void genairr_clear_locks(GenAIRRSimulator *sim, const char *segment);

/* ── Seed ────────────────────────────────────────────────────── */

/**
 * Set the random seed. 0 = use time-based seed.
 */
GENAIRR_EXPORT void genairr_set_seed(GenAIRRSimulator *sim, uint64_t seed);

/* ── Simulate ────────────────────────────────────────────────── */

/**
 * Run the simulation and write AIRR TSV to the given file path.
 *
 * @param sim          The configured simulator.
 * @param n            Number of sequences to generate.
 * @param output_path  Path for the output TSV file.
 * @return Number of sequences written, or -1 on error.
 */
GENAIRR_EXPORT int genairr_simulate(GenAIRRSimulator *sim, int n,
                      const char *output_path);

/**
 * Run the simulation and write AIRR TSV to a file descriptor.
 * Useful for writing to stdout or a pipe.
 *
 * @param sim   The configured simulator.
 * @param n     Number of sequences to generate.
 * @param fd    Open file descriptor to write to.
 * @return Number of sequences written, or -1 on error.
 */
GENAIRR_EXPORT int genairr_simulate_to_fd(GenAIRRSimulator *sim, int n, int fd);

/**
 * Simulate one sequence and write the AIRR TSV row to a buffer.
 *
 * No header or newline is written — just the tab-separated field values.
 * On the first call (or after genairr_set_seed), the RNG is seeded.
 * Subsequent calls continue from the current RNG state, enabling
 * streaming without re-seeding.
 *
 * @param sim          The configured simulator.
 * @param buffer       Output buffer for the TSV row.
 * @param buffer_size  Size of the buffer in bytes.
 * @return Number of bytes written (excl. null terminator),
 *         or -1 on error / buffer too small.
 */
GENAIRR_EXPORT int genairr_simulate_one(GenAIRRSimulator *sim,
                          char *buffer, int buffer_size);

/**
 * Write the AIRR TSV column names (tab-separated) to a buffer.
 *
 * @param buffer       Output buffer.
 * @param buffer_size  Size of the buffer in bytes.
 * @return Number of columns, or -1 on error.
 */
GENAIRR_EXPORT int genairr_get_header(char *buffer, int buffer_size);

/* ── Tracing ─────────────────────────────────────────────────── */

/**
 * Enable or disable simulation tracing.
 *
 * When enabled, simulate_one() records a detailed human-readable trace
 * of every internal decision (allele selection, trimming, assembly,
 * productivity, mutation, corruption, etc.).
 *
 * The trace buffer is cleared on each simulate_one() call.
 */
GENAIRR_EXPORT void genairr_set_trace(GenAIRRSimulator *sim, int enabled);

/**
 * Copy the trace log to the given output buffer.
 *
 * @return Number of bytes in the trace (may exceed buffer_size),
 *         or 0 if tracing is disabled.
 */
GENAIRR_EXPORT int genairr_get_trace(GenAIRRSimulator *sim,
                       char *buffer, int buffer_size);

/** Clear the trace log buffer without disabling tracing. */
GENAIRR_EXPORT void genairr_clear_trace(GenAIRRSimulator *sim);

/* ── Introspection hooks ─────────────────────────────────────── */

/**
 * Set which hook points to capture as snapshots during simulate_one().
 * Pass a bitmask of (1 << HookPoint) values. Pass 0 to disable all hooks.
 *
 * Hook points: post_assembly(0), post_functionality(1), post_d_inversion(2),
 * post_receptor_rev(3), post_mutation(4), post_selection(5), post_corrupt_5(6),
 * post_corrupt_3(7), post_indels(8), post_ns(9), post_pcr(10),
 * post_quality(11), final(12).
 */
GENAIRR_EXPORT void genairr_set_hooks(GenAIRRSimulator *sim, uint32_t hook_mask);

/** Get the number of captured snapshots from the last simulate_one() call. */
GENAIRR_EXPORT int genairr_get_snapshot_count(GenAIRRSimulator *sim);

/** Get the hook point name for a snapshot by index. */
GENAIRR_EXPORT const char *genairr_get_snapshot_name(GenAIRRSimulator *sim, int index);

/**
 * Dump a snapshot's ASeq as a JSON array of node objects.
 * Each node: {"pos":N,"cur":"A","germ":"A","seg":"V","flags":[...],"gp":N,"ph":N,"aa":"C","prod":bool}
 * @return Number of bytes written, or -1 on error/truncation.
 */
GENAIRR_EXPORT int genairr_dump_snapshot(GenAIRRSimulator *sim, int index,
                          char *buffer, int buffer_size);

/**
 * Dump a snapshot's codon rail as a JSON array.
 * Each codon: {"idx":N,"bases":"ATG","aa":"M","is_stop":false,"seg":"V"}
 * @return Number of bytes written, or -1 on error/truncation.
 */
GENAIRR_EXPORT int genairr_dump_snapshot_codon_rail(GenAIRRSimulator *sim, int index,
                                     char *buffer, int buffer_size);

/* ── Version info ────────────────────────────────────────────── */

/** Return version string (e.g. "0.1.0"). Static storage. */
GENAIRR_EXPORT const char *genairr_version(void);

#ifdef __cplusplus
}
#endif

#endif /* GENAIRR_API_H */
