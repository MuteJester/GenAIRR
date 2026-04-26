/**
 * airr.h — AIRR-format serialization from ASeq.
 *
 * Derives integer positions, boundary ambiguity adjustments, and
 * allele call corrections from the annotated sequence structure.
 * Positions are structural properties of the ASeq, computed at
 * serialization time. Allele-call corrections use lookup-table logic.
 */

#ifndef GENAIRR_AIRR_H
#define GENAIRR_AIRR_H

#include "aseq.h"
#include "pipeline.h"
#include "allele_bitmap.h"
#include <stdio.h>

/* ── AIRR position record (derived from ASeq at serialization) ── */

typedef struct {
    /* Sequence positions (0-based, in the assembled sequence) */
    int v_sequence_start;
    int v_sequence_end;        /* may be extended by boundary ambiguity */
    int d_sequence_start;
    int d_sequence_end;
    int j_sequence_start;
    int j_sequence_end;        /* may be extended by boundary ambiguity */

    /* Germline positions (positions within the germline allele) */
    int v_germline_start;
    int v_germline_end;
    int d_germline_start;
    int d_germline_end;
    int j_germline_start;
    int j_germline_end;

    /* Junction */
    int junction_start;
    int junction_end;
    int junction_length;

    /* Trim amounts (may be adjusted by boundary ambiguity) */
    int v_trim_3_adjusted;
    int d_trim_5_adjusted;
    int d_trim_3_adjusted;
    int j_trim_5_adjusted;
} AirrPositions;

/**
 * Derive AIRR positions from the ASeq structure.
 *
 * This single function replaces FixVPosition + FixDPosition + FixJPosition.
 * It walks the ASeq once, finds segment boundaries, then checks for
 * boundary ambiguity (where NP bases match the adjacent trimmed germline).
 *
 * @param seq   The annotated sequence.
 * @param rec   The simulation record (for allele pointers and trim amounts).
 * @param out   Output: filled with derived positions.
 */
void airr_derive_positions(const ASeq *seq, const SimRecord *rec,
                           AirrPositions *out);

/* ── Allele call correction ───────────────────────────────────── */
/*
 * Allele calls are now derived reactively via AlleleBitmap.
 * See allele_bitmap.h for the O(N·L) build + O(L·N/64) query approach
 * that replaces the old O(N²·L) correction maps + pairwise diff tables.
 */

/* ═══════════════════════════════════════════════════════════════
 * AirrRecord — Complete AIRR-format output record.
 *
 * This struct holds all fields needed to produce AIRR-compliant
 * output. It is filled by airr_serialize() which combines the
 * ASeq, SimRecord, positions, and allele corrections into a
 * single output record.
 * ═══════════════════════════════════════════════════════════════ */

#define AIRR_MAX_CALLS       2048   /* comma-separated allele names  */
#define AIRR_MAX_MUTATIONS   4096   /* mutation log string           */
#define AIRR_MAX_JUNCTION     256   /* junction nucleotide string    */
#define AIRR_MAX_NP            64   /* NP region nucleotide string   */

typedef struct {
    /* ── Assembled sequence ───────────────────────────────────── */
    char    sequence[GENAIRR_MAX_SEQ_LEN];
    char    germline_alignment[GENAIRR_MAX_SEQ_LEN];
    int     sequence_length;

    /* ── Allele calls (comma-separated strings) ───────────────── *
     * `v_call` / `d_call` / `j_call` are god-aligner-derived
     * (allele(s) whose germline best matches the survived sequence
     * — what an external aligner would report). Under heavy SHM
     * these can drift away from the truly-sampled allele.
     *
     * `v_call_true` / etc. expose the simulator's ground truth:
     * the name of the allele actually sampled at rearrangement.
     * Use these for AIRR self-consistency checks (seq + mutations
     * must reconstruct the *true* germline). */
    char    v_call[AIRR_MAX_CALLS];
    char    d_call[AIRR_MAX_CALLS];
    char    j_call[AIRR_MAX_CALLS];
    char    c_call[AIRR_MAX_CALLS];
    char    v_call_true[GENAIRR_MAX_ALLELE_NAME];
    char    d_call_true[GENAIRR_MAX_ALLELE_NAME];
    char    j_call_true[GENAIRR_MAX_ALLELE_NAME];
    bool    d_inverted;
    bool    receptor_revised;
    int     revision_footprint_length;
    char    original_v_allele_name[GENAIRR_MAX_ALLELE_NAME];

    /* ── Segment positions (0-based, in assembled sequence) ───── */
    int     v_sequence_start;
    int     v_sequence_end;
    int     v_germline_start;
    int     v_germline_end;

    int     d_sequence_start;
    int     d_sequence_end;
    int     d_germline_start;
    int     d_germline_end;

    int     j_sequence_start;
    int     j_sequence_end;
    int     j_germline_start;
    int     j_germline_end;

    /* ── Junction ─────────────────────────────────────────────── */
    int     junction_start;
    int     junction_end;
    char    junction_nt[AIRR_MAX_JUNCTION];
    char    junction_aa[AIRR_MAX_JUNCTION];  /* translated amino acids */
    int     junction_length;

    /* ── Trim amounts ─────────────────────────────────────────── */
    int     v_trim_5;
    int     v_trim_3;
    int     d_trim_5;
    int     d_trim_3;
    int     j_trim_5;
    int     j_trim_3;

    /* ── NP regions ───────────────────────────────────────────── */
    char    np1_region[AIRR_MAX_NP];
    int     np1_length;
    char    np2_region[AIRR_MAX_NP];
    int     np2_length;

    /* ── Mutation and error annotations ───────────────────────── */
    double  mutation_rate;
    int     n_mutations;
    int     n_sequencing_errors;
    int     n_pcr_errors;
    int     n_insertions;
    int     n_deletions;
    int     n_ns;

    /* Mutation logs: "pos:X>Y,pos:X>Y,..." format */
    char    mutations[AIRR_MAX_MUTATIONS];
    char    sequencing_errors[AIRR_MAX_MUTATIONS];
    char    pcr_errors[AIRR_MAX_MUTATIONS];

    /* ── Productivity ─────────────────────────────────────────── */
    bool    productive;
    bool    stop_codon;
    bool    vj_in_frame;
    char    note[256];

    /* ── Flags ────────────────────────────────────────────────── */
    bool    is_reverse_complement;
    bool    is_contaminant;
} AirrRecord;

/* ── Serialization ────────────────────────────────────────────── */

/**
 * Serialize an ASeq + SimRecord into a complete AirrRecord.
 *
 * Allele calls are derived reactively from the bitmap: for each
 * retained germline position, the bitmap records which alleles share
 * that base. The allele(s) with the most hits = the call. Trimming,
 * mutation, and deletion are automatically reflected because they
 * change which positions contribute votes.
 *
 * Pass NULL for `corr` to skip allele call correction (calls will
 * be the true allele name only).
 *
 * @param seq   The annotated sequence.
 * @param rec   The simulation record.
 * @param cfg   The simulation config (for allele pools).
 * @param corr  Pre-built allele bitmaps (NULL to skip corrections).
 * @param out   Output: filled AirrRecord.
 */
void  airr_serialize(const ASeq *seq, const SimRecord *rec,
                      const SimConfig *cfg, const AlleleCorrectionSet *corr,
                      AirrRecord *out);

/**
 * Write a TSV header line to a file.
 * Returns the number of columns written.
 */
int  airr_write_tsv_header(FILE *fp);

/**
 * Write one AirrRecord as a TSV row to a file.
 */
void  airr_write_tsv_row(FILE *fp, const AirrRecord *rec);

/**
 * Write one AirrRecord as a TSV row to a buffer (no header, no newline).
 * Returns bytes written (excl. null terminator), or -1 if buffer too small.
 */
int   airr_snprintf_tsv_row(char *buf, int size, const AirrRecord *rec);

/**
 * Write the TSV header (column names, tab-separated) to a buffer.
 * Returns number of columns, or -1 if buffer too small.
 */
int   airr_snprintf_tsv_header(char *buf, int size);

#endif /* GENAIRR_AIRR_H */
