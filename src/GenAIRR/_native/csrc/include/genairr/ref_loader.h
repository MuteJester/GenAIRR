/**
 * ref_loader.h — pluggable reference loader vtable + build report.
 *
 * Loaders read AIRR-format reference data (IMGT, OGRDB, AIRR-C, IgBLAST,
 * VDJbase, plain FASTA) and produce a uniform stream of
 * LoadedAlleleRecord values via a strategy-pattern vtable.
 *
 * OWNERSHIP CONTRACT
 * ------------------
 * 1. Each loader owns all char* memory referenced from any
 *    LoadedAlleleRecord it produces via next(). The pointers are
 *    valid only until the NEXT call to next() on the same loader,
 *    or until close().
 * 2. The caller MUST copy any strings they want to retain before
 *    calling next() again or close().
 * 3. The loader struct is heap-allocated by *_open and freed by
 *    close(). close() is idempotent and NULL-safe.
 * 4. err_msg pointers returned from open()/next() are static literals
 *    or loader-owned and remain valid until close().
 */

#ifndef GENAIRR_REF_LOADER_H
#define GENAIRR_REF_LOADER_H

#include "genairr/ref_record.h"
#include "genairr/anchor.h"
#include "genairr/export.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ── Loader vtable ────────────────────────────────────────────── */

typedef struct ReferenceLoader ReferenceLoader;

typedef struct {
    /* Get the next record. Returns:
     *    1  EOF (no more records — out unmodified)
     *    0  success (out filled, valid until next call)
     *   -1  error (*err_msg populated; loader is left in an undefined
     *              state and the caller should close() it)
     */
    int  (*next)(ReferenceLoader *self,
                 LoadedAlleleRecord *out,
                 const char **err_msg);

    /* Free all resources. Idempotent and NULL-safe. */
    void (*close)(ReferenceLoader *self);
} ReferenceLoaderVTable;

struct ReferenceLoader {
    const ReferenceLoaderVTable *vt;
    void *state;     /* loader-private; opaque to callers */
};

/* ── Concrete loader factories ───────────────────────────────── */

/* IMGT V-QUEST FASTA loader. Reads pipe-delimited headers and
 * IMGT-gapped sequences. Filters out partial alleles (`partial`
 * keyword in header) and non-functional (status not in {F, ORF})
 * automatically; pseudogenes are NOT preserved in Phase 1A.
 *
 * `segment_hint` is used when a per-allele segment cannot be inferred
 * from the gene name (rare in standard IMGT data). Pass SEG_UNKNOWN
 * to require name-based inference for every record.
 *
 * Returns NULL on file open failure or malformed data.
 * *err_msg points at a static string in the failure case. */
GENAIRR_EXPORT ReferenceLoader *imgt_vquest_loader_open(
    const char *fasta_path,
    Segment segment_hint,
    const char **err_msg);

/* Plain (bare) FASTA loader. Used for custom lab references with no
 * IMGT alignment, no pipe-delimited metadata. The caller must supply
 * a `locus_hint` so locus-aware J anchor resolution stays meaningful;
 * pass LOCUS_UNKNOWN if even the locus is uncertain (J anchor scan
 * downgrades to CONF_BEST_GUESS in that case). */
GENAIRR_EXPORT ReferenceLoader *plain_fasta_loader_open(
    const char *fasta_path,
    Locus locus_hint,
    Segment segment_hint,
    const char **err_msg);

/* AIRR-C GermlineSet JSON loader (Phase 2b). Reads the modern
 * source-of-truth schema with explicit V/J anchor coordinates. */
GENAIRR_EXPORT ReferenceLoader *airrc_germline_loader_open(
    const char *json_path,
    const char **err_msg);

/* OGRDB combined loader (Phase 2c). Reads a FASTA + AIRR-C-style JSON
 * sidecar; falls back to motif resolution if the sidecar is missing. */
GENAIRR_EXPORT ReferenceLoader *ogrdb_loader_open(
    const char *fasta_path,
    const char *json_sidecar_path,        /* may be NULL */
    Locus locus_hint,
    Segment segment_hint,
    const char **err_msg);

/* IgBLAST custom-DB loader (Phase 3a). Reads an IgBLAST-format FASTA
 * (plain headers, ungapped) plus optionally one of:
 *   - .aux  (J anchors, 0-based cdr3_end)
 *   - .ndm.imgt (V coords, 1-based)
 * Pass NULL for `aux_or_ndm_path` if no sidecar is available;
 * resolution falls back to motif search. */
GENAIRR_EXPORT ReferenceLoader *igblast_loader_open(
    const char *fasta_path,
    const char *aux_or_ndm_path,
    Locus locus_hint,
    Segment segment_hint,
    const char **err_msg);

/* Thin C-callable wrappers around the vtable for ctypes consumers
 * (Python bindings). Internal C code can call self->vt->next/close
 * directly; ctypes can't dereference function-pointer tables in a
 * struct, so we expose flat entry points. */
GENAIRR_EXPORT int  reference_loader_next(ReferenceLoader *self,
                                          LoadedAlleleRecord *out,
                                          const char **err_msg);
GENAIRR_EXPORT void reference_loader_close(ReferenceLoader *self);

/* Phase 2 (placeholders for future use; not implemented yet):
 *   airrc_germline_loader_open(json_path, ...)
 *   ogrdb_loader_open(fasta_path, json_sidecar_path, ...)
 *   plain_fasta_loader_open(fasta_path, locus_hint, ...)
 *
 * Phase 3:
 *   igblast_bundle_loader_open(bundle_dir, ...)
 *   vdjbase_loader_open(fasta_path, json_sidecar_path, ...)
 */

/* ── Loader auto-dispatch (Phase 2+) ─────────────────────────── */

typedef struct {
    Locus   locus_hint;       /* used by formats that don't carry locus */
    Segment segment_hint;
} LoaderHints;

/* Auto-detection — pick the right concrete loader by inspecting
 * the path's extension and the file's first few bytes. Falls back
 * to plain FASTA. Returns NULL on error with *err_msg populated. */
GENAIRR_EXPORT ReferenceLoader *reference_loader_open_auto(
    const char *path,
    const LoaderHints *hints,
    const char **err_msg);

/* ── Locus / segment inference ────────────────────────────────── */

/* Infer locus from an IMGT-style allele or gene name.
 *   "IGHV1-2*01" -> LOCUS_IGH
 *   "TRBV20-1*01" -> LOCUS_TRB
 *   etc.
 * Returns LOCUS_UNKNOWN if the name doesn't match a known prefix. */
Locus locus_from_gene_name(const char *name);

/* Infer locus from a filename's last path component (heuristic).
 *   "/data/IGH_V.fasta" -> LOCUS_IGH
 *   "human_TRBV.fa"     -> LOCUS_TRB
 * Returns LOCUS_UNKNOWN on no match. */
Locus locus_from_filename(const char *path);

/* Infer segment (V/D/J/C) from an IMGT-style gene name. */
Segment segment_from_gene_name(const char *name);

/* ── BuildReport ──────────────────────────────────────────────── */

#define BUILD_REPORT_REASON_LEN 200

typedef struct {
    char allele_name[64];
    char reason[BUILD_REPORT_REASON_LEN];
    int  position;
    char codon[4];
} BuildReportRejection;

/* Bucket array sizes — pick comfortably larger than the biggest enum
 * value so future additions don't silently overflow. The actual enum
 * upper bounds today are SEG_COUNT (=8 in types.h), FUNC_PARTIAL (=4),
 * CONF_CANONICAL (=3); we round up. */
#define BUILD_REPORT_SEG_BUCKETS    9
#define BUILD_REPORT_STATUS_BUCKETS 5
#define BUILD_REPORT_CONF_BUCKETS   4

typedef struct {
    int total_seen;
    int accepted;
    int by_segment[BUILD_REPORT_SEG_BUCKETS];
    int by_status[BUILD_REPORT_STATUS_BUCKETS];
    int by_confidence[BUILD_REPORT_CONF_BUCKETS];

    BuildReportRejection *rejections;
    int n_rejections;
    int cap_rejections;
} BuildReport;

void build_report_init(BuildReport *r);

void build_report_record_accepted(BuildReport *r,
                                  const LoadedAlleleRecord *rec,
                                  const AnchorResult *result);

void build_report_record_rejected(BuildReport *r,
                                  const LoadedAlleleRecord *rec,
                                  const AnchorResult *result);

/* Free dynamically-allocated rejection storage. NULL-safe. */
void build_report_destroy(BuildReport *r);

#ifdef __cplusplus
}
#endif

#endif /* GENAIRR_REF_LOADER_H */
