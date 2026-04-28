/**
 * ref_record.h — uniform input record for the anchor resolver.
 *
 * Each ReferenceLoader produces a stream of LoadedAlleleRecord values.
 * The resolver and downstream DataConfig builder consume them without
 * caring about the source format (IMGT / OGRDB / AIRR-C / IgBLAST / …).
 *
 * OWNERSHIP CONTRACT
 * ------------------
 * All char* fields in a LoadedAlleleRecord are loader-owned. The
 * pointers are valid only until the next call to next() on the same
 * loader, or until close(). Callers must copy any strings they want
 * to retain. See `ref_loader.h` for full ownership rules.
 */

#ifndef GENAIRR_REF_RECORD_H
#define GENAIRR_REF_RECORD_H

#include <stdbool.h>
#include "genairr/anchor.h"   /* Locus, Segment */
#include "genairr/export.h"   /* GENAIRR_EXPORT */

#ifdef __cplusplus
extern "C" {
#endif

/* ── Functional status (IMGT classification) ──────────────────── */

typedef enum {
    FUNC_UNKNOWN = 0,
    FUNC_F,          /* Functional */
    FUNC_ORF,        /* Open Reading Frame */
    FUNC_PSEUDO,     /* Pseudogene */
    FUNC_PARTIAL,    /* Sequence is incomplete */
} FunctionalStatus;

/* ── Loaded allele record ─────────────────────────────────────── */

#define REF_MAX_ALIASES 8

typedef struct LoadedAlleleRecord {
    /* Allele identification */
    const char       *name;                   /* e.g. "IGHV1-2*01" */
    const char       *aliases[REF_MAX_ALIASES];
    int               n_aliases;

    /* Classification */
    Segment           segment;                /* SEG_V / SEG_D / SEG_J / SEG_C */
    Locus             locus;
    const char       *species;                /* may be NULL */

    /* Sequence data */
    const char       *sequence;               /* ungapped, lowercase ACGT(+IUPAC) */
    int               sequence_length;
    const char       *gapped_sequence;        /* NULL if not provided */
    int               gapped_length;
    bool              gap_convention_imgt;    /* true → gapped uses IMGT positions */

    /* Provenance */
    FunctionalStatus  functional_status;
    int               explicit_anchor;        /* -1 if absent; else 0-based ungapped */
    const char       *source;                 /* "imgt-vquest", "ogrdb", … */
} LoadedAlleleRecord;

/* Reset all fields to safe defaults: zeroed, sequence_length=0,
 * explicit_anchor=-1, all enums to UNKNOWN/0. */
GENAIRR_EXPORT void loaded_allele_record_init(LoadedAlleleRecord *r);

/* Human-readable name for a functional status. */
GENAIRR_EXPORT const char *functional_status_name(FunctionalStatus s);

#ifdef __cplusplus
}
#endif

#endif /* GENAIRR_REF_RECORD_H */
