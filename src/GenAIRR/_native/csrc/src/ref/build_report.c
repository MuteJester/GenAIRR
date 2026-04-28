/**
 * build_report.c — aggregate per-allele resolution outcomes.
 *
 * The report is a Builder-style accumulator: callers feed it the
 * outcome of each LoadedAlleleRecord (accepted or rejected) and the
 * report tracks counts plus a bounded list of rejection details for
 * later inspection. Counts use plain int arrays indexed by enum
 * values — allocator-free, fast, easy to dump as JSON.
 *
 * Memory: rejections array grows on demand (doubling). Caller must
 * call build_report_destroy() when done.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "genairr/ref_loader.h"

/* Initial and growth policy for the rejection vector. */
#define BUILD_REPORT_INITIAL_CAP    16

void build_report_init(BuildReport *r) {
    if (!r) return;
    memset(r, 0, sizeof(*r));
    /* Count arrays are zeroed by memset; rejections start NULL. */
}

static int rejections_grow(BuildReport *r) {
    int new_cap = r->cap_rejections ? r->cap_rejections * 2
                                    : BUILD_REPORT_INITIAL_CAP;
    BuildReportRejection *grown = (BuildReportRejection *)realloc(
        r->rejections, sizeof(BuildReportRejection) * (size_t)new_cap);
    if (!grown) return -1;
    r->rejections = grown;
    r->cap_rejections = new_cap;
    return 0;
}

/* Safe copy with truncation. dst must be sized correctly. */
static void copy_truncated(char *dst, int dst_size, const char *src) {
    if (!dst || dst_size <= 0) return;
    if (!src) { dst[0] = '\0'; return; }
    int n = (int)strlen(src);
    if (n >= dst_size) n = dst_size - 1;
    memcpy(dst, src, (size_t)n);
    dst[n] = '\0';
}

/* Bound-check before indexing into the by_segment / by_status /
 * by_confidence count arrays — protects against future enum values
 * that exceed the static array sizes. */
static void bump_counter(int *counters, int max, int idx) {
    if (idx >= 0 && idx < max) counters[idx]++;
}

void build_report_record_accepted(BuildReport *r,
                                  const LoadedAlleleRecord *rec,
                                  const AnchorResult *result) {
    if (!r || !rec || !result) return;
    r->total_seen++;
    r->accepted++;
    bump_counter(r->by_segment,    BUILD_REPORT_SEG_BUCKETS,    (int)rec->segment);
    bump_counter(r->by_status,     BUILD_REPORT_STATUS_BUCKETS, (int)rec->functional_status);
    bump_counter(r->by_confidence, BUILD_REPORT_CONF_BUCKETS,   (int)result->confidence);
}

void build_report_record_rejected(BuildReport *r,
                                  const LoadedAlleleRecord *rec,
                                  const AnchorResult *result) {
    if (!r || !rec || !result) return;
    r->total_seen++;
    bump_counter(r->by_segment,    BUILD_REPORT_SEG_BUCKETS,    (int)rec->segment);
    bump_counter(r->by_status,     BUILD_REPORT_STATUS_BUCKETS, (int)rec->functional_status);
    bump_counter(r->by_confidence, BUILD_REPORT_CONF_BUCKETS,   (int)result->confidence);

    if (r->n_rejections >= r->cap_rejections) {
        if (rejections_grow(r) != 0) {
            /* OOM — silently drop the detail entry. The aggregate
             * count is still recorded. */
            return;
        }
    }
    BuildReportRejection *slot = &r->rejections[r->n_rejections++];
    copy_truncated(slot->allele_name, (int)sizeof(slot->allele_name),
                   rec->name);
    copy_truncated(slot->reason, BUILD_REPORT_REASON_LEN,
                   result->reason ? result->reason : "(no reason given)");
    slot->position = result->position;
    /* Codon is fixed-size; copy verbatim, ensure NUL. */
    memcpy(slot->codon, result->codon, sizeof(slot->codon));
    slot->codon[sizeof(slot->codon) - 1] = '\0';
}

void build_report_destroy(BuildReport *r) {
    if (!r) return;
    free(r->rejections);
    r->rejections = NULL;
    r->n_rejections = 0;
    r->cap_rejections = 0;
}
