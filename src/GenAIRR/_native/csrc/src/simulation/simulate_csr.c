/**
 * simulate_csr.c — Isotype-aware mutation rate adjustment for CSR.
 *
 * Provides csr_adjust_rates() which modifies mutation rate parameters
 * based on the selected C allele isotype. Called by the pipeline
 * before mutation to apply isotype-specific SHM rates.
 *
 * This is NOT a pipeline step (doesn't match StepFn signature) — it's
 * a utility called by higher-level code that has access to both the
 * SimConfig and the S5F model.
 *
 * Default isotype rates (from Python SimulateCSR):
 *   IGHM:  0.001 - 0.03
 *   IGHD:  0.001 - 0.03
 *   IGHG3: 0.05  - 0.12
 *   IGHG1: 0.05  - 0.12
 *   IGHA1: 0.04  - 0.10
 *   IGHG2: 0.05  - 0.12
 *   IGHG4: 0.05  - 0.12
 *   IGHE:  0.05  - 0.15
 *   IGHA2: 0.04  - 0.10
 */

#include "genairr/pipeline.h"
#include "genairr/trace.h"
#include <string.h>

typedef struct {
    const char *gene;
    double min_rate;
    double max_rate;
} IsotypeRate;

static const IsotypeRate ISOTYPE_RATES[] = {
    { "IGHM",  0.001, 0.03 },
    { "IGHD",  0.001, 0.03 },
    { "IGHG1", 0.05,  0.12 },
    { "IGHG2", 0.05,  0.12 },
    { "IGHG3", 0.05,  0.12 },
    { "IGHG4", 0.05,  0.12 },
    { "IGHA1", 0.04,  0.10 },
    { "IGHA2", 0.04,  0.10 },
    { "IGHE",  0.05,  0.15 },
    { NULL,    0.0,   0.0  },
};

void csr_adjust_rates(const SimRecord *rec,
                      double *min_rate, double *max_rate) {
    if (!rec->c_allele) return;

    /* Extract gene name (before '*') */
    char gene[64];
    const char *name = rec->c_allele->name;
    const char *star = strchr(name, '*');
    int len = star ? (int)(star - name) : (int)strlen(name);
    if (len >= 64) len = 63;
    memcpy(gene, name, len);
    gene[len] = '\0';

    /* Look up isotype-specific rates */
    double old_min = *min_rate, old_max = *max_rate;
    for (const IsotypeRate *r = ISOTYPE_RATES; r->gene; r++) {
        if (strcmp(gene, r->gene) == 0) {
            *min_rate = r->min_rate;
            *max_rate = r->max_rate;
            TRACE("[csr] isotype=%s, adjusted rates [%.4f, %.4f] → [%.4f, %.4f]",
                  gene, old_min, old_max, r->min_rate, r->max_rate);
            return;
        }
    }
    /* Unknown isotype: keep original rates */
    TRACE("[csr] isotype=%s (unknown), keeping rates [%.4f, %.4f]", gene, old_min, old_max);
}
