/**
 * strategy_motif_v.c — V-anchor motif fallback for ungapped references.
 *
 * Used when neither an explicit anchor nor an IMGT-gapped derivation
 * is available (plain custom FASTA, OGRDB without sidecar, etc.).
 *
 * Heuristic: V genes are conventionally trimmed to ~285 nt at the
 * 2nd-Cys; the conserved Cys anchor lies in the LAST ~30 codons of the
 * gene. Scanning that window in frame 0 for {TGT, TGC, TGG} and
 * accepting if exactly one candidate exists gives a reasonably reliable
 * best-guess.
 *
 * Multiple candidates → REJECTED with a descriptive reason rather than
 * picking arbitrarily — better to surface the ambiguity than silently
 * mis-anchor.
 */

#include <stddef.h>   /* NULL */
#include "strategies.h"
#include "codon.h"

/* Search window: last 30 codons (90 nt) of the ungapped sequence. */
#define MOTIF_V_WINDOW_CODONS  30
#define MOTIF_V_MIN_LEN        60   /* need ≥ 20 codons for the heuristic */

AnchorResult try_v_motif(const AnchorResolverConfig *cfg,
                         const LoadedAlleleRecord *rec) {
    if (!rec || !rec->sequence || rec->sequence_length < MOTIF_V_MIN_LEN) {
        return anchor_make_rejected(NULL, METHOD_MOTIF_SEARCH);
    }

    const int seq_len = rec->sequence_length;
    /* Frame 0 scan: codon starts at multiples of 3. Constrain to the
     * last MOTIF_V_WINDOW_CODONS codons. */
    int last_codon_start = (seq_len / 3) * 3 - 3;     /* index of last full codon */
    if (last_codon_start < 0) {
        return anchor_make_rejected(NULL, METHOD_MOTIF_SEARCH);
    }
    int window_start = last_codon_start - (MOTIF_V_WINDOW_CODONS - 1) * 3;
    if (window_start < 0) window_start = 0;
    /* Ensure window_start aligns to frame 0. */
    window_start -= window_start % 3;

    int best_pos = -1;
    int n_candidates = 0;
    AnchorResidue best_residue = ANCHOR_RES_UNKNOWN;
    char best_codon[4] = "???";

    for (int i = window_start; i + 3 <= seq_len; i += 3) {
        const char *codon = rec->sequence + i;
        if (codon_is_stop(codon)) continue;
        AnchorResidue r = codon_classify(codon);
        if (r == ANCHOR_CYS || r == ANCHOR_TRP) {
            n_candidates++;
            if (best_pos < 0) {
                best_pos = i;
                best_residue = r;
                codon_normalize(codon, best_codon);
            }
        }
    }

    if (n_candidates == 0) {
        return anchor_make_rejected(
            "no Cys/Trp codon found in last 30 codons of V sequence",
            METHOD_MOTIF_SEARCH);
    }
    if (n_candidates > 1) {
        return anchor_make_rejected(
            "ambiguous V anchor: multiple Cys/Trp codons in search window",
            METHOD_MOTIF_SEARCH);
    }
    if (best_residue == ANCHOR_TRP && cfg && !cfg->allow_trp_v) {
        return anchor_make_rejected(
            "Trp anchor disallowed by resolver config",
            METHOD_MOTIF_SEARCH);
    }
    /* Single canonical candidate found — return as best_guess. */
    return anchor_make_result(best_pos, best_codon, best_residue,
                              CONF_BEST_GUESS, METHOD_MOTIF_SEARCH);
}
