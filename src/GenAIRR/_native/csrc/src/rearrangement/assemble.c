/**
 * assemble.c — Sequence assembly and functionality assessment.
 */

#include "genairr/pipeline.h"
#include "genairr/functionality.h"
#include "genairr/trace.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ── Internal: generate a random NP region ────────────────────── */

static int generate_np(char *buf, int mean_len, int max_len) {
    static const char bases[] = "ACGT";

    int len;
    if (max_len <= 0) {
        len = 0;
    } else {
        int lo = mean_len > 3 ? mean_len - 3 : 0;
        int hi = mean_len + 3 < max_len ? mean_len + 3 : max_len;
        if (hi <= lo) hi = lo + 1;
        len = lo + rand() % (hi - lo + 1);
    }

    for (int i = 0; i < len; i++) {
        buf[i] = bases[rand() % 4];
    }
    buf[len] = '\0';
    return len;
}

/* ── step_assemble ────────────────────────────────────────────── */

void step_assemble(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    bool has_d = chain_has_d(cfg->chain_type) && rec->d_allele != NULL;

    /* ── Append V segment (trimmed) ─────────────────────────── */
    if (rec->v_allele) {
        const char *v_seq = rec->v_allele->seq;
        int v_start = rec->v_trim_5;
        int v_end   = rec->v_allele->length - rec->v_trim_3;
        int v_len   = v_end - v_start;

        if (v_len > 0) {
            aseq_append_segment(seq, v_seq + v_start, v_len,
                                SEG_V, v_start, rec->v_allele->anchor);
        }
        TRACE("[assemble] V: %dbp (germline[%d:%d], anchor_at=%d)",
              v_len, v_start, v_end, rec->v_allele->anchor);
    }

    /* ── Generate and append NP1 ────────────────────────────── */
    {
        char np1_buf[64];
        int np1_len = generate_np(np1_buf, cfg->np1_length_mean,
                                  cfg->np1_length_max);
        if (np1_len > 0) {
            aseq_append_np(seq, np1_buf, np1_len,
                           SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
        }
        rec->np1_length = np1_len;
        TRACE("[assemble] NP1: %dbp (%s)", np1_len, np1_len > 0 ? np1_buf : "");
    }

    /* ── Append D segment (trimmed) — VDJ only ─────────────── */
    if (has_d) {
        const char *d_seq = rec->d_allele->seq;
        int d_start = rec->d_trim_5;
        int d_end   = rec->d_allele->length - rec->d_trim_3;
        int d_len   = d_end - d_start;

        if (d_len > 0) {
            aseq_append_segment(seq, d_seq + d_start, d_len,
                                SEG_D, d_start, -1);
        }
        TRACE("[assemble] D: %dbp (germline[%d:%d])", d_len, d_start, d_end);

        /* ── Generate and append NP2 ────────────────────────── */
        char np2_buf[64];
        int np2_len = generate_np(np2_buf, cfg->np2_length_mean,
                                  cfg->np2_length_max);
        if (np2_len > 0) {
            aseq_append_np(seq, np2_buf, np2_len,
                           SEG_NP2, NUC_FLAG_N_NUCLEOTIDE);
        }
        rec->np2_length = np2_len;
        TRACE("[assemble] NP2: %dbp (%s)", np2_len, np2_len > 0 ? np2_buf : "");
    }

    /* ── Append J segment (trimmed) ─────────────────────────── */
    if (rec->j_allele) {
        const char *j_seq = rec->j_allele->seq;
        int j_start = rec->j_trim_5;
        int j_len   = rec->j_allele->length - rec->j_trim_5 - rec->j_trim_3;

        if (j_len > 0) {
            aseq_append_segment(seq, j_seq + j_start, j_len,
                                SEG_J, j_start, rec->j_allele->anchor);
        }
        TRACE("[assemble] J: %dbp (germline[%d:%d], anchor_at=%d)",
              j_len, j_start, j_start + j_len, rec->j_allele->anchor);
    }

    TRACE("[assemble] total: %dbp = V(%d) + NP1(%d) + D(%d) + NP2(%d) + J(%d)",
          seq->length,
          aseq_segment_length(seq, SEG_V),
          rec->np1_length,
          has_d ? aseq_segment_length(seq, SEG_D) : 0,
          has_d ? rec->np2_length : 0,
          aseq_segment_length(seq, SEG_J));
}

/* ── step_assess_functionality ────────────────────────────────── */

void step_assess_functionality(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    (void)cfg;

    if (!seq->codon_rail_valid) {
        aseq_build_codon_rail(seq);
    }

    static FunctionalityValidator validator;
    static bool initialized = false;
    if (!initialized) {
        validator = functionality_validator_default();
        initialized = true;
    }

    functionality_assess(&validator, seq, rec);

    TRACE("[assess] productive=%s, stop_codon=%s, vj_in_frame=%s%s%s",
          rec->productive ? "yes" : "no",
          rec->stop_codon ? "yes" : "no",
          rec->vj_in_frame ? "yes" : "no",
          rec->note[0] ? ", note=" : "",
          rec->note[0] ? rec->note : "");
}
