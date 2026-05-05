/**
 * assemble.c — Sequence assembly and functionality assessment.
 */

#include "genairr/pipeline.h"
#include "genairr/functionality.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ── Internal: NP region sampling ─────────────────────────────── */

#define GENAIRR_NP_BUF_SIZE 256
/* Lowercase to match V/D/J reference case stored by IMGT/AIRR loaders.
 * Mixing cases (uppercase NP + lowercase germline) silently breaks any
 * case-sensitive comparison between `sequence` and `germline_alignment`
 * — including at extension boundaries where NP bases are claimed as
 * V/D/J. AIRR convention emits sequences in lowercase. */
static const char NP_BASES[] = "acgt";

/* Sample an index in [0, n) from a categorical distribution.
 * Caller is responsible for ensuring `probs` sums to ~1.0; the
 * final-element fallback handles small numerical drift. */
static int sample_categorical(RngState *rng, const double *probs, int n) {
    double u = rng_uniform(rng);
    double cum = 0.0;
    for (int i = 0; i < n; i++) {
        cum += probs[i];
        if (u <= cum) return i;
    }
    return n - 1;
}

/* Uniform fallback used when the SimConfig was built without GDC NP
 * data (legacy paths, manual C tests). Generates a length in
 * [max(0, mean-3), min(max_len, mean+3)] and uniform A/C/G/T bases. */
static int generate_np_uniform(RngState *rng, char *buf, int max_buf) {
    /* Without empirical data, default to a small NP of 3-9 bases —
     * mid-range for human IGH. This path is only hit by manual
     * SimConfigs; production GDC-loaded configs always have np data. */
    const int mean_len = 6;
    const int spread   = 3;
    int hi = mean_len + spread < max_buf - 1 ? mean_len + spread
                                              : max_buf - 1;
    int lo = mean_len - spread > 0 ? mean_len - spread : 0;
    if (hi < lo) hi = lo;
    int len = lo + (int)rng_range(rng, (uint32_t)(hi - lo + 1));
    if (len > hi) len = hi;
    for (int i = 0; i < len; i++) {
        buf[i] = NP_BASES[rng_range(rng, 4)];
    }
    buf[len] = '\0';
    return len;
}

/* Sample an NP length from the empirical distribution, optionally
 * constrained to len % 3 == target_mod_3 (constraint-aware sampling
 * for V5 frame-respecting rearrangement, step 29).
 *
 * target_mod_3 < 0 → unconstrained (matches legacy behaviour).
 * target_mod_3 ∈ {0,1,2} → renormalise empirical probs over allowed
 *   lengths and draw from that filtered distribution.
 *
 * Returns -1 when the empirical distribution has zero mass over the
 * allowed lengths (caller should fall back to unconstrained sampling
 * or accept that no constraint-respecting NP exists for this allele
 * pair).
 */
static int sample_np_length(RngState *rng, const NpDist *np,
                            int target_mod_3) {
    if (np->max_length <= 0) return 0;
    if (target_mod_3 < 0) {
        return sample_categorical(rng, np->length_probs, np->max_length + 1);
    }

    double total = 0.0;
    for (int i = 0; i <= np->max_length; i++) {
        if ((i % 3) == target_mod_3) total += np->length_probs[i];
    }
    if (total < 1e-12) return -1;

    double u = rng_uniform(rng) * total;
    double cum = 0.0;
    int last_valid = -1;
    for (int i = 0; i <= np->max_length; i++) {
        if ((i % 3) != target_mod_3) continue;
        cum += np->length_probs[i];
        last_valid = i;
        if (u <= cum) return i;
    }
    return last_valid;
}

/* Sample an NP region from the empirical TdT Markov model. Length is
 * drawn from `np->length_probs`; the first base from `np->first_base`;
 * subsequent bases from `np->transitions[pos][prev]`. Position is
 * clamped to the last available transition row when the sampled
 * length exceeds n_positions. When a transition row is unobserved
 * (sums to zero), falls back to the marginal first_base distribution.
 *
 * `target_mod_3` constrains the sampled length: see `sample_np_length`.
 * Pass -1 for legacy behaviour. Returns -1 when the constraint cannot
 * be satisfied; caller decides whether to fall back. */
static int generate_np_markov(RngState *rng, char *buf, int max_buf,
                              const NpDist *np, int target_mod_3) {
    if (np->max_length <= 0) {
        buf[0] = '\0';
        return 0;
    }

    /* 1. Sample length from empirical distribution. */
    int len = sample_np_length(rng, np, target_mod_3);
    if (len < 0) {
        /* Constraint unsatisfiable: signal caller. */
        buf[0] = '\0';
        return -1;
    }
    if (len == 0) {
        buf[0] = '\0';
        return 0;
    }
    if (len > max_buf - 1) len = max_buf - 1;

    /* 2. Sample first base. */
    int prev = sample_categorical(rng, np->first_base, 4);
    buf[0] = NP_BASES[prev];

    /* 3. Sample subsequent bases via the Markov transition matrix. */
    for (int i = 1; i < len; i++) {
        int pos = (i - 1 < np->n_positions) ? (i - 1)
                                            : (np->n_positions - 1);
        int next;
        if (np->n_positions == 0 || np->transitions == NULL) {
            /* No transition data; fall back to first_base marginal. */
            next = sample_categorical(rng, np->first_base, 4);
        } else {
            const double *row = &np->transitions[pos * 16 + prev * 4];
            double s = row[0] + row[1] + row[2] + row[3];
            if (s < 1e-12) {
                /* Unobserved (pos, prev) state — back off to the
                 * first-base marginal rather than uniform, since it
                 * is the closest empirical proxy for "unknown TdT
                 * preference at this state". */
                next = sample_categorical(rng, np->first_base, 4);
            } else {
                next = sample_categorical(rng, row, 4);
            }
        }
        buf[i] = NP_BASES[next];
        prev = next;
    }
    buf[len] = '\0';
    return len;
}

/* Generate an NP region (length + bases). When `target_mod_3 >= 0`
 * the empirical Markov path constrains the length; the uniform fallback
 * silently ignores the constraint (it has no length distribution to
 * filter). Returns -1 on a hard constraint failure. */
static int generate_np(RngState *rng, char *buf, int max_buf,
                       const NpDist *np, int target_mod_3) {
    if (np && np->length_probs) {
        return generate_np_markov(rng, buf, max_buf, np, target_mod_3);
    }
    return generate_np_uniform(rng, buf, max_buf);
}

/* ── Internal: P-nucleotide sampling ──────────────────────────── */

#define GENAIRR_P_NUC_BUF 16

static char dna_complement(char b) {
    switch (b) {
        case 'A': return 'T'; case 'a': return 't';
        case 'T': return 'A'; case 't': return 'a';
        case 'C': return 'G'; case 'c': return 'g';
        case 'G': return 'C'; case 'g': return 'c';
        default:  return 'N';
    }
}

/* Sample a P-nucleotide length from the empirical distribution.
 * Returns 0 when the distribution is absent (no P-nuc data in
 * SimConfig) — the assemble step then emits no P-nuc at that end. */
static int sample_p_nuc_length(RngState *rng, const PNucDist *d) {
    if (!d || !d->length_probs || d->max_length <= 0) return 0;
    int k = sample_categorical(rng, d->length_probs, d->max_length + 1);
    if (k < 0) k = 0;
    if (k > GENAIRR_P_NUC_BUF) k = GENAIRR_P_NUC_BUF;
    return k;
}

/* Tail P-nuc (V 3', D 3'): RC of seg[seg_len-K..seg_len-1] read 5'→3'.
 *   p[0] = complement(seg[seg_len-1])
 *   p[1] = complement(seg[seg_len-2])
 *   ...
 * Caller guarantees k <= seg_len and out has at least k bytes. */
static void rc_tail(const char *seg, int seg_len, int k, char *out) {
    for (int i = 0; i < k; i++) {
        out[i] = dna_complement(seg[seg_len - 1 - i]);
    }
}

/* Head P-nuc (D 5', J 5'): RC of seg[0..K-1] read 5'→3'.
 *   p[0] = complement(seg[K-1])
 *   p[1] = complement(seg[K-2])
 *   ...
 *   p[K-1] = complement(seg[0])
 * The resulting K-base string is prepended to the segment, so
 * p[K-1] sits adjacent to seg[0]. */
static void rc_head(const char *seg, int k, char *out) {
    for (int i = 0; i < k; i++) {
        out[i] = dna_complement(seg[k - 1 - i]);
    }
}

/* Append a P-nucleotide of length k derived from `seg` at the given
 * end. Returns the length actually appended (clamped by seg_len and
 * buffer size). When k <= 0 or the segment is empty this is a no-op. */
static int append_p_nuc(ASeq *seq, const char *seg, int seg_len, int k,
                        Segment np_seg, bool head_side) {
    if (k <= 0 || seg_len <= 0) return 0;
    if (k > seg_len) k = seg_len;
    if (k > GENAIRR_P_NUC_BUF) k = GENAIRR_P_NUC_BUF;

    char p_buf[GENAIRR_P_NUC_BUF];
    if (head_side) {
        rc_head(seg, k, p_buf);
    } else {
        rc_tail(seg, seg_len, k, p_buf);
    }
    aseq_append_np(seq, p_buf, k, np_seg, NUC_FLAG_P_NUCLEOTIDE);
    return k;
}

/* ── step_assemble ────────────────────────────────────────────── */

void step_assemble(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    bool has_d = chain_has_d(cfg->chain_type) && rec->d_allele != NULL;

    /* P-nucleotides at the four eligible segment ends are tagged with
     * the surrounding NP segment so the AIRR np1_region / np2_region
     * naturally absorb them. NUC_FLAG_P_NUCLEOTIDE distinguishes them
     * from TdT-generated N-nucleotides for downstream introspection. */
    const PNucDist *p_dist = &cfg->p_nuc_dist;

    /* In a VJ chain (no D segment) all P-nucs sit in NP1; in a VDJ
     * chain D-distal ends (D 3', J 5') sit in NP2. */
    Segment j5_np_seg = has_d ? SEG_NP2 : SEG_NP1;

    int v_p_len = 0, d5_p_len = 0, d3_p_len = 0, j_p_len = 0;

    /* V5 step 29: frame-aware NP1 length sampling for VJ chains.
     *
     * In a VJ chain (no D) the junction length is fully determined by
     * V/J alleles, V/J trims, V3' P-nuc length, J5' P-nuc length, and
     * NP1 length. The first four are fixed at this point; only NP1 is
     * still being sampled. We pre-compute the deterministic part of
     * the junction modulo 3 and tell `sample_np_length` to draw NP1
     * only from lengths that bring the total to a multiple of 3.
     *
     * This eliminates frame-induced rearrangement retries on light
     * chains entirely (modulo stop codons in the junction, which are
     * rarer). The non-PRODUCTIVE_ONLY path passes target_mod_3=-1
     * and behaves exactly as before — empirical NP length distribution
     * is preserved untouched.
     *
     * VDJ chains stay on the legacy unconstrained path for now (NP1+
     * NP2+D interplay is the next step). */
    int np1_target_mod_3 = -1;
    int j5_p_len_pre = 0;          /* J5' P-nuc length, pre-sampled when
                                       constraining NP1 (see below). */
    bool j5_p_pre_sampled = false;
    if (!has_d &&
        simcfg_productive_only(cfg) &&
        rec->v_allele && allele_has_anchor(rec->v_allele) &&
        rec->j_allele && allele_has_anchor(rec->j_allele)) {

        /* V contribution to junction: bases from V Cys to end of
         * trimmed V. Cys codon starts at v_anchor in original allele
         * coordinates and the trimmed V ends at length - v_trim_3. */
        int v_part = (rec->v_allele->length - rec->v_trim_3)
                     - rec->v_allele->anchor;
        if (v_part < 0) v_part = 0;

        /* J contribution: from start of trimmed J to end of W/F codon
         * (anchor + 3 bases inclusive, exclusive upper bound). */
        int j_part = (rec->j_allele->anchor - rec->j_trim_5) + 3;
        if (j_part < 0) j_part = 0;

        /* V3' P-nuc length is determined here too — replicate the
         * "v_trim_3 == 0" gate of the actual append below. The append
         * reads `v_p_len` later; we set it now and skip the sample
         * call inside the V branch. */
        if (rec->v_trim_3 == 0 && (rec->v_allele->length - rec->v_trim_5) > 0) {
            v_p_len = sample_p_nuc_length(cfg->rng, p_dist);
        }

        /* J5' P-nuc is normally sampled in the J branch AFTER NP1.
         * Pre-sample it here so we know its contribution to the
         * fixed junction modulo before we draw NP1. */
        if (rec->j_trim_5 == 0 &&
            (rec->j_allele->length - rec->j_trim_5 - rec->j_trim_3) > 0) {
            j5_p_len_pre = sample_p_nuc_length(cfg->rng, p_dist);
        }
        j5_p_pre_sampled = true;

        int fixed_mod = (v_part + v_p_len + j5_p_len_pre + j_part) % 3;
        np1_target_mod_3 = (3 - fixed_mod) % 3;

        TRACE("[assemble] VJ frame-aware sampling: v_part=%d v3p=%d "
              "j5p=%d j_part=%d fixed_mod=%d => np1_target_mod_3=%d",
              v_part, v_p_len, j5_p_len_pre, j_part, fixed_mod,
              np1_target_mod_3);
    }

    /* ── Append V segment (trimmed) ─────────────────────────── */
    if (rec->v_allele) {
        const char *v_seq = rec->v_allele->seq;
        int v_start = rec->v_trim_5;
        int v_end   = rec->v_allele->length - rec->v_trim_3;
        int v_len   = v_end - v_start;

        if (v_len > 0) {
            /* Normalize the anchor sentinel: the helper rejects both
             * the GDC sentinel (-1) and the legacy embedded-data
             * sentinel (0) so anchorless V alleles never tag the V
             * first base as a Cys anchor. */
            int v_anchor = allele_has_anchor(rec->v_allele)
                           ? rec->v_allele->anchor : -1;
            aseq_append_segment(seq, v_seq + v_start, v_len,
                                SEG_V, v_start, v_anchor);

            /* V 3' P-nuc: only when v_trim_3 == 0 (untrimmed end).
             * In the VJ frame-aware path the length is pre-sampled
             * above so we don't double-draw RNG; just use it. */
            if (rec->v_trim_3 == 0) {
                int k = (np1_target_mod_3 >= 0)
                        ? v_p_len
                        : sample_p_nuc_length(cfg->rng, p_dist);
                v_p_len = append_p_nuc(seq, v_seq + v_start, v_len, k,
                                       SEG_NP1, /*head_side=*/false);
            }
        }
        TRACE("[assemble] V: %dbp (germline[%d:%d], anchor_at=%d) "
              "+ %dbp V3' P-nuc",
              v_len, v_start, v_end, (int)rec->v_allele->anchor, v_p_len);
    }

    /* ── Generate and append NP1 (TdT N-nucleotides) ───────── */
    {
        char np1_buf[GENAIRR_NP_BUF_SIZE];
        int np1_n_len = generate_np(cfg->rng, np1_buf,
                                    GENAIRR_NP_BUF_SIZE, &cfg->np[0],
                                    np1_target_mod_3);
        if (np1_n_len < 0) {
            /* Frame constraint unsatisfiable for this V/J pair; fall
             * back to unconstrained sampling. The retry loop will
             * still catch any non-productive output. */
            TRACE("[assemble] NP1 frame constraint mod=%d unsatisfiable, "
                  "falling back to unconstrained", np1_target_mod_3);
            np1_n_len = generate_np(cfg->rng, np1_buf,
                                    GENAIRR_NP_BUF_SIZE, &cfg->np[0], -1);
        }
        if (np1_n_len > 0) {
            aseq_append_np(seq, np1_buf, np1_n_len,
                           SEG_NP1, NUC_FLAG_N_NUCLEOTIDE);
        }
        TRACE("[assemble] NP1 N-nucs: %dbp (%s)",
              np1_n_len, np1_n_len > 0 ? np1_buf : "");
    }

    /* ── Append D segment (trimmed) — VDJ only ─────────────── */
    if (has_d) {
        const char *d_seq = rec->d_allele->seq;
        int d_start = rec->d_trim_5;
        int d_end   = rec->d_allele->length - rec->d_trim_3;
        int d_len   = d_end - d_start;

        if (d_len > 0) {
            /* D 5' P-nuc: only when d_trim_5 == 0 and D contributes
             * non-zero length. Goes BEFORE the D segment, tagged as
             * NP1 (between V and D). */
            if (rec->d_trim_5 == 0) {
                int k = sample_p_nuc_length(cfg->rng, p_dist);
                d5_p_len = append_p_nuc(seq, d_seq + d_start, d_len, k,
                                        SEG_NP1, /*head_side=*/true);
            }

            aseq_append_segment(seq, d_seq + d_start, d_len,
                                SEG_D, d_start, -1);

            /* D 3' P-nuc: only when d_trim_3 == 0. Tail-side, tagged
             * as NP2 (between D and J). */
            if (rec->d_trim_3 == 0) {
                int k = sample_p_nuc_length(cfg->rng, p_dist);
                d3_p_len = append_p_nuc(seq, d_seq + d_start, d_len, k,
                                        SEG_NP2, /*head_side=*/false);
            }
        }
        TRACE("[assemble] D: %dbp + D5' P-nuc=%dbp + D3' P-nuc=%dbp",
              d_len, d5_p_len, d3_p_len);

        /* ── Generate and append NP2 (TdT N-nucleotides) ───── */
        char np2_buf[GENAIRR_NP_BUF_SIZE];
        /* VDJ NP2 stays on the legacy unconstrained path; the joint
         * NP1+NP2 frame-aware path is a follow-up step. */
        int np2_n_len = generate_np(cfg->rng, np2_buf,
                                    GENAIRR_NP_BUF_SIZE, &cfg->np[1], -1);
        if (np2_n_len > 0) {
            aseq_append_np(seq, np2_buf, np2_n_len,
                           SEG_NP2, NUC_FLAG_N_NUCLEOTIDE);
        }
        TRACE("[assemble] NP2 N-nucs: %dbp (%s)",
              np2_n_len, np2_n_len > 0 ? np2_buf : "");
    }

    /* ── Append J segment (trimmed) ─────────────────────────── */
    if (rec->j_allele) {
        const char *j_seq = rec->j_allele->seq;
        int j_start = rec->j_trim_5;
        int j_len   = rec->j_allele->length - rec->j_trim_5 - rec->j_trim_3;

        if (j_len > 0) {
            /* J 5' P-nuc: only when j_trim_5 == 0. Goes BEFORE the J
             * segment, tagged with the segment between D and J (NP2
             * for VDJ, NP1 for VJ chains).
             *
             * V5 step 29: when the VJ frame-aware path pre-sampled
             * the J5'P length above, reuse it instead of drawing a
             * second time so RNG consumption stays consistent. */
            if (rec->j_trim_5 == 0) {
                int k = j5_p_pre_sampled
                        ? j5_p_len_pre
                        : sample_p_nuc_length(cfg->rng, p_dist);
                j_p_len = append_p_nuc(seq, j_seq + j_start, j_len, k,
                                       j5_np_seg, /*head_side=*/true);
            }

            /* Normalize anchor sentinel — see step_assemble's V branch. */
            int j_anchor = allele_has_anchor(rec->j_allele)
                           ? rec->j_allele->anchor : -1;
            aseq_append_segment(seq, j_seq + j_start, j_len,
                                SEG_J, j_start, j_anchor);
        }
        TRACE("[assemble] J: %dbp (germline[%d:%d], anchor_at=%d) "
              "+ %dbp J5' P-nuc",
              j_len, j_start, j_start + j_len, (int)rec->j_allele->anchor,
              j_p_len);
    }

    /* Record total NP1/NP2 lengths (P-nucs + N-nucs) so the trace and
     * AIRR np1_length/np2_length both report the full segment span. */
    rec->np1_length = aseq_segment_length(seq, SEG_NP1);
    rec->np2_length = has_d ? aseq_segment_length(seq, SEG_NP2) : 0;

    /* ── Junction marking ─────────────────────────────────────── */
    /* Tag every node in [V Cys, J W/F + 3) with NUC_FLAG_JUNCTION so
     * the junction can be tracked through the rest of the pipeline by
     * node flags rather than re-derived after each step. Capture the
     * original length so airr_serialize can detect truncation. */
    rec->original_junction_length = aseq_mark_junction(seq);

    TRACE("[assemble] total: %dbp = V(%d) + NP1(%d, P=%d) + D(%d) "
          "+ NP2(%d, P=%d) + J(%d)",
          seq->length,
          aseq_segment_length(seq, SEG_V),
          rec->np1_length, v_p_len + d5_p_len,
          has_d ? aseq_segment_length(seq, SEG_D) : 0,
          rec->np2_length, d3_p_len + (has_d ? j_p_len : 0),
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
