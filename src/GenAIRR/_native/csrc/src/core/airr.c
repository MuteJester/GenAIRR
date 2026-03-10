/**
 * airr.c — AIRR-format serialization from ASeq.
 *
 * Derives all AIRR positional fields (v_sequence_start, d_sequence_end, etc.)
 * directly from the ASeq linked list. Each nucleotide knows its segment tag
 * and germline position, so positions are always consistent.
 *
 * Boundary ambiguity check:
 *   When V is trimmed, the removed bases might be identical to the
 *   start of NP1 (by random chance or biological homology). An aligner
 *   seeing the flat sequence would extend V into NP1. We replicate this
 *   by checking if V's trimmed 3' bases match the adjacent NP1 bases.
 */

#include "genairr/airr.h"
#include "genairr/codon.h"
#include "genairr/trace.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

/* MSVC names POSIX strdup as _strdup */
#ifdef _MSC_VER
  #define strdup _strdup
#endif

/* ── Internal: count nodes from head to a target node ─────────── */

static int position_of_seg_boundary(const ASeq *seq, const Nuc *node) {
    if (!node) return -1;
    int pos = 0;
    for (Nuc *n = seq->head; n; n = n->next) {
        if (n == node) return pos;
        pos++;
    }
    return -1;
}

/* ── Internal: check V 3' boundary ambiguity ──────────────────── */

/**
 * Check how many bases of V's trimmed 3' end match the start of NP1.
 *
 * This replicates FixVPosition: walking the trimmed V reference
 * and comparing to the NP1 bases that follow.
 *
 * Returns the number of matching bases (0 if no ambiguity).
 */
static int check_v_3prime_ambiguity(const SimRecord *rec, const ASeq *seq) {
    if (!rec->v_allele || rec->v_trim_3 == 0) return 0;
    if (!seq->seg_first[SEG_NP1]) return 0;

    const char *v_ref = rec->v_allele->seq;
    int v_len = rec->v_allele->length;
    int trim_3 = rec->v_trim_3;

    /* The trimmed V 3' bases: v_ref[v_len - trim_3 .. v_len) */
    const char *trimmed = v_ref + (v_len - trim_3);

    /* Walk NP1 nodes and compare */
    int match = 0;
    Nuc *np = seq->seg_first[SEG_NP1];
    for (int i = 0; i < trim_3 && np && np->segment == SEG_NP1; i++, np = np->next) {
        if (trimmed[i] == np->current) {
            match++;
        } else {
            break;
        }
    }
    return match;
}

/* ── Internal: check J 5' boundary ambiguity ──────────────────── */

/**
 * Check how many bases of J's trimmed 5' end match the end of NP
 * (NP2 for VDJ chains, NP1 for VJ chains).
 *
 * This replicates FixJPosition: walking the trimmed J reference
 * backwards and comparing to the NP bases that precede J.
 */
static int check_j_5prime_ambiguity(const SimRecord *rec, const ASeq *seq) {
    if (!rec->j_allele || rec->j_trim_5 == 0) return 0;

    /* Find which NP segment precedes J */
    Segment np_seg = seq->seg_first[SEG_NP2] ? SEG_NP2 : SEG_NP1;
    if (!seq->seg_last[np_seg]) return 0;

    const char *j_ref = rec->j_allele->seq;
    int trim_5 = rec->j_trim_5;

    /* The trimmed J 5' bases: j_ref[0 .. trim_5), compared backwards */
    int match = 0;
    Nuc *np = seq->seg_last[np_seg];
    for (int i = trim_5 - 1; i >= 0 && np && np->segment == np_seg; i--, np = np->prev) {
        if (j_ref[i] == np->current) {
            match++;
        } else {
            break;
        }
    }
    return match;
}

/* ── Internal: check D 5' and 3' boundary ambiguity ───────────── */

static int check_d_5prime_ambiguity(const SimRecord *rec, const ASeq *seq) {
    if (!rec->d_allele || rec->d_trim_5 == 0) return 0;
    if (!seq->seg_last[SEG_NP1]) return 0;

    const char *d_ref = rec->d_allele->seq;
    int trim_5 = rec->d_trim_5;

    /* D's trimmed 5' bases: d_ref[0 .. trim_5), compared backwards against NP1 end */
    int match = 0;
    Nuc *np = seq->seg_last[SEG_NP1];
    for (int i = trim_5 - 1; i >= 0 && np && np->segment == SEG_NP1; i--, np = np->prev) {
        if (d_ref[i] == np->current) {
            match++;
        } else {
            break;
        }
    }
    return match;
}

static int check_d_3prime_ambiguity(const SimRecord *rec, const ASeq *seq) {
    if (!rec->d_allele || rec->d_trim_3 == 0) return 0;
    if (!seq->seg_first[SEG_NP2]) return 0;

    const char *d_ref = rec->d_allele->seq;
    int d_len = rec->d_allele->length;
    int trim_3 = rec->d_trim_3;

    /* D's trimmed 3' bases: d_ref[d_len - trim_3 .. d_len) */
    const char *trimmed = d_ref + (d_len - trim_3);

    int match = 0;
    Nuc *np = seq->seg_first[SEG_NP2];
    for (int i = 0; i < trim_3 && np && np->segment == SEG_NP2; i++, np = np->next) {
        if (trimmed[i] == np->current) {
            match++;
        } else {
            break;
        }
    }
    return match;
}

/* ── Main derivation function ─────────────────────────────────── */

void airr_derive_positions(const ASeq *seq, const SimRecord *rec,
                           AirrPositions *out) {
    memset(out, 0, sizeof(*out));

    bool has_d = rec->d_allele != NULL && seq->seg_first[SEG_D] != NULL;

    /* ── Ground truth positions (from segment boundaries) ────── */

    /* V positions */
    if (seq->seg_first[SEG_V]) {
        out->v_sequence_start = position_of_seg_boundary(seq, seq->seg_first[SEG_V]);
        out->v_sequence_end   = position_of_seg_boundary(seq, seq->seg_last[SEG_V]) + 1;
        uint16_t vg0 = seq->seg_first[SEG_V]->germline_pos;
        uint16_t vg1 = seq->seg_last[SEG_V]->germline_pos;
        out->v_germline_start = vg0 < vg1 ? vg0 : vg1;
        out->v_germline_end   = (vg0 > vg1 ? vg0 : vg1) + 1;
    }

    /* D positions */
    if (has_d) {
        out->d_sequence_start = position_of_seg_boundary(seq, seq->seg_first[SEG_D]);
        out->d_sequence_end   = position_of_seg_boundary(seq, seq->seg_last[SEG_D]) + 1;
        uint16_t dg0 = seq->seg_first[SEG_D]->germline_pos;
        uint16_t dg1 = seq->seg_last[SEG_D]->germline_pos;
        out->d_germline_start = dg0 < dg1 ? dg0 : dg1;
        out->d_germline_end   = (dg0 > dg1 ? dg0 : dg1) + 1;
    }

    /* J positions */
    if (seq->seg_first[SEG_J]) {
        out->j_sequence_start = position_of_seg_boundary(seq, seq->seg_first[SEG_J]);
        out->j_sequence_end   = position_of_seg_boundary(seq, seq->seg_last[SEG_J]) + 1;
        uint16_t jg0 = seq->seg_first[SEG_J]->germline_pos;
        uint16_t jg1 = seq->seg_last[SEG_J]->germline_pos;
        out->j_germline_start = jg0 < jg1 ? jg0 : jg1;
        out->j_germline_end   = (jg0 > jg1 ? jg0 : jg1) + 1;
    }

    /* ── Boundary ambiguity adjustments ──────────────────────── */
    /* These replicate FixV/D/JPosition: extend segment boundaries
     * where the trimmed germline bases match the adjacent NP region.
     * This is what an external aligner would report. */

    TRACE("[airr] ground truth: V=[%d:%d] germ[%d:%d], D=[%d:%d] germ[%d:%d], J=[%d:%d] germ[%d:%d]",
          out->v_sequence_start, out->v_sequence_end, out->v_germline_start, out->v_germline_end,
          out->d_sequence_start, out->d_sequence_end, out->d_germline_start, out->d_germline_end,
          out->j_sequence_start, out->j_sequence_end, out->j_germline_start, out->j_germline_end);

    int v_ambig = check_v_3prime_ambiguity(rec, seq);
    out->v_sequence_end   += v_ambig;
    out->v_germline_end   += v_ambig;
    out->v_trim_3_adjusted = rec->v_trim_3 - v_ambig;

    int j_ambig = check_j_5prime_ambiguity(rec, seq);
    out->j_sequence_start -= j_ambig;
    out->j_germline_start -= j_ambig;
    out->j_trim_5_adjusted = rec->j_trim_5 - j_ambig;

    if (has_d) {
        int d5_ambig = check_d_5prime_ambiguity(rec, seq);
        out->d_sequence_start -= d5_ambig;
        out->d_germline_start -= d5_ambig;
        out->d_trim_5_adjusted = rec->d_trim_5 - d5_ambig;

        int d3_ambig = check_d_3prime_ambiguity(rec, seq);
        out->d_sequence_end   += d3_ambig;
        out->d_germline_end   += d3_ambig;
        out->d_trim_3_adjusted = rec->d_trim_3 - d3_ambig;

        if (d5_ambig || d3_ambig) {
            TRACE("[airr] D boundary ambiguity: 5'=%d, 3'=%d", d5_ambig, d3_ambig);
        }
    }

    if (v_ambig || j_ambig) {
        TRACE("[airr] boundary ambiguity: V_3'=%d bases, J_5'=%d bases", v_ambig, j_ambig);
    }

    /* ── Junction ────────────────────────────────────────────── */
    /* IMGT convention: junction runs from conserved C (V anchor) to the
     * end of the conserved W/F codon (J anchor + 3 bases).
     * Must match functionality.c which uses junction_end = j_anchor_pos + 3. */
    Nuc *v_anchor = aseq_find_anchor(seq, SEG_V);
    Nuc *j_anchor = aseq_find_anchor(seq, SEG_J);
    if (v_anchor && j_anchor) {
        int seq_len = aseq_length(seq);
        out->junction_start = position_of_seg_boundary(seq, v_anchor);
        out->junction_end   = position_of_seg_boundary(seq, j_anchor) + 3;
        /* Cap to sequence length: 3' corruption can remove bases after
         * the J anchor, making junction_end exceed the actual length. */
        if (out->junction_end > seq_len)
            out->junction_end = seq_len;
        out->junction_length = out->junction_end - out->junction_start;
        TRACE("[airr] junction: [%d:%d] (%dbp)", out->junction_start, out->junction_end, out->junction_length);
    }
}

/* Correction maps replaced by AlleleBitmap (see allele_bitmap.h). */

#if 0

static bool seq_contains(const char *haystack, int haystack_len,
                          const char *needle, int needle_len) {
    if (needle_len == 0) return true;
    if (needle_len > haystack_len) return false;
    for (int i = 0; i <= haystack_len - needle_len; i++) {
        bool match = true;
        for (int j = 0; j < needle_len; j++) {
            if (haystack[i + j] != needle[j]) {
                match = false;
                break;
            }
        }
        if (match) return true;
    }
    return false;
}

/* ── Internal: find an allele index by name in a correction map ── */

static int find_allele_index(char **names, int count, const char *name) {
    for (int i = 0; i < count; i++) {
        if (strcmp(names[i], name) == 0) return i;
    }
    return -1;
}

/* ── build_v_end_correction_map (3' trim) ─────────────────────── */

TrimCorrectionMap build_v_end_correction_map(const AllelePool *pool) {
    TrimCorrectionMap map = {0};
    map.count = pool->count;
    map.names   = calloc(pool->count, sizeof(char *));
    map.entries = calloc(pool->count, sizeof(TrimCorrectionEntry));

    for (int i = 0; i < pool->count; i++) {
        map.names[i] = strdup(pool->alleles[i].name);
    }

    for (int i = 0; i < pool->count; i++) {
        const Allele *a = &pool->alleles[i];
        int seq_len = a->length;

        TrimCorrectionEntry *e = &map.entries[i];
        e->max_trim    = seq_len;
        e->max_alleles = pool->count;
        e->equiv_names  = calloc((seq_len + 1) * pool->count, sizeof(char *));
        e->equiv_counts = calloc(seq_len + 1, sizeof(int));

        for (int trim_3 = 0; trim_3 <= seq_len; trim_3++) {
            /* Trimmed subsequence: a->seq[0 .. seq_len - trim_3) */
            int trimmed_len = seq_len - trim_3;
            int count = 0;

            for (int j = 0; j < pool->count; j++) {
                const Allele *other = &pool->alleles[j];
                if (seq_contains(other->seq, other->length,
                                 a->seq, trimmed_len)) {
                    e->equiv_names[trim_3 * pool->count + count] =
                        map.names[j];
                    count++;
                }
            }
            e->equiv_counts[trim_3] = count;
        }
    }

    return map;
}

/* ── build_j_trim_correction_map (5' trim) ────────────────────── */

TrimCorrectionMap build_j_trim_correction_map(const AllelePool *pool) {
    TrimCorrectionMap map = {0};
    map.count = pool->count;
    map.names   = calloc(pool->count, sizeof(char *));
    map.entries = calloc(pool->count, sizeof(TrimCorrectionEntry));

    for (int i = 0; i < pool->count; i++) {
        map.names[i] = strdup(pool->alleles[i].name);
    }

    for (int i = 0; i < pool->count; i++) {
        const Allele *a = &pool->alleles[i];
        int seq_len = a->length;

        TrimCorrectionEntry *e = &map.entries[i];
        e->max_trim    = seq_len;
        e->max_alleles = pool->count;
        e->equiv_names  = calloc((seq_len + 1) * pool->count, sizeof(char *));
        e->equiv_counts = calloc(seq_len + 1, sizeof(int));

        for (int trim_5 = 0; trim_5 <= seq_len; trim_5++) {
            /* Trimmed subsequence: a->seq[trim_5 .. seq_len) */
            const char *trimmed = a->seq + trim_5;
            int trimmed_len = seq_len - trim_5;
            int count = 0;

            for (int j = 0; j < pool->count; j++) {
                const Allele *other = &pool->alleles[j];
                if (seq_contains(other->seq, other->length,
                                 trimmed, trimmed_len)) {
                    e->equiv_names[trim_5 * pool->count + count] =
                        map.names[j];
                    count++;
                }
            }
            e->equiv_counts[trim_5] = count;
        }
    }

    return map;
}

/* ── build_d_trim_correction_map (5' + 3' trim) ──────────────── */

DTrimCorrectionMap build_d_trim_correction_map(const AllelePool *pool) {
    DTrimCorrectionMap map = {0};
    map.count = pool->count;
    map.names   = calloc(pool->count, sizeof(char *));
    map.entries = calloc(pool->count, sizeof(DTrimCorrectionEntry));

    for (int i = 0; i < pool->count; i++) {
        map.names[i] = strdup(pool->alleles[i].name);
    }

    for (int i = 0; i < pool->count; i++) {
        const Allele *a = &pool->alleles[i];
        int seq_len = a->length;

        DTrimCorrectionEntry *e = &map.entries[i];
        e->max_trim_5  = seq_len;
        e->max_trim_3  = seq_len;
        e->max_alleles = pool->count;

        int dim_3 = seq_len + 1;  /* stride for trim_3 dimension */
        e->equiv_names  = calloc((seq_len + 1) * dim_3 * pool->count,
                                  sizeof(char *));
        e->equiv_counts = calloc((seq_len + 1) * dim_3, sizeof(int));

        for (int t5 = 0; t5 <= seq_len; t5++) {
            for (int t3 = 0; t3 <= seq_len - t5; t3++) {
                /* Trimmed subsequence: a->seq[t5 .. seq_len - t3) */
                const char *trimmed = a->seq + t5;
                int trimmed_len = seq_len - t5 - t3;
                int flat_2d = t5 * dim_3 + t3;
                int count = 0;

                for (int j = 0; j < pool->count; j++) {
                    const Allele *other = &pool->alleles[j];
                    if (seq_contains(other->seq, other->length,
                                     trimmed, trimmed_len)) {
                        e->equiv_names[flat_2d * pool->count + count] =
                            map.names[j];
                        count++;
                    }
                }
                e->equiv_counts[flat_2d] = count;
            }
        }
    }

    return map;
}

/* ── Free correction maps ─────────────────────────────────────── */

void trim_correction_map_destroy(TrimCorrectionMap *map) {
    if (!map) return;
    for (int i = 0; i < map->count; i++) {
        free(map->names[i]);
        free(map->entries[i].equiv_names);
        free(map->entries[i].equiv_counts);
    }
    free(map->names);
    free(map->entries);
    *map = (TrimCorrectionMap){0};
}

void dtrim_correction_map_destroy(DTrimCorrectionMap *map) {
    if (!map) return;
    for (int i = 0; i < map->count; i++) {
        free(map->names[i]);
        free(map->entries[i].equiv_names);
        free(map->entries[i].equiv_counts);
    }
    free(map->names);
    free(map->entries);
    *map = (DTrimCorrectionMap){0};
}

/* ═══════════════════════════════════════════════════════════════
 * Correction query functions
 *
 * These take a correction map + the allele name + trim amounts,
 * and return the list of equivalent allele names.
 * ═══════════════════════════════════════════════════════════════ */

int correct_v_end_cut(const TrimCorrectionMap *map,
                       const char *allele_name, int trim_3,
                       const char **out_names, int max_out) {
    int idx = find_allele_index(map->names, map->count, allele_name);
    if (idx < 0) return 0;

    const TrimCorrectionEntry *e = &map->entries[idx];
    int valid_trim = trim_3 < e->max_trim ? trim_3 : e->max_trim;
    int n = e->equiv_counts[valid_trim];
    int out_count = 0;

    for (int i = 0; i < n && out_count < max_out; i++) {
        out_names[out_count++] = e->equiv_names[valid_trim * e->max_alleles + i];
    }
    return out_count;
}

int correct_j_trims(const TrimCorrectionMap *map,
                     const char *allele_name, int trim_5,
                     const char **out_names, int max_out) {
    int idx = find_allele_index(map->names, map->count, allele_name);
    if (idx < 0) return 0;

    const TrimCorrectionEntry *e = &map->entries[idx];
    int valid_trim = trim_5 < e->max_trim ? trim_5 : e->max_trim;
    int n = e->equiv_counts[valid_trim];
    int out_count = 0;

    for (int i = 0; i < n && out_count < max_out; i++) {
        out_names[out_count++] = e->equiv_names[valid_trim * e->max_alleles + i];
    }
    return out_count;
}

int correct_d_trims(const DTrimCorrectionMap *map,
                     const char *allele_name, int trim_5, int trim_3,
                     const char **out_names, int max_out) {
    int idx = find_allele_index(map->names, map->count, allele_name);
    if (idx < 0) return 0;

    const DTrimCorrectionEntry *e = &map->entries[idx];
    if (trim_5 > e->max_trim_5) trim_5 = e->max_trim_5;
    if (trim_3 > e->max_trim_3) trim_3 = e->max_trim_3;

    int dim_3 = e->max_trim_3 + 1;
    int flat_2d = trim_5 * dim_3 + trim_3;
    int n = e->equiv_counts[flat_2d];
    int out_count = 0;

    for (int i = 0; i < n && out_count < max_out; i++) {
        out_names[out_count++] = e->equiv_names[flat_2d * e->max_alleles + i];
    }
    return out_count;
}

/* ═══════════════════════════════════════════════════════════════
 * Short-D threshold computation
 *
 * Mirrors Python compute_short_d_threshold(): for increasing L=1,2,...
 * checks if any substring of length L from a D allele is unique
 * (not found in any other D allele). Returns the first such L.
 * ═══════════════════════════════════════════════════════════════ */

int compute_short_d_threshold(const AllelePool *pool) {
    if (pool->count <= 1) return 1;

    /* Deduplicate allele sequences (Python uses a set comprehension).
     * Two alleles with identical sequences should count as one.
     * Without dedup, all-identical alleles would return shortest+1
     * instead of 1. */
    const char *unique_seqs[GENAIRR_MAX_ALLELES];
    int unique_lens[GENAIRR_MAX_ALLELES];
    int n_unique = 0;

    for (int i = 0; i < pool->count; i++) {
        bool dup = false;
        for (int j = 0; j < n_unique; j++) {
            if (pool->alleles[i].length == unique_lens[j] &&
                memcmp(pool->alleles[i].seq, unique_seqs[j],
                       pool->alleles[i].length) == 0) {
                dup = true;
                break;
            }
        }
        if (!dup && n_unique < GENAIRR_MAX_ALLELES) {
            unique_seqs[n_unique] = pool->alleles[i].seq;
            unique_lens[n_unique] = pool->alleles[i].length;
            n_unique++;
        }
    }

    if (n_unique <= 1) return 1;

    /* Find shortest unique allele length */
    int shortest = unique_lens[0];
    for (int i = 1; i < n_unique; i++) {
        if (unique_lens[i] < shortest)
            shortest = unique_lens[i];
    }

    for (int len = 1; len <= shortest; len++) {
        /* For each unique allele, extract all substrings of this length.
         * Check if any substring appears in only one allele. */
        for (int i = 0; i < n_unique; i++) {
            for (int s = 0; s <= unique_lens[i] - len; s++) {
                const char *sub = unique_seqs[i] + s;
                bool unique = true;

                for (int j = 0; j < n_unique && unique; j++) {
                    if (j == i) continue;
                    if (seq_contains(unique_seqs[j], unique_lens[j],
                                     sub, len)) {
                        unique = false;
                    }
                }

                if (unique) return len;
            }
        }
    }

    return shortest + 1;
}

/* ═══════════════════════════════════════════════════════════════
 * ReassessAlleleCalls — pairwise difference table approach
 *
 * Mirrors Python ReassessAlleleCalls._reassess_segment_fast().
 *
 * At config time: for every pair of alleles (A, B), record the
 * germline positions where they differ. This is O(N²·L) but done
 * once and cached.
 *
 * At runtime: score the true allele once (O(L)), then for each
 * candidate, compute delta using only distinguishing positions
 * (O(D) where D << L for same-gene alleles).
 * ═══════════════════════════════════════════════════════════════ */

PairwiseDiffTable build_pairwise_diff_table(const AllelePool *pool) {
    PairwiseDiffTable table = {0};
    table.count   = pool->count;
    table.names   = calloc(pool->count, sizeof(char *));
    table.lengths = calloc(pool->count, sizeof(int));
    table.entries = calloc(pool->count, sizeof(PairwiseDiffEntry));

    for (int i = 0; i < pool->count; i++) {
        table.names[i]   = strdup(pool->alleles[i].name);
        table.lengths[i] = pool->alleles[i].length;
    }

    for (int a = 0; a < pool->count; a++) {
        PairwiseDiffEntry *e = &table.entries[a];
        e->diffs   = calloc(pool->count, sizeof(AlleleDiff *));
        e->n_diffs = calloc(pool->count, sizeof(int));

        const char *a_seq = pool->alleles[a].seq;
        int a_len = pool->alleles[a].length;

        for (int b = 0; b < pool->count; b++) {
            if (b == a) continue;

            const char *b_seq = pool->alleles[b].seq;
            int b_len = pool->alleles[b].length;
            int min_len = a_len < b_len ? a_len : b_len;

            /* First pass: count differences */
            int n = 0;
            for (int p = 0; p < min_len; p++) {
                if (a_seq[p] != b_seq[p]) n++;
            }

            if (n == 0) continue;

            /* Second pass: record them */
            e->diffs[b] = calloc(n, sizeof(AlleleDiff));
            e->n_diffs[b] = n;
            int idx = 0;
            for (int p = 0; p < min_len; p++) {
                if (a_seq[p] != b_seq[p]) {
                    e->diffs[b][idx].germ_pos = p;
                    e->diffs[b][idx].a_base   = a_seq[p];
                    e->diffs[b][idx].b_base   = b_seq[p];
                    idx++;
                }
            }
        }
    }

    return table;
}

void pairwise_diff_table_destroy(PairwiseDiffTable *table) {
    if (!table) return;
    for (int i = 0; i < table->count; i++) {
        free(table->names[i]);
        if (table->entries[i].diffs) {
            for (int j = 0; j < table->count; j++) {
                free(table->entries[i].diffs[j]);
            }
            free(table->entries[i].diffs);
        }
        free(table->entries[i].n_diffs);
    }
    free(table->names);
    free(table->lengths);
    free(table->entries);
    *table = (PairwiseDiffTable){0};
}

int reassess_allele_calls(const PairwiseDiffTable *table,
                           const AllelePool *pool,
                           const char *seq,
                           const char *true_allele,
                           int seq_start, int seq_end,
                           int germ_start, int germ_end,
                           const char **out_names, int max_out) {
    /* Guard: empty or invalid region */
    int expected_len = seq_end - seq_start;
    int germ_len = germ_end - germ_start;
    if (expected_len <= 0 || germ_len <= 0) return 0;
    if (expected_len != germ_len) return 0;  /* indels → no 1:1 alignment */

    /* Find true allele index in table (table built from pool, same order) */
    int true_idx = find_allele_index(table->names, table->count, true_allele);
    if (true_idx < 0) return 0;

    /* Get the true allele's germline sequence from pool */
    const char *true_seq = pool->alleles[true_idx].seq;
    if (table->lengths[true_idx] < germ_end) return 0;

    int true_score = 0;
    int offset = germ_start - seq_start;  /* germ_pos = seq_pos + offset */

    /* Precompute: which germline positions have obs matching true allele.
     * Use a bitset for fast lookup (max 1024 positions). */
    bool obs_matches_true[GENAIRR_MAX_SEQ_LEN] = {false};

    for (int i = 0; i < expected_len; i++) {
        int gp = germ_start + i;
        char obs = seq[seq_start + i];
        if (obs == 'N' || obs == 'n') continue;
        if (obs == true_seq[gp]) {
            obs_matches_true[gp] = true;
        } else {
            true_score++;
        }
    }

    /* Delta scoring for each candidate */
    int out_count = 0;
    const PairwiseDiffEntry *true_entry = &table->entries[true_idx];

    for (int c = 0; c < table->count && out_count < max_out; c++) {
        if (c == true_idx) continue;

        /* Length check: candidate must cover [0, germ_end) */
        if (table->lengths[c] < germ_end) continue;

        int n_diffs = true_entry->n_diffs[c];
        if (n_diffs == 0) {
            /* Identical alleles → same score as true → include */
            out_names[out_count++] = table->names[c];
            continue;
        }

        const AlleleDiff *diffs = true_entry->diffs[c];
        int delta = 0;
        int remaining = n_diffs;

        for (int d = 0; d < n_diffs; d++) {
            remaining--;
            int gp = diffs[d].germ_pos;

            if (gp < germ_start || gp >= germ_end) continue;

            char obs = seq[gp - offset];
            if (obs == 'N' || obs == 'n') continue;

            if (obs_matches_true[gp]) {
                /* obs matches true but not candidate → candidate gets +1 */
                delta++;
                /* Early exit: if delta exceeds what remaining diffs
                 * could recover (each can subtract at most 1) */
                if (delta > remaining) break;
            } else if (obs == diffs[d].b_base) {
                /* obs matches candidate but not true → candidate gets -1 */
                delta--;
            }
            /* else: obs matches neither → both get 1, delta unchanged */
        }

        /* Include if candidate scores <= true */
        if (delta <= 0) {
            out_names[out_count++] = table->names[c];
        }
    }

    return out_count;
}

/* ═══════════════════════════════════════════════════════════════
 * CorrectionContext — build and destroy
 * ═══════════════════════════════════════════════════════════════ */

CorrectionContext correction_context_build(const SimConfig *cfg) {
    CorrectionContext ctx = {0};

    if (cfg->v_alleles.count > 0) {
        ctx.v_correction = build_v_end_correction_map(&cfg->v_alleles);
        ctx.v_diff_table = build_pairwise_diff_table(&cfg->v_alleles);
    }
    if (cfg->j_alleles.count > 0) {
        ctx.j_correction = build_j_trim_correction_map(&cfg->j_alleles);
        ctx.j_diff_table = build_pairwise_diff_table(&cfg->j_alleles);
    }
    if (cfg->d_alleles.count > 0) {
        ctx.d_correction = build_d_trim_correction_map(&cfg->d_alleles);
        ctx.d_diff_table = build_pairwise_diff_table(&cfg->d_alleles);
        ctx.short_d_threshold = compute_short_d_threshold(&cfg->d_alleles);
    }

    return ctx;
}

void correction_context_destroy(CorrectionContext *ctx) {
    trim_correction_map_destroy(&ctx->v_correction);
    trim_correction_map_destroy(&ctx->j_correction);
    dtrim_correction_map_destroy(&ctx->d_correction);
    pairwise_diff_table_destroy(&ctx->v_diff_table);
    pairwise_diff_table_destroy(&ctx->d_diff_table);
    pairwise_diff_table_destroy(&ctx->j_diff_table);
    memset(ctx, 0, sizeof(*ctx));
}

#endif  /* dead code: old correction maps */

/* ═══════════════════════════════════════════════════════════════
 * Internal: collect mutation/error annotations from ASeq
 *
 * Walks the list once, writing "pos:G>A,pos:C>T,..." for each
 * category of base change.
 * ═══════════════════════════════════════════════════════════════ */

static int append_annotation(char *buf, int pos, int cap,
                              int seq_pos, char germline, char current) {
    if (pos >= cap - 16) return pos;  /* leave room */
    if (pos > 0) buf[pos++] = ',';
    pos += snprintf(buf + pos, cap - pos, "%d:%c>%c",
                    seq_pos, germline, current);
    return pos;
}

static void collect_annotations(const ASeq *seq, AirrRecord *out) {
    int m_pos = 0, se_pos = 0, pcr_pos = 0;
    int n_mut = 0, n_se = 0, n_pcr = 0, n_ins = 0, n_ns = 0;
    int total_germline = 0;
    int i = 0;

    for (Nuc *n = seq->head; n; n = n->next, i++) {
        if (nuc_is_germline_segment(n)) total_germline++;

        if (n->flags & NUC_FLAG_MUTATED) {
            m_pos = append_annotation(out->mutations, m_pos,
                                       AIRR_MAX_MUTATIONS,
                                       i, n->germline, n->current);
            n_mut++;
        }
        if (n->flags & NUC_FLAG_SEQ_ERROR) {
            se_pos = append_annotation(out->sequencing_errors, se_pos,
                                        AIRR_MAX_MUTATIONS,
                                        i, n->germline, n->current);
            n_se++;
        }
        if (n->flags & NUC_FLAG_PCR_ERROR) {
            pcr_pos = append_annotation(out->pcr_errors, pcr_pos,
                                         AIRR_MAX_MUTATIONS,
                                         i, n->germline, n->current);
            n_pcr++;
        }
        if (n->flags & NUC_FLAG_INDEL_INS) n_ins++;
        if (n->flags & NUC_FLAG_IS_N)      n_ns++;
    }

    out->mutations[m_pos] = '\0';
    out->sequencing_errors[se_pos] = '\0';
    out->pcr_errors[pcr_pos] = '\0';

    out->n_mutations        = n_mut;
    out->n_sequencing_errors = n_se;
    out->n_pcr_errors        = n_pcr;
    out->n_insertions        = n_ins;
    out->n_ns                = n_ns;
    out->mutation_rate = total_germline > 0
                         ? (double)n_mut / total_germline : 0.0;
}

/* ── Internal: build germline alignment string ────────────────── */

static void build_germline_alignment(const ASeq *seq, char *buf, int maxlen) {
    int i = 0;
    for (Nuc *n = seq->head; n && i < maxlen - 1; n = n->next, i++) {
        if (nuc_is_germline_segment(n) && n->germline != '\0') {
            buf[i] = n->germline;
        } else {
            /* NP regions, adapter nodes, and indel-inserted nodes
             * (which have segment=V/D/J but germline='\0'): use 'N'. */
            buf[i] = 'N';
        }
    }
    buf[i] = '\0';
}

/* ── Internal: patch germline alignment for boundary extensions ── */

/**
 * When boundary ambiguity extends V/D/J coordinates into NP territory,
 * the germline alignment should show the allele's germline bases at
 * those positions, not 'N'. A god-aligner would reconstruct the
 * germline reference there.
 */
static void patch_boundary_extensions(const SimRecord *rec,
                                       const AirrPositions *pos,
                                       char *germ, int germ_len) {
    /* V 3' extension: NP1 bases claimed as V */
    int v_ambig = rec->v_trim_3 - pos->v_trim_3_adjusted;
    if (v_ambig > 0 && rec->v_allele) {
        int v_len = rec->v_allele->length;
        int patch_start = pos->v_sequence_end - v_ambig;
        int allele_offset = v_len - rec->v_trim_3;
        if (patch_start >= 0 && allele_offset >= 0) {
            for (int i = 0; i < v_ambig; i++) {
                if (patch_start + i >= germ_len) break;
                if (allele_offset + i >= v_len) break;
                germ[patch_start + i] = rec->v_allele->seq[allele_offset + i];
            }
        }
    }

    /* J 5' extension: NP bases claimed as J */
    int j_ambig = rec->j_trim_5 - pos->j_trim_5_adjusted;
    if (j_ambig > 0 && rec->j_allele) {
        int patch_start = pos->j_sequence_start;
        int allele_offset = rec->j_trim_5 - j_ambig;
        int j_len = rec->j_allele->length;
        if (patch_start >= 0 && allele_offset >= 0) {
            for (int i = 0; i < j_ambig; i++) {
                if (patch_start + i >= germ_len) break;
                if (allele_offset + i >= j_len) break;
                germ[patch_start + i] = rec->j_allele->seq[allele_offset + i];
            }
        }
    }

    /* D 5' extension: NP1 bases claimed as D */
    int d5_ambig = rec->d_trim_5 - pos->d_trim_5_adjusted;
    if (d5_ambig > 0 && rec->d_allele) {
        int patch_start = pos->d_sequence_start;
        int allele_offset = rec->d_trim_5 - d5_ambig;
        int d_len = rec->d_allele->length;
        if (patch_start >= 0 && allele_offset >= 0) {
            for (int i = 0; i < d5_ambig; i++) {
                if (patch_start + i >= germ_len) break;
                if (allele_offset + i >= d_len) break;
                germ[patch_start + i] = rec->d_allele->seq[allele_offset + i];
            }
        }
    }

    /* D 3' extension: NP2 bases claimed as D */
    int d3_ambig = rec->d_trim_3 - pos->d_trim_3_adjusted;
    if (d3_ambig > 0 && rec->d_allele) {
        int d_len = rec->d_allele->length;
        int patch_start = pos->d_sequence_end - d3_ambig;
        int allele_offset = d_len - rec->d_trim_3;
        if (patch_start >= 0 && allele_offset >= 0) {
            for (int i = 0; i < d3_ambig; i++) {
                if (patch_start + i >= germ_len) break;
                if (allele_offset + i >= d_len) break;
                germ[patch_start + i] = rec->d_allele->seq[allele_offset + i];
            }
        }
    }
}

/* ── Internal: collect mutations from boundary-extended positions ── */

/**
 * After patching germline_alignment for boundary extension, NP bases
 * absorbed into V/D/J may differ from the extended germline. A god-
 * aligner would count these as mutations. This function scans the
 * extended regions and appends any mismatches to the mutations string.
 *
 * Returns the number of extra mutations found.
 */
static int collect_boundary_mutations(const char *seq_str, const char *germ,
                                       int seq_len,
                                       const AirrPositions *pos,
                                       const SimRecord *rec,
                                       char *mutations_buf, int mut_pos,
                                       int bufsize) {
    int extra = 0;

    /* Helper macro to scan one extended region.
     * Guards against out-of-bounds access: after indels or corruption,
     * boundary extension can produce negative or too-large positions. */
    #define SCAN_EXTENSION(ambig, start)                                    \
        if ((ambig) > 0 && (start) >= 0) {                                 \
            for (int i = 0; i < (ambig); i++) {                            \
                int p = (start) + i;                                        \
                if (p < 0 || p >= seq_len) continue;                       \
                char s = seq_str[p], g = germ[p];                          \
                /* Case-insensitive, skip N */                              \
                char su = (s >= 'a' && s <= 'z') ? s - 32 : s;            \
                char gu = (g >= 'a' && g <= 'z') ? g - 32 : g;            \
                if (su != gu && su != 'N' && gu != 'N') {                  \
                    mut_pos = append_annotation(mutations_buf, mut_pos,     \
                                                 bufsize, p, g, s);        \
                    extra++;                                                \
                }                                                           \
            }                                                               \
        }

    /* V 3' extension */
    int v_ambig = rec->v_trim_3 - pos->v_trim_3_adjusted;
    SCAN_EXTENSION(v_ambig, pos->v_sequence_end - v_ambig);

    /* J 5' extension */
    int j_ambig = rec->j_trim_5 - pos->j_trim_5_adjusted;
    SCAN_EXTENSION(j_ambig, pos->j_sequence_start);

    /* D 5' extension */
    int d5_ambig = rec->d_trim_5 - pos->d_trim_5_adjusted;
    SCAN_EXTENSION(d5_ambig, pos->d_sequence_start);

    /* D 3' extension */
    int d3_ambig = rec->d_trim_3 - pos->d_trim_3_adjusted;
    SCAN_EXTENSION(d3_ambig, pos->d_sequence_end - d3_ambig);

    #undef SCAN_EXTENSION
    return extra;
}

/* ═══════════════════════════════════════════════════════════════
 * airr_serialize — the main serialization function
 * ═══════════════════════════════════════════════════════════════ */

void airr_serialize(const ASeq *seq, const SimRecord *rec,
                     const SimConfig *cfg, const AlleleCorrectionSet *corr,
                     AirrRecord *out) {
    memset(out, 0, sizeof(*out));

    /* ── 1. Extract flat sequence ────────────────────────────── */
    out->sequence_length = aseq_to_string(seq, out->sequence,
                                           GENAIRR_MAX_SEQ_LEN);

    /* ── 2. Build germline alignment ─────────────────────────── */
    build_germline_alignment(seq, out->germline_alignment,
                              GENAIRR_MAX_SEQ_LEN);

    /* ── 3. Derive positions (FixV/D/JPosition) ─────────────── */
    AirrPositions pos;
    airr_derive_positions(seq, rec, &pos);

    /* ── 3b. Patch germline alignment for boundary extension ── */
    /* NP bases claimed by V/D/J via boundary ambiguity should show
     * the allele's germline bases, not 'N'. */
    patch_boundary_extensions(rec, &pos, out->germline_alignment,
                              out->sequence_length);

    out->v_sequence_start = pos.v_sequence_start;
    out->v_sequence_end   = pos.v_sequence_end;
    out->v_germline_start = pos.v_germline_start;
    out->v_germline_end   = pos.v_germline_end;

    out->d_sequence_start = pos.d_sequence_start;
    out->d_sequence_end   = pos.d_sequence_end;
    out->d_germline_start = pos.d_germline_start;
    out->d_germline_end   = pos.d_germline_end;

    out->j_sequence_start = pos.j_sequence_start;
    out->j_sequence_end   = pos.j_sequence_end;
    out->j_germline_start = pos.j_germline_start;
    out->j_germline_end   = pos.j_germline_end;

    out->junction_start  = pos.junction_start;
    out->junction_end    = pos.junction_end;
    out->junction_length = pos.junction_length;

    /* ── 4. Extract junction nucleotides and translate ────────── */
    if (pos.junction_length > 0 && pos.junction_length < AIRR_MAX_JUNCTION) {
        memcpy(out->junction_nt, out->sequence + pos.junction_start,
               pos.junction_length);
        out->junction_nt[pos.junction_length] = '\0';

        /* Translate junction to amino acids */
        int aa_len = 0;
        for (int i = 0; i + 2 < pos.junction_length; i += 3) {
            out->junction_aa[aa_len++] = translate_codon(
                out->junction_nt[i],
                out->junction_nt[i + 1],
                out->junction_nt[i + 2]);
        }
        out->junction_aa[aa_len] = '\0';
    }

    /* ── 5. Trim amounts (use adjusted values from position derivation) */
    out->v_trim_5 = rec->v_trim_5;
    out->v_trim_3 = pos.v_trim_3_adjusted;
    out->d_trim_5 = pos.d_trim_5_adjusted;
    out->d_trim_3 = pos.d_trim_3_adjusted;
    out->j_trim_5 = pos.j_trim_5_adjusted;
    out->j_trim_3 = rec->j_trim_3;

    /* ── 6. NP regions ───────────────────────────────────────── */
    out->np1_length = aseq_segment_to_string(seq, SEG_NP1,
                                              out->np1_region, AIRR_MAX_NP);
    out->np2_length = aseq_segment_to_string(seq, SEG_NP2,
                                              out->np2_region, AIRR_MAX_NP);

    /* ── 7. Collect mutations / errors ───────────────────────── */
    collect_annotations(seq, out);

    /* ── 7b. Collect mutations from boundary-extended positions ── */
    /* NP bases absorbed into V/D/J by boundary extension may differ
     * from the extended germline — a god-aligner would count those
     * as mutations. */
    {
        int mut_len = (int)strlen(out->mutations);
        int extra = collect_boundary_mutations(
            out->sequence, out->germline_alignment, out->sequence_length,
            &pos, rec, out->mutations, mut_len, AIRR_MAX_MUTATIONS);
        out->n_mutations += extra;
        /* Recompute mutation rate with extra mutations */
        int total_germ = 0;
        for (Nuc *n = seq->head; n; n = n->next)
            if (nuc_is_germline_segment(n)) total_germ++;
        out->mutation_rate = total_germ > 0
            ? (double)out->n_mutations / total_germ : 0.0;
    }

    /* ── 8. Productivity flags ────────────────────────────── */
    /* When the codon rail has been built (functionality_assess ran),
     * derive productivity reactively from the ASeq.  Otherwise fall
     * back to whatever the SimRecord already recorded.              */
    if (seq->codon_rail_valid) {
        out->productive  = aseq_is_productive(seq);
        out->stop_codon  = seq->n_stop_codons > 0;
        out->vj_in_frame = seq->junction_in_frame;
    } else {
        out->productive  = rec->productive;
        out->stop_codon  = rec->stop_codon;
        out->vj_in_frame = rec->vj_in_frame;
    }
    strncpy(out->note, rec->note, sizeof(out->note) - 1);

    /* ── 9. Sequence-level flags ─────────────────────────────── */
    out->is_reverse_complement = rec->is_reverse_complement;
    out->is_contaminant        = rec->is_contaminant;

    /* ── 10. Allele calls (bitmap-based reactive derivation) ──── */
    bool has_d = rec->d_allele != NULL && aseq_has_segment(seq, SEG_D);

    /* C allele: no correction needed (constant region, single call) */
    if (rec->c_allele)
        strncpy(out->c_call, rec->c_allele->name, AIRR_MAX_CALLS - 1);

    if (corr) {
        /* ── V call: derive from bitmap ──────────────────────── */
        if (rec->v_allele && corr->v_bitmap.n_alleles > 0) {
            int true_idx = allele_pool_find_index(&cfg->v_alleles,
                                                    rec->v_allele->name);
            AlleleCallResult vr;
            allele_call_derive(&corr->v_bitmap, seq, SEG_V, &vr, rec);
            allele_call_format(&vr, &cfg->v_alleles, true_idx,
                                out->v_call, AIRR_MAX_CALLS);
        } else if (rec->v_allele) {
            strncpy(out->v_call, rec->v_allele->name, AIRR_MAX_CALLS - 1);
        }

        /* ── J call: derive from bitmap ──────────────────────── */
        if (rec->j_allele && corr->j_bitmap.n_alleles > 0) {
            int true_idx = allele_pool_find_index(&cfg->j_alleles,
                                                    rec->j_allele->name);
            AlleleCallResult jr;
            allele_call_derive(&corr->j_bitmap, seq, SEG_J, &jr, rec);
            allele_call_format(&jr, &cfg->j_alleles, true_idx,
                                out->j_call, AIRR_MAX_CALLS);
        } else if (rec->j_allele) {
            strncpy(out->j_call, rec->j_allele->name, AIRR_MAX_CALLS - 1);
        }

        /* ── D call: derive from bitmap + short-D check ──────── */
        if (has_d && rec->d_allele && corr->d_bitmap.n_alleles > 0) {
            AlleleCallResult dr;
            allele_call_derive(&corr->d_bitmap, seq, SEG_D, &dr, rec);

            if (allele_call_is_short_d(&dr, &corr->d_bitmap)) {
                strncpy(out->d_call, "Short-D", AIRR_MAX_CALLS - 1);
            } else {
                int true_idx = allele_pool_find_index(&cfg->d_alleles,
                                                        rec->d_allele->name);
                allele_call_format(&dr, &cfg->d_alleles, true_idx,
                                    out->d_call, AIRR_MAX_CALLS);
            }
        } else if (has_d && rec->d_allele) {
            strncpy(out->d_call, rec->d_allele->name, AIRR_MAX_CALLS - 1);
        }
    } else {
        /* No corrections: use true allele names */
        if (rec->v_allele)
            strncpy(out->v_call, rec->v_allele->name, AIRR_MAX_CALLS - 1);
        if (rec->j_allele)
            strncpy(out->j_call, rec->j_allele->name, AIRR_MAX_CALLS - 1);
        if (has_d && rec->d_allele)
            strncpy(out->d_call, rec->d_allele->name, AIRR_MAX_CALLS - 1);
    }
}

/* ═══════════════════════════════════════════════════════════════
 * TSV output
 * ═══════════════════════════════════════════════════════════════ */

static const char *AIRR_COLUMNS[] = {
    "sequence", "germline_alignment",
    "v_call", "d_call", "j_call", "c_call",
    "v_sequence_start", "v_sequence_end",
    "v_germline_start", "v_germline_end",
    "d_sequence_start", "d_sequence_end",
    "d_germline_start", "d_germline_end",
    "j_sequence_start", "j_sequence_end",
    "j_germline_start", "j_germline_end",
    "junction_start", "junction_end", "junction_nt", "junction_aa", "junction_length",
    "v_trim_5", "v_trim_3", "d_trim_5", "d_trim_3", "j_trim_5", "j_trim_3",
    "np1_region", "np1_length", "np2_region", "np2_length",
    "mutation_rate", "n_mutations", "mutations",
    "n_sequencing_errors", "sequencing_errors",
    "n_pcr_errors", "pcr_errors",
    "productive", "stop_codon", "vj_in_frame", "note",
    "is_reverse_complement", "is_contaminant",
    "sequence_length",
};
#define AIRR_N_COLUMNS (sizeof(AIRR_COLUMNS) / sizeof(AIRR_COLUMNS[0]))

int airr_write_tsv_header(FILE *fp) {
    for (int i = 0; i < (int)AIRR_N_COLUMNS; i++) {
        if (i > 0) fputc('\t', fp);
        fputs(AIRR_COLUMNS[i], fp);
    }
    fputc('\n', fp);
    return (int)AIRR_N_COLUMNS;
}

void airr_write_tsv_row(FILE *fp, const AirrRecord *r) {
    fprintf(fp,
        "%s\t%s\t"
        "%s\t%s\t%s\t%s\t"
        "%d\t%d\t%d\t%d\t"
        "%d\t%d\t%d\t%d\t"
        "%d\t%d\t%d\t%d\t"
        "%d\t%d\t%s\t%s\t%d\t"
        "%d\t%d\t%d\t%d\t%d\t%d\t"
        "%s\t%d\t%s\t%d\t"
        "%.6f\t%d\t%s\t"
        "%d\t%s\t"
        "%d\t%s\t"
        "%s\t%s\t%s\t%s\t"
        "%s\t%s\t"
        "%d\n",
        r->sequence, r->germline_alignment,
        r->v_call, r->d_call, r->j_call, r->c_call,
        r->v_sequence_start, r->v_sequence_end,
        r->v_germline_start, r->v_germline_end,
        r->d_sequence_start, r->d_sequence_end,
        r->d_germline_start, r->d_germline_end,
        r->j_sequence_start, r->j_sequence_end,
        r->j_germline_start, r->j_germline_end,
        r->junction_start, r->junction_end,
        r->junction_nt, r->junction_aa, r->junction_length,
        r->v_trim_5, r->v_trim_3, r->d_trim_5, r->d_trim_3,
        r->j_trim_5, r->j_trim_3,
        r->np1_region, r->np1_length,
        r->np2_region, r->np2_length,
        r->mutation_rate, r->n_mutations, r->mutations,
        r->n_sequencing_errors, r->sequencing_errors,
        r->n_pcr_errors, r->pcr_errors,
        r->productive ? "T" : "F",
        r->stop_codon ? "T" : "F",
        r->vj_in_frame ? "T" : "F",
        r->note,
        r->is_reverse_complement ? "T" : "F",
        r->is_contaminant ? "T" : "F",
        r->sequence_length);
}

int airr_snprintf_tsv_row(char *buf, int size, const AirrRecord *r) {
    if (!buf || size <= 0 || !r) return -1;

    int n = snprintf(buf, size,
        "%s\t%s\t"
        "%s\t%s\t%s\t%s\t"
        "%d\t%d\t%d\t%d\t"
        "%d\t%d\t%d\t%d\t"
        "%d\t%d\t%d\t%d\t"
        "%d\t%d\t%s\t%s\t%d\t"
        "%d\t%d\t%d\t%d\t%d\t%d\t"
        "%s\t%d\t%s\t%d\t"
        "%.6f\t%d\t%s\t"
        "%d\t%s\t"
        "%d\t%s\t"
        "%s\t%s\t%s\t%s\t"
        "%s\t%s\t"
        "%d",
        r->sequence, r->germline_alignment,
        r->v_call, r->d_call, r->j_call, r->c_call,
        r->v_sequence_start, r->v_sequence_end,
        r->v_germline_start, r->v_germline_end,
        r->d_sequence_start, r->d_sequence_end,
        r->d_germline_start, r->d_germline_end,
        r->j_sequence_start, r->j_sequence_end,
        r->j_germline_start, r->j_germline_end,
        r->junction_start, r->junction_end,
        r->junction_nt, r->junction_aa, r->junction_length,
        r->v_trim_5, r->v_trim_3, r->d_trim_5, r->d_trim_3,
        r->j_trim_5, r->j_trim_3,
        r->np1_region, r->np1_length,
        r->np2_region, r->np2_length,
        r->mutation_rate, r->n_mutations, r->mutations,
        r->n_sequencing_errors, r->sequencing_errors,
        r->n_pcr_errors, r->pcr_errors,
        r->productive ? "T" : "F",
        r->stop_codon ? "T" : "F",
        r->vj_in_frame ? "T" : "F",
        r->note,
        r->is_reverse_complement ? "T" : "F",
        r->is_contaminant ? "T" : "F",
        r->sequence_length);

    if (n < 0 || n >= size) return -1;  /* truncated or error */
    return n;
}

int airr_snprintf_tsv_header(char *buf, int size) {
    if (!buf || size <= 0) return -1;

    int pos = 0;
    for (int i = 0; i < (int)AIRR_N_COLUMNS; i++) {
        if (i > 0) {
            if (pos >= size - 1) return -1;
            buf[pos++] = '\t';
        }
        int len = (int)strlen(AIRR_COLUMNS[i]);
        if (pos + len >= size) return -1;
        memcpy(buf + pos, AIRR_COLUMNS[i], len);
        pos += len;
    }
    buf[pos] = '\0';
    return (int)AIRR_N_COLUMNS;
}
