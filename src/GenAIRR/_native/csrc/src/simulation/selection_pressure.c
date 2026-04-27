/**
 * selection_pressure.c — Antigen-driven selection after SHM.
 *
 * For every mutated nucleotide in the *coding* segments (V/NP1/D/NP2/J),
 * classify the mutation as:
 *
 *   - Silent (S): always kept.
 *   - Replacement (R) in the V/J anchor codon: kept with
 *     `anchor_r_acceptance` (default 0.0). The conserved Cys (V) and
 *     W/F (J) anchors are structural: disrupting them breaks the V/J
 *     domain fold, so anchor R-mutations are nearly always purged.
 *   - Replacement (R) in a CDR position: kept with `cdr_r_acceptance`.
 *   - Replacement (R) in a FWR position: kept with `fwr_r_acceptance`.
 *
 * Effective acceptance = 1.0 − strength × (1.0 − base_acceptance), so
 * `selection_strength == 0` is a no-op and `1.0` is the configured
 * floor.
 *
 * Region classification (anchor-aware):
 *   - V segment: positions before the V anchor (Cys, IMGT pos 104 →
 *     ungapped ≈ 285) use the IMGT-gapped boundaries 78/114/165/195
 *     (FR1/CDR1/FR2/CDR2/FR3). Positions ≥ V anchor are CDR3.
 *   - NP1 / D / NP2: always CDR3.
 *   - J segment: positions strictly before (J-anchor + 3) are CDR3
 *     (the W/F codon is the last CDR3 codon); positions from
 *     (J-anchor + 3) onward are FR4.
 *
 * Anchorless V or J alleles (anchor ≤ 0): selection on that segment is
 * skipped to avoid mis-classifying CDR3 vs FR. NP1/D/NP2 continue to
 * be selected as CDR3. This matches the 34/106 builtin configs that
 * have anchorless V alleles and TCRB's anchorless J alleles.
 *
 * R/S classification uses the codon rail when available (O(1) codon
 * head lookup, codon spans segment boundaries cleanly). Falls back to
 * `germline_pos % 3` only when rail is unavailable (manual SimConfigs
 * in C tests). NP-inserted nodes (germline == '\0') count as
 * replacement by definition: they have no germline amino-acid to
 * compare against, and reverting them would require deletion.
 */

#include "genairr/pipeline.h"
#include "genairr/rand_util.h"
#include "genairr/trace.h"
#include <stdlib.h>

/* ── IMGT-gapped boundaries for V (ungapped germline_pos) ─────────
 * These are the canonical V-segment IMGT region boundaries used
 * throughout GenAIRR. They are gapped boundaries applied to ungapped
 * positions, which works for the canonical alleles in IMGT but can
 * drift slightly for alleles with insertions/deletions vs IMGT
 * reference. Per-allele boundary tables are tracked as a follow-up
 * ticket. */
#define V_FR1_END   78
#define V_CDR1_END  114
#define V_FR2_END   165
#define V_CDR2_END  195

typedef enum {
    REGION_FWR     = 0,
    REGION_CDR     = 1,
    REGION_ANCHOR  = 2,   /* V Cys / J W or F anchor codon */
    REGION_SKIP    = 3,   /* anchorless allele — segment is unselected */
} RegionType;

/* ── Codon translation (compact) ────────────────────────────────── */

static int base_idx(char c) {
    switch (c) {
        case 'T': case 't': return 0;
        case 'C': case 'c': return 1;
        case 'A': case 'a': return 2;
        case 'G': case 'g': return 3;
        default: return -1;
    }
}

static const char CODON_TABLE[64] = {
    'F','F','L','L', 'S','S','S','S', 'Y','Y','*','*', 'C','C','*','W',
    'L','L','L','L', 'P','P','P','P', 'H','H','Q','Q', 'R','R','R','R',
    'I','I','I','M', 'T','T','T','T', 'N','N','K','K', 'S','S','R','R',
    'V','V','V','V', 'A','A','A','A', 'D','D','E','E', 'G','G','G','G',
};

static char translate3(char b1, char b2, char b3) {
    int i1 = base_idx(b1), i2 = base_idx(b2), i3 = base_idx(b3);
    if (i1 < 0 || i2 < 0 || i3 < 0) return '?';
    return CODON_TABLE[i1 * 16 + i2 * 4 + i3];
}

/* ── Region classification ───────────────────────────────────── */

static RegionType classify_node_region(const Nuc *n,
                                       const SimRecord *rec,
                                       Nuc *codon_head) {
    /* Anchor-codon protection: any node whose codon contains the V or
     * J anchor flag is treated as REGION_ANCHOR. The anchor flag is
     * set on the head of the V Cys / J W/F codon by the assemble step. */
    if (codon_head && (codon_head->flags & NUC_FLAG_ANCHOR)) {
        return REGION_ANCHOR;
    }
    /* Some legacy paths set FLAG_ANCHOR on the node itself rather than
     * the codon head; cover that case too. */
    if (n->flags & NUC_FLAG_ANCHOR) {
        return REGION_ANCHOR;
    }

    switch (n->segment) {
        case SEG_V: {
            int va = (rec && rec->v_allele) ? rec->v_allele->anchor : 0;
            if (va <= 0) return REGION_SKIP;   /* anchorless V */
            if (n->germline_pos >= (uint16_t)va) return REGION_CDR;  /* CDR3 */
            if (n->germline_pos < V_FR1_END)  return REGION_FWR;
            if (n->germline_pos < V_CDR1_END) return REGION_CDR;
            if (n->germline_pos < V_FR2_END)  return REGION_FWR;
            if (n->germline_pos < V_CDR2_END) return REGION_CDR;
            return REGION_FWR;                 /* FR3 (up to V anchor) */
        }
        case SEG_NP1:
        case SEG_D:
        case SEG_NP2:
            return REGION_CDR;                 /* always CDR3 */
        case SEG_J: {
            int ja = (rec && rec->j_allele) ? rec->j_allele->anchor : 0;
            if (ja <= 0) return REGION_SKIP;   /* anchorless J */
            /* CDR3 ends at the W/F codon end (anchor + 3, exclusive). */
            if (n->germline_pos < (uint16_t)(ja + 3)) return REGION_CDR;
            return REGION_FWR;                 /* FR4 */
        }
        default:
            return REGION_SKIP;                /* C / UMI / adapter */
    }
}

/* ── R/S classification via codon rail ──────────────────────────── */

/**
 * True if reverting `node->current` to `node->germline` changes the
 * codon's translated amino-acid (replacement). NP-inserted nodes have
 * germline == '\0' and are always counted as replacement.
 *
 * Uses the codon rail when valid: walks 3 nodes from the codon head.
 * The codon may span a segment boundary (e.g. last 1-2 bases of V
 * combine with first base of NP1 to form the V/CDR3 transition codon)
 * — the rail handles this transparently because frame_phase is global
 * across the linked list.
 *
 * Falls back to `germline_pos % 3` when the rail is invalid (manual
 * SimConfig path in C tests). The fallback assumes the node sits in
 * a same-segment same-allele neighborhood and is not codon-spanning
 * across a segment boundary; this is acceptable for the test path.
 */
static bool is_replacement(const ASeq *seq, Nuc *node) {
    if (node->germline == '\0') return true;   /* NP / inserted */

    Nuc *n0 = NULL, *n1 = NULL, *n2 = NULL;

    if (seq->codon_rail_valid) {
        Nuc *head = nuc_codon_head_of(node);
        if (!head) return true;                 /* fragmentary codon at edge */
        n0 = head;
        n1 = head->next;
        n2 = n1 ? n1->next : NULL;
        if (!n1 || !n2) return true;            /* incomplete tail codon */
    } else {
        /* Fallback: align via germline_pos %% 3 (single-segment heuristic). */
        int offset = node->germline_pos % 3;
        if (offset == 0) {
            n0 = node;
            n1 = node->next;
            n2 = n1 ? n1->next : NULL;
        } else if (offset == 1 && node->prev) {
            n0 = node->prev;
            n1 = node;
            n2 = node->next;
        } else if (offset == 2 && node->prev && node->prev->prev) {
            n0 = node->prev->prev;
            n1 = node->prev;
            n2 = node;
        }
        if (!n0 || !n1 || !n2) return true;
    }

    /* Translate with current bases */
    char aa_current = translate3(n0->current, n1->current, n2->current);

    /* Translate with this node's germline base substituted in. Each
     * mutation is evaluated independently against the current codon
     * state; if two mutations land in the same codon and both arrive
     * here, each is judged "would reverting *me* alone change the
     * amino acid?". Order-stability is verified by the MC test. */
    char c0 = n0->current, c1 = n1->current, c2 = n2->current;
    if (n0 == node) c0 = node->germline;
    if (n1 == node) c1 = node->germline;
    if (n2 == node) c2 = node->germline;
    char aa_germline = translate3(c0, c1, c2);

    return aa_current != aa_germline;
}

/* ── Main step ──────────────────────────────────────────────────── */

void step_selection_pressure(const SimConfig *cfg, ASeq *seq, SimRecord *rec) {
    double strength = cfg->selection_strength;
    if (strength <= 0.0) return;

    double cdr_accept    = 1.0 - strength * (1.0 - cfg->cdr_r_acceptance);
    double fwr_accept    = 1.0 - strength * (1.0 - cfg->fwr_r_acceptance);
    double anchor_accept = 1.0 - strength * (1.0 - cfg->anchor_r_acceptance);

    TRACE("[selection] strength=%.2f, CDR_R=%.2f, FWR_R=%.2f, ANCHOR_R=%.2f",
          strength, cdr_accept, fwr_accept, anchor_accept);

    int reverted_cdr = 0, reverted_fwr = 0, reverted_anchor = 0;
    int kept_r = 0, kept_s = 0, skipped_seg = 0;

    /* Walk the entire sequence; only act on coding-segment mutations. */
    for (Nuc *n = seq->head; n; n = n->next) {
        if (!(n->flags & NUC_FLAG_MUTATED)) continue;
        if (!nuc_is_coding_segment(n))      continue;

        /* Silent mutations always kept. */
        if (!is_replacement(seq, n)) { kept_s++; continue; }

        Nuc *head = seq->codon_rail_valid ? nuc_codon_head_of(n) : NULL;
        RegionType region = classify_node_region(n, rec, head);

        if (region == REGION_SKIP) { skipped_seg++; continue; }

        double accept;
        switch (region) {
            case REGION_ANCHOR: accept = anchor_accept; break;
            case REGION_CDR:    accept = cdr_accept;    break;
            case REGION_FWR:
            default:            accept = fwr_accept;    break;
        }

        if (rng_uniform(cfg->rng) < accept) { kept_r++; continue; }

        aseq_revert(seq, n);
        if      (region == REGION_ANCHOR) reverted_anchor++;
        else if (region == REGION_CDR)    reverted_cdr++;
        else                              reverted_fwr++;
    }

    TRACE("[selection] kept %d S + %d R, reverted %d FWR + %d CDR + %d ANCHOR, skipped %d unselectable",
          kept_s, kept_r, reverted_fwr, reverted_cdr, reverted_anchor, skipped_seg);
}
