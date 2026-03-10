/**
 * gdc_format.h — GenAIRR DataConfig binary format (.gdc).
 *
 * Self-contained binary format encoding all reference data needed for
 * BCR/TCR V(D)J simulation:
 *
 *   METADATA      Species, chain type, config name
 *   ALLELES       V/D/J/C gene pools with sequences and anchors
 *   GENE_USE      Per-segment gene usage probabilities
 *   TRIM_DISTS    Per-family/gene trim probability arrays
 *   NP_PARAMS     TdT model: length dists, first-base probs, Markov chains
 *   P_NUC_PROBS   P-nucleotide length distribution
 *   MUTATION_MODEL S5F 5-mer mutability + substitution tables
 *
 * Design goals:
 *   - Cross-language: writable by Python, readable by C (and vice versa)
 *   - Section-based with TOC: readers skip unknown sections (forward compat)
 *   - Compact: dense arrays for probability data, no JSON/text overhead
 *   - Little-endian throughout (x86/ARM native)
 *
 * String encoding: uint16_t length prefix + raw bytes (no null terminator).
 * Floats: IEEE 754 double (8 bytes).
 */

#ifndef GENAIRR_GDC_FORMAT_H
#define GENAIRR_GDC_FORMAT_H

#include <stdint.h>

/* ── Magic and version ────────────────────────────────────────── */

#define GDC_MAGIC_0  'G'
#define GDC_MAGIC_1  'D'
#define GDC_MAGIC_2  'C'
#define GDC_MAGIC_3  '\x01'

#define GDC_FORMAT_VERSION  1

/* ── Section IDs ──────────────────────────────────────────────── */

typedef enum {
    GDC_SECTION_METADATA       = 0,
    GDC_SECTION_ALLELES        = 1,
    GDC_SECTION_GENE_USE       = 2,
    GDC_SECTION_TRIM_DISTS     = 3,
    GDC_SECTION_NP_PARAMS      = 4,
    GDC_SECTION_P_NUC_PROBS    = 5,
    GDC_SECTION_MUTATION_MODEL = 6,
    GDC_SECTION_COUNT                 /* sentinel */
} GdcSectionId;

/* ── Chain type encoding (matches C ChainType enum) ───────────── */

#define GDC_CHAIN_IGH   0
#define GDC_CHAIN_IGK   1
#define GDC_CHAIN_IGL   2
#define GDC_CHAIN_TCRA  3
#define GDC_CHAIN_TCRB  4
#define GDC_CHAIN_TCRD  5
#define GDC_CHAIN_TCRG  6

/* ── Segment type encoding ────────────────────────────────────── */

#define GDC_SEG_V  0
#define GDC_SEG_D  1
#define GDC_SEG_J  2
#define GDC_SEG_C  3

/* ── Trim type encoding ───────────────────────────────────────── */

#define GDC_TRIM_V3  0
#define GDC_TRIM_D5  1
#define GDC_TRIM_D3  2
#define GDC_TRIM_J5  3

/* ── Mutation model type ──────────────────────────────────────── */

#define GDC_MUTMODEL_NONE     0
#define GDC_MUTMODEL_S5F      1
#define GDC_MUTMODEL_UNIFORM  2

/* ── S5F constants ────────────────────────────────────────────── */

#define GDC_S5F_KMER_SPACE  3125   /* 5^5 (A=0,C=1,G=2,T=3,N=4) */

/* ── On-disk structures (packed, little-endian) ───────────────── */

#pragma pack(push, 1)

typedef struct {
    uint8_t   magic[4];          /* "GDC\x01"                     */
    uint16_t  format_version;    /* GDC_FORMAT_VERSION             */
    uint16_t  n_sections;        /* number of section table entries */
    uint8_t   reserved[24];      /* zero-filled, future use        */
} GdcFileHeader;                 /* 32 bytes                       */

typedef struct {
    uint32_t  section_id;        /* GdcSectionId                   */
    uint32_t  size;              /* section payload size in bytes   */
    uint64_t  offset;            /* absolute file offset            */
} GdcSectionEntry;               /* 16 bytes                       */

#pragma pack(pop)

/* ── NP base index encoding ───────────────────────────────────── */

/* Row/col order for 4×4 Markov matrices: A=0, C=1, G=2, T=3 */
#define GDC_BASE_A  0
#define GDC_BASE_C  1
#define GDC_BASE_G  2
#define GDC_BASE_T  3

static inline int gdc_base_index(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return -1;
    }
}

#endif /* GENAIRR_GDC_FORMAT_H */
