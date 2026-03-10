/**
 * codon.c — Codon translation implementation.
 *
 * Contains the standard genetic code table and the base→index
 * encoding. Used by the codon rail, functionality validator,
 * and selection pressure.
 */

#include "genairr/codon.h"

/* Encoding: T=0, C=1, A=2, G=3.  Index = b1*16 + b2*4 + b3. */
static const char CODON_TABLE[64] = {
    /* TTT */ 'F', /* TTC */ 'F', /* TTA */ 'L', /* TTG */ 'L',
    /* TCT */ 'S', /* TCC */ 'S', /* TCA */ 'S', /* TCG */ 'S',
    /* TAT */ 'Y', /* TAC */ 'Y', /* TAA */ '*', /* TAG */ '*',
    /* TGT */ 'C', /* TGC */ 'C', /* TGA */ '*', /* TGG */ 'W',
    /* CTT */ 'L', /* CTC */ 'L', /* CTA */ 'L', /* CTG */ 'L',
    /* CCT */ 'P', /* CCC */ 'P', /* CCA */ 'P', /* CCG */ 'P',
    /* CAT */ 'H', /* CAC */ 'H', /* CAA */ 'Q', /* CAG */ 'Q',
    /* CGT */ 'R', /* CGC */ 'R', /* CGA */ 'R', /* CGG */ 'R',
    /* ATT */ 'I', /* ATC */ 'I', /* ATA */ 'I', /* ATG */ 'M',
    /* ACT */ 'T', /* ACC */ 'T', /* ACA */ 'T', /* ACG */ 'T',
    /* AAT */ 'N', /* AAC */ 'N', /* AAA */ 'K', /* AAG */ 'K',
    /* AGT */ 'S', /* AGC */ 'S', /* AGA */ 'R', /* AGG */ 'R',
    /* GTT */ 'V', /* GTC */ 'V', /* GTA */ 'V', /* GTG */ 'V',
    /* GCT */ 'A', /* GCC */ 'A', /* GCA */ 'A', /* GCG */ 'A',
    /* GAT */ 'D', /* GAC */ 'D', /* GAA */ 'E', /* GAG */ 'E',
    /* GGT */ 'G', /* GGC */ 'G', /* GGA */ 'G', /* GGG */ 'G',
};

int codon_base_to_idx(char c) {
    switch (c) {
        case 'T': case 't': return 0;
        case 'C': case 'c': return 1;
        case 'A': case 'a': return 2;
        case 'G': case 'g': return 3;
        default: return -1;
    }
}

char translate_codon(char b1, char b2, char b3) {
    int i1 = codon_base_to_idx(b1);
    int i2 = codon_base_to_idx(b2);
    int i3 = codon_base_to_idx(b3);
    if (i1 < 0 || i2 < 0 || i3 < 0) return '?';
    return CODON_TABLE[i1 * 16 + i2 * 4 + i3];
}
