/**
 * codon.h — Codon translation utility.
 *
 * Provides a fast base→index encoding and a 64-entry codon table
 * for translating nucleotide triplets to amino acids. Shared by
 * the codon rail (aseq.c), functionality validator, and any future
 * amino acid-aware operations (e.g. selection pressure).
 */

#ifndef GENAIRR_CODON_H
#define GENAIRR_CODON_H

/**
 * Translate a single codon (three nucleotide characters) to an amino acid.
 *
 * Encoding: T=0, C=1, A=2, G=3. Index = b1*16 + b2*4 + b3.
 *
 * Returns:
 *   - Standard single-letter amino acid code (e.g. 'M', 'F', 'L', ...)
 *   - '*' for stop codons (TAA, TAG, TGA)
 *   - '?' if any base is not in {A, C, G, T} (handles N, ambiguous bases)
 */
char translate_codon(char b1, char b2, char b3);

/**
 * Encode a single nucleotide base to its codon table index.
 * Returns 0-3 for T/C/A/G (case-insensitive), -1 for invalid bases.
 */
int codon_base_to_idx(char c);

#endif /* GENAIRR_CODON_H */
