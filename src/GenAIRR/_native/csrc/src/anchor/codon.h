/**
 * codon.h — internal codon-level helpers shared across anchor strategies.
 *
 * NOT a public API. Lives in src/anchor/ and is consumed only by the
 * strategy_*.c files and result.c.
 */

#ifndef GENAIRR_ANCHOR_CODON_H
#define GENAIRR_ANCHOR_CODON_H

#include <stdbool.h>
#include "genairr/anchor.h"

/* Lowercase the first 3 characters of `src` into a 4-byte NUL-terminated
 * buffer. `src` must point at ≥ 3 readable bytes. */
void codon_normalize(const char *src, char out[4]);

/* True if all 3 chars are unambiguous DNA (A/C/G/T, case-insensitive). */
bool codon_is_canonical(const char *codon);

/* True if any of the 3 chars is an IUPAC ambiguity code (Y/R/S/W/K/M/B/D/H/V/N).
 * Counter-example: 'Z' is invalid and returns false (caller treats as garbage). */
bool codon_has_iupac(const char *codon);

/* Classify a codon at the V/J anchor position into one of the
 * conserved residues used by the simulator. Returns:
 *
 *   ANCHOR_CYS    for tgt, tgc, and tg[y] (Y = C/T ambiguity)
 *   ANCHOR_TRP    for tgg
 *   ANCHOR_PHE    for ttt, ttc, and tt[y]
 *   ANCHOR_RES_UNKNOWN otherwise (includes stop codons and unrelated codons)
 *
 * Case-insensitive. Caller is responsible for first checking
 * `codon_is_stop()` if pseudogene rejection is desired — this
 * classifier alone does not distinguish "stop codon" from "non-anchor
 * codon".
 */
AnchorResidue codon_classify(const char *codon);

/* True if codon is a stop codon. Recognizes TAA, TAG, TGA, and the
 * IUPAC-aware TAR (R = A/G, covers TAA+TAG). Case-insensitive. */
bool codon_is_stop(const char *codon);

#endif /* GENAIRR_ANCHOR_CODON_H */
