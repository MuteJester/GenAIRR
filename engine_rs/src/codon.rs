//! Genetic-code arithmetic: standard codon-to-amino-acid translation
//! plus the canonical "stop" and "ambiguous" sentinel bytes.
//!
//! Lifted out of `ir` so that consumers that only need translation
//! (the contract bundle, the feasibility builder, the AIRR-record
//! column walker) don't have to depend on the persistent-IR types.

use crate::ir::encode_base;

/// Standard genetic code (NCBI table 1). Indexed by
/// `(b1<<4) | (b2<<2) | b3` with the 0..=3 base encoding from
/// `encode_base`. Stop codons (TAA, TAG, TGA) emit `b'*'`.
#[rustfmt::skip]
const GENETIC_CODE: [u8; 64] = [
    // AAA   AAC   AAG   AAT
    b'K', b'N', b'K', b'N',
    // ACA   ACC   ACG   ACT
    b'T', b'T', b'T', b'T',
    // AGA   AGC   AGG   AGT
    b'R', b'S', b'R', b'S',
    // ATA   ATC   ATG   ATT
    b'I', b'I', b'M', b'I',
    // CAA   CAC   CAG   CAT
    b'Q', b'H', b'Q', b'H',
    // CCA   CCC   CCG   CCT
    b'P', b'P', b'P', b'P',
    // CGA   CGC   CGG   CGT
    b'R', b'R', b'R', b'R',
    // CTA   CTC   CTG   CTT
    b'L', b'L', b'L', b'L',
    // GAA   GAC   GAG   GAT
    b'E', b'D', b'E', b'D',
    // GCA   GCC   GCG   GCT
    b'A', b'A', b'A', b'A',
    // GGA   GGC   GGG   GGT
    b'G', b'G', b'G', b'G',
    // GTA   GTC   GTG   GTT
    b'V', b'V', b'V', b'V',
    // TAA   TAC   TAG   TAT
    b'*', b'Y', b'*', b'Y',
    // TCA   TCC   TCG   TCT
    b'S', b'S', b'S', b'S',
    // TGA   TGC   TGG   TGT
    b'*', b'C', b'W', b'C',
    // TTA   TTC   TTG   TTT
    b'L', b'F', b'L', b'F',
];

/// Special amino-acid byte for a codon containing any non-{A,C,G,T,U}
/// base — e.g. an N from quality masking, or an indel-inserted byte
/// that hasn't been resolved.
pub const AMINO_AMBIGUOUS: u8 = b'X';

/// Special amino-acid byte for a stop codon (TAA, TAG, TGA).
pub const AMINO_STOP: u8 = b'*';

/// Translate one codon to an amino-acid byte. Returns `AMINO_AMBIGUOUS`
/// (`b'X'`) if any of the three input bases is not one of A/C/G/T/U
/// (case-insensitive).
pub fn translate_codon(b1: u8, b2: u8, b3: u8) -> u8 {
    match (encode_base(b1), encode_base(b2), encode_base(b3)) {
        (Some(i1), Some(i2), Some(i3)) => {
            let idx = ((i1 as usize) << 4) | ((i2 as usize) << 2) | (i3 as usize);
            GENETIC_CODE[idx]
        }
        _ => AMINO_AMBIGUOUS,
    }
}

/// Slice form of `translate_codon`. Returns `'X'` for any input that
/// isn't exactly three bytes, otherwise returns the amino-acid char
/// (one of the standard 20, `'*'` for stops, or `'X'` for ambiguous).
pub fn translate_codon_slice(codon: &[u8]) -> char {
    if codon.len() != 3 {
        return 'X';
    }
    translate_codon(codon[0], codon[1], codon[2]) as char
}

/// Translate an in-frame nucleotide string into amino-acid chars.
/// Reads codons from the start and stops at the last full triplet;
/// trailing partial codons are ignored.
pub fn translate_seq(seq: &str) -> String {
    let bytes = seq.as_bytes();
    let n_codons = bytes.len() / 3;
    let mut out = String::with_capacity(n_codons);
    for i in 0..n_codons {
        let start = i * 3;
        out.push(translate_codon_slice(&bytes[start..start + 3]));
    }
    out
}
