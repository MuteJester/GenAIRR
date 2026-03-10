"""
gdc_io.py — Read and write GenAIRR DataConfig binary files (.gdc).

The .gdc format is a self-contained binary encoding of all reference data
needed for BCR/TCR V(D)J simulation: alleles, gene usage, trim distributions,
NP parameters, P-nucleotide probs, and S5F mutation model tables.

Cross-language: this module writes files readable by the C genairr library,
and reads files written by either Python or C.

All integers are little-endian. Strings are uint16 length-prefixed.
Floats are IEEE 754 double (8 bytes).
"""

import struct
import math
import pickle
from pathlib import Path
from typing import Dict, Optional, Any
from datetime import date

from .data_config import DataConfig
from .enums import ChainType, Species
from .config_info import ConfigInfo

# ═════════════════════════════════════════════════════════════════════
# Constants — must match gdc_format.h exactly
# ═════════════════════════════════════════════════════════════════════

GDC_MAGIC = b'GDC\x01'
GDC_FORMAT_VERSION = 1

# Section IDs
SEC_METADATA       = 0
SEC_ALLELES        = 1
SEC_GENE_USE       = 2
SEC_TRIM_DISTS     = 3
SEC_NP_PARAMS      = 4
SEC_P_NUC_PROBS    = 5
SEC_MUTATION_MODEL = 6
SEC_COUNT          = 7

# Chain type encoding
_CHAIN_TO_INT = {
    ChainType.BCR_HEAVY:        0,
    ChainType.BCR_LIGHT_KAPPA:  1,
    ChainType.BCR_LIGHT_LAMBDA: 2,
    ChainType.TCR_ALPHA:        3,
    ChainType.TCR_BETA:         4,
    ChainType.TCR_DELTA:        5,
    ChainType.TCR_GAMMA:        6,
}
_INT_TO_CHAIN = {v: k for k, v in _CHAIN_TO_INT.items()}

# Segment type encoding
_SEG_V, _SEG_D, _SEG_J, _SEG_C = 0, 1, 2, 3

# Trim type encoding
_TRIM_V3, _TRIM_D5, _TRIM_D3, _TRIM_J5 = 0, 1, 2, 3
_TRIM_KEY_TO_INT = {"V_3": 0, "D_5": 1, "D_3": 2, "J_5": 3}

# Mutation model types
_MUTMODEL_NONE    = 0
_MUTMODEL_S5F     = 1
_MUTMODEL_UNIFORM = 2

# S5F key space
_S5F_KMER_SPACE = 3125  # 5^5

# S5F base encoding: A=0, C=1, G=2, T=3, N=4
_S5F_BASE = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}


def _s5f_str_to_key(kmer: str) -> int:
    """Encode a 5-mer string to an integer key (base-5 positional)."""
    b = _S5F_BASE
    try:
        return (b[kmer[0]] * 625 + b[kmer[1]] * 125 +
                b[kmer[2]] * 25 + b[kmer[3]] * 5 + b[kmer[4]])
    except (KeyError, IndexError):
        return -1


def _s5f_key_to_str(key: int) -> str:
    """Decode an integer key back to a 5-mer string."""
    bases = "ACGTN"
    s = []
    for _ in range(5):
        s.append(bases[key % 5])
        key //= 5
    return ''.join(reversed(s))


# NP base index
_NP_BASE_IDX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
_NP_IDX_BASE = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}


# ═════════════════════════════════════════════════════════════════════
# Low-level binary I/O
# ═════════════════════════════════════════════════════════════════════

class _BinaryWriter:
    """Buffered little-endian binary writer."""
    __slots__ = ('_buf',)

    def __init__(self):
        self._buf = bytearray()

    def tell(self) -> int:
        return len(self._buf)

    def u8(self, v: int):
        self._buf.append(v & 0xFF)

    def u16(self, v: int):
        self._buf.extend(struct.pack('<H', v))

    def u32(self, v: int):
        self._buf.extend(struct.pack('<I', v))

    def u64(self, v: int):
        self._buf.extend(struct.pack('<Q', v))

    def i16(self, v: int):
        self._buf.extend(struct.pack('<h', v))

    def f64(self, v: float):
        self._buf.extend(struct.pack('<d', v))

    def string(self, s: str):
        encoded = s.encode('utf-8')
        self.u16(len(encoded))
        self._buf.extend(encoded)

    def f64_array(self, arr):
        self._buf.extend(struct.pack(f'<{len(arr)}d', *arr))

    def raw(self, data: bytes):
        self._buf.extend(data)

    def patch_at(self, offset: int, data: bytes):
        self._buf[offset:offset + len(data)] = data

    def getvalue(self) -> bytes:
        return bytes(self._buf)


class _BinaryReader:
    """Buffered little-endian binary reader."""
    __slots__ = ('_data', '_pos')

    def __init__(self, data: bytes):
        self._data = data
        self._pos = 0

    def seek(self, pos: int):
        self._pos = pos

    def tell(self) -> int:
        return self._pos

    def _read(self, n: int) -> bytes:
        end = self._pos + n
        if end > len(self._data):
            raise ValueError(f"Read past end of data at offset {self._pos}")
        chunk = self._data[self._pos:end]
        self._pos = end
        return chunk

    def u8(self) -> int:
        return self._read(1)[0]

    def u16(self) -> int:
        return struct.unpack_from('<H', self._read(2))[0]

    def u32(self) -> int:
        return struct.unpack_from('<I', self._read(4))[0]

    def u64(self) -> int:
        return struct.unpack_from('<Q', self._read(8))[0]

    def i16(self) -> int:
        return struct.unpack_from('<h', self._read(2))[0]

    def f64(self) -> float:
        return struct.unpack_from('<d', self._read(8))[0]

    def string(self) -> str:
        length = self.u16()
        return self._read(length).decode('utf-8')

    def f64_array(self, n: int) -> list:
        data = self._read(n * 8)
        return list(struct.unpack_from(f'<{n}d', data))

    def raw(self, n: int) -> bytes:
        return self._read(n)


# ═════════════════════════════════════════════════════════════════════
# Section writers
# ═════════════════════════════════════════════════════════════════════

def _write_metadata(w: _BinaryWriter, config: DataConfig):
    w.string(config.name or "")

    chain_int = 0
    has_d = False
    species_str = ""
    ref_str = ""
    build_days = 0

    if config.metadata:
        chain_int = _CHAIN_TO_INT.get(config.metadata.chain_type, 0)
        has_d = config.metadata.has_d
        species_str = config.metadata.species.value if config.metadata.species else ""
        ref_str = config.metadata.reference_set or ""
        if config.metadata.last_updated:
            epoch = date(2000, 1, 1)
            build_days = (config.metadata.last_updated - epoch).days

    w.u8(chain_int)
    w.u8(1 if has_d else 0)
    w.string(species_str)
    w.string(ref_str)
    w.u32(build_days)


def _write_alleles(w: _BinaryWriter, config: DataConfig):
    for attr in ('v_alleles', 'd_alleles', 'j_alleles', 'c_alleles'):
        allele_dict = getattr(config, attr)
        if not allele_dict:
            w.u16(0)
            continue

        # Flatten gene groups into a flat allele list
        all_alleles = [a for gene_list in allele_dict.values() for a in gene_list]
        w.u16(len(all_alleles))

        for allele in all_alleles:
            w.string(allele.name)
            w.string(allele.gene)
            w.string(allele.family)
            seq = allele.ungapped_seq
            w.u16(len(seq))
            w.raw(seq.encode('ascii'))
            # anchor=0 means "not found" (set by keep_anchorless=True);
            # a real V/J anchor is never at position 0.  Write -1 so the
            # C engine treats the allele as anchorless → non-productive.
            anchor = allele.anchor if allele.anchor else -1
            w.i16(anchor)


def _write_gene_use(w: _BinaryWriter, config: DataConfig):
    seg_map = {'V': _SEG_V, 'D': _SEG_D, 'J': _SEG_J, 'C': _SEG_C}
    active = {k: v for k, v in config.gene_use_dict.items() if v and k in seg_map}

    w.u8(len(active))
    for key, usage_dict in active.items():
        w.u8(seg_map[key])
        w.u32(len(usage_dict))
        for gene_name, prob in usage_dict.items():
            w.string(gene_name)
            w.f64(prob)


def _write_trim_dists(w: _BinaryWriter, config: DataConfig):
    active_types = []
    for trim_key, trim_int in _TRIM_KEY_TO_INT.items():
        if trim_key in config.trim_dicts and config.trim_dicts[trim_key]:
            active_types.append((trim_int, config.trim_dicts[trim_key]))

    w.u8(len(active_types))

    for trim_int, family_dict in active_types:
        w.u8(trim_int)

        # Flatten: family → gene → probs becomes list of (family, gene, probs)
        entries = []
        for family_name, gene_dict in family_dict.items():
            for gene_name, prob_dict in gene_dict.items():
                entries.append((family_name, gene_name, prob_dict))

        w.u16(len(entries))

        for family_name, gene_name, prob_dict in entries:
            w.string(family_name)
            w.string(gene_name)

            if not prob_dict:
                w.u16(0)
                continue

            max_trim = max(prob_dict.keys())
            w.u16(max_trim)

            # Dense array [0..max_trim]
            probs = [prob_dict.get(i, 0.0) for i in range(max_trim + 1)]
            w.f64_array(probs)


def _write_np_params(w: _BinaryWriter, config: DataConfig):
    regions = []
    for key in ('NP1', 'NP2'):
        if key in config.NP_lengths and config.NP_lengths[key]:
            regions.append(key)

    w.u8(len(regions))

    for key in regions:
        # Length distribution
        length_dict = config.NP_lengths.get(key, {})
        max_length = max(length_dict.keys()) if length_dict else 0
        w.u16(max_length)
        length_probs = [length_dict.get(i, 0.0) for i in range(max_length + 1)]
        w.f64_array(length_probs)

        # First base probabilities (A, C, G, T order)
        first_bases = config.NP_first_bases.get(key, {})
        for base in ('A', 'C', 'G', 'T'):
            w.f64(first_bases.get(base, 0.25))

        # Markov transitions: position → from_base → {to_base: prob}
        trans = config.NP_transitions.get(key, {})
        n_positions = max(trans.keys()) + 1 if trans else 0
        w.u16(n_positions)

        for pos in range(n_positions):
            pos_dict = trans.get(pos, {})
            # Write 4×4 matrix: from_base × to_base
            for from_base in ('A', 'C', 'G', 'T'):
                from_dict = pos_dict.get(from_base, {})
                for to_base in ('A', 'C', 'G', 'T'):
                    w.f64(from_dict.get(to_base, 0.25))


def _write_p_nuc_probs(w: _BinaryWriter, config: DataConfig):
    probs = config.p_nucleotide_length_probs
    max_len = max(probs.keys()) if probs else 0
    w.u8(max_len)
    arr = [probs.get(i, 0.0) for i in range(max_len + 1)]
    w.f64_array(arr)


def _write_mutation_model(w: _BinaryWriter, config: DataConfig,
                          mutability: Optional[Dict] = None,
                          substitution: Optional[Dict] = None):
    if mutability is None or substitution is None:
        w.u8(_MUTMODEL_NONE)
        return

    w.u8(_MUTMODEL_S5F)

    # Mutability: 3125 doubles indexed by integer key
    mut_array = [0.0] * _S5F_KMER_SPACE
    for kmer, val in mutability.items():
        k = _s5f_str_to_key(kmer)
        if 0 <= k < _S5F_KMER_SPACE:
            mut_array[k] = val if isinstance(val, (int, float)) and not math.isnan(val) else 0.0
    w.f64_array(mut_array)

    # Substitution: per 5-mer, count + bases + weights
    for k in range(_S5F_KMER_SPACE):
        kmer = _s5f_key_to_str(k)
        subs = substitution.get(kmer, {})
        if not subs:
            w.u8(0)
            continue

        # Filter valid entries
        entries = [(base, weight) for base, weight in subs.items()
                   if isinstance(weight, (int, float)) and not math.isnan(weight)]
        if not entries:
            w.u8(0)
            continue

        w.u8(len(entries))
        for base, _ in entries:
            w.raw(base.encode('ascii'))
        for _, weight in entries:
            w.f64(weight)


# ═════════════════════════════════════════════════════════════════════
# Section readers
# ═════════════════════════════════════════════════════════════════════

def _read_metadata(r: _BinaryReader, config: DataConfig):
    name = r.string()
    config.name = name if name else None
    chain_int = r.u8()
    has_d = r.u8() != 0
    species_str = r.string()
    ref_str = r.string()
    build_days = r.u32()

    chain_type = _INT_TO_CHAIN.get(chain_int, ChainType.BCR_HEAVY)

    species = None
    for sp in Species:
        if sp.value == species_str:
            species = sp
            break

    epoch = date(2000, 1, 1)
    from datetime import timedelta
    last_updated = epoch + timedelta(days=build_days) if build_days else None

    config.metadata = ConfigInfo(
        species=species,
        chain_type=chain_type,
        reference_set=ref_str,
        last_updated=last_updated,
        has_d=has_d,
    )


def _read_alleles(r: _BinaryReader, config: DataConfig):
    from ..alleles.allele import VAllele, DAllele, JAllele, CAllele

    allele_classes = [VAllele, DAllele, JAllele, CAllele]
    attr_names = ['v_alleles', 'd_alleles', 'j_alleles', 'c_alleles']

    for cls, attr in zip(allele_classes, attr_names):
        count = r.u16()
        if count == 0:
            setattr(config, attr, None)
            continue

        gene_groups = {}
        for _ in range(count):
            name = r.string()
            gene = r.string()
            family = r.string()
            seq_len = r.u16()
            seq = r.raw(seq_len).decode('ascii')
            anchor = r.i16()

            # Reconstruct allele with gapped_seq = ungapped_seq (no gaps in binary)
            # The anchor is stored directly, bypassing _find_anchor
            allele = cls(name=name, gapped_sequence=seq,
                         length=len(seq), anchor_override=anchor)

            if gene not in gene_groups:
                gene_groups[gene] = []
            gene_groups[gene].append(allele)

        setattr(config, attr, gene_groups)


def _read_gene_use(r: _BinaryReader, config: DataConfig):
    int_to_seg = {_SEG_V: 'V', _SEG_D: 'D', _SEG_J: 'J', _SEG_C: 'C'}
    n_types = r.u8()

    config.gene_use_dict = {}
    for _ in range(n_types):
        seg = int_to_seg[r.u8()]
        count = r.u32()
        usage = {}
        for _ in range(count):
            gene_name = r.string()
            prob = r.f64()
            usage[gene_name] = prob
        config.gene_use_dict[seg] = usage


def _read_trim_dists(r: _BinaryReader, config: DataConfig):
    int_to_trim = {0: "V_3", 1: "D_5", 2: "D_3", 3: "J_5"}
    n_types = r.u8()

    config.trim_dicts = {}
    for _ in range(n_types):
        trim_key = int_to_trim[r.u8()]
        count = r.u16()

        family_dict = {}
        for _ in range(count):
            family_name = r.string()
            gene_name = r.string()
            max_trim = r.u16()
            probs = r.f64_array(max_trim + 1)

            if family_name not in family_dict:
                family_dict[family_name] = {}
            family_dict[family_name][gene_name] = {
                i: p for i, p in enumerate(probs)
            }

        config.trim_dicts[trim_key] = family_dict


def _read_np_params(r: _BinaryReader, config: DataConfig):
    n_regions = r.u8()
    np_keys = ['NP1', 'NP2']

    config.NP_lengths = {}
    config.NP_first_bases = {}
    config.NP_transitions = {}

    for i in range(n_regions):
        key = np_keys[i]

        # Length distribution
        max_length = r.u16()
        length_probs = r.f64_array(max_length + 1)
        config.NP_lengths[key] = {i: p for i, p in enumerate(length_probs)}

        # First base probs
        bases_arr = r.f64_array(4)
        config.NP_first_bases[key] = {
            'A': bases_arr[0], 'C': bases_arr[1],
            'G': bases_arr[2], 'T': bases_arr[3],
        }

        # Markov transitions
        n_positions = r.u16()
        transitions = {}
        base_order = ('A', 'C', 'G', 'T')
        for pos in range(n_positions):
            pos_dict = {}
            for from_idx, from_base in enumerate(base_order):
                row = r.f64_array(4)
                pos_dict[from_base] = {
                    to_base: row[to_idx]
                    for to_idx, to_base in enumerate(base_order)
                }
            transitions[pos] = pos_dict
        config.NP_transitions[key] = transitions


def _read_p_nuc_probs(r: _BinaryReader, config: DataConfig):
    max_len = r.u8()
    probs = r.f64_array(max_len + 1)
    config.p_nucleotide_length_probs = {i: p for i, p in enumerate(probs)}


def _read_mutation_model(r: _BinaryReader) -> tuple:
    """Returns (mutability_dict, substitution_dict) or (None, None)."""
    model_type = r.u8()

    if model_type != _MUTMODEL_S5F:
        return None, None

    # Mutability
    mut_array = r.f64_array(_S5F_KMER_SPACE)
    mutability = {}
    for k in range(_S5F_KMER_SPACE):
        kmer = _s5f_key_to_str(k)
        mutability[kmer] = mut_array[k]

    # Substitution
    substitution = {}
    for k in range(_S5F_KMER_SPACE):
        kmer = _s5f_key_to_str(k)
        cnt = r.u8()
        if cnt == 0:
            substitution[kmer] = {}
            continue

        bases = [chr(b) for b in r.raw(cnt)]
        weights = r.f64_array(cnt)
        substitution[kmer] = {base: weight for base, weight in zip(bases, weights)}

    return mutability, substitution


# ═════════════════════════════════════════════════════════════════════
# Public API
# ═════════════════════════════════════════════════════════════════════

def save_gdc(path: str, config: DataConfig, *,
             mutability: Optional[Dict] = None,
             substitution: Optional[Dict] = None,
             s5f_chain_category: Optional[str] = None) -> None:
    """
    Write a DataConfig to a .gdc binary file.

    Args:
        path: Output file path.
        config: The DataConfig to serialize.
        mutability: S5F mutability dict (5-mer string → float).
                    If None and s5f_chain_category is set, auto-loads from package data.
        substitution: S5F substitution dict (5-mer → {base → float}).
                      If None and s5f_chain_category is set, auto-loads from package data.
        s5f_chain_category: "heavy" or "light" — auto-loads S5F tables from
                            built-in package data. Ignored if mutability/substitution
                            are provided directly.
    """
    # Auto-load S5F if requested
    if mutability is None and substitution is None and s5f_chain_category:
        mutability, substitution, _ = _load_s5f_tables(s5f_chain_category)

    w = _BinaryWriter()

    # File header (32 bytes)
    w.raw(GDC_MAGIC)
    w.u16(GDC_FORMAT_VERSION)
    w.u16(SEC_COUNT)
    w.raw(b'\x00' * 24)  # reserved

    # Section table placeholder (SEC_COUNT × 16 bytes)
    toc_offset = w.tell()
    for _ in range(SEC_COUNT):
        w.u32(0)  # section_id
        w.u32(0)  # size
        w.u64(0)  # offset

    # Write sections, recording offsets and sizes
    section_writers = [
        (SEC_METADATA,       lambda: _write_metadata(w, config)),
        (SEC_ALLELES,        lambda: _write_alleles(w, config)),
        (SEC_GENE_USE,       lambda: _write_gene_use(w, config)),
        (SEC_TRIM_DISTS,     lambda: _write_trim_dists(w, config)),
        (SEC_NP_PARAMS,      lambda: _write_np_params(w, config)),
        (SEC_P_NUC_PROBS,    lambda: _write_p_nuc_probs(w, config)),
        (SEC_MUTATION_MODEL, lambda: _write_mutation_model(w, config, mutability, substitution)),
    ]

    entries = []
    for sec_id, writer_fn in section_writers:
        start = w.tell()
        writer_fn()
        size = w.tell() - start
        entries.append((sec_id, size, start))

    # Patch section table
    for i, (sec_id, size, offset) in enumerate(entries):
        entry_offset = toc_offset + i * 16
        w.patch_at(entry_offset, struct.pack('<IIQ', sec_id, size, offset))

    Path(path).write_bytes(w.getvalue())


def to_gdc_bytes(config: DataConfig, *,
                 mutability: Optional[Dict] = None,
                 substitution: Optional[Dict] = None,
                 s5f_chain_category: Optional[str] = None) -> bytes:
    """
    Serialize a DataConfig to GDC binary format in memory.

    Same as save_gdc() but returns bytes instead of writing to a file.
    Used by the C backend bridge to pass configs without temp files.
    """
    if mutability is None and substitution is None and s5f_chain_category:
        mutability, substitution, _ = _load_s5f_tables(s5f_chain_category)

    w = _BinaryWriter()

    w.raw(GDC_MAGIC)
    w.u16(GDC_FORMAT_VERSION)
    w.u16(SEC_COUNT)
    w.raw(b'\x00' * 24)

    toc_offset = w.tell()
    for _ in range(SEC_COUNT):
        w.u32(0)
        w.u32(0)
        w.u64(0)

    section_writers = [
        (SEC_METADATA,       lambda: _write_metadata(w, config)),
        (SEC_ALLELES,        lambda: _write_alleles(w, config)),
        (SEC_GENE_USE,       lambda: _write_gene_use(w, config)),
        (SEC_TRIM_DISTS,     lambda: _write_trim_dists(w, config)),
        (SEC_NP_PARAMS,      lambda: _write_np_params(w, config)),
        (SEC_P_NUC_PROBS,    lambda: _write_p_nuc_probs(w, config)),
        (SEC_MUTATION_MODEL, lambda: _write_mutation_model(w, config, mutability, substitution)),
    ]

    entries = []
    for sec_id, writer_fn in section_writers:
        start = w.tell()
        writer_fn()
        size = w.tell() - start
        entries.append((sec_id, size, start))

    for i, (sec_id, size, offset) in enumerate(entries):
        entry_offset = toc_offset + i * 16
        w.patch_at(entry_offset, struct.pack('<IIQ', sec_id, size, offset))

    return w.getvalue()


def load_gdc(path: str) -> tuple:
    """
    Load a .gdc binary file.

    Returns:
        Tuple of (DataConfig, mutability_dict, substitution_dict).
        mutability/substitution are None if no mutation model was stored.
    """
    data = Path(path).read_bytes()
    r = _BinaryReader(data)

    # Validate header
    magic = r.raw(4)
    if magic != GDC_MAGIC:
        raise ValueError(f"Invalid GDC magic: {magic!r}")

    version = r.u16()
    if version > GDC_FORMAT_VERSION:
        raise ValueError(f"Unsupported GDC version: {version}")

    n_sections = r.u16()
    r.raw(24)  # skip reserved

    # Read section table
    sections = {}
    for _ in range(n_sections):
        sec_id = r.u32()
        size = r.u32()
        offset = r.u64()
        sections[sec_id] = (offset, size)

    config = DataConfig()
    mutability = None
    substitution = None

    # Process known sections
    if SEC_METADATA in sections:
        r.seek(sections[SEC_METADATA][0])
        _read_metadata(r, config)

    if SEC_ALLELES in sections:
        r.seek(sections[SEC_ALLELES][0])
        _read_alleles(r, config)

    if SEC_GENE_USE in sections:
        r.seek(sections[SEC_GENE_USE][0])
        _read_gene_use(r, config)

    if SEC_TRIM_DISTS in sections:
        r.seek(sections[SEC_TRIM_DISTS][0])
        _read_trim_dists(r, config)

    if SEC_NP_PARAMS in sections:
        r.seek(sections[SEC_NP_PARAMS][0])
        _read_np_params(r, config)

    if SEC_P_NUC_PROBS in sections:
        r.seek(sections[SEC_P_NUC_PROBS][0])
        _read_p_nuc_probs(r, config)

    if SEC_MUTATION_MODEL in sections:
        r.seek(sections[SEC_MUTATION_MODEL][0])
        mutability, substitution = _read_mutation_model(r)

    return config, mutability, substitution


def _load_s5f_tables(chain_category: str) -> tuple:
    """Load S5F tables from built-in package data."""
    from importlib import resources

    pkg = "GenAIRR.data.mutation_model_parameters"
    filename = "HH_S5F_META.pkl" if chain_category == "heavy" else "HKL_S5F_META.pkl"

    with resources.path(pkg, filename) as data_path:
        with open(data_path, "rb") as f:
            mutability, substitution, targeting = pickle.load(f)

    # Clean NaN values
    mutability = {
        k: (v if isinstance(v, (int, float)) and not math.isnan(v) else 0.0)
        for k, v in mutability.items()
    }

    # Convert DataFrame if needed
    if hasattr(substitution, "to_dict"):
        substitution = substitution.to_dict(orient="dict")
        substitution = {
            outer_key: {
                inner_key: inner_value
                for inner_key, inner_value in outer_dict.items()
                if isinstance(inner_value, (int, float)) and not math.isnan(inner_value)
            }
            for outer_key, outer_dict in substitution.items()
        }

    return mutability, substitution, targeting
