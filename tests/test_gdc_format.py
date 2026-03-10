"""
Tests for the GDC binary format: Python write/read round-trip and
cross-language compatibility (write in Python, read in C and vice versa).
"""

import os
import struct
import tempfile
import pytest
from pathlib import Path

from GenAIRR.dataconfig.gdc_io import (
    save_gdc, load_gdc,
    GDC_MAGIC, GDC_FORMAT_VERSION,
    _s5f_str_to_key, _s5f_key_to_str,
)
from GenAIRR.dataconfig.data_config import DataConfig
from GenAIRR.dataconfig.config_info import ConfigInfo
from GenAIRR.dataconfig.enums import ChainType, Species
from GenAIRR.alleles.allele import VAllele, DAllele, JAllele

from datetime import date


# ═════════════════════════════════════════════════════════════════════
# Fixtures
# ═════════════════════════════════════════════════════════════════════

@pytest.fixture
def tmp_gdc(tmp_path):
    """Return a path for a temporary .gdc file."""
    return str(tmp_path / "test.gdc")


@pytest.fixture
def sample_config():
    """Build a minimal DataConfig for testing."""
    config = DataConfig()
    config.name = "TEST_IGH"
    config.metadata = ConfigInfo(
        species=Species.HUMAN,
        chain_type=ChainType.BCR_HEAVY,
        reference_set="IMGT",
        last_updated=date(2026, 3, 8),
        has_d=True,
    )

    # V alleles
    v1 = VAllele("IGHV1-2*01", "ATGATGATGATGATGATG", 18, anchor_override=15)
    v2 = VAllele("IGHV1-2*02", "ATGATGATGATGATGATC", 18, anchor_override=15)
    config.v_alleles = {"IGHV1-2": [v1, v2]}

    # D alleles
    d1 = DAllele("IGHD1-1*01", "GATTACA", 7, anchor_override=0)
    config.d_alleles = {"IGHD1-1": [d1]}

    # J alleles
    j1 = JAllele("IGHJ4*02", "ACTACTTTGACTACTGG", 17, anchor_override=14)
    config.j_alleles = {"IGHJ4": [j1]}

    config.c_alleles = None

    # Gene use
    config.gene_use_dict = {
        "V": {"IGHV1-2": 0.7, "IGHV3-11": 0.3},
        "D": {"IGHD1-1": 1.0},
        "J": {"IGHJ4": 1.0},
    }

    # Trim dists
    config.trim_dicts = {
        "V_3": {
            "IGHV1": {
                "IGHV1-2": {i: (1.0 / (1 + i)) for i in range(11)},
            }
        },
        "D_5": {"IGHD1": {"IGHD1-1": {i: (1.0 / (1 + i)) for i in range(6)}}},
        "D_3": {"IGHD1": {"IGHD1-1": {i: (1.0 / (1 + i)) for i in range(6)}}},
        "J_5": {"IGHJ4": {"IGHJ4": {i: (1.0 / (1 + i)) for i in range(9)}}},
    }
    # Normalize
    for trim_key in config.trim_dicts:
        for fam in config.trim_dicts[trim_key]:
            for gene in config.trim_dicts[trim_key][fam]:
                d = config.trim_dicts[trim_key][fam][gene]
                total = sum(d.values())
                for k in d:
                    d[k] /= total

    # NP params
    config.NP_lengths = {
        "NP1": {i: (1.0 / (1 + i)) if i < 15 else 0.0 for i in range(21)},
        "NP2": {i: (1.0 / (1 + i)) if i < 15 else 0.0 for i in range(21)},
    }
    for key in config.NP_lengths:
        total = sum(config.NP_lengths[key].values())
        for k in config.NP_lengths[key]:
            config.NP_lengths[key][k] /= total

    config.NP_first_bases = {
        "NP1": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
        "NP2": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
    }

    config.NP_transitions = {}
    for key in ("NP1", "NP2"):
        config.NP_transitions[key] = {}
        for pos in range(3):
            config.NP_transitions[key][pos] = {
                fb: {tb: 0.25 for tb in "ACGT"} for fb in "ACGT"
            }

    config.p_nucleotide_length_probs = {0: 0.50, 1: 0.25, 2: 0.15, 3: 0.07, 4: 0.03}

    return config


@pytest.fixture
def sample_s5f_data():
    """Build minimal S5F test data."""
    mutability = {"AAAAA": 0.005, "ACGTA": 0.0023}
    substitution = {
        "AAAAA": {"C": 0.354, "G": 0.415, "T": 0.231},
        "ACGTA": {"A": 0.6, "T": 0.4},
    }
    return mutability, substitution


# ═════════════════════════════════════════════════════════════════════
# S5F key encoding tests
# ═════════════════════════════════════════════════════════════════════

class TestS5FEncoding:
    def test_key_roundtrip(self):
        for kmer in ("AAAAA", "ACGTA", "TTTTT", "NNNNN", "GATCN"):
            key = _s5f_str_to_key(kmer)
            assert 0 <= key < 3125
            assert _s5f_key_to_str(key) == kmer

    def test_known_keys(self):
        assert _s5f_str_to_key("AAAAA") == 0
        assert _s5f_str_to_key("NNNNN") == 3124
        # ACGTA: A=0, C=1, G=2, T=3, A=0 → 0*625 + 1*125 + 2*25 + 3*5 + 0 = 190
        assert _s5f_str_to_key("ACGTA") == 190


# ═════════════════════════════════════════════════════════════════════
# Round-trip tests
# ═════════════════════════════════════════════════════════════════════

class TestGDCRoundTrip:
    def test_metadata_roundtrip(self, tmp_gdc, sample_config):
        save_gdc(tmp_gdc, sample_config)
        loaded, _, _ = load_gdc(tmp_gdc)

        assert loaded.name == "TEST_IGH"
        assert loaded.metadata.chain_type == ChainType.BCR_HEAVY
        assert loaded.metadata.has_d is True
        assert loaded.metadata.species == Species.HUMAN
        assert loaded.metadata.reference_set == "IMGT"
        assert loaded.metadata.last_updated == date(2026, 3, 8)

    def test_alleles_roundtrip(self, tmp_gdc, sample_config):
        save_gdc(tmp_gdc, sample_config)
        loaded, _, _ = load_gdc(tmp_gdc)

        assert "IGHV1-2" in loaded.v_alleles
        assert len(loaded.v_alleles["IGHV1-2"]) == 2
        v1 = loaded.v_alleles["IGHV1-2"][0]
        assert v1.name == "IGHV1-2*01"
        assert v1.ungapped_seq == "ATGATGATGATGATGATG"
        assert v1.anchor == 15
        assert v1.gene == "IGHV1-2"
        assert v1.family == "IGHV1"

        assert "IGHD1-1" in loaded.d_alleles
        assert "IGHJ4" in loaded.j_alleles
        assert loaded.c_alleles is None

    def test_gene_use_roundtrip(self, tmp_gdc, sample_config):
        save_gdc(tmp_gdc, sample_config)
        loaded, _, _ = load_gdc(tmp_gdc)

        assert abs(loaded.gene_use_dict["V"]["IGHV1-2"] - 0.7) < 1e-10
        assert abs(loaded.gene_use_dict["V"]["IGHV3-11"] - 0.3) < 1e-10
        assert abs(loaded.gene_use_dict["J"]["IGHJ4"] - 1.0) < 1e-10

    def test_trim_dists_roundtrip(self, tmp_gdc, sample_config):
        save_gdc(tmp_gdc, sample_config)
        loaded, _, _ = load_gdc(tmp_gdc)

        orig = sample_config.trim_dicts["V_3"]["IGHV1"]["IGHV1-2"]
        loaded_d = loaded.trim_dicts["V_3"]["IGHV1"]["IGHV1-2"]

        for k in orig:
            assert abs(orig[k] - loaded_d[k]) < 1e-10

    def test_np_params_roundtrip(self, tmp_gdc, sample_config):
        save_gdc(tmp_gdc, sample_config)
        loaded, _, _ = load_gdc(tmp_gdc)

        for key in ("NP1", "NP2"):
            for i in range(21):
                orig = sample_config.NP_lengths[key][i]
                got = loaded.NP_lengths[key][i]
                assert abs(orig - got) < 1e-10, f"NP_lengths[{key}][{i}]"

            for base in "ACGT":
                assert abs(loaded.NP_first_bases[key][base] - 0.25) < 1e-10

            for pos in range(3):
                for fb in "ACGT":
                    for tb in "ACGT":
                        assert abs(loaded.NP_transitions[key][pos][fb][tb] - 0.25) < 1e-10

    def test_p_nuc_probs_roundtrip(self, tmp_gdc, sample_config):
        save_gdc(tmp_gdc, sample_config)
        loaded, _, _ = load_gdc(tmp_gdc)

        orig = sample_config.p_nucleotide_length_probs
        for k, v in orig.items():
            assert abs(loaded.p_nucleotide_length_probs[k] - v) < 1e-10

    def test_mutation_model_roundtrip(self, tmp_gdc, sample_config, sample_s5f_data):
        mut, sub = sample_s5f_data
        save_gdc(tmp_gdc, sample_config, mutability=mut, substitution=sub)
        loaded, loaded_mut, loaded_sub = load_gdc(tmp_gdc)

        assert loaded_mut is not None
        assert loaded_sub is not None

        assert abs(loaded_mut["AAAAA"] - 0.005) < 1e-10
        assert abs(loaded_mut["ACGTA"] - 0.0023) < 1e-10

        assert abs(loaded_sub["AAAAA"]["C"] - 0.354) < 1e-10
        assert abs(loaded_sub["AAAAA"]["G"] - 0.415) < 1e-10
        assert abs(loaded_sub["ACGTA"]["A"] - 0.6) < 1e-10

    def test_no_mutation_model(self, tmp_gdc, sample_config):
        save_gdc(tmp_gdc, sample_config)
        _, mut, sub = load_gdc(tmp_gdc)
        assert mut is None
        assert sub is None


# ═════════════════════════════════════════════════════════════════════
# Edge cases and error handling
# ═════════════════════════════════════════════════════════════════════

class TestGDCEdgeCases:
    def test_invalid_magic(self, tmp_gdc):
        with open(tmp_gdc, "wb") as f:
            f.write(b"NOPE" + b"\x00" * 28)
        with pytest.raises(ValueError, match="Invalid GDC magic"):
            load_gdc(tmp_gdc)

    def test_empty_config(self, tmp_gdc):
        config = DataConfig()
        config.name = "EMPTY"
        config.gene_use_dict = {}
        config.trim_dicts = {}
        config.NP_lengths = {}
        config.NP_first_bases = {}
        config.NP_transitions = {}
        config.p_nucleotide_length_probs = {0: 1.0}

        save_gdc(tmp_gdc, config)
        loaded, _, _ = load_gdc(tmp_gdc)
        assert loaded.name == "EMPTY"

    def test_file_size_reasonable(self, tmp_gdc, sample_config, sample_s5f_data):
        mut, sub = sample_s5f_data
        save_gdc(tmp_gdc, sample_config, mutability=mut, substitution=sub)
        size = os.path.getsize(tmp_gdc)
        # Should be small for test data (no full 3125 entries)
        assert size > 100
        assert size < 200000

    def test_binary_header_structure(self, tmp_gdc, sample_config):
        save_gdc(tmp_gdc, sample_config)
        with open(tmp_gdc, "rb") as f:
            magic = f.read(4)
            assert magic == GDC_MAGIC
            version = struct.unpack("<H", f.read(2))[0]
            assert version == GDC_FORMAT_VERSION
            n_sections = struct.unpack("<H", f.read(2))[0]
            assert n_sections == 7  # SEC_COUNT


# ═════════════════════════════════════════════════════════════════════
# Integration: load from real DataConfig pickle
# ═════════════════════════════════════════════════════════════════════

class TestGDCIntegration:
    def test_real_dataconfig_roundtrip(self, tmp_gdc):
        """Load a real builtin DataConfig, convert to .gdc, read back."""
        try:
            from GenAIRR.data import HUMAN_IGH_IMGT as config
        except Exception:
            pytest.skip("HUMAN_IGH_IMGT config not available")

        # Determine chain category for S5F
        chain_cat = "heavy" if config.metadata.chain_type == ChainType.BCR_HEAVY else "light"

        save_gdc(tmp_gdc, config, s5f_chain_category=chain_cat)

        loaded, mut, sub = load_gdc(tmp_gdc)

        # Metadata
        assert loaded.name == config.name
        assert loaded.metadata.chain_type == config.metadata.chain_type
        assert loaded.metadata.has_d == config.metadata.has_d

        # Allele counts should match
        orig_v = sum(len(v) for v in config.v_alleles.values())
        loaded_v = sum(len(v) for v in loaded.v_alleles.values())
        assert loaded_v == orig_v

        # Gene use keys should match
        assert set(loaded.gene_use_dict["V"].keys()) == set(config.gene_use_dict["V"].keys())

        # S5F should be loaded
        assert mut is not None
        assert sub is not None
        assert len(mut) == 3125

        # File size for full config with S5F
        size = os.path.getsize(tmp_gdc)
        print(f"Real HUMAN_IGH_IMGT .gdc file size: {size:,} bytes")
        assert size < 500_000  # should be well under 500 KB
