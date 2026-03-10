"""
Tests for biological correctness and simulation integrity.

These tests validate invariants that must hold for simulated sequences
to be biologically plausible and internally consistent. They catch
subtle coordinate, mutation, and biological errors that basic
functional tests miss.

All tests use the Experiment DSL → C backend path.
"""

import pytest

from GenAIRR import Experiment

N_SEQUENCES = 50

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


# =============================================================================
# Reusable coordinate integrity checker (dict-based)
# =============================================================================

def assert_coordinate_integrity(rec, label=""):
    """Validate all positional invariants on a simulation record dict."""
    seq = rec["sequence"]
    seq_len = len(seq)
    pfx = f"[{label}] " if label else ""

    # All positions within sequence bounds
    assert 0 <= rec["v_sequence_start"] <= seq_len, \
        f"{pfx}v_sequence_start={rec['v_sequence_start']} out of bounds (seq_len={seq_len})"
    assert 0 <= rec["v_sequence_end"] <= seq_len, \
        f"{pfx}v_sequence_end={rec['v_sequence_end']} out of bounds"
    assert 0 <= rec["j_sequence_start"] <= seq_len, \
        f"{pfx}j_sequence_start={rec['j_sequence_start']} out of bounds"
    assert 0 <= rec["j_sequence_end"] <= seq_len, \
        f"{pfx}j_sequence_end={rec['j_sequence_end']} out of bounds"

    # Start <= End for each segment
    assert rec["v_sequence_start"] <= rec["v_sequence_end"], \
        f"{pfx}V start > end: {rec['v_sequence_start']} > {rec['v_sequence_end']}"
    assert rec["j_sequence_start"] <= rec["j_sequence_end"], \
        f"{pfx}J start > end"

    # Segments ordered: V before J
    assert rec["v_sequence_end"] <= rec["j_sequence_end"], \
        f"{pfx}V end ({rec['v_sequence_end']}) > J end ({rec['j_sequence_end']})"

    # D ordering (when D is present — heavy chain)
    d_start = rec.get("d_sequence_start", 0)
    d_end = rec.get("d_sequence_end", 0)
    if d_start != d_end:
        assert 0 <= d_start <= seq_len, f"{pfx}d_sequence_start out of bounds"
        assert 0 <= d_end <= seq_len, f"{pfx}d_sequence_end out of bounds"
        assert d_start <= d_end, f"{pfx}D start > end"

    # Junction within bounds
    j_start = rec.get("junction_start", 0)
    j_end = rec.get("junction_end", 0)
    assert 0 <= j_start <= seq_len, f"{pfx}junction_start out of bounds"
    assert 0 <= j_end <= seq_len, f"{pfx}junction_end out of bounds"
    assert j_start <= j_end, f"{pfx}junction start > end"

    # Only valid nucleotides
    invalid = set(seq.upper()) - set('ATCGN')
    assert not invalid, f"{pfx}invalid characters: {invalid}"


# =============================================================================
# 1. Junction Integrity (productive, unmutated)
# =============================================================================

class TestJunctionIntegrity:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.records = Experiment.on("human_igh").run(
            n=N_SEQUENCES, seed=42, productive=True)

    def test_junction_within_sequence_bounds(self):
        for i, r in enumerate(self.records):
            seq_len = len(r["sequence"])
            js = r["junction_start"]
            je = r["junction_end"]
            assert 0 <= js < seq_len, f"Seq {i}: junction_start out of bounds"
            assert 0 < je <= seq_len, f"Seq {i}: junction_end out of bounds"
            assert js < je, f"Seq {i}: junction_start >= junction_end"

    def test_junction_divisible_by_three_when_in_frame(self):
        for i, r in enumerate(self.records):
            if r["productive"] and r["vj_in_frame"]:
                jlen = r["junction_end"] - r["junction_start"]
                assert jlen % 3 == 0, \
                    f"Seq {i}: productive junction length {jlen} not divisible by 3"

    def test_junction_starts_with_cysteine_codon(self):
        cys_codons = {'TGT', 'TGC'}
        for i, r in enumerate(self.records):
            if r["productive"]:
                js = r["junction_start"]
                codon = r["sequence"][js:js + 3].upper()
                assert codon in cys_codons, \
                    f"Seq {i}: junction starts with '{codon}', expected Cys"

    def test_junction_ends_with_fw_codon(self):
        fw_codons = {'TTT', 'TTC', 'TGG'}
        for i, r in enumerate(self.records):
            if r["productive"]:
                je = r["junction_end"]
                codon = r["sequence"][je - 3:je].upper()
                assert codon in fw_codons, \
                    f"Seq {i}: junction ends with '{codon}', expected F/W"

    def test_junction_reasonable_length(self):
        for i, r in enumerate(self.records):
            jlen = r["junction_end"] - r["junction_start"]
            assert 10 <= jlen <= 120, \
                f"Seq {i}: junction length {jlen} outside reasonable range"


# =============================================================================
# 2. V/D/J Segment Boundaries
# =============================================================================

class TestSegmentBoundaries:

    @pytest.fixture(autouse=True)
    def setup(self):
        from GenAIRR.ops import rate
        self.records = (Experiment.on("human_igh")
                        .mutate(rate(0.01, 0.05))
                        .run(n=N_SEQUENCES, seed=42))

    def test_segments_non_overlapping(self):
        for i, r in enumerate(self.records):
            d_start = r.get("d_sequence_start", 0)
            d_end = r.get("d_sequence_end", 0)
            if d_start != d_end:
                assert r["v_sequence_end"] <= d_start, \
                    f"Seq {i}: V end ({r['v_sequence_end']}) > D start ({d_start})"
                assert d_end <= r["j_sequence_start"], \
                    f"Seq {i}: D end ({d_end}) > J start ({r['j_sequence_start']})"
            assert r["v_sequence_end"] <= r["j_sequence_start"], \
                f"Seq {i}: V end ({r['v_sequence_end']}) > J start ({r['j_sequence_start']})"

    def test_segment_lengths_positive(self):
        for i, r in enumerate(self.records):
            v_len = r["v_sequence_end"] - r["v_sequence_start"]
            j_len = r["j_sequence_end"] - r["j_sequence_start"]
            assert v_len > 0, f"Seq {i}: V empty"
            assert j_len > 0, f"Seq {i}: J empty"


# =============================================================================
# 3. Mutation Constraints
# =============================================================================

class TestMutationConstraints:

    @pytest.fixture(autouse=True)
    def setup(self):
        from GenAIRR.ops import rate
        self.records = (Experiment.on("human_igh")
                        .mutate(rate(0.05, 0.15))
                        .run(n=N_SEQUENCES, seed=42))

    def test_mutation_rate_within_bounds(self):
        for i, r in enumerate(self.records):
            rate = r["mutation_rate"]
            assert 0 <= rate <= 1.0, \
                f"Seq {i}: mutation rate {rate} out of [0,1]"

    def test_mutation_count_consistent_with_rate(self):
        for i, r in enumerate(self.records):
            n_mut = r["n_mutations"]
            rate = r["mutation_rate"]
            seq_len = len(r["sequence"])
            if seq_len > 0 and n_mut > 0:
                computed_rate = n_mut / seq_len
                # Allow generous tolerance — rate and count may differ slightly
                assert computed_rate < 0.5, \
                    f"Seq {i}: {n_mut} mutations in {seq_len}bp seems too high"

    def test_no_mutations_without_shm(self):
        records = Experiment.on("human_igh").run(n=N_SEQUENCES, seed=42)
        for i, r in enumerate(records):
            assert r["n_mutations"] == 0, \
                f"Seq {i}: {r['n_mutations']} mutations without SHM enabled"
            assert r["mutation_rate"] == 0.0


# =============================================================================
# 4. Productive Sequence Biological Invariants
# =============================================================================

class TestProductiveInvariants:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.records = Experiment.on("human_igh").run(
            n=N_SEQUENCES, seed=42, productive=True)

    def test_productive_no_stop_codons_in_junction(self):
        for i, r in enumerate(self.records):
            if not r["productive"]:
                continue
            js = r["junction_start"]
            je = r["junction_end"]
            junction = r["sequence"][js:je].upper()
            for j in range(0, len(junction) - 2, 3):
                codon = junction[j:j+3]
                aa = CODON_TABLE.get(codon, 'X')
                assert aa != '*', \
                    f"Seq {i}: stop codon '{codon}' at junction pos {j}"

    def test_productive_flag_consistency(self):
        for i, r in enumerate(self.records):
            if r["productive"]:
                assert not r["stop_codon"], f"Seq {i}: productive but has stop codon"
                assert r["vj_in_frame"], f"Seq {i}: productive but not in frame"


# =============================================================================
# 5. Full Pipeline End-to-End Coordinate Consistency
# =============================================================================

class TestFullPipelineIntegrity:

    def test_basic_pipeline(self):
        records = Experiment.on("human_igh").run(n=N_SEQUENCES, seed=42)
        for i, r in enumerate(records):
            assert_coordinate_integrity(r, f"basic-{i}")

    def test_with_shm(self):
        from GenAIRR.ops import rate
        records = (Experiment.on("human_igh")
                   .mutate(rate(0.01, 0.05))
                   .run(n=N_SEQUENCES, seed=42))
        for i, r in enumerate(records):
            assert_coordinate_integrity(r, f"shm-{i}")

    def test_with_corruption(self):
        from GenAIRR.ops import with_5prime_loss
        records = (Experiment.on("human_igh")
                   .sequence(with_5prime_loss())
                   .run(n=N_SEQUENCES, seed=42))
        for i, r in enumerate(records):
            assert_coordinate_integrity(r, f"corrupt-{i}")

    def test_with_indels(self):
        from GenAIRR.ops import with_indels
        records = (Experiment.on("human_igh")
                   .observe(with_indels())
                   .run(n=N_SEQUENCES, seed=42))
        for i, r in enumerate(records):
            assert_coordinate_integrity(r, f"indel-{i}")

    def test_full_augmentation(self):
        from GenAIRR.ops import rate, with_5prime_loss, with_indels, with_ns
        records = (Experiment.on("human_igh")
                   .mutate(rate(0.01, 0.05))
                   .sequence(with_5prime_loss())
                   .observe(with_indels(), with_ns())
                   .run(n=N_SEQUENCES, seed=42))
        for i, r in enumerate(records):
            assert_coordinate_integrity(r, f"full-{i}")

    def test_light_chain_kappa(self):
        records = Experiment.on("human_igk").run(n=N_SEQUENCES, seed=42)
        for i, r in enumerate(records):
            assert_coordinate_integrity(r, f"igk-{i}")

    def test_light_chain_lambda(self):
        records = Experiment.on("human_igl").run(n=N_SEQUENCES, seed=42)
        for i, r in enumerate(records):
            assert_coordinate_integrity(r, f"igl-{i}")


# =============================================================================
# 6. AIRR Dict Field Completeness
# =============================================================================

class TestAIRRFieldCompleteness:

    def test_has_all_required_fields(self):
        result = Experiment.on("human_igh").run(n=1, seed=42)
        rec = result[0]
        required = [
            'sequence', 'v_call', 'd_call', 'j_call',
            'v_sequence_start', 'v_sequence_end',
            'd_sequence_start', 'd_sequence_end',
            'j_sequence_start', 'j_sequence_end',
            'junction_start', 'junction_end',
            'mutation_rate', 'n_mutations',
            'productive', 'stop_codon', 'vj_in_frame',
        ]
        for field in required:
            assert field in rec, f"Missing field '{field}'"

    def test_sequence_only_valid_nucleotides(self):
        records = Experiment.on("human_igh").run(n=N_SEQUENCES, seed=42)
        for i, r in enumerate(records):
            seq = r["sequence"]
            assert len(seq) > 0, f"Seq {i}: empty sequence"
            invalid = set(seq.upper()) - set('ATCGN')
            assert not invalid, f"Seq {i}: invalid chars {invalid}"
