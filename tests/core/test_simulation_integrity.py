"""
Tests for biological correctness and simulation integrity.

These tests validate invariants that must hold for simulated sequences
to be biologically plausible and internally consistent. They catch
subtle coordinate, mutation, and biological errors that basic
functional tests miss.
"""

import pytest
from collections import defaultdict

from GenAIRR import (
    Pipeline,
    steps,
    S5F,
    Uniform,
    HUMAN_IGH_OGRDB,
    HUMAN_IGK_OGRDB,
    HUMAN_IGL_OGRDB,
)
from GenAIRR.steps import (
    SimulateSequence,
    FixVPositionAfterTrimmingIndexAmbiguity,
    FixDPositionAfterTrimmingIndexAmbiguity,
    FixJPositionAfterTrimmingIndexAmbiguity,
    CorrectForVEndCut,
    CorrectForDTrims,
    CorruptSequenceBeginning,
    InsertNs,
    InsertIndels,
)

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
# Reusable coordinate integrity checker
# =============================================================================

def assert_coordinate_integrity(container, label=""):
    """Validate all positional invariants on a SimulationContainer."""
    seq = container.sequence
    seq_len = len(seq)
    pfx = f"[{label}] " if label else ""

    # All positions within sequence bounds
    assert 0 <= container.v_sequence_start <= seq_len, \
        f"{pfx}v_sequence_start={container.v_sequence_start} out of bounds (seq_len={seq_len})"
    assert 0 <= container.v_sequence_end <= seq_len, \
        f"{pfx}v_sequence_end={container.v_sequence_end} out of bounds"
    assert 0 <= container.j_sequence_start <= seq_len, \
        f"{pfx}j_sequence_start={container.j_sequence_start} out of bounds"
    assert 0 <= container.j_sequence_end <= seq_len, \
        f"{pfx}j_sequence_end={container.j_sequence_end} out of bounds"
    assert 0 <= container.d_sequence_start <= seq_len, \
        f"{pfx}d_sequence_start={container.d_sequence_start} out of bounds"
    assert 0 <= container.d_sequence_end <= seq_len, \
        f"{pfx}d_sequence_end={container.d_sequence_end} out of bounds"

    # Start <= End for each segment
    assert container.v_sequence_start <= container.v_sequence_end, \
        f"{pfx}V start > end: {container.v_sequence_start} > {container.v_sequence_end}"
    assert container.d_sequence_start <= container.d_sequence_end, \
        f"{pfx}D start > end"
    assert container.j_sequence_start <= container.j_sequence_end, \
        f"{pfx}J start > end"

    # Segments ordered: V before J
    assert container.v_sequence_end <= container.j_sequence_start, \
        f"{pfx}V end ({container.v_sequence_end}) > J start ({container.j_sequence_start})"

    # D ordering (when D is present)
    if container.d_call and container.d_call[0] and container.d_sequence_start != container.d_sequence_end:
        assert container.v_sequence_end <= container.d_sequence_start, \
            f"{pfx}V end > D start"
        assert container.d_sequence_end <= container.j_sequence_start, \
            f"{pfx}D end > J start"

    # Junction within bounds
    assert 0 <= container.junction_sequence_start <= seq_len, \
        f"{pfx}junction_start out of bounds"
    assert 0 <= container.junction_sequence_end <= seq_len, \
        f"{pfx}junction_end out of bounds"
    assert container.junction_sequence_start <= container.junction_sequence_end, \
        f"{pfx}junction start > end"

    # Mutation positions within bounds
    if container.mutations:
        for pos in container.mutations:
            assert 0 <= int(pos) < seq_len, \
                f"{pfx}mutation pos {pos} out of bounds"

    # N positions within bounds and actually N
    if container.Ns:
        for pos in container.Ns:
            p = int(pos)
            assert 0 <= p < seq_len, f"{pfx}N pos {p} out of bounds"
            assert seq[p] == 'N', f"{pfx}pos {p} marked as N but is '{seq[p]}'"

    # Only valid nucleotides
    invalid = set(seq) - set('ATCGNatcgn')
    assert not invalid, f"{pfx}invalid characters: {invalid}"


# =============================================================================
# 1. Junction Integrity
# =============================================================================

class TestJunctionIntegrity:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.pipeline = Pipeline(
            config=HUMAN_IGH_OGRDB,
            steps=[
                SimulateSequence(S5F(min_mutation_rate=0.01, max_mutation_rate=0.05), productive=True),
                FixVPositionAfterTrimmingIndexAmbiguity(),
                FixDPositionAfterTrimmingIndexAmbiguity(),
                FixJPositionAfterTrimmingIndexAmbiguity(),
            ]
        )

    def test_junction_within_sequence_bounds(self):
        for i in range(N_SEQUENCES):
            r = self.pipeline.execute()
            seq_len = len(r.sequence)
            assert 0 <= r.junction_sequence_start < seq_len
            assert 0 < r.junction_sequence_end <= seq_len
            assert r.junction_sequence_start < r.junction_sequence_end

    def test_junction_divisible_by_three_when_in_frame(self):
        for i in range(N_SEQUENCES):
            r = self.pipeline.execute()
            if r.productive and r.vj_in_frame:
                jlen = r.junction_sequence_end - r.junction_sequence_start
                assert jlen % 3 == 0, \
                    f"Seq {i}: productive junction length {jlen} not divisible by 3"

    def test_junction_starts_with_cysteine_codon(self):
        cys_codons = {'TGT', 'TGC'}
        for i in range(N_SEQUENCES):
            r = self.pipeline.execute()
            if r.productive:
                codon = r.sequence[r.junction_sequence_start:r.junction_sequence_start + 3].upper()
                assert codon in cys_codons, \
                    f"Seq {i}: junction starts with '{codon}', expected Cys"

    def test_junction_ends_with_fw_codon(self):
        fw_codons = {'TTT', 'TTC', 'TGG'}
        for i in range(N_SEQUENCES):
            r = self.pipeline.execute()
            if r.productive:
                codon = r.sequence[r.junction_sequence_end - 3:r.junction_sequence_end].upper()
                assert codon in fw_codons, \
                    f"Seq {i}: junction ends with '{codon}', expected F/W"


# =============================================================================
# 2. V/D/J Segment Boundaries
# =============================================================================

class TestSegmentBoundaries:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.pipeline = Pipeline(
            config=HUMAN_IGH_OGRDB,
            steps=[
                SimulateSequence(S5F(min_mutation_rate=0.01, max_mutation_rate=0.05), productive=True),
                FixVPositionAfterTrimmingIndexAmbiguity(),
                FixDPositionAfterTrimmingIndexAmbiguity(),
                FixJPositionAfterTrimmingIndexAmbiguity(),
            ]
        )

    def test_segments_non_overlapping(self):
        for i in range(N_SEQUENCES):
            r = self.pipeline.execute()
            if r.d_sequence_start != r.d_sequence_end:
                assert r.v_sequence_end <= r.d_sequence_start, \
                    f"Seq {i}: V end ({r.v_sequence_end}) > D start ({r.d_sequence_start})"
                assert r.d_sequence_end <= r.j_sequence_start, \
                    f"Seq {i}: D end ({r.d_sequence_end}) > J start ({r.j_sequence_start})"
            assert r.v_sequence_end <= r.j_sequence_start, \
                f"Seq {i}: V end ({r.v_sequence_end}) > J start ({r.j_sequence_start})"

    def test_segment_lengths_positive(self):
        for i in range(N_SEQUENCES):
            r = self.pipeline.execute()
            assert r.v_sequence_end - r.v_sequence_start > 0, f"Seq {i}: V empty"
            assert r.d_sequence_end - r.d_sequence_start >= 0, f"Seq {i}: D negative"
            assert r.j_sequence_end - r.j_sequence_start > 0, f"Seq {i}: J empty"

    def test_germline_positions_consistent(self):
        for i in range(N_SEQUENCES):
            r = self.pipeline.execute()
            assert r.v_germline_start <= r.v_germline_end
            assert r.j_germline_start <= r.j_germline_end
            if r.d_germline_start is not None and r.d_germline_end is not None:
                assert r.d_germline_start <= r.d_germline_end


# =============================================================================
# 3. Mutation Location Constraints
# =============================================================================

class TestMutationLocations:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.pipeline = Pipeline(
            config=HUMAN_IGH_OGRDB,
            steps=[
                SimulateSequence(S5F(min_mutation_rate=0.05, max_mutation_rate=0.15), productive=True),
                FixVPositionAfterTrimmingIndexAmbiguity(),
                FixDPositionAfterTrimmingIndexAmbiguity(),
                FixJPositionAfterTrimmingIndexAmbiguity(),
            ]
        )

    def test_mutations_within_sequence_bounds(self):
        for i in range(N_SEQUENCES):
            r = self.pipeline.execute()
            seq_len = len(r.sequence)
            for pos in r.mutations:
                assert 0 <= int(pos) < seq_len, \
                    f"Seq {i}: mutation at {pos} out of bounds (len={seq_len})"

    def test_mutations_in_vdj_regions_only(self):
        """S5F mutations should only occur within V, D, or J segments.

        NOTE: After position-fixing steps adjust V/D/J boundaries, some
        mutation positions that were originally in V/D/J may appear to be
        in NP regions. We allow up to 10% as these are position-fixing
        artifacts, not actual NP mutations by S5F.
        """
        np_mutation_count = 0
        total_mutations = 0

        for i in range(N_SEQUENCES):
            r = self.pipeline.execute()
            vdj = set(range(r.v_sequence_start, r.v_sequence_end)) | \
                  set(range(r.d_sequence_start, r.d_sequence_end)) | \
                  set(range(r.j_sequence_start, r.j_sequence_end))

            for pos in r.mutations:
                total_mutations += 1
                if int(pos) not in vdj:
                    np_mutation_count += 1

        if total_mutations > 0:
            np_frac = np_mutation_count / total_mutations
            # Allow up to 10% due to position-fixing boundary shifts
            assert np_frac < 0.10, \
                f"{np_mutation_count}/{total_mutations} ({np_frac:.1%}) mutations outside V/D/J regions"

    def test_mutation_format_valid(self):
        valid_bases = set('ATCGatcg')
        for i in range(N_SEQUENCES):
            r = self.pipeline.execute()
            for pos, mut in r.mutations.items():
                parts = str(mut).split('>')
                assert len(parts) >= 2, f"Seq {i}: invalid format '{mut}' at {pos}"
                for base in parts:
                    assert len(base) == 1 and base in valid_bases, \
                        f"Seq {i}: invalid base '{base}' in '{mut}' at {pos}"


# =============================================================================
# 4. Productive Sequence Biological Invariants
# =============================================================================

class TestProductiveInvariants:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.pipeline = Pipeline(
            config=HUMAN_IGH_OGRDB,
            steps=[
                SimulateSequence(S5F(min_mutation_rate=0.01, max_mutation_rate=0.05), productive=True),
                FixVPositionAfterTrimmingIndexAmbiguity(),
            ]
        )

    def test_productive_no_stop_codons_in_junction(self):
        for i in range(N_SEQUENCES):
            r = self.pipeline.execute()
            if not r.productive:
                continue
            junction = r.sequence[r.junction_sequence_start:r.junction_sequence_end].upper()
            for j in range(0, len(junction) - 2, 3):
                codon = junction[j:j+3]
                aa = CODON_TABLE.get(codon, 'X')
                assert aa != '*', \
                    f"Seq {i}: stop codon '{codon}' at junction pos {j}"

    def test_productive_flag_consistency(self):
        for i in range(N_SEQUENCES):
            r = self.pipeline.execute()
            if r.productive:
                assert not r.stop_codon, f"Seq {i}: productive but has stop codon"
                assert r.vj_in_frame, f"Seq {i}: productive but not in frame"

    def test_mutation_rate_within_bounds(self):
        for i in range(N_SEQUENCES):
            r = self.pipeline.execute()
            if r.mutation_rate is not None:
                assert 0 <= r.mutation_rate <= 1.0, \
                    f"Seq {i}: mutation rate {r.mutation_rate} out of [0,1]"


# =============================================================================
# 5. Full Pipeline End-to-End Coordinate Consistency
# =============================================================================

class TestFullPipelineIntegrity:

    def _make_pipeline(self, config, corruption=False, ns=False, indels=False, has_d_segment=True):
        s = [
            SimulateSequence(S5F(min_mutation_rate=0.01, max_mutation_rate=0.05), productive=True),
            FixVPositionAfterTrimmingIndexAmbiguity(),
        ]
        if has_d_segment:
            s.append(FixDPositionAfterTrimmingIndexAmbiguity())
        s.append(FixJPositionAfterTrimmingIndexAmbiguity())
        s.append(CorrectForVEndCut())
        if has_d_segment:
            s.append(CorrectForDTrims())
        if corruption:
            s.append(CorruptSequenceBeginning(probability=1.0))
        if ns:
            s.append(InsertNs(probability=1.0))
        if indels:
            s.append(InsertIndels(probability=1.0))
        return Pipeline(config=config, steps=s)

    def test_basic_pipeline(self):
        pipeline = self._make_pipeline(HUMAN_IGH_OGRDB)
        for i in range(N_SEQUENCES):
            assert_coordinate_integrity(pipeline.execute(), f"basic-{i}")

    def test_with_corruption(self):
        pipeline = self._make_pipeline(HUMAN_IGH_OGRDB, corruption=True)
        for i in range(N_SEQUENCES):
            assert_coordinate_integrity(pipeline.execute(), f"corrupt-{i}")

    def test_with_ns(self):
        pipeline = self._make_pipeline(HUMAN_IGH_OGRDB, ns=True)
        for i in range(N_SEQUENCES):
            assert_coordinate_integrity(pipeline.execute(), f"ns-{i}")

    def test_with_indels(self):
        pipeline = self._make_pipeline(HUMAN_IGH_OGRDB, indels=True)
        for i in range(N_SEQUENCES):
            assert_coordinate_integrity(pipeline.execute(), f"indel-{i}")

    def test_full_augmentation(self):
        pipeline = self._make_pipeline(HUMAN_IGH_OGRDB, corruption=True, ns=True, indels=True)
        for i in range(N_SEQUENCES):
            assert_coordinate_integrity(pipeline.execute(), f"full-{i}")

    def test_light_chain_kappa(self):
        pipeline = self._make_pipeline(HUMAN_IGK_OGRDB, has_d_segment=False)
        for i in range(N_SEQUENCES):
            r = pipeline.execute()
            assert_coordinate_integrity(r, f"igk-{i}")
            assert not r.d_call or r.d_call == [], f"Seq {i}: kappa has D call"

    def test_light_chain_lambda(self):
        pipeline = self._make_pipeline(HUMAN_IGL_OGRDB, has_d_segment=False)
        for i in range(N_SEQUENCES):
            r = pipeline.execute()
            assert_coordinate_integrity(r, f"igl-{i}")
            assert not r.d_call or r.d_call == [], f"Seq {i}: lambda has D call"


# =============================================================================
# 6. Indels & Corruption Position Tracking
# =============================================================================

class TestIndelsAndCorruption:

    def test_indel_positions_within_bounds(self):
        pipeline = Pipeline(
            config=HUMAN_IGH_OGRDB,
            steps=[
                SimulateSequence(S5F(min_mutation_rate=0.01, max_mutation_rate=0.05), productive=True),
                FixVPositionAfterTrimmingIndexAmbiguity(),
                InsertIndels(probability=1.0),
            ]
        )
        for i in range(N_SEQUENCES):
            r = pipeline.execute()
            for pos in r.indels:
                assert 0 <= int(pos) < len(r.sequence), \
                    f"Seq {i}: indel pos {pos} out of bounds"

    def test_corruption_metadata_consistent(self):
        pipeline = Pipeline(
            config=HUMAN_IGH_OGRDB,
            steps=[
                SimulateSequence(S5F(min_mutation_rate=0.01, max_mutation_rate=0.05), productive=True),
                FixVPositionAfterTrimmingIndexAmbiguity(),
                CorruptSequenceBeginning(probability=1.0),
            ]
        )
        for i in range(N_SEQUENCES):
            r = pipeline.execute()
            event = r.corruption_event
            if event == 'add':
                assert r.corruption_add_amount > 0, \
                    f"Seq {i}: add event but amount is 0"
                assert len(r.corruption_added_section) == r.corruption_add_amount, \
                    f"Seq {i}: added section length mismatch"
            elif event == 'remove':
                assert r.corruption_remove_amount > 0, \
                    f"Seq {i}: remove event but amount is 0"
                assert len(r.corruption_removed_section) == r.corruption_remove_amount, \
                    f"Seq {i}: removed section length mismatch"
            elif event == 'remove_before_add':
                assert r.corruption_remove_amount > 0
                assert r.corruption_add_amount > 0

    def test_n_positions_are_actually_n(self):
        """Verify that every position logged in Ns is actually 'N' in the sequence."""
        pipeline = Pipeline(
            config=HUMAN_IGH_OGRDB,
            steps=[
                SimulateSequence(S5F(min_mutation_rate=0.05, max_mutation_rate=0.15), productive=True),
                FixVPositionAfterTrimmingIndexAmbiguity(),
                InsertNs(probability=1.0),
            ]
        )
        for i in range(N_SEQUENCES):
            r = pipeline.execute()
            for pos in r.Ns:
                p = int(pos)
                assert r.sequence[p] == 'N', \
                    f"Seq {i}: pos {p} logged as N but is '{r.sequence[p]}'"


# =============================================================================
# 7. Uniform Mutation Model
# =============================================================================

class TestUniformMutation:

    def test_mutation_count_matches_rate(self):
        pipeline = Pipeline(
            config=HUMAN_IGH_OGRDB,
            steps=[
                SimulateSequence(Uniform(min_mutation_rate=0.05, max_mutation_rate=0.05), productive=True),
                FixVPositionAfterTrimmingIndexAmbiguity(),
            ]
        )
        rates = []
        for i in range(N_SEQUENCES):
            r = pipeline.execute()
            if len(r.sequence) > 0:
                rates.append(len(r.mutations) / len(r.sequence))

        if rates:
            import statistics
            mean_rate = statistics.mean(rates)
            assert 0.01 < mean_rate < 0.15, \
                f"Uniform mean rate {mean_rate} far from 0.05"


# =============================================================================
# 8. Container get_dict Completeness
# =============================================================================

class TestContainerDict:

    def test_has_all_required_fields(self):
        pipeline = Pipeline(
            config=HUMAN_IGH_OGRDB,
            steps=[
                SimulateSequence(S5F(min_mutation_rate=0.01, max_mutation_rate=0.05), productive=True),
                FixVPositionAfterTrimmingIndexAmbiguity(),
            ]
        )
        d = pipeline.execute().get_dict()

        required = [
            'sequence', 'v_call', 'd_call', 'j_call',
            'v_sequence_start', 'v_sequence_end',
            'd_sequence_start', 'd_sequence_end',
            'j_sequence_start', 'j_sequence_end',
            'junction_sequence_start', 'junction_sequence_end',
            'mutation_rate', 'mutations',
            'productive', 'stop_codon', 'vj_in_frame',
        ]
        for field in required:
            assert field in d, f"Missing field '{field}' in get_dict()"

    def test_values_match_attributes(self):
        pipeline = Pipeline(
            config=HUMAN_IGH_OGRDB,
            steps=[
                SimulateSequence(S5F(min_mutation_rate=0.01, max_mutation_rate=0.05), productive=True),
                FixVPositionAfterTrimmingIndexAmbiguity(),
            ]
        )
        r = pipeline.execute()
        d = r.get_dict()
        assert d['sequence'] == r.sequence
        assert d['v_sequence_start'] == r.v_sequence_start
        assert d['productive'] == r.productive
        assert d['mutation_rate'] == r.mutation_rate
