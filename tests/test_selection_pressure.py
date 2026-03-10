"""Tests for IMGT region utilities used by selection pressure and clonal simulation."""

import pytest

from GenAIRR import Experiment
from GenAIRR.data import HUMAN_IGH_OGRDB
from GenAIRR.utilities.imgt_regions import (
    gapped_to_ungapped_pos,
    compute_v_region_boundaries,
    classify_position,
    IMGT_GAPPED_BOUNDARIES,
)
from GenAIRR.utilities.misc import translate


# ---------------------------------------------------------------------------
# IMGT region mapper tests
# ---------------------------------------------------------------------------

class TestGappedToUngapped:
    """Test the gapped-to-ungapped position mapping."""

    def test_no_gaps(self):
        seq = "atcgatcg"
        assert gapped_to_ungapped_pos(seq, 0) == 0
        assert gapped_to_ungapped_pos(seq, 4) == 4
        assert gapped_to_ungapped_pos(seq, 8) == 8

    def test_all_gaps(self):
        seq = "........"
        assert gapped_to_ungapped_pos(seq, 0) == 0
        assert gapped_to_ungapped_pos(seq, 4) == 0
        assert gapped_to_ungapped_pos(seq, 8) == 0

    def test_mixed_gaps(self):
        seq = "atc...gatcg"
        assert gapped_to_ungapped_pos(seq, 3) == 3
        assert gapped_to_ungapped_pos(seq, 6) == 3
        assert gapped_to_ungapped_pos(seq, 7) == 4

    def test_pos_beyond_length(self):
        seq = "atcg"
        assert gapped_to_ungapped_pos(seq, 10) == 4


class TestComputeVRegionBoundaries:
    """Test boundary computation for real V alleles."""

    def test_returns_all_five_regions(self):
        allele_list = list(HUMAN_IGH_OGRDB.v_alleles.values())[0]
        allele = allele_list[0]
        boundaries = compute_v_region_boundaries(allele)
        assert set(boundaries.keys()) == {'FWR1', 'CDR1', 'FWR2', 'CDR2', 'FWR3'}

    def test_boundaries_non_overlapping(self):
        for allele_list in HUMAN_IGH_OGRDB.v_alleles.values():
            for allele in allele_list:
                b = compute_v_region_boundaries(allele)
                assert b['FWR1'][0] == 0
                assert b['FWR1'][1] == b['CDR1'][0]
                assert b['CDR1'][1] == b['FWR2'][0]
                assert b['FWR2'][1] == b['CDR2'][0]
                assert b['CDR2'][1] == b['FWR3'][0]
                assert b['FWR3'][1] <= allele.ungapped_len

    def test_fwr3_end_near_anchor(self):
        allele_list = list(HUMAN_IGH_OGRDB.v_alleles.values())[0]
        allele = allele_list[0]
        b = compute_v_region_boundaries(allele)
        assert abs(b['FWR3'][1] - allele.anchor) <= 3


class TestClassifyPosition:
    """Test position classification into IMGT regions."""

    def test_cdr3_takes_priority(self):
        v_boundaries = {'FWR1': (0, 78), 'CDR1': (78, 100), 'FWR2': (100, 150),
                        'CDR2': (150, 180), 'FWR3': (180, 290)}
        result = classify_position(
            pos=300, v_seq_start=0, v_seq_end=290,
            v_boundaries=v_boundaries,
            junction_start=285, junction_end=320, j_seq_end=350,
        )
        assert result == 'CDR3'

    def test_fwr4(self):
        result = classify_position(
            pos=325, v_seq_start=0, v_seq_end=290,
            v_boundaries={},
            junction_start=285, junction_end=320, j_seq_end=350,
        )
        assert result == 'FWR4'

    def test_v_region_classification(self):
        v_boundaries = {'FWR1': (0, 78), 'CDR1': (78, 100), 'FWR2': (100, 150),
                        'CDR2': (150, 180), 'FWR3': (180, 290)}
        assert classify_position(10, 0, 290, v_boundaries, 285, 320, 350) == 'FWR1'
        assert classify_position(85, 0, 290, v_boundaries, 285, 320, 350) == 'CDR1'
        assert classify_position(120, 0, 290, v_boundaries, 285, 320, 350) == 'FWR2'
        assert classify_position(160, 0, 290, v_boundaries, 285, 320, 350) == 'CDR2'
        assert classify_position(200, 0, 290, v_boundaries, 285, 320, 350) == 'FWR3'

    def test_np_region(self):
        result = classify_position(
            pos=400, v_seq_start=0, v_seq_end=290,
            v_boundaries={},
            junction_start=285, junction_end=320, j_seq_end=350,
        )
        assert result == 'NP'

    def test_v_seq_start_offset(self):
        v_boundaries = {'FWR1': (0, 78)}
        result = classify_position(
            pos=10, v_seq_start=5, v_seq_end=100,
            v_boundaries=v_boundaries,
            junction_start=100, junction_end=120, j_seq_end=150,
        )
        assert result == 'FWR1'


# ---------------------------------------------------------------------------
# Integration: selection pressure via Experiment DSL
# ---------------------------------------------------------------------------

class TestSelectionIntegration:

    def test_selection_end_to_end(self):
        from GenAIRR.ops import rate, with_antigen_selection
        result = (Experiment.on("human_igh")
                  .mutate(rate(0.02, 0.08), with_antigen_selection(0.5))
                  .run(n=10, seed=42))
        df = result.to_dataframe()
        assert len(df) == 10
        assert 'mutation_rate' in df.columns

    def test_with_corruption(self):
        from GenAIRR.ops import rate, with_antigen_selection, with_5prime_loss, with_indels, with_ns
        result = (Experiment.on("human_igh")
                  .mutate(rate(0.02, 0.08), with_antigen_selection(0.5))
                  .sequence(with_5prime_loss())
                  .observe(with_indels(), with_ns())
                  .run(n=10, seed=42))
        df = result.to_dataframe()
        assert len(df) == 10

    def test_with_productive_mode(self):
        from GenAIRR.ops import rate, with_antigen_selection
        result = (Experiment.on("human_igh")
                  .mutate(rate(0.02, 0.08), with_antigen_selection(0.5))
                  .run(n=10, seed=42, productive=True))
        df = result.to_dataframe()
        assert len(df) == 10

    def test_light_chain(self):
        from GenAIRR.ops import rate, with_antigen_selection
        result = (Experiment.on("human_igk")
                  .mutate(rate(0.02, 0.08), with_antigen_selection(0.5))
                  .run(n=10, seed=42))
        assert len(result) == 10

    def test_tcr(self):
        from GenAIRR.ops import rate, with_antigen_selection
        result = (Experiment.on("human_tcrb")
                  .mutate(rate(0.02, 0.08), with_antigen_selection(0.5))
                  .run(n=10, seed=42))
        assert len(result) == 10

    def test_deterministic_with_seed(self):
        from GenAIRR.ops import rate, with_antigen_selection
        def _run():
            return (Experiment.on("human_igh")
                    .mutate(rate(0.05, 0.10), with_antigen_selection(0.7))
                    .run(n=20, seed=123))
        df1 = _run().to_dataframe()
        df2 = _run().to_dataframe()
        assert df1['mutation_rate'].tolist() == df2['mutation_rate'].tolist()
        assert df1['sequence'].tolist() == df2['sequence'].tolist()
