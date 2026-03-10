"""
Tests for GenAIRR.ops — clause construction, validation, immutability, labels.
"""

import dataclasses
import pytest

from GenAIRR.ops import (
    # Base types
    Clause, RecombineClause, MutateClause, PrepareClause,
    SequenceClause, ObserveClause,
    # Concrete clause types
    UsingClause, DInversionClause, ReceptorRevisionClause,
    ModelClause, RateClause, IsotypeRatesClause, AntigenSelectionClause,
    UMIClause, PCRClause, PrimerMaskClause,
    PairedEndClause, LongReadClause, FivePrimeLossClause, ThreePrimeLossClause,
    QualityProfileClause, ReverseComplementClause,
    ContaminantsClause, IndelsClause, NsClause,
    # Constructor functions
    using, with_d_inversion, with_receptor_revision,
    model, rate, with_antigen_selection, with_isotype_rates,
    with_umi, with_pcr, with_primer_mask,
    paired_end, long_read, with_5prime_loss, with_3prime_loss,
    with_quality_profile, with_reverse_complement,
    with_contaminants, with_indels, with_ns,
)


# =====================================================================
# using()
# =====================================================================

class TestUsing:

    def test_default_all_none(self):
        c = using()
        assert c.v is None
        assert c.d is None
        assert c.j is None
        assert c.c is None

    def test_single_string_normalized_to_tuple(self):
        c = using(v="IGHV1-2*01")
        assert c.v == ("IGHV1-2*01",)

    def test_list_normalized_to_tuple(self):
        c = using(v=["IGHV1-2*01", "IGHV1-2*02"])
        assert c.v == ("IGHV1-2*01", "IGHV1-2*02")

    def test_keyword_only(self):
        with pytest.raises(TypeError):
            using("IGHV1-2*01")  # positional arg not allowed

    def test_multiple_segments(self):
        c = using(v="IGHV1-2*01", j="IGHJ4*02", d=["IGHD3-10*01"])
        assert c.v == ("IGHV1-2*01",)
        assert c.j == ("IGHJ4*02",)
        assert c.d == ("IGHD3-10*01",)
        assert c.c is None

    def test_is_recombine_clause(self):
        assert isinstance(using(), RecombineClause)


# =====================================================================
# model()
# =====================================================================

class TestModel:

    def test_default_s5f(self):
        c = model()
        assert c.name == "s5f"

    def test_uniform_accepted(self):
        c = model("uniform")
        assert c.name == "uniform"

    def test_invalid_raises_valueerror(self):
        with pytest.raises(ValueError, match="Unknown mutation model"):
            model("xyz")

    def test_is_mutate_clause(self):
        assert isinstance(model(), MutateClause)


# =====================================================================
# rate()
# =====================================================================

class TestRate:

    def test_defaults(self):
        c = rate()
        assert c.min_rate == 0.01
        assert c.max_rate == 0.05

    def test_custom(self):
        c = rate(0.02, 0.08)
        assert c.min_rate == 0.02
        assert c.max_rate == 0.08

    def test_min_greater_than_max_raises(self):
        with pytest.raises(ValueError, match="min_rate.*<=.*max_rate"):
            rate(0.08, 0.02)

    def test_negative_min_raises(self):
        with pytest.raises(ValueError, match="min_rate.*non-negative"):
            rate(-0.01, 0.05)

    def test_negative_max_raises(self):
        with pytest.raises(ValueError, match="max_rate.*non-negative"):
            rate(0.01, -0.05)

    def test_zero_rate_valid(self):
        c = rate(0.0, 0.0)
        assert c.min_rate == 0.0
        assert c.max_rate == 0.0

    def test_is_mutate_clause(self):
        assert isinstance(rate(), MutateClause)


# =====================================================================
# with_d_inversion()
# =====================================================================

class TestDInversion:

    def test_default(self):
        c = with_d_inversion()
        assert c.prob == 0.15

    def test_custom(self):
        c = with_d_inversion(0.3)
        assert c.prob == 0.3

    def test_negative_prob_raises(self):
        with pytest.raises(ValueError, match="prob.*between 0 and 1"):
            with_d_inversion(-0.1)

    def test_prob_over_1_raises(self):
        with pytest.raises(ValueError, match="prob.*between 0 and 1"):
            with_d_inversion(1.5)

    def test_zero_valid(self):
        c = with_d_inversion(0.0)
        assert c.prob == 0.0

    def test_one_valid(self):
        c = with_d_inversion(1.0)
        assert c.prob == 1.0

    def test_is_recombine_clause(self):
        assert isinstance(with_d_inversion(), RecombineClause)


# =====================================================================
# with_receptor_revision()
# =====================================================================

class TestReceptorRevision:

    def test_default(self):
        c = with_receptor_revision()
        assert c.prob == 0.05
        assert c.footprint_min == 5
        assert c.footprint_max == 20

    def test_custom(self):
        c = with_receptor_revision(prob=0.1, footprint=(3, 15))
        assert c.prob == 0.1
        assert c.footprint_min == 3
        assert c.footprint_max == 15

    def test_inverted_footprint_raises(self):
        with pytest.raises(ValueError, match="footprint max.*>=.*footprint min"):
            with_receptor_revision(footprint=(20, 5))

    def test_negative_footprint_raises(self):
        with pytest.raises(ValueError, match="footprint min.*non-negative"):
            with_receptor_revision(footprint=(-1, 5))

    def test_negative_prob_raises(self):
        with pytest.raises(ValueError, match="prob.*between 0 and 1"):
            with_receptor_revision(prob=-0.1)

    def test_is_recombine_clause(self):
        assert isinstance(with_receptor_revision(), RecombineClause)


# =====================================================================
# with_antigen_selection()
# =====================================================================

class TestAntigenSelection:

    def test_defaults(self):
        c = with_antigen_selection()
        assert c.strength == 0.5
        assert c.cdr_r_acceptance == 0.85
        assert c.fwr_r_acceptance == 0.40

    def test_custom(self):
        c = with_antigen_selection(strength=0.8, cdr_r_acceptance=0.9,
                                   fwr_r_acceptance=0.3)
        assert c.strength == 0.8
        assert c.cdr_r_acceptance == 0.9
        assert c.fwr_r_acceptance == 0.3

    def test_negative_strength_raises(self):
        with pytest.raises(ValueError, match="strength.*between 0 and 1"):
            with_antigen_selection(strength=-0.1)

    def test_is_mutate_clause(self):
        assert isinstance(with_antigen_selection(), MutateClause)


# =====================================================================
# with_isotype_rates()
# =====================================================================

class TestIsotypeRates:

    def test_returns_clause(self):
        c = with_isotype_rates()
        assert isinstance(c, IsotypeRatesClause)
        assert isinstance(c, MutateClause)


# =====================================================================
# with_umi()
# =====================================================================

class TestUMI:

    def test_default(self):
        c = with_umi()
        assert c.length == 12

    def test_custom(self):
        c = with_umi(16)
        assert c.length == 16

    def test_zero_raises(self):
        with pytest.raises(ValueError, match="length.*at least 1"):
            with_umi(0)

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="length.*at least 1"):
            with_umi(-1)

    def test_is_prepare_clause(self):
        assert isinstance(with_umi(), PrepareClause)


# =====================================================================
# with_pcr()
# =====================================================================

class TestPCR:

    def test_defaults(self):
        c = with_pcr()
        assert c.error_rate == 1e-4
        assert c.cycles == 30

    def test_custom(self):
        c = with_pcr(error_rate=5e-5, cycles=25)
        assert c.error_rate == 5e-5
        assert c.cycles == 25

    def test_negative_cycles_raises(self):
        with pytest.raises(ValueError, match="cycles.*at least 1"):
            with_pcr(cycles=0)

    def test_negative_error_rate_raises(self):
        with pytest.raises(ValueError, match="error_rate.*non-negative"):
            with_pcr(error_rate=-0.1)

    def test_is_prepare_clause(self):
        assert isinstance(with_pcr(), PrepareClause)


# =====================================================================
# with_primer_mask()
# =====================================================================

class TestPrimerMask:

    def test_default(self):
        c = with_primer_mask()
        assert c.length == 0  # 0 = full FR1

    def test_custom(self):
        c = with_primer_mask(50)
        assert c.length == 50

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="length.*non-negative"):
            with_primer_mask(-1)

    def test_is_prepare_clause(self):
        assert isinstance(with_primer_mask(), PrepareClause)


# =====================================================================
# paired_end()
# =====================================================================

class TestPairedEnd:

    def test_default(self):
        c = paired_end()
        assert c.read_length == 300

    def test_custom(self):
        c = paired_end(150)
        assert c.read_length == 150

    def test_zero_raises(self):
        with pytest.raises(ValueError, match="read_length.*at least 1"):
            paired_end(0)

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="read_length.*at least 1"):
            paired_end(-100)

    def test_is_sequence_clause(self):
        assert isinstance(paired_end(), SequenceClause)


# =====================================================================
# long_read()
# =====================================================================

class TestLongRead:

    def test_defaults(self):
        c = long_read()
        assert c.error_rate == 0.03
        assert c.min_run_length == 3
        assert c.insertion_bias == 0.6

    def test_custom(self):
        c = long_read(error_rate=0.05, min_run_length=4, insertion_bias=0.7)
        assert c.error_rate == 0.05
        assert c.min_run_length == 4
        assert c.insertion_bias == 0.7

    def test_negative_error_rate_raises(self):
        with pytest.raises(ValueError, match="error_rate.*non-negative"):
            long_read(error_rate=-0.01)

    def test_invalid_insertion_bias_raises(self):
        with pytest.raises(ValueError, match="insertion_bias.*between 0 and 1"):
            long_read(insertion_bias=1.5)

    def test_is_sequence_clause(self):
        assert isinstance(long_read(), SequenceClause)


# =====================================================================
# with_5prime_loss() / with_3prime_loss()
# =====================================================================

class TestPrimeLoss:

    def test_5prime_returns_clause(self):
        c = with_5prime_loss()
        assert isinstance(c, FivePrimeLossClause)
        assert isinstance(c, SequenceClause)

    def test_3prime_returns_clause(self):
        c = with_3prime_loss()
        assert isinstance(c, ThreePrimeLossClause)
        assert isinstance(c, SequenceClause)


# =====================================================================
# with_quality_profile()
# =====================================================================

class TestQualityProfile:

    def test_defaults(self):
        c = with_quality_profile()
        assert c.base == 0.001
        assert c.peak == 0.02

    def test_custom(self):
        c = with_quality_profile(base=0.005, peak=0.05)
        assert c.base == 0.005
        assert c.peak == 0.05

    def test_negative_base_raises(self):
        with pytest.raises(ValueError, match="base.*non-negative"):
            with_quality_profile(base=-0.001)

    def test_negative_peak_raises(self):
        with pytest.raises(ValueError, match="peak.*non-negative"):
            with_quality_profile(peak=-0.02)

    def test_is_sequence_clause(self):
        assert isinstance(with_quality_profile(), SequenceClause)


# =====================================================================
# with_reverse_complement()
# =====================================================================

class TestReverseComplement:

    def test_default(self):
        c = with_reverse_complement()
        assert c.prob == 0.5

    def test_custom(self):
        c = with_reverse_complement(0.3)
        assert c.prob == 0.3

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="prob.*between 0 and 1"):
            with_reverse_complement(-0.1)

    def test_over_1_raises(self):
        with pytest.raises(ValueError, match="prob.*between 0 and 1"):
            with_reverse_complement(1.5)

    def test_is_sequence_clause(self):
        assert isinstance(with_reverse_complement(), SequenceClause)


# =====================================================================
# with_contaminants()
# =====================================================================

class TestContaminants:

    def test_default(self):
        c = with_contaminants()
        assert c.rate == 0.01
        assert c.source == "random"

    def test_phix(self):
        c = with_contaminants(source="phix")
        assert c.source == "phix"

    def test_invalid_source_raises(self):
        with pytest.raises(ValueError, match="Unknown contaminant source"):
            with_contaminants(source="xyz")

    def test_negative_rate_raises(self):
        with pytest.raises(ValueError, match="rate.*between 0 and 1"):
            with_contaminants(rate=-0.1)

    def test_rate_over_1_raises(self):
        with pytest.raises(ValueError, match="rate.*between 0 and 1"):
            with_contaminants(rate=1.5)

    def test_is_observe_clause(self):
        assert isinstance(with_contaminants(), ObserveClause)


# =====================================================================
# with_indels()
# =====================================================================

class TestIndels:

    def test_default(self):
        c = with_indels()
        assert c.prob == 0.01

    def test_custom(self):
        c = with_indels(0.005)
        assert c.prob == 0.005

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="prob.*between 0 and 1"):
            with_indels(-0.01)

    def test_over_1_raises(self):
        with pytest.raises(ValueError, match="prob.*between 0 and 1"):
            with_indels(1.5)

    def test_is_observe_clause(self):
        assert isinstance(with_indels(), ObserveClause)


# =====================================================================
# with_ns()
# =====================================================================

class TestNs:

    def test_default(self):
        c = with_ns()
        assert c.prob == 0.01

    def test_custom(self):
        c = with_ns(0.005)
        assert c.prob == 0.005

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="prob.*between 0 and 1"):
            with_ns(-0.01)

    def test_over_1_raises(self):
        with pytest.raises(ValueError, match="prob.*between 0 and 1"):
            with_ns(1.5)

    def test_is_observe_clause(self):
        assert isinstance(with_ns(), ObserveClause)


# =====================================================================
# Immutability — all clauses are frozen
# =====================================================================

class TestImmutability:

    @pytest.mark.parametrize("clause", [
        using(v="IGHV1-2*01"),
        with_d_inversion(),
        with_receptor_revision(),
        model(),
        rate(),
        with_isotype_rates(),
        with_antigen_selection(),
        with_umi(),
        with_pcr(),
        with_primer_mask(),
        paired_end(),
        long_read(),
        with_5prime_loss(),
        with_3prime_loss(),
        with_quality_profile(),
        with_reverse_complement(),
        with_contaminants(),
        with_indels(),
        with_ns(),
    ], ids=lambda c: type(c).__name__)
    def test_frozen(self, clause):
        with pytest.raises(dataclasses.FrozenInstanceError):
            # Normal assignment triggers frozen __setattr__
            clause._frozen_test = True


# =====================================================================
# Labels — every clause produces a non-empty label
# =====================================================================

class TestLabels:

    @pytest.mark.parametrize("clause", [
        using(v="IGHV1-2*01"),
        with_d_inversion(),
        with_receptor_revision(),
        model(),
        model("uniform"),
        rate(),
        with_isotype_rates(),
        with_antigen_selection(),
        with_umi(),
        with_pcr(),
        with_primer_mask(),
        with_primer_mask(50),
        paired_end(),
        long_read(),
        with_5prime_loss(),
        with_3prime_loss(),
        with_quality_profile(),
        with_reverse_complement(),
        with_contaminants(),
        with_contaminants(source="phix"),
        with_indels(),
        with_ns(),
    ], ids=lambda c: f"{type(c).__name__}_{id(c)}")
    def test_label_not_empty(self, clause):
        label = clause.__label__()
        assert isinstance(label, str)
        assert len(label) > 0

    def test_using_label_shows_segments(self):
        c = using(v="IGHV1-2*01", j="IGHJ4*02")
        label = c.__label__()
        assert "V=" in label
        assert "J=" in label

    def test_model_label_shows_name(self):
        assert "S5F" in model("s5f").__label__()
        assert "UNIFORM" in model("uniform").__label__()

    def test_rate_label_shows_range(self):
        label = rate(0.02, 0.08).__label__()
        assert "0.02" in label
        assert "0.08" in label

    def test_primer_mask_label_full_fr1(self):
        assert "FR1" in with_primer_mask().__label__()

    def test_primer_mask_label_custom(self):
        assert "50" in with_primer_mask(50).__label__()


# =====================================================================
# Type hierarchy
# =====================================================================

class TestTypeHierarchy:

    def test_all_inherit_from_clause(self):
        clauses = [
            using(), with_d_inversion(), with_receptor_revision(),
            model(), rate(), with_isotype_rates(), with_antigen_selection(),
            with_umi(), with_pcr(), with_primer_mask(),
            paired_end(), long_read(), with_5prime_loss(), with_3prime_loss(),
            with_quality_profile(), with_reverse_complement(),
            with_contaminants(), with_indels(), with_ns(),
        ]
        for c in clauses:
            assert isinstance(c, Clause), f"{type(c).__name__} is not a Clause"

    def test_recombine_clauses(self):
        assert isinstance(using(), RecombineClause)
        assert isinstance(with_d_inversion(), RecombineClause)
        assert isinstance(with_receptor_revision(), RecombineClause)

    def test_mutate_clauses(self):
        assert isinstance(model(), MutateClause)
        assert isinstance(rate(), MutateClause)
        assert isinstance(with_isotype_rates(), MutateClause)
        assert isinstance(with_antigen_selection(), MutateClause)

    def test_prepare_clauses(self):
        assert isinstance(with_umi(), PrepareClause)
        assert isinstance(with_pcr(), PrepareClause)
        assert isinstance(with_primer_mask(), PrepareClause)

    def test_sequence_clauses(self):
        assert isinstance(paired_end(), SequenceClause)
        assert isinstance(long_read(), SequenceClause)
        assert isinstance(with_5prime_loss(), SequenceClause)
        assert isinstance(with_3prime_loss(), SequenceClause)
        assert isinstance(with_quality_profile(), SequenceClause)
        assert isinstance(with_reverse_complement(), SequenceClause)

    def test_observe_clauses(self):
        assert isinstance(with_contaminants(), ObserveClause)
        assert isinstance(with_indels(), ObserveClause)
        assert isinstance(with_ns(), ObserveClause)

    def test_no_cross_phase(self):
        """Recombine clauses are NOT mutate clauses, etc."""
        assert not isinstance(using(), MutateClause)
        assert not isinstance(model(), RecombineClause)
        assert not isinstance(with_umi(), SequenceClause)
        assert not isinstance(paired_end(), ObserveClause)
        assert not isinstance(with_contaminants(), PrepareClause)
