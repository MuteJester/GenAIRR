"""
Tests for clause-to-step lowering in GenAIRR.protocol._clauses_to_steps.
"""

import warnings
import pytest

from GenAIRR.ops import (
    using, with_d_inversion, with_receptor_revision,
    model, rate, with_antigen_selection, with_isotype_rates,
    with_umi, with_pcr, with_primer_mask,
    paired_end, long_read, with_5prime_loss, with_3prime_loss,
    with_quality_profile, with_reverse_complement,
    with_contaminants, with_indels, with_ns,
)
from GenAIRR.protocol import _clauses_to_steps
from GenAIRR.steps import (
    Rearrange, LockAlleles, SimulateDGeneInversion,
    SimulateReceptorRevision, Mutate, SimulateCSR, SelectionPressure,
    PrimerMask, SimulateUMI, PCRAmplification,
    Corrupt5Prime, Corrupt3Prime, CorruptQuality, SimulatePairedEnd,
    SkewBaseComposition, ReverseComplement,
    SpikeContaminants, InsertIndels, InsertNs,
)


def _lower(*,
           recombine=None, mutate=None, prepare=None,
           sequence=None, observe=None, config=None):
    """Helper: call _clauses_to_steps with keyword lists."""
    return _clauses_to_steps(
        recombine_clauses=recombine or [],
        mutate_clauses=mutate or [],
        prepare_clauses=prepare or [],
        sequence_clauses=sequence or [],
        observe_clauses=observe or [],
        resolved_config=config,
    )


# =====================================================================
# Empty experiment
# =====================================================================

class TestEmptyExperiment:

    def test_empty_produces_rearrange_only(self):
        steps = _lower()
        assert len(steps) == 1
        assert isinstance(steps[0], Rearrange)


# =====================================================================
# Recombine phase
# =====================================================================

class TestRecombineLowering:

    def test_using_produces_lock_alleles(self):
        steps = _lower(recombine=[using(v="IGHV1-2*01")])
        lock_steps = [s for s in steps if isinstance(s, LockAlleles)]
        assert len(lock_steps) == 1
        assert lock_steps[0].locks == {"v": ["IGHV1-2*01"]}

    def test_using_multiple_segments(self):
        steps = _lower(recombine=[using(v="IGHV1-2*01", j="IGHJ4*02")])
        lock = [s for s in steps if isinstance(s, LockAlleles)][0]
        assert lock.locks["v"] == ["IGHV1-2*01"]
        assert lock.locks["j"] == ["IGHJ4*02"]

    def test_using_segment_merge_across_calls(self):
        """Two using() clauses with different segments merge."""
        steps = _lower(recombine=[
            using(v="IGHV1-2*01"),
            using(j="IGHJ4*02"),
        ])
        lock = [s for s in steps if isinstance(s, LockAlleles)][0]
        assert lock.locks["v"] == ["IGHV1-2*01"]
        assert lock.locks["j"] == ["IGHJ4*02"]

    def test_using_same_segment_last_write_wins(self):
        """Two using() clauses for the same segment: last wins."""
        steps = _lower(recombine=[
            using(v="IGHV1-2*01"),
            using(v="IGHV3-30*01"),
        ])
        lock = [s for s in steps if isinstance(s, LockAlleles)][0]
        assert lock.locks["v"] == ["IGHV3-30*01"]

    def test_d_inversion_produces_step(self):
        steps = _lower(recombine=[with_d_inversion(0.2)])
        dinv = [s for s in steps if isinstance(s, SimulateDGeneInversion)]
        assert len(dinv) == 1
        assert dinv[0].probability == 0.2

    def test_d_inversion_last_write_wins(self):
        steps = _lower(recombine=[
            with_d_inversion(0.1),
            with_d_inversion(0.3),
        ])
        dinv = [s for s in steps if isinstance(s, SimulateDGeneInversion)]
        assert len(dinv) == 1
        assert dinv[0].probability == 0.3

    def test_receptor_revision_produces_step(self):
        steps = _lower(recombine=[
            with_receptor_revision(prob=0.1, footprint=(3, 15))
        ])
        rev = [s for s in steps if isinstance(s, SimulateReceptorRevision)]
        assert len(rev) == 1
        assert rev[0].probability == 0.1
        assert rev[0].footprint_min == 3
        assert rev[0].footprint_max == 15

    def test_all_recombine_clauses_together(self):
        steps = _lower(recombine=[
            using(v="IGHV1-2*01"),
            with_d_inversion(0.15),
            with_receptor_revision(),
        ])
        assert any(isinstance(s, LockAlleles) for s in steps)
        assert any(isinstance(s, SimulateDGeneInversion) for s in steps)
        assert any(isinstance(s, SimulateReceptorRevision) for s in steps)

    def test_empty_using_produces_no_lock(self):
        """using() with all None segments produces no LockAlleles step."""
        steps = _lower(recombine=[using()])
        assert not any(isinstance(s, LockAlleles) for s in steps)


# =====================================================================
# Mutate phase
# =====================================================================

class TestMutateLowering:

    def test_rate_only_produces_mutate(self):
        steps = _lower(mutate=[rate(0.02, 0.08)])
        m = [s for s in steps if isinstance(s, Mutate)]
        assert len(m) == 1
        assert m[0].min_rate == 0.02
        assert m[0].max_rate == 0.08

    def test_defaults_when_no_rate(self):
        """If only model() specified, default rate is used."""
        steps = _lower(mutate=[model("s5f")])
        m = [s for s in steps if isinstance(s, Mutate)]
        assert len(m) == 1
        assert m[0].min_rate == 0.01
        assert m[0].max_rate == 0.05

    def test_rate_last_write_wins(self):
        steps = _lower(mutate=[rate(0.01, 0.05), rate(0.02, 0.08)])
        m = [s for s in steps if isinstance(s, Mutate)][0]
        assert m.min_rate == 0.02
        assert m.max_rate == 0.08

    def test_model_last_write_wins(self):
        """model() is last-write-wins — no error for calling twice."""
        # Both accepted; last wins (both produce Mutate, just label differs)
        steps = _lower(mutate=[model("s5f"), model("uniform")])
        # Should still produce one Mutate step
        m = [s for s in steps if isinstance(s, (Mutate, SimulateCSR))]
        assert len(m) == 1

    def test_isotype_rates_produces_csr(self):
        steps = _lower(mutate=[rate(0.01, 0.05), with_isotype_rates()])
        csr = [s for s in steps if isinstance(s, SimulateCSR)]
        assert len(csr) == 1
        assert csr[0].min_rate == 0.01
        assert csr[0].max_rate == 0.05
        # Should NOT also have a Mutate
        assert not any(isinstance(s, Mutate) for s in steps)

    def test_selection_produces_step(self):
        steps = _lower(mutate=[rate(0.02, 0.08), with_antigen_selection(0.7)])
        sel = [s for s in steps if isinstance(s, SelectionPressure)]
        assert len(sel) == 1
        assert sel[0].strength == 0.7

    def test_selection_last_write_wins(self):
        steps = _lower(mutate=[
            with_antigen_selection(0.3),
            with_antigen_selection(0.8),
        ])
        sel = [s for s in steps if isinstance(s, SelectionPressure)]
        assert len(sel) == 1
        assert sel[0].strength == 0.8

    def test_full_mutate_phase(self):
        steps = _lower(mutate=[
            model("s5f"),
            rate(0.01, 0.05),
            with_isotype_rates(),
            with_antigen_selection(0.5),
        ])
        assert any(isinstance(s, SimulateCSR) for s in steps)
        assert any(isinstance(s, SelectionPressure) for s in steps)
        assert not any(isinstance(s, Mutate) for s in steps)

    def test_uniform_model_warns(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _lower(mutate=[model("uniform"), rate(0.01, 0.05)])
        runtime_warnings = [x for x in w if issubclass(x.category, RuntimeWarning)]
        assert len(runtime_warnings) == 1
        assert "uniform" in str(runtime_warnings[0].message).lower()

    def test_s5f_model_no_warning(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _lower(mutate=[model("s5f"), rate(0.01, 0.05)])
        runtime_warnings = [x for x in w if issubclass(x.category, RuntimeWarning)]
        assert len(runtime_warnings) == 0


# =====================================================================
# Prepare phase
# =====================================================================

class TestPrepareLowering:

    def test_umi_produces_step(self):
        steps = _lower(prepare=[with_umi(16)])
        umi = [s for s in steps if isinstance(s, SimulateUMI)]
        assert len(umi) == 1
        assert umi[0].umi_length == 16

    def test_pcr_produces_step(self):
        steps = _lower(prepare=[with_pcr(error_rate=5e-5, cycles=25)])
        pcr = [s for s in steps if isinstance(s, PCRAmplification)]
        assert len(pcr) == 1
        assert pcr[0].error_rate == 5e-5
        assert pcr[0].n_cycles == 25

    def test_primer_mask_produces_step(self):
        steps = _lower(prepare=[with_primer_mask(50)])
        pm = [s for s in steps if isinstance(s, PrimerMask)]
        assert len(pm) == 1
        assert pm[0].mask_length == 50

    def test_all_prepare_clauses(self):
        steps = _lower(prepare=[
            with_primer_mask(),
            with_umi(12),
            with_pcr(),
        ])
        assert any(isinstance(s, PrimerMask) for s in steps)
        assert any(isinstance(s, SimulateUMI) for s in steps)
        assert any(isinstance(s, PCRAmplification) for s in steps)

    def test_last_write_wins_for_umi(self):
        steps = _lower(prepare=[with_umi(8), with_umi(16)])
        umi = [s for s in steps if isinstance(s, SimulateUMI)]
        assert len(umi) == 1
        assert umi[0].umi_length == 16


# =====================================================================
# Sequence phase
# =====================================================================

class TestSequenceLowering:

    def test_paired_end_produces_step(self):
        steps = _lower(sequence=[paired_end(250)])
        pe = [s for s in steps if isinstance(s, SimulatePairedEnd)]
        assert len(pe) == 1
        assert pe[0].read_length == 250

    def test_long_read_produces_step(self):
        steps = _lower(sequence=[long_read(error_rate=0.05)])
        lr = [s for s in steps if isinstance(s, SkewBaseComposition)]
        assert len(lr) == 1
        assert lr[0].error_rate == 0.05

    def test_5prime_loss_produces_step(self):
        steps = _lower(sequence=[with_5prime_loss()])
        assert any(isinstance(s, Corrupt5Prime) for s in steps)

    def test_3prime_loss_produces_step(self):
        steps = _lower(sequence=[with_3prime_loss()])
        assert any(isinstance(s, Corrupt3Prime) for s in steps)

    def test_quality_profile_produces_step(self):
        steps = _lower(sequence=[with_quality_profile(base=0.005, peak=0.05)])
        q = [s for s in steps if isinstance(s, CorruptQuality)]
        assert len(q) == 1
        assert q[0].base_error_rate == 0.005
        assert q[0].peak_error_rate == 0.05

    def test_reverse_complement_produces_step(self):
        steps = _lower(sequence=[with_reverse_complement(0.3)])
        rc = [s for s in steps if isinstance(s, ReverseComplement)]
        assert len(rc) == 1
        assert rc[0].probability == 0.3

    def test_paired_end_plus_long_read_warns(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _lower(sequence=[paired_end(300), long_read()])
        user_warnings = [x for x in w if issubclass(x.category, UserWarning)]
        assert len(user_warnings) == 1
        assert "conflict" in str(user_warnings[0].message).lower()

    def test_full_sequence_phase(self):
        steps = _lower(sequence=[
            paired_end(300),
            with_5prime_loss(),
            with_3prime_loss(),
            with_quality_profile(),
            with_reverse_complement(),
        ])
        assert any(isinstance(s, SimulatePairedEnd) for s in steps)
        assert any(isinstance(s, Corrupt5Prime) for s in steps)
        assert any(isinstance(s, Corrupt3Prime) for s in steps)
        assert any(isinstance(s, CorruptQuality) for s in steps)
        assert any(isinstance(s, ReverseComplement) for s in steps)


# =====================================================================
# Observe phase
# =====================================================================

class TestObserveLowering:

    def test_contaminants_produces_step(self):
        steps = _lower(observe=[with_contaminants(0.05, source="phix")])
        c = [s for s in steps if isinstance(s, SpikeContaminants)]
        assert len(c) == 1
        assert c[0].probability == 0.05
        assert c[0].contaminant_type == "phix"

    def test_indels_produces_step(self):
        steps = _lower(observe=[with_indels(0.005)])
        i = [s for s in steps if isinstance(s, InsertIndels)]
        assert len(i) == 1
        assert i[0].probability == 0.005

    def test_ns_produces_step(self):
        steps = _lower(observe=[with_ns(0.005)])
        n = [s for s in steps if isinstance(s, InsertNs)]
        assert len(n) == 1
        assert n[0].probability == 0.005

    def test_all_observe_clauses(self):
        steps = _lower(observe=[
            with_contaminants(0.01),
            with_indels(0.005),
            with_ns(0.005),
        ])
        assert any(isinstance(s, SpikeContaminants) for s in steps)
        assert any(isinstance(s, InsertIndels) for s in steps)
        assert any(isinstance(s, InsertNs) for s in steps)

    def test_last_write_wins_for_contaminants(self):
        steps = _lower(observe=[
            with_contaminants(0.01),
            with_contaminants(0.05, source="phix"),
        ])
        c = [s for s in steps if isinstance(s, SpikeContaminants)]
        assert len(c) == 1
        assert c[0].probability == 0.05
        assert c[0].contaminant_type == "phix"


# =====================================================================
# Full pipeline
# =====================================================================

class TestFullPipeline:

    def test_all_phases_produce_correct_step_count(self):
        steps = _lower(
            recombine=[
                using(v="IGHV1-2*01"),
                with_d_inversion(0.15),
                with_receptor_revision(),
            ],
            mutate=[
                rate(0.01, 0.05),
                with_isotype_rates(),
                with_antigen_selection(0.5),
            ],
            prepare=[
                with_primer_mask(),
                with_umi(12),
                with_pcr(),
            ],
            sequence=[
                paired_end(300),
                with_5prime_loss(),
                with_3prime_loss(),
                with_quality_profile(),
                with_reverse_complement(),
            ],
            observe=[
                with_contaminants(0.01),
                with_indels(0.005),
                with_ns(0.005),
            ],
        )
        # Rearrange(1) + Lock(1) + DInv(1) + RecRev(1) +
        # CSR(1) + Selection(1) +
        # PrimerMask(1) + UMI(1) + PCR(1) +
        # 5'(1) + 3'(1) + Quality(1) + PairedEnd(1) + RC(1) +
        # Contaminants(1) + Indels(1) + Ns(1) = 17
        assert len(steps) == 17
        assert isinstance(steps[0], Rearrange)

    def test_step_ordering(self):
        """Steps are emitted in biological order regardless of clause order."""
        steps = _lower(
            observe=[with_ns(0.01)],
            mutate=[rate(0.02, 0.08)],
            recombine=[with_d_inversion(0.15)],
            sequence=[with_5prime_loss()],
        )
        types = [type(s).__name__ for s in steps]
        # Should be: Rearrange, DInversion, Mutate, 5'Corrupt, Ns
        assert types == [
            "Rearrange",
            "SimulateDGeneInversion",
            "Mutate",
            "Corrupt5Prime",
            "InsertNs",
        ]


# =====================================================================
# Chain compatibility warnings
# =====================================================================

class TestChainCompatibilityWarnings:

    def _make_config_with_chain(self, chain_type):
        """Create a minimal mock config with a chain type."""
        class MockMetadata:
            def __init__(self, ct):
                self.chain_type = ct
        class MockConfig:
            def __init__(self, ct):
                self.metadata = MockMetadata(ct)
        return MockConfig(chain_type)

    def test_d_inversion_on_kappa_warns(self):
        from GenAIRR.dataconfig.enums import ChainType
        config = self._make_config_with_chain(ChainType.BCR_LIGHT_KAPPA)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _lower(
                recombine=[with_d_inversion(0.15)],
                config=config,
            )
        user_warnings = [x for x in w if issubclass(x.category, UserWarning)]
        assert len(user_warnings) == 1
        assert "no D segment" in str(user_warnings[0].message)

    def test_d_inversion_on_lambda_warns(self):
        from GenAIRR.dataconfig.enums import ChainType
        config = self._make_config_with_chain(ChainType.BCR_LIGHT_LAMBDA)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _lower(
                recombine=[with_d_inversion(0.15)],
                config=config,
            )
        user_warnings = [x for x in w if issubclass(x.category, UserWarning)]
        assert len(user_warnings) == 1

    def test_d_inversion_on_heavy_no_warning(self):
        from GenAIRR.dataconfig.enums import ChainType
        config = self._make_config_with_chain(ChainType.BCR_HEAVY)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _lower(
                recombine=[with_d_inversion(0.15)],
                config=config,
            )
        user_warnings = [x for x in w if issubclass(x.category, UserWarning)]
        assert len(user_warnings) == 0

    def test_d_inversion_on_tcr_beta_no_warning(self):
        from GenAIRR.dataconfig.enums import ChainType
        config = self._make_config_with_chain(ChainType.TCR_BETA)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _lower(
                recombine=[with_d_inversion(0.15)],
                config=config,
            )
        user_warnings = [x for x in w if issubclass(x.category, UserWarning)]
        assert len(user_warnings) == 0

    def test_d_inversion_on_tcr_alpha_warns(self):
        from GenAIRR.dataconfig.enums import ChainType
        config = self._make_config_with_chain(ChainType.TCR_ALPHA)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _lower(
                recombine=[with_d_inversion(0.15)],
                config=config,
            )
        user_warnings = [x for x in w if issubclass(x.category, UserWarning)]
        assert len(user_warnings) == 1

    def test_no_config_no_warning(self):
        """Without a resolved config, no chain compatibility warning."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _lower(recombine=[with_d_inversion(0.15)])
        user_warnings = [x for x in w if issubclass(x.category, UserWarning)]
        assert len(user_warnings) == 0
