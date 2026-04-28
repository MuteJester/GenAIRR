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
        """T2-2: duplicate rate() emits UserWarning and last-write-wins."""
        with pytest.warns(UserWarning, match=r"rate\(\) specified more than once"):
            steps = _lower(mutate=[rate(0.01, 0.05), rate(0.02, 0.08)])
        m = [s for s in steps if isinstance(s, Mutate)][0]
        assert m.min_rate == 0.02
        assert m.max_rate == 0.08

    def test_model_last_write_wins(self):
        """T2-2: duplicate model() emits UserWarning and last-write-wins."""
        with pytest.warns(UserWarning,
                          match=r"model\(\) specified more than once"):
            steps = _lower(mutate=[model("s5f"), model("uniform")])
        # Both accepted; last wins. The compiled step now carries the
        # last-written model name (T2-4: uniform is no longer silently
        # rewritten to S5F).
        m = [s for s in steps if isinstance(s, (Mutate, SimulateCSR))]
        assert len(m) == 1
        assert m[0].model == "uniform"

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
        """T2-2: duplicate with_antigen_selection() emits UserWarning."""
        with pytest.warns(UserWarning,
                          match=r"with_antigen_selection\(\) specified more than once"):
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

    def test_uniform_model_no_warning(self):
        """T2-4: uniform is now wired through to the C kernel — no
        RuntimeWarning fires (the silent S5F-fallback warning has been
        deleted; the kernel actually runs)."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            steps = _lower(mutate=[model("uniform"), rate(0.01, 0.05)])
        runtime_warnings = [x for x in w if issubclass(x.category, RuntimeWarning)]
        assert len(runtime_warnings) == 0
        # Step carries the model name end-to-end.
        m = [s for s in steps if isinstance(s, Mutate)][0]
        assert m.model == "uniform"

    def test_s5f_model_no_warning(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            steps = _lower(mutate=[model("s5f"), rate(0.01, 0.05)])
        runtime_warnings = [x for x in w if issubclass(x.category, RuntimeWarning)]
        assert len(runtime_warnings) == 0
        m = [s for s in steps if isinstance(s, Mutate)][0]
        assert m.model == "s5f"


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
        """T2-2: duplicate with_umi() emits UserWarning."""
        with pytest.warns(UserWarning,
                          match=r"with_umi\(\) specified more than once"):
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
        """T2-2: duplicate with_contaminants() emits UserWarning."""
        with pytest.warns(UserWarning,
                          match=r"with_contaminants\(\) specified more than once"):
            steps = _lower(observe=[
                with_contaminants(0.01),
                with_contaminants(0.05, source="phix"),
            ])
        c = [s for s in steps if isinstance(s, SpikeContaminants)]
        assert len(c) == 1
        assert c[0].probability == 0.05
        assert c[0].contaminant_type == "phix"


# =====================================================================
# T2-2: duplicate-clause warning across all phases
# =====================================================================

class TestDuplicateClauseWarnings:
    """T2-2: every silent last-write-wins overwrite must surface as a
    ``UserWarning`` so accidental over-specification is visible. The
    behavior of the compiled steps is unchanged (last-write-wins)."""

    # ── Mutate phase clauses without a dedicated last-write-wins test ──

    def test_isotype_rates_duplicate_warns(self):
        """Two with_isotype_rates() are merged silently otherwise.
        Not strictly a value-overwrite (clause is a flag), but a
        duplicate is still a likely typo."""
        with pytest.warns(UserWarning,
                          match=r"with_isotype_rates\(\) specified more than once"):
            _lower(mutate=[
                rate(0.01, 0.05),
                with_isotype_rates(),
                with_isotype_rates(),
            ])

    # ── Prepare phase clauses ──────────────────────────────────────────

    def test_primer_mask_duplicate_warns(self):
        with pytest.warns(UserWarning,
                          match=r"with_primer_mask\(\) specified more than once"):
            _lower(prepare=[with_primer_mask(), with_primer_mask()])

    def test_pcr_duplicate_warns(self):
        with pytest.warns(UserWarning,
                          match=r"with_pcr\(\) specified more than once"):
            _lower(prepare=[with_pcr(), with_pcr()])

    # ── Sequence phase clauses ─────────────────────────────────────────

    def test_paired_end_duplicate_warns(self):
        with pytest.warns(UserWarning,
                          match=r"paired_end\(\) specified more than once"):
            _lower(sequence=[paired_end(250), paired_end(300)])

    def test_long_read_duplicate_warns(self):
        with pytest.warns(UserWarning,
                          match=r"long_read\(\) specified more than once"):
            _lower(sequence=[long_read(), long_read(error_rate=0.05)])

    def test_5prime_loss_duplicate_warns(self):
        with pytest.warns(UserWarning,
                          match=r"with_5prime_loss\(\) specified more than once"):
            _lower(sequence=[with_5prime_loss(), with_5prime_loss()])

    def test_3prime_loss_duplicate_warns(self):
        with pytest.warns(UserWarning,
                          match=r"with_3prime_loss\(\) specified more than once"):
            _lower(sequence=[with_3prime_loss(), with_3prime_loss()])

    def test_quality_profile_duplicate_warns(self):
        with pytest.warns(UserWarning,
                          match=r"with_quality_profile\(\) specified more than once"):
            _lower(sequence=[
                with_quality_profile(base=0.001, peak=0.01),
                with_quality_profile(base=0.005, peak=0.05),
            ])

    def test_reverse_complement_duplicate_warns(self):
        with pytest.warns(UserWarning,
                          match=r"with_reverse_complement\(\) specified more than once"):
            _lower(sequence=[
                with_reverse_complement(0.3),
                with_reverse_complement(0.5),
            ])

    # ── Observe phase clauses ──────────────────────────────────────────

    def test_indels_duplicate_warns(self):
        with pytest.warns(UserWarning,
                          match=r"with_indels\(\) specified more than once"):
            _lower(observe=[with_indels(0.005), with_indels(0.01)])

    def test_ns_duplicate_warns(self):
        with pytest.warns(UserWarning,
                          match=r"with_ns\(\) specified more than once"):
            _lower(observe=[with_ns(0.005), with_ns(0.01)])

    # ── Negative cases: NO warning when clauses are not duplicates ─────

    def test_no_warning_for_distinct_clauses_in_same_phase(self):
        """A phase with one of each clause type should produce no
        duplicate-clause warnings."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _lower(
                mutate=[rate(0.01, 0.05), model("s5f"), with_antigen_selection()],
                prepare=[with_primer_mask(), with_umi(12), with_pcr()],
                sequence=[
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
        dup_warnings = [
            x for x in w
            if issubclass(x.category, UserWarning)
            and "specified more than once" in str(x.message)
        ]
        assert len(dup_warnings) == 0, \
            f"Expected no duplicate-clause warnings, got: " \
            f"{[str(x.message) for x in dup_warnings]}"

    def test_warning_count_matches_duplicate_count(self):
        """Three rate() clauses → exactly two warnings (one per
        overwrite). One rate() = no warning, two rate() = one warning."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _lower(mutate=[
                rate(0.01, 0.05),
                rate(0.02, 0.06),
                rate(0.03, 0.07),
            ])
        dup_warnings = [
            x for x in w
            if issubclass(x.category, UserWarning)
            and "rate() specified more than once" in str(x.message)
        ]
        assert len(dup_warnings) == 2

    def test_each_phase_warns_independently(self):
        """A duplicate in one phase should not suppress the warning
        in another phase."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _lower(
                mutate=[rate(0.01, 0.05), rate(0.02, 0.06)],
                prepare=[with_umi(8), with_umi(16)],
            )
        dup_warnings = [
            x for x in w
            if issubclass(x.category, UserWarning)
            and "specified more than once" in str(x.message)
        ]
        # One for rate, one for with_umi.
        assert len(dup_warnings) == 2
        messages = " ".join(str(x.message) for x in dup_warnings)
        assert "rate()" in messages
        assert "with_umi()" in messages


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


# =====================================================================
# T2-3: using() allele typo validation with did-you-mean suggestions
# =====================================================================

class TestUsingAlleleValidation:
    """T2-3: typos in ``using()`` clauses must surface at compile time
    with a clear ``ValueError`` and ``difflib`` suggestions, not as an
    opaque C-side ``Unknown segment 'v' or allele 'X' not found in pool``
    error from deep in the call stack."""

    @pytest.fixture
    def imgt_config(self):
        # IMGT-style names (e.g. IGHV1-2*01) — covers the common case.
        from GenAIRR.data import HUMAN_IGH_IMGT
        return HUMAN_IGH_IMGT

    @pytest.fixture
    def ogrdb_config(self):
        # OGRDB-style names (e.g. IGHVF1-G2*02) — different convention.
        from GenAIRR.data import HUMAN_IGH_OGRDB
        return HUMAN_IGH_OGRDB

    # ── Happy path: no error for valid alleles ──────────────────────

    def test_valid_v_allele_compiles(self, imgt_config):
        """A correctly-spelled allele compiles without raising."""
        steps = _lower(recombine=[using(v="IGHV1-2*01")], config=imgt_config)
        assert any(isinstance(s, LockAlleles) for s in steps)

    def test_valid_alleles_all_segments_compile(self, imgt_config):
        steps = _lower(
            recombine=[using(
                v="IGHV1-2*01",
                d="IGHD1-1*01",
                j="IGHJ4*02",
            )],
            config=imgt_config,
        )
        assert any(isinstance(s, LockAlleles) for s in steps)

    # ── Bad allele names raise ValueError with suggestions ──────────

    def test_unknown_v_allele_raises_with_suggestion(self, imgt_config):
        """Typo in V allele suffix → ValueError with did-you-mean."""
        with pytest.raises(ValueError) as excinfo:
            _lower(
                recombine=[using(v="IGHV1-2*99")],
                config=imgt_config,
            )
        msg = str(excinfo.value)
        assert "Unknown V allele" in msg
        assert "'IGHV1-2*99'" in msg
        assert "Did you mean" in msg
        # The suggestion should include the same gene's real alleles.
        assert "IGHV1-2" in msg

    def test_unknown_j_allele_raises_with_suggestion(self, imgt_config):
        with pytest.raises(ValueError) as excinfo:
            _lower(
                recombine=[using(j="IGHJ4*99")],
                config=imgt_config,
            )
        msg = str(excinfo.value)
        assert "Unknown J allele" in msg
        assert "'IGHJ4*99'" in msg
        assert "IGHJ4" in msg  # close-match suggestion exists

    def test_typo_with_missing_dash_corrected(self, imgt_config):
        """Common typo: ``IGHV12*01`` (missing dash) → suggests
        ``IGHV1-2*01`` (the right answer)."""
        with pytest.raises(ValueError) as excinfo:
            _lower(
                recombine=[using(v="IGHV12*01")],
                config=imgt_config,
            )
        assert "IGHV1-2*01" in str(excinfo.value)

    def test_completely_garbage_name_no_suggestions(self, imgt_config):
        """Very dissimilar input → 'No close matches' message."""
        with pytest.raises(ValueError) as excinfo:
            _lower(
                recombine=[using(v="totally_made_up_xyz")],
                config=imgt_config,
            )
        msg = str(excinfo.value)
        assert "Unknown V allele" in msg
        assert "No close matches" in msg

    # ── Multiple typos collected into one error ─────────────────────

    def test_multiple_typos_all_reported(self, imgt_config):
        """Two bad alleles in one using() → one ValueError listing both,
        not two separate exceptions."""
        with pytest.raises(ValueError) as excinfo:
            _lower(
                recombine=[using(v="IGHV1-2*99", j="IGHJ4*99")],
                config=imgt_config,
            )
        msg = str(excinfo.value)
        assert "Unknown V allele" in msg
        assert "Unknown J allele" in msg
        assert "IGHV1-2*99" in msg
        assert "IGHJ4*99" in msg

    def test_multiple_typos_in_same_segment_all_reported(self, imgt_config):
        """List of V alleles where some are bad → all bad ones reported."""
        with pytest.raises(ValueError) as excinfo:
            _lower(
                recombine=[using(v=["IGHV1-2*01", "IGHV1-2*99", "IGHV9-9*99"])],
                config=imgt_config,
            )
        msg = str(excinfo.value)
        assert "IGHV1-2*99" in msg
        assert "IGHV9-9*99" in msg
        # The valid one should not appear as an error.
        assert "Unknown V allele 'IGHV1-2*01'" not in msg

    # ── Cross-config typo (IMGT name on OGRDB config) ───────────────

    def test_imgt_name_on_ogrdb_config_rejected(self, ogrdb_config):
        """User passes an IMGT-style ``IGHV1-2*02`` to an OGRDB config
        whose names use ``IGHVF1-G2*02``. Validator must catch."""
        with pytest.raises(ValueError) as excinfo:
            _lower(
                recombine=[using(v="IGHV1-2*02")],
                config=ogrdb_config,
            )
        assert "Unknown V allele" in str(excinfo.value)

    # ── Skip path: no config → no validation ────────────────────────

    def test_no_config_no_validation(self):
        """When ``resolved_config`` is None (some test paths do this),
        the validator must be a no-op so unrelated tests don't break."""
        # Should NOT raise.
        steps = _lower(recombine=[using(v="absolutely_not_real_xyz")])
        assert any(isinstance(s, LockAlleles) for s in steps)

    # ── Compile-time, not run-time ──────────────────────────────────

    def test_error_fires_at_compile_time_not_runtime(self):
        """Verifies the audit's central UX claim: typos surface during
        ``Experiment.compile()``, not when a sequence is generated.

        Uses the public Experiment API end-to-end so the test catches
        regressions where validation is moved out of the compile path.
        """
        from GenAIRR import Experiment
        with pytest.raises(ValueError) as excinfo:
            (Experiment.on("human_igh_imgt")
             .recombine(using(v="IGHV1-2*99"))
             .compile(seed=42))
        assert "Did you mean" in str(excinfo.value)
