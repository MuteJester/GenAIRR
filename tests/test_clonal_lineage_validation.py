"""Tests for clonal_lineage DSL validation (Fix 1–4)."""
import pytest
import GenAIRR as ga


def _base_exp(**kw):
    defaults = dict(n_clones=2, max_generations=4, n_max=100, n_sample=10,
                    rate=0.03, lambda_base=1.5)
    defaults.update(kw)
    return ga.Experiment.on("human_igh").recombine().clonal_lineage(**defaults)


# ---------------------------------------------------------------------------
# Fix 1: lambda_mut is removed — passing it must be a TypeError
# ---------------------------------------------------------------------------

def test_lambda_mut_removed_raises_typeerror():
    """lambda_mut is no longer a parameter; passing it must raise TypeError."""
    with pytest.raises(TypeError):
        _base_exp(lambda_mut=0.1)


# ---------------------------------------------------------------------------
# Fix 3: target_aa validation
# ---------------------------------------------------------------------------

def test_target_aa_empty_string_raises():
    with pytest.raises(ValueError, match="non-empty amino-acid"):
        _base_exp(target_aa="")


def test_target_aa_invalid_chars_raises():
    with pytest.raises(ValueError, match="invalid characters"):
        _base_exp(target_aa="not-an-aa-$$")


def test_target_aa_invalid_chars_mixed_raises():
    with pytest.raises(ValueError, match="invalid characters"):
        _base_exp(target_aa="ACDEFG1HIKLM")


def test_target_aa_non_aa_chars_raises():
    """Characters outside the 20 standard AAs (X, B, Z, U, O, numbers) must be rejected."""
    with pytest.raises(ValueError):
        _base_exp(target_aa="XBZUOACDEFGHIKL")
    with pytest.raises(ValueError):
        _base_exp(target_aa="ACDEFG0HIKLMN")


def test_target_aa_none_is_valid():
    """None (auto-target) must be accepted without error."""
    exp = _base_exp(target_aa=None)
    assert exp is not None


def test_target_aa_valid_long_string():
    """A realistic receptor-length AA string must be accepted."""
    aa = "ACDEFGHIKLMNPQRSTVWY" * 20   # 400 aa — long enough for no warning
    exp = _base_exp(target_aa=aa)
    assert exp is not None


def test_target_aa_short_triggers_warning(recwarn):
    """A very short target_aa (< 30 aa) must emit a UserWarning."""
    aa = "ACDEFGHIKLMN"  # 12 aa — below the 30-aa soft threshold
    import warnings
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _base_exp(target_aa=aa)
    user_warnings = [w for w in caught if issubclass(w.category, UserWarning)]
    assert len(user_warnings) == 1
    assert "receptor" in str(user_warnings[0].message).lower()


# ---------------------------------------------------------------------------
# Fix 3: s5f_model validation at call time
# ---------------------------------------------------------------------------

def test_s5f_model_bogus_raises_at_call_time():
    """A typo in s5f_model must error at clonal_lineage() call time."""
    with pytest.raises(ValueError, match="Unknown s5f_model"):
        _base_exp(s5f_model="bogus")


def test_s5f_model_valid_accepted():
    for model in ("hh_s5f", "hkl_s5f", "hh_s5f_60", "hh_s5f_opposite"):
        exp = _base_exp(s5f_model=model)
        assert exp is not None


# ---------------------------------------------------------------------------
# Fix 4: run_records kwargs on the lineage path
# ---------------------------------------------------------------------------

def _compiled_lineage():
    return _base_exp()


def test_n_not_none_raises_for_lineage():
    """Record count is not a fixed product for lineage; passing n is rejected."""
    with pytest.raises(ValueError):
        _compiled_lineage().run_records(n=5)


def test_seed_and_strict_still_work_for_lineage():
    """seed and strict must work."""
    result = _compiled_lineage().run_records(seed=42, strict=False)
    assert len(result.records) > 0


# ---------------------------------------------------------------------------
# validate_records / expose_provenance now SUPPORTED on the lineage path
# ---------------------------------------------------------------------------

def _lineage_exp(**kw):
    """A clonal_lineage experiment that reliably yields observed records."""
    defaults = dict(
        n_clones=2, max_generations=6, n_max=300, n_sample=15,
        rate=0.05, lambda_base=1.5, selection_strength=0.0, target_aa=None,
    )
    defaults.update(kw)
    return ga.Experiment.on("human_igh").recombine().clonal_lineage(**defaults)


# seed chosen so the (single-founder) families survive: sampling now draws from
# the living final generation, so an all-extinct seed yields zero records.
_SURVIVING_SEED = 1


def test_validate_records_true_no_corruption_does_not_raise():
    result = _lineage_exp().run_records(seed=_SURVIVING_SEED, validate_records=True)
    assert len(result.records) > 0


def test_validate_records_true_with_corruption_does_not_raise():
    exp = _lineage_exp().sequencing_errors(rate=0.02)
    result = exp.run_records(seed=_SURVIVING_SEED, validate_records=True)
    assert len(result.records) > 0


def test_outcomes_are_retained_and_index_aligned():
    result = _lineage_exp().run_records(seed=0)
    assert result.outcomes is not None
    assert len(result.outcomes) == len(result.records)


def test_outcomes_retained_with_corruption():
    result = _lineage_exp().sequencing_errors(rate=0.02).run_records(seed=0)
    assert result.outcomes is not None
    assert len(result.outcomes) == len(result.records)


def test_validate_records_called_directly_is_clean_no_corruption():
    exp = _lineage_exp()
    refdata = exp.compile().refdata
    result = exp.run_records(seed=0)
    report = result.validate_records(refdata)
    assert report.ok, report.summary()


def test_validate_records_called_directly_is_clean_with_corruption():
    exp = _lineage_exp().sequencing_errors(rate=0.02)
    refdata = exp.compile().refdata
    result = exp.run_records(seed=0)
    report = result.validate_records(refdata)
    assert report.ok, report.summary()


def test_expose_provenance_adds_truth_v_call_no_corruption():
    result = _lineage_exp().run_records(seed=_SURVIVING_SEED, expose_provenance=True)
    assert len(result.records) > 0
    for rec in result.records:
        assert "truth_v_call" in rec
        assert rec["truth_v_call"]


def test_expose_provenance_adds_truth_v_call_with_corruption():
    exp = _lineage_exp().sequencing_errors(rate=0.02)
    result = exp.run_records(seed=_SURVIVING_SEED, expose_provenance=True)
    assert len(result.records) > 0
    for rec in result.records:
        assert "truth_v_call" in rec
        assert rec["truth_v_call"]


# ---------------------------------------------------------------------------
# TCR guard: clonal_lineage models B-cell SHM => BCR/Ig only
# ---------------------------------------------------------------------------

def test_clonal_lineage_rejects_tcr_locus():
    """clonal_lineage applies S5F SHM (a B-cell process) so it must reject
    TCR loci with a clear BCR-only message — mirroring mutate()'s TCR guard.
    The guard must fire at clonal_lineage() call time (before compile) so the
    error is the BCR-only message, not a cartridge/compile error.
    """
    with pytest.raises(ValueError) as excinfo:
        (
            ga.Experiment.on("human_tcrb")
            .allow_curatable_refdata()
            .recombine()
            .clonal_lineage(n_clones=2, n_sample=10)
        )
    msg = str(excinfo.value)
    assert "BCR" in msg
    assert "TCR" in msg


def test_clonal_lineage_allows_bcr_locus():
    """Sanity: the BCR-only guard does not fire on an Ig locus."""
    exp = ga.Experiment.on("human_igh").recombine().clonal_lineage(
        n_clones=1, n_sample=5
    )
    assert exp is not None


# ---------------------------------------------------------------------------
# Founder-survival guard: every requested clone yields a surviving family
# ---------------------------------------------------------------------------

def test_survival_guard_yields_all_clones_by_default():
    """By default (allow_extinction=False) every requested clone survives via
    deterministic retry, so clone_ids == {0..n_clones-1} — even at a
    lambda_base where a single founder goes extinct ~25% of the time.
    """
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .clonal_lineage(n_clones=5, max_generations=6, n_max=200,
                        n_sample=15, rate=0.05, lambda_base=1.5)
        .run_records(seed=0)
    )
    cids = {r["clone_id"] for r in result.records}
    assert cids == {0, 1, 2, 3, 4}


def test_survival_guard_is_deterministic():
    """Same top-level seed => same result, even with retries."""
    def run():
        return (
            ga.Experiment.on("human_igh")
            .recombine()
            .clonal_lineage(n_clones=5, max_generations=6, n_max=200,
                            n_sample=15, rate=0.05, lambda_base=1.5)
            .run_records(seed=0)
        )
    a = run().records
    b = run().records
    assert len(a) == len(b)
    assert [r["sequence"] for r in a] == [r["sequence"] for r in b]
    assert [r["clone_id"] for r in a] == [r["clone_id"] for r in b]


def test_allow_extinction_true_allows_missing_clones():
    """allow_extinction=True accepts extinction and skips extinct clones, so
    clone_ids is a (possibly proper) subset of the requested range.
    """
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .clonal_lineage(n_clones=5, max_generations=6, n_max=200,
                        n_sample=15, rate=0.05, lambda_base=1.5,
                        allow_extinction=True)
        .run_records(seed=0)
    )
    cids = {r["clone_id"] for r in result.records}
    assert cids <= {0, 1, 2, 3, 4}


# ---------------------------------------------------------------------------
# AIRR-standard duplicate_count mirrors lineage_abundance
# ---------------------------------------------------------------------------

def test_records_emit_duplicate_count_equal_to_abundance():
    result = _base_exp(n_clones=3, lambda_base=1.6).run_records(seed=0)
    assert len(result.records) > 0
    for rec in result.records:
        assert "duplicate_count" in rec
        assert rec["duplicate_count"] == rec["lineage_abundance"]


# ---------------------------------------------------------------------------
# Smoke test: a valid call still works end-to-end
# ---------------------------------------------------------------------------

def test_valid_clonal_lineage_call_works():
    result = _base_exp(
        n_clones=2,
        selection_strength=0.0,
        target_aa=None,
    ).run_records(seed=7)
    assert len(result.records) > 0
    for rec in result.records:
        assert "lineage_affinity" in rec
        assert "clone_id" in rec
