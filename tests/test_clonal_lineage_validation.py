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


def test_validate_records_true_no_corruption_does_not_raise():
    result = _lineage_exp().run_records(seed=0, validate_records=True)
    assert len(result.records) > 0


def test_validate_records_true_with_corruption_does_not_raise():
    exp = _lineage_exp().sequencing_errors(rate=0.02)
    result = exp.run_records(seed=0, validate_records=True)
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
    result = _lineage_exp().run_records(seed=0, expose_provenance=True)
    assert len(result.records) > 0
    for rec in result.records:
        assert "truth_v_call" in rec
        assert rec["truth_v_call"]


def test_expose_provenance_adds_truth_v_call_with_corruption():
    exp = _lineage_exp().sequencing_errors(rate=0.02)
    result = exp.run_records(seed=0, expose_provenance=True)
    assert len(result.records) > 0
    for rec in result.records:
        assert "truth_v_call" in rec
        assert rec["truth_v_call"]


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
