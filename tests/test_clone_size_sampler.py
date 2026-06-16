from GenAIRR import _engine


def test_sample_clone_sizes_basic():
    sizes = _engine.sample_clone_sizes(500, 0)           # defaults: power_law exp 2
    assert len(sizes) == 500
    assert all(s >= 1 for s in sizes)
    assert max(sizes) > 1                                 # heavy tail reaches >1


def test_sample_clone_sizes_deterministic():
    a = _engine.sample_clone_sizes(200, 7, kind="power_law", exponent=2.0)
    b = _engine.sample_clone_sizes(200, 7, kind="power_law", exponent=2.0)
    assert a == b


def test_unexpanded_fraction_forces_singletons():
    sizes = _engine.sample_clone_sizes(1000, 3, unexpanded_fraction=0.4)
    assert sum(1 for s in sizes if s == 1) >= 400


def test_lognormal_and_validation():
    import pytest
    sizes = _engine.sample_clone_sizes(100, 1, kind="lognormal", mu=1.0, sigma=1.0, max_size=10000)
    assert len(sizes) == 100 and all(1 <= s <= 10000 for s in sizes)
    with pytest.raises(ValueError):
        _engine.sample_clone_sizes(10, 0, kind="bogus")
    with pytest.raises(ValueError):
        _engine.sample_clone_sizes(10, 0, unexpanded_fraction=2.0)
