"""The keyed seed derivation, as Python sees it.

A clonal run splits into sub-simulations — a clone's founder, a family's
retry attempt, a read, an observed node — and each needs its own seed.
Deriving those by adding an offset to the run seed makes the whole
schedule one lattice shared by every run seed, so a different seed
reuses streams instead of drawing new ones. The engine already has the
primitive that does not do that; this is the surface Python derives
through, so that there is one keyed derivation rather than two.
"""

import pytest

import GenAIRR._engine as _engine

pytestmark = pytest.mark.behavior


def test_a_role_and_an_index_name_one_seed():
    """Same arguments, same word; different role or index, different word."""
    assert _engine.substream_seed(7, "clone", 3) == _engine.substream_seed(7, "clone", 3)
    assert _engine.substream_seed(7, "clone", 3) != _engine.substream_seed(7, "read", 3)
    assert _engine.substream_seed(7, "clone", 3) != _engine.substream_seed(7, "clone", 4)
    assert _engine.substream_seed(7, "clone", 3) != _engine.substream_seed(8, "clone", 3)


def test_the_run_seed_is_mixed_in_rather_than_added_to():
    """The defect this exists to remove, stated as a property.

    Under ``seed + clone_idx * 1_000_000`` the seed of clone *c* at run
    seed *S* is the seed of clone 0 at run seed ``S + c * 1_000_000``, so
    two runs share streams. Keep the old collision in the test: a
    regression test for a scheme nobody can still write has to show what
    it is protecting against.
    """
    old = lambda seed, clone: seed + clone * 1_000_000
    assert old(0, 1) == old(1_000_000, 0)
    assert old(0, 2) == old(2_000_000, 0)

    new = lambda seed, clone: _engine.substream_seed(seed, "clone", clone)
    assert new(0, 1) != new(1_000_000, 0)
    assert new(0, 2) != new(2_000_000, 0)

    # Not just those two pairs: no seed drawn for one run is drawn for
    # another, over a grid wide enough that a translation could not hide.
    grid = {}
    for run_seed in (0, 1, 2, 1_000_000, 2_000_000, 2**63):
        for clone in range(64):
            word = new(run_seed, clone)
            assert word not in grid, f"{(run_seed, clone)} collides with {grid[word]}"
            grid[word] = (run_seed, clone)


def test_two_axes_under_one_clone_do_not_meet():
    """A read index and a retry index must not name the same stream.

    Under the old schedule a node took ``clone_seed + 1 + id`` and an
    attempt took ``clone_seed + attempt * 7919``, so node 7918 and
    attempt 1 drew the same words.
    """
    clone_seed = 4242
    assert clone_seed + 1 + 7918 == clone_seed + 1 * 7919

    seen = {}
    for role in ("attempt", "read", "node"):
        for index in list(range(256)) + [7918, 7919, 15837]:
            word = _engine.substream_seed(clone_seed, role, index)
            assert word not in seen, f"{(role, index)} collides with {seen[word]}"
            seen[word] = (role, index)


def test_an_index_the_label_cannot_hold_is_refused():
    """Truncating silently would alias two sub-simulations."""
    assert _engine.substream_seed(1, "read", 2**32 - 1) != _engine.substream_seed(1, "read", 0)
    with pytest.raises((OverflowError, ValueError)):
        _engine.substream_seed(1, "read", 2**32)
    with pytest.raises((OverflowError, ValueError)):
        _engine.substream_seed(2**64, "read", 0)
    with pytest.raises((OverflowError, ValueError)):
        _engine.substream_seed(-1, "read", 0)


def test_a_stream_with_no_role_is_refused():
    """An empty role would make every axis the same axis."""
    with pytest.raises(ValueError):
        _engine.substream_seed(1, "", 0)


def test_the_boundary_computes_what_the_crate_computes():
    """One pinned vector, asserted on both sides of the binding.

    The properties above would all hold for a wrapper that quietly
    derived something else of its own. ``engine_rs/src/python/rng.rs``
    asserts this same triple, so a wrapper that stops calling the
    crate's primitive fails here rather than drifting.
    """
    assert _engine.substream_seed(0, "clone", 0) == 0xE6DC_F4CF_2CC6_3646
