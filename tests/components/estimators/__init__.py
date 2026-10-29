"""Estimating model parameters from observed repertoires.

One directory per estimator. Inside each, test_pins.py pins the live surfaces
the estimator reuses, test_estimation.py exercises its own behavior, and
test_composition.py runs it against the full release stack -- the distinction
the _contract / _implementation / _release filenames used to carry.
"""
