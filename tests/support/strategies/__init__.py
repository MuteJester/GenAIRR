"""Hypothesis strategies for the suite.

Oracle verification is naturally property-based: the claim is about any
simulated sequence, not three hand-picked ones. Strategies live here so a
property test states what it is asserting rather than how to build inputs.
"""
