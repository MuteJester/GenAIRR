"""Shared machinery for the test suite.

Test modules build on this package instead of importing from each other.

    fixtures/   pytest fixtures, including session-scoped heavy builds
    builders/   construct inputs; importable without pytest
    oracles/    derive expected truth independently of the code under test
    asserts/    domain assertions with legible failure messages
    cases/      shared parametrization tables
    strategies/ Hypothesis strategies
    paths       resolve fixture-data locations
"""
