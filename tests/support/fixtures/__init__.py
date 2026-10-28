"""Fixtures and fixture-building helpers shared across test modules.

Each module here is the single definitional home for setup that several test
modules need.  Before this package existed, a test module that needed another's
setup imported it directly, which coupled their collection order and made the
lender impossible to relocate.
"""
