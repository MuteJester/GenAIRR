"""Verification that GenAIRR's recorded metadata is true.

GenAIRR is used to evaluate other tools, and everything downstream treats its
ground-truth metadata as a perfect oracle: when it says a read's np1 is five
untemplated bases, every benchmark built on it assumes that is what actually
happened. So the obligation here is not "does the code run" but "is the
recorded truth actually true".

Which mechanism a test verifies is declared by the test, wherever it lives --
see catalog. Coverage therefore aggregates from components/, from here, and
from the Rust suites alike, instead of depending on a directory.
"""
