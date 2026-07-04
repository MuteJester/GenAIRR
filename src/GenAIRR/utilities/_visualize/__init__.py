"""Internal implementation of the HTML "exploding view" renderer.

Split into cohesive submodules: parse / styles / components / alignment /
render. The public entry point is :func:`visualize_sequence`, re-exported
here and by :mod:`GenAIRR.utilities.visualize`.
"""

from .render import visualize_sequence

__all__ = ["visualize_sequence"]
