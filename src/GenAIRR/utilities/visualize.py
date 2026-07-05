"""Standalone HTML "exploding view" of a simulated AIRR record.

Public entry point: :func:`visualize_sequence`. Implementation lives in
the `_visualize/` package (parse / styles / components / alignment / render).
"""
from ._visualize.render import visualize_sequence

__all__ = ["visualize_sequence"]
