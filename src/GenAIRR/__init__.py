from importlib.metadata import version as _pkg_version, PackageNotFoundError as _PackageNotFoundError

# Expose package version without hardcoding (single source of truth in pyproject)
try:
	__version__ = _pkg_version("GenAIRR")
except _PackageNotFoundError:
	# Fallback for editable installs during development when metadata is missing
	__version__ = "0.0.0"

# Keep imports minimal to avoid heavy side-effects at top-level import
# Users can import submodules directly, e.g., `from GenAIRR.sequence import HeavyChainSequence`
# from . import utilities  # optional convenience re-exports
# from . import alleles
# from . import sequence
# from . import mutation
# from . import simulation
