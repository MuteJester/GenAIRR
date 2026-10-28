"""Published germline sources: lock an explicit OGRDB release and build a portable reference from it.

``OgrdbSource`` acquires one explicit germline-set release (``fetch``, the only operation that
uses the network), records local copies of its exports (``from_files``) or reads an existing
lock with a verified cache (``read``).  Its lock pins the exact bytes.  ``inspect()`` reports
what the source holds, and ``build(recipe)`` projects an explicit recipe into a
``ReferenceBuild``: a runtime ``Reference`` plus the complete build report.

A build is written as a portable dossier with ``write_dossier``.  ``ReferenceBuild.read_dossier``
checks one fully offline: it regenerates every derived file from the locked bytes and the
recipe, and compares the result with what is stored.

Failures raise ``ReferenceSourceError`` with a stable ``kind``, a ``path`` locating the problem
and a message naming the affected asset, record or allele.
"""
from GenAIRR._reference_sources.api import OgrdbSource
from GenAIRR._reference_sources.build import ReferenceBuild
from GenAIRR._reference_sources.errors import ReferenceSourceError

__all__ = ("OgrdbSource", "ReferenceBuild", "ReferenceSourceError")
