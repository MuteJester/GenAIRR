"""The closed registry of reference build dossiers packaged with GenAIRR.

Each entry maps an exact identifier to a dossier directory packaged as a resource of a Python
package, and to the manifest hash that dossier must have.  There are no aliases, no
pattern matching and no downloads: an identifier is either listed here or unknown.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
import importlib.resources
from pathlib import Path
import tempfile
from types import MappingProxyType
from typing import Iterator, Mapping

from GenAIRR._reference_sources._fields import bad, check_digest, check_member_path
from GenAIRR._reference_sources.manifest import MAX_MEMBERS

MAX_RESOURCE = 255
# A packaged dossier is copied out of a non-file package (for example a zip archive) within the
# same bounds the dossier reader enforces: entries and total bytes.
_MAX_ENTRIES = 2 * MAX_MEMBERS
_MAX_BYTES = 1 << 30


@dataclass(frozen=True)
class BundledDossier:
    """A packaged dossier: the directory ``resource`` ('/'-separated) inside ``package``, and the
    ``genairr_reference_source_build_manifest_v1`` hash its manifest must have."""

    package: str
    resource: str
    manifest_hash: str

    def __post_init__(self):
        if type(self.package) is not str or not all(part.isidentifier() for part in self.package.split(".")):
            raise bad("/package", f"expected an importable package name, got {self.package!r}")
        check_member_path(self.resource, "/resource", MAX_RESOURCE)
        check_digest(self.manifest_hash, "/manifest_hash")


# Exact identifier -> packaged dossier.
REGISTRY: Mapping[str, BundledDossier] = MappingProxyType({})


def entry(exact_id: str) -> BundledDossier:
    """The registry entry of exactly ``exact_id``; an unknown identifier lists the available ones."""
    if type(exact_id) is not str:
        raise TypeError(f"a bundled build identifier is a string, got {type(exact_id).__name__}")
    try:
        return REGISTRY[exact_id]
    except KeyError:
        available = ", ".join(sorted(REGISTRY)) or "none"
        raise bad("", f"no bundled reference build is named exactly {exact_id!r}; available: {available}",
                  "missing_asset") from None


def _copy(source, target: Path, budget: list) -> None:
    """Copy a packaged directory tree (any importlib resource) into ``target``, within bounds."""
    target.mkdir()
    for item in sorted(source.iterdir(), key=lambda x: x.name):
        budget[0] -= 1
        if budget[0] < 0:
            raise bad("", f"the packaged dossier holds more than {_MAX_ENTRIES} entries", "resource_limit")
        if item.is_dir():
            _copy(item, target / item.name, budget)
        else:
            data = item.read_bytes()
            budget[1] -= len(data)
            if budget[1] < 0:
                raise bad("", f"the packaged dossier holds more than {_MAX_BYTES} bytes", "resource_limit")
            (target / item.name).write_bytes(data)


@contextmanager
def dossier_directory(dossier: BundledDossier) -> Iterator[Path]:
    """A filesystem directory holding the packaged dossier, for reading only.

    Where the package's resources are plain files this is the packaged directory itself;
    otherwise the dossier is copied into a private temporary directory, removed afterwards.
    Nothing is ever written into the package.
    """
    try:
        root = importlib.resources.files(dossier.package)
    except (ImportError, TypeError) as exc:
        raise bad("", f"the package {dossier.package!r} of the bundled dossier cannot be found: {exc}",
                  "missing_asset") from exc
    resource = root.joinpath(*dossier.resource.split("/"))
    if isinstance(resource, Path):
        yield resource
        return
    if not resource.is_dir():
        raise bad("", f"the bundled dossier {dossier.resource} is not packaged in {dossier.package}", "missing_asset")
    with tempfile.TemporaryDirectory(prefix="genairr_bundled_") as temporary:
        target = Path(temporary) / "dossier"
        _copy(resource, target, [_MAX_ENTRIES, _MAX_BYTES])
        yield target
