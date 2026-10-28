#!/usr/bin/env python3
"""Generate the packaged archive of a reference build shipped with GenAIRR.

Maintainer tool (run from a source checkout; not part of the shipped wheel).  For one exact
build identifier it records the pinned provider exports under
``tests/data/reference_sources/ogrdb/`` as local files, together with the provider's data-terms
statement pages (the licence and the requested citation, recorded in the lock), and builds the
reviewed recipe next to them with ``OgrdbSource.build``.  The dossier is written into a private
temporary directory and checked there with ``ReferenceBuild.read_dossier``; its members are then
written as one deterministic zip archive (entries sorted by name, stored uncompressed, fixed
times and modes, nothing else), which is extracted again with the packaged-dossier loader and
compared byte for byte.  It never uses the network.  Acquisition receipts are machine-local
records and are not written.

The dossier records the GenAIRR version that runs this script as its builder.  Regenerate the
packaged archive whenever that version changes (after a release version bump) and pin the
printed manifest hash in ``src/GenAIRR/_reference_sources/bundled.py``.

Usage::

    # Write the starter archive into a scratch location:
    python tools/build_published_reference.py --output ./human_igh_ogrdb_10_observed_v1.zip

    # Regenerate the packaged starter (remove the old archive first):
    python tools/build_published_reference.py \\
        --output src/GenAIRR/data/published_references/human_igh_ogrdb_10_observed_v1.zip
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import re
import sys
import tempfile
from typing import Mapping, NamedTuple
import zipfile

import GenAIRR
from GenAIRR._reference_sources import bundled, ogrdb
from GenAIRR._reference_sources.acquisition import SuppliedFile, ingest_files
from GenAIRR._reference_sources.digests import manifest_hash
from GenAIRR._reference_sources.lock import Terms
from GenAIRR.reference_sources import OgrdbSource, ReferenceBuild, ReferenceSourceError

REPO = Path(__file__).resolve().parent.parent
PYPROJECT = REPO / "pyproject.toml"
SOURCE_DATA = REPO / "tests" / "data" / "reference_sources"
EXPORT_FILES = {"airr": "airr.json", "airr_ex": "airr_ex.json", "gapped": "gapped.fasta",
                "gapped_ex": "gapped_ex.fasta"}
# OGRDB's data terms: each statement page kept as a fixture, and the section of the fixture
# manifest (sources.json) that records the URL it was retrieved from and what it states.
OGRDB_STATEMENTS = {"licence_statement": ("ogrdb/statements/licensing_statement.html", "data_licence"),
                    "citation_statement": ("ogrdb/statements/help.html", "citation")}

MANIFEST = "manifest.json"
ARCHIVE_SUFFIX = bundled.ARCHIVE_SUFFIX
# Every archive entry is the same apart from its name and bytes: stored, at the earliest zip
# time, made on POSIX as a regular 0644 file, with no extra field and no comment.
ZIP_TIME = (1980, 1, 1, 0, 0, 0)
ZIP_FILE_ATTR = 0o100644 << 16


def _sources() -> dict:
    return json.loads((SOURCE_DATA / "sources.json").read_bytes().decode("utf-8"))


class Published(NamedTuple):
    """A packaged build: the provider set release, its interpretation profile and its recipe."""

    set_id: str
    release: int
    profile: str
    recipe: str

    @property
    def exports(self) -> dict:
        directory = SOURCE_DATA / "ogrdb" / self.set_id / str(self.release)
        return {fmt: directory / name for fmt, name in EXPORT_FILES.items()}

    @property
    def terms(self) -> Terms:
        """The provider's data terms, as ``sources.json`` records them: the licence, and the citation
        the provider requests (verbatim) followed by its DOI, so the citation stands on its own."""
        sources = _sources()
        citation = f"{sources['citation']['requested_by_provider']} doi:{sources['citation']['doi']}"
        return Terms(status="recorded", assets=tuple(OGRDB_STATEMENTS), licence=sources["data_licence"]["licence"],
                     citation=citation, reason=None)

    @property
    def recipe_path(self) -> Path:
        return SOURCE_DATA / self.recipe


BUILDS = {
    "human_igh_ogrdb_10_observed_v1": Published(
        set_id="9606.IGH_VDJ", release=10, profile="ogrdb_api_v2_human_ig_observed_v1",
        recipe="human_igh_ogrdb_10_observed_v1.recipe.json"),
}
DEFAULT_ID = "human_igh_ogrdb_10_observed_v1"


class GeneratorError(Exception):
    """A target the generator will not write into."""


def source(exact_id: str, *, cache) -> OgrdbSource:
    """The pinned exports of ``exact_id`` and the provider's statement pages, recorded as local files
    with the terms they state; ``cache`` is an empty scratch directory."""
    published = BUILDS[exact_id]
    files = {fmt: SuppliedFile(path) for fmt, path in published.exports.items()}
    sources = _sources()
    files.update({fmt: SuppliedFile(SOURCE_DATA / name, declared_url=sources[section]["statement_url"])
                  for fmt, (name, section) in OGRDB_STATEMENTS.items()})
    acquired = ingest_files(set_id=published.set_id, release=published.release, profile=published.profile,
                            files=files, cache=cache, identify=ogrdb.identify, terms=published.terms)
    return OgrdbSource._of(acquired)


def build(exact_id: str, *, cache) -> ReferenceBuild:
    """Build ``exact_id`` from its pinned exports, statements and recipe; ``cache`` is an empty scratch
    directory."""
    return source(exact_id, cache=cache).build(BUILDS[exact_id].recipe_path)


def declared_version(text: str) -> str:
    """``[project].version`` of a pyproject.toml text, read without tomllib (so on every Python)
    and with either line ending.  A missing table or version, or a dynamic version, is a
    ``GeneratorError`` saying so.  The tests read the declared version with this same function."""
    text = text.replace("\r\n", "\n")
    table = re.search(r"^\[project\][ \t]*(?:#.*)?$(.*?)(?=^\[|\Z)", text, re.M | re.S)
    if table is None:
        raise GeneratorError("pyproject.toml has no [project] table, so it declares no version to compare with")
    project = table.group(1)
    dynamic = re.search(r"^dynamic\s*=\s*\[([^\]]*)\]", project, re.M)
    if dynamic and re.search(r"[\"']version[\"']", dynamic.group(1)):
        raise GeneratorError("pyproject.toml declares the version as dynamic; the packaged dossier's builder can only "
                             "be compared with a version written in [project]")
    version = re.search(r"""^version\s*=\s*(["'])([^"']+)\1\s*(?:#.*)?$""", project, re.M)
    if version is None:
        raise GeneratorError("pyproject.toml has no version in its [project] table")
    return version.group(2)


def check_version() -> None:
    """The installed GenAIRR is the version this checkout declares: the dossier records it."""
    if not PYPROJECT.is_file():
        raise GeneratorError(f"{PYPROJECT} does not exist: run this script from a source checkout whose "
                             "pyproject.toml declares the version being packaged")
    declared = declared_version(PYPROJECT.read_bytes().decode("utf-8"))
    if GenAIRR.__version__ != declared:
        raise GeneratorError(f"the installed GenAIRR is {GenAIRR.__version__}, but pyproject.toml declares {declared}; "
                             "reinstall this checkout (for example with maturin develop --release) so that the "
                             "packaged dossier records the version being released")


def check_target(target: Path) -> None:
    """``target`` must be a new ``.zip`` file in an existing directory."""
    if not target.name.endswith(ARCHIVE_SUFFIX):
        raise GeneratorError(f"{target} must name a .zip file")
    if target.is_symlink() or target.exists():
        raise GeneratorError(f"{target} already exists; the archive is written as a new file")
    if not target.parent.is_dir():
        raise GeneratorError(f"the parent directory of {target} does not exist")


def archive(members: Mapping[str, bytes]) -> bytes:
    """The deterministic zip of ``members`` (name -> bytes): entries sorted by name, stored, with
    every entry field fixed, so equal members give equal bytes on every platform."""
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_STORED, allowZip64=False) as out:
        for name in sorted(members):
            entry = zipfile.ZipInfo(name, date_time=ZIP_TIME)
            entry.compress_type = zipfile.ZIP_STORED
            entry.create_system = 3
            entry.create_version = entry.extract_version = 20
            entry.external_attr = ZIP_FILE_ATTR
            entry.internal_attr = 0
            entry.flag_bits = 0
            entry.extra = b""
            entry.comment = b""
            out.writestr(entry, members[name])
        out.comment = b""
    return buffer.getvalue()


def write_archive(members: Mapping[str, bytes], target: Path) -> bytes:
    """Write the archive of ``members`` as the new file ``target`` and read it back as a package
    would; an archive that the loader refuses, or that does not give back exactly ``members``, is
    removed again."""
    check_target(target)
    data = archive(members)
    with open(target, "xb") as handle:
        handle.write(data)
    try:
        with tempfile.TemporaryDirectory(prefix="genairr_archive_") as scratch:
            unpacked = Path(scratch) / "dossier"
            bundled._extract(target, unpacked, str(target), None)
            found = {p.relative_to(unpacked).as_posix(): p.read_bytes() for p in unpacked.rglob("*") if p.is_file()}
        if found != dict(members):
            raise GeneratorError(f"the archive {target} does not extract to the dossier members it was written from")
    except BaseException:
        target.unlink()
        raise
    return data


def publish(result: ReferenceBuild, target: Path) -> ReferenceBuild:
    """Write ``result`` as the archive ``target`` (a new ``.zip`` file) and return the dossier checked
    before archiving.  Receipts are not written."""
    check_target(target)
    names = [member.name for member in result.manifest.members] + [MANIFEST]
    with tempfile.TemporaryDirectory(prefix="genairr_dossier_") as scratch:
        staged = Path(scratch) / "dossier"
        result.write_dossier(staged)
        checked = ReferenceBuild.read_dossier(staged)
        if checked.manifest != result.manifest:
            raise GeneratorError("the written dossier does not check as the build that was written")
        members = {name: staged.joinpath(*name.split("/")).read_bytes() for name in names}
    write_archive(members, target)
    return checked


def summary(exact_id: str, result: ReferenceBuild, data: bytes) -> str:
    """The identifiers to pin: the typed manifest hash, the runtime reference hashes and the archive digest."""
    reference = result.reference
    return "\n".join([
        f"build {exact_id}",
        f"builder GenAIRR {result.manifest.builder.version}",
        f"manifest {manifest_hash(result.manifest)}",
        f"runtime execution {reference.execution_hash}",
        f"runtime catalogue {reference.catalogue_hash}",
        f"runtime models {reference.models_hash}",
        f"runtime artifact {reference.artifact_hash}",
        f"archive sha256:{hashlib.sha256(data).hexdigest()} ({len(data)} bytes)",
        f"lock {result.lock.commitment}",
    ])


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--id", default=DEFAULT_ID, choices=sorted(BUILDS), help="exact build identifier")
    parser.add_argument("--output", required=True, type=Path, help="new .zip file for the dossier archive")
    args = parser.parse_args(argv)
    try:
        check_version()
        check_target(args.output)
        with tempfile.TemporaryDirectory(prefix="genairr_sources_") as scratch:
            result = build(args.id, cache=Path(scratch) / "cache")
        checked = publish(result, args.output)
    except (GeneratorError, ReferenceSourceError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(summary(args.id, checked, args.output.read_bytes()))
    print(f"checked with GenAIRR {GenAIRR.__version__}. Pin what changed: the manifest hash in the registry "
          "(src/GenAIRR/_reference_sources/bundled.py) and in MANIFEST_PIN, the archive digest in ARCHIVE_PIN, the "
          "runtime hashes in RUNTIME_PINS and the lock in LOCK_PIN (tests/test_published_reference_resources.py)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
