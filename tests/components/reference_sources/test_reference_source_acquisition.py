"""Bounded, explicit acquisition of published germline sources into the blob cache.

``fetch_source`` is the only operation that touches the network, and only through the
transport it is handed: these tests hand it a controlled fake, so nothing here needs a
live network.  Offline operations (reading a lock with a cache, verifying its blobs,
importing local files) take no transport at all and must not reach the network on any
path, including a cache miss or malformed input.  Receipts record each acquisition
event and never enter the lock: identical bytes always seal the same lock.
"""
from __future__ import annotations

import copy
import datetime
import errno
import hashlib
import http.client
import inspect
import io
import os
import pathlib
import socket
import ssl
import threading

import pytest

from GenAIRR import _artifact_io
from GenAIRR._artifact_io import canonical_json, strict_loads
from GenAIRR._reference_sources import acquisition as acq
from GenAIRR._reference_sources import model as m
from GenAIRR._reference_sources.acquisition import (
    ACQUISITION_DEADLINE_SECONDS, FETCH_FORMATS, MAX_FOLLOWED_REDIRECTS, OGRDB_ORIGIN, SOCKET_TIMEOUT_SECONDS,
    AcquiredSource, DeadlineExceeded, HttpsTransport, SuppliedFile, fetch_source, ingest_files, provider_url,
    read_locked, resolve_profile,
)
from GenAIRR._reference_sources.cache import BlobCache
from GenAIRR._reference_sources.errors import ReferenceSourceError
from tests.support.fixtures.source_acquisition import ALL, FakeResponse, FakeTransport, IGH, IGK, SOURCES, _fixture, _sources_asset, _url, ok, serve

pytestmark = pytest.mark.behavior

PROFILE = "ogrdb_api_v2_human_ig_observed_v1"
T0 = datetime.datetime(2026, 10, 27, 9, 30, tzinfo=datetime.timezone.utc)


def _sha(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def identify(*, set_id, release, assets):
    """Test-local identity reader: the set, release, taxon and locus of a record export."""
    for fmt in ("airr", "airr_ex"):
        if fmt in assets:
            (germline_set,) = strict_loads(assets[fmt].decode("utf-8"), max_bytes=1 << 25)["GermlineSet"]
            taxonomy_id = int(germline_set["species"]["id"].split(":")[1])
            served = f"{taxonomy_id}.{germline_set['germline_set_name'].removesuffix('_extended')}"
            if served != set_id or germline_set["release_version"] != release:
                raise ReferenceSourceError("identity_mismatch", "", f"served {served} release "
                                           f"{germline_set['release_version']}, requested {set_id} release {release}")
            locus = germline_set["locus"]
            return m.SourceIdentity(m.Taxon(germline_set["species"]["label"], taxonomy_id), locus,
                                    m.LOCUS_CHAIN_TYPES[locus])
    raise ReferenceSourceError("identity_mismatch", "", "no record export to identify the set from")


def _error(kind, call, *args, **kwargs) -> ReferenceSourceError:
    with pytest.raises(ReferenceSourceError) as info:
        call(*args, **kwargs)
    assert info.value.kind == kind, f"expected {kind!r}, got {info.value.kind!r}: {info.value}"
    return info.value


def _snapshot(directory: pathlib.Path) -> dict:
    if not directory.exists():
        return {}
    return {p.relative_to(directory).as_posix(): (p.read_bytes() if p.is_file() else None)
            for p in sorted(directory.rglob("*"))}


def _temporaries(directory: pathlib.Path) -> list:
    return sorted(p.name for p in directory.rglob("*.tmp")) if directory.exists() else []


# ── controlled transport ─────────────────────────────────────────────


def redirect(location, status=302):
    return lambda: FakeResponse(status, b"", headers={"Location": location} if location else {})


class FailingTransport:
    def get(self, url, *, timeout, deadline):
        pytest.fail(f"the transport was called for {url}")


def fetch(tmp_path, transport, *, source=IGH, formats=ALL, cache="cache", **kwargs):
    kwargs.setdefault("now", lambda: T0)
    return fetch_source(set_id=source[0], release=source[1], profile=PROFILE, formats=formats,
                        cache=tmp_path / cache, transport=transport, identify=identify, **kwargs)


@pytest.fixture
def no_network(monkeypatch):
    """Any attempt to open a connection fails the test."""
    def refuse(*args, **kwargs):
        pytest.fail("an offline operation tried to reach the network")

    monkeypatch.setattr(socket, "create_connection", refuse)
    monkeypatch.setattr(socket.socket, "connect", refuse)
    monkeypatch.setattr(http.client.HTTPSConnection, "connect", refuse)
    monkeypatch.setattr(http.client.HTTPConnection, "connect", refuse)
    monkeypatch.setattr(HttpsTransport, "get", refuse)


# ── constants and request validation ─────────────────────────────────


def test_limits_are_explicit_constants_shared_with_the_source_model():
    assert SOCKET_TIMEOUT_SECONDS == 30
    assert ACQUISITION_DEADLINE_SECONDS == 120
    assert MAX_FOLLOWED_REDIRECTS == 3
    assert acq.MAX_ASSET_BYTES == m.MAX_ASSET_BYTES == 16 * 1024 * 1024
    assert acq.MAX_SOURCE_PAYLOAD_BYTES == m.MAX_SOURCE_PAYLOAD_BYTES == 128 * 1024 * 1024
    assert FETCH_FORMATS == ("airr", "airr_ex", "gapped", "gapped_ex")
    assert OGRDB_ORIGIN == "https://ogrdb.airr-community.org"


def test_provider_urls_are_generated_from_validated_components():
    for entry in SOURCES["sets"]:
        for asset in entry["assets"]:
            if asset["format"] in FETCH_FORMATS:
                assert provider_url(entry["set_id"], entry["release"], asset["format"]) == asset["url"]


@pytest.mark.parametrize("set_id,release,fmt,kind", [
    ("9606.IGH_VDJ", "latest", "airr", "malformed"),
    ("9606.IGH_VDJ", True, "airr", "malformed"),
    ("9606.IGH_VDJ", 0, "airr", "malformed"),
    ("9606.IGH_VDJ", -3, "airr", "malformed"),
    ("9606.IGH_VDJ", 10.0, "airr", "malformed"),
    ("9606.IGH_VDJ", "10", "airr", "malformed"),
    ("9606.IGH_VDJ", 2 ** 31, "airr", "malformed"),
    ("../admin", 10, "airr", "malformed"),
    ("9606.IGH_VDJ/10/airr?x=", 10, "airr", "malformed"),
    ("9606 IGH", 10, "airr", "malformed"),
    ("%2e%2e", 10, "airr", "malformed"),
    ("", 10, "airr", "malformed"),
    (None, 10, "airr", "malformed"),
    ("9606.IGH_VDJ", 10, "ungapped", "unsupported_profile"),
    ("9606.IGH_VDJ", 10, "licence_statement", "unsupported_profile"),
    ("9606.IGH_VDJ", 10, "latest", "unsupported_profile"),
    ("9606.IGH_VDJ", 10, None, "malformed"),
])
def test_provider_urls_reject_unvalidated_components(set_id, release, fmt, kind):
    _error(kind, provider_url, set_id, release, fmt)


def test_profiles_resolve_to_a_supported_adapter():
    adapter = resolve_profile(PROFILE)
    assert adapter == m.AdapterRef("ogrdb_germline_set", 1, m.CoordinateProfileRef(PROFILE, 1))
    assert resolve_profile(adapter) is adapter
    _error("unsupported_profile", resolve_profile, "ogrdb_api_v2_human_ig_observed_v2")
    _error("unsupported_profile", resolve_profile, "latest")
    _error("malformed", resolve_profile, True)
    _error("malformed", resolve_profile, None)


BAD_REQUESTS = {
    "release-latest": (dict(release="latest"), "malformed", "/release"),
    "release-boolean": (dict(release=True), "malformed", "/release"),
    "release-zero": (dict(release=0), "malformed", "/release"),
    "release-negative": (dict(release=-1), "malformed", "/release"),
    "release-float": (dict(release=10.0), "malformed", "/release"),
    "set-path": (dict(set_id="../9606.IGH_VDJ"), "malformed", "/set_id"),
    "set-slash": (dict(set_id="9606.IGH_VDJ/10"), "malformed", "/set_id"),
    "profile-unknown": (dict(profile="ogrdb_api_v2_human_ig_observed_v9"), "unsupported_profile", "/profile"),
    "profile-boolean": (dict(profile=True), "malformed", "/profile"),
    "format-unknown": (dict(formats=("airr", "ungapped")), "unsupported_profile", "/formats/1"),
    "format-latest": (dict(formats=("latest",)), "unsupported_profile", "/formats/0"),
    "format-boolean": (dict(formats=(True,)), "malformed", "/formats/0"),
    "format-duplicate": (dict(formats=("airr", "airr")), "malformed", "/formats"),
    "format-empty": (dict(formats=()), "malformed", "/formats"),
    "format-string": (dict(formats="airr"), "malformed", "/formats"),
    "format-no-record-export": (dict(formats=("gapped", "gapped_ex")), "malformed", "/formats"),
}


@pytest.mark.parametrize("override,kind,path", list(BAD_REQUESTS.values()), ids=list(BAD_REQUESTS))
def test_bad_requests_are_rejected_before_any_transport_call_or_file_io(tmp_path, override, kind, path):
    request = dict(set_id=IGH[0], release=IGH[1], profile=PROFILE, formats=("airr", "airr_ex"))
    request.update(override)
    error = _error(kind, fetch_source, cache=tmp_path / "cache", transport=FailingTransport(), identify=identify,
                   **request)
    assert error.path == path
    assert not (tmp_path / "cache").exists()


def test_arguments_of_the_wrong_type_are_programming_errors(tmp_path):
    with pytest.raises(TypeError):
        fetch(tmp_path, FailingTransport(), identify="not callable")
    with pytest.raises(TypeError):
        fetch(tmp_path, FailingTransport(), expected_lock={"not": "a lock"})
    with pytest.raises(TypeError):
        fetch(tmp_path, FailingTransport(), terms={"status": "unresolved"})
    with pytest.raises(TypeError):
        fetch(tmp_path, object())
    assert not (tmp_path / "cache").exists()


# ── a clean fetch ────────────────────────────────────────────────────


def test_fetch_stores_every_requested_format_and_seals_a_candidate_lock(tmp_path):
    transport = FakeTransport(serve())
    seen = []

    def spy(**kwargs):
        seen.append(dict(kwargs["assets"]))
        return identify(**kwargs)

    result = fetch_source(set_id=IGH[0], release=IGH[1], profile=PROFILE, formats=ALL, cache=tmp_path / "cache",
                          transport=transport, identify=spy, now=lambda: T0)
    assert isinstance(result, AcquiredSource)
    lock = result.lock
    assert transport.calls == [(_url(IGH, fmt), SOCKET_TIMEOUT_SECONDS) for fmt in ALL]
    assert all(response.closed for response in transport.responses)
    assert seen == [{fmt: _fixture(IGH, fmt) for fmt in ALL}]
    assert (lock.provider, lock.set_id, lock.release) == ("ogrdb", IGH[0], IGH[1])
    assert lock.identity.to_dict() == {"taxon": {"scientific_name": "Homo sapiens", "taxonomy_id": 9606},
                                       "locus": "IGH", "chain_type": "vdj"}
    assert lock.adapter == resolve_profile(PROFILE)
    assert [a.to_dict() for a in lock.assets] == [
        {"asset_id": fmt, "role": m.FORMAT_ROLES[fmt], "format": fmt, "sha256": _sha(_fixture(IGH, fmt)),
         "bytes": len(_fixture(IGH, fmt)), "declared_url": _url(IGH, fmt), "origin_kind": "provider_download"}
        for fmt in ALL]
    assert lock.missing_assets == ()
    assert lock.terms.status == "unresolved" and lock.terms.reason
    m.check_profile_scope(lock)          # exactly the validated bytes
    assert dict(result.assets) == {fmt: _fixture(IGH, fmt) for fmt in ALL}
    cache = BlobCache(tmp_path / "cache")
    for asset in lock.assets:
        assert cache.read(asset.sha256, asset.bytes) == _fixture(IGH, asset.format)
    assert _temporaries(tmp_path) == []


def test_receipts_record_the_event_and_only_selected_headers(tmp_path):
    result = fetch(tmp_path, FakeTransport(serve(formats=("airr",))), formats=("airr",))
    (receipt,) = result.receipts
    asset = _sources_asset(IGH, "airr")
    assert receipt.to_dict() == {
        "artifact": "genairr_reference_source_receipt", "schema_version": 1,
        "event_time_utc": "2026-10-27T09:30:00.000000+00:00", "origin": "provider_download", "format": "airr",
        "sha256": "sha256:" + asset["sha256"], "bytes": asset["bytes"],
        "requested_url": asset["url"], "final_url": asset["url"], "redirects": [], "http_status": 200,
        "headers": {"content-disposition": asset["content_disposition"], "content-length": str(asset["bytes"]),
                    "content-type": "application/octet-stream"}}
    m.AcquisitionReceipt.from_dict(receipt.to_dict())


def test_unrecordable_header_values_are_left_out_of_the_receipt(tmp_path):
    routes = serve(formats=("airr",))
    body = _fixture(IGH, "airr")
    routes[_url(IGH, "airr")] = ok(body, headers={"ETag": "x" * 5000, "Date": " ", "Content-Type": "text/plain"})
    (receipt,) = fetch(tmp_path, FakeTransport(routes), formats=("airr",)).receipts
    assert dict(receipt.headers) == {"content-length": str(len(body)), "content-type": "text/plain"}


def test_default_formats_record_the_other_evidence_as_not_requested(tmp_path):
    transport = FakeTransport(serve())
    result = fetch_source(set_id=IGH[0], release=IGH[1], profile=PROFILE, cache=tmp_path / "c",
                          transport=transport, identify=identify)
    assert transport.urls == [_url(IGH, "airr"), _url(IGH, "airr_ex")]
    assert result.lock.asset_ids == ("airr", "airr_ex")
    assert [(a.format, a.reason) for a in result.lock.missing_assets] == [("gapped", "not_requested"),
                                                                          ("gapped_ex", "not_requested")]
    assert len(result.receipts) == 2


def test_requests_in_any_order_produce_the_canonical_lock(tmp_path):
    a = fetch(tmp_path, FakeTransport(serve()), formats=("gapped_ex", "airr", "gapped", "airr_ex"), cache="a")
    b = fetch(tmp_path, FakeTransport(serve()), formats=ALL, cache="b")
    assert a.lock == b.lock


def test_independent_clean_acquisitions_of_the_same_bytes_have_equal_lock_identity(tmp_path):
    first = fetch(tmp_path, FakeTransport(serve()), cache="first", now=lambda: T0)
    routes = {}
    for fmt in ALL:   # the second acquisition is redirected once and sees other headers and times
        hop = _url(IGH, fmt) + "/mirror"
        routes[_url(IGH, fmt)] = redirect(hop, status=307)
        routes[hop] = ok(_fixture(IGH, fmt), headers={"Content-Type": "application/json", "ETag": '"abc"',
                                                      "Date": "Tue, 27 Oct 2026 10:00:00 GMT"},
                         content_length=False, chunk=777)
    later = T0 + datetime.timedelta(days=3)
    second = fetch(tmp_path, FakeTransport(routes), cache="second", now=lambda: later)
    assert second.lock == first.lock
    assert second.lock.commitment == first.lock.commitment == m.lock_hash(first.lock)
    assert canonical_json(second.lock.to_dict()) == canonical_json(first.lock.to_dict())
    assert second.receipts != first.receipts
    assert all(r.redirects for r in second.receipts) and not any(r.redirects for r in first.receipts)
    assert {r.event_time_utc for r in second.receipts} != {r.event_time_utc for r in first.receipts}


def test_a_fetch_reacquires_every_requested_asset_even_when_cached(tmp_path):
    fetch(tmp_path, FakeTransport(serve()))
    transport = FakeTransport(serve())
    again = fetch(tmp_path, transport)
    assert transport.urls == [_url(IGH, fmt) for fmt in ALL]
    assert {r.origin for r in again.receipts} == {"provider_download"}


def test_a_fetch_is_attempted_once_and_never_retried(tmp_path):
    class Refused(FakeTransport):
        def get(self, url, *, timeout, deadline):
            self.calls.append((url, timeout))
            raise ConnectionRefusedError(errno.ECONNREFUSED, "Connection refused")

    transport = Refused()
    error = _error("network", fetch, tmp_path, transport)
    assert transport.urls == [_url(IGH, "airr")]
    assert _url(IGH, "airr") in error.message


def test_bytes_without_a_matching_identity_are_never_published(tmp_path):
    transport = FakeTransport(serve(source=IGK, at=IGH))     # another set behind the requested URLs
    error = _error("identity_mismatch", fetch, tmp_path, transport)
    assert "9606.IGKappa_VJ" in error.message
    assert _snapshot(tmp_path / "cache") == {"blobs": None, "blobs/sha256": None}


def test_a_wrong_release_is_rejected_by_the_identity_reader(tmp_path):
    transport = FakeTransport(serve(at=("9606.IGH_VDJ", 9)))  # release 10 bytes at the release 9 URLs
    _error("identity_mismatch", fetch, tmp_path, transport, source=("9606.IGH_VDJ", 9))
    assert _snapshot(tmp_path / "cache") == {"blobs": None, "blobs/sha256": None}


def test_an_identity_reader_must_return_a_source_identity(tmp_path):
    with pytest.raises(TypeError):
        fetch_source(set_id=IGH[0], release=IGH[1], profile=PROFILE, cache=tmp_path / "c",
                     transport=FakeTransport(serve()), identify=lambda **kwargs: {"locus": "IGH"})
    assert _temporaries(tmp_path) == [] and _snapshot(tmp_path / "c" / "blobs" / "sha256") == {}


# ── redirects ────────────────────────────────────────────────────────


@pytest.mark.parametrize("hops", [1, 2, 3])
def test_same_origin_https_redirects_are_followed_and_recorded(tmp_path, hops):
    start, body = _url(IGH, "airr"), _fixture(IGH, "airr")
    chain = [f"https://ogrdb.airr-community.org/mirror/{i}/airr" for i in range(hops)]
    routes = {start: redirect("/mirror/0/airr", status=301)}             # a relative Location
    for i in range(hops - 1):
        routes[chain[i]] = redirect(chain[i + 1], status=(302, 303, 308)[i % 3])
    routes[chain[-1]] = ok(body)
    result = fetch(tmp_path, FakeTransport(routes), formats=("airr",))
    (receipt,) = result.receipts
    assert receipt.requested_url == start
    assert list(receipt.redirects) == chain and receipt.final_url == chain[-1]
    assert result.lock.asset("airr").declared_url == start
    assert result.lock == fetch(tmp_path, FakeTransport(serve(formats=("airr",))), formats=("airr",),
                                cache="direct").lock


def test_more_than_three_redirects_are_refused(tmp_path):
    start = _url(IGH, "airr")
    chain = [f"https://ogrdb.airr-community.org/hop/{i}" for i in range(4)]
    routes = {start: redirect(chain[0])}
    for i in range(3):
        routes[chain[i]] = redirect(chain[i + 1])
    routes[chain[3]] = ok(_fixture(IGH, "airr"))
    transport = FakeTransport(routes)
    error = _error("network", fetch, tmp_path, transport, formats=("airr",))
    assert "redirect" in error.message and "3" in error.message
    assert transport.urls == [start] + chain[:3]
    assert _snapshot(tmp_path / "cache" / "blobs" / "sha256") == {}


@pytest.mark.parametrize("location", [
    "https://example.org/api_v2/germline/set/9606.IGH_VDJ/10/airr",
    "https://ogrdb.airr-community.org.example.org/airr",
    "https://mirror.ogrdb.airr-community.org/airr",
    "https://ogrdb.airr-community.org:8443/airr",
    "https://user:secret@ogrdb.airr-community.org/airr",
    "//example.org/airr",
    "ftp://ogrdb.airr-community.org/airr",
])
def test_a_redirect_to_another_origin_is_refused(tmp_path, location):
    transport = FakeTransport({_url(IGH, "airr"): redirect(location)})
    error = _error("network", fetch, tmp_path, transport, formats=("airr",))
    assert "secret" not in error.message
    assert transport.urls == [_url(IGH, "airr")]
    assert _snapshot(tmp_path / "cache" / "blobs" / "sha256") == {}


@pytest.mark.parametrize("location", ["http://ogrdb.airr-community.org/api_v2/germline/set/9606.IGH_VDJ/10/airr",
                                      "HTTP://ogrdb.airr-community.org/airr"])
def test_a_redirect_downgrade_to_http_is_refused(tmp_path, location):
    transport = FakeTransport({_url(IGH, "airr"): redirect(location)})
    error = _error("network", fetch, tmp_path, transport, formats=("airr",))
    assert "http" in error.message.lower() and "downgrade" in error.message
    assert transport.urls == [_url(IGH, "airr")]


@pytest.mark.parametrize("location", [
    "https://[::1/hop",                                    # urlsplit: invalid IPv6 URL
    "https://ogrdb.airr-community.org\uff0fevil/hop",      # netloc changes under NFKC normalization
    "https://ogrdb.airr-community.org/\u00e9",              # non-ASCII path
    "/\u00c3\u00a9",                                         # UTF-8 bytes as http.client decodes them
    "https://ogrdb.airr-community.org:abc/hop",
    "https://ogrdb.airr-community.org:99999/hop",
])
def test_a_malformed_or_non_ascii_redirect_target_is_a_network_error(tmp_path, location):
    transport = FakeTransport({_url(IGH, "airr"): redirect(location)})
    _error("network", fetch, tmp_path, transport, formats=("airr",))
    assert transport.urls == [_url(IGH, "airr")]
    assert _snapshot(tmp_path / "cache" / "blobs" / "sha256") == {}


def test_a_redirect_without_a_location_is_a_network_error(tmp_path):
    transport = FakeTransport({_url(IGH, "airr"): redirect(None)})
    _error("network", fetch, tmp_path, transport, formats=("airr",))


@pytest.mark.parametrize("status", [404, 500, 503, 204, 206, 304, 100])
def test_a_status_other_than_200_is_a_network_error(tmp_path, status):
    transport = FakeTransport({_url(IGH, "airr"): lambda: FakeResponse(status, b"body")})
    error = _error("network", fetch, tmp_path, transport, formats=("airr",))
    assert str(status) in error.message
    assert all(response.closed for response in transport.responses)


# ── truncated, oversized and slow responses ──────────────────────────


TRUNCATIONS = {
    "short-body": dict(content_length=len(_fixture(IGH, "airr")) + 10),
    "connection-reset": dict(fail_after=5000, error=ConnectionResetError(errno.ECONNRESET, "reset")),
    "socket-timeout": dict(fail_after=5000, error=socket.timeout("timed out")),
    "incomplete-read": dict(fail_after=5000, error=http.client.IncompleteRead(b"", 10)),
    "reset-without-length": dict(content_length=False, fail_after=5000,
                                 error=ConnectionResetError(errno.ECONNRESET, "reset")),
}


@pytest.mark.parametrize("options", list(TRUNCATIONS.values()), ids=list(TRUNCATIONS))
def test_a_truncated_download_is_a_network_error_and_publishes_nothing(tmp_path, options):
    routes = serve(formats=("airr", "airr_ex"))
    routes[_url(IGH, "airr_ex")] = ok(_fixture(IGH, "airr_ex"), **options)
    error = _error("network", fetch, tmp_path, FakeTransport(routes), formats=("airr", "airr_ex"))
    assert "airr_ex" in error.message
    assert _snapshot(tmp_path / "cache" / "blobs" / "sha256") == {}
    assert _temporaries(tmp_path) == []


def test_a_body_longer_than_its_content_length_is_refused(tmp_path):
    body = _fixture(IGH, "airr")
    routes = {_url(IGH, "airr"): ok(body, content_length=len(body) - 1)}
    error = _error("network", fetch, tmp_path, FakeTransport(routes), formats=("airr",))
    assert "Content-Length" in error.message
    assert _snapshot(tmp_path / "cache" / "blobs" / "sha256") == {}


@pytest.mark.parametrize("value", ["abc", "-1", "12, 12", "", "1e3", "0x10", "+5", "\u0661\u0662",
                                   pytest.param("9" * 20, id="20-digits"), pytest.param("9" * 5000, id="5000-digits")])
def test_an_invalid_content_length_is_a_network_error(tmp_path, value):
    routes = {_url(IGH, "airr"): ok(_fixture(IGH, "airr"), content_length=value)}
    _error("network", fetch, tmp_path, FakeTransport(routes), formats=("airr",))


def test_a_nineteen_digit_content_length_is_parsed_and_bounded(tmp_path):
    response = FakeResponse(200, b"x", content_length="9" * 19)
    error = _error("resource_limit", fetch, tmp_path, FakeTransport({_url(IGH, "airr"): lambda: response}),
                   formats=("airr",))
    assert "9" * 19 in error.message and response.reads == 0


class CloseFails(FakeResponse):
    def close(self):
        super().close()
        raise ConnectionResetError(errno.ECONNRESET, "reset while closing")


def test_an_error_closing_a_complete_response_is_a_network_error_and_leaves_nothing(tmp_path):
    routes = {_url(IGH, "airr"): lambda: CloseFails(200, _fixture(IGH, "airr"))}
    error = _error("network", fetch, tmp_path, FakeTransport(routes), formats=("airr",))
    assert "clos" in error.message
    assert _snapshot(tmp_path / "cache" / "blobs" / "sha256") == {}
    assert _temporaries(tmp_path) == []


@pytest.mark.parametrize("response,kind,text", [
    (lambda: CloseFails(500, b"oops"), "network", "500"),
    (lambda: CloseFails(200, filler=(b"x", m.MAX_ASSET_BYTES + 1)), "resource_limit", "Content-Length"),
    (lambda: CloseFails(302, b"", headers={"Location": "https://example.org/x"}), "network", "another origin"),
], ids=["status", "oversized", "redirect"])
def test_an_error_closing_a_failed_response_never_hides_the_failure(tmp_path, response, kind, text):
    error = _error(kind, fetch, tmp_path, FakeTransport({_url(IGH, "airr"): response}), formats=("airr",))
    assert text in error.message


def test_a_failure_after_the_body_is_staged_leaves_no_temporary_file(tmp_path):
    calls = []

    def now():
        calls.append(None)
        if len(calls) > 1:     # valid when checked up front, failing once the body is staged
            raise RuntimeError("clock failure")
        return T0

    with pytest.raises(RuntimeError, match="clock failure"):
        fetch(tmp_path, FakeTransport(serve(formats=("airr",))), formats=("airr",), now=now)
    assert len(calls) == 2
    assert _snapshot(tmp_path / "cache" / "blobs" / "sha256") == {}
    assert _temporaries(tmp_path) == []


def test_cleanup_failures_never_hide_the_original_error(tmp_path, monkeypatch):
    real_unlink = os.unlink

    def unlink(path, *args, **kwargs):
        if os.fspath(path).endswith(".tmp"):
            raise PermissionError(errno.EACCES, "The process cannot access the file")
        return real_unlink(path, *args, **kwargs)

    monkeypatch.setattr(os, "unlink", unlink)
    _error("identity_mismatch", fetch, tmp_path, FakeTransport(serve(source=IGK, at=IGH)))
    _error("network", fetch, tmp_path, FakeTransport({_url(IGH, "airr"): ok(_fixture(IGH, "airr"),
                                                                             content_length=10 ** 6)}),
           formats=("airr",))
    # only the temporaries that could not be removed remain; nothing was published
    assert all(name.endswith(".tmp") for name in _snapshot(tmp_path / "cache" / "blobs" / "sha256"))


def test_an_encoded_body_is_refused_because_raw_bytes_are_locked(tmp_path):
    routes = {_url(IGH, "airr"): ok(_fixture(IGH, "airr"), headers={"Content-Encoding": "gzip"})}
    error = _error("network", fetch, tmp_path, FakeTransport(routes), formats=("airr",))
    assert "gzip" in error.message


def test_an_empty_body_is_malformed(tmp_path):
    routes = {_url(IGH, "airr"): ok(b"")}
    _error("malformed", fetch, tmp_path, FakeTransport(routes), formats=("airr",))
    assert _temporaries(tmp_path) == []


def test_an_oversized_content_length_is_refused_before_reading_the_body(tmp_path):
    response = FakeResponse(200, filler=(b"x", m.MAX_ASSET_BYTES + 1))
    transport = FakeTransport({_url(IGH, "airr"): lambda: response})
    error = _error("resource_limit", fetch, tmp_path, transport, formats=("airr",))
    assert str(m.MAX_ASSET_BYTES + 1) in error.message and str(m.MAX_ASSET_BYTES) in error.message
    assert response.reads == 0 and response.closed
    assert _temporaries(tmp_path) == []


def test_an_oversized_body_without_content_length_is_refused_while_streaming(tmp_path):
    response = FakeResponse(200, filler=(b"x", m.MAX_ASSET_BYTES + 5_000_000), content_length=False,
                            chunk=1 << 20)
    transport = FakeTransport({_url(IGH, "airr"): lambda: response})
    error = _error("resource_limit", fetch, tmp_path, transport, formats=("airr",))
    assert str(m.MAX_ASSET_BYTES) in error.message
    assert response.served == m.MAX_ASSET_BYTES + 1          # never read beyond one byte past the limit
    assert _snapshot(tmp_path / "cache" / "blobs" / "sha256") == {}
    assert _temporaries(tmp_path) == []


def test_a_content_length_is_never_trusted_alone(tmp_path):
    # A small declared size does not raise the limit: the body is cut at the declared size.
    response = FakeResponse(200, filler=(b"x", m.MAX_ASSET_BYTES + 100), content_length=10)
    transport = FakeTransport({_url(IGH, "airr"): lambda: response})
    _error("network", fetch, tmp_path, transport, formats=("airr",))
    assert response.served == 11


@pytest.mark.parametrize("declared", [True, False], ids=["with-content-length", "without-content-length"])
def test_the_aggregate_payload_limit_spans_the_whole_operation(tmp_path, monkeypatch, declared):
    first, second = _fixture(IGH, "airr"), _fixture(IGH, "airr_ex")
    allowed = len(first) + len(second) - 1
    monkeypatch.setattr(acq, "MAX_SOURCE_PAYLOAD_BYTES", allowed)
    routes = serve(formats=("airr", "airr_ex"), content_length=declared)
    error = _error("resource_limit", fetch, tmp_path, FakeTransport(routes), formats=("airr", "airr_ex"))
    assert "aggregate" in error.message and str(allowed) in error.message
    assert _snapshot(tmp_path / "cache" / "blobs" / "sha256") == {}
    assert _temporaries(tmp_path) == []


class FakeClock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now


def test_the_acquisition_deadline_is_checked_while_streaming(tmp_path):
    clock = FakeClock()

    def slow():
        clock.now += 50.0

    routes = {_url(IGH, "airr"): ok(_fixture(IGH, "airr"), chunk=1000, on_read=slow)}
    error = _error("resource_limit", fetch, tmp_path, FakeTransport(routes), formats=("airr",), clock=clock)
    assert "deadline" in error.message and "150" in error.message and "120" in error.message
    assert _snapshot(tmp_path / "cache" / "blobs" / "sha256") == {}
    assert _temporaries(tmp_path) == []


def test_the_absolute_deadline_is_passed_to_the_transport(tmp_path):
    clock = FakeClock()
    transport = FakeTransport(serve(formats=("airr", "airr_ex")))
    fetch(tmp_path, transport, formats=("airr", "airr_ex"), clock=clock)
    assert transport.deadlines == [clock.now + ACQUISITION_DEADLINE_SECONDS] * 2


def test_the_deadline_covers_the_whole_operation(tmp_path):
    clock = FakeClock()

    class Slow(FakeTransport):
        def get(self, url, *, timeout, deadline):
            clock.now += 70.0
            return super().get(url, timeout=timeout, deadline=deadline)

    transport = Slow(serve(formats=("airr", "airr_ex")))
    _error("resource_limit", fetch, tmp_path, transport, formats=("airr", "airr_ex"), clock=clock)
    assert _snapshot(tmp_path / "cache" / "blobs" / "sha256") == {}


def test_cancellation_mid_stream_leaves_no_temporary_file_and_publishes_nothing(tmp_path):
    routes = serve(formats=("airr", "airr_ex"))
    routes[_url(IGH, "airr_ex")] = ok(_fixture(IGH, "airr_ex"), fail_after=4096, error=KeyboardInterrupt())
    with pytest.raises(KeyboardInterrupt):
        fetch(tmp_path, FakeTransport(routes), formats=("airr", "airr_ex"))
    assert _snapshot(tmp_path / "cache" / "blobs" / "sha256") == {}
    assert _temporaries(tmp_path) == []


# ── expected locks: reacquisition and drift ──────────────────────────


def _lock_file(tmp_path, lock) -> pathlib.Path:
    path = tmp_path / "source.lock.json"
    path.write_bytes(canonical_json(lock.to_dict()).encode("utf-8"))
    return path


def test_an_expected_lock_reacquires_cached_assets_and_reproduces_the_lock(tmp_path):
    lock = fetch(tmp_path, FakeTransport(serve())).lock
    transport = FakeTransport(serve())
    again = fetch(tmp_path, transport, expected_lock=lock, now=lambda: T0 + datetime.timedelta(hours=1))
    assert again.lock == lock
    assert transport.urls == [_url(IGH, fmt) for fmt in ALL]
    assert [r.origin for r in again.receipts] == ["provider_download"] * 4


def test_unrequested_pinned_assets_are_read_from_the_verified_cache(tmp_path):
    lock = fetch(tmp_path, FakeTransport(serve())).lock
    transport = FakeTransport(serve())
    again = fetch(tmp_path, transport, formats=("airr", "airr_ex"), expected_lock=lock)
    assert again.lock == lock
    assert transport.urls == [_url(IGH, "airr"), _url(IGH, "airr_ex")]
    assert [(r.format, r.origin) for r in again.receipts] == [
        ("airr", "provider_download"), ("airr_ex", "provider_download"),
        ("gapped", "verified_cache"), ("gapped_ex", "verified_cache")]
    assert again.receipts[2].requested_url is None and again.receipts[2].headers == ()


def test_a_pinned_asset_missing_from_the_cache_fails_before_any_http(tmp_path):
    lock = fetch(tmp_path, FakeTransport(serve())).lock
    BlobCache(tmp_path / "cache").blob_path(lock.asset("gapped").sha256).unlink()
    error = _error("missing_asset", fetch, tmp_path, FailingTransport(), formats=("airr", "airr_ex"),
                   expected_lock=lock)
    assert "gapped" in error.message and lock.asset("gapped").sha256 in error.message


def _with(lock, **changes):
    fields = dict(provider=lock.provider, set_id=lock.set_id, release=lock.release, identity=lock.identity,
                  adapter=lock.adapter, assets=lock.assets, missing_assets=lock.missing_assets, terms=lock.terms)
    fields.update(changes)
    return m.SourceLock.seal(**fields)


def test_an_expected_lock_must_match_the_request_before_any_http(tmp_path):
    lock = fetch(tmp_path, FakeTransport(serve(formats=("airr", "airr_ex"))), formats=("airr", "airr_ex")).lock
    mismatches = [
        dict(source=IGK),                               # another set
        dict(source=("9606.IGH_VDJ", 9)),               # another release
        dict(formats=("airr", "airr_ex", "gapped")),    # a format the lock does not pin
    ]
    for request in mismatches:
        before = _snapshot(tmp_path / "cache")
        error = _error("identity_mismatch", fetch, tmp_path, FailingTransport(), expected_lock=lock,
                       **{"formats": ("airr", "airr_ex"), **request})
        assert error.path.startswith("/expected_lock")
        assert _snapshot(tmp_path / "cache") == before


def test_a_local_or_relocated_expected_asset_cannot_be_reacquired(tmp_path):
    lock = fetch(tmp_path, FakeTransport(serve(formats=("airr",))), formats=("airr",)).lock
    local = _with(lock, assets=(m.Asset(**{**lock.asset("airr").__dict__, "origin_kind": "local_file",
                                           "declared_url": None}),))
    moved = _with(lock, assets=(m.Asset(**{**lock.asset("airr").__dict__,
                                           "declared_url": "https://ogrdb.airr-community.org/elsewhere"}),))
    for expected in (local, moved):
        _error("identity_mismatch", fetch, tmp_path, FailingTransport(), formats=("airr",), expected_lock=expected)


def test_an_expected_lock_requires_the_same_identity(tmp_path):
    lock = fetch(tmp_path, FakeTransport(serve(formats=("airr",))), formats=("airr",)).lock
    other = _with(lock, identity=m.SourceIdentity(m.Taxon("Homo sapiens", 9606), "IGK", "vj"))
    before = _snapshot(tmp_path / "cache")
    error = _error("identity_mismatch", fetch, tmp_path, FakeTransport(serve(formats=("airr",))),
                   formats=("airr",), expected_lock=other)
    assert "IGK" in error.message and "IGH" in error.message
    assert _snapshot(tmp_path / "cache") == before


def test_terms_of_an_expected_lock_come_from_the_lock(tmp_path):
    lock = fetch(tmp_path, FakeTransport(serve(formats=("airr",))), formats=("airr",)).lock
    other_terms = m.Terms("unresolved", (), None, None, "another reason")
    _error("malformed", fetch, tmp_path, FailingTransport(), formats=("airr",), expected_lock=lock,
           terms=other_terms)
    assert fetch(tmp_path, FakeTransport(serve(formats=("airr",))), formats=("airr",), expected_lock=lock,
                 terms=lock.terms).lock == lock


DRIFTS = {
    "same-size": lambda data: data[:-2] + b"X\n",
    "longer": lambda data: data + b"\n",
    "shorter": lambda data: data[:-1],
}


@pytest.mark.parametrize("change", list(DRIFTS.values()), ids=list(DRIFTS))
def test_changing_upstream_bytes_are_reported_as_drift_and_change_nothing(tmp_path, change):
    lock = fetch(tmp_path, FakeTransport(serve())).lock
    lock_path = _lock_file(tmp_path, lock)
    before = _snapshot(tmp_path)
    original = _fixture(IGH, "airr_ex")
    drifted = change(original)
    routes = serve()
    routes[_url(IGH, "airr_ex")] = ok(drifted)
    error = _error("digest_mismatch", fetch, tmp_path, FakeTransport(routes), expected_lock=lock)
    assert "airr_ex" in error.message and "drift" in error.message
    assert _sha(drifted) in error.message and str(len(drifted)) in error.message
    assert _sha(original) in error.message and str(len(original)) in error.message
    assert _snapshot(tmp_path) == before                     # old lock and cache intact
    assert not BlobCache(tmp_path / "cache").blob_path(_sha(drifted)).exists()
    assert m.SourceLock.read(lock_path) == lock


def test_without_an_expected_lock_changed_bytes_make_a_new_candidate_lock(tmp_path):
    lock = fetch(tmp_path, FakeTransport(serve())).lock
    routes = serve()
    routes[_url(IGH, "gapped")] = ok(_fixture(IGH, "gapped") + b"\n")
    candidate = fetch(tmp_path, FakeTransport(routes)).lock
    assert candidate.commitment != lock.commitment
    assert m.profile_scope_problem(candidate) is not None and m.profile_scope_problem(lock) is None


# ── failures never replace a valid cache blob or lock ────────────────


def _failing_routes():
    def drift():
        routes = serve()
        routes[_url(IGH, "gapped")] = ok(_fixture(IGH, "gapped")[:-1] + b"N")
        return routes

    def truncated():
        routes = serve()
        routes[_url(IGH, "gapped_ex")] = ok(_fixture(IGH, "gapped_ex"), content_length=10 ** 6)
        return routes

    def oversized():
        routes = serve()
        routes[_url(IGH, "airr_ex")] = lambda: FakeResponse(200, filler=(b"x", m.MAX_ASSET_BYTES + 1))
        return routes

    def other_origin():
        routes = serve()
        routes[_url(IGH, "gapped")] = redirect("https://example.org/gapped")
        return routes

    def server_error():
        routes = serve()
        routes[_url(IGH, "gapped_ex")] = lambda: FakeResponse(500, b"oops")
        return routes

    return {"drift": (drift, "digest_mismatch"), "truncated": (truncated, "network"),
            "oversized": (oversized, "resource_limit"), "other-origin": (other_origin, "network"),
            "server-error": (server_error, "network"),
            "other-set": (lambda: serve(source=IGK, at=IGH), "identity_mismatch")}


FAILURES = _failing_routes()
# Changed bytes are a failure only against an expected lock; without one they seal a new candidate.
FAILURE_CASES = [pytest.param(case, pinned, id=f"{case}-{'expected-lock' if pinned else 'candidate'}")
                 for case in FAILURES for pinned in (True, False) if pinned or case != "drift"]


@pytest.mark.parametrize("case,with_expected_lock", FAILURE_CASES)
def test_failures_never_replace_a_valid_cache_blob_or_lock(tmp_path, case, with_expected_lock):
    routes, kind = FAILURES[case]
    lock = fetch(tmp_path, FakeTransport(serve())).lock
    _lock_file(tmp_path, lock)
    before = _snapshot(tmp_path)
    expected = dict(expected_lock=lock) if with_expected_lock else {}
    if kind == "identity_mismatch" and with_expected_lock:
        kind = "digest_mismatch"          # the other set's bytes are drift against the pins
    _error(kind, fetch, tmp_path, FakeTransport(routes()), **expected)
    assert _snapshot(tmp_path) == before


def test_a_conflicting_cache_file_fails_the_fetch_and_publishes_nothing(tmp_path):
    cache = BlobCache(tmp_path / "cache")
    cache.blob_dir.mkdir(parents=True)
    conflicting = cache.blob_path(_sha(_fixture(IGH, "gapped_ex")))
    conflicting.write_bytes(b"not the gapped_ex bytes")
    error = _error("digest_mismatch", fetch, tmp_path, FakeTransport(serve()))
    assert str(conflicting) in error.message
    assert sorted(p.name for p in cache.blob_dir.iterdir()) == [conflicting.name]
    assert conflicting.read_bytes() == b"not the gapped_ex bytes"


def test_concurrent_fetches_into_one_cache_agree(tmp_path):
    barrier = threading.Barrier(2)
    results, errors = [], []

    class Synchronized(FakeTransport):
        def get(self, url, *, timeout, deadline):
            if url == _url(IGH, "airr"):
                barrier.wait(timeout=10)
            return super().get(url, timeout=timeout, deadline=deadline)

    def run():
        try:
            results.append(fetch(tmp_path, Synchronized(serve())))
        except BaseException as exc:  # reported below
            errors.append(exc)

    threads = [threading.Thread(target=run) for _ in range(2)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(timeout=60)
    assert errors == []
    assert results[0].lock == results[1].lock
    assert len(list((tmp_path / "cache" / "blobs" / "sha256").iterdir())) == 4
    assert _temporaries(tmp_path) == []


def test_a_cache_without_hard_links_fails_typed_and_creates_no_blob(tmp_path, monkeypatch):
    monkeypatch.setattr(_artifact_io, "_WINDOWS", False)

    def link(source, destination, *args, **kwargs):
        raise PermissionError(errno.EPERM, "Operation not permitted")

    monkeypatch.setattr(os, "link", link)
    error = _error("io", fetch, tmp_path, FakeTransport(serve()))
    assert "atomic no-overwrite publication" in error.message
    assert _snapshot(tmp_path / "cache" / "blobs" / "sha256") == {}


# ── data terms ───────────────────────────────────────────────────────


LICENCE = b"Information provided on this site is published under a CC0 (Creative Commons 0) license.\n"


def _recorded_terms():
    return m.Terms("recorded", ("licence_statement",), "CC0 1.0 Universal", None, None)


def test_recorded_terms_come_with_their_supplied_statement(tmp_path):
    statement = tmp_path / "statement.txt"
    statement.write_bytes(LICENCE)
    transport = FakeTransport(serve(formats=("airr",)))
    result = fetch(tmp_path, transport, formats=("airr",), terms=_recorded_terms(),
                   statements={"licence_statement": statement})
    assert transport.urls == [_url(IGH, "airr")]               # statements are never fetched
    assert result.lock.terms == _recorded_terms()
    asset = result.lock.asset("licence_statement")
    assert (asset.origin_kind, asset.declared_url, asset.sha256) == ("local_file", None, _sha(LICENCE))
    assert [(r.format, r.origin) for r in result.receipts] == [("airr", "provider_download"),
                                                                ("licence_statement", "local_file")]
    assert BlobCache(tmp_path / "cache").read(asset.sha256, asset.bytes) == LICENCE


def test_terms_and_statements_must_agree_before_any_transport_call(tmp_path):
    statement = tmp_path / "statement.txt"
    statement.write_bytes(LICENCE)
    _error("malformed", fetch, tmp_path, FailingTransport(), terms=_recorded_terms())
    _error("malformed", fetch, tmp_path, FailingTransport(), statements={"licence_statement": statement})
    _error("malformed", fetch, tmp_path, FailingTransport(), terms=_recorded_terms(),
           statements={"airr": statement})
    assert not (tmp_path / "cache").exists()


# ── local files ──────────────────────────────────────────────────────


def _local_copies(directory: pathlib.Path, source=IGH, formats=ALL) -> dict:
    directory.mkdir(parents=True, exist_ok=True)
    files = {}
    for fmt in formats:
        path = directory / f"download-{fmt}.bin"
        path.write_bytes(_fixture(source, fmt))
        files[fmt] = path
    return files


def ingest(tmp_path, files, *, source=IGH, cache="cache", **kwargs):
    return ingest_files(set_id=source[0], release=source[1], profile=PROFILE, files=files,
                        cache=tmp_path / cache, identify=identify, **kwargs)


def test_local_files_are_recorded_as_supplied_bytes(tmp_path, no_network):
    result = ingest(tmp_path, _local_copies(tmp_path / "exports"), now=lambda: T0)
    assert [(a.format, a.origin_kind, a.declared_url, a.sha256) for a in result.lock.assets] == [
        (fmt, "local_file", None, _sha(_fixture(IGH, fmt))) for fmt in ALL]
    assert [r.to_dict() for r in result.receipts] == [
        {"artifact": "genairr_reference_source_receipt", "schema_version": 1,
         "event_time_utc": "2026-10-27T09:30:00.000000+00:00", "origin": "local_file", "format": fmt,
         "sha256": _sha(_fixture(IGH, fmt)), "bytes": len(_fixture(IGH, fmt)), "requested_url": None,
         "final_url": None, "redirects": [], "http_status": None, "headers": {}} for fmt in ALL]
    cache = BlobCache(tmp_path / "cache")
    assert all(cache.read(a.sha256, a.bytes) == _fixture(IGH, a.format) for a in result.lock.assets)
    m.check_profile_scope(result.lock)


def test_identical_local_bytes_moved_between_directories_yield_an_identical_lock(tmp_path, no_network):
    first = ingest(tmp_path, _local_copies(tmp_path / "a" / "deep"), cache="c1", now=lambda: T0)
    second = ingest(tmp_path, {fmt: str(path) for fmt, path in _local_copies(tmp_path / "b").items()},
                    cache="c2", now=lambda: T0 + datetime.timedelta(days=1))
    assert second.lock == first.lock
    assert canonical_json(second.lock.to_dict()) == canonical_json(first.lock.to_dict())
    assert "download-" not in canonical_json(first.lock.to_dict())
    assert str(tmp_path) not in canonical_json(first.lock.to_dict())


def test_local_and_provider_locks_over_the_same_bytes_differ(tmp_path):
    local = ingest(tmp_path, _local_copies(tmp_path / "exports"), cache="local").lock
    provider = fetch(tmp_path, FakeTransport(serve()), cache="provider").lock
    assert [a.sha256 for a in local.assets] == [a.sha256 for a in provider.assets]
    assert local.commitment != provider.commitment


def test_a_declared_provider_url_is_recorded_without_claiming_a_fetch(tmp_path, no_network):
    files = _local_copies(tmp_path / "exports", formats=("airr",))
    result = ingest(tmp_path, {"airr": SuppliedFile(files["airr"], declared_url=_url(IGH, "airr"))})
    asset = result.lock.asset("airr")
    assert (asset.origin_kind, asset.declared_url) == ("local_file", _url(IGH, "airr"))
    (receipt,) = result.receipts
    assert (receipt.origin, receipt.requested_url, receipt.final_url) == ("local_file", None, None)
    assert [(a.format, a.reason) for a in result.lock.missing_assets] == [
        ("airr_ex", "not_supplied"), ("gapped", "not_supplied"), ("gapped_ex", "not_supplied")]


def test_local_ingestion_accepts_statements_with_recorded_terms(tmp_path, no_network):
    files = _local_copies(tmp_path / "exports", formats=("airr",))
    statement = tmp_path / "exports" / "licence.txt"
    statement.write_bytes(LICENCE)
    files["licence_statement"] = statement
    result = ingest(tmp_path, files, terms=_recorded_terms())
    assert result.lock.asset_ids == ("airr", "licence_statement")
    assert result.lock.terms == _recorded_terms()


def _bad_local_inputs(tmp_path):
    exports = _local_copies(tmp_path / "exports", formats=("airr", "gapped"))
    empty = tmp_path / "empty.json"
    empty.write_bytes(b"")
    return {
        "missing-file": ({"airr": tmp_path / "absent.json"}, "io"),
        "directory": ({"airr": tmp_path / "exports"}, "io"),
        "empty-file": ({"airr": empty}, "malformed"),
        "unknown-format": ({"airr": exports["airr"], "fasta": exports["gapped"]}, "malformed"),
        "no-record-export": ({"gapped": exports["gapped"]}, "malformed"),
        "no-files": ({}, "malformed"),
        "other-release-url": ({"airr": SuppliedFile(exports["airr"], declared_url=_url(("9606.IGH_VDJ", 9), "airr"))},
                              "identity_mismatch"),
        "other-format-url": ({"airr": SuppliedFile(exports["airr"], declared_url=_url(IGH, "airr_ex"))},
                             "identity_mismatch"),
        "off-origin-url": ({"airr": SuppliedFile(exports["airr"], declared_url="https://example.org/airr.json")},
                           "identity_mismatch"),
        "terms-without-statement": ({"airr": exports["airr"]}, "malformed"),
        "wrong-set": ({"airr": _local_copies(tmp_path / "igk", source=IGK, formats=("airr",))["airr"]},
                      "identity_mismatch"),
    }


@pytest.mark.parametrize("case", ["missing-file", "directory", "empty-file", "unknown-format", "no-record-export",
                                  "no-files", "other-release-url", "other-format-url", "off-origin-url",
                                  "terms-without-statement", "wrong-set"])
def test_bad_local_input_fails_typed_offline_and_publishes_nothing(tmp_path, no_network, case):
    files, kind = _bad_local_inputs(tmp_path)[case]
    terms = _recorded_terms() if case == "terms-without-statement" else None
    _error(kind, ingest, tmp_path, files, terms=terms)
    assert _snapshot(tmp_path / "cache" / "blobs" / "sha256") == {}
    assert _temporaries(tmp_path) == []


@pytest.mark.parametrize("url", ["http://ogrdb.airr-community.org/api_v2/germline/set/9606.IGH_VDJ/10/airr",
                                 "https://user:pw@ogrdb.airr-community.org/x", "not a url"])
def test_a_declared_url_is_an_https_url_without_credentials(tmp_path, url):
    error = _error("malformed", SuppliedFile, tmp_path / "airr.json", declared_url=url)
    assert error.path == "/declared_url"
    with pytest.raises(TypeError):
        SuppliedFile(b"airr.json")


def test_a_declared_statement_url_must_be_on_the_provider_origin(tmp_path, no_network):
    files = _local_copies(tmp_path / "exports", formats=("airr",))
    statement = tmp_path / "licence.txt"
    statement.write_bytes(LICENCE)
    page = "https://ogrdb.airr-community.org/render_page/licensing_statement.html"
    result = ingest(tmp_path, {**files, "licence_statement": SuppliedFile(statement, declared_url=page)},
                    terms=_recorded_terms())
    assert result.lock.asset("licence_statement").declared_url == page
    _error("identity_mismatch", ingest, tmp_path,
           {**files, "licence_statement": SuppliedFile(statement, declared_url="https://example.org/terms")},
           terms=_recorded_terms(), cache="other")


@pytest.mark.parametrize("port", ["abc", "99999"])
def test_an_invalid_port_in_a_declared_statement_url_is_malformed(tmp_path, no_network, port):
    files = _local_copies(tmp_path / "exports", formats=("airr",))
    statement = tmp_path / "licence.txt"
    statement.write_bytes(LICENCE)
    supplied = SuppliedFile(statement, declared_url=f"https://ogrdb.airr-community.org:{port}/terms")
    error = _error("malformed", ingest, tmp_path, {**files, "licence_statement": supplied}, terms=_recorded_terms())
    assert error.path == "/files/licence_statement/declared_url"
    error = _error("malformed", fetch, tmp_path, FailingTransport(), formats=("airr",), terms=_recorded_terms(),
                   statements={"licence_statement": supplied})
    assert error.path == "/statements/licence_statement/declared_url"
    assert not (tmp_path / "cache").exists()


def test_an_oversized_local_file_is_a_resource_limit(tmp_path, no_network):
    big = tmp_path / "big.json"
    with open(big, "wb") as handle:
        handle.truncate(m.MAX_ASSET_BYTES + 1)
    error = _error("resource_limit", ingest, tmp_path, {"airr": big})
    assert str(m.MAX_ASSET_BYTES + 1) in error.message and str(m.MAX_ASSET_BYTES) in error.message
    assert _temporaries(tmp_path) == []


def test_local_ingestion_has_the_same_aggregate_budget(tmp_path, monkeypatch, no_network):
    files = _local_copies(tmp_path / "exports", formats=("airr", "airr_ex"))
    monkeypatch.setattr(acq, "MAX_SOURCE_PAYLOAD_BYTES", len(_fixture(IGH, "airr")) + 10)
    error = _error("resource_limit", ingest, tmp_path, files)
    assert "aggregate" in error.message
    assert _snapshot(tmp_path / "cache" / "blobs" / "sha256") == {}


# ── offline reads ────────────────────────────────────────────────────


NAIVE = datetime.datetime(2026, 10, 27, 9, 30)


@pytest.mark.parametrize("now", [lambda: NAIVE, lambda: "2026-10-27T09:30:00Z", "not callable"],
                         ids=["naive", "text", "not-callable"])
def test_now_is_validated_before_any_io(tmp_path, now):
    with pytest.raises(TypeError):
        fetch(tmp_path, FailingTransport(), now=now)
    with pytest.raises(TypeError):
        ingest(tmp_path, {"airr": tmp_path / "absent.json"}, now=now)
    with pytest.raises(TypeError):
        read_locked(tmp_path / "absent.lock.json", cache=tmp_path / "cache", now=now)
    assert not (tmp_path / "cache").exists()


def test_offline_operations_take_no_transport():
    for operation in (read_locked, ingest_files):
        parameters = inspect.signature(operation).parameters
        assert "transport" not in parameters and "clock" not in parameters


def test_reading_a_lock_verifies_every_blob_offline(tmp_path, monkeypatch):
    lock = fetch(tmp_path, FakeTransport(serve())).lock
    path = _lock_file(tmp_path, lock)
    monkeypatch.setattr(acq, "fetch_source", lambda **kwargs: pytest.fail("read tried to fetch"))
    for given in (lock, path, str(path)):
        result = read_locked(given, cache=tmp_path / "cache", now=lambda: T0)
        assert result.lock == lock
        assert dict(result.assets) == {fmt: _fixture(IGH, fmt) for fmt in ALL}
        assert [(r.format, r.origin, r.requested_url) for r in result.receipts] == [
            (fmt, "verified_cache", None) for fmt in ALL]
    assert read_locked(lock, cache=BlobCache(tmp_path / "cache"), identify=identify).lock == lock


def test_offline_reads_never_touch_the_network(tmp_path, no_network):
    lock = ingest(tmp_path, _local_copies(tmp_path / "exports")).lock
    path = _lock_file(tmp_path, lock)
    assert read_locked(path, cache=tmp_path / "cache").lock == lock
    # a cache miss
    error = _error("missing_asset", read_locked, path, cache=tmp_path / "empty-cache")
    for asset in lock.assets:
        assert asset.format in error.message and asset.sha256 in error.message
    assert not (tmp_path / "empty-cache").exists()
    # malformed input
    garbage = tmp_path / "garbage.lock.json"
    garbage.write_bytes(b'{"artifact": "genairr_reference_source_lock", "schema_version": 1}')
    _error("malformed", read_locked, garbage, cache=tmp_path / "cache")
    _error("io", read_locked, tmp_path / "absent.lock.json", cache=tmp_path / "cache")
    _error("malformed", ingest, tmp_path, {"airr": tmp_path / "exports"}, source=("9606.IGH_VDJ", "latest"))


def test_a_missing_blob_names_its_asset_and_digest_and_never_downloads(tmp_path, no_network):
    lock = ingest(tmp_path, _local_copies(tmp_path / "exports")).lock
    gapped = lock.asset("gapped")
    BlobCache(tmp_path / "cache").blob_path(gapped.sha256).unlink()
    error = _error("missing_asset", read_locked, lock, cache=tmp_path / "cache")
    assert error.path == "/assets/2"
    assert "gapped" in error.message and gapped.sha256 in error.message and str(gapped.bytes) in error.message
    assert "airr_ex" not in error.message


def test_a_corrupt_cache_blob_fails_an_offline_read(tmp_path, no_network):
    lock = ingest(tmp_path, _local_copies(tmp_path / "exports")).lock
    asset = lock.asset("airr_ex")
    path = BlobCache(tmp_path / "cache").blob_path(asset.sha256)
    corrupted = bytearray(path.read_bytes())
    corrupted[100] ^= 1
    path.write_bytes(bytes(corrupted))
    error = _error("digest_mismatch", read_locked, lock, cache=tmp_path / "cache")
    assert error.path == "/assets/1" and asset.sha256 in error.message
    assert path.read_bytes() == bytes(corrupted)


def test_an_offline_read_can_check_the_identity_again(tmp_path, no_network):
    lock = ingest(tmp_path, _local_copies(tmp_path / "exports")).lock
    other = _with(lock, identity=m.SourceIdentity(m.Taxon("Homo sapiens", 9606), "IGK", "vj"))
    _error("identity_mismatch", read_locked, other, cache=tmp_path / "cache", identify=identify)


def test_assets_with_identical_bytes_share_one_verified_blob(tmp_path, no_network):
    exports = _local_copies(tmp_path / "exports", formats=("airr",))
    statement = tmp_path / "terms.txt"
    statement.write_bytes(LICENCE)
    terms = m.Terms("recorded", ("licence_statement", "citation_statement"), "CC0 1.0 Universal", None, None)
    result = ingest(tmp_path, {**exports, "licence_statement": statement, "citation_statement": statement},
                    terms=terms)
    assert result.lock.asset("licence_statement").sha256 == result.lock.asset("citation_statement").sha256
    assert read_locked(result.lock, cache=tmp_path / "cache").assets["citation_statement"] == LICENCE


# ── the standard-library HTTPS transport ─────────────────────────────
#
# A stand-in for a connected TLS socket drives the real http.client request and response
# code, and a simulated monotonic clock advances as the "server" sends bytes: no network
# and no real waiting.


class SimulatedSocket:
    """A connected socket whose receives cost simulated time.

    Each receive yields at most ``chunk`` bytes and costs ``head_step`` seconds while the
    status line and headers are sent, then ``body_step``.  A receive slower than the
    current socket timeout advances the clock by the timeout and raises
    ``socket.timeout``, as a real socket would.
    """

    def __init__(self, clock, data, head_size, *, head_step=0.0, body_step=0.0, chunk=1 << 16):
        self.clock, self.data, self.head_size = clock, data, head_size
        self.head_step, self.body_step, self.chunk = head_step, body_step, chunk
        self.pos = self.recvs = 0
        self.timeout = None
        self.timeouts = []        # (simulated time, timeout) at each settimeout
        self.sent = b""
        self.closed = False

    def settimeout(self, value):
        self.timeout = value
        self.timeouts.append((self.clock.now, value))

    def gettimeout(self):
        return self.timeout

    def sendall(self, data, *args):
        self.sent += bytes(data)

    def makefile(self, mode, buffering=None, **kwargs):
        assert (mode, buffering) == ("rb", 0)     # the transport adds its own buffering
        return SimulatedReader(self)

    def recv_into(self, buffer):
        if self.pos >= len(self.data):
            return 0
        in_head = self.pos < self.head_size
        cost = self.head_step if in_head else self.body_step
        if self.timeout is not None and cost > self.timeout:
            self.clock.now += self.timeout
            raise socket.timeout("timed out")
        self.clock.now += cost
        end = min(len(self.data), self.pos + self.chunk, self.pos + len(buffer),
                  self.head_size if in_head else len(self.data))
        n = end - self.pos
        buffer[:n] = self.data[self.pos:end]
        self.pos, self.recvs = end, self.recvs + 1
        return n

    def close(self):
        self.closed = True


class SimulatedReader(io.RawIOBase):
    def __init__(self, sock):
        self.sock = sock

    def readable(self):
        return True

    def readinto(self, buffer):
        return self.sock.recv_into(buffer)


def http_response(body: bytes, *, status="200 OK", headers=()):
    head = (f"HTTP/1.1 {status}\r\nContent-Length: {len(body)}\r\n"
            + "".join(f"{name}: {value}\r\n" for name, value in headers) + "\r\n").encode("ascii")
    return head + body, len(head)


def simulate_connections(monkeypatch, clock, make_socket, *, connect_cost=0.0):
    """Replace ``http.client.HTTPSConnection`` by one whose ``connect`` yields simulated sockets."""
    opened = []

    # A plain HTTPConnection underneath: the simulation needs no TLS, and the patched module
    # name must not be the base class (HTTPSConnection.__init__ looks itself up by name).
    class Connection(http.client.HTTPConnection):
        default_port = http.client.HTTPS_PORT

        def __init__(self, host, port=None, *, timeout=None, context=None):
            self.given_context = context
            super().__init__(host, port, timeout=timeout)

        def connect(self):
            opened.append(self)
            if connect_cost > self.timeout:
                clock.now += self.timeout
                raise socket.timeout("timed out")
            clock.now += connect_cost
            self.simulated = make_socket()
            self.sock = self.simulated

    monkeypatch.setattr(http.client, "HTTPSConnection", Connection)
    return opened


def test_the_https_transport_drives_the_standard_library_client(monkeypatch):
    clock = FakeClock()
    data, head = http_response(b"payload", headers=[("Content-Type", "text/plain")])
    opened = simulate_connections(monkeypatch, clock, lambda: SimulatedSocket(clock, data, head))
    response = HttpsTransport(clock=clock).get(_url(IGH, "airr") + "?q=1", timeout=SOCKET_TIMEOUT_SECONDS,
                                               deadline=clock.now + ACQUISITION_DEADLINE_SECONDS)
    (connection,) = opened
    assert (connection.host, connection.port, connection.timeout) == ("ogrdb.airr-community.org", 443, 30)
    lines = connection.simulated.sent.decode("ascii").split("\r\n")
    assert lines[0] == "GET /api_v2/germline/set/9606.IGH_VDJ/10/airr?q=1 HTTP/1.1"
    headers = {line.split(":", 1)[0].lower(): line.split(":", 1)[1].strip() for line in lines[1:] if line}
    assert headers["host"] == "ogrdb.airr-community.org" and headers["accept-encoding"] == "identity"
    assert not {"authorization", "cookie", "proxy-authorization"} & set(headers)
    assert response.status == 200 and response.getheader("content-type") == "text/plain"
    assert response.read(100) == b"payload" and response.read(100) == b""
    response.close()
    assert connection.simulated.closed


def test_the_https_transport_verifies_certificates_and_host_names(monkeypatch):
    clock = FakeClock()
    data, head = http_response(b"x")
    opened = simulate_connections(monkeypatch, clock, lambda: SimulatedSocket(clock, data, head))
    HttpsTransport(clock=clock).get(_url(IGH, "airr"), timeout=30, deadline=clock.now + 120).close()
    context = opened[0].given_context
    assert context.verify_mode == ssl.CERT_REQUIRED and context.check_hostname is True
    unverified = ssl.create_default_context()
    unverified.check_hostname = False
    with pytest.raises(ValueError):
        HttpsTransport(context=unverified)
    unverified.verify_mode = ssl.CERT_NONE
    with pytest.raises(ValueError):
        HttpsTransport(context=unverified)


def test_the_https_transport_returns_after_one_receive(monkeypatch):
    clock = FakeClock()
    data, head = http_response(b"0123456789")
    sockets = []
    simulate_connections(monkeypatch, clock, lambda: sockets.append(SimulatedSocket(clock, data, head, chunk=1))
                         or sockets[-1])
    response = HttpsTransport(clock=clock).get(_url(IGH, "airr"), timeout=30, deadline=clock.now + 120)
    before = sockets[0].recvs
    assert response.read(1 << 16) == b"0"
    assert sockets[0].recvs == before + 1


def _trickle(monkeypatch, clock, body, **steps):
    data, head = http_response(body)
    sockets = []
    simulate_connections(monkeypatch, clock,
                         lambda: sockets.append(SimulatedSocket(clock, data, head, chunk=1, **steps)) or sockets[-1])
    return sockets, head


def _bounded_by_the_deadline(sock, deadline):
    assert sock.timeouts
    for at, timeout in sock.timeouts:
        assert 0 < timeout <= SOCKET_TIMEOUT_SECONDS and at + timeout <= deadline + 1e-9


def test_a_trickled_body_hits_the_deadline_promptly(tmp_path, monkeypatch):
    clock = FakeClock()
    start = clock.now
    sockets, head = _trickle(monkeypatch, clock, _fixture(IGH, "airr"), body_step=29.0)
    error = _error("resource_limit", fetch, tmp_path, HttpsTransport(clock=clock), formats=("airr",), clock=clock)
    assert "deadline" in error.message and "120" in error.message
    assert clock.now - start == pytest.approx(ACQUISITION_DEADLINE_SECONDS)
    assert sockets[0].pos - head == 4                    # four bytes at 29 s; the fifth is cut at the deadline
    _bounded_by_the_deadline(sockets[0], start + ACQUISITION_DEADLINE_SECONDS)
    assert _snapshot(tmp_path / "cache" / "blobs" / "sha256") == {}
    assert _temporaries(tmp_path) == []


def test_a_trickled_status_line_and_headers_hit_the_deadline_promptly(tmp_path, monkeypatch):
    clock = FakeClock()
    start = clock.now
    sockets, head = _trickle(monkeypatch, clock, _fixture(IGH, "airr"), head_step=29.0)
    error = _error("resource_limit", fetch, tmp_path, HttpsTransport(clock=clock), formats=("airr",), clock=clock)
    assert "deadline" in error.message
    assert clock.now - start == pytest.approx(ACQUISITION_DEADLINE_SECONDS)
    assert sockets[0].pos == 4 < head
    _bounded_by_the_deadline(sockets[0], start + ACQUISITION_DEADLINE_SECONDS)
    assert _temporaries(tmp_path) == []


def test_a_receive_slower_than_the_socket_timeout_is_a_network_error(tmp_path, monkeypatch):
    clock = FakeClock()
    start = clock.now
    _trickle(monkeypatch, clock, _fixture(IGH, "airr"), body_step=40.0)
    _error("network", fetch, tmp_path, HttpsTransport(clock=clock), formats=("airr",), clock=clock)
    assert clock.now - start == pytest.approx(SOCKET_TIMEOUT_SECONDS)


def test_a_complete_download_through_the_https_transport(tmp_path, monkeypatch):
    clock = FakeClock()
    body = _fixture(IGH, "airr")
    data, head = http_response(body, headers=[("Content-Type", "application/octet-stream")])
    simulate_connections(monkeypatch, clock, lambda: SimulatedSocket(clock, data, head, body_step=0.001))
    result = fetch(tmp_path, HttpsTransport(clock=clock), formats=("airr",), clock=clock)
    assert result.assets["airr"] == body
    assert result.lock.asset("airr").sha256 == _sha(body)
    assert dict(result.receipts[0].headers) == {"content-length": str(len(body)),
                                                "content-type": "application/octet-stream"}


def test_a_slow_connect_is_bounded_by_the_remaining_deadline(monkeypatch):
    clock = FakeClock()
    data, head = http_response(b"x")
    opened = simulate_connections(monkeypatch, clock, lambda: SimulatedSocket(clock, data, head), connect_cost=20.0)
    start = clock.now
    with pytest.raises(DeadlineExceeded):
        HttpsTransport(clock=clock).get(_url(IGH, "airr"), timeout=30, deadline=start + 10)
    assert opened[0].timeout == pytest.approx(10) and clock.now - start == pytest.approx(10)
    assert issubclass(DeadlineExceeded, socket.timeout)


def test_a_passed_deadline_opens_no_connection(monkeypatch):
    clock = FakeClock()
    opened = simulate_connections(monkeypatch, clock, lambda: pytest.fail("a socket was created"))
    with pytest.raises(DeadlineExceeded):
        HttpsTransport(clock=clock).get(_url(IGH, "airr"), timeout=30, deadline=clock.now)
    assert opened == []


@pytest.mark.parametrize("url", ["http://ogrdb.airr-community.org/x", "ftp://ogrdb.airr-community.org/x",
                                 "https://user:pw@ogrdb.airr-community.org/x",
                                 "https://ogrdb.airr-community.org/\u00e9"])
def test_the_https_transport_refuses_anything_but_plain_https(monkeypatch, url):
    monkeypatch.setattr(http.client, "HTTPSConnection",
                        lambda *a, **k: pytest.fail("a connection was opened"))
    with pytest.raises(ValueError):
        HttpsTransport().get(url, timeout=SOCKET_TIMEOUT_SECONDS, deadline=float("inf"))


def test_the_source_model_is_not_changed_by_acquisition(tmp_path):
    lock = fetch(tmp_path, FakeTransport(serve())).lock
    assert m.SourceLock.from_dict(copy.deepcopy(lock.to_dict())) == lock
    for receipt in fetch(tmp_path, FakeTransport(serve()), cache="again").receipts:
        with pytest.raises(TypeError):
            m.lock_hash(receipt)
