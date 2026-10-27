"""Explicit, bounded acquisition of published germline-source bytes.

Three operations produce a sealed :class:`~.lock.SourceLock` together with the exact
verified bytes of every locked asset and one acquisition receipt per asset:

* :func:`fetch_source` downloads the requested OGRDB exports through the transport it is
  handed.  It is the only operation that can reach the network.
* :func:`ingest_files` records bytes supplied as local files.  It never claims that
  GenAIRR fetched them.
* :func:`read_locked` checks every blob a lock requires in the cache.  A missing blob is
  an actionable error, never a download trigger.

The last two take no transport at all.

Interpretation stays outside this module.  The caller supplies an ``identify`` callable
that reads the set identity (taxon, locus, chain type) from the verified bytes by format.
The caller also supplies the data terms, or they are recorded as unresolved.

Provider URLs are generated from validated components only.  Only HTTPS to the known
OGRDB origin is accepted, and at most three same-origin HTTPS redirects are followed.
Each response streams under a hard byte bound, checked against Content-Length when one
is present but never trusting it alone.  Every blocking step (connecting, sending, the
status line and headers, and each receive of the body) waits at most the socket timeout
and never beyond the whole-operation deadline.  Nothing is retried automatically.

Every asset of an operation is staged in an owned temporary file and validated: digest,
pins of an expected lock, and identity.  Only then is anything published, so a failure
never replaces a valid cache blob or lock.

Receipts are audit events.  They never enter the lock, so identical bytes always seal
the same lock.
"""
from __future__ import annotations

from dataclasses import dataclass
import datetime
import http.client
import io
import os
import re
import socket
import ssl
import stat
import time
from types import MappingProxyType
from typing import Callable, Mapping, Optional, Protocol, Union
from urllib.parse import urldefrag, urljoin, urlsplit

from GenAIRR._reference_sources._fields import (
    MAX_ASSET_BYTES, MAX_ASSETS, MAX_RELEASE, MAX_SOURCE_PAYLOAD_BYTES, MAX_TEXT, bad, check_text, check_url, ptr,
    show,
)
from GenAIRR._reference_sources.cache import BlobCache
from GenAIRR._reference_sources.common import (
    EVIDENCE_FORMATS, FORMAT_ROLES, FORMATS, RECORD_EXPORT_FORMATS, SUPPORTED_ADAPTERS, AdapterRef,
    CoordinateProfileRef, SourceIdentity,
)
from GenAIRR._reference_sources.errors import ReferenceSourceError
from GenAIRR._reference_sources.lock import _SET_ID_RE, Asset, MissingAsset, SourceLock, Terms
from GenAIRR._reference_sources.receipt import RECEIPT_HEADERS, AcquisitionReceipt

# ── Policy constants ─────────────────────────────────────────────────

OGRDB_HOST = "ogrdb.airr-community.org"
OGRDB_ORIGIN = "https://" + OGRDB_HOST
_API_PATH = "/api_v2/germline/set"
# Provider exports that can be downloaded.  Statement pages are never fetched.
FETCH_FORMATS = ("airr", "airr_ex", "gapped", "gapped_ex")
# Provider exports a local file may stand for (their endpoints share one URL pattern).
_EXPORT_FORMATS = ("airr", "airr_ex", "gapped", "gapped_ex", "ungapped", "ungapped_ex")
DATA_TERMS_FORMATS = tuple(f for f in FORMATS if FORMAT_ROLES[f] == "data_terms")
SOCKET_TIMEOUT_SECONDS = 30
ACQUISITION_DEADLINE_SECONDS = 120
MAX_FOLLOWED_REDIRECTS = 3
REDIRECT_STATUSES = (301, 302, 303, 307, 308)
_CHUNK = 1 << 16
_DIGITS_RE = re.compile(r"[0-9]+")
_FETCH_TERMS_REASON = ("no data-terms statement was supplied with this acquisition; GenAIRR does not fetch "
                       "provider statement pages")
_LOCAL_TERMS_REASON = "no data-terms statement was supplied with these files"


# ── Request components ──────────────────────────────────────────────


def _check_set_id(value, path: str) -> str:
    if type(value) is not str or not _SET_ID_RE.fullmatch(value):
        raise bad(path, f"expected a provider set id matching [A-Za-z0-9][A-Za-z0-9._-]* (at most 128 "
                        f"characters), got {show(value)}")
    return value


def _check_release(value, path: str) -> int:
    if type(value) is not int or not 1 <= value <= MAX_RELEASE:
        raise bad(path, f"expected an explicit positive integer release (1..{MAX_RELEASE}), got {show(value)}; "
                        "'latest' and other implicit releases are not supported")
    return value


def _check_fetch_format(value, path: str) -> str:
    if type(value) is not str:
        raise bad(path, f"expected a format name, got {show(value)}")
    if value not in FETCH_FORMATS:
        raise bad(path, f"format {value!r} cannot be fetched; fetchable formats: {list(FETCH_FORMATS)}",
                  "unsupported_profile")
    return value


def _endpoint(set_id: str, release: int, fmt: str) -> str:
    return f"{OGRDB_ORIGIN}{_API_PATH}/{set_id}/{release}/{fmt}"


def provider_url(set_id, release, fmt) -> str:
    """The OGRDB endpoint of one export, generated from validated components only."""
    return _endpoint(_check_set_id(set_id, "/set_id"), _check_release(release, "/release"),
                     _check_fetch_format(fmt, "/format"))


def resolve_profile(profile) -> AdapterRef:
    """The adapter reading an interpretation profile, from its exact id or an :class:`AdapterRef`."""
    if isinstance(profile, AdapterRef):
        return profile
    if type(profile) is not str:
        raise bad("", f"expected an interpretation profile id, got {show(profile)}")
    matches = [AdapterRef(adapter_id, version, CoordinateProfileRef(profile_id, profile_version))
               for (adapter_id, version), (_, profiles) in SUPPORTED_ADAPTERS.items()
               for profile_id, profile_version in profiles if profile_id == profile]
    if len(matches) != 1:
        known = sorted(p for _, profiles in SUPPORTED_ADAPTERS.values() for p, _ in profiles)
        raise bad("", f"unknown interpretation profile {profile!r}; supported: {known}", "unsupported_profile")
    return matches[0]


def _profile(value, path: str) -> AdapterRef:
    try:
        return resolve_profile(value)
    except ReferenceSourceError as exc:
        raise exc.under(path) from None


def _check_formats(formats, path: str) -> tuple:
    if isinstance(formats, (str, bytes)) or not isinstance(formats, (tuple, list)):
        raise bad(path, f"expected a list of formats, got {show(formats)}")
    if not formats:
        raise bad(path, "no formats were requested")
    for i, fmt in enumerate(formats):
        _check_fetch_format(fmt, ptr(path, i))
    if len(set(formats)) != len(formats):
        raise bad(path, f"duplicate formats in {list(formats)}")
    if not any(fmt in RECORD_EXPORT_FORMATS for fmt in formats):
        raise bad(path, f"at least one record export {list(RECORD_EXPORT_FORMATS)} must be requested")
    return tuple(fmt for fmt in FETCH_FORMATS if fmt in formats)


# ── Transport ───────────────────────────────────────────────────────


class DeadlineExceeded(socket.timeout):
    """The operation's acquisition deadline passed during a transport step.

    It is a socket timeout, so code unaware of it still sees a timeout; the acquisition
    reports it as ``resource_limit`` with the elapsed and allowed time.
    """


class TransportResponse(Protocol):
    """One HTTP response: ``status``, header lookup, bounded reads and ``close``."""

    status: int

    def getheader(self, name: str, default: Optional[str] = None) -> Optional[str]: ...

    def read(self, amt: int) -> bytes:
        """At most ``amt`` bytes, returning as soon as some arrive; ``b""`` at the end."""

    def close(self) -> None: ...


class Transport(Protocol):
    """One HTTP GET that never follows redirects.

    ``timeout`` bounds each blocking socket step.  ``deadline`` is the absolute time, on
    the monotonic clock of the operation, by which the whole operation must end: once it
    passes, every step (connecting, sending, the status line and headers, and each body
    read) raises :class:`DeadlineExceeded`.  Other failures raise ``OSError`` or
    ``http.client.HTTPException``.
    """

    def get(self, url: str, *, timeout: float, deadline: float) -> TransportResponse: ...


def _user_agent() -> str:
    try:
        from GenAIRR import __version__ as version
    except ImportError:  # pragma: no cover - partially initialized package
        version = "unknown"
    return f"GenAIRR/{version} (reference source import)"


class _StepBound:
    """Socket timeout of each blocking step: ``min(timeout, time left before the deadline)``."""

    def __init__(self, clock: Callable[[], float], timeout: float, deadline: float):
        self._clock, self._timeout, self._deadline = clock, timeout, deadline

    def next(self) -> tuple:
        """``(seconds, bounded_by_deadline)``; :class:`DeadlineExceeded` once the deadline has passed."""
        remaining = self._deadline - self._clock()
        if remaining <= 0:
            raise DeadlineExceeded("the acquisition deadline has passed")
        if remaining < self._timeout:
            return remaining, True
        return self._timeout, False

    def run(self, sock, step, *args):
        """Run one blocking socket step under the next bound."""
        seconds, by_deadline = self.next()
        sock.settimeout(seconds)
        try:
            return step(*args)
        except socket.timeout as exc:
            if by_deadline and not isinstance(exc, DeadlineExceeded):
                raise DeadlineExceeded("the acquisition deadline passed while waiting for the server") from exc
            raise


class _BoundedReader(io.RawIOBase):
    """Raw receives from a connected socket, each bounded by :class:`_StepBound`."""

    def __init__(self, sock, raw, bound: _StepBound):
        self._sock, self._raw, self._bound = sock, raw, bound

    def readable(self) -> bool:
        return True

    def readinto(self, buffer) -> int:
        return self._bound.run(self._sock, self._raw.readinto, buffer)

    def close(self) -> None:
        try:
            self._raw.close()
        finally:
            super().close()


class _BoundedSocket:
    """The connected socket as ``http.client`` uses it, with every send and receive bounded."""

    def __init__(self, sock, bound: _StepBound):
        self._sock, self._bound = sock, bound

    def makefile(self, mode="r", buffering=None, **kwargs):
        if mode != "rb":
            raise ValueError("only binary reads are supported")
        raw = self._sock.makefile("rb", buffering=0)
        return io.BufferedReader(_BoundedReader(self._sock, raw, self._bound))

    def sendall(self, data, *args):
        return self._bound.run(self._sock, self._sock.sendall, data, *args)

    def close(self) -> None:
        self._sock.close()

    def __getattr__(self, name):
        return getattr(self._sock, name)


class _HttpsResponse:
    def __init__(self, connection, response):
        self._connection = connection
        self._response = response
        self.status = response.status

    def getheader(self, name: str, default: Optional[str] = None) -> Optional[str]:
        return self._response.getheader(name, default)

    def read(self, amt: int) -> bytes:
        return self._response.read1(amt)    # at most one underlying receive

    def close(self) -> None:
        try:
            self._response.close()
        finally:
            self._connection.close()


class HttpsTransport:
    """Standard-library HTTPS GET.

    Redirects are never followed here; the caller applies the redirect policy.  The
    request carries no credentials and asks for the unencoded bytes.  Certificates and host
    names are always verified: a custom ``context`` that does not verify both is refused.
    Connecting (with the TLS handshake), sending, receiving the status line and headers,
    and every body receive each run under ``min(timeout, time left before deadline)`` on
    ``clock``; body reads return after at most one receive.
    """

    def __init__(self, *, context: Optional[ssl.SSLContext] = None,
                 clock: Callable[[], float] = time.monotonic):
        if context is not None:
            if not isinstance(context, ssl.SSLContext):
                raise TypeError(f"context must be an ssl.SSLContext, got {type(context).__name__}")
            if context.verify_mode != ssl.CERT_REQUIRED or not context.check_hostname:
                raise ValueError("the TLS context must verify certificates (CERT_REQUIRED) and host names")
        self._context = context
        self._clock = clock

    def get(self, url: str, *, timeout: float, deadline: float) -> TransportResponse:
        parts = urlsplit(url)
        if parts.scheme != "https" or not parts.hostname or parts.username is not None or parts.password is not None:
            raise ValueError("only https URLs without credentials can be requested")
        bound = _StepBound(self._clock, timeout, deadline)
        seconds, by_deadline = bound.next()
        context = self._context if self._context is not None else ssl.create_default_context()
        connection = http.client.HTTPSConnection(parts.hostname, parts.port or 443, timeout=seconds,
                                                 context=context)
        target = (parts.path or "/") + (f"?{parts.query}" if parts.query else "")
        try:
            try:
                connection.connect()
            except socket.timeout as exc:
                if by_deadline and not isinstance(exc, DeadlineExceeded):
                    raise DeadlineExceeded("the acquisition deadline passed while connecting") from exc
                raise
            connection.sock = _BoundedSocket(connection.sock, bound)
            connection.request("GET", target, headers={"User-Agent": _user_agent(), "Accept": "*/*",
                                                       "Accept-Encoding": "identity", "Connection": "close"})
            response = connection.getresponse()
        except BaseException:
            connection.close()
            raise
        return _HttpsResponse(connection, response)


_TRANSPORT_ERRORS = (OSError, http.client.HTTPException)


# ── Results and inputs ──────────────────────────────────────────────


@dataclass(frozen=True)
class SuppliedFile:
    """A local file standing for one asset, optionally declaring the provider URL it came from.

    ``declared_url`` records the caller's claim; GenAIRR never claims it fetched the file.
    """

    path: Union[str, "os.PathLike[str]"]
    declared_url: Optional[str] = None

    def __post_init__(self):
        if isinstance(self.path, bool) or not isinstance(self.path, (str, os.PathLike)) \
                or not isinstance(os.fspath(self.path), str):
            raise TypeError(f"path must be a text filesystem path, got {type(self.path).__name__}")
        check_url(self.declared_url, "/declared_url", nullable=True)


@dataclass(frozen=True)
class AcquiredSource:
    """A sealed lock, the exact verified bytes of each locked asset, and one receipt per asset.

    ``assets`` is a read-only mapping from asset id to bytes in lock order; ``receipts``
    follow the same order.  Receipts are audit events and are not part of the lock.
    """

    lock: SourceLock
    assets: Mapping[str, bytes]
    receipts: tuple


# ── Budgets and clocks ──────────────────────────────────────────────


class _Budget:
    """Per-asset and per-operation byte bounds (read at use, from the shared budget)."""

    def __init__(self):
        self.used = 0

    def allowance(self) -> tuple:
        remaining = max(MAX_SOURCE_PAYLOAD_BYTES - self.used, 0)
        if remaining < MAX_ASSET_BYTES:
            return remaining, (f"the aggregate limit of {MAX_SOURCE_PAYLOAD_BYTES} bytes per operation "
                               f"({self.used} bytes already acquired, {remaining} remaining)")
        return MAX_ASSET_BYTES, f"the per-asset limit of {MAX_ASSET_BYTES} bytes"

    def consume(self, size: int) -> None:
        self.used += size


class _Deadline:
    """The whole-operation deadline; ``at`` is its absolute time on ``clock``."""

    def __init__(self, clock: Callable[[], float]):
        self._clock = clock
        self._start = clock()
        self.at = self._start + ACQUISITION_DEADLINE_SECONDS

    def exceeded(self, activity: str) -> ReferenceSourceError:
        elapsed = self._clock() - self._start
        return bad("", f"the acquisition deadline was exceeded while {activity}: {elapsed:.1f} s elapsed, "
                       f"allowed {ACQUISITION_DEADLINE_SECONDS} s", "resource_limit")

    def check(self, activity: str) -> None:
        if self._clock() - self._start > ACQUISITION_DEADLINE_SECONDS:
            raise self.exceeded(activity)


def _timestamp(now: Optional[Callable[[], datetime.datetime]]) -> str:
    moment = now() if now is not None else datetime.datetime.now(datetime.timezone.utc)
    if not isinstance(moment, datetime.datetime) or moment.utcoffset() is None:
        raise TypeError("now() must return a timezone-aware datetime")
    return moment.astimezone(datetime.timezone.utc).isoformat(timespec="microseconds")


def _offline_receipt(origin: str, fmt: str, sha256: str, size: int, now) -> AcquisitionReceipt:
    return AcquisitionReceipt(event_time_utc=_timestamp(now), origin=origin, format=fmt, sha256=sha256,
                              bytes=size, requested_url=None, final_url=None, redirects=(), http_status=None,
                              headers=())


# ── Argument guards ─────────────────────────────────────────────────


def _require_callable(value, name: str) -> None:
    if not callable(value):
        raise TypeError(f"{name} must be callable, got {type(value).__name__}")


def _as_cache(cache) -> BlobCache:
    return cache if isinstance(cache, BlobCache) else BlobCache(cache)


def _check_files(files, path: str, allowed: tuple) -> dict:
    """``{format: SuppliedFile}`` in canonical format order."""
    if not isinstance(files, Mapping):
        raise bad(path, f"expected a mapping of format to file, got {type(files).__name__}")
    if len(files) > MAX_ASSETS:
        raise bad(path, f"{len(files)} files exceed the limit of {MAX_ASSETS} assets", "resource_limit")
    supplied = {}
    for fmt, value in files.items():
        if type(fmt) is not str or fmt not in allowed:
            raise bad(ptr(path, fmt), f"unknown or unsupported format {show(fmt)} here; expected one of "
                                      f"{list(allowed)}")
        supplied[fmt] = value if isinstance(value, SuppliedFile) else SuppliedFile(value)
    return {fmt: supplied[fmt] for fmt in FORMATS if fmt in supplied}


def _check_declared(fmt: str, supplied: SuppliedFile, set_id: str, release: int, path: str) -> None:
    url = supplied.declared_url
    if url is None:
        return
    if fmt in _EXPORT_FORMATS:
        expected = _endpoint(set_id, release, fmt)
        if url != expected:
            raise bad(ptr(path, "declared_url"), f"declared URL {url} is not the provider endpoint {expected} of "
                                                 "this set, release and format", "identity_mismatch")
        return
    parts = urlsplit(url)
    if (parts.hostname or "").lower() != OGRDB_HOST or parts.port not in (None, 443):
        raise bad(ptr(path, "declared_url"), f"declared URL {url} is not on the provider origin {OGRDB_ORIGIN}",
                  "identity_mismatch")


def _check_terms(terms, statement_formats: tuple, default_reason: str, path: str) -> Terms:
    if terms is None:
        if statement_formats:
            raise bad(path, f"data-terms statements {list(statement_formats)} were supplied without recorded "
                            "terms")
        return Terms(status="unresolved", assets=(), licence=None, citation=None, reason=default_reason)
    if not isinstance(terms, Terms):
        raise TypeError(f"terms must be a Terms, got {type(terms).__name__}")
    if tuple(terms.assets) != tuple(statement_formats):
        raise bad(path, f"terms reference the statements {list(terms.assets)}, but the supplied statements are "
                        f"{list(statement_formats)}")
    return terms


# ── Staging ─────────────────────────────────────────────────────────


def _stage_file(fmt: str, supplied: SuppliedFile, cache: BlobCache, budget: _Budget, path: str):
    """Stream a local file into an owned temporary blob; return ``(staged, bytes)``."""
    name = os.fspath(supplied.path)
    try:
        status = os.stat(name)
    except FileNotFoundError:
        raise bad(path, f"{fmt}: the supplied file {name} does not exist", "io") from None
    except OSError as exc:
        raise bad(path, f"{fmt}: cannot inspect the supplied file {name} ({exc})", "io") from exc
    if not stat.S_ISREG(status.st_mode):
        raise bad(path, f"{fmt}: the supplied path {name} is not a regular file", "io")
    allowed, limit = budget.allowance()
    if status.st_size > allowed:
        raise bad(path, f"{fmt}: the supplied file {name} holds {status.st_size} bytes, more than {limit}",
                  "resource_limit")
    chunks, total = [], 0
    try:
        with open(name, "rb") as handle, cache.writer() as writer:
            if not stat.S_ISREG(os.fstat(handle.fileno()).st_mode):
                raise bad(path, f"{fmt}: the supplied path {name} is not a regular file", "io")
            while True:
                chunk = handle.read(min(_CHUNK, allowed + 1 - total))
                if not chunk:
                    break
                total += len(chunk)
                if total > allowed:
                    raise bad(path, f"{fmt}: the supplied file {name} grew beyond {limit} while being read "
                                    f"({total} bytes read, allowed {allowed})", "resource_limit")
                writer.write(chunk)
                chunks.append(chunk)
            if total == 0:
                raise bad(path, f"{fmt}: the supplied file {name} is empty")
            staged = writer.finish()
    except OSError as exc:
        raise bad(path, f"{fmt}: cannot read the supplied file {name} ({exc})", "io") from exc
    budget.consume(total)
    return staged, b"".join(chunks)


def _network(fmt: str, message: str) -> ReferenceSourceError:
    return bad("", f"{fmt}: {message}", "network")


def _content_length(fmt: str, url: str, value: Optional[str]) -> Optional[int]:
    if value is None:
        return None
    if not _DIGITS_RE.fullmatch(value.strip()):
        raise _network(fmt, f"the response from {url} has an invalid Content-Length {value!r}")
    return int(value.strip())


def _redirect_target(fmt: str, current: str, location: Optional[str], status: int) -> str:
    """The validated absolute target of a redirect: HTTPS on the OGRDB origin, no credentials."""
    if location is None or not location.strip():
        raise _network(fmt, f"{current} answered redirect status {status} without a Location header")
    target = urldefrag(urljoin(current, location.strip())).url
    try:
        parts = urlsplit(target)
        port = parts.port
    except ValueError:
        raise _network(fmt, f"refused a redirect from {current} to an invalid URL") from None
    scheme = parts.scheme.lower()
    if scheme == "http":
        raise _network(fmt, f"refused a redirect downgrade from https to http ({current} redirected to an "
                            "http URL); only https is accepted")
    if scheme != "https":
        raise _network(fmt, f"refused a redirect from {current} to a non-https URL (scheme {scheme!r})")
    if parts.username is not None or parts.password is not None:
        raise _network(fmt, f"refused a redirect from {current} to a URL carrying credentials")
    host = (parts.hostname or "").lower()
    if host != OGRDB_HOST or port not in (None, 443):
        origin = f"https://{host}" + (f":{port}" if port not in (None, 443) else "")
        raise _network(fmt, f"refused a redirect from {current} to another origin {origin}; only "
                            f"{OGRDB_ORIGIN} is accepted")
    try:
        check_url(target, "")
    except ReferenceSourceError as exc:
        raise _network(fmt, f"refused a redirect from {current} to an invalid URL ({exc.message})") from None
    return target


def _receipt_headers(response: TransportResponse) -> tuple:
    """The allowlisted response headers whose values a receipt can record."""
    kept = []
    for name in RECEIPT_HEADERS:
        value = response.getheader(name)
        if value is None:
            continue
        try:
            check_text(value, "", max_len=MAX_TEXT, multiline=False)
        except ReferenceSourceError:
            continue
        kept.append((name, value))
    return tuple(sorted(kept))


def _stream_body(fmt: str, url: str, response: TransportResponse, cache: BlobCache, deadline: _Deadline,
                 budget: _Budget):
    """Stream one 200 response into an owned temporary blob; return ``(staged, bytes)``."""
    encoding = response.getheader("content-encoding")
    if encoding is not None and encoding.strip().lower() not in ("", "identity"):
        raise _network(fmt, f"the response from {url} is encoded ({encoding}); only the raw bytes can be locked")
    declared = _content_length(fmt, url, response.getheader("content-length"))
    allowed, limit = budget.allowance()
    if declared is not None and declared > allowed:
        raise bad("", f"{fmt}: {url} declares a Content-Length of {declared} bytes, more than {limit}",
                  "resource_limit")
    cap = allowed if declared is None else declared
    chunks, total = [], 0
    with cache.writer() as writer:
        while True:
            try:
                chunk = response.read(min(_CHUNK, cap + 1 - total))
            except DeadlineExceeded as exc:
                raise deadline.exceeded(f"downloading {fmt}") from exc
            except _TRANSPORT_ERRORS as exc:
                raise _network(fmt, f"the connection to {url} failed after {total} bytes ({exc!r})") from exc
            deadline.check(f"downloading {fmt}")
            if not chunk:
                break
            total += len(chunk)
            if total > cap:
                if declared is not None:
                    raise _network(fmt, f"the response from {url} is longer than its Content-Length of "
                                        f"{declared} bytes")
                raise bad("", f"{fmt}: the response from {url} exceeds {limit}: {total} bytes received, "
                              f"allowed {allowed}", "resource_limit")
            writer.write(chunk)
            chunks.append(chunk)
        if declared is not None and total < declared:
            raise _network(fmt, f"truncated response from {url}: received {total} of the {declared} bytes declared "
                                "by Content-Length")
        if total == 0:
            raise bad("", f"{fmt}: the response from {url} is empty")
        staged = writer.finish()
    budget.consume(total)
    return staged, b"".join(chunks)


def _download(fmt: str, url: str, *, transport: Transport, cache: BlobCache, deadline: _Deadline,
              budget: _Budget, now):
    """GET one export, applying the redirect policy; return ``(staged, bytes, receipt)``."""
    current, redirects = url, []
    while True:
        deadline.check(f"requesting {fmt}")
        try:
            response = transport.get(current, timeout=SOCKET_TIMEOUT_SECONDS, deadline=deadline.at)
        except DeadlineExceeded as exc:
            raise deadline.exceeded(f"requesting {fmt}") from exc
        except _TRANSPORT_ERRORS as exc:
            raise _network(fmt, f"the request for {current} failed ({exc!r}); it is not retried") from exc
        staged = None
        try:
            try:
                status = response.status
                if status in REDIRECT_STATUSES:
                    if len(redirects) >= MAX_FOLLOWED_REDIRECTS:
                        raise _network(fmt, f"refused to follow more than {MAX_FOLLOWED_REDIRECTS} redirects "
                                            f"(requested {url}, redirected again by {current})")
                    current = _redirect_target(fmt, current, response.getheader("location"), status)
                    redirects.append(current)
                    continue
                if status != 200:
                    raise _network(fmt, f"{current} answered HTTP status {status}; only 200 with the complete "
                                        "body is accepted")
                staged, data = _stream_body(fmt, current, response, cache, deadline, budget)
                headers = _receipt_headers(response)
            finally:
                response.close()
            receipt = AcquisitionReceipt(
                event_time_utc=_timestamp(now), origin="provider_download", format=fmt, sha256=staged.sha256,
                bytes=staged.bytes, requested_url=url, final_url=current, redirects=tuple(redirects),
                http_status=200, headers=headers)
        except BaseException:
            if staged is not None:
                staged.discard()
            raise
        return staged, data, receipt


# ── Shared steps ────────────────────────────────────────────────────


def _identify(identify, set_id: str, release: int, data: dict) -> SourceIdentity:
    assets = MappingProxyType({fmt: data[fmt] for fmt in FORMATS if fmt in data})
    identity = identify(set_id=set_id, release=release, assets=assets)
    if not isinstance(identity, SourceIdentity):
        raise TypeError(f"identify must return a SourceIdentity, got {type(identity).__name__}")
    return identity


def _publish_all(cache: BlobCache, staged: list) -> None:
    """Publish every staged blob, after checking every destination first.

    An existing identical blob is accepted after full verification and a conflicting file
    fails the whole operation before anything new is published.
    """
    pending, seen = [], set()
    for item in staged:
        if item.sha256 in seen or cache.contains(item.sha256, item.bytes):
            item.discard()
        else:
            pending.append(item)
        seen.add(item.sha256)
    for item in pending:
        cache.publish(item)


def _result(lock: SourceLock, data: dict, receipts: dict) -> AcquiredSource:
    return AcquiredSource(lock=lock, assets=MappingProxyType({a.asset_id: data[a.format] for a in lock.assets}),
                          receipts=tuple(receipts[a.format] for a in lock.assets))


def _check_expected(lock, set_id: str, release: int, adapter: AdapterRef, requested: tuple) -> None:
    """An expected lock must pin exactly this request before any HTTP."""
    if not isinstance(lock, SourceLock):
        raise TypeError(f"expected_lock must be a SourceLock, got {type(lock).__name__}")
    if lock.set_id != set_id:
        raise bad("/expected_lock/set_id", f"the expected lock pins set {lock.set_id!r}, not the requested "
                                           f"{set_id!r}", "identity_mismatch")
    if lock.release != release:
        raise bad("/expected_lock/release", f"the expected lock pins release {lock.release}, not the requested "
                                            f"release {release}", "identity_mismatch")
    if lock.adapter != adapter:
        raise bad("/expected_lock/adapter", f"the expected lock was read with {lock.adapter.to_dict()}, not the "
                                            f"requested {adapter.to_dict()}", "identity_mismatch")
    formats = [a.format for a in lock.assets]
    for fmt in requested:
        if fmt not in formats:
            raise bad("/expected_lock/assets", f"the expected lock does not pin the requested format {fmt!r}; it "
                                               f"pins {formats}", "identity_mismatch")
        i = formats.index(fmt)
        asset = lock.assets[i]
        url = _endpoint(set_id, release, fmt)
        if asset.origin_kind != "provider_download" or asset.declared_url != url:
            raise bad(ptr("/expected_lock/assets", i), f"the expected lock pins {fmt} as {asset.origin_kind} from "
                                                       f"{asset.declared_url}, so it cannot be reacquired from {url}",
                      "identity_mismatch")


# ── Operations ──────────────────────────────────────────────────────


def fetch_source(*, set_id, release, profile, cache, transport: Transport, identify,
                 formats=("airr", "airr_ex"), terms: Optional[Terms] = None, statements=None,
                 expected_lock: Optional[SourceLock] = None, clock: Callable[[], float] = time.monotonic,
                 now: Optional[Callable[[], datetime.datetime]] = None) -> AcquiredSource:
    """Download the requested exports of one explicit set release and seal a lock.

    ``identify(set_id=, release=, assets=)`` receives the verified bytes by format and
    returns the :class:`SourceIdentity` or raises :class:`ReferenceSourceError`.  Recorded
    ``terms`` need their statement files in ``statements`` (``{format: path or
    SuppliedFile}``); without terms the lock records them as unresolved.

    Without ``expected_lock`` the result is a new candidate lock and makes no historical
    claim.  With one, the request must match its pins before any HTTP, every requested
    asset is downloaded again even if cached, and the bytes and identity must equal the
    pins.  Changed bytes are reported as drift (``digest_mismatch``) and are not stored.
    Pinned assets that were not requested are read from the verified cache.
    """
    _require_callable(identify, "identify")
    if not callable(getattr(transport, "get", None)):
        raise TypeError(f"transport must provide get(url, *, timeout, deadline), got {type(transport).__name__}")
    _require_callable(clock, "clock")
    if expected_lock is not None and not isinstance(expected_lock, SourceLock):
        raise TypeError(f"expected_lock must be a SourceLock, got {type(expected_lock).__name__}")
    set_id = _check_set_id(set_id, "/set_id")
    release = _check_release(release, "/release")
    adapter = _profile(profile, "/profile")
    requested = _check_formats(formats, "/formats")
    supplied = _check_files({} if statements is None else statements, "/statements", DATA_TERMS_FORMATS)
    for fmt, item in supplied.items():
        _check_declared(fmt, item, set_id, release, ptr("/statements", fmt))
    if expected_lock is not None:
        if terms is not None and terms != expected_lock.terms:
            raise bad("/terms", "the terms of an expected lock come from that lock; the supplied terms differ")
        if supplied:
            raise bad("/statements", "the statements of an expected lock are read from the cache, not supplied")
        _check_expected(expected_lock, set_id, release, adapter, requested)
        terms = expected_lock.terms
    else:
        terms = _check_terms(terms, tuple(supplied), _FETCH_TERMS_REASON, "/terms")
    cache = _as_cache(cache)

    deadline, budget = _Deadline(clock), _Budget()
    staged, data, receipts, assets = [], {}, {}, {}
    try:
        if expected_lock is not None:
            blobs = {}
            for i, asset in enumerate(expected_lock.assets):
                if asset.format in requested:
                    continue
                if asset.sha256 not in blobs:
                    try:
                        blobs[asset.sha256] = cache.read(asset.sha256, asset.bytes, asset=asset.asset_id)
                    except ReferenceSourceError as exc:
                        raise exc.under(ptr("/expected_lock/assets", i)) from None
                    budget.consume(asset.bytes)
                data[asset.format] = blobs[asset.sha256]
                receipts[asset.format] = _offline_receipt("verified_cache", asset.format, asset.sha256,
                                                          asset.bytes, now)
                assets[asset.format] = asset
        for fmt, item in supplied.items():
            blob, data[fmt] = _stage_file(fmt, item, cache, budget, ptr("/statements", fmt))
            staged.append(blob)
            receipts[fmt] = _offline_receipt("local_file", fmt, blob.sha256, blob.bytes, now)
            assets[fmt] = Asset(asset_id=fmt, role=FORMAT_ROLES[fmt], format=fmt, sha256=blob.sha256,
                                bytes=blob.bytes, declared_url=item.declared_url, origin_kind="local_file")
        for fmt in requested:
            url = _endpoint(set_id, release, fmt)
            blob, data[fmt], receipts[fmt] = _download(fmt, url, transport=transport, cache=cache,
                                                       deadline=deadline, budget=budget, now=now)
            staged.append(blob)
            if expected_lock is not None:
                pinned = expected_lock.asset(fmt)
                if (blob.sha256, blob.bytes) != (pinned.sha256, pinned.bytes):
                    i = expected_lock.asset_ids.index(fmt)
                    raise bad(ptr("/expected_lock/assets", i),
                              f"drift: {fmt} of {set_id} release {release} changed upstream: {url} served "
                              f"{blob.sha256} ({blob.bytes} bytes), but the expected lock pins {pinned.sha256} "
                              f"({pinned.bytes} bytes); the changed bytes were not stored and the lock and cache "
                              "are unchanged", "digest_mismatch")
            assets[fmt] = Asset(asset_id=fmt, role=FORMAT_ROLES[fmt], format=fmt, sha256=blob.sha256,
                                bytes=blob.bytes, declared_url=url, origin_kind="provider_download")
        identity = _identify(identify, set_id, release, data)
        if expected_lock is not None and identity != expected_lock.identity:
            raise bad("/expected_lock/identity", f"the bytes identify as {identity.to_dict()}, but the expected lock "
                                                 f"pins {expected_lock.identity.to_dict()}", "identity_mismatch")
        ordered = tuple(assets[fmt] for fmt in FORMATS if fmt in assets)
        if expected_lock is not None:
            missing = expected_lock.missing_assets
        else:
            missing = tuple(MissingAsset(format=fmt, role=FORMAT_ROLES[fmt], reason="not_requested", detail=None)
                            for fmt in EVIDENCE_FORMATS if fmt not in assets)
        lock = SourceLock.seal(provider=adapter.provider, set_id=set_id, release=release, identity=identity,
                               adapter=adapter, assets=ordered, missing_assets=missing, terms=terms)
        if expected_lock is not None and lock != expected_lock:
            raise bad("/expected_lock", f"the reacquired lock {lock.commitment} differs from the expected lock "
                                        f"{expected_lock.commitment}", "identity_mismatch")
        _publish_all(cache, staged)
    finally:
        for item in staged:
            item.discard()
    return _result(lock, data, receipts)


def ingest_files(*, set_id, release, profile, files, cache, identify, terms: Optional[Terms] = None,
                 now: Optional[Callable[[], datetime.datetime]] = None) -> AcquiredSource:
    """Record supplied local files (``{format: path or SuppliedFile}``) and seal a local lock.

    Never touches the network.  Assets have origin ``local_file`` and a null or
    caller-declared provider URL; receipts have origin ``local_file``.  Lock identity
    depends on the bytes only, never on paths or file names.
    """
    _require_callable(identify, "identify")
    set_id = _check_set_id(set_id, "/set_id")
    release = _check_release(release, "/release")
    adapter = _profile(profile, "/profile")
    supplied = _check_files(files, "/files", FORMATS)
    if not supplied:
        raise bad("/files", "no files were supplied")
    if not any(fmt in RECORD_EXPORT_FORMATS for fmt in supplied):
        raise bad("/files", f"at least one record export {list(RECORD_EXPORT_FORMATS)} must be supplied")
    for fmt, item in supplied.items():
        _check_declared(fmt, item, set_id, release, ptr("/files", fmt))
    statements = tuple(fmt for fmt in supplied if FORMAT_ROLES[fmt] == "data_terms")
    terms = _check_terms(terms, statements, _LOCAL_TERMS_REASON, "/terms")
    cache = _as_cache(cache)

    budget = _Budget()
    staged, data, receipts, assets = [], {}, {}, []
    try:
        for fmt, item in supplied.items():
            blob, data[fmt] = _stage_file(fmt, item, cache, budget, ptr("/files", fmt))
            staged.append(blob)
            receipts[fmt] = _offline_receipt("local_file", fmt, blob.sha256, blob.bytes, now)
            assets.append(Asset(asset_id=fmt, role=FORMAT_ROLES[fmt], format=fmt, sha256=blob.sha256,
                                bytes=blob.bytes, declared_url=item.declared_url, origin_kind="local_file"))
        identity = _identify(identify, set_id, release, data)
        missing = tuple(MissingAsset(format=fmt, role=FORMAT_ROLES[fmt], reason="not_supplied", detail=None)
                        for fmt in EVIDENCE_FORMATS if fmt not in supplied)
        lock = SourceLock.seal(provider=adapter.provider, set_id=set_id, release=release, identity=identity,
                               adapter=adapter, assets=tuple(assets), missing_assets=missing, terms=terms)
        _publish_all(cache, staged)
    finally:
        for item in staged:
            item.discard()
    return _result(lock, data, receipts)


def read_locked(lock, *, cache, identify=None,
                now: Optional[Callable[[], datetime.datetime]] = None) -> AcquiredSource:
    """Read a lock (a :class:`SourceLock` or a lock file path) and verify every blob it needs.

    Never touches the network.  Blobs missing from the cache are ``missing_asset``,
    naming every such asset and its expected digest.  A corrupt blob is
    ``digest_mismatch``.  With ``identify``, the identity read from the bytes must equal
    the lock's.  Receipts have origin ``verified_cache``: a cached read is not a fresh
    acquisition.
    """
    if isinstance(lock, (str, os.PathLike)) and not isinstance(lock, bool):
        lock = SourceLock.read(lock)
    elif not isinstance(lock, SourceLock):
        raise TypeError(f"lock must be a SourceLock or a lock file path, got {type(lock).__name__}")
    if identify is not None:
        _require_callable(identify, "identify")
    cache = _as_cache(cache)
    data, blobs, missing = {}, {}, []
    for i, asset in enumerate(lock.assets):
        if asset.sha256 not in blobs:
            try:
                blobs[asset.sha256] = cache.read(asset.sha256, asset.bytes, asset=asset.asset_id)
            except ReferenceSourceError as exc:
                if exc.kind != "missing_asset":
                    raise exc.under(ptr("/assets", i)) from None
                missing.append((i, asset))
                continue
        data[asset.format] = blobs[asset.sha256]
    if missing:
        names = "; ".join(f"asset {a.asset_id} ({a.sha256}, {a.bytes} bytes)" for _, a in missing)
        path = ptr("/assets", missing[0][0]) if len(missing) == 1 else "/assets"
        raise bad(path, f"not in the cache {cache.root}: {names}; reading never downloads: fetch the source "
                        "again with this lock as the expected lock, or import the files again", "missing_asset")
    if identify is not None:
        identity = _identify(identify, lock.set_id, lock.release, data)
        if identity != lock.identity:
            raise bad("/identity", f"the cached bytes identify as {identity.to_dict()}, but the lock pins "
                                   f"{lock.identity.to_dict()}", "identity_mismatch")
    receipts = {a.format: _offline_receipt("verified_cache", a.format, a.sha256, a.bytes, now) for a in lock.assets}
    return _result(lock, data, receipts)
