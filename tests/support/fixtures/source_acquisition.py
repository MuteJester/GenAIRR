"""A fake HTTP transport serving the pinned export bytes.

Lets acquisition be exercised without a network: responses stream in bounded chunks and
can be made to fail partway.
"""
from __future__ import annotations

from tests.support import paths

from GenAIRR._artifact_io import strict_loads

DATA = paths.data("reference_sources")

SOURCES = strict_loads((DATA / "sources.json").read_bytes().decode("utf-8"), max_bytes=1 << 22)

SUFFIX = {"airr": "json", "airr_ex": "json", "gapped": "fasta", "gapped_ex": "fasta"}

IGH, IGK = ("9606.IGH_VDJ", 10), ("9606.IGKappa_VJ", 5)

ALL = ("airr", "airr_ex", "gapped", "gapped_ex")

def _fixture(source, fmt) -> bytes:
    set_id, release = source
    return (DATA / "ogrdb" / set_id / str(release) / f"{fmt}.{SUFFIX[fmt]}").read_bytes()

def _sources_asset(source, fmt) -> dict:
    (entry,) = [s for s in SOURCES["sets"] if (s["set_id"], s["release"]) == source]
    (asset,) = [a for a in entry["assets"] if a["format"] == fmt]
    return asset

def _url(source, fmt) -> str:
    return f"https://ogrdb.airr-community.org/api_v2/germline/set/{source[0]}/{source[1]}/{fmt}"

class FakeResponse:
    """One HTTP response: status, headers and a body streamed in bounded chunks."""

    def __init__(self, status=200, body=b"", *, headers=None, content_length=True, chunk=1 << 16,
                 fail_after=None, error=None, filler=None, on_read=None):
        self.status = status
        self._body = body
        self._filler = filler                   # (byte, length): a lazily produced body
        length = filler[1] if filler else len(body)
        self._headers = {k.lower(): v for k, v in (headers or {}).items()}
        if content_length is True:
            self._headers.setdefault("content-length", str(length))
        elif content_length is not False:
            self._headers["content-length"] = str(content_length)
        self._chunk = chunk
        self._fail_after = fail_after
        self._error = error
        self._on_read = on_read
        self.served = 0
        self.reads = 0
        self.closed = False

    def getheader(self, name, default=None):
        return self._headers.get(name.lower(), default)

    def read(self, amt):
        self.reads += 1
        if self._on_read:
            self._on_read()
        if self._fail_after is not None and self.served >= self._fail_after:
            raise self._error
        size = self._filler[1] if self._filler else len(self._body)
        n = min(amt, self._chunk, size - self.served)
        if self._fail_after is not None:
            n = min(n, self._fail_after - self.served)
        if self._filler:
            data = self._filler[0] * n
        else:
            data = self._body[self.served:self.served + n]
        self.served += n
        return data

    def close(self):
        self.closed = True

def ok(body, **kwargs):
    return lambda: FakeResponse(200, body, **kwargs)

class FakeTransport:
    def __init__(self, routes=None):
        self.routes = dict(routes or {})
        self.calls = []
        self.deadlines = []
        self.responses = []

    def get(self, url, *, timeout, deadline):
        self.calls.append((url, timeout))
        self.deadlines.append(deadline)
        route = self.routes.get(url)
        response = route() if route else FakeResponse(404, b"not found")
        self.responses.append(response)
        return response

    @property
    def urls(self):
        return [url for url, _ in self.calls]

def headers_for(source, fmt, **extra):
    asset = _sources_asset(source, fmt)
    return {"Content-Type": asset["content_type"], "Content-Disposition": asset["content_disposition"],
            "Set-Cookie": "session=secret", "Server": "nginx", **extra}

def serve(source=IGH, formats=ALL, *, at=None, **kwargs):
    """Routes serving the pinned fixture bytes of ``source`` (at the URLs of ``at``)."""
    at = at or source
    return {_url(at, fmt): ok(_fixture(source, fmt), headers=headers_for(source, fmt), **kwargs)
            for fmt in formats}
