"""Acquisition receipt ``genairr_reference_source_receipt`` version 1.

A receipt is an audit event: when bytes were obtained, from where (the provider, a
supplied local file, or an already verified cache), the URLs and actual redirects, the
HTTP status and a few selected response headers, and the digest and size of the bytes.
Receipts carry no credentials, are not biological evidence and enter no semantic digest:
two downloads of identical bytes may have different receipts and the same source lock.
"""
from __future__ import annotations

from dataclasses import dataclass
import datetime
import re
from typing import Optional

from GenAIRR._reference_sources._fields import (
    MAX_ASSET_BYTES, bad, check_digest, check_enum, check_int, check_text, check_tuple, check_url, construct, parse_array,
    parse_json, ptr, read_json, require_object, check_header, show,
)
from GenAIRR._reference_sources.common import FORMATS

ARTIFACT = "genairr_reference_source_receipt"
SCHEMA_VERSION = 1
MAX_RECEIPT_BYTES = 64 * 1024
MAX_REDIRECTS = 16
RECEIPT_ORIGINS = ("provider_download", "local_file", "verified_cache")
# Response headers a receipt may keep; anything else (cookies, authorization) is refused.
RECEIPT_HEADERS = ("content-type", "content-length", "content-disposition", "etag", "last-modified", "date")

_TIME_RE = re.compile(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(\.\d{1,6})?(Z|\+00:00)")
_KEYS = ("artifact", "schema_version", "event_time_utc", "origin", "format", "sha256", "bytes", "requested_url",
         "final_url", "redirects", "http_status", "headers")


def check_utc_time(value, path: str) -> str:
    if type(value) is str and _TIME_RE.fullmatch(value):
        try:
            datetime.datetime.fromisoformat(value)
            return value
        except ValueError:
            pass
    raise bad(path, f"expected an ISO 8601 UTC time (YYYY-MM-DDTHH:MM:SS[.ffffff] with Z or +00:00), "
                    f"got {show(value)}")


@dataclass(frozen=True)
class AcquisitionReceipt:
    """One acquisition event; ``headers`` is a name-sorted tuple of ``(name, value)`` pairs."""

    event_time_utc: str
    origin: str
    format: str
    sha256: str
    bytes: int
    requested_url: Optional[str]
    final_url: Optional[str]
    redirects: tuple
    http_status: Optional[int]
    headers: tuple

    def __post_init__(self):
        check_utc_time(self.event_time_utc, "/event_time_utc")
        check_enum(self.origin, RECEIPT_ORIGINS, "/origin")
        check_enum(self.format, FORMATS, "/format")
        check_digest(self.sha256, "/sha256")
        check_int(self.bytes, "/bytes", minimum=1, maximum=MAX_ASSET_BYTES, limit_kind="resource_limit")
        check_tuple(self.redirects, "/redirects", MAX_REDIRECTS, str)
        check_tuple(self.headers, "/headers", len(RECEIPT_HEADERS), tuple)
        names = []
        for i, pair in enumerate(self.headers):
            if len(pair) != 2:
                raise bad(ptr("/headers", i), "expected a (name, value) pair")
            name, value = pair
            if name not in RECEIPT_HEADERS:
                raise bad("/headers", f"header {show(name)} is not a recorded header; receipts keep only "
                                      f"{list(RECEIPT_HEADERS)} (lower case) and never credentials")
            check_text(value, ptr("/headers", name), multiline=False)
            names.append(name)
        if names != sorted(set(names)):
            raise bad("/headers", "header names must be unique and sorted")
        if self.origin == "provider_download":
            check_url(self.requested_url, "/requested_url")
            check_url(self.final_url, "/final_url")
            for i, url in enumerate(self.redirects):
                check_url(url, ptr("/redirects", i))
            if self.redirects and self.redirects[-1] != self.final_url:
                raise bad("/redirects", "the last redirect must be the final URL")
            if not self.redirects and self.final_url != self.requested_url:
                raise bad("/redirects", "a final URL different from the requested URL needs its redirects")
            check_int(self.http_status, "/http_status", minimum=200, maximum=299)
        else:
            for key in ("requested_url", "final_url", "http_status"):
                if getattr(self, key) is not None:
                    raise bad("/" + key, f"a {self.origin} receipt has no {key}")
            if self.redirects:
                raise bad("/redirects", f"a {self.origin} receipt has no redirects")
            if self.headers:
                raise bad("/headers", f"a {self.origin} receipt has no response headers")

    def to_dict(self) -> dict:
        return {"artifact": ARTIFACT, "schema_version": SCHEMA_VERSION, "event_time_utc": self.event_time_utc,
                "origin": self.origin, "format": self.format, "sha256": self.sha256, "bytes": self.bytes,
                "requested_url": self.requested_url, "final_url": self.final_url,
                "redirects": list(self.redirects), "http_status": self.http_status,
                "headers": dict(self.headers)}

    @classmethod
    def from_dict(cls, d, path: str = "") -> "AcquisitionReceipt":
        check_header(d, ARTIFACT, SCHEMA_VERSION, path)
        require_object(d, _KEYS, path)
        headers = d["headers"]
        if not isinstance(headers, dict):
            raise bad(ptr(path, "headers"), "expected an object of header names to values")
        return construct(
            cls, path, event_time_utc=d["event_time_utc"], origin=d["origin"], format=d["format"],
            sha256=d["sha256"], bytes=d["bytes"], requested_url=d["requested_url"], final_url=d["final_url"],
            redirects=parse_array(d["redirects"], ptr(path, "redirects"), MAX_REDIRECTS, lambda v, p: v),
            http_status=d["http_status"], headers=tuple(sorted(headers.items(), key=lambda kv: str(kv[0]))))

    @classmethod
    def from_json(cls, text: str) -> "AcquisitionReceipt":
        return cls.from_dict(parse_json(text, max_bytes=MAX_RECEIPT_BYTES))

    @classmethod
    def read(cls, path) -> "AcquisitionReceipt":
        return cls.from_dict(read_json(path, max_bytes=MAX_RECEIPT_BYTES))
