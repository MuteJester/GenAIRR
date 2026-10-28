"""Strict field validators, size limits and JSON helpers for the published-source model.

Every validator raises :class:`ReferenceSourceError` with a JSON-pointer-like path
relative to the object being validated; ``construct`` and the ``from_dict`` readers
prefix that path with the object's location in its document.  Exceeding a named size
limit is a ``resource_limit`` error; every other violation is ``malformed`` unless a
caller passes a more specific kind.
"""
from __future__ import annotations

from contextlib import contextmanager
import datetime
import json
import math
from pathlib import Path
import re
import tempfile
import unicodedata
from typing import Any, Callable, Iterable, Iterator, Optional
from urllib.parse import urlsplit

from GenAIRR._artifact_io import (
    ArtifactIOError, canonical_json, read_text_bounded, strict_loads, typed_digest, typed_digest_of_canonical,
)
from GenAIRR._reference.model import check_label as _runtime_label_problem
from GenAIRR._reference_sources.errors import ReferenceSourceError

# ── Resource budget shared by every format (acquisition and cache reuse these) ──

MAX_ASSET_BYTES = 16 * 1024 * 1024             # one source asset
MAX_SOURCE_PAYLOAD_BYTES = 128 * 1024 * 1024   # all source assets of one lock together
MAX_ASSETS = 16                                # assets per lock
MAX_SOURCE_RECORDS = 100_000                   # source record occurrences per inventory
MAX_JSON_DEPTH = 64                            # nesting of owned JSON documents and raw provider payloads
# The deepest place a preserved provider value sits in any owned document: a build report's
# embedded inventory -> observations -> one observation -> its value.  Claim values and
# finding observations sit at the same depth.
VALUE_EMBEDDING_PATH = ("inventory", "observations", "<index>", "value")
VALUE_EMBEDDING_DEPTH = len(VALUE_EMBEDDING_PATH)
# Nesting allowed inside a preserved value, so that every owned document embedding it stays
# within MAX_JSON_DEPTH and is accepted by the strict parser and the typed digest alike.
MAX_EMBEDDED_VALUE_DEPTH = MAX_JSON_DEPTH - VALUE_EMBEDDING_DEPTH

# ── Size limits ─────────────────────────────────────────────────────

MAX_TEXT = 4096                 # free text: reasons, details, problems, citations, hints
MAX_NAME = 256                  # labels, reference names/versions, scientific names
MAX_URL = 2048                  # declared and receipt URLs
MAX_POINTER = 1024              # JSON pointers inside a source asset
MAX_SEQUENCE = 65535            # one normalized coding sequence
MAX_VERSION = 1_000_000         # adapter, coordinate-profile and projection versions
MAX_RELEASE = 2 ** 31 - 1       # provider release numbers
MAX_TAXONOMY_ID = 2 ** 31 - 1
MAX_REFS = 4096                 # locator or identifier lists inside one item
MAX_VALUE_NODES = 65536         # JSON values inside one preserved provider value
MAX_VALUE_STRING = 1 << 20      # characters of one string inside a preserved provider value

_DIGEST_RE = re.compile(r"sha256:[0-9a-f]{64}")
_DATE_RE = re.compile(r"\d{4}-\d{2}-\d{2}")
_INDEX_RE = re.compile(r"0|[1-9][0-9]*")
_SHOW = 80


def bad(path: str, message: str, kind: str = "malformed") -> ReferenceSourceError:
    return ReferenceSourceError(kind, path, message)


def show(value: Any) -> str:
    """A short ``repr`` for error messages."""
    text = repr(value)
    return text if len(text) <= _SHOW else text[:_SHOW - 3] + "..."


@contextmanager
def private_directory(prefix: str, purpose: str) -> Iterator[Path]:
    """A new private temporary directory for ``purpose``, removed afterwards.

    Failing to create it is a typed ``io`` error.  Removing it is best effort: a failure (such as
    a Windows sharing violation) never replaces an error being reported and never fails an
    operation that succeeded.
    """
    try:
        holder = tempfile.TemporaryDirectory(prefix=prefix)
    except OSError as exc:
        raise bad("", f"cannot create a private temporary directory for {purpose}: {exc}", "io") from exc
    try:
        yield Path(holder.name)
    finally:
        try:
            holder.cleanup()
        except OSError:
            pass


def short_digest(value: Any) -> str:
    """A digest shortened for display: its algorithm and first 12 hex digits."""
    text = str(value)
    return text[:19] + "..." if len(text) > 19 else text


def ptr(path: str, token: Any) -> str:
    """Append one JSON-pointer token (escaping ``~`` and ``/``)."""
    return path + "/" + str(token).replace("~", "~0").replace("/", "~1")


def construct(cls, path: str, **kwargs):
    """Build ``cls(**kwargs)``, locating any validation error below ``path``."""
    try:
        return cls(**kwargs)
    except ReferenceSourceError as exc:
        raise exc.under(path) from None


def closed(cls, d: Any, path: str, **nested):
    """Generic ``from_dict`` for a flat shape: exactly ``cls.KEYS``; ``nested`` parses listed keys."""
    require_object(d, cls.KEYS, path)
    return construct(cls, path, **{k: (nested[k](d[k], ptr(path, k)) if k in nested else d[k]) for k in cls.KEYS})


def require_type(value: Any, cls: type, name: str) -> None:
    """Programming-error guard for function arguments (``TypeError``, not a document error)."""
    if not isinstance(value, cls):
        raise TypeError(f"{name} must be a {cls.__name__}, got {type(value).__name__}")


# ── Containers ──────────────────────────────────────────────────────


def require_object(value: Any, keys: Iterable[str], path: str) -> dict:
    """A JSON object with exactly ``keys`` (missing and unknown keys are both errors)."""
    if not isinstance(value, dict):
        raise bad(path, "expected an object")
    keys = tuple(keys)
    missing = [k for k in keys if k not in value]
    if missing:
        raise bad(path, f"missing keys {missing}")
    unknown = [k for k in value if k not in keys]
    if unknown:
        raise bad(path, f"unknown keys {[show(k) for k in unknown]}")
    return value


def require_array(value: Any, path: str, max_items: int) -> list:
    if not isinstance(value, list):
        raise bad(path, "expected an array")
    if len(value) > max_items:
        raise bad(path, f"{len(value)} items exceed the limit of {max_items}", "resource_limit")
    return value


def parse_array(value: Any, path: str, max_items: int, parse: Callable[[Any, str], Any]) -> tuple:
    return tuple(parse(item, ptr(path, i)) for i, item in enumerate(require_array(value, path, max_items)))


def check_tuple(value: Any, path: str, max_items: int, item_type: Optional[type] = None) -> tuple:
    if not isinstance(value, tuple):
        raise bad(path, "expected a tuple")
    if len(value) > max_items:
        raise bad(path, f"{len(value)} items exceed the limit of {max_items}", "resource_limit")
    if item_type is not None:
        for i, item in enumerate(value):
            if not isinstance(item, item_type):
                raise bad(ptr(path, i), f"expected {item_type.__name__}, got {type(item).__name__}")
    return value


def check_instance(value: Any, cls: type, path: str, *, nullable: bool = False):
    if value is None and nullable:
        return None
    if not isinstance(value, cls):
        raise bad(path, f"expected {cls.__name__}" + (" or null" if nullable else "")
                  + f", got {type(value).__name__}")
    return value


def check_canonical_order(keys: list, path: str, what: str) -> None:
    """``keys`` strictly increase: canonical order without duplicates."""
    for i in range(1, len(keys)):
        if keys[i] == keys[i - 1]:
            raise bad(path, f"duplicate {what} at index {i}")
        if keys[i] < keys[i - 1]:
            raise bad(path, f"{what} at index {i} is out of canonical order")


def check_unique(keys: list, path: str, what: str, kind: str = "malformed") -> None:
    seen = {}
    for i, key in enumerate(keys):
        if key in seen:
            raise bad(path, f"duplicate {what} at indexes {seen[key]} and {i}", kind)
        seen[key] = i


# ── Scalars ─────────────────────────────────────────────────────────


def _utf8(value: str, path: str) -> None:
    try:
        value.encode("utf-8")
    except UnicodeEncodeError:
        raise bad(path, "text is not valid UTF-8 (lone surrogate)") from None


def check_text(value: Any, path: str, *, nullable: bool = False, max_len: int = MAX_TEXT,
               multiline: bool = True) -> Optional[str]:
    """Non-blank text without control characters (newline and tab allowed when ``multiline``)."""
    if value is None and nullable:
        return None
    if type(value) is not str:
        raise bad(path, "expected a string" + (" or null" if nullable else ""))
    if len(value) > max_len:
        raise bad(path, f"{len(value)} characters exceed the limit of {max_len}", "resource_limit")
    if not value.strip():
        raise bad(path, "must not be blank")
    _utf8(value, path)
    for position, char in enumerate(value):
        if multiline and char in "\n\t":
            continue
        if unicodedata.category(char) == "Cc":
            raise bad(path, f"contains control character U+{ord(char):04X} at position {position}")
    return value


def check_label(value: Any, path: str) -> str:
    """An executable label: 1-256 characters, no whitespace, commas, control or format characters."""
    problem = _runtime_label_problem(value)
    if problem:
        raise bad(path, f"label {show(value)} {problem}")
    _utf8(value, path)
    return value


def check_token(value: Any, path: str, pattern: "re.Pattern", description: str, *,
                nullable: bool = False) -> Optional[str]:
    if value is None and nullable:
        return None
    if type(value) is not str or not pattern.fullmatch(value):
        raise bad(path, f"expected {description}, got {show(value)}")
    return value


def check_int(value: Any, path: str, *, minimum: Optional[int] = None, maximum: Optional[int] = None,
              nullable: bool = False, limit_kind: str = "malformed") -> Optional[int]:
    """An exact integer (``true`` and ``1.0`` are rejected)."""
    if value is None and nullable:
        return None
    if type(value) is not int:
        raise bad(path, "expected an integer" + (" or null" if nullable else "") + f", got {show(value)}")
    if minimum is not None and value < minimum:
        raise bad(path, f"must be >= {minimum}, got {value}")
    if maximum is not None and value > maximum:
        raise bad(path, f"must be <= {maximum}, got {value}", limit_kind)
    return value


def check_bool(value: Any, path: str) -> bool:
    if type(value) is not bool:
        raise bad(path, f"expected a Boolean, got {show(value)}")
    return value


def check_enum(value: Any, allowed: Iterable[str], path: str, *, nullable: bool = False) -> Optional[str]:
    if value is None and nullable:
        return None
    allowed = tuple(allowed)
    if type(value) is not str or value not in allowed:
        raise bad(path, f"must be one of {list(allowed)}" + (" or null" if nullable else "")
                  + f", got {show(value)}")
    return value


def check_digest(value: Any, path: str, *, nullable: bool = False) -> Optional[str]:
    if value is None and nullable:
        return None
    if type(value) is not str or not _DIGEST_RE.fullmatch(value):
        raise bad(path, f"expected sha256:<64 lowercase hex>, got {show(value)}")
    return value


def check_date(value: Any, path: str, *, nullable: bool = False) -> Optional[str]:
    """An ISO calendar date ``YYYY-MM-DD``."""
    if value is None and nullable:
        return None
    if type(value) is str and _DATE_RE.fullmatch(value):
        try:
            datetime.date.fromisoformat(value)
            return value
        except ValueError:
            pass
    raise bad(path, "expected an ISO date YYYY-MM-DD" + (" or null" if nullable else "") + f", got {show(value)}")


def check_url(value: Any, path: str, *, nullable: bool = False) -> Optional[str]:
    """An ``https`` URL without credentials, whitespace or control characters."""
    if value is None and nullable:
        return None
    if type(value) is not str:
        raise bad(path, "expected a URL string" + (" or null" if nullable else ""))
    if len(value) > MAX_URL:
        raise bad(path, f"URL of {len(value)} characters exceeds the limit of {MAX_URL}", "resource_limit")
    _utf8(value, path)
    if any(char.isspace() or unicodedata.category(char) == "Cc" for char in value):
        raise bad(path, "URL contains whitespace or control characters")
    try:
        parts = urlsplit(value)
    except ValueError as exc:
        raise bad(path, f"not a valid URL: {exc}") from None
    if parts.scheme != "https" or not parts.netloc:
        raise bad(path, f"expected an https URL, got {show(value)}")
    if "@" in parts.netloc:
        raise bad(path, "URL must not carry credentials")
    return value


def check_pointer(value: Any, path: str) -> str:
    """A non-root RFC 6901 JSON pointer into a source asset."""
    if type(value) is not str:
        raise bad(path, "expected a JSON pointer string")
    if len(value) > MAX_POINTER:
        raise bad(path, f"pointer of {len(value)} characters exceeds the limit of {MAX_POINTER}",
                  "resource_limit")
    if not value.startswith("/"):
        raise bad(path, f"expected a JSON pointer starting with '/', got {show(value)}")
    _utf8(value, path)
    for position, char in enumerate(value):
        if char == "~" and value[position + 1:position + 2] not in ("0", "1"):
            raise bad(path, f"invalid '~' escape at position {position}")
        if unicodedata.category(char) == "Cc":
            raise bad(path, f"contains control character U+{ord(char):04X} at position {position}")
    return value


def pointer_tokens(pointer: str) -> tuple:
    return tuple(t.replace("~1", "/").replace("~0", "~") for t in pointer.split("/")[1:])


def pointer_sort_key(pointer: str) -> tuple:
    """Document order: array-index tokens compare as integers, other tokens as text."""
    return tuple((0, int(t), "") if _INDEX_RE.fullmatch(t) else (1, 0, t) for t in pointer_tokens(pointer))


def check_member_path(value: Any, path: str, max_len: int) -> str:
    """A relative, forward-slash path: no absolute, drive, backslash, empty, '.' or '..' forms."""
    if type(value) is not str:
        raise bad(path, f"expected a name string, got {show(value)}")
    if not value:
        raise bad(path, "name must not be empty")
    if len(value) > max_len:
        raise bad(path, f"name of {len(value)} characters exceeds the limit of {max_len}", "resource_limit")
    _utf8(value, path)
    if "\\" in value:
        raise bad(path, f"name {show(value)} contains a backslash")
    if value.startswith("/"):
        raise bad(path, f"name {show(value)} is absolute")
    if re.match(r"[A-Za-z]:", value):
        raise bad(path, f"name {show(value)} has a drive form")
    if any(unicodedata.category(char) == "Cc" for char in value):
        raise bad(path, f"name {show(value)} contains a control character")
    for part in value.split("/"):
        if part in ("", ".", ".."):
            raise bad(path, f"name {show(value)} has an empty, '.' or '..' component")
    return value


# ── Preserved JSON values ───────────────────────────────────────────


def _check_json(value: Any, path: str) -> None:
    nodes = 0
    stack = [(value, path, 0)]
    while stack:
        item, where, depth = stack.pop()
        nodes += 1
        if nodes > MAX_VALUE_NODES:
            raise bad(path, f"value has more than {MAX_VALUE_NODES} JSON values", "resource_limit")
        if depth > MAX_EMBEDDED_VALUE_DEPTH:
            raise bad(where, f"value nests {depth} levels deep; the limit for a preserved provider value is "
                             f"{MAX_EMBEDDED_VALUE_DEPTH} ({MAX_JSON_DEPTH} for a whole document, less "
                             f"{VALUE_EMBEDDING_DEPTH} levels of embedding)", "resource_limit")
        kind = type(item)
        if item is None or kind is bool:
            continue
        if kind is int:
            if not -2 ** 63 <= item < 2 ** 64:
                raise bad(where, "integer outside [-2**63, 2**64)")
        elif kind is float:
            if not math.isfinite(item):
                raise bad(where, "non-finite number")
        elif kind is str:
            if len(item) > MAX_VALUE_STRING:
                raise bad(where, f"string of {len(item)} characters exceeds the limit of {MAX_VALUE_STRING}",
                          "resource_limit")
            _utf8(item, where)
        elif kind is list:
            stack.extend((x, ptr(where, i), depth + 1) for i, x in enumerate(item))
        elif kind is dict:
            for key, x in item.items():
                if type(key) is not str:
                    raise bad(where, f"object key {show(key)} is not a string")
                if len(key) > MAX_VALUE_STRING:
                    raise bad(where, "object key exceeds the string limit", "resource_limit")
                _utf8(key, where)
                stack.append((x, ptr(where, key), depth + 1))
        else:
            raise bad(where, f"not a JSON value ({kind.__name__})")


class JsonValue:
    """An immutable JSON value: finite numbers, string keys, bounded depth and size.

    Stored as canonical JSON text, so equality ignores object key order; ``get()``
    returns a fresh copy each time.
    """

    __slots__ = ("_text",)

    def __init__(self, value: Any, path: str = ""):
        if isinstance(value, JsonValue):
            text = value._text
        else:
            _check_json(value, path)
            text = canonical_json(value)
        object.__setattr__(self, "_text", text)

    def __setattr__(self, name, value):
        raise AttributeError("JsonValue is immutable")

    def __delattr__(self, name):
        raise AttributeError("JsonValue is immutable")

    @property
    def text(self) -> str:
        return self._text

    def get(self) -> Any:
        return json.loads(self._text)

    def __eq__(self, other) -> bool:
        return isinstance(other, JsonValue) and other._text == self._text

    def __hash__(self) -> int:
        return hash(self._text)

    def __repr__(self) -> str:
        text = self._text if len(self._text) <= _SHOW else self._text[:_SHOW - 3] + "..."
        return f"JsonValue({text})"


# ── Documents ───────────────────────────────────────────────────────


def check_header(document: Any, artifact: str, version: int, path: str = "") -> None:
    """Artifact discriminator, then schema version, before any key-set check.

    A wrong discriminator or a non-integer version is ``malformed``; an integer version
    other than ``version`` is ``unsupported_schema``.
    """
    if not isinstance(document, dict):
        raise bad(path, "expected a JSON object")
    if "artifact" not in document:
        raise bad(path, "missing keys ['artifact']")
    if document["artifact"] != artifact:
        raise bad(ptr(path, "artifact"), f"expected artifact {artifact!r}, got {show(document['artifact'])}")
    if "schema_version" not in document:
        raise bad(path, "missing keys ['schema_version']")
    found = document["schema_version"]
    if type(found) is not int:
        raise bad(ptr(path, "schema_version"), f"expected an integer, got {show(found)}")
    if found != version:
        raise bad(ptr(path, "schema_version"),
                  f"unsupported {artifact} schema version {found}; this GenAIRR reads version {version}",
                  "unsupported_schema")


def parse_json(text: Any, *, max_bytes: int) -> Any:
    """Strict JSON: duplicate keys, non-finite numbers and trailing data are ``malformed``."""
    if type(text) is not str:
        raise bad("", f"expected JSON text (str), got {type(text).__name__}")
    try:
        return strict_loads(text, max_bytes=max_bytes)
    except ArtifactIOError as exc:
        raise bad("", str(exc), "resource_limit" if exc.kind == "resource" else "malformed") from exc


def read_json(path, *, max_bytes: int) -> Any:
    """Read and strictly parse a UTF-8 JSON file of at most ``max_bytes``."""
    try:
        text = read_text_bounded(path, max_bytes=max_bytes)
    except ArtifactIOError as exc:
        kind = {"io": "io", "resource": "resource_limit"}.get(exc.kind, "malformed")
        raise bad("", f"{path}: {exc}", kind) from exc
    return parse_json(text, max_bytes=max_bytes)


def digest(value: Any, domain: str) -> str:
    """``typed_digest``; a value beyond the native parser's limits is a ``resource_limit`` error."""
    try:
        return typed_digest(value, domain)
    except (ValueError, TypeError, RuntimeError) as exc:
        kind = "resource_limit" if "resource limit" in str(exc).lower() else "malformed"
        raise bad("", f"cannot compute the {domain} digest: {exc}", kind) from exc


def document_bytes(document: dict) -> bytes:
    """A generated JSON document's bytes: canonical JSON, one LF, UTF-8."""
    return (canonical_json(document) + "\n").encode("utf-8")


def digest_of_canonical(text: str, domain: str) -> str:
    """:func:`digest` of a value given as its canonical JSON ``text`` (``canonical_json(value)``),
    without serialising it again: ``typed_digest`` is the same native digest over that text."""
    try:
        return typed_digest_of_canonical(text, domain)
    except (ValueError, TypeError, RuntimeError) as exc:
        kind = "resource_limit" if "resource limit" in str(exc).lower() else "malformed"
        raise bad("", f"cannot compute the {domain} digest: {exc}", kind) from exc
