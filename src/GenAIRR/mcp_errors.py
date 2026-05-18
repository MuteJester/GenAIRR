"""Error codes and response-envelope decorator for the GenAIRR MCP server.

Every MCP tool wraps its return through the `envelope` decorator, which
produces a uniform success/failure shape:

    success: {"ok": True,  "tool": "<name>", "elapsed_ms": int, "result": {...}}
    failure: {"ok": False, "tool": "<name>", "elapsed_ms": int, "error": {...}}

Tools surface domain-specific failures by raising `MCPError(code, message, hint)`
inside the wrapped function -- the decorator catches it and produces the
failure envelope. Any other exception is mapped to error.code == "internal_error".

Stable error-code constants live here so agents and tests can branch on them
across releases. Adding a new code is fine; renaming an existing code is a
breaking change to the MCP surface.
"""
from __future__ import annotations

import functools
import time
from typing import Any, Callable, Dict, Optional, TypeVar


# -- Stable error-code tokens ---------------------------------------

CONFIG_NOT_FOUND = "config_not_found"
ALLELE_NOT_FOUND = "allele_not_found"
INVALID_PRESET = "invalid_preset"
INVALID_PARAMETER = "invalid_parameter"
MALFORMED_RECORD = "malformed_record"
SEED_REPLAY_MISMATCH = "seed_replay_mismatch"


# -- Domain exception type ------------------------------------------


class MCPError(Exception):
    """Raised by a tool implementation to surface a structured error.

    The `envelope` decorator catches this and produces the standard
    failure envelope with the supplied code, message, and optional hint.
    Any unhandled exception (not an MCPError) is mapped to
    error.code == "internal_error" with the exception's class + str as
    the message -- a fallback so tools never silently 500.
    """

    def __init__(self, code: str, message: str, hint: Optional[str] = None) -> None:
        super().__init__(message)
        self.code = code
        self.message = message
        self.hint = hint


# -- Envelope decorator ---------------------------------------------

F = TypeVar("F", bound=Callable[..., Dict[str, Any]])


def envelope(tool_name: str) -> Callable[[F], F]:
    """Wrap a tool implementation in the standard response envelope.

    Usage:

        @mcp.tool()
        @envelope("list_configs")
        def list_configs(...) -> Dict[str, Any]:
            ...
            return result_dict

    The decorated function should return ONLY the inner `result` dict --
    the envelope wraps it with ok/tool/elapsed_ms keys. Raise MCPError
    inside the function for domain failures; the decorator translates
    them to the failure envelope.
    """

    def decorator(fn: F) -> F:
        @functools.wraps(fn)
        def wrapper(*args: Any, **kwargs: Any) -> Dict[str, Any]:
            start = time.perf_counter()
            try:
                result = fn(*args, **kwargs)
                elapsed_ms = int((time.perf_counter() - start) * 1000)
                return {
                    "ok": True,
                    "tool": tool_name,
                    "elapsed_ms": elapsed_ms,
                    "result": result,
                }
            except MCPError as exc:
                elapsed_ms = int((time.perf_counter() - start) * 1000)
                error: Dict[str, Any] = {"code": exc.code, "message": exc.message}
                if exc.hint is not None:
                    error["hint"] = exc.hint
                return {
                    "ok": False,
                    "tool": tool_name,
                    "elapsed_ms": elapsed_ms,
                    "error": error,
                }
            except Exception as exc:  # noqa: BLE001 -- last-resort safety net
                elapsed_ms = int((time.perf_counter() - start) * 1000)
                return {
                    "ok": False,
                    "tool": tool_name,
                    "elapsed_ms": elapsed_ms,
                    "error": {
                        "code": "internal_error",
                        "message": f"{type(exc).__name__}: {exc}",
                    },
                }

        return wrapper  # type: ignore[return-value]

    return decorator
