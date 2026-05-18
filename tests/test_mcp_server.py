"""Tests for the redesigned GenAIRR MCP server.

Layered per the spec at docs/superpowers/specs/2026-05-18-mcp-redesign-v2-design.md:
- Tier 1: happy path per tool (14 tests)
- Tier 2: error envelope per error code (~6 tests)
- Tier 3: one end-to-end agent-style chain (1 test)

Plus a foundation block here for the envelope decorator itself.
"""
from __future__ import annotations

import pytest

from GenAIRR.mcp_errors import (
    CONFIG_NOT_FOUND,
    ALLELE_NOT_FOUND,
    INVALID_PRESET,
    INVALID_PARAMETER,
    MALFORMED_RECORD,
    SEED_REPLAY_MISMATCH,
    MCPError,
    envelope,
)


# -- Foundation: envelope decorator ---------------------------------


def test_envelope_wraps_success_with_ok_true_and_elapsed_ms():
    @envelope("test_tool")
    def my_tool() -> dict:
        return {"x": 1}

    resp = my_tool()
    assert resp["ok"] is True
    assert resp["tool"] == "test_tool"
    assert isinstance(resp["elapsed_ms"], int) and resp["elapsed_ms"] >= 0
    assert resp["result"] == {"x": 1}
    assert "error" not in resp


def test_envelope_catches_mcp_error_into_failure_shape():
    @envelope("test_tool")
    def my_tool() -> dict:
        raise MCPError(CONFIG_NOT_FOUND, "Config 'foo' missing.", hint="Try list_configs.")

    resp = my_tool()
    assert resp["ok"] is False
    assert resp["tool"] == "test_tool"
    assert isinstance(resp["elapsed_ms"], int) and resp["elapsed_ms"] >= 0
    assert resp["error"]["code"] == CONFIG_NOT_FOUND
    assert resp["error"]["message"] == "Config 'foo' missing."
    assert resp["error"]["hint"] == "Try list_configs."
    assert "result" not in resp


def test_envelope_catches_unexpected_exception_into_internal_error():
    @envelope("test_tool")
    def my_tool() -> dict:
        raise ValueError("something went wrong")

    resp = my_tool()
    assert resp["ok"] is False
    assert resp["error"]["code"] == "internal_error"
    assert "ValueError" in resp["error"]["message"]
    assert "something went wrong" in resp["error"]["message"]


def test_envelope_preserves_function_metadata():
    @envelope("test_tool")
    def my_tool(x: int) -> dict:
        """My docstring."""
        return {"x": x}

    assert my_tool.__name__ == "my_tool"
    assert my_tool.__doc__ == "My docstring."


def test_error_code_constants_are_stable_strings():
    # The agent branches on these -- they must be plain strings, not enums.
    assert CONFIG_NOT_FOUND == "config_not_found"
    assert ALLELE_NOT_FOUND == "allele_not_found"
    assert INVALID_PRESET == "invalid_preset"
    assert INVALID_PARAMETER == "invalid_parameter"
    assert MALFORMED_RECORD == "malformed_record"
    assert SEED_REPLAY_MISMATCH == "seed_replay_mismatch"
