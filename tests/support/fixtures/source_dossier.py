"""A child-process script that proves a dossier loads with no network.

Shared by the dossier tests and the packaged-starter tests, which both
need to assert that an offline load touches no socket.
"""
from __future__ import annotations


OFFLINE_SCRIPT = r'''
import _socket
import json
import socket
import sys
import warnings
from pathlib import Path

GUARD = "network access is disabled in this offline check"
ATTEMPTS = []


def _refuse(*args, **kwargs):
    ATTEMPTS.append(repr(args)[:200])          # counted even if a caller swallows the error
    raise OSError(GUARD)


def _audit(event, args):
    # Below the Python socket surface: every connection, datagram send and name lookup raises an audit event.
    if event in ("socket.connect", "socket.sendto", "socket.sendmsg", "socket.getaddrinfo", "socket.gethostbyname",
                 "socket.gethostbyaddr", "socket.getnameinfo"):
        _refuse(event, *args)


class _NoSocket(socket.socket):
    def __init__(self, *args, **kwargs):
        _refuse(*args)


sys.addaudithook(_audit)
socket.socket = _NoSocket
socket.create_connection = _refuse
socket.getaddrinfo = _refuse
warnings.simplefilter("error")

def _raw_connect():
    sock = _socket.socket()
    try:
        sock.connect(("127.0.0.1", 9))
    finally:
        sock.close()


def _datagram():
    sock = _socket.socket(_socket.AF_INET, _socket.SOCK_DGRAM)
    try:
        sock.sendto(b"x", ("127.0.0.1", 9))
    finally:
        sock.close()


# Positive controls: each route fails with the guard's own error, and is counted.
CONTROLS = []
for name, attempt in (("raw socket connect", _raw_connect),
                      ("datagram send", _datagram),
                      ("name lookup", lambda: socket.gethostbyname("localhost")),
                      ("reverse lookup", lambda: _socket.gethostbyaddr("127.0.0.1")),
                      ("address lookup", lambda: _socket.getaddrinfo("localhost", 443)),
                      ("connection", lambda: socket.create_connection(("127.0.0.1", 9), timeout=1))):
    try:
        attempt()
        CONTROLS.append((name, "connected"))
    except OSError as exc:
        CONTROLS.append((name, str(exc)))
CONTROL_ATTEMPTS = len(ATTEMPTS)
del ATTEMPTS[:]

from GenAIRR import Experiment
from GenAIRR.reference_sources import OgrdbSource, ReferenceBuild

dossier, out = Path(sys.argv[1]), Path(sys.argv[2])
checked = ReferenceBuild.read_dossier(dossier)
locked = OgrdbSource.read(dossier / "source.lock.json", cache=dossier / "cache")
rebuilt = locked.build(dossier / "recipe.json")
rebuilt.write_dossier(out)
reference = ReferenceBuild.read_dossier(out).reference
result = Experiment.on(reference).recombine().run_records(n=25, seed=7, validate_records=True)
print(json.dumps({
    "checked_execution_hash": checked.reference.execution_hash,
    "rebuilt_execution_hash": rebuilt.reference.execution_hash,
    "reference_content_hash": result.reference_content_hash,
    "execution_hash": reference.execution_hash,
    "records": len(result),
    "alleles": sorted({record["v_call"] for record in result} | {record["j_call"] for record in result}),
    "network_attempts": ATTEMPTS,
    "controls": CONTROLS,
    "control_attempts": CONTROL_ATTEMPTS,
}))
'''
