from __future__ import annotations

import time
import uuid
from dataclasses import dataclass, field
from typing import Any

_cache: dict[str, "CacheEntry"] = {}
_TTL = 1800.0  # 30 minutes


@dataclass
class CacheEntry:
    intermediates: list
    intermediates_all: list
    node_mols: dict
    graph: Any | None
    score_matrix: Any | None
    name_a: str
    name_b: str
    warnings: list[str] = field(default_factory=list)
    created_at: float = field(default_factory=time.time)
    last_accessed: float = field(default_factory=time.time)


def store(entry: CacheEntry) -> str:
    """Store a cache entry and return a new session_id (12-char hex)."""
    cleanup()
    session_id = uuid.uuid4().hex[:12]
    _cache[session_id] = entry
    return session_id


def get(session_id: str) -> CacheEntry | None:
    """Retrieve a cache entry, extending TTL on access. Returns None if expired/missing."""
    entry = _cache.get(session_id)
    if entry is None:
        return None
    now = time.time()
    if now - entry.last_accessed > _TTL:
        del _cache[session_id]
        return None
    entry.last_accessed = now
    return entry


def cleanup() -> None:
    """Remove entries that have exceeded the TTL."""
    now = time.time()
    expired = [sid for sid, e in list(_cache.items()) if now - e.last_accessed > _TTL]
    for sid in expired:
        del _cache[sid]
