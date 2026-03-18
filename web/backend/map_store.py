"""SQLite job store for Map Mode v1 (PairMap engine) jobs."""
from __future__ import annotations

import json
import sqlite3
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Optional

from .models import MapJobStatus

MAP_JOBS_DB = Path(__file__).parent.parent / "jobs" / "map_jobs.db"


def _parse_dt(s: Optional[str]) -> Optional[str]:
    return s  # stored as ISO string, returned as-is for the model


@contextmanager
def _conn():
    MAP_JOBS_DB.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(str(MAP_JOBS_DB))
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA journal_mode=WAL")
    try:
        yield con
        con.commit()
    finally:
        con.close()


def init_db() -> None:
    with _conn() as con:
        con.execute("""
            CREATE TABLE IF NOT EXISTS map_jobs (
                id TEXT PRIMARY KEY,
                status TEXT NOT NULL DEFAULT 'queued',
                engine TEXT NOT NULL,
                config_json TEXT NOT NULL DEFAULT '{}',
                created_at TEXT NOT NULL,
                started_at TEXT,
                completed_at TEXT,
                error TEXT,
                progress TEXT
            )
        """)


def create_job(job_id: str, engine: str, config: dict) -> MapJobStatus:
    now = datetime.utcnow().isoformat()
    with _conn() as con:
        con.execute(
            "INSERT INTO map_jobs (id, status, engine, config_json, created_at) VALUES (?, ?, ?, ?, ?)",
            (job_id, "queued", engine, json.dumps(config), now),
        )
    return get_job(job_id)


def get_job(job_id: str) -> Optional[MapJobStatus]:
    with _conn() as con:
        row = con.execute("SELECT * FROM map_jobs WHERE id = ?", (job_id,)).fetchone()
    if row is None:
        return None
    return _row_to_status(row)


def list_jobs() -> list[MapJobStatus]:
    with _conn() as con:
        rows = con.execute("SELECT * FROM map_jobs ORDER BY created_at DESC").fetchall()
    return [_row_to_status(r) for r in rows]


def update_job(job_id: str, **kwargs) -> None:
    allowed = {"status", "started_at", "completed_at", "error", "progress"}
    updates = {k: v for k, v in kwargs.items() if k in allowed}
    if not updates:
        return
    set_clause = ", ".join(f"{k} = ?" for k in updates)
    values = list(updates.values()) + [job_id]
    with _conn() as con:
        con.execute(f"UPDATE map_jobs SET {set_clause} WHERE id = ?", values)


def _row_to_status(row: sqlite3.Row) -> MapJobStatus:
    return MapJobStatus(
        id=row["id"],
        status=row["status"],
        engine=row["engine"],
        config=json.loads(row["config_json"]),
        created_at=row["created_at"],
        started_at=row["started_at"],
        completed_at=row["completed_at"],
        error=row["error"],
        progress=row["progress"],
    )
