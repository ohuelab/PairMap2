"""SQLite job store for Map Mode v1 (PairMap engine) jobs."""
from __future__ import annotations

import json
import os
import sqlite3
from contextlib import contextmanager
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional

from .models import MapJobStatus

_JOBS_DIR = Path(os.environ.get("PAIRMAP_JOBS_DIR", str(Path(__file__).parent.parent / "jobs")))
MAP_JOBS_DB = _JOBS_DIR / "map_jobs.db"


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
                progress TEXT,
                session_id TEXT
            )
        """)
        try:
            con.execute("ALTER TABLE map_jobs ADD COLUMN session_id TEXT")
        except sqlite3.OperationalError:
            pass


def create_job(job_id: str, engine: str, config: dict, session_id: Optional[str] = None) -> MapJobStatus:
    now = datetime.utcnow().isoformat()
    with _conn() as con:
        con.execute(
            "INSERT INTO map_jobs (id, status, engine, config_json, created_at, session_id) VALUES (?, ?, ?, ?, ?, ?)",
            (job_id, "queued", engine, json.dumps(config), now, session_id),
        )
    return get_job(job_id, session_id=session_id)


def get_job(job_id: str, session_id: Optional[str] = None) -> Optional[MapJobStatus]:
    with _conn() as con:
        if session_id is not None:
            row = con.execute(
                "SELECT * FROM map_jobs WHERE id = ? AND session_id = ?", (job_id, session_id)
            ).fetchone()
        else:
            row = con.execute("SELECT * FROM map_jobs WHERE id = ?", (job_id,)).fetchone()
    if row is None:
        return None
    return _row_to_status(row)


def list_jobs(session_id: Optional[str] = None) -> list[MapJobStatus]:
    with _conn() as con:
        if session_id is not None:
            rows = con.execute(
                "SELECT * FROM map_jobs WHERE session_id = ? ORDER BY created_at DESC", (session_id,)
            ).fetchall()
        else:
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


def purge_old_jobs(max_age_days: int = 90) -> list[str]:
    """Delete map_jobs records and artifacts older than max_age_days."""
    import shutil
    from .map_worker import MAP_JOBS_DIR

    cutoff = (datetime.utcnow() - timedelta(days=max_age_days)).isoformat()
    with _conn() as con:
        rows = con.execute(
            "SELECT id FROM map_jobs WHERE created_at < ?", (cutoff,)
        ).fetchall()
        purged = [r["id"] for r in rows]
        if purged:
            con.execute("DELETE FROM map_jobs WHERE created_at < ?", (cutoff,))
    for job_id in purged:
        job_dir = MAP_JOBS_DIR / job_id
        if job_dir.exists():
            shutil.rmtree(job_dir)
    return purged


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
