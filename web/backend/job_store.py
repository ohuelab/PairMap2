"""SQLite-based job state management for PairMap2 WebUI."""
from __future__ import annotations

import json
import os
import sqlite3
import time
from contextlib import contextmanager
from pathlib import Path
from typing import List, Optional

from .models import JobResponse, JobStatus, StageTimingModel

# Jobs directory: configurable via env var; default is ./jobs relative to project root
JOBS_DIR = Path(os.environ.get("PAIRMAP_JOBS_DIR", str(Path(__file__).parent.parent / "jobs")))
DB_PATH = JOBS_DIR / "jobs.db"


@contextmanager
def _conn():
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(str(DB_PATH))
    con.row_factory = sqlite3.Row
    # WAL mode for safe concurrent access from worker threads
    con.execute("PRAGMA journal_mode=WAL")
    try:
        yield con
        con.commit()
    finally:
        con.close()


def init_db() -> None:
    """Create the jobs table if it does not exist."""
    with _conn() as con:
        con.execute("""
            CREATE TABLE IF NOT EXISTS jobs (
                job_id       TEXT PRIMARY KEY,
                status       TEXT NOT NULL DEFAULT 'pending',
                created_at   REAL NOT NULL,
                updated_at   REAL NOT NULL,
                config_json  TEXT NOT NULL DEFAULT '{}',
                timings_json TEXT NOT NULL DEFAULT '[]',
                error        TEXT,
                n_nodes      INTEGER,
                n_edges      INTEGER
            )
        """)


def create_job(job_id: str, config: dict) -> JobResponse:
    now = time.time()
    with _conn() as con:
        con.execute(
            """INSERT INTO jobs
               (job_id, status, created_at, updated_at, config_json, timings_json)
               VALUES (?, ?, ?, ?, ?, ?)""",
            (job_id, JobStatus.PENDING.value, now, now, json.dumps(config), "[]"),
        )
    return get_job(job_id)


def get_job(job_id: str) -> Optional[JobResponse]:
    with _conn() as con:
        row = con.execute("SELECT * FROM jobs WHERE job_id = ?", (job_id,)).fetchone()
    if row is None:
        return None
    return _row_to_response(row)


def list_jobs() -> List[JobResponse]:
    with _conn() as con:
        rows = con.execute("SELECT * FROM jobs ORDER BY created_at DESC").fetchall()
    return [_row_to_response(r) for r in rows]


def update_job(job_id: str, **kwargs) -> None:
    """Update allowed fields on a job row.

    Accepted kwargs: status, timings_json, error, n_nodes, n_edges.
    updated_at is always refreshed automatically.
    """
    allowed = {"status", "timings_json", "error", "n_nodes", "n_edges"}
    updates = {k: v for k, v in kwargs.items() if k in allowed}
    updates["updated_at"] = time.time()
    set_clause = ", ".join(f"{k} = ?" for k in updates)
    values = list(updates.values()) + [job_id]
    with _conn() as con:
        con.execute(f"UPDATE jobs SET {set_clause} WHERE job_id = ?", values)


def cancel_job(job_id: str) -> Optional[JobResponse]:
    job = get_job(job_id)
    if job is None:
        return None
    if job.status not in (JobStatus.PENDING, JobStatus.RUNNING):
        return job
    update_job(job_id, status=JobStatus.CANCELLED.value)
    return get_job(job_id)


def _row_to_response(row: sqlite3.Row) -> JobResponse:
    raw_timings = json.loads(row["timings_json"] or "[]")
    timings: list[StageTimingModel] = []
    for t in raw_timings:
        if isinstance(t, dict):
            try:
                timings.append(StageTimingModel(**t))
            except Exception:
                pass
    return JobResponse(
        job_id=row["job_id"],
        status=JobStatus(row["status"]),
        created_at=row["created_at"],
        updated_at=row["updated_at"],
        config=json.loads(row["config_json"] or "{}"),
        timings=timings,
        error=row["error"],
        n_nodes=row["n_nodes"],
        n_edges=row["n_edges"],
    )
