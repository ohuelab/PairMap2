"""Thread-safe SQLite-backed (or in-memory) cache for pairwise molecule scores."""
import sqlite3
import hashlib
import json
import threading
import time
from typing import Optional

import numpy as np


class ScoreCache:
    """Thread-safe SQLite-backed cache for pairwise molecule similarity scores.

    When ``db_path`` is ``None`` the cache lives entirely in a dict in memory
    and is not persisted across process restarts.  Pass a file path to enable
    an SQLite backing store that survives between runs.
    """

    def __init__(self, db_path: Optional[str] = None):
        self._lock = threading.Lock()
        self._db_path = db_path
        if db_path is not None:
            self._conn = sqlite3.connect(db_path, check_same_thread=False)
            self._init_db()
        else:
            self._conn = None
            self._mem_cache: dict = {}

    def _init_db(self):
        with self._conn:
            self._conn.execute(
                """
                CREATE TABLE IF NOT EXISTS scores (
                    smiles_a     TEXT    NOT NULL,
                    smiles_b     TEXT    NOT NULL,
                    options_hash TEXT    NOT NULL,
                    score        REAL    NOT NULL,
                    created_at   REAL    NOT NULL,
                    PRIMARY KEY (smiles_a, smiles_b, options_hash)
                )
                """
            )
            self._conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_pair ON scores(smiles_a, smiles_b)"
            )

    @staticmethod
    def _normalize_key(smiles_a: str, smiles_b: str):
        """Ensure smiles_a <= smiles_b for a canonical pair key."""
        if smiles_a <= smiles_b:
            return smiles_a, smiles_b
        return smiles_b, smiles_a

    @staticmethod
    def _options_hash(options: dict) -> str:
        return hashlib.md5(
            json.dumps(sorted(options.items())).encode()
        ).hexdigest()[:8]

    def get(
        self, smiles_a: str, smiles_b: str, options: Optional[dict] = None
    ) -> Optional[float]:
        """Return cached score or ``None`` if not present."""
        key_a, key_b = self._normalize_key(smiles_a, smiles_b)
        opts_hash = self._options_hash(options or {})
        with self._lock:
            if self._conn is not None:
                cur = self._conn.execute(
                    "SELECT score FROM scores "
                    "WHERE smiles_a=? AND smiles_b=? AND options_hash=?",
                    (key_a, key_b, opts_hash),
                )
                row = cur.fetchone()
                return row[0] if row else None
            else:
                return self._mem_cache.get((key_a, key_b, opts_hash))

    def put(
        self,
        smiles_a: str,
        smiles_b: str,
        score: float,
        options: Optional[dict] = None,
    ):
        """Store a score in the cache."""
        key_a, key_b = self._normalize_key(smiles_a, smiles_b)
        opts_hash = self._options_hash(options or {})
        with self._lock:
            if self._conn is not None:
                with self._conn:
                    self._conn.execute(
                        "INSERT OR REPLACE INTO scores VALUES (?, ?, ?, ?, ?)",
                        (key_a, key_b, opts_hash, score, time.time()),
                    )
            else:
                self._mem_cache[(key_a, key_b, opts_hash)] = score

    def get_matrix(
        self, smiles_list: list, options: Optional[dict] = None
    ) -> Optional[np.ndarray]:
        """Return a full N×N score matrix if every pair is cached, else ``None``."""
        n = len(smiles_list)
        matrix = np.zeros((n, n))
        opts_hash = self._options_hash(options or {})
        with self._lock:
            for i in range(n):
                for j in range(i + 1, n):
                    key_a, key_b = self._normalize_key(smiles_list[i], smiles_list[j])
                    if self._conn is not None:
                        cur = self._conn.execute(
                            "SELECT score FROM scores "
                            "WHERE smiles_a=? AND smiles_b=? AND options_hash=?",
                            (key_a, key_b, opts_hash),
                        )
                        row = cur.fetchone()
                        if row is None:
                            return None
                        score = row[0]
                    else:
                        score = self._mem_cache.get((key_a, key_b, opts_hash))
                        if score is None:
                            return None
                    matrix[i][j] = score
                    matrix[j][i] = score
        return matrix

    def close(self):
        if self._conn is not None:
            self._conn.close()
