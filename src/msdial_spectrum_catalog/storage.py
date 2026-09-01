from __future__ import annotations

import sqlite3
from contextlib import contextmanager
from pathlib import Path

from .schema import SCHEMA_SQL, SCHEMA_VERSION


def connect(path: str | Path) -> sqlite3.Connection:
    connection = sqlite3.connect(Path(path))
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys = ON")
    return connection


def initialize(path: str | Path) -> None:
    database = Path(path)
    database.parent.mkdir(parents=True, exist_ok=True)
    connection = connect(database)
    try:
        connection.executescript(SCHEMA_SQL)
        _ensure_column(connection, "artifact", "source_path", "TEXT NOT NULL DEFAULT ''")
        _ensure_column(connection, "sample", "raw_file_path", "TEXT")
        connection.execute(
            "INSERT OR REPLACE INTO catalog_meta(key, value) VALUES ('schema_version', ?)",
            (str(SCHEMA_VERSION),),
        )
        connection.commit()
    finally:
        connection.close()


def _ensure_column(connection: sqlite3.Connection, table: str, column: str, definition: str) -> None:
    columns = {row[1] for row in connection.execute(f"PRAGMA table_info({table})")}
    if column not in columns:
        connection.execute(f"ALTER TABLE {table} ADD COLUMN {column} {definition}")


@contextmanager
def transaction(path: str | Path):
    connection = connect(path)
    try:
        connection.execute("BEGIN")
        yield connection
        connection.commit()
    except Exception:
        connection.rollback()
        raise
    finally:
        connection.close()
