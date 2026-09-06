from __future__ import annotations

import sqlite3
from contextlib import contextmanager
from pathlib import Path

from .schema import MIGRATIONS, SCHEMA_SQL, SCHEMA_VERSION, AddColumn, Migration, RunSql


class SchemaVersionError(RuntimeError):
    """Raised when a database cannot be handled by this build of the catalog."""


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
        current = read_schema_version(connection)
        # Refuse before writing anything: a newer database opened by older code would otherwise be
        # stamped back down to SCHEMA_VERSION and silently lose its version marker.
        if current is not None and current > SCHEMA_VERSION:
            raise SchemaVersionError(
                f"database {database} reports schema_version {current}, "
                f"but this build supports at most {SCHEMA_VERSION}; upgrade msdial_spectrum_catalog"
            )
        connection.executescript(SCHEMA_SQL)
        applied = apply_migrations(connection, from_version=current)
        descriptions = {migration.version: migration.description for migration in MIGRATIONS}
        for version in applied:
            connection.execute(
                "INSERT OR IGNORE INTO schema_migration(version, description) VALUES (?, ?)",
                (version, descriptions[version]),
            )
        connection.execute(
            "INSERT OR REPLACE INTO catalog_meta(key, value) VALUES ('schema_version', ?)",
            (str(SCHEMA_VERSION),),
        )
        connection.commit()
    finally:
        connection.close()


def read_schema_version(connection: sqlite3.Connection) -> int | None:
    if not _table_exists(connection, "catalog_meta"):
        return None
    row = connection.execute(
        "SELECT value FROM catalog_meta WHERE key = 'schema_version'"
    ).fetchone()
    if row is None:
        return None
    try:
        return int(str(row[0]).strip())
    except ValueError:
        return None


def apply_migrations(connection: sqlite3.Connection, *, from_version: int | None) -> tuple[int, ...]:
    baseline = 0 if from_version is None else from_version
    applied: list[int] = []
    for migration in MIGRATIONS:
        if migration.version <= baseline:
            continue
        _apply_migration(connection, migration)
        applied.append(migration.version)
    return tuple(applied)


def _apply_migration(connection: sqlite3.Connection, migration: Migration) -> None:
    for step in migration.steps:
        if isinstance(step, AddColumn):
            _ensure_column(connection, step.table, step.column, step.definition)
        elif isinstance(step, RunSql):
            connection.execute(step.statement)
        else:
            raise SchemaVersionError(
                f"migration {migration.version} declares unsupported step {step!r}"
            )


def _table_exists(connection: sqlite3.Connection, table: str) -> bool:
    row = connection.execute(
        "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?", (table,)
    ).fetchone()
    return row is not None


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
