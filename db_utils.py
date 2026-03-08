import re
import sqlite3
from contextlib import contextmanager


IDENTIFIER_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


@contextmanager
def connect(db_path):
    conn = sqlite3.connect(db_path)
    try:
        yield conn
        conn.commit()
    finally:
        conn.close()


def validate_identifier(name):
    if not IDENTIFIER_RE.match(name):
        raise ValueError(f"Invalid SQL identifier: {name}")
    return name


def execute(db_path, sql, params=()):
    with connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute(sql, params)
        return cursor.rowcount


def executemany(db_path, sql, params_seq):
    with connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.executemany(sql, params_seq)
        return cursor.rowcount


def query_all(db_path, sql, params=()):
    with connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute(sql, params)
        return cursor.fetchall()


def query_one(db_path, sql, params=()):
    with connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute(sql, params)
        return cursor.fetchone()
