"""Central VAXINTAIC database configuration."""

from __future__ import annotations

import os
from urllib.parse import quote_plus


DB_HOST = os.environ.get("VAXINTAIC_DB_HOST", "localhost")
DB_PORT = int(os.environ.get("VAXINTAIC_DB_PORT", "3306"))
DB_USER = os.environ.get("VAXINTAIC_DB_USER", "vaxuser")
DB_NAME = os.environ.get("VAXINTAIC_DB_NAME", "vaxintaic")
DB_PASSWORD = os.environ.get("VAXINTAIC_DB_PASSWORD")

if not DB_PASSWORD:
    raise RuntimeError(
        "VAXINTAIC_DB_PASSWORD is not set; refusing database access."
    )


def mysql_connector_config() -> dict:
    return {
        "host": DB_HOST,
        "port": DB_PORT,
        "user": DB_USER,
        "password": DB_PASSWORD,
        "database": DB_NAME,
    }


def sqlalchemy_url(driver: str = "mysqlconnector") -> str:
    return (
        f"mysql+{driver}://{DB_USER}:{quote_plus(DB_PASSWORD)}"
        f"@{DB_HOST}:{DB_PORT}/{DB_NAME}"
    )
