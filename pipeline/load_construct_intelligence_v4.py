#!/usr/bin/env python3
from pathlib import Path
from datetime import datetime, timezone
import os
import json
import pandas as pd
import mysql.connector

BASE_DIR = Path(__file__).resolve().parents[1]
INPUT_CSV = BASE_DIR / "data" / "reports" / "model_scores_v3.csv"

DB_HOST = os.environ.get("VAXINTAIC_DB_HOST", "localhost")
DB_PORT = int(os.environ.get("VAXINTAIC_DB_PORT", "3306"))
DB_USER = os.environ.get("VAXINTAIC_DB_USER", "vaxuser")
DB_PASSWORD = os.environ.get("VAXINTAIC_DB_PASSWORD")
DB_NAME = os.environ.get("VAXINTAIC_DB_NAME", "vaxintaic")

if not DB_PASSWORD:
    raise SystemExit("VAXINTAIC_DB_PASSWORD is not set; refusing database write.")

if not INPUT_CSV.exists():
    raise SystemExit(f"Missing corrected Stage 07 artifact: {INPUT_CSV}")

df = pd.read_csv(INPUT_CSV)

required = {
    "construct_id", "prrsv_type", "num_epitopes", "length_nt",
    "GC_percent", "length_kb", "stability_index",
    "encapsulation_efficiency", "delivery_stability",
    "delivery_index", "mean_immunogenicity",
    "predicted_model_score", "vaccine_index",
}

missing = required - set(df.columns)
if missing:
    raise SystemExit(f"Missing required corrected columns: {sorted(missing)}")

for forbidden in ("CAI", "predicted_efficacy"):
    if forbidden in df.columns:
        raise SystemExit(f"Forbidden legacy field present: {forbidden}")

if df.empty:
    raise SystemExit("Corrected Stage 07 artifact contains no constructs.")

timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")

conn = mysql.connector.connect(
    host=DB_HOST,
    port=DB_PORT,
    user=DB_USER,
    password=DB_PASSWORD,
    database=DB_NAME,
)
cur = conn.cursor()

try:
    for row in df.to_dict(orient="records"):
        payload = {
            key: row[key]
            for key in required
        }
        payload["codon_adaptation_status"] = "unavailable"
        payload["codon_adaptation_reference"] = None

        cur.execute(
            """
            INSERT INTO construct_intelligence_v3
            (construct_id, prrsv_type, data, timestamp_utc)
            VALUES (%s, %s, %s, %s)
            ON DUPLICATE KEY UPDATE
              data = VALUES(data),
              timestamp_utc = VALUES(timestamp_utc)
            """,
            (
                str(row["construct_id"]),
                str(row["prrsv_type"]),
                json.dumps(payload, sort_keys=True),
                timestamp,
            ),
        )

    conn.commit()
except Exception:
    conn.rollback()
    raise
finally:
    cur.close()
    conn.close()

print(f"Committed {len(df)} corrected construct records.")
