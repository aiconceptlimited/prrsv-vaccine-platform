import os
import pandas as pd
from pathlib import Path
import mysql.connector
BASE_DIR = Path(__file__).resolve().parents[1]
from datetime import datetime, timezone
# VAXINTAIC database configuration.
# Credentials are supplied through the environment and must never be committed.
DB_HOST = os.environ.get("VAXINTAIC_DB_HOST", "localhost")
DB_PORT = int(os.environ.get("VAXINTAIC_DB_PORT", "3306"))
DB_USER = os.environ.get("VAXINTAIC_DB_USER", "vaxuser")
DB_PASSWORD = os.environ.get("VAXINTAIC_DB_PASSWORD")
DB_NAME = os.environ.get("VAXINTAIC_DB_NAME", "vaxintaic")

if not DB_PASSWORD:
    raise RuntimeError(
        "VAXINTAIC_DB_PASSWORD is not set; refusing database connection."
    )

DB_URL = (
    f"mysql+mysqlconnector://{DB_USER}:{DB_PASSWORD}"
    f"@{DB_HOST}:{DB_PORT}/{DB_NAME}"
)

CSV_PATH = BASE_DIR / "data" / "models" / "nanoparticle_model_v3.csv"

df = pd.read_csv(CSV_PATH)

conn = mysql.connector.connect(
    host=DB_HOST,
    port=DB_PORT,
    user=DB_USER,
    password=DB_PASSWORD,
    database=DB_NAME
)
cur = conn.cursor()

try:
    for _, r in df.iterrows():
        cur.execute("""
            INSERT INTO nanoparticle_delivery_v3
            (construct_id, prrsv_type, delivery_index,
             encapsulation_efficiency, delivery_stability, timestamp_utc)
            VALUES (%s,%s,%s,%s,%s,%s)
            ON DUPLICATE KEY UPDATE
              delivery_index = VALUES(delivery_index),
              encapsulation_efficiency = VALUES(encapsulation_efficiency),
              delivery_stability = VALUES(delivery_stability),
              timestamp_utc = VALUES(timestamp_utc)
        """, (
            r["construct_id"],
            r["prrsv_type"],
            float(r["delivery_index"]),
            float(r.get("encapsulation_efficiency", 0)),
            float(r.get("delivery_stability", 0)),
            datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
        ))

    conn.commit()
except Exception:
    conn.rollback()
    raise
finally:
    cur.close()
    conn.close()

print("✅ Nanoparticle delivery metrics loaded.")

