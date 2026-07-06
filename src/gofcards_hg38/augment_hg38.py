from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

import pandas as pd

from .api import append_jsonl, fetch_summary, read_jsonl
from .io_utils import read_excel, write_excel
from .schema import ALLELE_COLUMNS, allele_key_series, normalize_backend_frame

SUMMARY_KEEP_FIELDS = [
    "hg19_start",
    "hg19_end",
    "hg38_start",
    "hg38_end",
    "rsID",
    "AAChange_refGene",
    "rare",
    "Accession",
    "CLNDN",
    "CLNSIG",
    "CLNREVSTAT",
]


def _read_backend_records(path: str | Path) -> pd.DataFrame:
    try:
        return read_excel(path, "backend_records")
    except ValueError:
        return read_excel(path, 0)


def _cache_key(record: dict[str, Any]) -> str:
    return "|".join(str(record.get(col, "")) for col in ALLELE_COLUMNS)


def _load_cache(path: str | Path) -> dict[str, dict[str, Any]]:
    cache: dict[str, dict[str, Any]] = {}
    for row in read_jsonl(path):
        key = row.get("allele_key")
        if key:
            cache[str(key)] = row
    return cache


def _flatten_summary(cache_row: dict[str, Any]) -> dict[str, Any]:
    flat: dict[str, Any] = {
        "summary_status": cache_row.get("status", ""),
        "summary_error": cache_row.get("error", ""),
    }
    payload = cache_row.get("response") or {}
    data = payload.get("data") if isinstance(payload, dict) else {}
    data = data or {}
    for field in SUMMARY_KEEP_FIELDS:
        flat[f"summary_{field}"] = data.get(field, "")
    for list_name in ("codingList", "nonCodingList", "splicingList"):
        value = data.get(list_name)
        flat[f"summary_{list_name}_json"] = json.dumps(value, ensure_ascii=False) if value else ""
    return flat


def augment_hg38(
    input_xlsx: str | Path,
    out_xlsx: str | Path,
    cache_jsonl: str | Path,
    sleep_seconds: float = 0.15,
) -> None:
    df = normalize_backend_frame(_read_backend_records(input_xlsx))
    df["allele_key"] = allele_key_series(df)
    unique = df.drop_duplicates("allele_key")[["allele_key", *ALLELE_COLUMNS]].copy()

    cache = _load_cache(cache_jsonl)
    for _, row in unique.iterrows():
        key = str(row["allele_key"])
        if key in cache:
            continue
        record = {col: row[col] for col in ALLELE_COLUMNS}
        try:
            response = fetch_summary(record)
            cache_row = {"allele_key": key, "request": record, "status": "ok", "response": response}
        except Exception as exc:
            cache_row = {"allele_key": key, "request": record, "status": "error", "error": str(exc)}
        append_jsonl(cache_jsonl, cache_row)
        cache[key] = cache_row
        time.sleep(sleep_seconds)

    flat_rows = []
    for _, row in unique.iterrows():
        key = str(row["allele_key"])
        flat = {"allele_key": key}
        flat.update(_flatten_summary(cache.get(key, {"status": "missing_cache"})))
        flat_rows.append(flat)
    flat_df = pd.DataFrame(flat_rows)
    augmented = df.merge(flat_df, on="allele_key", how="left")

    summary = pd.DataFrame(
        [
            {"metric": "input_records", "value": len(df)},
            {"metric": "unique_alleles", "value": unique["allele_key"].nunique()},
            {"metric": "summary_ok", "value": int((flat_df["summary_status"] == "ok").sum())},
            {"metric": "summary_error", "value": int((flat_df["summary_status"] == "error").sum())},
            {"metric": "summary_cache_jsonl", "value": str(cache_jsonl)},
        ]
    )
    write_excel(out_xlsx, {"augmented_records": augmented, "unique_summary": flat_df, "summary": summary})

