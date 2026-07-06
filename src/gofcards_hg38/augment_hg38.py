from __future__ import annotations

import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
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
    workers: int = 1,
) -> None:
    df = normalize_backend_frame(_read_backend_records(input_xlsx))
    df["allele_key"] = allele_key_series(df)
    unique = df.drop_duplicates("allele_key")[["allele_key", *ALLELE_COLUMNS]].copy()

    cache = _load_cache(cache_jsonl)
    pending = []
    for index, (_, row) in enumerate(unique.iterrows(), start=1):
        key = str(row["allele_key"])
        if key not in cache:
            pending.append((index, key, {col: row[col] for col in ALLELE_COLUMNS}))

    total = len(unique)
    if pending:
        print(
            f"Querying GoFCards summary for {len(pending)} uncached alleles "
            f"({total} total unique alleles) with workers={workers}",
            file=sys.stderr,
        )

    def query_one(item: tuple[int, str, dict[str, Any]]) -> dict[str, Any]:
        index, key, record = item
        try:
            response = fetch_summary(record)
            return {"allele_key": key, "request": record, "status": "ok", "response": response, "index": index}
        except Exception as exc:
            return {"allele_key": key, "request": record, "status": "error", "error": str(exc), "index": index}

    if workers <= 1:
        for done_count, item in enumerate(pending, start=1):
            index, key, _record = item
            if done_count == 1 or done_count % 100 == 0:
                print(f"Querying GoFCards summary {done_count}/{len(pending)}: {key}", file=sys.stderr)
            cache_row = query_one(item)
            append_jsonl(cache_jsonl, cache_row)
            cache[key] = cache_row
            time.sleep(sleep_seconds)
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(query_one, item) for item in pending]
            for done_count, future in enumerate(as_completed(futures), start=1):
                cache_row = future.result()
                key = str(cache_row["allele_key"])
                if done_count == 1 or done_count % 100 == 0:
                    print(f"Completed GoFCards summary {done_count}/{len(pending)}: {key}", file=sys.stderr)
                append_jsonl(cache_jsonl, cache_row)
                cache[key] = cache_row

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
