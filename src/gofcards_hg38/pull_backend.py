from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

from .api import PUBLIC_EXCEL_URL, fetch_table_records, write_jsonl
from .io_utils import ensure_parent, write_excel
from .schema import normalize_backend_frame


def pull_backend(out_xlsx: str | Path, raw_jsonl: str | Path, page_size: int = 5000) -> None:
    all_records: list[dict] = []
    raw_rows: list[dict] = []
    for variant_type in ("SNV", "Indel"):
        print(f"Pulling GoFCards backend {variant_type} table...", file=sys.stderr)
        records, pages = fetch_table_records(variant_type, page_size=page_size)
        print(f"Pulled {len(records)} {variant_type} records.", file=sys.stderr)
        all_records.extend(records)
        for page in pages:
            raw_rows.append({"variant_type": variant_type, "response": page})
    write_jsonl(raw_jsonl, raw_rows)

    df = normalize_backend_frame(pd.DataFrame(all_records))
    summary = pd.DataFrame(
        [
            {"metric": "backend_total_records", "value": len(df)},
            {"metric": "backend_snv_records", "value": int((df["variant_type"] == "SNV").sum())},
            {"metric": "backend_indel_records", "value": int((df["variant_type"] == "Indel").sum())},
            {"metric": "raw_jsonl", "value": str(raw_jsonl)},
        ]
    )
    write_excel(out_xlsx, {"backend_records": df, "summary": summary})


def download_public_excel(url: str | None, out_xlsx: str | Path) -> None:
    import requests

    url = url or PUBLIC_EXCEL_URL
    out = ensure_parent(out_xlsx)
    headers = {"User-Agent": "Mozilla/5.0", "Accept": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet,*/*"}
    resp = requests.get(url, headers=headers, timeout=120)
    resp.raise_for_status()
    out.write_bytes(resp.content)
    meta = {
        "url": url,
        "status_code": resp.status_code,
        "content_length": resp.headers.get("Content-Length"),
        "last_modified": resp.headers.get("Last-Modified"),
    }
    out.with_suffix(out.suffix + ".metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
