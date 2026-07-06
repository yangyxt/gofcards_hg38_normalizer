from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import requests

BASE = "https://java.genemed.tech/admin-api/backend/data/hg19"
PUBLIC_EXCEL_URL = "https://download.genemed.tech/upload/GainFunCards/gofcards_data_download.xlsx"

HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
        "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/125 Safari/537.36"
    ),
    "Accept": "application/json, text/plain, */*",
    "Referer": "https://www.genemed.tech/gofcards/",
    "Origin": "https://www.genemed.tech",
    "Admin-Work-Version": "x",
    "Content-Type": "application/json; charset=UTF-8",
}

TABLE_ENDPOINTS = {
    "SNV": f"{BASE}/GainFunCards_SNV/geneSymbol/page",
    "Indel": f"{BASE}/GainFunCards_Indel/geneSymbol/page",
}

SUMMARY_ENDPOINT = f"{BASE}/variantLevel/summary"


@dataclass(frozen=True)
class RetryConfig:
    attempts: int = 4
    sleep_seconds: float = 1.5
    timeout_seconds: int = 60


def _json_or_error(resp: requests.Response) -> dict[str, Any]:
    try:
        payload = resp.json()
    except ValueError as exc:
        text = resp.text[:500].replace("\n", " ")
        raise RuntimeError(f"Non-JSON response from {resp.url}: HTTP {resp.status_code}: {text}") from exc
    if resp.status_code >= 400:
        raise RuntimeError(f"HTTP {resp.status_code} from {resp.url}: {payload}")
    if payload.get("code") not in (0, "0", None):
        raise RuntimeError(f"GoFCards API error from {resp.url}: {payload}")
    return payload


def request_json(
    session: requests.Session,
    method: str,
    url: str,
    retry: RetryConfig | None = None,
    **kwargs: Any,
) -> dict[str, Any]:
    retry = retry or RetryConfig()
    last_error: Exception | None = None
    for attempt in range(1, retry.attempts + 1):
        try:
            resp = session.request(
                method,
                url,
                headers=HEADERS,
                timeout=retry.timeout_seconds,
                **kwargs,
            )
            return _json_or_error(resp)
        except Exception as exc:  # requests raises several transient subclasses.
            last_error = exc
            if attempt == retry.attempts:
                break
            time.sleep(retry.sleep_seconds * attempt)
    raise RuntimeError(f"Failed {method} {url} after {retry.attempts} attempts") from last_error


def table_payload(page: int, page_size: int) -> dict[str, Any]:
    return {
        "reference": "hg19",
        "queryby": "geneSymbol",
        "terms": "",
        "page": page,
        "pageNo": page,
        "currentPage": page,
        "pageSize": page_size,
    }


def fetch_table_records(
    variant_type: str,
    page_size: int = 5000,
    retry: RetryConfig | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    variant_type = "Indel" if variant_type.lower() == "indel" else variant_type.upper()
    if variant_type not in TABLE_ENDPOINTS:
        raise ValueError(f"variant_type must be one of {sorted(TABLE_ENDPOINTS)}")
    session = requests.Session()
    url = TABLE_ENDPOINTS[variant_type]
    first = request_json(session, "POST", url, retry=retry, json=table_payload(1, page_size))
    data = first.get("data") or {}
    pages = int(data.get("pages") or 1)
    raw_pages = [first]
    records = list(data.get("records") or [])
    for page in range(2, pages + 1):
        payload = request_json(session, "POST", url, retry=retry, json=table_payload(page, page_size))
        raw_pages.append(payload)
        records.extend((payload.get("data") or {}).get("records") or [])
    for rec in records:
        rec["variant_type"] = variant_type
    return records, raw_pages


def summary_params(record: dict[str, Any]) -> dict[str, str]:
    variant_type = str(record.get("variant_type") or record.get("Type") or "SNV")
    variant_type = "Indel" if variant_type.lower() == "indel" else "SNV"
    return {
        "projectCode": "GoFCards",
        "variantLevelType": variant_type,
        "chr": str(record.get("Chr", "")).replace("chr", ""),
        "start": str(record.get("Start", "")),
        "end": str(record.get("End", "")),
        "ref": str(record.get("Ref", "")),
        "alt": str(record.get("Alt", "")),
    }


def fetch_summary(record: dict[str, Any], retry: RetryConfig | None = None) -> dict[str, Any]:
    session = requests.Session()
    return request_json(session, "GET", SUMMARY_ENDPOINT, retry=retry, params=summary_params(record))


def write_jsonl(path: str | Path, rows: list[dict[str, Any]]) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")


def append_jsonl(path: str | Path, row: dict[str, Any]) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")


def read_jsonl(path: str | Path) -> list[dict[str, Any]]:
    in_path = Path(path)
    if not in_path.exists():
        return []
    rows: list[dict[str, Any]] = []
    with in_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows
