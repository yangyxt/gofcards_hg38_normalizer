from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .io_utils import ensure_parent, read_excel, write_excel


def _read_refalt(path: str | Path) -> pd.DataFrame:
    try:
        return read_excel(path, "refalt_checked")
    except ValueError:
        return read_excel(path, 0)


def _parse_aachange(value: object) -> list[dict[str, str]]:
    text = "" if value is None or pd.isna(value) else str(value)
    rows: list[dict[str, str]] = []
    for item in text.split(","):
        parts = item.strip().split(":")
        if len(parts) < 5:
            continue
        gene, transcript = parts[0], parts[1]
        cdna = next((p for p in parts if p.startswith("c.")), "")
        protein = next((p for p in parts if p.startswith("p.")), "")
        if cdna:
            cdna = _annovar_cdna_to_hgvs(cdna)
        rows.append({"gene": gene, "transcript": transcript, "cdna": cdna, "protein": protein})
    return rows


def _annovar_cdna_to_hgvs(cdna: str) -> str:
    match = re.fullmatch(r"c\.([ACGTN])(\d+)([ACGTN])", cdna, flags=re.IGNORECASE)
    if match:
        ref, pos, alt = match.groups()
        return f"c.{pos}{ref.upper()}>{alt.upper()}"
    return cdna


def write_transvar_queries(input_xlsx: str | Path, out_dir: str | Path) -> None:
    df = _read_refalt(input_xlsx)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    query_rows: list[dict[str, str]] = []
    canno: list[str] = []
    panno: list[str] = []
    for _, row in df.iterrows():
        allele_key = str(row.get("allele_key", ""))
        source = row.get("AAChange_refGene", "") or row.get("summary_AAChange_refGene", "")
        for parsed in _parse_aachange(source):
            if parsed["cdna"]:
                query = f"{parsed['gene']}:{parsed['cdna']}"
                canno.append(query)
                query_rows.append({**parsed, "query_type": "canno", "query": query, "allele_key": allele_key})
            if parsed["protein"]:
                query = f"{parsed['gene']}:{parsed['protein']}"
                panno.append(query)
                query_rows.append({**parsed, "query_type": "panno", "query": query, "allele_key": allele_key})

    (out / "transvar_canno_queries.txt").write_text("\n".join(sorted(set(canno))) + "\n", encoding="utf-8")
    (out / "transvar_panno_queries.txt").write_text("\n".join(sorted(set(panno))) + "\n", encoding="utf-8")
    runner = out / "run_transvar.sh"
    runner.write_text(
        """#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
transvar canno -l transvar_canno_queries.txt --refversion hg19 --ensembl > transvar.canno.hg19.txt
transvar panno -l transvar_panno_queries.txt --refversion hg19 --ensembl > transvar.panno.hg19.txt
""",
        encoding="utf-8",
    )
    runner.chmod(0o755)
    write_excel(out / "transvar_query_map.xlsx", {"transvar_query_map": pd.DataFrame(query_rows)})


def read_transvar_outputs(transvar_dir: str | Path) -> dict[str, pd.DataFrame]:
    out: dict[str, pd.DataFrame] = {}
    base = Path(transvar_dir)
    for name in ("transvar.canno.hg19.txt", "transvar.panno.hg19.txt"):
        path = base / name
        if path.exists():
            rows = [{"line": line.rstrip("\n")} for line in path.read_text(encoding="utf-8", errors="replace").splitlines()]
            out[name[:31]] = pd.DataFrame(rows)
    query_map = base / "transvar_query_map.xlsx"
    if query_map.exists():
        out["transvar_query_map"] = read_excel(query_map, 0)
    return out

