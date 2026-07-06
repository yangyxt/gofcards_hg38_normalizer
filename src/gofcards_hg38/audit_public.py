from __future__ import annotations

from pathlib import Path

import pandas as pd

from .io_utils import read_excel, write_excel
from .schema import ALLELE_COLUMNS, allele_key_series, normalize_backend_frame, normalize_public_frame


def _read_backend_records(path: str | Path) -> pd.DataFrame:
    try:
        return read_excel(path, "backend_records")
    except ValueError:
        return read_excel(path, 0)


def _read_public_records(path: str | Path) -> pd.DataFrame:
    sheets = pd.read_excel(path, sheet_name=None, dtype=object)
    if "total3161" in sheets:
        df = sheets["total3161"]
    else:
        first_name = next(iter(sheets))
        df = sheets[first_name]
    return df.dropna(how="all").copy()


def audit_public_excel(backend_xlsx: str | Path, public_xlsx: str | Path, out_xlsx: str | Path) -> None:
    backend = normalize_backend_frame(_read_backend_records(backend_xlsx))
    public = normalize_public_frame(_read_public_records(public_xlsx))

    backend = backend.copy()
    public = public.copy()
    backend["allele_key"] = allele_key_series(backend)
    public["allele_key"] = allele_key_series(public)
    backend["backend_row_number"] = range(2, len(backend) + 2)
    public["public_row_number"] = range(2, len(public) + 2)

    backend_keys = backend[["allele_key", *ALLELE_COLUMNS, "Gene_Symbol", "Transcript", "backend_row_number"]].drop_duplicates()
    public_keys = public[["allele_key", *ALLELE_COLUMNS, "Gene_Symbol", "Transcript", "public_row_number"]].drop_duplicates()

    merged = backend_keys.merge(
        public_keys,
        on="allele_key",
        how="outer",
        suffixes=("_backend", "_public"),
        indicator=True,
    )
    merged["audit_status"] = merged["_merge"].map(
        {"both": "matched_by_allele", "left_only": "backend_only", "right_only": "public_only"}
    )

    comparable_cols = ["Gene_Symbol", "Transcript", "Chr", "Start", "End", "Ref", "Alt"]
    for col in comparable_cols:
        b = f"{col}_backend"
        p = f"{col}_public"
        if b in merged.columns and p in merged.columns:
            merged[f"{col}_matches"] = merged[b].fillna("").astype(str) == merged[p].fillna("").astype(str)

    backend_only = merged.loc[merged["audit_status"] == "backend_only"].copy()
    public_only = merged.loc[merged["audit_status"] == "public_only"].copy()
    mismatched = merged.loc[
        (merged["audit_status"] == "matched_by_allele")
        & (
            (~merged.get("Gene_Symbol_matches", pd.Series(True, index=merged.index)))
            | (~merged.get("Transcript_matches", pd.Series(True, index=merged.index)))
        )
    ].copy()
    summary = pd.DataFrame(
        [
            {"metric": "backend_records", "value": len(backend)},
            {"metric": "public_records", "value": len(public)},
            {"metric": "backend_unique_alleles", "value": backend["allele_key"].nunique()},
            {"metric": "public_unique_alleles", "value": public["allele_key"].nunique()},
            {"metric": "matched_unique_alleles", "value": int((merged["audit_status"] == "matched_by_allele").sum())},
            {"metric": "backend_only_unique_alleles", "value": len(backend_only)},
            {"metric": "public_only_unique_alleles", "value": len(public_only)},
            {"metric": "matched_alleles_gene_or_transcript_mismatch", "value": len(mismatched)},
        ]
    )

    write_excel(
        out_xlsx,
        {
            "summary": summary,
            "allele_audit": merged.drop(columns=["_merge"]),
            "backend_only": backend_only.drop(columns=["_merge"]),
            "public_only": public_only.drop(columns=["_merge"]),
            "gene_tx_mismatches": mismatched.drop(columns=["_merge"]),
        },
    )

