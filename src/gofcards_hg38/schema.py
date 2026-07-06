from __future__ import annotations

import re
from typing import Any

import pandas as pd

ALLELE_COLUMNS = ["variant_type", "Chr", "Start", "End", "Ref", "Alt"]

PUBLIC_COLUMN_MAP = {
    "genesymbol": "Gene_Symbol",
    "transcript": "Transcript",
    "chr": "Chr",
    "hg19start": "Start",
    "hg19end": "End",
    "ref": "Ref",
    "alt": "Alt",
    "function": "Function",
    "pathways proteins involved": "Pathways_proteins_involved",
    "disorder involved": "Phenotype",
    "pmid": "PMID",
    "animal model": "Animal_model",
    "cell model": "Cell_model",
    "pscore": "Pscore",
}


def normalize_chrom(value: Any) -> str:
    text = "" if value is None or pd.isna(value) else str(value).strip()
    text = re.sub(r"^chr", "", text, flags=re.IGNORECASE)
    return text


def normalize_allele(value: Any) -> str:
    text = "" if value is None or pd.isna(value) else str(value).strip()
    if text in {".", "-", "nan", "None"}:
        return ""
    return text.upper()


def normalize_int_text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    text = str(value).strip()
    if text.endswith(".0"):
        return text[:-2]
    return text


def add_variant_type(df: pd.DataFrame, variant_type: str | None = None) -> pd.DataFrame:
    out = df.copy()
    if "variant_type" not in out.columns:
        out.insert(0, "variant_type", variant_type or "")
    elif variant_type:
        out["variant_type"] = out["variant_type"].fillna("").replace("", variant_type)
    return out


def normalize_backend_frame(df: pd.DataFrame, variant_type: str | None = None) -> pd.DataFrame:
    out = add_variant_type(df, variant_type)
    for col in ["Chr", "Start", "End", "Ref", "Alt"]:
        if col not in out.columns:
            out[col] = ""
    out["variant_type"] = out["variant_type"].astype(str).str.strip()
    out["Chr"] = out["Chr"].map(normalize_chrom)
    out["Start"] = out["Start"].map(normalize_int_text)
    out["End"] = out["End"].map(normalize_int_text)
    out["Ref"] = out["Ref"].map(normalize_allele)
    out["Alt"] = out["Alt"].map(normalize_allele)
    return out


def normalize_public_frame(df: pd.DataFrame) -> pd.DataFrame:
    rename = {}
    for col in df.columns:
        key = str(col).strip().lower()
        if key in PUBLIC_COLUMN_MAP:
            rename[col] = PUBLIC_COLUMN_MAP[key]
    out = df.rename(columns=rename).copy()
    out = normalize_backend_frame(out, None)
    out["variant_type"] = out.apply(
        lambda row: "SNV"
        if len(str(row.get("Ref", ""))) == 1 and len(str(row.get("Alt", ""))) == 1
        else "Indel",
        axis=1,
    )
    return out


def allele_key_series(df: pd.DataFrame) -> pd.Series:
    tmp = normalize_backend_frame(df)
    return tmp[ALLELE_COLUMNS].astype(str).agg("|".join, axis=1)

