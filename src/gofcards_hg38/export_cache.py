from __future__ import annotations

import gzip
import re
from datetime import date
from pathlib import Path

import pandas as pd

from .io_utils import ensure_parent, read_excel


def _clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() == "nan" else text


def _normalize_hgvsp(value: object) -> str:
    text = _clean(value)
    if not text or text == "-":
        return ""
    text = text.split(":")[-1]
    text = re.sub(r"^p\.", "", text)
    text = text.replace("%3D", "=").replace("Ter", "*")
    return text


def _hgvsp_key(value: object) -> str:
    return _normalize_hgvsp(value).upper()


def _normalize_symbol(value: object) -> str:
    return _clean(value).upper()


def _read_sheet_or_empty(workbook_xlsx: str | Path, sheet_name: str) -> pd.DataFrame:
    try:
        df = read_excel(workbook_xlsx, sheet_name)
        return df.where(pd.notna(df), "")
    except Exception:
        return pd.DataFrame()


def _core_metadata_by_allele(workbook_xlsx: str | Path) -> dict[str, dict[str, str]]:
    core = _read_sheet_or_empty(workbook_xlsx, "normalized_core")
    if core.empty or "allele_key" not in core.columns:
        return {}
    for col in [
        "Transcript",
        "Function",
        "PMID",
        "Pscore",
        "Pathways_proteins_involved",
        "Phenotype",
    ]:
        if col not in core.columns:
            core[col] = ""
    out: dict[str, dict[str, str]] = {}
    for _, row in core.iterrows():
        key = _clean(row.get("allele_key"))
        if not key or key in out:
            continue
        out[key] = {
            "transcript": _clean(row.get("Transcript")),
            "function": _clean(row.get("Function")),
            "pmids": _clean(row.get("PMID")),
            "pscore": _clean(row.get("Pscore")),
            "pathway": _clean(row.get("Pathways_proteins_involved")),
            "disease": _clean(row.get("Phenotype")),
        }
    return out


def export_priva_gof_tsv(
    workbook_xlsx: str | Path,
    out_tsv: str | Path,
    preferred_only: bool = False,
) -> None:
    sheet = "preferred_transcript_table" if preferred_only else "variant_transcript_table"
    df = read_excel(workbook_xlsx, sheet)
    if df.empty:
        raise ValueError(f"{sheet} is empty in {workbook_xlsx}")

    for col in [
        "gofcards_symbol",
        "gofcards_symbol_resolved",
        "vep_symbol",
        "vep_symbol_resolved",
        "source_refseq_transcript",
        "vep_transcript",
        "HGVSc",
        "HGVSp",
        "gofcards_hgvsc_key",
        "gofcards_hgvsp_key",
        "vep_hgvsc_key",
        "vep_hgvsp_key",
        "gofcards_gene_match",
        "gofcards_hgvsc_match",
        "gofcards_hgvsp_match",
        "gofcards_hgvs_match_status",
        "CANONICAL",
        "hg19_chrom",
        "hg19_start",
        "hg19_end",
        "hg19_ref",
        "hg19_alt",
        "hg38_chrom",
        "hg38_start",
        "hg38_end",
        "hg38_ref",
        "hg38_alt",
        "hg38_refalt_status",
        "gofcards_AAChange_refGene",
        "summary_rsID",
        "allele_key",
    ]:
        if col not in df.columns:
            df[col] = ""

    gofcards_symbol = df["gofcards_symbol_resolved"].where(
        df["gofcards_symbol_resolved"].fillna("").astype(str).str.strip() != "",
        df["gofcards_symbol"],
    )
    vep_symbol = df["vep_symbol_resolved"].where(
        df["vep_symbol_resolved"].fillna("").astype(str).str.strip() != "",
        df["vep_symbol"],
    )
    match_symbol = vep_symbol.where(
        vep_symbol.fillna("").astype(str).str.strip() != "",
        gofcards_symbol,
    )
    normalized_hgvsp = df["HGVSp"].map(_normalize_hgvsp)
    metadata = _core_metadata_by_allele(workbook_xlsx)

    def meta_value(allele_key: object, field: str) -> str:
        return metadata.get(_clean(allele_key), {}).get(field, "")

    out = pd.DataFrame(
        {
            "source": "GoFCards",
            "mechanism": "GOF",
            "build": "hg19_and_hg38",
            "symbol": match_symbol.map(_normalize_symbol),
            "hgvsp_key": df["HGVSp"].map(_hgvsp_key),
            "gofcards_symbol": df["gofcards_symbol"].map(_clean),
            "gofcards_symbol_resolved": gofcards_symbol.map(_clean),
            "vep_symbol": df["vep_symbol"].map(_clean),
            "vep_symbol_resolved": vep_symbol.map(_clean),
            "match_symbol": match_symbol.map(_normalize_symbol),
            "source_refseq_transcript": df["source_refseq_transcript"].map(_clean),
            "vep_transcript": df["vep_transcript"].map(_clean),
            "HGVSc": df["HGVSc"].map(_clean),
            "HGVSp": df["HGVSp"].map(_clean),
            "normalized_hgvsp": normalized_hgvsp,
            "gofcards_hgvsc_key": df["gofcards_hgvsc_key"].map(_clean),
            "gofcards_hgvsp_key": df["gofcards_hgvsp_key"].map(_clean),
            "vep_hgvsc_key": df["vep_hgvsc_key"].map(_clean),
            "vep_hgvsp_key": df["vep_hgvsp_key"].map(_clean),
            "gofcards_gene_match": df["gofcards_gene_match"].map(_clean),
            "gofcards_hgvsc_match": df["gofcards_hgvsc_match"].map(_clean),
            "gofcards_hgvsp_match": df["gofcards_hgvsp_match"].map(_clean),
            "gofcards_hgvs_match_status": df["gofcards_hgvs_match_status"].map(_clean),
            "canonical_transcript": df["CANONICAL"].map(_clean),
            "hg19_chrom": df["hg19_chrom"].map(_clean),
            "hg19_start": df["hg19_start"].map(_clean),
            "hg19_end": df["hg19_end"].map(_clean),
            "hg19_ref": df["hg19_ref"].map(_clean),
            "hg19_alt": df["hg19_alt"].map(_clean),
            "chrom": df["hg19_chrom"].map(lambda value: f"chr{_clean(value)}" if _clean(value) and not _clean(value).startswith("chr") else _clean(value)),
            "pos": df["hg19_start"].map(_clean),
            "ref": df["hg19_ref"].map(_clean),
            "alt": df["hg19_alt"].map(_clean),
            "hg38_chrom": df["hg38_chrom"].map(_clean),
            "hg38_start": df["hg38_start"].map(_clean),
            "hg38_end": df["hg38_end"].map(_clean),
            "hg38_ref": df["hg38_ref"].map(_clean),
            "hg38_alt": df["hg38_alt"].map(_clean),
            "hg38_refalt_status": df["hg38_refalt_status"].map(_clean),
            "gofcards_AAChange_refGene": df["gofcards_AAChange_refGene"].map(_clean),
            "gofcards_accession_id": df["summary_rsID"].map(_clean),
            "gofcards_variant_id": df["allele_key"].map(_clean),
            "disease": df["allele_key"].map(lambda value: meta_value(value, "disease")),
            "pmids": df["allele_key"].map(lambda value: meta_value(value, "pmids")),
            "pscore": df["allele_key"].map(lambda value: meta_value(value, "pscore")),
            "function": df["allele_key"].map(lambda value: meta_value(value, "function")),
            "pathway": df["allele_key"].map(lambda value: meta_value(value, "pathway")),
            "transcript": df["allele_key"].map(lambda value: meta_value(value, "transcript")),
            "derived_on": date.today().isoformat(),
            "allele_key": df["allele_key"].map(_clean),
        }
    )
    out = out.loc[(out["match_symbol"] != "") & (out["hgvsp_key"] != "")].drop_duplicates()
    out_path = ensure_parent(out_tsv)
    if str(out_path).endswith(".gz"):
        with gzip.open(out_path, "wt", encoding="utf-8", newline="") as handle:
            out.to_csv(handle, sep="\t", index=False)
    else:
        out.to_csv(out_path, sep="\t", index=False)
