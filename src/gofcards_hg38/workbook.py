from __future__ import annotations

from pathlib import Path

import pandas as pd

from .io_utils import read_excel, write_excel
from .schema import allele_key_series
from .transvar_io import read_transvar_outputs


def _optional_sheet(path: str | Path, sheet_name: str | int = 0) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        return pd.DataFrame()
    try:
        return read_excel(p, sheet_name)
    except Exception:
        return pd.DataFrame()


def _optional_first_available_sheet(path: str | Path, sheet_names: list[str]) -> pd.DataFrame:
    for sheet_name in sheet_names:
        df = _optional_sheet(path, sheet_name)
        if not df.empty:
            return df
    return _optional_sheet(path, 0)


CORE_COLUMNS = [
    "allele_key",
    "Gene_Symbol",
    "Transcript",
    "AAChange_refGene",
    "Chr",
    "Start",
    "End",
    "Ref",
    "Alt",
    "hg38_Chr",
    "hg38_Start",
    "hg38_End",
    "hg38_Ref_for_vep",
    "hg38_Alt_for_vep",
    "hg38_REF_fasta",
    "hg38_refalt_status",
    "hg38_refalt_needs_review",
    "summary_hg38_start",
    "summary_hg38_end",
    "summary_rsID",
]

TRANSCRIPT_TABLE_COLUMNS = [
    "gofcards_symbol",
    "source_refseq_transcript",
    "assembly",
    "vep_symbol",
    "vep_transcript",
    "feature_type",
    "consequence",
    "HGVSc",
    "HGVSp",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "CANONICAL",
    "Existing_variation",
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
    "hg38_fasta_ref",
    "hg38_refalt_status",
    "hg38_refalt_needs_review",
    "gofcards_AAChange_refGene",
    "summary_rsID",
    "allele_key",
    "Uploaded_variation",
]


def _prepare_core(refalt: pd.DataFrame) -> pd.DataFrame:
    if refalt.empty:
        return refalt
    core = refalt.copy()
    if "allele_key" not in core.columns:
        core["allele_key"] = allele_key_series(core)
    for col in CORE_COLUMNS:
        if col not in core.columns:
            core[col] = ""
    if (core["hg38_Start"].astype(str).str.strip() == "").all() and "summary_hg38_start" in core.columns:
        core["hg38_Start"] = core["summary_hg38_start"]
    if (core["hg38_End"].astype(str).str.strip() == "").all() and "summary_hg38_end" in core.columns:
        core["hg38_End"] = core["summary_hg38_end"]
    if (core["hg38_Chr"].astype(str).str.strip() == "").all() and "Chr" in core.columns:
        core["hg38_Chr"] = core["Chr"]
    return core[CORE_COLUMNS].drop_duplicates("allele_key")


def _core_only_transcript_table(core: pd.DataFrame) -> pd.DataFrame:
    if core.empty:
        return pd.DataFrame(columns=TRANSCRIPT_TABLE_COLUMNS)
    out = pd.DataFrame(
        {
            "gofcards_symbol": core["Gene_Symbol"],
            "source_refseq_transcript": core["Transcript"],
            "assembly": "",
            "vep_symbol": "",
            "vep_transcript": "",
            "feature_type": "",
            "consequence": "",
            "HGVSc": "",
            "HGVSp": "",
            "MANE_SELECT": "",
            "MANE_PLUS_CLINICAL": "",
            "CANONICAL": "",
            "Existing_variation": "",
            "hg19_chrom": core["Chr"],
            "hg19_start": core["Start"],
            "hg19_end": core["End"],
            "hg19_ref": core["Ref"],
            "hg19_alt": core["Alt"],
            "hg38_chrom": core["hg38_Chr"],
            "hg38_start": core["hg38_Start"],
            "hg38_end": core["hg38_End"],
            "hg38_ref": core["hg38_Ref_for_vep"],
            "hg38_alt": core["hg38_Alt_for_vep"],
            "hg38_fasta_ref": core["hg38_REF_fasta"],
            "hg38_refalt_status": core["hg38_refalt_status"],
            "hg38_refalt_needs_review": core["hg38_refalt_needs_review"],
            "gofcards_AAChange_refGene": core["AAChange_refGene"],
            "summary_rsID": core["summary_rsID"],
            "allele_key": core["allele_key"],
            "Uploaded_variation": "",
        }
    )
    return out[TRANSCRIPT_TABLE_COLUMNS]


def _build_transcript_table(core: pd.DataFrame, vep_all: pd.DataFrame) -> pd.DataFrame:
    if core.empty:
        return pd.DataFrame(columns=TRANSCRIPT_TABLE_COLUMNS)
    if vep_all.empty or "allele_key" not in vep_all.columns:
        return _core_only_transcript_table(core)

    vep = vep_all.copy()
    for col in [
        "assembly",
        "allele_key",
        "SYMBOL",
        "Feature",
        "Feature_type",
        "Consequence",
        "HGVSc",
        "HGVSp",
        "MANE_SELECT",
        "MANE_PLUS_CLINICAL",
        "CANONICAL",
        "Existing_variation",
        "Uploaded_variation",
    ]:
        if col not in vep.columns:
            vep[col] = ""

    merged = vep.merge(core, on="allele_key", how="left", suffixes=("_vep", "_core"))
    out = pd.DataFrame(
        {
            "gofcards_symbol": merged["Gene_Symbol"],
            "source_refseq_transcript": merged["Transcript"],
            "assembly": merged["assembly"],
            "vep_symbol": merged["SYMBOL"],
            "vep_transcript": merged["Feature"],
            "feature_type": merged["Feature_type"],
            "consequence": merged["Consequence"],
            "HGVSc": merged["HGVSc"],
            "HGVSp": merged["HGVSp"],
            "MANE_SELECT": merged["MANE_SELECT"],
            "MANE_PLUS_CLINICAL": merged["MANE_PLUS_CLINICAL"],
            "CANONICAL": merged["CANONICAL"],
            "Existing_variation": merged["Existing_variation"],
            "hg19_chrom": merged["Chr"],
            "hg19_start": merged["Start"],
            "hg19_end": merged["End"],
            "hg19_ref": merged["Ref"],
            "hg19_alt": merged["Alt"],
            "hg38_chrom": merged["hg38_Chr"],
            "hg38_start": merged["hg38_Start"],
            "hg38_end": merged["hg38_End"],
            "hg38_ref": merged["hg38_Ref_for_vep"],
            "hg38_alt": merged["hg38_Alt_for_vep"],
            "hg38_fasta_ref": merged["hg38_REF_fasta"],
            "hg38_refalt_status": merged["hg38_refalt_status"],
            "hg38_refalt_needs_review": merged["hg38_refalt_needs_review"],
            "gofcards_AAChange_refGene": merged["AAChange_refGene"],
            "summary_rsID": merged["summary_rsID"],
            "allele_key": merged["allele_key"],
            "Uploaded_variation": merged["Uploaded_variation"],
        }
    )
    return out[TRANSCRIPT_TABLE_COLUMNS]


def _preferred_rank(row: pd.Series) -> tuple[int, str]:
    if str(row.get("MANE_SELECT", "") or "").strip():
        return (0, str(row.get("vep_transcript", "")))
    if str(row.get("MANE_PLUS_CLINICAL", "") or "").strip():
        return (1, str(row.get("vep_transcript", "")))
    if str(row.get("CANONICAL", "") or "").upper() == "YES":
        return (2, str(row.get("vep_transcript", "")))
    if str(row.get("HGVSc", "") or "").strip() or str(row.get("HGVSp", "") or "").strip():
        return (3, str(row.get("vep_transcript", "")))
    return (4, str(row.get("vep_transcript", "")))


def _build_preferred_table(transcript_table: pd.DataFrame) -> pd.DataFrame:
    if transcript_table.empty:
        return transcript_table
    ranked = transcript_table.copy()
    ranks = ranked.apply(_preferred_rank, axis=1)
    ranked["_preferred_rank"] = [r[0] for r in ranks]
    ranked["_preferred_transcript_sort"] = [r[1] for r in ranks]
    ranked = ranked.sort_values(["allele_key", "assembly", "_preferred_rank", "_preferred_transcript_sort"])
    preferred = ranked.drop_duplicates(["allele_key", "assembly"], keep="first")
    return preferred.drop(columns=["_preferred_rank", "_preferred_transcript_sort"])


def build_workbook(
    refalt_xlsx: str | Path,
    audit_xlsx: str | Path,
    vep_xlsx: str | Path,
    transvar_dir: str | Path,
    out_xlsx: str | Path,
) -> None:
    refalt = _optional_first_available_sheet(refalt_xlsx, ["refalt_checked", "augmented_records"])
    refalt_summary = _optional_sheet(refalt_xlsx, "summary")
    audit_summary = _optional_sheet(audit_xlsx, "summary")
    audit = _optional_sheet(audit_xlsx, "allele_audit")
    vep_all = _optional_sheet(vep_xlsx, "vep_all")
    core = _prepare_core(refalt)
    transcript_table = _build_transcript_table(core, vep_all)
    preferred_table = _build_preferred_table(transcript_table)
    sheets: dict[str, pd.DataFrame] = {
        "variant_transcript_table": transcript_table,
        "preferred_transcript_table": preferred_table,
        "normalized_core": refalt,
        "refalt_summary": refalt_summary,
        "public_audit_summary": audit_summary,
        "public_allele_audit": audit,
        "vep_all": vep_all,
    }
    sheets.update(read_transvar_outputs(transvar_dir))
    summary = pd.DataFrame(
        [
            {"artifact": "refalt_xlsx", "path": str(refalt_xlsx), "exists": Path(refalt_xlsx).exists()},
            {"artifact": "audit_xlsx", "path": str(audit_xlsx), "exists": Path(audit_xlsx).exists()},
            {"artifact": "vep_xlsx", "path": str(vep_xlsx), "exists": Path(vep_xlsx).exists()},
            {"artifact": "transvar_dir", "path": str(transvar_dir), "exists": Path(transvar_dir).exists()},
        ]
    )
    sheets = {"run_summary": summary, **sheets}
    write_excel(out_xlsx, sheets)
