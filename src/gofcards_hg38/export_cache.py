from __future__ import annotations

import gzip
import re
from datetime import date
from pathlib import Path

import pandas as pd

from .io_utils import ensure_parent, read_excel

HGVS_MATCH_STATUSES = {
    "both_cdna_protein_match",
    "protein_only_match",
    "cdna_only_match",
}
GENOMIC_ONLY_STATUSES = {
    "no_parseable_gofcards_hgvs",
    "same_gene_no_hgvs_match",
    "no_vep_hgvs",
    "no_same_gene_vep_row",
    "other_no_match",
}
INVALID_SYMBOLS = {"", "-", ".", "NA", "N/A", "UNKNOWN"}


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


def _valid_symbol_mask(values: pd.Series) -> pd.Series:
    return ~values.map(_normalize_symbol).isin(INVALID_SYMBOLS)


def _first_nonblank(*values: object) -> str:
    for value in values:
        text = _clean(value)
        if text:
            return text
    return ""


def _genomic_key(chrom: object, pos: object, ref: object, alt: object) -> str:
    chrom_text = _clean(chrom)
    pos_text = _clean(pos)
    ref_text = _clean(ref)
    alt_text = _clean(alt)
    if not chrom_text or not pos_text or not ref_text or not alt_text:
        return ""
    return "|".join([chrom_text, pos_text, ref_text, alt_text])


def _read_sheet_or_empty(workbook_xlsx: str | Path, sheet_name: str) -> pd.DataFrame:
    try:
        df = read_excel(workbook_xlsx, sheet_name)
        return df.where(pd.notna(df), "")
    except Exception:
        return pd.DataFrame()


def _read_normalized_vcf_by_allele(workbook_xlsx: str | Path, assembly: str) -> dict[str, dict[str, str]]:
    workbook_path = Path(workbook_xlsx)
    key_path = workbook_path.parent / "vep_inputs" / "gofcards_vep_input_key.xlsx"
    vcf_path = workbook_path.parent / "vep_inputs" / f"gofcards.{assembly}.norm.vcf"
    if not key_path.exists() or not vcf_path.exists():
        return {}
    key_df = read_excel(key_path, 0).where(lambda df: pd.notna(df), "")
    required = {"assembly", "vcf_id", "allele_key"}
    if not required.issubset(set(key_df.columns)):
        return {}
    id_to_key = {
        _clean(row.get("vcf_id")): _clean(row.get("allele_key"))
        for _, row in key_df[key_df["assembly"].map(_clean).eq(assembly)].iterrows()
    }
    out: dict[str, dict[str, str]] = {}
    with open(vcf_path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, record_id, ref, alt = fields[:5]
            allele_key = id_to_key.get(record_id)
            if not allele_key or allele_key in out:
                continue
            out[allele_key] = {
                "chrom": _clean(chrom),
                "pos": _clean(pos),
                "ref": _clean(ref),
                "alt": _clean(alt),
            }
    return out


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
        "hg19_vcf_pos",
        "hg19_vcf_ref",
        "hg19_vcf_alt",
        "hg19_vcf_status",
        "hg38_chrom",
        "hg38_start",
        "hg38_end",
        "hg38_ref",
        "hg38_alt",
        "hg38_vcf_pos",
        "hg38_vcf_ref",
        "hg38_vcf_alt",
        "hg38_refalt_status",
        "gofcards_AAChange_refGene",
        "summary_rsID",
        "allele_key",
    ]:
        if col not in df.columns:
            df[col] = ""

    df = df.where(pd.notna(df), "")
    df["gofcards_hgvs_match_status"] = df["gofcards_hgvs_match_status"].map(_clean)
    keep_status = HGVS_MATCH_STATUSES | GENOMIC_ONLY_STATUSES
    df = df[df["gofcards_hgvs_match_status"].isin(keep_status)].copy()
    if df.empty:
        raise ValueError(f"{sheet} has no GoFCards rows usable for exact variant matching")

    gofcards_symbol = df["gofcards_symbol_resolved"].where(
        _valid_symbol_mask(df["gofcards_symbol_resolved"]),
        df["gofcards_symbol"],
    )
    vep_symbol = df["vep_symbol_resolved"].where(
        _valid_symbol_mask(df["vep_symbol_resolved"]),
        df["vep_symbol"],
    )
    match_symbol = vep_symbol.where(_valid_symbol_mask(vep_symbol), gofcards_symbol)
    is_genomic_only = ~df["gofcards_hgvs_match_status"].isin(HGVS_MATCH_STATUSES)
    hgvsc = df["HGVSc"].where(~is_genomic_only, "")
    hgvsp = df["HGVSp"].where(~is_genomic_only, "")
    hgvsp_key = hgvsp.map(_hgvsp_key)
    vep_assembly = df["assembly"].where(~is_genomic_only, "")
    vep_transcript = df["vep_transcript"].where(~is_genomic_only, "")
    feature_type = df["feature_type"].where(~is_genomic_only, "")
    consequence = df["consequence"].where(~is_genomic_only, "")
    canonical_transcript = df["CANONICAL"].where(~is_genomic_only, "")

    hg19_pos = [
        _first_nonblank(vcf_pos, start)
        for vcf_pos, start in zip(df["hg19_vcf_pos"], df["hg19_start"])
    ]
    hg19_ref = [
        _first_nonblank(vcf_ref, raw_ref)
        for vcf_ref, raw_ref in zip(df["hg19_vcf_ref"], df["hg19_ref"])
    ]
    hg19_alt = [
        _first_nonblank(vcf_alt, raw_alt)
        for vcf_alt, raw_alt in zip(df["hg19_vcf_alt"], df["hg19_alt"])
    ]
    hg38_pos = [
        _first_nonblank(vcf_pos, start)
        for vcf_pos, start in zip(df["hg38_vcf_pos"], df["hg38_start"])
    ]
    hg38_ref = [
        _first_nonblank(vcf_ref, raw_ref)
        for vcf_ref, raw_ref in zip(df["hg38_vcf_ref"], df["hg38_ref"])
    ]
    hg38_alt = [
        _first_nonblank(vcf_alt, raw_alt)
        for vcf_alt, raw_alt in zip(df["hg38_vcf_alt"], df["hg38_alt"])
    ]

    metadata = _core_metadata_by_allele(workbook_xlsx)
    normalized_vcf = {
        "hg19": _read_normalized_vcf_by_allele(workbook_xlsx, "hg19"),
        "hg38": _read_normalized_vcf_by_allele(workbook_xlsx, "hg38"),
    }

    def meta_value(allele_key: object, field: str) -> str:
        return metadata.get(_clean(allele_key), {}).get(field, "")

    out = pd.DataFrame(
        {
            "source": "GoFCards",
            "mechanism": "GOF",
            "build": "hg19_and_hg38",
            "HGNC_Symbol": match_symbol.map(_normalize_symbol),
            "VEP_assembly": vep_assembly.map(_clean),
            "VEP_transcript": vep_transcript.map(_clean),
            "feature_type": feature_type.map(_clean),
            "consequence": consequence.map(_clean),
            "HGVSc": hgvsc.map(_clean),
            "HGVSp": hgvsp.map(_clean),
            "hgvsp_key": hgvsp_key,
            "match_status": df["gofcards_hgvs_match_status"].map(_clean),
            "raw_GoFCards_HGVS": df["gofcards_AAChange_refGene"].map(_clean),
            "GoFCards_transcript": df["source_refseq_transcript"].map(_clean),
            "canonical_transcript": canonical_transcript.map(_clean),
            "hg19_chrom": df["hg19_chrom"].map(_clean),
            "hg19_pos": pd.Series(hg19_pos, index=df.index).map(_clean),
            "hg19_ref": pd.Series(hg19_ref, index=df.index).map(_clean),
            "hg19_alt": pd.Series(hg19_alt, index=df.index).map(_clean),
            "hg19_vcf_status": df["hg19_vcf_status"].map(_clean),
            "hg38_chrom": df["hg38_chrom"].map(_clean),
            "hg38_pos": pd.Series(hg38_pos, index=df.index).map(_clean),
            "hg38_ref": pd.Series(hg38_ref, index=df.index).map(_clean),
            "hg38_alt": pd.Series(hg38_alt, index=df.index).map(_clean),
            "hg38_refalt_status": df["hg38_refalt_status"].map(_clean),
            "gofcards_accession_id": df["summary_rsID"].map(_clean),
            "gofcards_variant_id": df["allele_key"].map(_clean),
            "disease": df["allele_key"].map(lambda value: meta_value(value, "disease")),
            "pmids": df["allele_key"].map(lambda value: meta_value(value, "pmids")),
            "pscore": df["allele_key"].map(lambda value: meta_value(value, "pscore")),
            "function": df["allele_key"].map(lambda value: meta_value(value, "function")),
            "pathway": df["allele_key"].map(lambda value: meta_value(value, "pathway")),
            "derived_on": date.today().isoformat(),
            "allele_key": df["allele_key"].map(_clean),
        }
    )
    out["hg19_genomic_key"] = [
        _genomic_key(chrom, pos, ref, alt)
        for chrom, pos, ref, alt in zip(out["hg19_chrom"], out["hg19_pos"], out["hg19_ref"], out["hg19_alt"])
    ]
    for assembly in ("hg19", "hg38"):
        norm_rows = normalized_vcf[assembly]
        out[f"{assembly}_vcf_pos"] = [
            norm_rows.get(allele_key, {}).get("pos", raw_pos)
            for allele_key, raw_pos in zip(out["allele_key"], out[f"{assembly}_pos"])
        ]
        out[f"{assembly}_vcf_ref"] = [
            norm_rows.get(allele_key, {}).get("ref", raw_ref)
            for allele_key, raw_ref in zip(out["allele_key"], out[f"{assembly}_ref"])
        ]
        out[f"{assembly}_vcf_alt"] = [
            norm_rows.get(allele_key, {}).get("alt", raw_alt)
            for allele_key, raw_alt in zip(out["allele_key"], out[f"{assembly}_alt"])
        ]
    out["hg19_vcf_key"] = [
        _genomic_key(chrom, pos, ref, alt)
        for chrom, pos, ref, alt in zip(out["hg19_chrom"], out["hg19_vcf_pos"], out["hg19_vcf_ref"], out["hg19_vcf_alt"])
    ]
    out["hg38_genomic_key"] = [
        _genomic_key(chrom, pos, ref, alt)
        for chrom, pos, ref, alt in zip(out["hg38_chrom"], out["hg38_pos"], out["hg38_ref"], out["hg38_alt"])
    ]
    out["hg38_vcf_key"] = [
        _genomic_key(chrom, pos, ref, alt)
        for chrom, pos, ref, alt in zip(out["hg38_chrom"], out["hg38_vcf_pos"], out["hg38_vcf_ref"], out["hg38_vcf_alt"])
    ]
    key_types: list[str] = []
    for _, row in out.iterrows():
        keys = []
        if row["hgvsp_key"]:
            keys.append("hgvsp")
        if row["hg19_genomic_key"]:
            keys.append("hg19_genomic")
        if row["hg19_vcf_key"]:
            keys.append("hg19_vcf")
        if row["hg38_genomic_key"]:
            keys.append("hg38_genomic")
        if row["hg38_vcf_key"]:
            keys.append("hg38_vcf")
        key_types.append(";".join(keys))
    out["match_key_types"] = key_types
    required = [
        "HGNC_Symbol",
        "hg19_chrom",
        "hg19_pos",
        "hg19_ref",
        "hg19_alt",
        "hg38_chrom",
        "hg38_pos",
        "hg38_ref",
        "hg38_alt",
        "hg19_genomic_key",
        "hg38_genomic_key",
    ]
    missing_required = out[required].astype(str).apply(lambda col: col.str.strip().eq("")).any(axis=1)
    valid_symbol = ~out["HGNC_Symbol"].astype(str).str.strip().str.upper().isin(INVALID_SYMBOLS)
    out = out.loc[valid_symbol & ~missing_required & (out["match_key_types"] != "")].drop_duplicates()
    out_path = ensure_parent(out_tsv)
    if str(out_path).endswith(".gz"):
        with gzip.open(out_path, "wt", encoding="utf-8", newline="") as handle:
            out.to_csv(handle, sep="\t", index=False)
    else:
        out.to_csv(out_path, sep="\t", index=False)
