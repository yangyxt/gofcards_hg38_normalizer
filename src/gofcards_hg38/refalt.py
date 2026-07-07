from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
from pyfaidx import Fasta

from .io_utils import read_excel, write_excel
from .schema import allele_key_series, normalize_backend_frame


def _read_augmented(path: str | Path) -> pd.DataFrame:
    try:
        return read_excel(path, "augmented_records")
    except ValueError:
        return read_excel(path, 0)


def _chrom_name(fasta: Fasta, chrom: str) -> str | None:
    candidates = [chrom, f"chr{chrom}"]
    if chrom.upper() == "M":
        candidates.extend(["MT", "chrM", "chrMT"])
    if chrom.upper() == "MT":
        candidates.extend(["M", "chrM", "chrMT"])
    for candidate in candidates:
        if candidate in fasta:
            return candidate
    return None


def _to_int(value: Any) -> int | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def _fasta_seq(fasta: Fasta, chrom: str, start_1based: int, length: int) -> str:
    name = _chrom_name(fasta, chrom)
    if name is None:
        raise KeyError(f"Chromosome {chrom!r} not found in FASTA")
    if start_1based < 1:
        raise ValueError(f"Invalid 1-based FASTA coordinate: {start_1based}")
    start0 = start_1based - 1
    end0 = start0 + length
    return str(fasta[name][start0:end0]).upper()


def _vcf_fields(
    prefix: str,
    *,
    chrom: str,
    start: int | None,
    ref: str,
    alt: str,
    fasta: Fasta | None,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        f"{prefix}_REF_fasta": "",
        f"{prefix}_VCF_Pos": "",
        f"{prefix}_Ref_for_vep": "",
        f"{prefix}_Alt_for_vep": "",
        f"{prefix}_vcf_status": "",
        f"{prefix}_vcf_needs_review": "Y",
    }
    if start is None:
        result[f"{prefix}_vcf_status"] = "missing_coordinate"
        return result

    if ref and alt:
        result[f"{prefix}_VCF_Pos"] = start
        result[f"{prefix}_Ref_for_vep"] = ref
        result[f"{prefix}_Alt_for_vep"] = alt
        result[f"{prefix}_vcf_status"] = "raw_ref_alt"
        result[f"{prefix}_vcf_needs_review"] = "N"
        if fasta is not None:
            try:
                fasta_ref = _fasta_seq(fasta, chrom, start, max(len(ref), 1))
                result[f"{prefix}_REF_fasta"] = fasta_ref
                if fasta_ref != ref:
                    result[f"{prefix}_vcf_status"] = "raw_ref_alt_ref_mismatch"
                    result[f"{prefix}_vcf_needs_review"] = "Y"
            except Exception as exc:
                result[f"{prefix}_vcf_status"] = f"fasta_lookup_error:{exc}"
        return result

    if fasta is None:
        result[f"{prefix}_vcf_status"] = "missing_ref_or_alt_no_fasta"
        return result

    try:
        if ref and not alt:
            anchor_pos = start - 1
            anchor = _fasta_seq(fasta, chrom, anchor_pos, 1)
            fasta_ref = _fasta_seq(fasta, chrom, start, max(len(ref), 1))
            result[f"{prefix}_REF_fasta"] = fasta_ref
            result[f"{prefix}_VCF_Pos"] = anchor_pos
            result[f"{prefix}_Ref_for_vep"] = anchor + (ref if fasta_ref == ref else fasta_ref)
            result[f"{prefix}_Alt_for_vep"] = anchor
            result[f"{prefix}_vcf_status"] = (
                "deletion_padded_ref_match" if fasta_ref == ref else "deletion_padded_ref_mismatch_using_fasta"
            )
            result[f"{prefix}_vcf_needs_review"] = "N" if fasta_ref == ref else "Y"
            return result

        if alt and not ref:
            anchor = _fasta_seq(fasta, chrom, start, 1)
            result[f"{prefix}_REF_fasta"] = anchor
            result[f"{prefix}_VCF_Pos"] = start
            result[f"{prefix}_Ref_for_vep"] = anchor
            result[f"{prefix}_Alt_for_vep"] = anchor + alt
            result[f"{prefix}_vcf_status"] = "insertion_padded"
            result[f"{prefix}_vcf_needs_review"] = "N"
            return result
    except Exception as exc:
        result[f"{prefix}_vcf_status"] = f"fasta_lookup_error:{exc}"
        return result

    result[f"{prefix}_vcf_status"] = "missing_ref_and_alt"
    return result


def _source_vcf_row(row: pd.Series, fasta: Fasta | None) -> dict[str, Any]:
    chrom = str(row.get("Chr", "")).replace("chr", "")
    start = _to_int(row.get("Start"))
    ref = str(row.get("Ref", "") or "").upper()
    alt = str(row.get("Alt", "") or "").upper()
    return _vcf_fields("hg19", chrom=chrom, start=start, ref=ref, alt=alt, fasta=fasta)


def _check_row(row: pd.Series, fasta: Fasta) -> dict[str, Any]:
    chrom = str(row.get("Chr", "")).replace("chr", "")
    start = _to_int(row.get("summary_hg38_start"))
    end = _to_int(row.get("summary_hg38_end"))
    ref = str(row.get("Ref", "") or "").upper()
    alt = str(row.get("Alt", "") or "").upper()
    variant_type = str(row.get("variant_type", "") or "")

    result: dict[str, Any] = {
        "hg38_Chr": chrom,
        "hg38_Start": start or "",
        "hg38_End": end or "",
        "hg38_REF_fasta": "",
        "hg38_VCF_Pos": "",
        "hg38_Ref_for_vep": "",
        "hg38_Alt_for_vep": "",
        "hg38_refalt_status": "",
        "hg38_refalt_needs_review": "Y",
    }
    if start is None or end is None:
        result["hg38_refalt_status"] = "missing_hg38_coordinate"
        return result
    if not ref or not alt:
        vcf = _vcf_fields("hg38", chrom=chrom, start=start, ref=ref, alt=alt, fasta=fasta)
        result.update(
            {
                "hg38_REF_fasta": vcf["hg38_REF_fasta"],
                "hg38_VCF_Pos": vcf["hg38_VCF_Pos"],
                "hg38_Ref_for_vep": vcf["hg38_Ref_for_vep"],
                "hg38_Alt_for_vep": vcf["hg38_Alt_for_vep"],
                "hg38_refalt_status": vcf["hg38_vcf_status"],
                "hg38_refalt_needs_review": vcf["hg38_vcf_needs_review"],
            }
        )
        return result
    try:
        fasta_ref = _fasta_seq(fasta, chrom, start, max(len(ref), 1))
    except Exception as exc:
        result["hg38_refalt_status"] = f"fasta_lookup_error:{exc}"
        return result

    result["hg38_REF_fasta"] = fasta_ref
    is_snv = variant_type.upper() == "SNV" or (len(ref) == 1 and len(alt) == 1)
    if is_snv:
        if len(ref) != 1 or len(alt) != 1:
            result["hg38_refalt_status"] = "variant_type_snv_but_allele_length_not_1"
            return result
        result["hg38_VCF_Pos"] = start
        result["hg38_Ref_for_vep"] = fasta_ref[:1]
        if fasta_ref[:1] == ref:
            result["hg38_Alt_for_vep"] = alt
            result["hg38_refalt_status"] = "snv_ref_match"
            result["hg38_refalt_needs_review"] = "N"
        else:
            result["hg38_refalt_status"] = "snv_ref_mismatch_needs_review"
        return result

    if fasta_ref == ref:
        result["hg38_VCF_Pos"] = start
        result["hg38_Ref_for_vep"] = ref
        result["hg38_Alt_for_vep"] = alt
        result["hg38_refalt_status"] = "indel_ref_match_not_left_normalized"
        result["hg38_refalt_needs_review"] = "N"
    else:
        result["hg38_refalt_status"] = "indel_ref_mismatch_needs_normalization"
    return result


def validate_refalt(
    input_xlsx: str | Path,
    hg38_fasta: str | Path,
    out_xlsx: str | Path,
    hg19_fasta: str | Path | None = None,
) -> None:
    df = normalize_backend_frame(_read_augmented(input_xlsx))
    if "allele_key" not in df.columns:
        df["allele_key"] = allele_key_series(df)
    source_fasta = Fasta(str(hg19_fasta), rebuild=False) if hg19_fasta else None
    fasta = Fasta(str(hg38_fasta), rebuild=False)
    source_checks = pd.DataFrame([_source_vcf_row(row, source_fasta) for _, row in df.iterrows()])
    checks = pd.DataFrame([_check_row(row, fasta) for _, row in df.iterrows()])
    out = pd.concat([df.reset_index(drop=True), source_checks, checks], axis=1)
    summary = (
        out["hg38_refalt_status"]
        .value_counts(dropna=False)
        .rename_axis("hg38_refalt_status")
        .reset_index(name="records")
    )
    write_excel(out_xlsx, {"refalt_checked": out, "summary": summary})
