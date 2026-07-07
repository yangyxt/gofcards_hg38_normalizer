from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .hgnc import resolve_symbol
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
    "summary_AAChange_refGene",
]

TRANSCRIPT_TABLE_COLUMNS = [
    "gofcards_symbol",
    "gofcards_symbol_resolved",
    "source_refseq_transcript",
    "assembly",
    "vep_symbol",
    "vep_symbol_resolved",
    "vep_transcript",
    "feature_type",
    "consequence",
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

AA3_TO_1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
    "Sec": "U",
    "Pyl": "O",
    "Xaa": "X",
}
AA3_PATTERN = re.compile("|".join(map(re.escape, sorted(AA3_TO_1, key=len, reverse=True))))


def _clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() == "nan" else text


def _normalize_hgvsc(value: object) -> str:
    text = _clean(value)
    if not text or text in {"-", "."}:
        return ""
    text = text.split(":")[-1]
    match = re.fullmatch(r"c\.([ACGTN])(\d+)([ACGTN])", text, flags=re.IGNORECASE)
    if match:
        ref, pos, alt = match.groups()
        text = f"c.{pos}{ref.upper()}>{alt.upper()}"
    return text.upper().replace(" ", "")


def _normalize_hgvsp(value: object) -> str:
    text = _clean(value).replace("%3D", "=")
    if not text or text in {"-", "."}:
        return ""
    text = text.split(":")[-1]
    text = re.sub(r"^p\.", "", text)
    text = AA3_PATTERN.sub(lambda match: AA3_TO_1[match.group(0)], text)
    return text.replace("Ter", "*").replace("Stop", "*").upper().replace(" ", "")


def _parse_aachange(value: object) -> list[dict[str, str]]:
    entries: list[dict[str, str]] = []
    for item in _clean(value).split(","):
        parts = [part.strip() for part in item.split(":") if part.strip()]
        if len(parts) < 2:
            continue
        gene = parts[0]
        cdna = next((part for part in parts if part.startswith("c.")), "")
        protein = next((part for part in parts if part.startswith("p.")), "")
        entries.append(
            {
                "gene": resolve_symbol(gene),
                "hgvsc": _normalize_hgvsc(cdna),
                "hgvsp": _normalize_hgvsp(protein),
            }
        )
    return entries


def _join_keys(values: list[str]) -> str:
    return ";".join(sorted({value for value in values if value}))


def _has_parseable_aachange(value: object) -> bool:
    return any(entry["hgvsc"] or entry["hgvsp"] for entry in _parse_aachange(value))


def _select_aachange(raw_value: object, summary_value: object) -> str:
    raw = _clean(raw_value)
    summary = _clean(summary_value)
    if _has_parseable_aachange(raw):
        return raw
    if _has_parseable_aachange(summary):
        return summary
    return raw or summary


def _annotate_hgvs_matches(table: pd.DataFrame) -> pd.DataFrame:
    if table.empty:
        for col in TRANSCRIPT_TABLE_COLUMNS:
            if col not in table.columns:
                table[col] = ""
        return table[TRANSCRIPT_TABLE_COLUMNS]

    out = table.copy()
    annotations: list[dict[str, str]] = []
    for _, row in out.iterrows():
        gofcards_symbol_resolved = resolve_symbol(row.get("gofcards_symbol"))
        vep_symbol_resolved = resolve_symbol(row.get("vep_symbol"))
        source_entries = _parse_aachange(row.get("gofcards_AAChange_refGene"))
        vep_hgvsc_key = _normalize_hgvsc(row.get("HGVSc"))
        vep_hgvsp_key = _normalize_hgvsp(row.get("HGVSp"))
        source_hgvsc_keys = [entry["hgvsc"] for entry in source_entries]
        source_hgvsp_keys = [entry["hgvsp"] for entry in source_entries]

        gene_match = any(entry["gene"] and entry["gene"] == vep_symbol_resolved for entry in source_entries)
        hgvsc_match = any(
            gene_match
            and entry["gene"] == vep_symbol_resolved
            and entry["hgvsc"]
            and entry["hgvsc"] == vep_hgvsc_key
            for entry in source_entries
        )
        hgvsp_match = any(
            gene_match
            and entry["gene"] == vep_symbol_resolved
            and entry["hgvsp"]
            and entry["hgvsp"] == vep_hgvsp_key
            for entry in source_entries
        )
        both_match = any(
            entry["gene"] == vep_symbol_resolved
            and entry["hgvsc"]
            and entry["hgvsp"]
            and entry["hgvsc"] == vep_hgvsc_key
            and entry["hgvsp"] == vep_hgvsp_key
            for entry in source_entries
        )
        has_source_hgvs = any(entry["hgvsc"] or entry["hgvsp"] for entry in source_entries)
        has_vep_hgvs = bool(vep_hgvsc_key or vep_hgvsp_key)

        if both_match:
            status = "both_cdna_protein_match"
        elif hgvsp_match:
            status = "protein_only_match"
        elif hgvsc_match:
            status = "cdna_only_match"
        elif gene_match and has_source_hgvs and has_vep_hgvs:
            status = "same_gene_no_hgvs_match"
        elif not has_source_hgvs:
            status = "no_parseable_gofcards_hgvs"
        elif not has_vep_hgvs:
            status = "no_vep_hgvs"
        elif not gene_match:
            status = "no_same_gene_vep_row"
        else:
            status = "other_no_match"

        annotations.append(
            {
                "gofcards_symbol_resolved": gofcards_symbol_resolved,
                "vep_symbol_resolved": vep_symbol_resolved,
                "gofcards_hgvsc_key": _join_keys(source_hgvsc_keys),
                "gofcards_hgvsp_key": _join_keys(source_hgvsp_keys),
                "vep_hgvsc_key": vep_hgvsc_key,
                "vep_hgvsp_key": vep_hgvsp_key,
                "gofcards_gene_match": "Y" if gene_match else "N",
                "gofcards_hgvsc_match": "Y" if hgvsc_match else "N",
                "gofcards_hgvsp_match": "Y" if hgvsp_match else "N",
                "gofcards_hgvs_match_status": status,
            }
        )

    annotation_df = pd.DataFrame(annotations, index=out.index)
    for col in annotation_df.columns:
        out[col] = annotation_df[col]
    return out[TRANSCRIPT_TABLE_COLUMNS]


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
            "gofcards_symbol_resolved": "",
            "source_refseq_transcript": core["Transcript"],
            "assembly": "",
            "vep_symbol": "",
            "vep_symbol_resolved": "",
            "vep_transcript": "",
            "feature_type": "",
            "consequence": "",
            "HGVSc": "",
            "HGVSp": "",
            "gofcards_hgvsc_key": "",
            "gofcards_hgvsp_key": "",
            "vep_hgvsc_key": "",
            "vep_hgvsp_key": "",
            "gofcards_gene_match": "",
            "gofcards_hgvsc_match": "",
            "gofcards_hgvsp_match": "",
            "gofcards_hgvs_match_status": "",
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
            "gofcards_AAChange_refGene": [
                _select_aachange(raw, summary)
                for raw, summary in zip(core["AAChange_refGene"], core["summary_AAChange_refGene"])
            ],
            "summary_rsID": core["summary_rsID"],
            "allele_key": core["allele_key"],
            "Uploaded_variation": "",
        }
    )
    return _annotate_hgvs_matches(out)


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
            "gofcards_symbol_resolved": "",
            "source_refseq_transcript": merged["Transcript"],
            "assembly": merged["assembly"],
            "vep_symbol": merged["SYMBOL"],
            "vep_symbol_resolved": "",
            "vep_transcript": merged["Feature"],
            "feature_type": merged["Feature_type"],
            "consequence": merged["Consequence"],
            "HGVSc": merged["HGVSc"],
            "HGVSp": merged["HGVSp"],
            "gofcards_hgvsc_key": "",
            "gofcards_hgvsp_key": "",
            "vep_hgvsc_key": "",
            "vep_hgvsp_key": "",
            "gofcards_gene_match": "",
            "gofcards_hgvsc_match": "",
            "gofcards_hgvsp_match": "",
            "gofcards_hgvs_match_status": "",
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
            "gofcards_AAChange_refGene": [
                _select_aachange(raw, summary)
                for raw, summary in zip(merged["AAChange_refGene"], merged["summary_AAChange_refGene"])
            ],
            "summary_rsID": merged["summary_rsID"],
            "allele_key": merged["allele_key"],
            "Uploaded_variation": merged["Uploaded_variation"],
        }
    )
    out = out[TRANSCRIPT_TABLE_COLUMNS]

    # Preserve alleles that VEP skipped, e.g. source records with incomplete
    # indel alleles. The final table is also a coordinate/ref-alt audit table,
    # so every GoFCards allele must remain visible even without transcript HGVS.
    vep_keys = set(out["allele_key"].dropna().astype(str))
    missing_core = core.loc[~core["allele_key"].astype(str).isin(vep_keys)].copy()
    if not missing_core.empty:
        out = pd.concat([out, _core_only_transcript_table(missing_core)], ignore_index=True)
    return _annotate_hgvs_matches(out)


def _preferred_rank(row: pd.Series) -> tuple[int, int, str]:
    match_rank = {
        "both_cdna_protein_match": 0,
        "protein_only_match": 1,
        "cdna_only_match": 2,
        "same_gene_no_hgvs_match": 3,
        "no_parseable_gofcards_hgvs": 4,
        "no_vep_hgvs": 5,
        "no_same_gene_vep_row": 6,
        "other_no_match": 7,
    }.get(str(row.get("gofcards_hgvs_match_status", "")), 8)
    if str(row.get("MANE_SELECT", "") or "").strip():
        mane_rank = 0
    elif str(row.get("MANE_PLUS_CLINICAL", "") or "").strip():
        mane_rank = 1
    elif str(row.get("CANONICAL", "") or "").upper() == "YES":
        mane_rank = 2
    elif str(row.get("HGVSc", "") or "").strip() or str(row.get("HGVSp", "") or "").strip():
        mane_rank = 3
    else:
        mane_rank = 4
    return (match_rank, mane_rank, str(row.get("vep_transcript", "")))


def _build_preferred_table(transcript_table: pd.DataFrame) -> pd.DataFrame:
    if transcript_table.empty:
        return transcript_table
    ranked = transcript_table.copy()
    ranks = ranked.apply(_preferred_rank, axis=1)
    ranked["_preferred_rank"] = [r[0] for r in ranks]
    ranked["_preferred_mane_rank"] = [r[1] for r in ranks]
    ranked["_preferred_transcript_sort"] = [r[2] for r in ranks]
    ranked = ranked.sort_values(
        ["allele_key", "assembly", "_preferred_rank", "_preferred_mane_rank", "_preferred_transcript_sort"]
    )
    preferred = ranked.drop_duplicates(["allele_key", "assembly"], keep="first")
    return preferred.drop(columns=["_preferred_rank", "_preferred_mane_rank", "_preferred_transcript_sort"])


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
