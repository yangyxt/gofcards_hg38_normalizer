from __future__ import annotations

from io import StringIO
from pathlib import Path
from urllib.parse import quote

import pandas as pd

from .io_utils import ensure_parent, read_excel, write_excel
from .schema import allele_key_series, normalize_backend_frame


def _read_refalt(path: str | Path) -> pd.DataFrame:
    try:
        return read_excel(path, "refalt_checked")
    except ValueError:
        return read_excel(path, 0)


def _vcf_escape(value: object) -> str:
    text = "" if value is None or pd.isna(value) else str(value)
    return quote(text.replace(";", ","), safe="._:-")


def _vcf_id(prefix: str, key: str, index: int) -> str:
    safe = quote(key, safe="").replace("%", "_")
    return f"{prefix}_{index}_{safe[:80]}"


def _write_vcf(rows: list[dict], path: str | Path) -> None:
    out = ensure_parent(path)
    with out.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##INFO=<ID=ALLELE_KEY,Number=1,Type=String,Description=\"GoFCards allele key\">\n")
        handle.write("##INFO=<ID=GENE,Number=1,Type=String,Description=\"GoFCards gene symbol\">\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in rows:
            info = f"ALLELE_KEY={_vcf_escape(row['allele_key'])};GENE={_vcf_escape(row.get('Gene_Symbol', ''))}"
            handle.write(
                "\t".join(
                    [
                        str(row["chrom"]),
                        str(row["pos"]),
                        str(row["id"]),
                        str(row["ref"]),
                        str(row["alt"]),
                        ".",
                        "PASS",
                        info,
                    ]
                )
                + "\n"
            )


def write_vep_inputs(input_xlsx: str | Path, out_dir: str | Path) -> None:
    df = normalize_backend_frame(_read_refalt(input_xlsx))
    if "allele_key" not in df.columns:
        df["allele_key"] = allele_key_series(df)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    hg19_rows: list[dict] = []
    hg38_rows: list[dict] = []
    key_rows: list[dict] = []
    for idx, row in df.iterrows():
        key = str(row["allele_key"])
        hg19_id = _vcf_id("gofcards_hg19", key, idx + 1)
        hg19_rows.append(
            {
                "chrom": row.get("Chr"),
                "pos": row.get("Start"),
                "id": hg19_id,
                "ref": row.get("Ref"),
                "alt": row.get("Alt"),
                "allele_key": key,
                "Gene_Symbol": row.get("Gene_Symbol", ""),
            }
        )
        key_rows.append({"assembly": "hg19", "vcf_id": hg19_id, "allele_key": key})

        ref38 = str(row.get("hg38_Ref_for_vep", "") or "")
        alt38 = str(row.get("hg38_Alt_for_vep", "") or "")
        start38 = str(row.get("hg38_Start", "") or "")
        if ref38 and alt38 and start38:
            hg38_id = _vcf_id("gofcards_hg38", key, idx + 1)
            hg38_rows.append(
                {
                    "chrom": row.get("hg38_Chr") or row.get("Chr"),
                    "pos": start38,
                    "id": hg38_id,
                    "ref": ref38,
                    "alt": alt38,
                    "allele_key": key,
                    "Gene_Symbol": row.get("Gene_Symbol", ""),
                }
            )
            key_rows.append({"assembly": "hg38", "vcf_id": hg38_id, "allele_key": key})

    _write_vcf(hg19_rows, out_dir / "gofcards.hg19.vcf")
    _write_vcf(hg38_rows, out_dir / "gofcards.hg38.vcf")
    write_excel(out_dir / "gofcards_vep_input_key.xlsx", {"vep_input_key": pd.DataFrame(key_rows)})


def _read_vep_tsv(path: str | Path) -> pd.DataFrame:
    in_path = Path(path)
    if not in_path.exists():
        return pd.DataFrame()
    kept: list[str] = []
    with in_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                line = line[1:]
            kept.append(line)
    if not kept:
        return pd.DataFrame()
    return pd.read_csv(StringIO("".join(kept)), sep="\t", dtype=object)


def parse_vep(hg19_vep_tsv: str | Path, hg38_vep_tsv: str | Path, out_xlsx: str | Path) -> None:
    hg19 = _read_vep_tsv(hg19_vep_tsv)
    hg38 = _read_vep_tsv(hg38_vep_tsv)
    if not hg19.empty:
        hg19.insert(0, "assembly", "hg19")
    if not hg38.empty:
        hg38.insert(0, "assembly", "hg38")
    both = pd.concat([hg19, hg38], ignore_index=True) if not (hg19.empty and hg38.empty) else pd.DataFrame()
    write_excel(out_xlsx, {"vep_all": both, "vep_hg19": hg19, "vep_hg38": hg38})

