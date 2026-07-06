from __future__ import annotations

from pathlib import Path

import pandas as pd

from .io_utils import read_excel, write_excel
from .transvar_io import read_transvar_outputs


def _optional_sheet(path: str | Path, sheet_name: str | int = 0) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        return pd.DataFrame()
    try:
        return read_excel(p, sheet_name)
    except Exception:
        return pd.DataFrame()


def build_workbook(
    refalt_xlsx: str | Path,
    audit_xlsx: str | Path,
    vep_xlsx: str | Path,
    transvar_dir: str | Path,
    out_xlsx: str | Path,
) -> None:
    refalt = _optional_sheet(refalt_xlsx, "refalt_checked")
    refalt_summary = _optional_sheet(refalt_xlsx, "summary")
    audit_summary = _optional_sheet(audit_xlsx, "summary")
    audit = _optional_sheet(audit_xlsx, "allele_audit")
    vep_all = _optional_sheet(vep_xlsx, "vep_all")
    sheets: dict[str, pd.DataFrame] = {
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

