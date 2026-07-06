from __future__ import annotations

from pathlib import Path
from typing import Mapping

import pandas as pd
from openpyxl.cell.cell import ILLEGAL_CHARACTERS_RE
from openpyxl.styles import Alignment


def _excel_safe_value(value: object) -> object:
    if value is None:
        return ""
    if isinstance(value, (int, float)):
        return value
    return ILLEGAL_CHARACTERS_RE.sub("", str(value))


def ensure_parent(path: str | Path) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    return out


def stringify_for_excel(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in out.columns:
        if out[col].dtype == "object":
            out[col] = out[col].map(_excel_safe_value)
    return out


def write_excel(path: str | Path, sheets: Mapping[str, pd.DataFrame]) -> None:
    out = ensure_parent(path)
    with pd.ExcelWriter(out, engine="openpyxl") as writer:
        for name, df in sheets.items():
            sheet_name = name[:31]
            clean = stringify_for_excel(df)
            clean.to_excel(writer, sheet_name=sheet_name, index=False)
            ws = writer.sheets[sheet_name]
            ws.freeze_panes = "A2"
            ws.auto_filter.ref = ws.dimensions
            for column_cells in ws.columns:
                header = str(column_cells[0].value or "")
                max_len = len(header)
                for cell in column_cells[1:200]:
                    value = "" if cell.value is None else str(cell.value)
                    max_len = max(max_len, min(len(value), 80))
                    cell.alignment = Alignment(wrap_text=True, vertical="top")
                width = min(max(max_len + 2, 10), 48)
                ws.column_dimensions[column_cells[0].column_letter].width = width


def read_excel(path: str | Path, sheet_name: str | int = 0) -> pd.DataFrame:
    return pd.read_excel(path, sheet_name=sheet_name, dtype=object)
