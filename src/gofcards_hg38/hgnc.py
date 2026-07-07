from __future__ import annotations

import csv
import os
import re
from functools import lru_cache
from pathlib import Path

import pandas as pd


def _clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() == "nan" else text


def _split_symbol_list(value: object) -> list[str]:
    text = _clean(value)
    if not text:
        return []
    return [item.strip() for item in re.split(r"[|,]", text) if item.strip()]


def _candidate_hgnc_tables() -> list[Path]:
    paths: list[Path] = []
    if os.environ.get("HGNC_COMPLETE_SET_TSV"):
        paths.append(Path(os.environ["HGNC_COMPLETE_SET_TSV"]).expanduser())
    if os.environ.get("HGNC_SYMBOL_RESOLVER_CACHE"):
        paths.append(Path(os.environ["HGNC_SYMBOL_RESOLVER_CACHE"]).expanduser() / "hgnc_complete_set.txt")
    paths.append(Path("~/.cache/hgnc-symbol-resolver/hgnc_complete_set.txt").expanduser())
    return paths


class HgncSymbolResolver:
    """Small offline HGNC alias/previous-symbol resolver.

    The standalone hgnc-symbol-resolver package remains the canonical resolver,
    but this workflow only needs deterministic offline symbol normalization.
    Set HGNC_COMPLETE_SET_TSV to the HGNC complete-set TSV to enable it.
    """

    def __init__(self, table: str | Path | None = None) -> None:
        self.table = Path(table).expanduser() if table else self._first_existing_table()
        self._approved: dict[str, str] = {}
        self._aliases: dict[str, set[str]] = {}
        if self.table and self.table.exists():
            self._load(self.table)

    @staticmethod
    def _first_existing_table() -> Path | None:
        for path in _candidate_hgnc_tables():
            if path.exists():
                return path
        return None

    def _load(self, table: Path) -> None:
        with table.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                symbol = _clean(row.get("symbol")).upper()
                if not symbol:
                    continue
                self._approved[symbol] = symbol
                for value in [row.get("alias_symbol"), row.get("prev_symbol")]:
                    for alias in _split_symbol_list(value):
                        self._aliases.setdefault(alias.upper(), set()).add(symbol)

    def resolve(self, value: object) -> str:
        query = _clean(value).upper()
        if not query:
            return ""
        if query in self._approved:
            return self._approved[query]
        matches = self._aliases.get(query, set())
        if len(matches) == 1:
            return next(iter(matches))
        return query


@lru_cache(maxsize=1)
def default_resolver() -> HgncSymbolResolver:
    return HgncSymbolResolver()


def resolve_symbol(value: object) -> str:
    return default_resolver().resolve(value)
