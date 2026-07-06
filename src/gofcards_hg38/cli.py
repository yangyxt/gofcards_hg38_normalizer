from __future__ import annotations

import argparse


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gofcards-hg38",
        description="Normalize GoFCards variant records across backend, public Excel, hg38, VEP, and TransVar.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("pull-backend", help="Pull GoFCards backend SNV and Indel tables.")
    sub.add_parser("download-public-excel", help="Download the public GoFCards Excel export.")
    sub.add_parser("audit-public-excel", help="Audit backend table against public Excel.")
    sub.add_parser("augment-hg38", help="Query GoFCards summary endpoint for hg38 coordinates.")
    sub.add_parser("validate-refalt", help="Validate/refill hg38 REF/ALT from FASTA.")
    sub.add_parser("write-vep-inputs", help="Write hg19 and hg38 VCF inputs for VEP.")
    sub.add_parser("parse-vep", help="Parse VEP tab outputs.")
    sub.add_parser("write-transvar-queries", help="Write TransVar query files and runner.")
    sub.add_parser("build-workbook", help="Build final integrated workbook.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    parser.parse_args(argv)
    parser.error("implementation not installed yet")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())

