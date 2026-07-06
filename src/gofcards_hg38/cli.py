from __future__ import annotations

import argparse

from .audit_public import audit_public_excel
from .augment_hg38 import augment_hg38
from .pull_backend import download_public_excel, pull_backend


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gofcards-hg38",
        description="Normalize GoFCards variant records across backend, public Excel, hg38, VEP, and TransVar.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    pull = sub.add_parser("pull-backend", help="Pull GoFCards backend SNV and Indel tables.")
    pull.add_argument("--out-xlsx", required=True)
    pull.add_argument("--raw-jsonl", required=True)
    pull.add_argument("--page-size", type=int, default=5000)
    pull.set_defaults(func=lambda args: pull_backend(args.out_xlsx, args.raw_jsonl, args.page_size))

    public = sub.add_parser("download-public-excel", help="Download the public GoFCards Excel export.")
    public.add_argument("--url", default=None)
    public.add_argument("--out-xlsx", required=True)
    public.set_defaults(func=lambda args: download_public_excel(args.url, args.out_xlsx))

    audit = sub.add_parser("audit-public-excel", help="Audit backend table against public Excel.")
    audit.add_argument("--backend-xlsx", required=True)
    audit.add_argument("--public-xlsx", required=True)
    audit.add_argument("--out-xlsx", required=True)
    audit.set_defaults(func=lambda args: audit_public_excel(args.backend_xlsx, args.public_xlsx, args.out_xlsx))

    augment = sub.add_parser("augment-hg38", help="Query GoFCards summary endpoint for hg38 coordinates.")
    augment.add_argument("--input-xlsx", required=True)
    augment.add_argument("--out-xlsx", required=True)
    augment.add_argument("--cache-jsonl", required=True)
    augment.add_argument("--sleep-seconds", type=float, default=0.15)
    augment.set_defaults(
        func=lambda args: augment_hg38(args.input_xlsx, args.out_xlsx, args.cache_jsonl, args.sleep_seconds)
    )

    sub.add_parser("validate-refalt", help="Validate/refill hg38 REF/ALT from FASTA.")
    sub.add_parser("write-vep-inputs", help="Write hg19 and hg38 VCF inputs for VEP.")
    sub.add_parser("parse-vep", help="Parse VEP tab outputs.")
    sub.add_parser("write-transvar-queries", help="Write TransVar query files and runner.")
    sub.add_parser("build-workbook", help="Build final integrated workbook.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.error(f"{args.command!r} is not implemented yet")
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
