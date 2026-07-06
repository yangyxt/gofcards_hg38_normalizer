# GoFCards HG38 Normalizer

Standalone workflow to rebuild a current GoFCards variant table around stable
variant-level evidence instead of trusting one exported coordinate table.

The workflow:

1. Pulls the GoFCards backend SNV and Indel tables.
2. Downloads and audits the public Excel export.
3. Queries the GoFCards variant summary endpoint for hg38 coordinates.
4. Validates or refills hg38 REF/ALT from a user-supplied hg38 FASTA.
5. Builds VEP input for hg19 and hg38 and parses VEP MANE/HGVS output.
6. Runs TransVar as an optional cross-check, not as the final transcript authority.

## Quick Start

```bash
cd /paedyl01/disk1/yangyxt/gofcards_hg38_normalizer
bin/gofcards_workflow.sh create_env
source ~/.bashrc && mamba activate gofcards_hg38

export HG38_FASTA=/path/to/GRCh38.fa
export HG19_FASTA=/path/to/hg19.fa
export VEP=/path/to/vep
export VEP_CACHE=/path/to/vep/cache

bin/gofcards_workflow.sh run_all
```

Outputs are written to `work/` by default. Set `WORKDIR=/path/to/output` to
use another directory.

## Modular Commands

Each workflow step is callable on its own:

```bash
bin/gofcards_workflow.sh pull_backend_tables
bin/gofcards_workflow.sh download_public_excel
bin/gofcards_workflow.sh join_public_excel
bin/gofcards_workflow.sh augment_hg38
bin/gofcards_workflow.sh validate_hg38_refalt
bin/gofcards_workflow.sh write_vep_inputs
bin/gofcards_workflow.sh run_vep_hg19
bin/gofcards_workflow.sh run_vep_hg38
bin/gofcards_workflow.sh parse_vep_outputs
bin/gofcards_workflow.sh run_transvar_crosscheck
bin/gofcards_workflow.sh build_workbook
```

The Python CLI is installed into the mamba environment by `create_env`:

```bash
source ~/.bashrc && mamba activate gofcards_hg38
gofcards-hg38 --help
```

## Main Artifacts

Default output paths under `work/`:

- `gofcards_backend.xlsx`: backend SNV + Indel table.
- `gofcards_backend.raw.jsonl`: raw backend page responses.
- `gofcards_public.xlsx`: public Excel export.
- `gofcards_backend_vs_public_audit.xlsx`: allele-level backend/public audit.
- `gofcards_augmented_hg38.xlsx`: backend table plus summary endpoint hg38 fields.
- `gofcards_summary_cache.jsonl`: resumable per-allele summary endpoint cache.
- `gofcards_hg38_refalt_checked.xlsx`: hg38 FASTA REF/ALT validation result.
- `vep_inputs/gofcards.hg19.vcf` and `vep_inputs/gofcards.hg38.vcf`: VEP inputs.
- `vep_outputs/gofcards_vep_parsed.xlsx`: parsed VEP MANE/HGVS output.
- `transvar/`: TransVar query files, runner, and optional raw outputs.
- `gofcards_hg38_normalized_workbook.xlsx`: integrated final workbook.

## Smoke Test Status

On this HPC, the backend/public/audit path completed with:

- backend records: 3,160
- public Excel records: 3,161
- backend unique alleles: 2,033
- public unique alleles: 2,034
- matched unique alleles: 1,983
- backend-only unique alleles: 50
- public-only unique alleles: 51

A three-allele summary endpoint test also returned hg38 coordinates
successfully. Full hg38 REF/ALT validation still requires an explicit
`HG38_FASTA` path.

## Design Notes

GoFCards is treated as a curated GoF/DN variant source, but transcript truth is
rebuilt with current VEP/MANE. TransVar is useful for reconciliation and
additional HGVS checks; it is deliberately not used as the final authority.

The backend API currently requires browser-like request headers. The client
sets those headers and retries transient failures. If the API blocks a host, the
public Excel audit can still run, but hg38 augmentation from the summary
endpoint cannot be completed from that host.
