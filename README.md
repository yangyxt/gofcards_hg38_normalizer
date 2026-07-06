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

## Design Notes

GoFCards is treated as a curated GoF/DN variant source, but transcript truth is
rebuilt with current VEP/MANE. TransVar is useful for reconciliation and
additional HGVS checks; it is deliberately not used as the final authority.

The backend API currently requires browser-like request headers. The client
sets those headers and retries transient failures. If the API blocks a host, the
public Excel audit can still run, but hg38 augmentation from the summary
endpoint cannot be completed from that host.

