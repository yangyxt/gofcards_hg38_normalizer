# GoFCards HG38 Normalizer

Standalone workflow to rebuild a current GoFCards variant table around stable
variant-level evidence instead of trusting one exported coordinate table.

The ultimate deliverable is `variant_transcript_table` in the final workbook:
one row per GoFCards allele, reference assembly, and VEP transcript annotation,
with clear symbol, transcript ID, HGVSc, HGVSp, hg19 genomic position/ref/alt,
and hg38 genomic position/ref/alt. The companion `preferred_transcript_table`
keeps one preferred transcript row per allele and assembly, ranked by MANE
Select, MANE Plus Clinical, canonical transcript, then any transcript with HGVS.
The workflow also exports `gofcards_priva_exact_gof_hgvsp.tsv.gz`, a compact
runtime cache for PriVA exact variant-level GoF matching.

The workflow:

1. Pulls the GoFCards backend SNV and Indel tables.
2. Downloads and audits the public Excel export.
3. Queries the GoFCards variant summary endpoint for hg38 coordinates.
4. Validates or refills hg38 REF/ALT from a user-supplied hg38 FASTA.
5. VCF-pads deletion/insertion records with blank source REF or ALT using the
   supplied FASTA files, while preserving original GoFCards coordinates.
6. Builds VEP input for hg19 and hg38 and parses VEP MANE/HGVS output.
7. Runs TransVar as an optional cross-check, not as the final transcript authority.

## Quick Start

```bash
cd /paedyl01/disk1/yangyxt/gofcards_hg38_normalizer
bin/gofcards_workflow.sh create_env
source ~/.bashrc && mamba activate gofcards_hg38

export HG38_FASTA=/path/to/GRCh38.fa
export HG19_FASTA=/path/to/hg19.fa
export VEP=/path/to/vep
export VEP_CACHE_HG19=/path/to/GRCh37/vep/cache
export VEP_CACHE_VERSION_HG19=112
# Optional: set these only when a local GRCh38 VEP cache is available.
export VEP_CACHE_HG38=/path/to/GRCh38/vep/cache
export VEP_CACHE_VERSION_HG38=112
# Optional but recommended: enables HGNC alias/previous-symbol normalization
# before GoFCards-to-VEP HGVS concordance ranking.
export HGNC_COMPLETE_SET_TSV=/path/to/hgnc_complete_set.txt

bin/gofcards_workflow.sh run_all
```

Outputs are written to `work/` by default. Set `WORKDIR=/path/to/output` to
use another directory. `run_all` refreshes the integrated workbook and the
compact PriVA TSV cache. If `VEP_CACHE_HG38` is not set, `run_all` skips GRCh38
VEP transcript annotation, while still validating hg38 genomic coordinates and
REF/ALT from the GoFCards summary endpoint plus the supplied hg38 FASTA.

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
bin/gofcards_workflow.sh export_priva_gof_tsv
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
  - `variant_transcript_table`: full transcript-level target table.
  - `preferred_transcript_table`: one preferred transcript row per allele and assembly.
    Representative transcript selection first prioritizes concordance between
    GoFCards RefSeq-style cDNA/protein HGVS and VEP ENST/ENSP HGVS after HGNC
    symbol normalization, then falls back to MANE/canonical status.
- `gofcards_priva_exact_gof_hgvsp.tsv.gz`: compact PriVA cache keyed by
  normalized HGNC symbol plus exact protein change when HGVSp is available, and
  retaining genomic-only rows for future hg19/hg38 position/ref/alt matching.
  This file is for variant-level GoF matching only; it must not be used as
  gene-level GoF evidence.

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
PriVA should consume the compact TSV cache for exact HGVSp matching; a match
means the candidate protein change is represented in GoFCards, not that every
variant in the same gene is GoF.
Rows without HGVSp are retained with `hg19_genomic_key`, `hg19_vcf_key`,
`hg38_genomic_key`, and `hg38_vcf_key` so a later genomic allele matching
function can use them.
Set `HGNC_COMPLETE_SET_TSV` to the official HGNC complete-set TSV when building
the workbook so previous symbols and aliases such as `TMEM173/STING1` and
`PARK2/PRKN` are reconciled before transcript ranking.

The backend API currently requires browser-like request headers. The client
sets those headers and retries transient failures. If the API blocks a host, the
public Excel audit can still run, but hg38 augmentation from the summary
endpoint cannot be completed from that host.
