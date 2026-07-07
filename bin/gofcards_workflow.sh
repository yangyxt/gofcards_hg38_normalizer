#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORKDIR="${WORKDIR:-${REPO_ROOT}/work}"
CONDA_ENV="${CONDA_ENV:-gofcards_hg38}"

PUBLIC_EXCEL_URL="${PUBLIC_EXCEL_URL:-https://download.genemed.tech/upload/GainFunCards/gofcards_data_download.xlsx}"
BACKEND_XLSX="${BACKEND_XLSX:-${WORKDIR}/gofcards_backend.xlsx}"
BACKEND_JSONL="${BACKEND_JSONL:-${WORKDIR}/gofcards_backend.raw.jsonl}"
PUBLIC_XLSX="${PUBLIC_XLSX:-${WORKDIR}/gofcards_public.xlsx}"
AUDIT_XLSX="${AUDIT_XLSX:-${WORKDIR}/gofcards_backend_vs_public_audit.xlsx}"
AUGMENTED_XLSX="${AUGMENTED_XLSX:-${WORKDIR}/gofcards_augmented_hg38.xlsx}"
SUMMARY_CACHE="${SUMMARY_CACHE:-${WORKDIR}/gofcards_summary_cache.jsonl}"
REFALT_XLSX="${REFALT_XLSX:-${WORKDIR}/gofcards_hg38_refalt_checked.xlsx}"
VEP_INPUT_DIR="${VEP_INPUT_DIR:-${WORKDIR}/vep_inputs}"
VEP_OUTPUT_DIR="${VEP_OUTPUT_DIR:-${WORKDIR}/vep_outputs}"
TRANSVAR_DIR="${TRANSVAR_DIR:-${WORKDIR}/transvar}"
FINAL_XLSX="${FINAL_XLSX:-${WORKDIR}/gofcards_hg38_normalized_workbook.xlsx}"
PRIVA_GOF_TSV="${PRIVA_GOF_TSV:-${WORKDIR}/gofcards_priva_exact_gof_hgvsp.tsv.gz}"

mkdir -p "${WORKDIR}"

usage() {
  cat <<'USAGE'
Usage: bin/gofcards_workflow.sh <function>

Functions:
  create_env              Create/update the mamba environment.
  pull_backend_tables     Pull GoFCards backend SNV and Indel records.
  download_public_excel   Download the public GoFCards Excel export.
  join_public_excel       Compare backend records with public Excel.
  augment_hg38            Query summary endpoint for hg38 coordinates.
  validate_hg38_refalt    Validate/refill hg38 REF/ALT from HG38_FASTA.
  write_vep_inputs        Write hg19 and hg38 VCFs for VEP.
  run_vep_hg19            Run VEP on the hg19 input VCF.
  run_vep_hg38            Run VEP on the hg38 input VCF.
  parse_vep_outputs       Parse VEP tab outputs into one workbook.
  run_transvar_crosscheck Generate and optionally run TransVar checks.
  build_workbook          Merge core, VEP, and TransVar outputs.
  export_priva_gof_tsv    Export compact TSV for PriVA exact GoF variant matching.
  run_all                 Run all non-optional steps, plus VEP/TransVar when configured.

Required env vars for selected steps:
  HG38_FASTA              Required by validate_hg38_refalt and hg38 VEP.
  HG19_FASTA              Required by hg19 VEP; used by validate_hg38_refalt to VCF-pad hg19 indels when set.
  VEP                     Optional; defaults to "vep" if available.
  VEP_CACHE               Optional shared VEP cache directory.
  VEP_CACHE_HG19          Optional GRCh37/hg19 VEP cache directory.
  VEP_CACHE_HG38          Optional GRCh38 VEP cache directory. If unset, run_all skips hg38 VEP.
  HGNC_COMPLETE_SET_TSV   Optional HGNC complete-set TSV for symbol alias normalization.
  PRIVA_GOF_TSV           Optional compact PriVA exact GoF TSV output path.
USAGE
}

create_env() {
  source ~/.bashrc
  if mamba env list | awk '{print $1}' | grep -qx "${CONDA_ENV}"; then
    mamba env update -y -n "${CONDA_ENV}" -f "${REPO_ROOT}/environment.yml" --prune
  else
    mamba env create -y -n "${CONDA_ENV}" -f "${REPO_ROOT}/environment.yml"
  fi
  mamba run -n "${CONDA_ENV}" python -m pip install -e "${REPO_ROOT}"
}

run_py() {
  PYTHONNOUSERSITE=1 PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}" python -m gofcards_hg38.cli "$@"
}

pull_backend_tables() {
  mkdir -p "${WORKDIR}"
  run_py pull-backend --out-xlsx "${BACKEND_XLSX}" --raw-jsonl "${BACKEND_JSONL}"
}

download_public_excel() {
  mkdir -p "$(dirname "${PUBLIC_XLSX}")"
  run_py download-public-excel --url "${PUBLIC_EXCEL_URL}" --out-xlsx "${PUBLIC_XLSX}"
}

join_public_excel() {
  run_py audit-public-excel \
    --backend-xlsx "${BACKEND_XLSX}" \
    --public-xlsx "${PUBLIC_XLSX}" \
    --out-xlsx "${AUDIT_XLSX}"
}

augment_hg38() {
  run_py augment-hg38 \
    --input-xlsx "${BACKEND_XLSX}" \
    --out-xlsx "${AUGMENTED_XLSX}" \
    --cache-jsonl "${SUMMARY_CACHE}" \
    --workers "${SUMMARY_WORKERS:-1}"
}

validate_hg38_refalt() {
  : "${HG38_FASTA:?Set HG38_FASTA before validate_hg38_refalt}"
  local args=(
    --input-xlsx "${AUGMENTED_XLSX}"
    --hg38-fasta "${HG38_FASTA}"
    --out-xlsx "${REFALT_XLSX}"
  )
  if [[ -n "${HG19_FASTA:-}" ]]; then
    args+=(--hg19-fasta "${HG19_FASTA}")
  fi
  run_py validate-refalt "${args[@]}"
}

write_vep_inputs() {
  mkdir -p "${VEP_INPUT_DIR}"
  run_py write-vep-inputs \
    --input-xlsx "${REFALT_XLSX}" \
    --out-dir "${VEP_INPUT_DIR}"
}

vep_common_args() {
  local assembly="$1"
  local input_vcf="$2"
  local output_tsv="$3"
  local fasta="$4"
  local cache_dir="${5:-${VEP_CACHE:-}}"
  local cache_version="${6:-${VEP_CACHE_VERSION:-}}"
  local cache_merged="${7:-${VEP_CACHE_MERGED:-}}"
  local vep_bin="${VEP:-vep}"
  local cache_args=()
  if [[ -n "${cache_dir}" ]]; then
    cache_args=(--cache --offline --dir_cache "${cache_dir}")
    if [[ -n "${cache_version}" ]]; then
      cache_args+=(--cache_version "${cache_version}")
    fi
    if [[ "${cache_merged}" == "1" || "${cache_merged}" == "true" || "${cache_merged}" == "TRUE" || "${cache_merged}" == "yes" || "${cache_merged}" == "YES" ]]; then
      cache_args+=(--merged)
    fi
  fi
  "${vep_bin}" \
    --input_file "${input_vcf}" \
    --output_file "${output_tsv}" \
    --format vcf \
    --tab \
    --force_overwrite \
    --assembly "${assembly}" \
    --fasta "${fasta}" \
    --symbol --hgvs --mane --canonical --transcript_version \
    --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,HGVSc,HGVSp,MANE_SELECT,MANE_PLUS_CLINICAL,CANONICAL,SYMBOL,Existing_variation" \
    "${cache_args[@]}"
}

run_vep_hg19() {
  : "${HG19_FASTA:?Set HG19_FASTA before run_vep_hg19}"
  mkdir -p "${VEP_OUTPUT_DIR}"
  vep_common_args GRCh37 "${VEP_INPUT_DIR}/gofcards.hg19.vcf" "${VEP_OUTPUT_DIR}/gofcards.hg19.vep.tsv" "${HG19_FASTA}" "${VEP_CACHE_HG19:-${VEP_CACHE:-}}" "${VEP_CACHE_VERSION_HG19:-${VEP_CACHE_VERSION:-}}" "${VEP_CACHE_MERGED_HG19:-${VEP_CACHE_MERGED:-}}"
}

run_vep_hg38() {
  : "${HG38_FASTA:?Set HG38_FASTA before run_vep_hg38}"
  mkdir -p "${VEP_OUTPUT_DIR}"
  vep_common_args GRCh38 "${VEP_INPUT_DIR}/gofcards.hg38.vcf" "${VEP_OUTPUT_DIR}/gofcards.hg38.vep.tsv" "${HG38_FASTA}" "${VEP_CACHE_HG38:-${VEP_CACHE:-}}" "${VEP_CACHE_VERSION_HG38:-${VEP_CACHE_VERSION:-}}" "${VEP_CACHE_MERGED_HG38:-${VEP_CACHE_MERGED:-}}"
}

parse_vep_outputs() {
  run_py parse-vep \
    --hg19-vep-tsv "${VEP_OUTPUT_DIR}/gofcards.hg19.vep.tsv" \
    --hg38-vep-tsv "${VEP_OUTPUT_DIR}/gofcards.hg38.vep.tsv" \
    --vep-input-key-xlsx "${VEP_INPUT_DIR}/gofcards_vep_input_key.xlsx" \
    --out-xlsx "${VEP_OUTPUT_DIR}/gofcards_vep_parsed.xlsx"
}

run_transvar_crosscheck() {
  mkdir -p "${TRANSVAR_DIR}"
  run_py write-transvar-queries \
    --input-xlsx "${REFALT_XLSX}" \
    --out-dir "${TRANSVAR_DIR}"
  if command -v transvar >/dev/null 2>&1; then
    bash "${TRANSVAR_DIR}/run_transvar.sh" || echo "TransVar cross-check failed; continuing because it is not the final transcript authority." >&2
  else
    echo "transvar is not on PATH; generated ${TRANSVAR_DIR}/run_transvar.sh only." >&2
  fi
}

build_workbook() {
  run_py build-workbook \
    --refalt-xlsx "${REFALT_XLSX}" \
    --audit-xlsx "${AUDIT_XLSX}" \
    --vep-xlsx "${VEP_OUTPUT_DIR}/gofcards_vep_parsed.xlsx" \
    --transvar-dir "${TRANSVAR_DIR}" \
    --out-xlsx "${FINAL_XLSX}"
}

export_priva_gof_tsv() {
  run_py export-priva-gof-tsv \
    --workbook-xlsx "${FINAL_XLSX}" \
    --out-tsv "${PRIVA_GOF_TSV}"
}

run_all() {
  pull_backend_tables
  download_public_excel
  join_public_excel
  augment_hg38
  validate_hg38_refalt
  write_vep_inputs
  local ran_any_vep=0
  if command -v "${VEP:-vep}" >/dev/null 2>&1 && [[ -n "${HG19_FASTA:-}" ]]; then
    run_vep_hg19
    ran_any_vep=1
  else
    echo "Skipping hg19 VEP: set VEP/HG19_FASTA and ensure VEP is on PATH." >&2
  fi
  if command -v "${VEP:-vep}" >/dev/null 2>&1 && [[ -n "${HG38_FASTA:-}" && -n "${VEP_CACHE_HG38:-${VEP_CACHE:-}}" ]]; then
    run_vep_hg38
    ran_any_vep=1
  else
    echo "Skipping hg38 VEP: set VEP/HG38_FASTA and a local VEP_CACHE_HG38. hg38 genomic REF/ALT validation still runs." >&2
  fi
  if [[ "${ran_any_vep}" -eq 1 ]]; then
    parse_vep_outputs
  else
    echo "Skipping VEP parsing: no VEP run completed." >&2
  fi
  run_transvar_crosscheck
  build_workbook
  export_priva_gof_tsv
}

main() {
  local cmd="${1:-}"
  if [[ -z "${cmd}" || "${cmd}" == "-h" || "${cmd}" == "--help" ]]; then
    usage
    exit 0
  fi
  shift || true
  case "${cmd}" in
    create_env|pull_backend_tables|download_public_excel|join_public_excel|augment_hg38|validate_hg38_refalt|write_vep_inputs|run_vep_hg19|run_vep_hg38|parse_vep_outputs|run_transvar_crosscheck|build_workbook|export_priva_gof_tsv|run_all)
      "${cmd}" "$@"
      ;;
    *)
      usage >&2
      exit 2
      ;;
  esac
}

main "$@"
