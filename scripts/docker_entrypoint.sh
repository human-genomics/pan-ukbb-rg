#!/usr/bin/env bash
set -euo pipefail

cd /app

OUTPUT_DIR="${PAN_UKBB_RG_OUTPUT_DIR:-/app/pipeline-output}"

link_runtime_dir() {
  local name="$1"
  mkdir -p "${OUTPUT_DIR}/${name}"
  if [[ -e "${name}" && ! -L "${name}" ]]; then
    echo "Refusing to replace existing /app/${name}; set PAN_UKBB_RG_OUTPUT_DIR or run from a clean image." >&2
    exit 2
  fi
  ln -sfn "${OUTPUT_DIR}/${name}" "${name}"
}

prepare_runtime_dirs() {
  mkdir -p "${OUTPUT_DIR}"
  link_runtime_dir data
  link_runtime_dir results
  link_runtime_dir logs
}

print_help() {
  cat <<'EOF'
Pan-UKBB EUR LDSC one-vs-all container

Mount a host directory at /app/pipeline-output to persist generated data:

  docker run --rm -v $(pwd)/pan-ukbb-rg-work:/app/pipeline-output \
    ghcr.io/human-genomics/pan-ukbb-rg:latest setup --jobs 8

Commands:
  setup [--jobs N]
      Download manifests/LD scores, prepare all EUR sumstats, and verify LDSC.

  dry-run [run_one_vs_all.py args]
      Resolve a lead phenotype and check whether prepared sumstats are present.
      Example: dry-run --phenocode 20016

  one-vs-all [run_one_vs_all.py args]
      Run one selected GWAS against all EUR GWAS.
      Example: one-vs-all --phenocode 20016 --jobs 8

  make [make args]
      Run a Makefile target directly.

  bash
      Open a shell in /app.

Runtime state is written under /app/pipeline-output/{data,results,logs}.
The image contains software only; the full setup cache is generated at runtime.
EOF
}

make_with_common_args() {
  local target="$1"
  shift
  local make_args=("${target}")
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --jobs)
        if [[ $# -lt 2 ]]; then
          echo "Missing value for --jobs" >&2
          exit 2
        fi
        make_args+=("JOBS=$2")
        shift 2
        ;;
      --jobs=*)
        make_args+=("JOBS=${1#--jobs=}")
        shift
        ;;
      *)
        make_args+=("$1")
        shift
        ;;
    esac
  done
  exec make "${make_args[@]}"
}

run_one_vs_all() {
  local extra_args=()
  if [[ "${1:-}" == "--dry-run" ]]; then
    extra_args+=("--dry-run")
    shift
  fi
  exec python3 scripts/run_one_vs_all.py \
    --manifest data/catalog/eur_gwas_manifest.tsv \
    --sumstats-dir data/sumstats/eur \
    --ldsc-dir "${LDSC_DIR}" \
    --ldsc-python "${LDSC_PYTHON}" \
    --ld-prefix data/ld/UKBB.EUR \
    "${extra_args[@]}" \
    "$@"
}

prepare_runtime_dirs

case "${1:-help}" in
  help|--help|-h)
    print_help
    ;;
  setup|prepare-all-sumstats|prepare-sumstats|benchmark90|validate-catalog|catalog|prepare-ldscores)
    cmd="$1"
    shift
    make_with_common_args "${cmd}" "$@"
    ;;
  dry-run|one-vs-all-dry-run)
    shift
    run_one_vs_all --dry-run "$@"
    ;;
  one-vs-all)
    shift
    run_one_vs_all "$@"
    ;;
  make)
    shift
    exec make "$@"
    ;;
  bash|sh|python3|python|make)
    exec "$@"
    ;;
  *)
    exec "$@"
    ;;
esac
