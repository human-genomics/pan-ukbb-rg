SHELL := /bin/bash

PYTHON ?= python3
LDSC_ENV_PREFIX ?= .envs/ldsc-neale
LDSC_PYTHON ?= $(LDSC_ENV_PREFIX)/bin/python
ENV_MANAGER ?= micromamba
JOBS ?= 16
BENCH_N ?= 90
RAYON_THREADS ?= 50
MAX_PARALLEL_SHARDS ?= 1
TRAIT_BLOCK_SIZE ?= 256
EXTERNAL_RG_TRAIT_BLOCK_SIZE ?= 64
CARGO_BUILD_JOBS ?= 16
VALIDATION_PAIRS ?= 100
VALIDATION_SEED ?= 20260527
VALIDATION_JOBS ?= 8
VALIDATION_RAYON_THREADS ?= 8
EXTERNAL_GWAS_TARGETS ?= config/external_gwas_targets.tsv
EXTERNAL_GWAS_DIR ?= data/external_gwas
EXTERNAL_GWAS_SOURCES ?= gwas_catalog,zenodo_indep107,figshare_file,direct_url,nextcloud_file
EXTERNAL_GWAS_INCLUDE ?=
EXTERNAL_GWAS_LIMIT ?=
EXTERNAL_GWAS_MAX_ROWS ?=
EXTERNAL_GWAS_ZENODO_ARCHIVE ?=
EXTERNAL_GWAS_RSID_REFERENCE ?=
EXTERNAL_GWAS_URL_OVERRIDES ?= config/external_gwas_restricted_urls.tsv
EXTERNAL_GWAS_STRICT ?=
PAN_GWAS_MANIFEST ?= data/catalog/eur_gwas_manifest.tsv
PAN_SUMSTATS_DIR ?= data/sumstats/eur
EXTERNAL_RG_MANIFEST ?= $(EXTERNAL_GWAS_DIR)/catalog/external_gwas_manifest.tsv
EXTERNAL_RG_DIR ?= results/external_rg
EXTERNAL_RG_SUMSTATS_DIR ?= data/external_rg/sumstats
EXTERNAL_RG_INCREMENTAL_DIR ?= results/external_rg_neale_incremental
EXTERNAL_RG_INCREMENTAL_SUMSTATS_DIR ?= data/external_rg_neale_incremental/sumstats
EXTERNAL_RG_INCREMENTAL_TRAIT_PREFIX ?= NealeLab.
EXTERNAL_RG_INCREMENTAL_TRAIT_IDS ?=
EXTERNAL_RG_INCREMENTAL_RAYON_THREADS ?= 8
EXTERNAL_RG_INCREMENTAL_MAX_PARALLEL_SHARDS ?= 1

PAN_BASE := https://pan-ukb-us-east-1.s3.amazonaws.com
MANIFEST := data/manifests/phenotype_manifest.tsv.bgz
H2_MANIFEST := data/manifests/h2_manifest.tsv.bgz
BENCH_PHENOS := data/benchmark90/phenotypes.tsv
CATALOG := data/catalog/eur_gwas_manifest.tsv
SUMSTATS_DIR := data/sumstats/eur
LD_PREFIX := data/ld/UKBB.EUR
LDSC_DIR ?= external/ldsc-neale
LDSC_RS_DIR ?= external/ldsc-rs-rg-batch
LDSC_RS_TARGET ?= external/ldsc-rs-rg-batch-target
LDSC_RS_BIN ?= $(LDSC_RS_TARGET)/release/ldsc
ALL_RG_DIR ?= results/all_rg
ALL_RG_VALIDATION_DIR ?= results/all_rg_validation

.PHONY: all init setup fetch-manifests catalog validate-catalog select-benchmark prepare-ldscores setup-ldsc setup-ldsc-env prepare-sumstats prepare-all-sumstats external-gwas-manifest external-gwas external-gwas-smoke external-rg-prepare external-rg external-rg-dry-run external-rg-progress external-rg-collect external-rg-incremental-prepare external-rg-incremental external-rg-incremental-dry-run external-rg-incremental-progress external-rg-incremental-collect one-vs-all one-vs-all-dry-run setup-ldsc-rs-rg-batch all-rg-prepare all-rg all-rg-dry-run all-rg-progress all-rg-collect all-rg-validation run-benchmark summarize benchmark90 hardware clean-small

all: setup

init:
	@git init

fetch-manifests:
	@mkdir -p data/manifests
	curl -fL --retry 5 --retry-delay 5 -o $(MANIFEST) $(PAN_BASE)/sumstats_release/phenotype_manifest.tsv.bgz
	curl -fL --retry 5 --retry-delay 5 -o $(H2_MANIFEST) $(PAN_BASE)/sumstats_release/h2_manifest.tsv.bgz

catalog: fetch-manifests
	@mkdir -p data/catalog
	$(PYTHON) scripts/build_eur_gwas_catalog.py \
		--phenotype-manifest $(MANIFEST) \
		--out $(CATALOG)

validate-catalog: catalog
	$(PYTHON) scripts/validate_eur_gwas_catalog.py \
		--phenotype-manifest $(MANIFEST) \
		--catalog $(CATALOG)

select-benchmark: fetch-manifests
	@mkdir -p data/benchmark90
	$(PYTHON) scripts/select_benchmark_phenotypes.py \
		--manifest $(MANIFEST) \
		--config config/benchmark90.yaml \
		--out $(BENCH_PHENOS)

prepare-ldscores:
	@mkdir -p data/ld
	curl -fL --retry 5 --retry-delay 5 -o $(LD_PREFIX).l2.ldscore.gz $(PAN_BASE)/ld_release/UKBB.EUR.l2.ldscore.gz
	curl -fL --retry 5 --retry-delay 5 -o $(LD_PREFIX).l2.M $(PAN_BASE)/ld_release/UKBB.EUR.l2.M
	curl -fL --retry 5 --retry-delay 5 -o $(LD_PREFIX).l2.M_5_50 $(PAN_BASE)/ld_release/UKBB.EUR.l2.M_5_50
	gzip -dc $(LD_PREFIX).l2.ldscore.gz | awk 'NR>1 {print $$2}' > data/ld/UKBB.EUR.snps

setup-ldsc:
	bash scripts/setup_neale_ldsc.sh

setup-ldsc-env:
	@if [ -x "$(LDSC_ENV_PREFIX)/bin/python" ] && bash scripts/check_ldsc_env.sh "$(LDSC_ENV_PREFIX)/bin/python" >/dev/null; then \
		echo "LDSC environment already exists at $(LDSC_ENV_PREFIX)"; \
	else \
		if [ -e "$(LDSC_ENV_PREFIX)" ]; then \
			echo "Removing incomplete or invalid LDSC environment at $(LDSC_ENV_PREFIX)"; \
			rm -rf "$(LDSC_ENV_PREFIX)"; \
		fi; \
		if [ "$(ENV_MANAGER)" = "conda" ] || [ "$(ENV_MANAGER)" = "mamba" ]; then \
			$(ENV_MANAGER) env create -y -p $(LDSC_ENV_PREFIX) -f envs/ldsc-neale.yml; \
		else \
			$(ENV_MANAGER) create -y -p $(LDSC_ENV_PREFIX) -f envs/ldsc-neale.yml; \
		fi; \
		bash scripts/check_ldsc_env.sh "$(LDSC_ENV_PREFIX)/bin/python"; \
	fi

setup: setup-ldsc-env validate-catalog prepare-ldscores setup-ldsc prepare-all-sumstats

prepare-sumstats: $(BENCH_PHENOS) prepare-ldscores
	mkdir -p data/benchmark90/sumstats results/benchmark90/prepare_stats logs/prepare_sumstats
	$(PYTHON) scripts/prepare_sumstats_batch.py \
		--phenotypes $(BENCH_PHENOS) \
		--ld-snps data/ld/UKBB.EUR.snps \
		--out-dir data/benchmark90/sumstats \
		--stats-dir results/benchmark90/prepare_stats \
		--jobs $(JOBS)

prepare-all-sumstats: catalog prepare-ldscores
	mkdir -p $(SUMSTATS_DIR) results/prepare_all/prepare_stats logs/prepare_all_sumstats
	$(PYTHON) scripts/prepare_sumstats_batch.py \
		--phenotypes $(CATALOG) \
		--ld-snps data/ld/UKBB.EUR.snps \
		--out-dir $(SUMSTATS_DIR) \
		--stats-dir results/prepare_all/prepare_stats \
		--log-dir logs/prepare_all_sumstats \
		--jobs $(JOBS)

external-gwas-manifest:
	mkdir -p $(EXTERNAL_GWAS_DIR)
	$(PYTHON) scripts/prepare_external_gwas.py \
		--targets $(EXTERNAL_GWAS_TARGETS) \
		--ld-snps data/ld/UKBB.EUR.snps \
		--out-dir $(EXTERNAL_GWAS_DIR) \
		--sources "$(EXTERNAL_GWAS_SOURCES)" \
		--url-overrides "$(EXTERNAL_GWAS_URL_OVERRIDES)" \
		$(if $(EXTERNAL_GWAS_INCLUDE),--include "$(EXTERNAL_GWAS_INCLUDE)",) \
		$(if $(EXTERNAL_GWAS_LIMIT),--limit $(EXTERNAL_GWAS_LIMIT),) \
		--manifest-only

external-gwas:
	mkdir -p $(EXTERNAL_GWAS_DIR)
	$(PYTHON) scripts/prepare_external_gwas.py \
		--targets $(EXTERNAL_GWAS_TARGETS) \
		--ld-snps data/ld/UKBB.EUR.snps \
		--out-dir $(EXTERNAL_GWAS_DIR) \
		--sources "$(EXTERNAL_GWAS_SOURCES)" \
		--url-overrides "$(EXTERNAL_GWAS_URL_OVERRIDES)" \
		$(if $(EXTERNAL_GWAS_INCLUDE),--include "$(EXTERNAL_GWAS_INCLUDE)",) \
		$(if $(EXTERNAL_GWAS_LIMIT),--limit $(EXTERNAL_GWAS_LIMIT),) \
		$(if $(EXTERNAL_GWAS_MAX_ROWS),--max-rows $(EXTERNAL_GWAS_MAX_ROWS),) \
		$(if $(EXTERNAL_GWAS_ZENODO_ARCHIVE),--zenodo-archive "$(EXTERNAL_GWAS_ZENODO_ARCHIVE)",) \
		$(if $(EXTERNAL_GWAS_RSID_REFERENCE),--rsid-reference "$(EXTERNAL_GWAS_RSID_REFERENCE)",) \
		$(if $(EXTERNAL_GWAS_STRICT),--strict,)

external-gwas-smoke:
	mkdir -p tmp/external_gwas_smoke
	$(PYTHON) scripts/prepare_external_gwas.py \
		--targets $(EXTERNAL_GWAS_TARGETS) \
		--ld-snps data/ld/UKBB.EUR.snps \
		--out-dir tmp/external_gwas_smoke \
		--sources gwas_catalog \
		--include GCST90011874 \
		--max-rows 200000 \
		--strict

external-rg-prepare:
	$(PYTHON) scripts/run_external_rg_hybrid.py \
		--pan-manifest $(PAN_GWAS_MANIFEST) \
		--pan-sumstats-dir $(PAN_SUMSTATS_DIR) \
		--external-manifest $(EXTERNAL_RG_MANIFEST) \
		--combined-sumstats-dir $(EXTERNAL_RG_SUMSTATS_DIR) \
		--ld-prefix $(LD_PREFIX) \
		--ldsc-bin $(LDSC_RS_BIN) \
		--out-dir $(EXTERNAL_RG_DIR) \
		--trait-block-size $(EXTERNAL_RG_TRAIT_BLOCK_SIZE) \
		--rayon-threads $(RAYON_THREADS) \
		--max-parallel-shards $(MAX_PARALLEL_SHARDS) \
		--force-shards \
		--force-symlinks \
		--prepare-only

external-rg:
	$(PYTHON) scripts/run_external_rg_hybrid.py \
		--pan-manifest $(PAN_GWAS_MANIFEST) \
		--pan-sumstats-dir $(PAN_SUMSTATS_DIR) \
		--external-manifest $(EXTERNAL_RG_MANIFEST) \
		--combined-sumstats-dir $(EXTERNAL_RG_SUMSTATS_DIR) \
		--ld-prefix $(LD_PREFIX) \
		--ldsc-bin $(LDSC_RS_BIN) \
		--out-dir $(EXTERNAL_RG_DIR) \
		--trait-block-size $(EXTERNAL_RG_TRAIT_BLOCK_SIZE) \
		--rayon-threads $(RAYON_THREADS) \
		--max-parallel-shards $(MAX_PARALLEL_SHARDS) \
		--force-shards \
		--force-symlinks

external-rg-dry-run:
	$(PYTHON) scripts/run_external_rg_hybrid.py \
		--pan-manifest $(PAN_GWAS_MANIFEST) \
		--pan-sumstats-dir $(PAN_SUMSTATS_DIR) \
		--external-manifest $(EXTERNAL_RG_MANIFEST) \
		--combined-sumstats-dir $(EXTERNAL_RG_SUMSTATS_DIR) \
		--ld-prefix $(LD_PREFIX) \
		--ldsc-bin $(LDSC_RS_BIN) \
		--out-dir $(EXTERNAL_RG_DIR) \
		--trait-block-size $(EXTERNAL_RG_TRAIT_BLOCK_SIZE) \
		--rayon-threads $(RAYON_THREADS) \
		--max-parallel-shards $(MAX_PARALLEL_SHARDS) \
		--force-shards \
		--force-symlinks \
		--dry-run

external-rg-progress:
	$(PYTHON) scripts/run_external_rg_hybrid.py \
		--combined-sumstats-dir $(EXTERNAL_RG_SUMSTATS_DIR) \
		--ld-prefix $(LD_PREFIX) \
		--ldsc-bin $(LDSC_RS_BIN) \
		--out-dir $(EXTERNAL_RG_DIR) \
		--rayon-threads $(RAYON_THREADS) \
		--max-parallel-shards $(MAX_PARALLEL_SHARDS) \
		--progress

external-rg-collect:
	$(PYTHON) scripts/run_external_rg_hybrid.py \
		--combined-sumstats-dir $(EXTERNAL_RG_SUMSTATS_DIR) \
		--ld-prefix $(LD_PREFIX) \
		--ldsc-bin $(LDSC_RS_BIN) \
		--out-dir $(EXTERNAL_RG_DIR) \
		--rayon-threads $(RAYON_THREADS) \
		--max-parallel-shards $(MAX_PARALLEL_SHARDS) \
		--collect

external-rg-incremental-prepare:
	$(PYTHON) scripts/run_external_rg_hybrid.py \
		--pan-manifest $(PAN_GWAS_MANIFEST) \
		--pan-sumstats-dir $(PAN_SUMSTATS_DIR) \
		--external-manifest $(EXTERNAL_RG_MANIFEST) \
		--combined-sumstats-dir $(EXTERNAL_RG_INCREMENTAL_SUMSTATS_DIR) \
		--ld-prefix $(LD_PREFIX) \
		--ldsc-bin $(LDSC_RS_BIN) \
		--out-dir $(EXTERNAL_RG_INCREMENTAL_DIR) \
		--trait-block-size $(EXTERNAL_RG_TRAIT_BLOCK_SIZE) \
		--rayon-threads $(EXTERNAL_RG_INCREMENTAL_RAYON_THREADS) \
		--max-parallel-shards $(EXTERNAL_RG_INCREMENTAL_MAX_PARALLEL_SHARDS) \
		$(if $(EXTERNAL_RG_INCREMENTAL_TRAIT_PREFIX),--pair-include-trait-prefix "$(EXTERNAL_RG_INCREMENTAL_TRAIT_PREFIX)",) \
		$(if $(EXTERNAL_RG_INCREMENTAL_TRAIT_IDS),--pair-include-trait-id "$(EXTERNAL_RG_INCREMENTAL_TRAIT_IDS)",) \
		--force-shards \
		--force-symlinks \
		--prepare-only

external-rg-incremental:
	$(PYTHON) scripts/run_external_rg_hybrid.py \
		--pan-manifest $(PAN_GWAS_MANIFEST) \
		--pan-sumstats-dir $(PAN_SUMSTATS_DIR) \
		--external-manifest $(EXTERNAL_RG_MANIFEST) \
		--combined-sumstats-dir $(EXTERNAL_RG_INCREMENTAL_SUMSTATS_DIR) \
		--ld-prefix $(LD_PREFIX) \
		--ldsc-bin $(LDSC_RS_BIN) \
		--out-dir $(EXTERNAL_RG_INCREMENTAL_DIR) \
		--trait-block-size $(EXTERNAL_RG_TRAIT_BLOCK_SIZE) \
		--rayon-threads $(EXTERNAL_RG_INCREMENTAL_RAYON_THREADS) \
		--max-parallel-shards $(EXTERNAL_RG_INCREMENTAL_MAX_PARALLEL_SHARDS) \
		$(if $(EXTERNAL_RG_INCREMENTAL_TRAIT_PREFIX),--pair-include-trait-prefix "$(EXTERNAL_RG_INCREMENTAL_TRAIT_PREFIX)",) \
		$(if $(EXTERNAL_RG_INCREMENTAL_TRAIT_IDS),--pair-include-trait-id "$(EXTERNAL_RG_INCREMENTAL_TRAIT_IDS)",) \
		--force-shards \
		--force-symlinks

external-rg-incremental-dry-run:
	$(PYTHON) scripts/run_external_rg_hybrid.py \
		--pan-manifest $(PAN_GWAS_MANIFEST) \
		--pan-sumstats-dir $(PAN_SUMSTATS_DIR) \
		--external-manifest $(EXTERNAL_RG_MANIFEST) \
		--combined-sumstats-dir $(EXTERNAL_RG_INCREMENTAL_SUMSTATS_DIR) \
		--ld-prefix $(LD_PREFIX) \
		--ldsc-bin $(LDSC_RS_BIN) \
		--out-dir $(EXTERNAL_RG_INCREMENTAL_DIR) \
		--trait-block-size $(EXTERNAL_RG_TRAIT_BLOCK_SIZE) \
		--rayon-threads $(EXTERNAL_RG_INCREMENTAL_RAYON_THREADS) \
		--max-parallel-shards $(EXTERNAL_RG_INCREMENTAL_MAX_PARALLEL_SHARDS) \
		$(if $(EXTERNAL_RG_INCREMENTAL_TRAIT_PREFIX),--pair-include-trait-prefix "$(EXTERNAL_RG_INCREMENTAL_TRAIT_PREFIX)",) \
		$(if $(EXTERNAL_RG_INCREMENTAL_TRAIT_IDS),--pair-include-trait-id "$(EXTERNAL_RG_INCREMENTAL_TRAIT_IDS)",) \
		--force-shards \
		--force-symlinks \
		--dry-run

external-rg-incremental-progress:
	$(PYTHON) scripts/run_external_rg_hybrid.py \
		--combined-sumstats-dir $(EXTERNAL_RG_INCREMENTAL_SUMSTATS_DIR) \
		--ld-prefix $(LD_PREFIX) \
		--ldsc-bin $(LDSC_RS_BIN) \
		--out-dir $(EXTERNAL_RG_INCREMENTAL_DIR) \
		--rayon-threads $(EXTERNAL_RG_INCREMENTAL_RAYON_THREADS) \
		--max-parallel-shards $(EXTERNAL_RG_INCREMENTAL_MAX_PARALLEL_SHARDS) \
		--progress

external-rg-incremental-collect:
	$(PYTHON) scripts/run_external_rg_hybrid.py \
		--combined-sumstats-dir $(EXTERNAL_RG_INCREMENTAL_SUMSTATS_DIR) \
		--ld-prefix $(LD_PREFIX) \
		--ldsc-bin $(LDSC_RS_BIN) \
		--out-dir $(EXTERNAL_RG_INCREMENTAL_DIR) \
		--rayon-threads $(EXTERNAL_RG_INCREMENTAL_RAYON_THREADS) \
		--max-parallel-shards $(EXTERNAL_RG_INCREMENTAL_MAX_PARALLEL_SHARDS) \
		--collect

one-vs-all: setup-ldsc setup-ldsc-env catalog prepare-ldscores
	@if [ -z "$(PHENOCODE)$(PHENOTYPE_ID)$(QUERY)" ]; then \
		echo "Set one selector, e.g. make one-vs-all PHENOCODE=20016 JOBS=16"; \
		exit 2; \
	fi
	$(PYTHON) scripts/run_one_vs_all.py \
		--manifest $(CATALOG) \
		--sumstats-dir $(SUMSTATS_DIR) \
		--ldsc-dir $(LDSC_DIR) \
		--ldsc-python "$(LDSC_PYTHON)" \
		--ld-prefix $(LD_PREFIX) \
		--jobs $(JOBS) \
		$(if $(PHENOCODE),--phenocode $(PHENOCODE),) \
		$(if $(PHENOTYPE_ID),--phenotype-id $(PHENOTYPE_ID),) \
		$(if $(QUERY),--query "$(QUERY)",)

one-vs-all-dry-run: setup-ldsc setup-ldsc-env catalog prepare-ldscores
	@if [ -z "$(PHENOCODE)$(PHENOTYPE_ID)$(QUERY)" ]; then \
		echo "Set one selector, e.g. make one-vs-all-dry-run PHENOCODE=20016"; \
		exit 2; \
	fi
	$(PYTHON) scripts/run_one_vs_all.py \
		--manifest $(CATALOG) \
		--sumstats-dir $(SUMSTATS_DIR) \
		--ldsc-dir $(LDSC_DIR) \
		--ldsc-python "$(LDSC_PYTHON)" \
		--ld-prefix $(LD_PREFIX) \
		--jobs $(JOBS) \
		--dry-run \
		$(if $(PHENOCODE),--phenocode $(PHENOCODE),) \
		$(if $(PHENOTYPE_ID),--phenotype-id $(PHENOTYPE_ID),) \
		$(if $(QUERY),--query "$(QUERY)",)

setup-ldsc-rs-rg-batch:
	LDSC_RS_DIR="$(LDSC_RS_DIR)" \
	LDSC_RS_TARGET="$(LDSC_RS_TARGET)" \
	CARGO_BUILD_JOBS="$(CARGO_BUILD_JOBS)" \
		bash scripts/setup_ldsc_rs_rg_batch.sh

all-rg-prepare:
	$(PYTHON) scripts/run_all_rg_hybrid.py \
		--manifest $(CATALOG) \
		--sumstats-dir $(SUMSTATS_DIR) \
		--ld-prefix $(LD_PREFIX) \
		--ldsc-bin $(LDSC_RS_BIN) \
		--out-dir $(ALL_RG_DIR) \
		--trait-block-size $(TRAIT_BLOCK_SIZE) \
		--prepare-only

all-rg: setup-ldsc-rs-rg-batch
	$(PYTHON) scripts/run_all_rg_hybrid.py \
		--manifest $(CATALOG) \
		--sumstats-dir $(SUMSTATS_DIR) \
		--ld-prefix $(LD_PREFIX) \
		--ldsc-bin $(LDSC_RS_BIN) \
		--out-dir $(ALL_RG_DIR) \
		--trait-block-size $(TRAIT_BLOCK_SIZE) \
		--rayon-threads $(RAYON_THREADS) \
		--max-parallel-shards $(MAX_PARALLEL_SHARDS)

all-rg-dry-run: setup-ldsc-rs-rg-batch
	$(PYTHON) scripts/run_all_rg_hybrid.py \
		--manifest $(CATALOG) \
		--sumstats-dir $(SUMSTATS_DIR) \
		--ld-prefix $(LD_PREFIX) \
		--ldsc-bin $(LDSC_RS_BIN) \
		--out-dir $(ALL_RG_DIR) \
		--trait-block-size $(TRAIT_BLOCK_SIZE) \
		--rayon-threads $(RAYON_THREADS) \
		--max-parallel-shards $(MAX_PARALLEL_SHARDS) \
		--dry-run

all-rg-progress:
	$(PYTHON) scripts/run_all_rg_hybrid.py \
		--out-dir $(ALL_RG_DIR) \
		--progress

all-rg-collect:
	$(PYTHON) scripts/run_all_rg_hybrid.py \
		--out-dir $(ALL_RG_DIR) \
		--collect

all-rg-validation: setup-ldsc-rs-rg-batch
	$(PYTHON) scripts/run_all_rg_validation.py \
		--manifest $(CATALOG) \
		--sumstats-dir $(SUMSTATS_DIR) \
		--ld-prefix $(LD_PREFIX) \
		--ldsc-bin $(LDSC_RS_BIN) \
		--out-dir $(ALL_RG_VALIDATION_DIR) \
		--n-pairs $(VALIDATION_PAIRS) \
		--seed $(VALIDATION_SEED) \
		--jobs $(VALIDATION_JOBS) \
		--hybrid-threads $(VALIDATION_RAYON_THREADS)

run-benchmark: setup-ldsc setup-ldsc-env prepare-sumstats
	mkdir -p results/benchmark90/rg logs/ldsc
	$(PYTHON) scripts/run_ldsc_triangle.py \
		--phenotypes $(BENCH_PHENOS) \
		--sumstats-dir data/benchmark90/sumstats \
		--ldsc-dir $(LDSC_DIR) \
		--ldsc-python "$(LDSC_PYTHON)" \
		--ld-prefix $(LD_PREFIX) \
		--out-dir results/benchmark90/rg \
		--jobs $(JOBS)

summarize:
	mkdir -p benchmarks/benchmark90
	$(PYTHON) scripts/summarize_benchmark.py \
		--phenotypes $(BENCH_PHENOS) \
		--prepare-stats results/benchmark90/prepare_stats \
		--rg-dir results/benchmark90/rg \
		--out-md benchmarks/benchmark90/summary.md \
		--out-tsv benchmarks/benchmark90/summary.tsv

hardware:
	mkdir -p results/benchmark90
	bash scripts/hardware_report.sh > results/benchmark90/hardware.txt

benchmark90: fetch-manifests select-benchmark run-benchmark summarize

clean-small:
	rm -rf benchmarks/benchmark90/*.tmp logs/*.tmp
