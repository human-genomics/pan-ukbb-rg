# Pan-UKBB EUR Genetic Correlations

Compute genome-wide genetic correlations (r<sub>g</sub>) across the Pan-UK Biobank
European-ancestry GWAS — one trait against all others, or the full all-pairs matrix —
and align external published GWAS to the same SNP set. The results power the
**[Genetic Correlation Catalog](https://rgcatalog.org)**.

**[🔎 Explore the results](https://rgcatalog.org)**  ·  **[⬇ Download the full export (v2.0)](https://github.com/jesseICR/pan-ukbb-rg/releases/tag/v2.0)**  ·  MIT-licensed

---

## What it does

Starting from public Pan-UKBB summary statistics, the pipeline:

1. **Prepares** every Pan-UKBB EUR GWAS (7,160 traits) as compact LDSC `.sumstats.gz` files (`make setup`).
2. **Computes** genetic correlations in any of three modes:
   - **One vs all** — one lead GWAS against all 7,159 others (Neale Lab LDSC fork, chunked).
   - **All pairs** — the full 25.6M-pair Pan-UKBB matrix (fast Rust `ldsc-rs` engine).
   - **External GWAS** — download and align ~270 published GWAS (PGC, GIANT, ENIGMA, EGG, …), then correlate them against Pan-UKBB and each other.

The published catalog combines these into **~28 million genetic correlations across 7,521 GWAS** (Pan-UKBB + external + Neale Lab sex-specific). No AWS account is required — Pan-UKBB data are streamed from public HTTPS URLs.

## Just want the data?

- **Browse interactively:** <https://rgcatalog.org>
- **Download the full r<sub>g</sub> export:** [release v2.0](https://github.com/jesseICR/pan-ukbb-rg/releases/tag/v2.0)

You only need this repository if you want to recompute or extend the correlations.

## Quick start

### Docker (recommended)

```bash
# 1. One-time setup into a persistent host directory (downloads + prepares the cache).
mkdir -p pan-ukbb-rg-work
docker run --rm -v "$(pwd)/pan-ukbb-rg-work:/app/pipeline-output" \
  ghcr.io/jesseicr/pan-ukbb-rg:latest setup --jobs 8

# 2. All correlations for one trait (phenocode 20016 = fluid intelligence score).
docker run --rm -v "$(pwd)/pan-ukbb-rg-work:/app/pipeline-output" \
  ghcr.io/jesseicr/pan-ukbb-rg:latest one-vs-all --phenocode 20016 --jobs 8
```

### Local

```bash
git clone https://github.com/jesseICR/pan-ukbb-rg.git
cd pan-ukbb-rg
make setup                                    # one-time data prep (long — see Runtime)
python3 scripts/run_one_vs_all.py --phenocode 20016
```

Select a trait by phenocode, phenotype ID, or free text:

```bash
python3 scripts/run_one_vs_all.py --phenotype-id continuous-20016-both_sexes-irnt
python3 scripts/run_one_vs_all.py --query "fluid intelligence"
```

`make setup` is fully **resumable** — rerun it after any interruption and completed traits are skipped.

## Workflows

### 1. One trait vs all

`make setup`, then `run_one_vs_all.py` (above).
**Output:** `results/one_vs_all/<phenotype_id>/rg.tsv` — the lead GWAS against every other EUR GWAS.

### 2. All pairs (Rust engine)

The full Pan-UKBB matrix (25,629,220 pairs), using the patched `ldsc-rs` `rg-batch` engine on the same `make setup` cache:

```bash
make all-rg-dry-run RAYON_THREADS=50 MAX_PARALLEL_SHARDS=1
make all-rg         RAYON_THREADS=50 MAX_PARALLEL_SHARDS=1
make all-rg-collect          # → results/all_rg/rg.tsv.gz
```

Traits are split into blocks and computed as resumable upper-triangle block-pair shards, so every pair is computed exactly once and runs can be stopped and resumed. Design details: [`docs/design.md`](docs/design.md).

### 3. External / published GWAS

```bash
make external-gwas           # download + align ~270 external GWAS to the Pan-UKBB SNP array
make external-rg RAYON_THREADS=9 MAX_PARALLEL_SHARDS=1   # external × Pan-UKBB and external × external
make external-rg-collect
```

The target list lives in `config/external_gwas_targets.tsv` (GWAS Catalog, Zenodo S-LDSC `indep107`, PGC, CNCR/CTG, EGG, GIANT, Levey Lab, and more). Restricted sources (ENIGMA, SSGAC/EA4) are tracked but their URLs are not stored in git — see [External data access](#external-data-access--terms).

## Requirements

- Unix-like shell with `bash`, `make`, `curl`, `gzip`, `awk`, `git`
- Python ≥ 3.10
- `micromamba` / `mamba` / `conda` — creates the isolated Python 2.7 LDSC environment
- ≥ 150 GiB free disk; ample RAM/CPU for parallel LDSC jobs
- Public internet access to Pan-UKBB

For the all-pairs and external Rust engine, additionally either **Rust/Cargo** (to build the `sharifhsn/ldsc` checkout) or **Docker** (used to build it).

```bash
make setup ENV_MANAGER=mamba   # or conda; default is micromamba
```

## Runtime & hardware

Reference machine: AMD Ryzen Threadripper PRO 5995WX (64 cores / 128 threads), 503 GiB RAM, NVMe.

| Step | Scope | Approx. wall time |
|---|---|---|
| `make setup` | stream ~6.7 TiB, prepare 7,160 sumstats | ~36–37 h (`JOBS=16`) |
| One vs all | 1 trait × 7,159 | ~33 min (`JOBS=16`) |
| All pairs | 25.6M pairs | ~10–16 days (~50 threads) |

GPUs don't help — LDSC is CPU-bound and setup is I/O-bound. A small `make benchmark90` is included to sanity-check LDSC on a new machine.

## Outputs

| Path | What it is |
|---|---|
| `data/catalog/eur_gwas_manifest.tsv` | searchable EUR GWAS catalog (phenocodes, IDs, QC, N, source) |
| `data/sumstats/eur/<id>.sumstats.gz` | LDSC-ready per-trait sumstats |
| `results/one_vs_all/<id>/rg.tsv` | one-vs-all correlations |
| `results/all_rg/rg.tsv.gz` | combined all-pairs matrix |
| `results/external_rg/` | external-GWAS correlations |

Large inputs and generated outputs under `data/`, `external/`, `results/`, `logs/`, and `.envs/` are gitignored.

**Sumstats conversion.** Pan-UKBB flat files are converted to minimal LDSC sumstats:
`SNP = chr:pos:ref:alt`, `A1 = alt`, `A2 = ref`, `Z = beta_EUR / se_EUR`,
`N = cases + controls` (binary) or `cases` (continuous). Rows outside the EUR
LD-score SNP set, flagged `low_confidence_EUR`, or with missing/invalid statistics are dropped.

## External data access & terms

Running the external download step makes **you** responsible for complying with each source's terms (GWAS Catalog, Zenodo, **PGC**, Figshare, and any linked source). PGC summary statistics are restricted: running `make external-gwas` means you accept the PGC Data Access Terms (research use only, no reposting, appropriate citation, etc.).

Restricted GWAS with no public URL (ENIGMA cortical/subcortical, SSGAC EA4) are tracked in `config/external_gwas_targets.tsv` but their URLs are not committed. After obtaining access, copy `config/external_gwas_restricted_urls.tsv.example` → `config/external_gwas_restricted_urls.tsv` and add local URLs or `file://` paths; that override file is gitignored.

## Credits & citation

This pipeline stands on the work of others:

- **ldsc-rs** — Haason, S. & Khan, Y. (2026). *ldsc-rs: Exact and approximate LD Score Regression at biobank scale.* bioRxiv. <https://github.com/sharifhsn/ldsc>
- **LD Score Regression** — Bulik-Sullivan et al. (2015), *Nature Genetics* 47:291–295 (heritability) and 47:1236–1241 (genetic correlation); run here via the [Neale Lab LDSC fork](https://github.com/astheeggeggs/ldsc).
- **Pan-UK Biobank** and the **Neale Lab** UK Biobank GWAS, plus external consortia (PGC, GIANT, SSGAC, ENIGMA, EGG, …). Please also cite the primary GWAS for any specific trait you rely on.

If you use the catalog or this pipeline, please credit the above and link <https://rgcatalog.org>.

## License

This repository's code is released under the **MIT License** (see [`LICENSE`](LICENSE)). It builds on GPL-3.0 LDSC engines — the Neale Lab fork and `ldsc-rs` — which are fetched and built at setup time rather than redistributed here.
