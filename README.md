# Pan-UKBB EUR Genetic Correlations

This repository prepares public Pan-UKBB EUR GWAS summary statistics for LDSC and computes genetic correlations either between one selected GWAS and every other EUR GWAS, or across all EUR GWAS pairs.

The workflow uses the Neale Lab batching idea: convert Pan-UKBB flat files into compact LDSC `.sumstats.gz` files, then run the patched Neale LDSC fork with `--rg-file` and `--write-rg` so one lead phenotype can be compared against target chunks in parallel.

The optional all-pairs workflow reuses the same `make setup` data cache but runs a patched Rust LDSC engine with an `rg-batch` command. It computes the same pair-specific LDSC rg model, but organizes work into resumable block-pair shards for throughput.

## Quick Start

### Docker

The Docker image bundles this repo, Python 3, the pinned Neale LDSC checkout,
and the Python 2.7 LDSC environment. It does not include the large Pan-UKBB
data cache; generated files are written to the mounted directory.

Pull the pre-built image:

```bash
docker pull ghcr.io/jesseicr/pan-ukbb-rg:latest
```

Run the one-time setup into a persistent host directory:

```bash
mkdir -p pan-ukbb-rg-work
docker run --rm \
  -v "$(pwd)/pan-ukbb-rg-work:/app/pipeline-output" \
  ghcr.io/jesseicr/pan-ukbb-rg:latest setup --jobs 8
```

Then compute all genetic correlations involving one GWAS:

```bash
docker run --rm \
  -v "$(pwd)/pan-ukbb-rg-work:/app/pipeline-output" \
  ghcr.io/jesseicr/pan-ukbb-rg:latest one-vs-all --phenocode 20016 --jobs 8
```

Dry-run phenotype resolution and cache completeness checks:

```bash
docker run --rm \
  -v "$(pwd)/pan-ukbb-rg-work:/app/pipeline-output" \
  ghcr.io/jesseicr/pan-ukbb-rg:latest dry-run --phenocode 20016
```

Runtime state is written under:

```text
pan-ukbb-rg-work/data/
pan-ukbb-rg-work/results/
pan-ukbb-rg-work/logs/
```

Build locally:

```bash
docker build -t pan-ukbb-rg .
docker run --rm -v "$(pwd)/pan-ukbb-rg-work:/app/pipeline-output" pan-ukbb-rg help
```

### Local

Clone the repo, then run the one-time setup:

```bash
git clone https://github.com/jesseICR/pan-ukbb-rg.git
cd pan-ukbb-rg
```

Run setup:

```bash
make setup
```

Optionally, after `make setup`, download and align external GWAS summary statistics
to the same Pan-UKBB EUR SNP array:

```bash
make external-gwas
```

This writes LDSC-ready files under:

```text
data/external_gwas/aligned/*.sumstats.gz
data/external_gwas/catalog/external_gwas_manifest.tsv
data/external_gwas/prepare_stats/*.json
```

The external target list is `config/external_gwas_targets.tsv`. It currently covers
the requested GWAS Catalog education traits, the Zenodo S-LDSC `sumstats_indep107`
traits, static Figshare download URLs for the best/latest PGC row for each
requested unique PGC phenotype, explicit CNCR/CTG Nextcloud rows, the PGC
Alzheimer 2021 CNCR file
`PGCALZ2sumstatsExcluding23andMe.txt.gz`, the CNCR Jansen/Nagel 2020 brain-volume
files, Tielbeek et al. 2022 BroadABC summary statistics, static EGG Consortium
releases, Levey Lab public EUR files, GIANT anthropometric releases, All of Us
EUR height/BMI summary statistics, and a small number of direct public URLs for
requested traits.

PGC rows are limited to one best/latest row for each unique phenotype value in the
PGC table, except where the table exposes distinct phenotypes within a working
group such as panic disorder, quantitative anxiety symptoms, case-control anxiety,
obsessive-compulsive symptoms, hoarding symptoms, and the specific substance-use
phenotypes. The `cdg2025` PGC release is expanded into its six factor GWAS files.
`alz2021` is included via the CNCR/CTG Nextcloud file approved by the PGC download
page instead of as a duplicate unsupported PGC table row. Maciel et al. 2026
`META.txt.gz` is excluded because it has SNP annotation columns but no association
`Z`, `BETA`, `SE`, or `P` columns. The 2024 EGG trans-ancestry pubertal-growth
rows are also excluded from the default target list. `PASS.Multiple_Sclerosis.IMSGC2019`
is still recorded as missing because the Zenodo readme states that file requires
separate approval; the public older MS GWAS Catalog row `GCST001198` is included
as the automatic fallback.

External GWAS not auto-added:

- `Risktaking-Auto speeding` and `Risktaking-Nr. Sexual Partners`: the Karlsson
  Linner et al. 2019 paper points to the SSGAC data portal, and the current
  portal requires an account; OpenGWAS mirrors for the close UKB traits require
  generated JWT/signed links for automated download.
- IMSGC 2019 multiple sclerosis (`GCST009597` / `ieu-b-18`): GWAS Catalog marks
  full summary statistics unavailable, and IMSGC requires a data-access request.
- PGC suicide attempt (`sui2023`): the PGC download table points to a data-access
  form rather than a public Figshare/Nextcloud file, so it is tracked as a
  non-automated source.
- Maciel et al. 2026 `META.txt.gz`: metadata only; not an RG input.
- 2024 EGG trans-ancestry pubertal-growth rows: excluded by ancestry policy for
  this external-GWAS run.

External GWAS to keep/add:

- `PGC.alz2021_PGCALZ2_no23andMe`
- `broadabc2022_combined`, `broadabc2022_females`, `broadabc2022_males`
- `GCST001198` public multiple sclerosis fallback
- `PASS.General_Risk_Tolerance.KarlssonLinner2019`, already present from the
  Zenodo S-LDSC archive, for the SSGAC risk-tolerance/risk-taking signal.

The aligned files use the Pan-UKBB SNP identifier (`chr:pos:ref:alt`), Pan-UKBB
ALT as `A1`, Pan-UKBB REF as `A2`, and flip `Z` when the external effect allele is
the Pan-UKBB REF allele. Coordinate/allele matching is preferred. rsID+allele
matching is used as a fallback when a source only has rsIDs or when it helps with
legacy formats. Sources are restricted to European/EUR files where the source
offers population-specific downloads; AFR, EAS, SAS, and other non-EUR files are
excluded. Explicit EGG trans-ancestry birth-weight and childhood-obesity rows are
retained, while the 2024 EGG trans-ancestry pubertal-growth rows are excluded.

After `make external-gwas`, compute all genetic correlations involving at least
one external GWAS:

```bash
make external-rg RAYON_THREADS=9 MAX_PARALLEL_SHARDS=1
```

This prepares a combined symlinked sumstats directory, writes pair shards only for
external x Pan-UKBB and external x external pairs, and runs the same patched Rust
`ldsc rg-batch` engine used by `make all-rg`. It deliberately does not recompute
Pan-UKBB x Pan-UKBB pairs. External rg defaults to
`EXTERNAL_RG_TRAIT_BLOCK_SIZE=64`, separately from the Pan-UKBB all-rg
`TRAIT_BLOCK_SIZE=256`, so each external shard stays small enough to run
alongside an existing Pan-UKBB rg job. Runtime state is written under:

```text
data/external_rg/sumstats/
results/external_rg/metadata/
results/external_rg/pair_shards/
results/external_rg/rg_shards/
```

If the Pan-UKBB cache is outside this checkout, point the target at it:

```bash
make external-rg \
  PAN_GWAS_MANIFEST=/path/to/data/catalog/eur_gwas_manifest.tsv \
  PAN_SUMSTATS_DIR=/path/to/data/sumstats/eur \
  LD_PREFIX=/path/to/data/ld/UKBB.EUR \
  LDSC_RS_BIN=/path/to/external/ldsc-rs-rg-batch-target/release/ldsc \
  RAYON_THREADS=9 \
  MAX_PARALLEL_SHARDS=1 \
  EXTERNAL_RG_TRAIT_BLOCK_SIZE=64
```

Progress and collection commands:

```bash
make external-rg-progress
make external-rg-collect
```

Useful variants:

```bash
make external-gwas-manifest
make external-gwas EXTERNAL_GWAS_SOURCES=gwas_catalog
make external-gwas EXTERNAL_GWAS_SOURCES=cncr EXTERNAL_GWAS_LIMIT=20
make external-gwas EXTERNAL_GWAS_SOURCES=egg EXTERNAL_GWAS_LIMIT=20
make external-gwas EXTERNAL_GWAS_SOURCES=zenodo_face_cgwas EXTERNAL_GWAS_LIMIT=1 EXTERNAL_GWAS_MAX_ROWS=200000
make external-gwas EXTERNAL_GWAS_INCLUDE=GCST90011874
make external-gwas-smoke
make external-rg-prepare
make external-rg-dry-run
```

The default `make external-gwas` path uses explicit rows in
`config/external_gwas_targets.tsv`. The `pgc`, `cncr`, and `egg` source modes are
kept as discovery tools for refreshing the static target list, not as the default
runtime path. The optional `zenodo_face_cgwas` source supports the 2024
facial-shape Zenodo record by downloading `SnpInfo.tsv.bz2`, joining it to the
946 `Beta/P` files inside the `.tar.bz2` archives, and aligning each facial
distance GWAS to Pan-UKBB SNPs. It is not part of the default source list because
the raw record is about 23.5 GB and the current converter processes those 946
phenotypes serially.

Restricted ENIGMA rows are tracked in `config/external_gwas_targets.tsv`, but
their direct download URLs are intentionally not stored in git. Users can request
ENIGMA GWAS access at <https://enigma.ini.usc.edu/>. After access is granted,
copy `config/external_gwas_restricted_urls.tsv.example` to
`config/external_gwas_restricted_urls.tsv` and fill `source_url` locally. That
override file is gitignored and is read by `make external-gwas` through
`EXTERNAL_GWAS_URL_OVERRIDES`. Without the local override file, those ENIGMA rows
are reported as `restricted_url_missing` instead of downloaded.

SSGAC rows that require an account use the same pattern. The tracked target rows
describe the GWAS but omit direct URLs; users should request access through
<https://www.thessgac.org/data>, then add local URLs or `file://` paths to
`config/external_gwas_restricted_urls.tsv`. The local override file on this
machine points `SSGAC.EA4_additive_excl_23andMe` at
`/home/jesse/Downloads/EA4_additive_excl_23andMe.txt.gz`.

If you already have the large Zenodo files locally, point the pipeline at them:

```bash
make external-gwas \
  EXTERNAL_GWAS_SOURCES=zenodo_indep107 \
  EXTERNAL_GWAS_ZENODO_ARCHIVE=/path/to/sumstats_indep107.tgz \
  EXTERNAL_GWAS_RSID_REFERENCE=/path/to/reference.1000G.maf.0.005.txt.gz
```

By running the external download step, you are responsible for complying with
the access terms of GWAS Catalog, Zenodo, PGC, Figshare, and any linked source.
For PGC summary statistics specifically, running `make external-gwas` means you
agree to the PGC Data Access Terms and Conditions published on the PGC download
page: the data are not unrestricted, are provided as-is, must not be reposted,
must not be used to identify participants, are for scientific research only
unless separate commercial permission is obtained, must be used in compliance
with applicable policies, must not be used for prenatal risk or predictive
testing, require appropriate PGC citation, and require respect for PGC scientific
priorities when data are released before publication.

Then compute all genetic correlations involving one GWAS:

```bash
python3 scripts/run_one_vs_all.py --phenocode 20016
```

`20016` is the Pan-UKBB phenocode for fluid intelligence score. Replace it with any phenocode in `data/catalog/eur_gwas_manifest.tsv`.

You can also select by phenotype ID or text query:

```bash
python3 scripts/run_one_vs_all.py --phenotype-id continuous-20016-both_sexes-irnt
python3 scripts/run_one_vs_all.py --query "fluid intelligence"
```

If a selector matches multiple GWAS, the script prints the matching rows and asks you to use a more specific selector.

To compute genetic correlations for **all EUR GWAS pairs** after the same `make setup` step, use the Rust hybrid all-pairs add-on:

```bash
make all-rg-dry-run RAYON_THREADS=50 MAX_PARALLEL_SHARDS=1
make all-rg RAYON_THREADS=50 MAX_PARALLEL_SHARDS=1
make all-rg-collect
```

This does not change the existing one-vs-all workflow. It consumes the same prepared files from `make setup`:

```text
data/catalog/eur_gwas_manifest.tsv
data/sumstats/eur/*.sumstats.gz
data/ld/UKBB.EUR.l2.ldscore.gz
data/ld/UKBB.EUR.l2.M_5_50
```

The all-pairs outputs are written under `results/all_rg/`. Progress and final collection:

```bash
make all-rg-progress
make all-rg-collect
```

To validate the Rust batch path against the existing Rust `rg` path on 100 random disjoint trait pairs:

```bash
make all-rg-validation
```

## Requirements

You need:

- a Unix-like shell with `bash` and `make`
- `curl`, `gzip`, `awk`, and `git`
- Python 3.10 or newer for this repo's scripts
- `micromamba`, `mamba`, or `conda` to create the isolated LDSC Python 2.7 environment
- public internet access to Pan-UKBB HTTPS URLs
- at least 150 GiB free disk for the prepared data cache
- enough RAM/CPU for parallel LDSC jobs

For the optional all-pairs Rust hybrid add-on, you also need one of:

- Rust/Cargo capable of building the patched `sharifhsn/ldsc` checkout, or
- Docker, which the setup script uses to build the patched Rust binary via `rust:1.91-bookworm`

The all-pairs add-on writes many shard outputs and the final combined rg table. Plan for additional result-space beyond the setup cache; compressed shard outputs for 25.6M pairs are expected to be several GiB, depending on formatting and compression.

No AWS account, AWS credentials, or AWS CLI are required. Pan-UKBB data are downloaded from public HTTPS URLs such as:

```text
https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz
```

By default, `make setup` uses `micromamba`. If you use another environment manager:

```bash
make setup ENV_MANAGER=mamba
make setup ENV_MANAGER=conda
```

The LDSC environment is created under `.envs/ldsc-neale/`, which is local generated state and is not meant to be committed.

## Runtime and Hardware

The default is `JOBS=16` for both setup and one-vs-all rg. `JOBS` means concurrent download/conversion jobs during setup or concurrent LDSC chunk jobs during one-vs-all analysis. Use fewer jobs if the machine becomes unresponsive or if network/disk throughput is saturated:

```bash
make setup JOBS=8
python3 scripts/run_one_vs_all.py --phenocode 20016 --jobs 8
```

The full run used `JOBS=16` on a Linux workstation with an AMD Ryzen Threadripper PRO 5995WX, 64 physical cores / 128 threads, 503 GiB RAM, and NVMe storage. Runtime on other systems will mainly depend on internet throughput, remote S3 throttling, gzip decompression, disk write speed, and available RAM for parallel LDSC jobs.

The runtimes below are approximate observed wall times from the completed run artifacts and per-GWAS conversion metadata. They are meant as planning guidance, not a guarantee for every network or storage environment.

Observed setup runtime:

| Step | What It Does | Observed Runtime |
|------|--------------|------------------|
| `setup-ldsc-env` | Create the isolated Python 2.7 LDSC environment | Short setup step; exact wall time was not recoverable from preserved artifacts |
| `validate-catalog` | Download Pan-UKBB manifests, build `eur_gwas_manifest.tsv`, and validate 7,160 EUR GWAS rows | Less than 1 minute when rerun locally |
| `prepare-ldscores` | Download Pan-UKBB EUR LD-score reference files and extract `UKBB.EUR.snps` | Less than 1 minute |
| `setup-ldsc` | Clone the patched Neale LDSC fork and check out the pinned commit | Less than 1 minute |
| `prepare-all-sumstats` | Stream every EUR GWAS, filter to EUR LD-score SNPs, and write LDSC-ready `.sumstats.gz` files | About 36.3 hours at `JOBS=16` |
| Full `make setup` | All setup steps above | About 36-37 hours, dominated by `prepare-all-sumstats` |

The full setup streamed roughly 6.7 TiB from Pan-UKBB and wrote 7,160 prepared sumstats files. The per-GWAS conversion diagnostics in `results/prepare_all/prepare_stats/` recorded about 598 job-hours of conversion work. At `JOBS=16`, the observed wall time for the conversion step was about 36.3 hours. Median per-GWAS conversion time was about 2.8 minutes, with the slowest GWAS taking about 18.9 minutes.

Observed one-vs-all runtime:

| Analysis | Targets | Chunks / Jobs | Observed Runtime |
|----------|---------|---------------|------------------|
| Fluid intelligence score, phenocode `20016` | 7,159 EUR GWAS | 16 chunks at `JOBS=16` | About 33 minutes |

For the fluid-intelligence example, all 16 LDSC chunks completed successfully. Individual chunk runtimes ranged from about 29 to 33 minutes. The combined output is `results/one_vs_all/continuous-20016-both_sexes-irnt/rg.tsv`.

Observed all-pairs Rust hybrid benchmark:

| Engine | Pairs | Threads / Jobs | Observed Runtime |
|--------|-------|----------------|------------------|
| Original Neale/Python LDSC | 100 | 8 workers | 161.54 seconds |
| Existing Rust grouped rg | 100 | 8 workers | 48.77 seconds |
| New Rust `rg-batch` hybrid | 100 | 8 Rayon threads | 25.74 seconds |
| New Rust `rg-batch` hybrid | 100 | 16 Rayon threads | 23.41 seconds |

The all-pairs add-on computes all `7160 * 7159 / 2 = 25,629,220` EUR GWAS pairs. On a 50-core machine similar to the Threadripper Pro system above, the current practical estimate is roughly **10-16 days** for the rg step after `make setup` has completed, with **about two weeks** as the planning number. The estimate is based on the 100-pair hybrid benchmark and validation run, then scaled to 50 useful compute threads. It may vary with memory bandwidth, allocator overhead, disk speed, and how well the machine handles either one 50-thread shard or several smaller shards.

Start with:

```bash
make all-rg RAYON_THREADS=50 MAX_PARALLEL_SHARDS=1
```

If CPU utilization is poor or the machine is more stable with smaller jobs, try the same total core budget split across shards:

```bash
make all-rg RAYON_THREADS=10 MAX_PARALLEL_SHARDS=5
```

High-level all-pairs design:

1. Traits are split into blocks, defaulting to 256 traits per block.
2. The runner creates one shard for each upper-triangle block pair: block 0 vs block 0, block 0 vs block 1, block 0 vs block 2, and so on.
3. A diagonal shard computes only within-block pairs where `i < j`.
4. An off-diagonal shard computes every pair between the two different blocks.
5. Because only the upper triangle is used, every GWAS pair is computed exactly once. The workflow does not compute A-vs-B and then B-vs-A.
6. For each shard, `rg-batch` loads LD scores once and loads the unique traits in that shard once. With the default block size, a shard loads at most 512 traits.
7. Within the shard, `rg-batch` computes pair-specific LDSC h2/gencov/rg for every pair, preserving pair-specific SNP intersections and allele checks.

So it is not literally running 7,160 one-vs-all jobs. It is closer to running many block-vs-block batches. This captures the useful idea behind one-vs-all batching, but avoids duplicate pairs and keeps each job resumable.

From the 90-phenotype benchmark in this repo:

| Metric | Value |
|--------|-------|
| Benchmark phenotypes | 90 |
| Genetic-correlation pairs | 4,005 |
| LDSC wall time | About 1.04 hours |
| Sum of LDSC worker time | About 16.4 job-hours |
| Effective runtime per pair | About 14.8 seconds |
| Sumstats preparation worker time | About 5.5 job-hours |

The benchmark is useful for checking that LDSC runs correctly on a new machine, but the full setup is more network- and disk-sensitive because it streams thousands of large Pan-UKBB files.

GPUs are not useful for this workflow. LDSC is CPU-bound Python, and setup is mostly streaming, decompression, filtering, and writing.

## Is Setup Resumable?

Yes. `make setup` is designed to be rerun.

The large setup step is `make prepare-all-sumstats`. For each phenotype, the converter writes temporary files first and only marks that phenotype complete after both of these final files exist:

```text
data/sumstats/eur/<phenotype_id>.sumstats.gz
results/prepare_all/prepare_stats/<phenotype_id>.json
```

If the computer shuts down mid-run, rerun:

```bash
make setup
```

Completed phenotypes are skipped. Interrupted or incomplete phenotypes are streamed and converted again. The raw Pan-UKBB flat files are not stored locally.

## Important Files

After `make setup`, the most useful files are:

```text
data/catalog/eur_gwas_manifest.tsv
```

The searchable EUR GWAS catalog. Use this to find phenocodes, phenotype IDs, descriptions, QC labels, sample sizes, and source URLs.

```text
data/sumstats/eur/<phenotype_id>.sumstats.gz
```

The compact LDSC-ready EUR sumstats cache. These are the files used by `run_one_vs_all.py`.

```text
results/prepare_all/prepare_stats/<phenotype_id>.json
```

Per-GWAS conversion diagnostics: input rows, LD-score SNPs seen, written rows, skipped low-confidence rows, skipped missing EUR stats, and elapsed time.

```text
data/ld/UKBB.EUR.l2.ldscore.gz
data/ld/UKBB.EUR.l2.M
data/ld/UKBB.EUR.l2.M_5_50
data/ld/UKBB.EUR.snps
```

EUR LD-score reference files used by LDSC. `UKBB.EUR.snps` is the extracted list of SNP IDs used to filter the Pan-UKBB flat files.

```text
external/ldsc-neale/
.envs/ldsc-neale/
```

The patched LDSC code and its isolated Python 2.7 environment.

After `run_one_vs_all.py`, the main output is:

```text
results/one_vs_all/<lead_phenotype_id>/rg.tsv
```

This is the combined table of genetic correlations between the selected lead GWAS and all other EUR GWAS.
The `p1` and `p2` columns are phenotype IDs from `data/catalog/eur_gwas_manifest.tsv`.

Other useful one-vs-all outputs:

```text
results/one_vs_all/<lead_phenotype_id>/lead.tsv
results/one_vs_all/<lead_phenotype_id>/driver_summary.tsv
results/one_vs_all/<lead_phenotype_id>/logs/
results/one_vs_all/<lead_phenotype_id>/chunks/
```

`lead.tsv` records the selected GWAS metadata. `driver_summary.tsv` records per-chunk runtime/status. `logs/` is for troubleshooting. `chunks/` contains the raw Neale LDSC `.r2` chunk outputs that are combined into `rg.tsv`.

After `make all-rg`, the all-pairs outputs are:

```text
results/all_rg/metadata/traits.tsv
results/all_rg/metadata/shards.tsv
results/all_rg/pair_shards/
results/all_rg/rg_shards/
results/all_rg/logs/
```

`metadata/shards.tsv` records every block-pair shard. `rg_shards/` contains compressed per-shard rg outputs plus `.done` markers. Rerunning `make all-rg` skips shards whose output exists and whose row count matches the expected pair count. `make all-rg-collect` concatenates completed shards into:

```text
results/all_rg/rg.tsv.gz
```

After `make all-rg-validation`, the validation outputs are:

```text
results/all_rg_validation/random_disjoint_pairs.tsv
results/all_rg_validation/hybrid_rg.tsv
results/all_rg_validation/baseline_rg.tsv
results/all_rg_validation/comparison_summary.json
```

The validation chooses 100 random pairs using 200 distinct traits when possible. It compares the new Rust `rg-batch` path against the existing Rust `rg` command, which exercises the same LDSC regression code through the older one-pair/one-lead interface.

Large public inputs and generated outputs under `data/`, `external/`, `results/`, `logs/`, and `.envs/` are gitignored.

## Dry Runs

A dry run is optional. It does not stream GWAS data and does not run LDSC. It only resolves the phenotype and reports whether the prepared sumstats cache exists.

```bash
python3 scripts/run_one_vs_all.py --phenocode 20016 --dry-run
```

This is useful before launching a long one-vs-all run or when checking whether `make setup` completed.

## What Setup Does

`make setup` runs these steps:

```text
make setup-ldsc-env
make validate-catalog
make prepare-ldscores
make setup-ldsc
make prepare-all-sumstats
```

`make setup-ldsc-env` creates `.envs/ldsc-neale/` from `envs/ldsc-neale.yml`.

`make validate-catalog` downloads the public Pan-UKBB manifests, builds `data/catalog/eur_gwas_manifest.tsv`, and validates it against the source manifest.

`make prepare-ldscores` downloads the Pan-UKBB EUR LD-score files and extracts `data/ld/UKBB.EUR.snps`.

`make setup-ldsc` checks out the patched Neale LDSC fork from `https://github.com/astheeggeggs/ldsc.git` at commit `a4ee4c8aa065a1c9a586c3b678e9b3040bbebafc`.

`make prepare-all-sumstats` streams every Pan-UKBB EUR GWAS flat file and writes compact LDSC-ready sumstats.

The optional all-pairs add-on is intentionally separate from `make setup`. After setup has produced the shared data cache, `make all-rg` checks out and patches the Rust LDSC rewrite under `external/ldsc-rs-rg-batch/`, builds it under `external/ldsc-rs-rg-batch-target/`, creates block-pair shards, and computes all pairwise rg results.

Expected data volume:

```text
streamed from Pan-UKBB: ~6.7 TiB
stored after setup: ~70-100 GiB in data/sumstats/eur/
recommended free disk: >=150 GiB
```

## Sumstats Conversion

The Pan-UKBB per-phenotype flat files are not pre-munged LDSC files. This repo converts EUR columns to minimal LDSC-compatible files with:

- `SNP = chr:pos:ref:alt`
- `A1 = alt`
- `A2 = ref`
- `Z = beta_EUR / se_EUR`
- `N = n_cases_EUR + n_controls_EUR` for binary traits, or `n_cases_EUR` otherwise

Rows are skipped if they are outside the EUR LD-score SNP list, have `low_confidence_EUR == true`, have missing EUR beta/se, or have nonpositive/nonfinite standard error.

## Optional Benchmark

The benchmark is separate from the main one-vs-all workflow. It estimates local LDSC runtime using 90 selected phenotypes: 30 `PASS`, 30 `h2_z_insignificant`, and 30 `not_EUR_plus_1`.

```bash
make benchmark90
```

Benchmark outputs are written under:

```text
benchmarks/benchmark90/
results/benchmark90/
```
