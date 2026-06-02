#!/usr/bin/env python3
"""Run external GWAS genetic correlations against Pan-UKBB and each other."""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

import run_all_rg_hybrid as rg


def symlink_sumstats(src: Path, dst: Path, force: bool = False) -> None:
    src = src.resolve()
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        if dst.is_symlink() and Path(os.readlink(dst)).resolve() == src:
            return
        if not force:
            raise SystemExit(f"Refusing to replace existing sumstats path: {dst}")
        dst.unlink()
    os.symlink(src, dst)


def read_pan_traits(manifest: Path, sumstats_dir: Path, allow_missing: bool) -> list[dict[str, str]]:
    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "phenotype_id" not in (reader.fieldnames or []):
            raise SystemExit(f"{manifest} does not contain a phenotype_id column")
        rows = [row for row in reader if row.get("phenotype_id")]

    seen: set[str] = set()
    traits: list[dict[str, str]] = []
    missing: list[str] = []
    for row in rows:
        trait_id = row["phenotype_id"]
        if trait_id in seen:
            continue
        seen.add(trait_id)
        path = sumstats_dir / f"{trait_id}.sumstats.gz"
        if not path.exists():
            missing.append(trait_id)
            continue
        traits.append(
            {
                "trait_id": trait_id,
                "source_group": "panukbb",
                "source_path": str(path.resolve()),
                "description": row.get("description", ""),
                "category": row.get("category", ""),
            }
        )
    if missing and not allow_missing:
        preview = ", ".join(missing[:10])
        raise SystemExit(
            f"{len(missing)} Pan-UKBB traits are missing prepared sumstats. "
            f"First missing: {preview}"
        )
    return traits


def read_external_traits(manifest: Path, allow_missing: bool) -> list[dict[str, str]]:
    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        for required in ["external_id", "status", "aligned_path"]:
            if required not in fieldnames:
                raise SystemExit(f"{manifest} does not contain required column {required}")
        rows = list(reader)

    traits: list[dict[str, str]] = []
    missing: list[str] = []
    seen: set[str] = set()
    for row in rows:
        if row.get("status") != "done":
            continue
        trait_id = row["external_id"]
        if trait_id in seen:
            continue
        seen.add(trait_id)
        path = Path(row.get("aligned_path") or "")
        if not path.is_absolute():
            path = (manifest.parent / path).resolve() if path.parts else path
            if not path.exists():
                path = Path(row.get("aligned_path") or "")
        if not path.exists():
            missing.append(trait_id)
            continue
        traits.append(
            {
                "trait_id": trait_id,
                "source_group": "external",
                "source_path": str(path.resolve()),
                "description": row.get("trait", ""),
                "category": row.get("category", ""),
            }
        )
    if missing and not allow_missing:
        preview = ", ".join(missing[:10])
        raise SystemExit(
            f"{len(missing)} external GWAS are missing aligned sumstats. "
            f"First missing: {preview}"
        )
    return traits


def write_external_traits(path: Path, traits: list[dict[str, str]], sumstats_dir: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["trait_index", "trait_id", "source_group", "sumstats_path", "description", "category"])
        for idx, trait in enumerate(traits):
            writer.writerow(
                [
                    idx,
                    trait["trait_id"],
                    trait["source_group"],
                    sumstats_dir / f"{trait['trait_id']}.sumstats.gz",
                    trait.get("description", ""),
                    trait.get("category", ""),
                ]
            )
    tmp.replace(path)


def write_combined_manifest(path: Path, traits: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["phenotype_id", "source_group", "description", "category", "source_path"])
        for trait in traits:
            writer.writerow(
                [
                    trait["trait_id"],
                    trait["source_group"],
                    trait.get("description", ""),
                    trait.get("category", ""),
                    trait["source_path"],
                ]
            )
    tmp.replace(path)


def pair_matches_filter(
    p1: str,
    p2: str,
    include_ids: set[str],
    include_prefixes: list[str],
) -> bool:
    if not include_ids and not include_prefixes:
        return True
    if p1 in include_ids or p2 in include_ids:
        return True
    return any(p1.startswith(prefix) or p2.startswith(prefix) for prefix in include_prefixes)


def pair_rows_for_block(
    traits: list[dict[str, str]],
    n_pan: int,
    block_i: int,
    block_j: int,
    block_size: int,
    include_ids: set[str],
    include_prefixes: list[str],
) -> list[tuple[int, int, int, str, str, str, str]]:
    n_traits = len(traits)
    i0 = block_i * block_size
    i1 = min(i0 + block_size, n_traits)
    j0 = block_j * block_size
    j1 = min(j0 + block_size, n_traits)
    rows: list[tuple[int, int, int, str, str, str, str]] = []
    for i in range(i0, i1):
        start_j = max(i + 1, j0)
        for j in range(start_j, j1):
            if i < n_pan and j < n_pan:
                continue
            p1 = traits[i]["trait_id"]
            p2 = traits[j]["trait_id"]
            if not pair_matches_filter(p1, p2, include_ids, include_prefixes):
                continue
            rows.append((rg.pair_id(n_traits, i, j), i, j, p1, p2, p1, p2))
    return rows


def split_csv_values(values: list[str]) -> list[str]:
    out: list[str] = []
    for value in values:
        out.extend(part.strip() for part in value.split(",") if part.strip())
    return out


def prepare_external_shards(args: argparse.Namespace) -> list[dict[str, str]]:
    pan_traits = read_pan_traits(args.pan_manifest, args.pan_sumstats_dir, args.allow_missing_sumstats)
    external_traits = read_external_traits(args.external_manifest, args.allow_missing_sumstats)
    if not external_traits:
        raise SystemExit("No done external GWAS found")

    pan_ids = {trait["trait_id"] for trait in pan_traits}
    collisions = [trait["trait_id"] for trait in external_traits if trait["trait_id"] in pan_ids]
    if collisions:
        raise SystemExit(f"External IDs collide with Pan-UKBB phenotype IDs: {', '.join(collisions[:10])}")

    traits = pan_traits + external_traits
    for trait in traits:
        symlink_sumstats(
            Path(trait["source_path"]),
            args.combined_sumstats_dir / f"{trait['trait_id']}.sumstats.gz",
            args.force_symlinks,
        )

    meta_dir = args.out_dir / "metadata"
    pair_dir = args.out_dir / "pair_shards"
    meta_dir.mkdir(parents=True, exist_ok=True)
    pair_dir.mkdir(parents=True, exist_ok=True)
    write_external_traits(meta_dir / "traits.tsv", traits, args.combined_sumstats_dir)
    write_combined_manifest(meta_dir / "combined_manifest.tsv", traits)

    n_traits = len(traits)
    n_pan = len(pan_traits)
    n_external = len(external_traits)
    n_blocks = (n_traits + args.trait_block_size - 1) // args.trait_block_size
    include_ids = set(split_csv_values(args.pair_include_trait_id))
    include_prefixes = split_csv_values(args.pair_include_trait_prefix)
    manifest_rows: list[dict[str, str]] = []
    total_pairs = 0
    external_pan_pairs = 0
    external_external_pairs = 0
    shard_id = 0
    for block_i in range(n_blocks):
        for block_j in range(block_i, n_blocks):
            rows = pair_rows_for_block(
                traits,
                n_pan,
                block_i,
                block_j,
                args.trait_block_size,
                include_ids,
                include_prefixes,
            )
            if not rows:
                continue
            shard = pair_dir / f"pairs.block_{block_i:03d}_{block_j:03d}.tsv"
            if args.force_shards or not shard.exists():
                rg.write_shard(shard, rows)
            total_pairs += len(rows)
            external_pan_pairs += sum(1 for _, i, j, *_ in rows if i < n_pan <= j)
            external_external_pairs += sum(1 for _, i, j, *_ in rows if n_pan <= i and n_pan <= j)
            manifest_rows.append(
                {
                    "shard_id": str(shard_id),
                    "block_i": str(block_i),
                    "block_j": str(block_j),
                    "n_pairs": str(len(rows)),
                    "path": str(shard.resolve()),
                }
            )
            shard_id += 1

    if include_ids or include_prefixes:
        expected = total_pairs
    else:
        expected = n_pan * n_external + n_external * (n_external - 1) // 2
        if total_pairs != expected:
            raise SystemExit(f"internal error: planned {total_pairs} pairs, expected {expected}")

    shard_manifest = meta_dir / "shards.tsv"
    tmp = shard_manifest.with_suffix(shard_manifest.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["shard_id", "block_i", "block_j", "n_pairs", "path"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    tmp.replace(shard_manifest)

    summary = {
        "pan_traits": n_pan,
        "external_traits": n_external,
        "total_traits": n_traits,
        "external_pan_pairs": external_pan_pairs,
        "external_external_pairs": external_external_pairs,
        "total_pairs": total_pairs,
        "shards": len(manifest_rows),
        "trait_block_size": args.trait_block_size,
        "combined_sumstats_dir": str(args.combined_sumstats_dir),
        "pair_include_trait_ids": sorted(include_ids),
        "pair_include_trait_prefixes": include_prefixes,
    }
    (meta_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"prepared {len(manifest_rows)} external-rg shards for "
        f"{n_pan} Pan-UKBB traits, {n_external} external traits, and {total_pairs} pairs",
        flush=True,
    )
    return manifest_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pan-manifest", default=Path("data/catalog/eur_gwas_manifest.tsv"), type=Path)
    parser.add_argument("--pan-sumstats-dir", default=Path("data/sumstats/eur"), type=Path)
    parser.add_argument(
        "--external-manifest",
        default=Path("data/external_gwas/catalog/external_gwas_manifest.tsv"),
        type=Path,
    )
    parser.add_argument("--combined-sumstats-dir", default=Path("data/external_rg/sumstats"), type=Path)
    parser.add_argument("--ld-prefix", default=Path("data/ld/UKBB.EUR"), type=Path)
    parser.add_argument(
        "--ldsc-bin",
        default=Path("external/ldsc-rs-rg-batch-target/release/ldsc"),
        type=Path,
    )
    parser.add_argument("--out-dir", default=Path("results/external_rg"), type=Path)
    parser.add_argument("--trait-block-size", default=64, type=int)
    parser.add_argument("--rayon-threads", default=9, type=int)
    parser.add_argument("--max-parallel-shards", default=1, type=int)
    parser.add_argument("--m-snps", default="auto")
    parser.add_argument("--allow-missing-sumstats", action="store_true")
    parser.add_argument("--force-shards", action="store_true")
    parser.add_argument("--force-symlinks", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--progress", action="store_true")
    parser.add_argument("--collect", action="store_true")
    parser.add_argument("--collect-out", type=Path)
    parser.add_argument("--allow-incomplete-collect", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--max-shards", type=int)
    parser.add_argument(
        "--pair-include-trait-id",
        action="append",
        default=[],
        help=(
            "Only emit pairs where at least one trait ID matches one of these comma-separated IDs. "
            "Can be passed multiple times."
        ),
    )
    parser.add_argument(
        "--pair-include-trait-prefix",
        action="append",
        default=[],
        help=(
            "Only emit pairs where at least one trait ID starts with one of these comma-separated prefixes. "
            "Can be passed multiple times."
        ),
    )
    parser.add_argument("--no-compress-output", dest="compress_output", action="store_false")
    parser.add_argument("--no-check-alleles", action="store_true")
    parser.set_defaults(compress_output=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.trait_block_size < 2:
        raise SystemExit("--trait-block-size must be >= 2")
    if args.max_parallel_shards < 1:
        raise SystemExit("--max-parallel-shards must be >= 1")
    if args.rayon_threads < 1:
        raise SystemExit("--rayon-threads must be >= 1")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.combined_sumstats_dir.mkdir(parents=True, exist_ok=True)
    # run_all_rg_hybrid.run_shards expects the Pan-UKBB runner argument name.
    args.sumstats_dir = args.combined_sumstats_dir

    if args.progress or args.collect:
        rows = rg.read_shard_manifest(args.out_dir)
    else:
        rows = prepare_external_shards(args)

    if args.prepare_only:
        return
    if args.progress:
        rg.print_progress(args, rows)
        return
    if args.collect:
        rg.collect_results(args, rows)
        return
    if not args.ldsc_bin.exists():
        raise SystemExit(f"Patched ldsc-rs binary not found: {args.ldsc_bin}")
    rg.run_shards(args, rows)
    rg.print_progress(args, rg.read_shard_manifest(args.out_dir))


if __name__ == "__main__":
    main()
