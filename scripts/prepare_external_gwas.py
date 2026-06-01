#!/usr/bin/env python3
"""Download external GWAS summary statistics and align them to Pan-UKBB SNP IDs."""

from __future__ import annotations

import argparse
import base64
import bz2
import csv
import gzip
import html
import io
import json
import math
import os
import re
import shutil
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
import zlib
from contextlib import contextmanager
from dataclasses import dataclass
from html.parser import HTMLParser
from pathlib import Path
from statistics import NormalDist
from typing import Iterator


GWAS_FTP_ROOT = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
PGC_DOWNLOADS_URL = "https://pgc.unc.edu/for-researchers/download-results/"
CNCR_SUMMARY_STATS_URL = "https://cncr.nl/research/summary_statistics/"
EGG_SUMMARY_STATS_URL = "https://egg-consortium.org/"
ZENODO_RECORD_URL = "https://zenodo.org/records/10515792"
ZENODO_SUMSTATS_URL = f"{ZENODO_RECORD_URL}/files/sumstats_indep107.tgz?download=1"
ZENODO_1000G_URL = f"{ZENODO_RECORD_URL}/files/1000G_Phase3_plinkfiles.tgz?download=1"
ZENODO_FACE_CGWAS_RECORD_URL = "https://zenodo.org/records/13730680"
USER_AGENT = "pan-ukbb-rg external-gwas downloader"
NA = {"", "NA", "NaN", "nan", "None", "none", "."}
COMPLEMENTS = str.maketrans("ACGT", "TGCA")
COMPRESSED_READ_ERRORS = (EOFError, gzip.BadGzipFile, tarfile.ReadError, zipfile.BadZipFile, zlib.error)
REQUESTED_PGC_PUBLICATIONS = {
    "adhd2022",
    "an2019",
    "anx2016",
    "anx2026",
    "anx2026_GADsymptsQuant",
    "asd2019",
    "bip2024",
    "bpd2025",
    "cdg2025",
    "ciac",
    "hoarding2022",
    "mdd2025",
    "ocs2024",
    "ocd2025",
    "panic2019",
    "ptsd2024",
    "scz2022",
    "sud2018-alc",
    "sud2019-alcuse",
    "sud2020-cud",
    "sud2020-op",
    "SUD2023",
    "sui2023",
    "ts2019",
    "alz2021",
}
CDG2025_FILES = {
    "PFactor_2025.tsv.gz": ("pfactor_grotzinger_2025", "Cross-disorder p-factor"),
    "F1_CompulsiveDisorders_2025.tsv.gz": (
        "compulsive_disorders_grotzinger_2025",
        "Cross-disorder compulsive disorders factor",
    ),
    "F2_SchizophreniaBipolar_2025.tsv.gz": (
        "schizophrenia_bipolar_grotzinger_2025",
        "Cross-disorder schizophrenia-bipolar factor",
    ),
    "F3_Neurodevelopmental_2025.tsv.gz": (
        "neurodevelopmental_grotzinger_2025",
        "Cross-disorder neurodevelopmental factor",
    ),
    "F4_Internalizing_2025.tsv.gz": (
        "internalizing_grotzinger_2025",
        "Cross-disorder internalizing factor",
    ),
    "F5_SubstanceUse_2025.tsv.gz": (
        "substance_use_grotzinger_2025",
        "Cross-disorder substance use factor",
    ),
}


@dataclass(frozen=True)
class PanSnp:
    snp: str
    chrom: str
    pos: str
    ref: str
    alt: str


@dataclass
class PanPanel:
    by_snp: dict[str, PanSnp]
    by_coord: dict[tuple[str, str, str, str], PanSnp]


class TableParser(HTMLParser):
    """Small HTML table parser for the PGC TablePress download table."""

    def __init__(self) -> None:
        super().__init__()
        self.in_table = False
        self.in_row = False
        self.in_cell = False
        self.current_cell: list[str] = []
        self.current_links: list[str] = []
        self.current_row: list[dict[str, object]] = []
        self.rows: list[list[dict[str, object]]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attr = {k: v or "" for k, v in attrs}
        if tag == "table" and attr.get("id") == "tablepress-5":
            self.in_table = True
            return
        if not self.in_table:
            return
        if tag == "tr":
            self.in_row = True
            self.current_row = []
        elif self.in_row and tag in {"td", "th"}:
            self.in_cell = True
            self.current_cell = []
            self.current_links = []
        elif self.in_cell and tag == "a" and attr.get("href"):
            self.current_links.append(attr["href"])

    def handle_endtag(self, tag: str) -> None:
        if tag == "table" and self.in_table:
            self.in_table = False
            return
        if not self.in_table:
            return
        if tag in {"td", "th"} and self.in_cell:
            text = html.unescape(" ".join("".join(self.current_cell).split()))
            self.current_row.append({"text": text, "links": list(self.current_links)})
            self.in_cell = False
        elif tag == "tr" and self.in_row:
            if self.current_row:
                self.rows.append(self.current_row)
            self.in_row = False

    def handle_data(self, data: str) -> None:
        if self.in_cell:
            self.current_cell.append(data)


def request_url(url: str, headers: dict[str, str] | None = None) -> urllib.request.Request:
    request_headers = {"User-Agent": USER_AGENT}
    if headers:
        request_headers.update(headers)
    return urllib.request.Request(url, headers=request_headers)


def safe_id(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return cleaned.strip("_") or "unknown"


def norm_chrom(value: str) -> str:
    value = value.strip()
    value = re.sub(r"^chr", "", value, flags=re.IGNORECASE)
    return value.upper()


def norm_allele(value: str) -> str:
    return value.strip().upper()


def is_missing(value: str | None) -> bool:
    return value is None or value.strip() in NA


def as_float(value: str | None) -> float | None:
    if is_missing(value):
        return None
    try:
        out = float(str(value).strip())
    except ValueError:
        return None
    if not math.isfinite(out):
        return None
    return out


def split_list(value: str | None) -> set[str]:
    if not value:
        return set()
    return {item.strip() for item in value.split(",") if item.strip()}


def parse_notes_kv(notes: str | None) -> dict[str, str]:
    out: dict[str, str] = {}
    if not notes:
        return out
    for part in notes.split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        out[key.strip().lower()] = value.strip()
    return out


def load_targets(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def apply_url_overrides(rows: list[dict[str, str]], path: Path | None) -> None:
    if path is None or not path.exists():
        return
    with path.open(newline="") as f:
        overrides = {
            row["external_id"]: row
            for row in csv.DictReader(f, delimiter="\t")
            if row.get("external_id")
        }
    for row in rows:
        override = overrides.get(row.get("external_id", ""))
        if not override:
            continue
        source_url = (override.get("source_url") or "").strip()
        if source_url:
            row["source_url"] = source_url
        source_file = (override.get("source_file") or "").strip()
        if source_file:
            row["source_file"] = source_file


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    with tmp.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def parse_pan_snp(value: str) -> PanSnp | None:
    parts = value.strip().split(":")
    if len(parts) != 4:
        return None
    chrom, pos, ref, alt = parts
    return PanSnp(
        snp=value.strip(),
        chrom=norm_chrom(chrom),
        pos=pos.strip(),
        ref=norm_allele(ref),
        alt=norm_allele(alt),
    )


def load_pan_panel(path: Path) -> PanPanel:
    by_snp: dict[str, PanSnp] = {}
    by_coord: dict[tuple[str, str, str, str], PanSnp] = {}
    with path.open() as f:
        for line in f:
            snp = parse_pan_snp(line.strip())
            if snp is None:
                continue
            by_snp[snp.snp] = snp
            by_coord[(snp.chrom, snp.pos, snp.ref, snp.alt)] = snp
    if not by_snp:
        raise ValueError(f"No Pan-UKBB SNP IDs loaded from {path}")
    return PanPanel(by_snp=by_snp, by_coord=by_coord)


def download_file(url: str, out: Path, force: bool = False) -> Path:
    if out.exists() and out.stat().st_size > 0 and not force:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(f"{out.name}.tmp.{os.getpid()}")
    last_error: Exception | None = None
    for attempt in range(1, 6):
        try:
            with urllib.request.urlopen(request_url(url), timeout=120) as resp, tmp.open("wb") as f:
                shutil.copyfileobj(resp, f, length=1024 * 1024)
            os.replace(tmp, out)
            return out
        except Exception as exc:  # noqa: BLE001 - record URL/download errors.
            last_error = exc
            if tmp.exists():
                tmp.unlink()
            time.sleep(min(30, attempt * 5))
    raise RuntimeError(f"Failed to download {url}: {last_error}")


@contextmanager
def open_text_path(path: Path, tolerate_compressed_eof: bool = False) -> Iterator[io.TextIOBase]:
    raw = path.open("rb")
    try:
        if path.name.endswith(".gz") or path.name.endswith(".bgz"):
            stream: io.BufferedIOBase = gzip.GzipFile(fileobj=raw)
        elif path.name.endswith(".bz2"):
            stream = bz2.BZ2File(raw)
        else:
            stream = raw
        text = io.TextIOWrapper(stream)
        try:
            yield text
        finally:
            try:
                text.close()
            except COMPRESSED_READ_ERRORS:
                if not tolerate_compressed_eof:
                    raise
    finally:
        if not raw.closed:
            raw.close()


@contextmanager
def open_text_url(url: str, name: str | None = None) -> Iterator[io.TextIOBase]:
    resp = urllib.request.urlopen(request_url(url), timeout=120)
    try:
        stream_name = name or url.split("?", 1)[0]
        if stream_name.endswith((".gz", ".bgz")):
            stream: io.BufferedIOBase = gzip.GzipFile(fileobj=resp)
        else:
            stream = resp
        text = io.TextIOWrapper(stream)
        try:
            yield text
        finally:
            text.close()
    finally:
        resp.close()


@contextmanager
def open_text_request(req: urllib.request.Request, name: str | None = None) -> Iterator[io.TextIOBase]:
    resp = urllib.request.urlopen(req, timeout=120)
    try:
        if name and name.endswith((".gz", ".bgz")):
            stream: io.BufferedIOBase = gzip.GzipFile(fileobj=resp)
        else:
            stream = resp
        text = io.TextIOWrapper(stream)
        try:
            yield text
        finally:
            text.close()
    finally:
        resp.close()


@contextmanager
def open_tar_member_text(archive: Path, member_name: str) -> Iterator[io.TextIOBase]:
    tar = tarfile.open(archive, "r:*")
    try:
        try:
            member = tar.getmember(member_name)
        except KeyError:
            names = {m.name for m in tar.getmembers()}
            if member_name not in names:
                raise FileNotFoundError(f"{member_name} is not present in {archive}")
            member = tar.getmember(member_name)
        extracted = tar.extractfile(member)
        if extracted is None:
            raise FileNotFoundError(f"{member_name} is not a regular file in {archive}")
        try:
            stream: io.BufferedIOBase
            if member_name.endswith(".gz"):
                stream = gzip.GzipFile(fileobj=extracted)
            else:
                stream = extracted
            text = io.TextIOWrapper(stream)
            try:
                yield text
            finally:
                text.close()
        finally:
            extracted.close()
    finally:
        tar.close()


def is_tar_archive_name(name: str) -> bool:
    lower = name.lower()
    return lower.endswith((".tar", ".tar.gz", ".tgz", ".tar.bz2", ".tbz2"))


def is_tar_archive_path(path: Path) -> bool:
    if is_tar_archive_name(path.name):
        return True
    try:
        return tarfile.is_tarfile(path)
    except OSError:
        return False


def is_zip_archive_name(name: str) -> bool:
    return name.lower().endswith(".zip")


def is_zip_archive_path(path: Path) -> bool:
    if is_zip_archive_name(path.name):
        return True
    try:
        return zipfile.is_zipfile(path)
    except OSError:
        return False


@contextmanager
def open_tar_data_text(archive: Path, row: dict[str, str]) -> Iterator[tuple[io.TextIOBase, str]]:
    notes = parse_notes_kv(row.get("notes", ""))
    requested = notes.get("tar_member")
    tar = tarfile.open(archive, "r:*")
    try:
        members = [member for member in tar.getmembers() if member.isfile()]
        chosen: tarfile.TarInfo | None = None
        if requested:
            for member in members:
                if member.name == requested or Path(member.name).name == requested:
                    chosen = member
                    break
            if chosen is None:
                raise FileNotFoundError(f"{requested} is not present in {archive}")
        else:
            candidates = [
                member
                for member in members
                if is_data_file_name(member.name) and not is_non_eur_name(member.name)
            ]
            if not candidates:
                raise FileNotFoundError(f"No data-like member found in {archive}")
            candidates.sort(
                key=lambda member: (
                    "readme" in member.name.lower(),
                    not is_eur_name(member.name),
                    -int(member.size or 0),
                    member.name,
                )
            )
            chosen = candidates[0]
        extracted = tar.extractfile(chosen)
        if extracted is None:
            raise FileNotFoundError(f"{chosen.name} is not a regular file in {archive}")
        try:
            stream: io.BufferedIOBase
            if chosen.name.endswith((".gz", ".bgz")):
                stream = gzip.GzipFile(fileobj=extracted)
            else:
                stream = extracted
            text = io.TextIOWrapper(stream)
            try:
                yield text, chosen.name
            finally:
                text.close()
        finally:
            extracted.close()
    finally:
        tar.close()


@contextmanager
def open_zip_data_text(archive: Path, row: dict[str, str]) -> Iterator[tuple[io.TextIOBase, str]]:
    notes = parse_notes_kv(row.get("notes", ""))
    requested = notes.get("zip_member") or notes.get("archive_member")
    zf = zipfile.ZipFile(archive)
    try:
        members = [member for member in zf.infolist() if not member.is_dir()]
        chosen: zipfile.ZipInfo | None = None
        if requested:
            for member in members:
                if member.filename == requested or Path(member.filename).name == requested:
                    chosen = member
                    break
            if chosen is None:
                raise FileNotFoundError(f"{requested} is not present in {archive}")
        else:
            candidates = [
                member
                for member in members
                if is_data_file_name(member.filename) and not is_non_eur_name(member.filename)
            ]
            if not candidates:
                raise FileNotFoundError(f"No data-like member found in {archive}")
            candidates.sort(
                key=lambda member: (
                    "readme" in member.filename.lower(),
                    not is_eur_name(member.filename),
                    -int(member.file_size or 0),
                    member.filename,
                )
            )
            chosen = candidates[0]
        extracted = zf.open(chosen)
        try:
            stream: io.BufferedIOBase
            if chosen.filename.endswith((".gz", ".bgz")):
                stream = gzip.GzipFile(fileobj=extracted)
            else:
                stream = extracted
            text = io.TextIOWrapper(stream)
            try:
                yield text, chosen.filename
            finally:
                text.close()
        finally:
            extracted.close()
    finally:
        zf.close()


def iter_text_lines(text: io.TextIOBase, tolerate_compressed_eof: bool = False) -> Iterator[str]:
    while True:
        try:
            yield next(text)
        except StopIteration:
            return
        except COMPRESSED_READ_ERRORS:
            if tolerate_compressed_eof:
                return
            raise


def read_header_and_rows(
    text: io.TextIOBase, tolerate_compressed_eof: bool = False
) -> tuple[list[str], Iterator[dict[str, str]]]:
    lines = iter_text_lines(text, tolerate_compressed_eof)
    header = ""
    for line in lines:
        if line.startswith("##"):
            continue
        if line.strip():
            header = line
            break
    if not header:
        raise ValueError("Input summary-statistics file is empty")
    header = header.rstrip("\n\r")
    if header.startswith("#"):
        header = header.lstrip("#")
    if "\t" in header:
        delimiter: str | None = "\t"
    elif "," in header:
        delimiter = ","
    else:
        delimiter = None
    if delimiter is None:
        fieldnames = header.split()

        def iter_split() -> Iterator[dict[str, str]]:
            for line in lines:
                parts = line.strip().split()
                if len(parts) < len(fieldnames):
                    continue
                yield dict(zip(fieldnames, parts))

        return fieldnames, iter_split()

    fieldnames = next(csv.reader([header], delimiter=delimiter))

    def iter_delimited() -> Iterator[dict[str, str]]:
        for line in lines:
            if line.startswith("##") or not line.strip():
                continue
            parts = next(csv.reader([line.rstrip("\n\r")], delimiter=delimiter))
            if len(parts) < len(fieldnames):
                whitespace_parts = line.strip().split()
                if len(whitespace_parts) >= len(fieldnames):
                    parts = whitespace_parts
                else:
                    continue
            yield {
                name: parts[index] if index < len(parts) else ""
                for index, name in enumerate(fieldnames)
            }

    return fieldnames, iter_delimited()


def normalize_col(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def column_map(fieldnames: list[str]) -> dict[str, str]:
    return {normalize_col(name): name for name in fieldnames}


def pick_col(cols: dict[str, str], aliases: list[str]) -> str | None:
    for alias in aliases:
        name = cols.get(normalize_col(alias))
        if name is not None:
            return name
    return None


def resolve_requested_col(fieldnames: list[str], requested: str | None) -> str | None:
    if not requested:
        return None
    if requested in fieldnames:
        return requested
    cols = column_map(fieldnames)
    return cols.get(normalize_col(requested))


def detect_columns(fieldnames: list[str], notes: str | None = None) -> dict[str, str | None]:
    cols = column_map(fieldnames)
    detected = {
        "chrom": pick_col(cols, ["chromosome", "chrom", "chr", "chr_hg19", "hm_chrom", "hm_chromosome"]),
        "pos": pick_col(
            cols,
            [
                "base_pair_location",
                "base_pair_position",
                "pos37",
                "position37",
                "bp37",
                "pos_hg19",
                "pos_grch37",
                "pos_hg19",
                "position",
                "pos",
                "bp",
                "hm_pos",
                "hm_coordinate",
            ],
        ),
        "effect": pick_col(
            cols,
            ["hm_effect_allele", "effect_allele", "ea", "eff", "a1", "a_1", "allele1", "tested_allele", "alt"],
        ),
        "other": pick_col(
            cols,
            [
                "hm_other_allele",
                "other_allele",
                "nea",
                "oa",
                "a2",
                "a_0",
                "allele2",
                "non_effect_allele",
                "ref",
            ],
        ),
        "beta": pick_col(cols, ["hm_beta", "beta", "eur_beta", "effect", "est", "b", "beta_0", "beta_t"]),
        "or": pick_col(cols, ["odds_ratio", "or"]),
        "se": pick_col(cols, ["standard_error", "stderr", "se", "eur_se", "sebeta", "se_0", "se_t"]),
        "z": pick_col(cols, ["hm_z", "z", "zscore", "z_score"]),
        "p": pick_col(cols, ["p_value", "pvalue", "pval", "p", "eur_p", "pval", "pvalue_association", "p_t"]),
        "rsid": pick_col(
            cols, ["hm_rsid", "rsid", "rs_id", "snp", "snpid", "markername", "marker_name", "id"]
        ),
        "n": pick_col(
            cols,
            [
                "n",
                "n_total",
                "ntotal",
                "totaln",
                "nsample",
                "n_sample",
                "total_sample_size",
                "totalsamplesize",
                "sample_size",
                "samplesize",
                "eur_n",
                "effective_n",
                "neff",
                "neffdiv2",
                "neff_half",
                "neffhalf",
                "weight",
            ],
        ),
        "n_cases": pick_col(cols, ["n_cases", "ncases", "nca", "n_case", "nca1"]),
        "n_controls": pick_col(cols, ["n_controls", "ncontrols", "nco", "n_control", "nco1"]),
    }
    overrides = parse_notes_kv(notes)
    for key in [
        "chrom",
        "pos",
        "effect",
        "other",
        "beta",
        "or",
        "se",
        "z",
        "p",
        "rsid",
        "n",
        "n_cases",
        "n_controls",
    ]:
        override = resolve_requested_col(fieldnames, overrides.get(f"{key}_col"))
        if override is not None:
            detected[key] = override
    return detected


def row_value(row: dict[str, str], column: str | None) -> str | None:
    if column is None:
        return None
    return row.get(column)


def signed_z(row: dict[str, str], cols: dict[str, str | None]) -> float | None:
    z = as_float(row_value(row, cols["z"]))
    if z is not None:
        return z
    beta = as_float(row_value(row, cols["beta"]))
    se = as_float(row_value(row, cols["se"]))
    if beta is not None and se is not None and se > 0:
        return beta / se
    odds_ratio = as_float(row_value(row, cols["or"]))
    if odds_ratio is not None and odds_ratio > 0 and se is not None and se > 0:
        return math.log(odds_ratio) / se
    p = as_float(row_value(row, cols["p"]))
    if p is not None and 0 < p <= 1 and beta is not None:
        return signed_z_from_beta_p(beta, p)
    return None


def signed_z_from_beta_p(beta: float, p: float) -> float | None:
    if not math.isfinite(beta) or not math.isfinite(p) or p <= 0 or p > 1:
        return None
    # Lower-tail inversion stays finite for very small p-values where 1 - p / 2
    # rounds to 1 in double precision.
    tail = min(max(p / 2, 1e-323), 0.5)
    z_abs = -NormalDist().inv_cdf(tail)
    if not math.isfinite(z_abs):
        return None
    return math.copysign(z_abs, beta)


def row_n(row: dict[str, str], cols: dict[str, str | None], fallback: str | None) -> float | None:
    n_col = cols["n"]
    value = as_float(row_value(row, n_col))
    if value is not None and value > 0:
        if n_col and normalize_col(n_col) in {"neffdiv2", "neffhalf"}:
            return 2 * value
        return value
    n_cases = as_float(row_value(row, cols.get("n_cases")))
    n_controls = as_float(row_value(row, cols.get("n_controls")))
    if n_cases is not None and n_controls is not None and n_cases > 0 and n_controls > 0:
        return n_cases + n_controls
    value = as_float(fallback)
    if value is not None and value > 0:
        return value
    return None


def orient_to_pan_alt(
    panel: PanPanel,
    chrom: str,
    pos: str,
    effect: str,
    other: str,
    z: float,
) -> tuple[PanSnp, float] | None:
    chrom = norm_chrom(chrom)
    pos = pos.strip()
    effect = norm_allele(effect)
    other = norm_allele(other)
    direct = panel.by_coord.get((chrom, pos, other, effect))
    reverse = panel.by_coord.get((chrom, pos, effect, other))
    if direct and reverse:
        return None
    if direct:
        return direct, z
    if reverse:
        return reverse, -z
    return None


def orient_rsid_to_pan_alt(
    rsid_map: dict[str, PanSnp],
    rsid: str,
    effect: str,
    other: str,
    z: float,
) -> tuple[PanSnp, float] | None:
    pan = rsid_map.get(rsid)
    if pan is None:
        return None
    effect = norm_allele(effect)
    other = norm_allele(other)
    if effect == pan.alt and other == pan.ref:
        return pan, z
    if effect == pan.ref and other == pan.alt:
        return pan, -z
    return None


def reverse_complement(allele: str) -> str:
    allele = norm_allele(allele)
    if not allele or any(base not in "ACGT" for base in allele):
        return allele
    return allele.translate(COMPLEMENTS)[::-1]


def convert_external_stream(
    text: io.TextIOBase,
    row: dict[str, str],
    panel: PanPanel,
    rsid_map: dict[str, PanSnp] | None,
    out: Path,
    stats_out: Path,
    max_rows: int | None,
    tolerate_compressed_eof: bool = False,
) -> dict[str, object]:
    started = time.time()
    fieldnames, rows = read_header_and_rows(text, tolerate_compressed_eof)
    cols = detect_columns(fieldnames, row.get("notes", ""))
    has_coord = all(cols[key] for key in ["chrom", "pos", "effect", "other"])
    has_rsid = rsid_map is not None and all(cols[key] for key in ["rsid", "effect", "other"])
    if not has_coord and not has_rsid:
        raise ValueError(
            "Could not find coordinate or rsID columns with effect/other alleles; "
            f"header={fieldnames}"
        )

    out.parent.mkdir(parents=True, exist_ok=True)
    stats_out.parent.mkdir(parents=True, exist_ok=True)
    out_tmp = out.with_name(f"{out.name}.tmp.{os.getpid()}")
    stats_tmp = stats_out.with_name(f"{stats_out.name}.tmp.{os.getpid()}")

    n_in = 0
    n_written = 0
    n_duplicate = 0
    n_missing = 0
    n_no_z = 0
    n_no_n = 0
    n_no_pan_match = 0
    seen: set[str] = set()

    try:
        with gzip.open(out_tmp, "wt", compresslevel=6) as dst:
            dst.write("SNP\tA1\tA2\tZ\tN\n")
            for data in rows:
                n_in += 1
                if max_rows is not None and n_in > max_rows:
                    break
                z = signed_z(data, cols)
                if z is None:
                    n_no_z += 1
                    continue
                n = row_n(data, cols, row.get("sample_size"))
                if n is None:
                    n_no_n += 1
                    continue
                effect = row_value(data, cols["effect"])
                other = row_value(data, cols["other"])
                if is_missing(effect) or is_missing(other):
                    n_missing += 1
                    continue

                oriented: tuple[PanSnp, float] | None = None
                if has_coord:
                    chrom = row_value(data, cols["chrom"])
                    pos = row_value(data, cols["pos"])
                    if not is_missing(chrom) and not is_missing(pos):
                        oriented = orient_to_pan_alt(panel, chrom or "", pos or "", effect or "", other or "", z)
                if oriented is None and has_rsid and rsid_map is not None:
                    rsid = row_value(data, cols["rsid"])
                    if not is_missing(rsid):
                        oriented = orient_rsid_to_pan_alt(rsid_map, rsid or "", effect or "", other or "", z)
                        if oriented is None:
                            rc_effect = reverse_complement(effect or "")
                            rc_other = reverse_complement(other or "")
                            if (rc_effect, rc_other) != (effect, other):
                                oriented = orient_rsid_to_pan_alt(
                                    rsid_map, rsid or "", rc_effect, rc_other, z
                                )
                if oriented is None:
                    n_no_pan_match += 1
                    continue

                pan, z_out = oriented
                if pan.snp in seen:
                    n_duplicate += 1
                    continue
                seen.add(pan.snp)
                dst.write(f"{pan.snp}\t{pan.alt}\t{pan.ref}\t{z_out:.8g}\t{n:.0f}\n")
                n_written += 1

        stats: dict[str, object] = {
            "external_id": row["external_id"],
            "source_type": row["source_type"],
            "trait": row.get("trait", ""),
            "status": "done",
            "population": row.get("population", ""),
            "source_url": row.get("source_url", ""),
            "source_file": row.get("source_file", ""),
            "genome_build": parse_notes_kv(row.get("notes", "")).get("genome_build", "GRCh37"),
            "input_rows": n_in,
            "written_rows": n_written,
            "total_panukbb_snps": len(panel.by_snp),
            "panukbb_overlap_pct": (100 * n_written / len(panel.by_snp)) if panel.by_snp else None,
            "input_overlap_pct": (100 * n_written / n_in) if n_in else None,
            "duplicate_snps": n_duplicate,
            "missing_alleles": n_missing,
            "missing_z": n_no_z,
            "missing_n": n_no_n,
            "no_panukbb_match": n_no_pan_match,
            "elapsed_seconds": time.time() - started,
            "out": str(out),
            "header": fieldnames,
        }
        if tolerate_compressed_eof:
            stats["tolerated_truncated_compressed_input"] = True
        with stats_tmp.open("w") as f:
            json.dump(stats, f, indent=2, sort_keys=True)
            f.write("\n")
        os.replace(out_tmp, out)
        os.replace(stats_tmp, stats_out)
        return stats
    finally:
        for tmp in [out_tmp, stats_tmp]:
            if tmp.exists():
                tmp.unlink()


def write_skip_stats(row: dict[str, str], stats_out: Path, status: str, message: str) -> None:
    stats_out.parent.mkdir(parents=True, exist_ok=True)
    tmp = stats_out.with_name(f"{stats_out.name}.tmp.{os.getpid()}")
    stats = {
        "external_id": row["external_id"],
        "source_type": row["source_type"],
        "trait": row.get("trait", ""),
        "population": row.get("population", ""),
        "source_url": row.get("source_url", ""),
        "source_file": row.get("source_file", ""),
        "genome_build": parse_notes_kv(row.get("notes", "")).get("genome_build", "GRCh37"),
        "status": status,
        "message": message,
        "written_rows": 0,
    }
    with tmp.open("w") as f:
        json.dump(stats, f, indent=2, sort_keys=True)
        f.write("\n")
    os.replace(tmp, stats_out)


def write_metadata_only_stats(
    row: dict[str, str],
    stats_out: Path,
    raw_path: Path | None,
    message: str,
) -> dict[str, object]:
    stats_out.parent.mkdir(parents=True, exist_ok=True)
    tmp = stats_out.with_name(f"{stats_out.name}.tmp.{os.getpid()}")
    stats: dict[str, object] = {
        "external_id": row["external_id"],
        "source_type": row["source_type"],
        "trait": row.get("trait", ""),
        "status": "metadata_only",
        "population": row.get("population", ""),
        "source_url": row.get("source_url", ""),
        "source_file": row.get("source_file", ""),
        "genome_build": parse_notes_kv(row.get("notes", "")).get("genome_build", "GRCh37"),
        "message": message,
        "written_rows": 0,
    }
    if raw_path is not None:
        stats["raw_path"] = str(raw_path)
    with tmp.open("w") as f:
        json.dump(stats, f, indent=2, sort_keys=True)
        f.write("\n")
    os.replace(tmp, stats_out)
    return stats


def existing_skip_result(
    row: dict[str, str],
    out: Path,
    stats: Path,
    force: bool,
) -> dict[str, object] | None:
    if force or not (out.exists() and stats.exists()):
        return None
    try:
        with stats.open() as f:
            existing = json.load(f)
    except Exception:  # noqa: BLE001 - stale/broken stats should be regenerated.
        return None
    if (
        existing.get("status") == "done"
        and int(existing.get("written_rows") or 0) == 0
    ):
        return None
    return {"external_id": row["external_id"], "status": "skip", "out": str(out)}


def gwas_catalog_dir(accession: str) -> str:
    suffix = accession.removeprefix("GCST")
    number = int(suffix)
    start = ((number - 1) // 1000) * 1000 + 1
    end = start + 999
    width = max(6, len(suffix))
    return f"{GWAS_FTP_ROOT}/GCST{start:0{width}d}-GCST{end:0{width}d}/{accession}/"


def hrefs_from_directory(url: str) -> list[str]:
    with urllib.request.urlopen(request_url(url), timeout=120) as resp:
        body = resp.read().decode("utf-8", errors="replace")
    return [urllib.parse.urljoin(url, html.unescape(h)) for h in re.findall(r'href="([^"]+)"', body)]


def gwas_catalog_raw_url(accession: str) -> str:
    root = gwas_catalog_dir(accession)
    hrefs = hrefs_from_directory(root)
    candidates = []
    for href in hrefs:
        name = href.rsplit("/", 1)[-1]
        lower = name.lower()
        if "harmonised" in href.lower() or "metadata" in lower:
            continue
        if lower.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz", ".csv", ".csv.gz")):
            candidates.append(href)
    if not candidates:
        raise FileNotFoundError(f"No raw summary-statistics file found in {root}")
    candidates.sort(
        key=lambda url: (
            not url.rsplit("/", 1)[-1].startswith(accession),
            "buildgrch37" not in url.lower(),
            url,
        )
    )
    return candidates[0]


def fetch_pgc_rows(page_url: str) -> list[dict[str, str]]:
    with urllib.request.urlopen(request_url(page_url), timeout=120) as resp:
        body = resp.read().decode("utf-8", errors="replace")
    parser = TableParser()
    parser.feed(body)
    rows: list[dict[str, str]] = []
    for cells in parser.rows:
        values = [str(cell["text"]) for cell in cells]
        if len(values) != 7 or values[0] == "pgcGroup":
            continue
        links = cells[6]["links"] if cells[6]["links"] else []
        download_url = str(links[0]) if links else values[6]
        data_doi = values[5].strip()
        publication = values[2].strip()
        publication_id = safe_id(publication)
        if publication_id not in REQUESTED_PGC_PUBLICATIONS:
            continue
        if publication_id == "alz2021" and not data_doi.startswith("10.6084/m9.figshare."):
            # This public PGC Alzheimer row is served through CNCR/CTG rather than Figshare.
            # Keep the explicit Nextcloud target row to avoid a duplicate unsupported PGC row.
            continue
        external_id = safe_id(f"PGC.{publication}")
        if not publication_id:
            external_id = safe_id(f"PGC.{values[0]}.{values[1]}")
        base_row = {
            "external_id": external_id,
            "source_type": "pgc",
            "trait": values[1],
            "source_url": download_url,
            "population": "European",
            "category": values[0],
            "sample_size": "NA",
            "notes": (
                f"pgc_publication={publication}; journal={values[3]}; PMID={values[4]}; "
                f"data_doi={data_doi}; genome_build=GRCh37"
            ),
            "data_doi": data_doi,
        }
        if not data_doi.startswith("10.6084/m9.figshare."):
            rows.append(base_row)
            continue
        try:
            article = figshare_article_from_doi(data_doi)
            if publication_id == "cdg2025":
                rows.extend(expand_cdg2025_rows(base_row, article))
                continue
            file_info = choose_figshare_file(article)
            base_row["source_file"] = str(file_info.get("name") or "")
            base_row["notes"] += (
                f"; file_name={file_info.get('name')}; "
                f"download_url={file_info.get('download_url')}"
            )
        except Exception as exc:  # noqa: BLE001 - leave row processable/skippable later.
            base_row["notes"] += f"; figshare_discovery_error={safe_id(str(exc))}"
        rows.append(base_row)
    return rows


def expand_cdg2025_rows(base_row: dict[str, str], article: dict[str, object]) -> list[dict[str, str]]:
    files = article.get("files") or []
    by_name = {
        str(item.get("name", "")): item
        for item in files
        if isinstance(item, dict) and str(item.get("name", "")) in CDG2025_FILES
    }
    rows: list[dict[str, str]] = []
    for file_name, (external_id, trait) in CDG2025_FILES.items():
        item = by_name.get(file_name)
        row = dict(base_row)
        row["external_id"] = external_id
        row["trait"] = trait
        row["source_file"] = file_name
        if item:
            row["notes"] += (
                f"; cdg2025_factor_file={file_name}; file_name={file_name}; "
                f"download_url={item.get('download_url')}"
            )
        else:
            row["notes"] += f"; cdg2025_factor_file={file_name}; figshare_file_missing=true"
        rows.append(row)
    return rows


def is_non_eur_name(name: str) -> bool:
    lower = name.lower()
    tokens = re.split(r"[^a-z0-9]+", lower)
    token_set = set(tokens)
    non_eur_tokens = {
        "afr",
        "afram",
        "african",
        "eas",
        "asian",
        "asn",
        "sas",
        "latino",
        "latinx",
        "amr",
        "hispanic",
        "aam",
        "hna",
        "native",
        "multiancestry",
        "multi",
    }
    if token_set & non_eur_tokens:
        return True
    return any(
        flag in lower
        for flag in [
            "multi-ancestry",
            "all_ancestry",
            "all-ancestry",
            "transethnic",
            "trans-ethnic",
            "trans_ancestral",
            "trans-ancestral",
        ]
    )


def is_eur_name(name: str) -> bool:
    lower = name.lower()
    tokens = set(re.split(r"[^a-z0-9]+", lower))
    return bool(tokens & {"eur", "euro", "european"}) or "wo23andme" in lower


def figshare_article_from_doi(doi: str) -> dict[str, object]:
    url = "https://api.figshare.com/v2/articles?doi=" + urllib.parse.quote(doi, safe="")
    with urllib.request.urlopen(request_url(url), timeout=120) as resp:
        matches = json.load(resp)
    if not matches:
        raise FileNotFoundError(f"No Figshare article found for DOI {doi}")
    article_id = matches[0]["id"]
    with urllib.request.urlopen(
        request_url(f"https://api.figshare.com/v2/articles/{article_id}"), timeout=120
    ) as resp:
        return json.load(resp)


def choose_figshare_file(article: dict[str, object]) -> dict[str, object]:
    files = article.get("files") or []
    if not isinstance(files, list):
        raise ValueError("Figshare article has no file list")
    bad_ext = (".pdf", ".doc", ".docx", ".xls", ".xlsx", ".png", ".jpg", ".jpeg")
    good_ext = (".gz", ".bgz", ".tsv", ".txt", ".csv", ".sumstats")
    candidates = []
    for item in files:
        if not isinstance(item, dict):
            continue
        name = str(item.get("name", ""))
        lower = name.lower()
        if lower.endswith(bad_ext) or "readme" in lower:
            continue
        if any(flag in lower for flag in ["top10k", "pgs", "prs", "hits"]):
            continue
        if is_non_eur_name(name):
            continue
        if lower.endswith(good_ext):
            candidates.append(item)
    if not candidates:
        raise FileNotFoundError("No data-like Figshare file found")
    candidates.sort(
        key=lambda item: (
            not is_eur_name(str(item.get("name", ""))),
            str(item.get("name", "")).lower().startswith("sensitivity"),
            not any(
                word in str(item.get("name", "")).lower()
                for word in ["noukb", "noukbb", "no_ukb", "no_ukbb", "without_ukb", "without_ukbb"]
            ),
            not any(word in str(item.get("name", "")).lower() for word in ["no23andme", "wo23andme"]),
            bool(re.search(r"(^|[_\-.])(male|female|mal|fem)([_\-.]|$)", str(item.get("name", "")).lower())),
            not any(
                word in str(item.get("name", "")).lower()
                for word in ["looukbb", "meta", "european", "eur"]
            ),
            -int(item.get("size") or 0),
        )
    )
    return candidates[0]


def nextcloud_share_token(url: str) -> str:
    path = urllib.parse.urlparse(url).path.rstrip("/")
    token = path.rsplit("/", 1)[-1]
    if not token:
        raise ValueError(f"Could not infer Nextcloud share token from {url}")
    return token


def nextcloud_base(url: str) -> str:
    parsed = urllib.parse.urlparse(url)
    netloc = "vu.data.surf.nl" if parsed.netloc == "vu.data.surfsara.nl" else parsed.netloc
    return f"{parsed.scheme or 'https'}://{netloc}"


def nextcloud_auth_header(token: str) -> str:
    return "Basic " + base64.b64encode(f"{token}:".encode()).decode()


def nextcloud_webdav_url(row: dict[str, str]) -> str:
    notes = parse_notes_kv(row.get("notes", ""))
    token = notes.get("share_token") or nextcloud_share_token(row["source_url"])
    file_name = row.get("source_file") or notes.get("file")
    if not file_name:
        raise ValueError(f"Nextcloud row is missing file name: {row['external_id']}")
    return f"{nextcloud_base(row['source_url'])}/public.php/webdav/{urllib.parse.quote(file_name, safe='/')}"


def nextcloud_request(row: dict[str, str], method: str = "GET") -> urllib.request.Request:
    notes = parse_notes_kv(row.get("notes", ""))
    token = notes.get("share_token") or nextcloud_share_token(row["source_url"])
    return request_url(
        nextcloud_webdav_url(row),
        headers={"Authorization": nextcloud_auth_header(token)},
    ) if method == "GET" else urllib.request.Request(
        f"{nextcloud_base(row['source_url'])}/public.php/webdav/",
        method=method,
        data=b"",
        headers={
            "User-Agent": USER_AGENT,
            "Depth": "1",
            "Authorization": nextcloud_auth_header(token),
        },
    )


def list_nextcloud_share(share_url: str) -> list[dict[str, object]]:
    token = nextcloud_share_token(share_url)
    base = nextcloud_base(share_url)
    ns = {"d": "DAV:"}
    files: list[dict[str, object]] = []
    pending = [""]
    seen_dirs: set[str] = set()
    while pending:
        prefix = pending.pop()
        if prefix in seen_dirs:
            continue
        seen_dirs.add(prefix)
        url = f"{base}/public.php/webdav/{urllib.parse.quote(prefix, safe='/')}"
        req = urllib.request.Request(
            url,
            method="PROPFIND",
            data=b"",
            headers={
                "User-Agent": USER_AGENT,
                "Depth": "1",
                "Authorization": nextcloud_auth_header(token),
            },
        )
        with urllib.request.urlopen(req, timeout=120) as resp:
            tree = ET.fromstring(resp.read())
        for response in tree.findall("d:response", ns):
            href = response.findtext("d:href", default="", namespaces=ns)
            rel = urllib.parse.unquote(href)
            marker = "/public.php/webdav/"
            if marker in rel:
                rel = rel.split(marker, 1)[1]
            rel = rel.strip("/")
            if rel == prefix.strip("/"):
                continue
            if not rel:
                continue
            if response.find(".//d:collection", ns) is not None:
                pending.append(rel)
                continue
            size_text = response.findtext(".//d:getcontentlength", default="0", namespaces=ns)
            files.append({"name": rel, "size": int(size_text or 0)})
    return files


def is_data_file_name(name: str) -> bool:
    lower = name.lower()
    if lower.endswith((".pdf", ".doc", ".docx", ".xls", ".xlsx", ".png", ".jpg", ".jpeg")):
        return False
    if any(
        flag in lower
        for flag in [
            "readme",
            "license",
            "clump",
            "manhattan",
            "qq.",
            "header_explanation",
            "file_header",
            "explanation",
        ]
    ):
        return False
    return lower.endswith((".gz", ".bgz", ".tsv", ".txt", ".csv", ".sumstats", ".tbl"))


def fetch_cncr_rows(page_url: str) -> list[dict[str, str]]:
    with urllib.request.urlopen(request_url(page_url), timeout=120) as resp:
        body = resp.read().decode("utf-8", errors="replace")
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for match in re.finditer(r'href=["\']([^"\']+)["\']', body, re.IGNORECASE):
        href = html.unescape(match.group(1))
        if "vu.data.surf" not in href and "vu.data.surfsara" not in href:
            continue
        share_url = urllib.parse.urljoin(page_url, href)
        token = nextcloud_share_token(share_url)
        try:
            files = list_nextcloud_share(share_url)
        except Exception:
            continue
        prefix = re.sub("<[^>]+>", " ", body[max(0, match.start() - 1600) : match.start()])
        prefix = " ".join(html.unescape(prefix).split())
        year = re.search(r"\b(20[0-9]{2}|19[0-9]{2})\b", prefix)
        title = prefix[-220:] if prefix else "CNCR summary statistics"
        for file_info in files:
            name = str(file_info["name"])
            if not is_data_file_name(name) or is_non_eur_name(name):
                continue
            basename = Path(name).name
            if basename.lower() == "meta.txt.gz":
                continue
            key = (token, name)
            if key in seen:
                continue
            seen.add(key)
            sample_size = "24442" if "_FC_sumstats" in basename else "NA"
            external_id = safe_id(basename)
            for suffix in [".sumstats.txt.gz", ".txt.gz", ".tsv.gz", ".gz", ".txt", ".tsv"]:
                external_id = external_id.removesuffix(safe_id(suffix))
            rows.append(
                {
                    "external_id": f"CNCR.{external_id}",
                    "source_type": "nextcloud_file",
                    "trait": basename.removesuffix(".gz").removesuffix(".txt").removesuffix(".tsv"),
                    "source_url": share_url,
                    "source_file": name,
                    "population": "European",
                    "category": "CNCR",
                    "sample_size": sample_size,
                    "notes": (
                        f"share_token={token}; file={name}; genome_build=GRCh37; "
                        f"cncr_year={year.group(1) if year else 'NA'}; cncr_context={title}"
                    ),
                }
            )
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["external_id"]] = counts.get(row["external_id"], 0) + 1
    for row in rows:
        if counts[row["external_id"]] > 1:
            token = parse_notes_kv(row.get("notes", "")).get("share_token", "share")
            row["external_id"] = f"{row['external_id']}.{token[:8]}"
    return rows


def fetch_egg_rows(page_url: str) -> list[dict[str, str]]:
    with urllib.request.urlopen(request_url(page_url), timeout=120) as resp:
        body = resp.read().decode("utf-8", errors="replace")
    page_links: list[tuple[str, str]] = []
    for match in re.finditer(
        r'<a[^>]+href=["\']([^"\']+\.html)["\'][^>]*>(.*?)</a>',
        body,
        re.IGNORECASE | re.DOTALL,
    ):
        href = html.unescape(match.group(1))
        if href.startswith(("http://", "https://")) and "egg-consortium.org" not in href:
            continue
        label = " ".join(html.unescape(re.sub("<[^>]+>", " ", match.group(2))).split())
        page_links.append((label or Path(href).stem, urllib.parse.urljoin(page_url, href)))

    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for page_label, source_page in page_links:
        try:
            with urllib.request.urlopen(request_url(source_page), timeout=120) as resp:
                page_body = resp.read().decode("utf-8", errors="replace")
        except Exception:
            continue
        for match in re.finditer(
            r'<a[^>]+href=["\']([^"\']+)["\'][^>]*>(.*?)</a>',
            page_body,
            re.IGNORECASE | re.DOTALL,
        ):
            href = html.unescape(match.group(1))
            data_url = urllib.parse.urljoin(source_page, href)
            parsed = urllib.parse.urlparse(data_url)
            if "egg-consortium.org" not in parsed.netloc:
                continue
            file_name = Path(urllib.parse.unquote(parsed.path)).name
            if not file_name or not is_data_file_name(file_name):
                continue
            if is_non_eur_name(file_name):
                continue
            label = " ".join(html.unescape(re.sub("<[^>]+>", " ", match.group(2))).split())
            external_id = f"EGG.{safe_id(file_name)}"
            for suffix in [
                ".txt.gz",
                ".tsv.gz",
                ".csv.gz",
                ".sumstats.gz",
                ".gz",
                ".txt",
                ".tsv",
                ".csv",
            ]:
                external_id = external_id.removesuffix(safe_id(suffix))
            if external_id in seen:
                continue
            seen.add(external_id)
            trait = file_name.removesuffix(".gz").removesuffix(".txt").removesuffix(".tsv").removesuffix(".csv")
            rows.append(
                {
                    "external_id": external_id,
                    "source_type": "egg",
                    "trait": trait,
                    "source_url": data_url,
                    "source_file": file_name,
                    "population": "European",
                    "category": "EGG",
                    "sample_size": "NA",
                    "notes": (
                        f"genome_build=GRCh37; egg_page={source_page}; "
                        f"egg_page_label={page_label}; egg_link_label={label}"
                    ),
                }
            )
    return rows


def zenodo_api_record_url(record_url: str) -> str:
    if "/api/records/" in record_url:
        return record_url
    parsed = urllib.parse.urlparse(record_url)
    record_id = parsed.path.rstrip("/").rsplit("/", 1)[-1]
    if not record_id:
        raise ValueError(f"Could not infer Zenodo record id from {record_url}")
    return f"https://zenodo.org/api/records/{record_id}"


def zenodo_record_files(record_url: str) -> list[dict[str, object]]:
    with urllib.request.urlopen(request_url(zenodo_api_record_url(record_url)), timeout=120) as resp:
        record = json.load(resp)
    files = record.get("files") or []
    if not isinstance(files, list):
        raise ValueError(f"Zenodo record has no file list: {record_url}")
    return [item for item in files if isinstance(item, dict)]


def zenodo_face_record_id(record_url: str) -> str:
    return zenodo_api_record_url(record_url).rstrip("/").rsplit("/", 1)[-1]


def zenodo_face_raw_dir(out_dir: Path, record_url: str) -> Path:
    return out_dir / "raw" / "zenodo_face_cgwas" / zenodo_face_record_id(record_url)


def zenodo_file_download_url(file_info: dict[str, object]) -> str:
    links = file_info.get("links") or {}
    if isinstance(links, dict) and links.get("self"):
        return str(links["self"])
    key = str(file_info.get("key") or "")
    if not key:
        raise ValueError(f"Zenodo file is missing a key/download link: {file_info}")
    return f"{zenodo_api_record_url(ZENODO_FACE_CGWAS_RECORD_URL)}/files/{urllib.parse.quote(key)}/content"


def download_zenodo_face_file(
    file_info: dict[str, object],
    out_dir: Path,
    record_url: str,
    force: bool,
) -> Path:
    key = str(file_info.get("key") or "")
    if not key:
        raise ValueError(f"Zenodo file is missing key: {file_info}")
    return download_file(
        zenodo_file_download_url(file_info),
        zenodo_face_raw_dir(out_dir, record_url) / key,
        force,
    )


def strip_archive_suffix(name: str) -> str:
    for suffix in [".tar.bz2", ".tar.gz", ".tbz2", ".tgz", ".tar", ".bz2", ".gz"]:
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def fetch_zenodo_face_rows(
    record_url: str,
    out_dir: Path,
    force: bool = False,
    limit: int | None = None,
) -> list[dict[str, str]]:
    record_id = zenodo_face_record_id(record_url)
    index_path = out_dir / "catalog" / f"zenodo_{record_id}_face_members.tsv"
    fields = [
        "external_id",
        "source_type",
        "trait",
        "source_url",
        "population",
        "category",
        "sample_size",
        "notes",
    ]
    if index_path.exists() and not force:
        with index_path.open(newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    files = zenodo_record_files(record_url)
    by_key = {str(item.get("key") or ""): item for item in files}
    snp_info = by_key.get("SnpInfo.tsv.bz2")
    if snp_info is None:
        raise FileNotFoundError(f"SnpInfo.tsv.bz2 is missing from Zenodo record {record_url}")
    snp_info_url = zenodo_file_download_url(snp_info)
    rows: list[dict[str, str]] = []
    archive_files = [
        item for item in files if str(item.get("key") or "").endswith(".tar.bz2")
    ]
    archive_files.sort(key=lambda item: (int(item.get("size") or 0), str(item.get("key") or "")))
    for file_info in archive_files:
        key = str(file_info.get("key") or "")
        archive = download_zenodo_face_file(file_info, out_dir, record_url, force)
        archive_stem = strip_archive_suffix(key)
        with tarfile.open(archive, "r:*") as tar:
            members = [
                member
                for member in tar.getmembers()
                if member.isfile() and member.name.lower().endswith(".tsv")
            ]
        for member in sorted(members, key=lambda item: item.name):
            member_stem = Path(member.name).name.removesuffix(".tsv")
            external_id = safe_id(f"face_cgwas_2024_{archive_stem}_{member_stem}")
            trait = f"Facial shape distance {archive_stem} {member_stem}"
            rows.append(
                {
                    "external_id": external_id,
                    "source_type": "zenodo_face_cgwas",
                    "trait": trait,
                    "source_url": zenodo_file_download_url(file_info),
                    "population": "European",
                    "category": "Facial shape",
                    "sample_size": "10115",
                    "notes": (
                        f"genome_build=GRCh37; zenodo_record={record_id}; archive_name={key}; "
                        f"tar_member={member.name}; snp_info_file=SnpInfo.tsv.bz2; "
                        f"snp_info_url={snp_info_url}; beta_aligned_to=A1; "
                        "publication=Shaffer et al. 2019 / C-GWAS facial-shape source data; "
                        "phenotype_description=Distance between facial landmarks from the "
                        "C-GWAS human facial-shape Zenodo source-data archive."
                    ),
                }
            )
            if limit is not None and len(rows) >= limit:
                return rows
    write_tsv(index_path, rows, fields)
    return rows


def ensure_zenodo_archive(args: argparse.Namespace) -> Path:
    if args.zenodo_archive:
        path = Path(args.zenodo_archive)
        if not path.exists():
            raise FileNotFoundError(path)
        return path
    return download_file(
        ZENODO_SUMSTATS_URL,
        args.out_dir / "raw" / "zenodo_indep107" / "sumstats_indep107.tgz",
        args.force,
    )


def ensure_1000g_archive(args: argparse.Namespace) -> Path:
    return download_file(
        ZENODO_1000G_URL,
        args.out_dir / "raw" / "zenodo_indep107" / "1000G_Phase3_plinkfiles.tgz",
        args.force,
    )


def match_reference_variant(
    panel: PanPanel, chrom: str, pos: str, allele1: str, allele2: str
) -> PanSnp | None:
    chrom = norm_chrom(chrom)
    pos = pos.strip()
    a1 = norm_allele(allele1)
    a2 = norm_allele(allele2)
    return panel.by_coord.get((chrom, pos, a1, a2)) or panel.by_coord.get((chrom, pos, a2, a1))


def write_rsid_map_entry(
    out: io.TextIOBase,
    seen: set[str],
    panel: PanPanel,
    rsid: str,
    chrom: str,
    pos: str,
    allele1: str,
    allele2: str,
) -> int:
    if not rsid.startswith("rs") or rsid in seen:
        return 0
    pan = match_reference_variant(panel, chrom, pos, allele1, allele2)
    if pan is None:
        return 0
    seen.add(rsid)
    out.write(f"{rsid}\t{pan.snp}\t{pan.ref}\t{pan.alt}\n")
    return 1


def build_rsid_map_from_reference(reference: Path, panel: PanPanel, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(f"{out.name}.tmp.{os.getpid()}")
    seen: set[str] = set()
    matched = 0
    with open_text_path(reference) as inp, gzip.open(tmp, "wt", compresslevel=6) as dst:
        dst.write("rsid\tSNP\tREF\tALT\n")
        first = inp.readline()
        if not first:
            raise ValueError(f"Empty rsID reference file: {reference}")
        columns = first.strip().split()
        has_header = {"SNP", "CHR", "BP", "A1", "A2"}.issubset(set(columns))
        if has_header:
            idx = {name: i for i, name in enumerate(columns)}
        else:
            parts = columns
            if len(parts) >= 6:
                matched += write_rsid_map_entry(dst, seen, panel, parts[1], parts[0], parts[3], parts[4], parts[5])
            idx = {}
        for line in inp:
            parts = line.strip().split()
            if len(parts) < 6:
                continue
            if has_header:
                matched += write_rsid_map_entry(
                    dst,
                    seen,
                    panel,
                    parts[idx["SNP"]],
                    parts[idx["CHR"]],
                    parts[idx["BP"]],
                    parts[idx["A1"]],
                    parts[idx["A2"]],
                )
            else:
                matched += write_rsid_map_entry(dst, seen, panel, parts[1], parts[0], parts[3], parts[4], parts[5])
    if matched == 0:
        tmp.unlink(missing_ok=True)
        raise ValueError(f"No rsIDs in {reference} matched Pan-UKBB SNP IDs")
    os.replace(tmp, out)


def build_rsid_map_from_1000g_archive(archive: Path, panel: PanPanel, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(f"{out.name}.tmp.{os.getpid()}")
    seen: set[str] = set()
    matched = 0
    with tarfile.open(archive, "r:gz") as tar, gzip.open(tmp, "wt", compresslevel=6) as dst:
        dst.write("rsid\tSNP\tREF\tALT\n")
        for member in tar.getmembers():
            if not member.name.endswith(".bim") or not member.isfile():
                continue
            extracted = tar.extractfile(member)
            if extracted is None:
                continue
            with io.TextIOWrapper(extracted) as inp:
                for line in inp:
                    parts = line.strip().split()
                    if len(parts) < 6:
                        continue
                    matched += write_rsid_map_entry(
                        dst, seen, panel, parts[1], parts[0], parts[3], parts[4], parts[5]
                    )
    if matched == 0:
        tmp.unlink(missing_ok=True)
        raise ValueError(f"No rsIDs in {archive} matched Pan-UKBB SNP IDs")
    os.replace(tmp, out)


def load_rsid_map(path: Path) -> dict[str, PanSnp]:
    out: dict[str, PanSnp] = {}
    with gzip.open(path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pan = parse_pan_snp(row["SNP"])
            if pan is not None:
                out[row["rsid"]] = pan
    if not out:
        raise ValueError(f"No rsID mappings loaded from {path}")
    return out


def ensure_rsid_map(args: argparse.Namespace, panel: PanPanel) -> dict[str, PanSnp]:
    out = args.out_dir / "maps" / "panukbb_rsids.tsv.gz"
    if not out.exists() or args.force:
        if args.rsid_reference:
            build_rsid_map_from_reference(Path(args.rsid_reference), panel, out)
        else:
            build_rsid_map_from_1000g_archive(ensure_1000g_archive(args), panel, out)
    return load_rsid_map(out)


def selected_rows(rows: list[dict[str, str]], args: argparse.Namespace) -> list[dict[str, str]]:
    sources = split_list(args.sources)
    include = split_list(args.include)
    source_aliases = set(sources)
    if "cncr" in sources:
        source_aliases.add("nextcloud_file")
    if "pgc" in sources:
        source_aliases.add("figshare_file")
    selected = [
        row
        for row in rows
        if (not sources or row["source_type"] in source_aliases) and is_european_row(row)
    ]
    if include:
        selected = [row for row in selected if row["external_id"] in include]
    deduped: list[dict[str, str]] = []
    seen: set[str] = set()
    for row in selected:
        external_id = row["external_id"]
        if external_id in seen:
            continue
        seen.add(external_id)
        deduped.append(row)
    selected = deduped
    if args.limit:
        selected = selected[: args.limit]
    return selected


def is_european_row(row: dict[str, str]) -> bool:
    population = row.get("population", "").lower()
    if "east asian" in population or "african" in population or "south asian" in population:
        return False
    if "trans-ancestry" in population or "trans ancestry" in population:
        return True
    if "european" in population or population == "eur":
        return True
    return row["source_type"] in {"pgc", "figshare_file", "egg"}


MANIFEST_FIELDS = [
    "external_id",
    "source_type",
    "trait",
    "population",
    "category",
    "sample_size",
    "source_url",
    "source_file",
    "genome_build",
    "status",
    "aligned_path",
    "input_rows",
    "written_rows",
    "total_panukbb_snps",
    "panukbb_overlap_pct",
    "input_overlap_pct",
    "no_panukbb_match",
    "duplicate_snps",
    "missing_z",
    "missing_n",
    "message",
    "notes",
]


def load_stats(stats_dir: Path, external_id: str) -> dict[str, object]:
    path = stats_dir / f"{external_id}.json"
    if not path.exists():
        return {}
    with path.open() as f:
        return json.load(f)


def manifest_rows(rows: list[dict[str, str]], stats_dir: Path | None = None) -> list[dict[str, str]]:
    out_rows: list[dict[str, str]] = []
    for row in rows:
        stats = load_stats(stats_dir, row["external_id"]) if stats_dir else {}
        notes = parse_notes_kv(row.get("notes", ""))
        out_rows.append(
            {
                "external_id": row.get("external_id", ""),
                "source_type": row.get("source_type", ""),
                "trait": row.get("trait", ""),
                "population": row.get("population", ""),
                "category": row.get("category", ""),
                "sample_size": row.get("sample_size", ""),
                "source_url": row.get("source_url", ""),
                "source_file": str(stats.get("source_file") or row.get("source_file", "")),
                "genome_build": str(stats.get("genome_build") or notes.get("genome_build", "GRCh37")),
                "status": str(stats.get("status", "pending")),
                "aligned_path": str(stats.get("out", "")),
                "input_rows": str(stats.get("input_rows", "")),
                "written_rows": str(stats.get("written_rows", "")),
                "total_panukbb_snps": str(stats.get("total_panukbb_snps", "")),
                "panukbb_overlap_pct": format_float(stats.get("panukbb_overlap_pct")),
                "input_overlap_pct": format_float(stats.get("input_overlap_pct")),
                "no_panukbb_match": str(stats.get("no_panukbb_match", "")),
                "duplicate_snps": str(stats.get("duplicate_snps", "")),
                "missing_z": str(stats.get("missing_z", "")),
                "missing_n": str(stats.get("missing_n", "")),
                "message": str(stats.get("message", "")),
                "notes": row.get("notes", ""),
            }
        )
    return out_rows


def format_float(value: object) -> str:
    if value is None or value == "":
        return ""
    try:
        return f"{float(value):.6g}"
    except (TypeError, ValueError):
        return str(value)


def write_external_manifest(path: Path, rows: list[dict[str, str]], stats_dir: Path | None = None) -> None:
    write_tsv(path, manifest_rows(rows, stats_dir), MANIFEST_FIELDS)


def legacy_source_manifest_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    fields = [
        "external_id",
        "source_type",
        "trait",
        "source_url",
        "population",
        "category",
        "sample_size",
        "notes",
    ]
    return [{field: row.get(field, "") for field in fields} for row in rows]


def process_gwas_catalog(row: dict[str, str], args: argparse.Namespace, panel: PanPanel) -> dict[str, object]:
    out = args.out_dir / "aligned" / f"{row['external_id']}.sumstats.gz"
    stats = args.out_dir / "prepare_stats" / f"{row['external_id']}.json"
    if skip := existing_skip_result(row, out, stats, args.force):
        return skip
    raw_url = gwas_catalog_raw_url(row["external_id"])
    row["raw_url"] = raw_url
    if args.max_rows is not None:
        with open_text_url(raw_url) as text:
            return convert_external_stream(text, row, panel, None, out, stats, args.max_rows)
    raw_name = raw_url.rsplit("/", 1)[-1]
    raw = download_file(raw_url, args.out_dir / "raw" / "gwas_catalog" / row["external_id"] / raw_name, args.force)
    with open_text_path(raw) as text:
        return convert_external_stream(text, row, panel, None, out, stats, args.max_rows)


def process_zenodo(row: dict[str, str], args: argparse.Namespace, panel: PanPanel, rsid_map: dict[str, PanSnp]) -> dict[str, object]:
    out = args.out_dir / "aligned" / f"{row['external_id']}.sumstats.gz"
    stats = args.out_dir / "prepare_stats" / f"{row['external_id']}.json"
    if skip := existing_skip_result(row, out, stats, args.force):
        return skip
    archive = ensure_zenodo_archive(args)
    member = f"sumstats_107/{row['external_id']}.sumstats.gz"
    try:
        with open_tar_member_text(archive, member) as text:
            return convert_external_stream(text, row, panel, rsid_map, out, stats, args.max_rows)
    except FileNotFoundError as exc:
        write_skip_stats(row, stats, "missing_source", str(exc))
        return {"external_id": row["external_id"], "status": "missing_source", "message": str(exc)}


def raw_file_name(row: dict[str, str], default_suffix: str = ".dat") -> str:
    notes = parse_notes_kv(row.get("notes", ""))
    if row.get("source_file"):
        return Path(row["source_file"]).name
    if notes.get("file"):
        return Path(notes["file"]).name
    if notes.get("file_name"):
        return Path(notes["file_name"]).name
    source_url = row.get("source_url", "")
    if source_url:
        parsed = urllib.parse.urlparse(source_url)
        name = Path(urllib.parse.unquote(parsed.path)).name
        if name:
            return name
    return f"{row['external_id']}{default_suffix}"


def process_direct_url(row: dict[str, str], args: argparse.Namespace, panel: PanPanel, rsid_map: dict[str, PanSnp] | None) -> dict[str, object]:
    out = args.out_dir / "aligned" / f"{row['external_id']}.sumstats.gz"
    stats = args.out_dir / "prepare_stats" / f"{row['external_id']}.json"
    if skip := existing_skip_result(row, out, stats, args.force):
        return skip
    raw_name = raw_file_name(row)
    row["source_file"] = raw_name
    if not row.get("source_url"):
        message = (
            "Missing restricted source_url. Add it locally to "
            f"{args.url_overrides} with columns external_id and source_url."
        )
        write_skip_stats(row, stats, "restricted_url_missing", message)
        return {
            "external_id": row["external_id"],
            "status": "restricted_url_missing",
            "source_file": raw_name,
            "message": message,
        }
    if args.max_rows is not None:
        if is_tar_archive_name(raw_name) or is_zip_archive_name(raw_name):
            raw = download_file(
                row["source_url"],
                args.out_dir / "raw" / row["source_type"] / row["external_id"] / raw_name,
                args.force,
            )
            opener = open_zip_data_text if is_zip_archive_path(raw) else open_tar_data_text
            with opener(raw, row) as (text, member_name):
                row["source_file"] = f"{raw_name}:{member_name}"
                return convert_external_stream(text, row, panel, rsid_map, out, stats, args.max_rows)
        with open_text_url(row["source_url"], raw_name) as text:
            return convert_external_stream(text, row, panel, rsid_map, out, stats, args.max_rows)
    raw = download_file(
        row["source_url"],
        args.out_dir / "raw" / row["source_type"] / row["external_id"] / raw_name,
        args.force,
    )

    def convert_raw(raw_path: Path, tolerate_compressed_eof: bool = False) -> dict[str, object]:
        if is_tar_archive_path(raw_path):
            with open_tar_data_text(raw_path, row) as (text, member_name):
                row["source_file"] = f"{raw_name}:{member_name}"
                return convert_external_stream(
                    text,
                    row,
                    panel,
                    rsid_map,
                    out,
                    stats,
                    args.max_rows,
                    tolerate_compressed_eof,
                )
        if is_zip_archive_path(raw_path):
            with open_zip_data_text(raw_path, row) as (text, member_name):
                row["source_file"] = f"{raw_name}:{member_name}"
                return convert_external_stream(
                    text,
                    row,
                    panel,
                    rsid_map,
                    out,
                    stats,
                    args.max_rows,
                    tolerate_compressed_eof,
                )
        row["source_file"] = raw_name
        with open_text_path(raw_path, tolerate_compressed_eof) as text:
            return convert_external_stream(
                text,
                row,
                panel,
                rsid_map,
                out,
                stats,
                args.max_rows,
                tolerate_compressed_eof,
            )

    try:
        return convert_raw(raw)
    except COMPRESSED_READ_ERRORS:
        if args.force:
            raise
        raw.unlink(missing_ok=True)
        raw = download_file(
            row["source_url"],
            args.out_dir / "raw" / row["source_type"] / row["external_id"] / raw_name,
            force=True,
        )
        try:
            return convert_raw(raw)
        except COMPRESSED_READ_ERRORS:
            return convert_raw(raw, tolerate_compressed_eof=True)


def process_zenodo_face_cgwas(row: dict[str, str], args: argparse.Namespace, panel: PanPanel) -> dict[str, object]:
    out = args.out_dir / "aligned" / f"{row['external_id']}.sumstats.gz"
    stats = args.out_dir / "prepare_stats" / f"{row['external_id']}.json"
    if skip := existing_skip_result(row, out, stats, args.force):
        return skip

    notes = parse_notes_kv(row.get("notes", ""))
    archive_name = notes.get("archive_name")
    member_name = notes.get("tar_member")
    snp_info_file = notes.get("snp_info_file", "SnpInfo.tsv.bz2")
    snp_info_url = notes.get("snp_info_url")
    record_url = getattr(args, "zenodo_face_record_url", ZENODO_FACE_CGWAS_RECORD_URL)
    if not archive_name or not member_name or not snp_info_url:
        raise ValueError(f"Zenodo facial-shape row is missing archive/member/SNP metadata: {row['external_id']}")

    archive = download_file(
        row["source_url"],
        zenodo_face_raw_dir(args.out_dir, record_url) / archive_name,
        args.force,
    )
    snp_info = download_file(
        snp_info_url,
        zenodo_face_raw_dir(args.out_dir, record_url) / snp_info_file,
        args.force,
    )

    started = time.time()
    out.parent.mkdir(parents=True, exist_ok=True)
    stats.parent.mkdir(parents=True, exist_ok=True)
    out_tmp = out.with_name(f"{out.name}.tmp.{os.getpid()}")
    stats_tmp = stats.with_name(f"{stats.name}.tmp.{os.getpid()}")

    n_in = 0
    n_written = 0
    n_duplicate = 0
    n_missing = 0
    n_no_z = 0
    n_no_n = 0
    n_no_pan_match = 0
    seen: set[str] = set()
    sample_size = row.get("sample_size", "10115")
    sample_n = as_float(sample_size)
    if sample_n is None or sample_n <= 0:
        raise ValueError(f"Zenodo facial-shape row has invalid sample_size={sample_size}")

    try:
        with (
            bz2.open(snp_info, "rt") as snp_text,
            tarfile.open(archive, "r:*") as tar,
            gzip.open(out_tmp, "wt", compresslevel=6) as dst,
        ):
            member = tar.getmember(member_name)
            extracted = tar.extractfile(member)
            if extracted is None:
                raise FileNotFoundError(f"{member_name} is not a regular file in {archive}")
            try:
                assoc_text = io.TextIOWrapper(extracted)
                snp_header = snp_text.readline().strip().split()
                assoc_header = assoc_text.readline().strip().split()
                snp_index = {name: idx for idx, name in enumerate(snp_header)}
                assoc_index = {name: idx for idx, name in enumerate(assoc_header)}
                for required in ["CHR", "BP", "A1", "A2"]:
                    if required not in snp_index:
                        raise ValueError(f"{snp_info} is missing required column {required}")
                if "Beta" not in assoc_index or "P" not in assoc_index:
                    raise ValueError(f"{member_name} is missing required Beta/P columns")

                dst.write("SNP\tA1\tA2\tZ\tN\n")
                for snp_line, assoc_line in zip(snp_text, assoc_text):
                    n_in += 1
                    if args.max_rows is not None and n_in > args.max_rows:
                        break
                    snp_parts = snp_line.strip().split()
                    assoc_parts = assoc_line.strip().split()
                    if len(snp_parts) < len(snp_header) or len(assoc_parts) < len(assoc_header):
                        n_missing += 1
                        continue
                    beta = as_float(assoc_parts[assoc_index["Beta"]])
                    p = as_float(assoc_parts[assoc_index["P"]])
                    if beta is None or p is None:
                        n_no_z += 1
                        continue
                    z = signed_z_from_beta_p(beta, p)
                    if z is None:
                        n_no_z += 1
                        continue
                    oriented = orient_to_pan_alt(
                        panel,
                        snp_parts[snp_index["CHR"]],
                        snp_parts[snp_index["BP"]],
                        snp_parts[snp_index["A1"]],
                        snp_parts[snp_index["A2"]],
                        z,
                    )
                    if oriented is None:
                        n_no_pan_match += 1
                        continue
                    pan, z_out = oriented
                    if pan.snp in seen:
                        n_duplicate += 1
                        continue
                    seen.add(pan.snp)
                    dst.write(f"{pan.snp}\t{pan.alt}\t{pan.ref}\t{z_out:.8g}\t{sample_n:.0f}\n")
                    n_written += 1
            finally:
                extracted.close()

        stats_obj: dict[str, object] = {
            "external_id": row["external_id"],
            "source_type": row["source_type"],
            "trait": row.get("trait", ""),
            "status": "done",
            "population": row.get("population", ""),
            "source_url": row.get("source_url", ""),
            "source_file": f"{archive_name}:{member_name}",
            "genome_build": notes.get("genome_build", "GRCh37"),
            "input_rows": n_in,
            "written_rows": n_written,
            "total_panukbb_snps": len(panel.by_snp),
            "panukbb_overlap_pct": (100 * n_written / len(panel.by_snp)) if panel.by_snp else None,
            "input_overlap_pct": (100 * n_written / n_in) if n_in else None,
            "duplicate_snps": n_duplicate,
            "missing_alleles": n_missing,
            "missing_z": n_no_z,
            "missing_n": n_no_n,
            "no_panukbb_match": n_no_pan_match,
            "elapsed_seconds": time.time() - started,
            "out": str(out),
            "header": ["SnpInfo.tsv.bz2: CHR BP SNP A1 A2", f"{member_name}: Beta P"],
        }
        with stats_tmp.open("w") as f:
            json.dump(stats_obj, f, indent=2, sort_keys=True)
            f.write("\n")
        os.replace(out_tmp, out)
        os.replace(stats_tmp, stats)
        return stats_obj
    finally:
        for tmp in [out_tmp, stats_tmp]:
            if tmp.exists():
                tmp.unlink()


def process_nextcloud_file(row: dict[str, str], args: argparse.Namespace, panel: PanPanel, rsid_map: dict[str, PanSnp] | None) -> dict[str, object]:
    out = args.out_dir / "aligned" / f"{row['external_id']}.sumstats.gz"
    stats = args.out_dir / "prepare_stats" / f"{row['external_id']}.json"
    if skip := existing_skip_result(row, out, stats, args.force):
        return skip
    notes = parse_notes_kv(row.get("notes", ""))
    raw_name = row.get("source_file") or notes.get("file") or raw_file_name(row)
    row["source_file"] = raw_name
    if notes.get("metadata_only", "").lower() in {"1", "true", "yes"}:
        raw = args.out_dir / "raw" / "nextcloud" / row["external_id"] / raw_name
        if args.max_rows is None:
            if not raw.exists() or raw.stat().st_size == 0 or args.force:
                raw.parent.mkdir(parents=True, exist_ok=True)
                tmp = raw.with_name(f"{raw.name}.tmp.{os.getpid()}")
                with urllib.request.urlopen(nextcloud_request(row), timeout=120) as resp, tmp.open("wb") as f:
                    shutil.copyfileobj(resp, f, length=1024 * 1024)
                os.replace(tmp, raw)
            message = "Downloaded metadata-only file; it has no association Z/BETA/SE/P columns to align."
            result = write_metadata_only_stats(row, stats, raw, message)
        else:
            message = "Metadata-only file; raw download skipped because --max-rows was set."
            result = write_metadata_only_stats(row, stats, None, message)
        return {
            "external_id": row["external_id"],
            "status": result["status"],
            "source_file": raw_name,
            "message": result["message"],
        }
    if args.max_rows is not None:
        with open_text_request(nextcloud_request(row), raw_name) as text:
            return convert_external_stream(text, row, panel, rsid_map, out, stats, args.max_rows)
    raw = args.out_dir / "raw" / "nextcloud" / row["external_id"] / raw_name
    if not raw.exists() or raw.stat().st_size == 0 or args.force:
        raw.parent.mkdir(parents=True, exist_ok=True)
        tmp = raw.with_name(f"{raw.name}.tmp.{os.getpid()}")
        with urllib.request.urlopen(nextcloud_request(row), timeout=120) as resp, tmp.open("wb") as f:
            shutil.copyfileobj(resp, f, length=1024 * 1024)
        os.replace(tmp, raw)
    if is_tar_archive_path(raw):
        with open_tar_data_text(raw, row) as (text, member_name):
            row["source_file"] = f"{raw_name}:{member_name}"
            return convert_external_stream(text, row, panel, rsid_map, out, stats, args.max_rows)
    with open_text_path(raw) as text:
        return convert_external_stream(text, row, panel, rsid_map, out, stats, args.max_rows)


def process_figshare_file(row: dict[str, str], args: argparse.Namespace, panel: PanPanel, rsid_map: dict[str, PanSnp] | None) -> dict[str, object]:
    out = args.out_dir / "aligned" / f"{row['external_id']}.sumstats.gz"
    stats = args.out_dir / "prepare_stats" / f"{row['external_id']}.json"
    if skip := existing_skip_result(row, out, stats, args.force):
        return skip
    notes = parse_notes_kv(row.get("notes", ""))
    download_url = notes.get("download_url")
    file_name = notes.get("file_name") or row.get("source_file")
    if not download_url:
        doi = row.get("data_doi") or notes.get("data_doi")
        if not doi or not doi.startswith("10.6084/m9.figshare."):
            write_skip_stats(row, stats, "manual_or_unsupported", "No public Figshare DOI or download_url")
            return {"external_id": row["external_id"], "status": "manual_or_unsupported"}
        article = figshare_article_from_doi(doi)
        file_info = choose_figshare_file(article)
        download_url = str(file_info["download_url"])
        file_name = str(file_info.get("name") or file_name or f"{row['external_id']}.dat")
    row["source_file"] = file_name or raw_file_name(row)
    if args.max_rows is not None:
        if is_tar_archive_name(row["source_file"]):
            raw = download_file(
                download_url,
                args.out_dir / "raw" / "figshare" / row["external_id"] / row["source_file"],
                args.force,
            )
            with open_tar_data_text(raw, row) as (text, member_name):
                row["source_file"] = f"{row['source_file']}:{member_name}"
                return convert_external_stream(text, row, panel, rsid_map, out, stats, args.max_rows)
        with open_text_url(download_url, row["source_file"]) as text:
            return convert_external_stream(text, row, panel, rsid_map, out, stats, args.max_rows)
    raw = download_file(
        download_url,
        args.out_dir / "raw" / "figshare" / row["external_id"] / row["source_file"],
        args.force,
    )
    if is_tar_archive_path(raw):
        with open_tar_data_text(raw, row) as (text, member_name):
            row["source_file"] = f"{row['source_file']}:{member_name}"
            return convert_external_stream(text, row, panel, rsid_map, out, stats, args.max_rows)
    with open_text_path(raw) as text:
        return convert_external_stream(text, row, panel, rsid_map, out, stats, args.max_rows)


def process_pgc(row: dict[str, str], args: argparse.Namespace, panel: PanPanel, rsid_map: dict[str, PanSnp]) -> dict[str, object]:
    out = args.out_dir / "aligned" / f"{row['external_id']}.sumstats.gz"
    stats = args.out_dir / "prepare_stats" / f"{row['external_id']}.json"
    if skip := existing_skip_result(row, out, stats, args.force):
        return skip
    doi = row.get("data_doi", "")
    notes = parse_notes_kv(row.get("notes", ""))
    download_url = notes.get("download_url")
    name = notes.get("file_name") or row.get("source_file") or f"{row['external_id']}.dat"
    if not download_url:
        if not doi.startswith("10.6084/m9.figshare."):
            write_skip_stats(row, stats, "manual_or_unsupported", "PGC row does not point to a public Figshare DOI")
            return {"external_id": row["external_id"], "status": "manual_or_unsupported"}
        article = figshare_article_from_doi(doi)
        file_info = choose_figshare_file(article)
        download_url = str(file_info["download_url"])
        name = str(file_info.get("name") or name)
    row["source_file"] = name
    if args.max_rows is not None:
        if is_tar_archive_name(name):
            raw = download_file(download_url, args.out_dir / "raw" / "pgc" / row["external_id"] / name, args.force)
            with open_tar_data_text(raw, row) as (text, member_name):
                row["source_file"] = f"{name}:{member_name}"
                return convert_external_stream(text, row, panel, rsid_map, out, stats, args.max_rows)
        with open_text_url(download_url, name) as text:
            return convert_external_stream(text, row, panel, rsid_map, out, stats, args.max_rows)
    raw = download_file(download_url, args.out_dir / "raw" / "pgc" / row["external_id"] / name, args.force)
    if is_tar_archive_path(raw):
        with open_tar_data_text(raw, row) as (text, member_name):
            row["source_file"] = f"{name}:{member_name}"
            return convert_external_stream(text, row, panel, rsid_map, out, stats, args.max_rows)
    with open_text_path(raw) as text:
        return convert_external_stream(text, row, panel, rsid_map, out, stats, args.max_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", required=True, type=Path)
    parser.add_argument("--ld-snps", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--sources", default="gwas_catalog,zenodo_indep107,pgc,direct_url,cncr,egg")
    parser.add_argument("--include", default="")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--max-rows", type=int)
    parser.add_argument("--zenodo-archive")
    parser.add_argument("--rsid-reference")
    parser.add_argument("--url-overrides", type=Path, default=Path("config/external_gwas_restricted_urls.tsv"))
    parser.add_argument("--pgc-page-url", default=PGC_DOWNLOADS_URL)
    parser.add_argument("--cncr-page-url", default=CNCR_SUMMARY_STATS_URL)
    parser.add_argument("--egg-page-url", default=EGG_SUMMARY_STATS_URL)
    parser.add_argument("--zenodo-face-record-url", default=ZENODO_FACE_CGWAS_RECORD_URL)
    parser.add_argument("--manifest-only", action="store_true")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    static_rows = load_targets(args.targets)
    all_rows = list(static_rows)
    sources = split_list(args.sources)
    if "pgc" in sources:
        all_rows.extend(fetch_pgc_rows(args.pgc_page_url))
    if "cncr" in sources:
        all_rows.extend(fetch_cncr_rows(args.cncr_page_url))
    if "egg" in sources:
        all_rows.extend(fetch_egg_rows(args.egg_page_url))
    if "zenodo_face_cgwas" in sources:
        all_rows.extend(
            fetch_zenodo_face_rows(
                args.zenodo_face_record_url,
                args.out_dir,
                args.force,
                args.limit if not args.include else None,
            )
        )
    apply_url_overrides(all_rows, args.url_overrides)

    selected = selected_rows(all_rows, args)
    manifest = args.out_dir / "catalog" / "external_gwas_manifest.tsv"
    write_external_manifest(manifest, selected, args.out_dir / "prepare_stats")
    print(f"wrote {len(selected)} rows to {manifest}", flush=True)
    if args.manifest_only:
        return

    if not args.ld_snps.exists():
        raise FileNotFoundError(
            f"{args.ld_snps} is missing. Run make setup or make prepare-ldscores first."
        )
    panel = load_pan_panel(args.ld_snps)
    needs_rsid_map = any(row["source_type"] not in {"gwas_catalog", "zenodo_face_cgwas"} for row in selected)
    rsid_map = ensure_rsid_map(args, panel) if needs_rsid_map else None

    failures: list[str] = []
    for row in selected:
        try:
            source = row["source_type"]
            if source == "gwas_catalog":
                result = process_gwas_catalog(row, args, panel)
            elif source == "zenodo_indep107":
                if rsid_map is None:
                    raise RuntimeError("rsID map was not initialized")
                result = process_zenodo(row, args, panel, rsid_map)
            elif source == "pgc":
                if rsid_map is None:
                    raise RuntimeError("rsID map was not initialized")
                result = process_pgc(row, args, panel, rsid_map)
            elif source == "direct_url":
                result = process_direct_url(row, args, panel, rsid_map)
            elif source == "zenodo_face_cgwas":
                result = process_zenodo_face_cgwas(row, args, panel)
            elif source == "figshare_file":
                result = process_figshare_file(row, args, panel, rsid_map)
            elif source == "nextcloud_file":
                result = process_nextcloud_file(row, args, panel, rsid_map)
            elif source == "egg":
                result = process_direct_url(row, args, panel, rsid_map)
            else:
                raise ValueError(f"Unsupported source_type={source}")
            print(json.dumps(result, sort_keys=True), flush=True)
            write_external_manifest(manifest, selected, args.out_dir / "prepare_stats")
            if args.strict and result.get("status") not in {"done", "skip", "metadata_only"}:
                failures.append(f"{row['external_id']}: {result.get('status')}")
            if args.strict and int(result.get("written_rows", 1) or 0) == 0 and result.get("status") == "done":
                failures.append(f"{row['external_id']}: wrote zero rows")
        except Exception as exc:  # noqa: BLE001 - keep batch resumable and write diagnostics.
            stats = args.out_dir / "prepare_stats" / f"{row['external_id']}.json"
            write_skip_stats(row, stats, "failed", str(exc))
            write_external_manifest(manifest, selected, args.out_dir / "prepare_stats")
            message = f"{row['external_id']}: failed: {exc}"
            print(message, file=sys.stderr, flush=True)
            failures.append(message)
            if args.strict:
                break

    if failures and args.strict:
        raise SystemExit("External GWAS preparation failed:\n" + "\n".join(failures))


if __name__ == "__main__":
    main()
