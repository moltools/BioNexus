#!/usr/bin/env python3

"""Stream GenBank files to map nuccore accessions to GCF accessions via NCBI E-utilities."""

from __future__ import annotations

import argparse
import csv
import logging
import os
import re
import sys
import time
from collections import OrderedDict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass

import requests
from tqdm import tqdm

from bionexus.utils.logging import setup_logging


log = logging.getLogger(__name__)

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

ACCESSION_RE = re.compile(r"^[A-Z]{2}_[A-Z0-9]+\d+(?:\.\d+)?$")
HEADER_STOP_MARKERS = ("FEATURES", "ORIGIN")
EUTILS_TERM_CHUNK_SIZE = 200
EUTILS_ID_CHUNK_SIZE = 200


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Map nuccore accessions in GenBank files to GCF accessions."
    )
    parser.add_argument("--gbk-dir", required=True, help="directory with .gbk/.gbff files")
    parser.add_argument("--output", required=True, help="output CSV path")
    parser.add_argument(
        "--extensions",
        default=".gbk,.gbff,.gb",
        help="comma-separated file extensions to scan (default: .gbk,.gbff,.gb)",
    )
    parser.add_argument(
        "--max-header-lines",
        type=int,
        default=2000,
        help="max header lines to scan per file (default: 2000)",
    )
    parser.add_argument(
        "--keep-version",
        action="store_true",
        help="keep version suffix (e.g., .1) on accessions",
    )
    parser.add_argument(
        "--skip-missing",
        action="store_true",
        help="skip writing rows when a GCF accession is not found",
    )
    parser.add_argument(
        "--append",
        action="store_true",
        help="append to output CSV instead of overwriting",
    )
    parser.add_argument(
        "--log-every",
        type=int,
        default=5000,
        help="log progress every N files (default: 5000)",
    )
    parser.add_argument(
        "--cache-size",
        type=int,
        default=10000,
        help="LRU cache size for accession lookups (0 disables; default: 10000)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=200,
        help="accessions per NCBI batch (default: 200)",
    )
    parser.add_argument(
        "--rate",
        type=float,
        default=3.0,
        help="max NCBI requests per second (default: 3)",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=20,
        help="NCBI request timeout in seconds (default: 20)",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="NCBI request retries on failure (default: 3)",
    )
    parser.add_argument(
        "--email",
        default=os.getenv("NCBI_EMAIL", ""),
        help="NCBI contact email (or set NCBI_EMAIL)",
    )
    parser.add_argument(
        "--api-key",
        default=os.getenv("NCBI_API_KEY", ""),
        help="NCBI API key (or set NCBI_API_KEY)",
    )
    parser.add_argument(
        "--tool",
        default="bionexus-nuccore-gcf",
        help="NCBI tool name (default: bionexus-nuccore-gcf)",
    )
    parser.add_argument(
        "--allow-gca",
        action="store_true",
        help="allow returning GCA accessions when no GCF is available",
    )
    parser.add_argument(
        "--progress",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="show progress bar (default: enabled when stderr is a TTY)",
    )
    return parser.parse_args()


def iter_gbk_files(root: str, extensions: tuple[str, ...]) -> Iterator[str]:
    for dirpath, _, filenames in os.walk(root):
        for filename in filenames:
            if filename.lower().endswith(extensions):
                yield os.path.join(dirpath, filename)


def _strip_version(token: str) -> str:
    if "." in token:
        base, suffix = token.rsplit(".", 1)
        if suffix.isdigit():
            return base
    return token


def _unique_preserve(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


def _chunked(values: list[str], size: int) -> Iterator[list[str]]:
    if size <= 0:
        raise ValueError("chunk size must be positive")
    for idx in range(0, len(values), size):
        yield values[idx : idx + size]


def _find_accession_in_line(line: str, *, strip_version: bool) -> str | None:
    for raw in line.split():
        if raw in {"LOCUS", "ACCESSION", "VERSION"}:
            continue
        token = raw.strip().strip(";").strip(",")
        token = _strip_version(token) if strip_version else token
        if ACCESSION_RE.match(token):
            return token
    return None


def extract_nuccore_accession(
    path: str,
    *,
    max_header_lines: int,
    strip_version: bool,
) -> str | None:
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for idx, line in enumerate(handle):
                if idx >= max_header_lines:
                    break

                line = line.rstrip("\n")
                if not line:
                    continue

                if line.startswith("VERSION"):
                    acc = _find_accession_in_line(line, strip_version=strip_version)
                    if acc:
                        return acc

                if line.startswith("ACCESSION") or line.startswith("LOCUS"):
                    acc = _find_accession_in_line(line, strip_version=strip_version)
                    if acc:
                        return acc

                if line.startswith(HEADER_STOP_MARKERS):
                    break
    except OSError as exc:
        log.warning("Failed to read %s: %s", path, exc)
    return None


class RateLimiter:
    def __init__(self, rate_per_sec: float) -> None:
        self.min_interval = 1.0 / rate_per_sec if rate_per_sec > 0 else 0.0
        self._last_ts = 0.0

    def wait(self) -> None:
        if self.min_interval <= 0:
            return
        now = time.monotonic()
        if self._last_ts:
            delay = self.min_interval - (now - self._last_ts)
            if delay > 0:
                time.sleep(delay)
        self._last_ts = time.monotonic()


class LRUCache:
    def __init__(self, maxsize: int) -> None:
        self.maxsize = maxsize
        self._data: OrderedDict[str, str | None] = OrderedDict()

    def get(self, key: str) -> object:
        if self.maxsize <= 0:
            return _MISSING
        if key not in self._data:
            return _MISSING
        value = self._data.pop(key)
        self._data[key] = value
        return value

    def set(self, key: str, value: str | None) -> None:
        if self.maxsize <= 0:
            return
        if key in self._data:
            self._data.pop(key)
        self._data[key] = value
        if len(self._data) > self.maxsize:
            self._data.popitem(last=False)


_MISSING = object()


@dataclass
class NcbiClient:
    api_key: str
    email: str
    tool: str
    rate: float
    timeout: int
    retries: int

    def __post_init__(self) -> None:
        self.session = requests.Session()
        ua = self.tool
        if self.email:
            ua = f"{ua} ({self.email})"
        self.session.headers.update({"User-Agent": ua})
        self.limiter = RateLimiter(self.rate)

    def _request_json(
        self,
        endpoint: str,
        params: dict[str, str],
        *,
        method: str = "GET",
    ) -> dict | None:
        merged = {"retmode": "json", "tool": self.tool}
        if self.email:
            merged["email"] = self.email
        if self.api_key:
            merged["api_key"] = self.api_key
        merged.update(params)
        url = f"{EUTILS_BASE}/{endpoint}"
        method = method.upper()
        if method not in {"GET", "POST"}:
            raise ValueError(f"Unsupported method: {method}")

        backoff = 1.0
        last_err: Exception | None = None
        for _ in range(self.retries + 1):
            self.limiter.wait()
            try:
                if method == "POST":
                    resp = self.session.post(url, data=merged, timeout=self.timeout)
                else:
                    resp = self.session.get(url, params=merged, timeout=self.timeout)
            except requests.RequestException as exc:
                last_err = exc
            else:
                if resp.status_code == 200:
                    try:
                        return resp.json()
                    except ValueError as exc:
                        last_err = exc
                elif resp.status_code == 429:
                    retry_after = resp.headers.get("Retry-After")
                    if retry_after and retry_after.isdigit():
                        time.sleep(float(retry_after))
                    else:
                        time.sleep(backoff)
                        backoff = min(backoff * 2, 60.0)
                    continue
                elif 500 <= resp.status_code < 600:
                    last_err = RuntimeError(f"NCBI {endpoint} {resp.status_code}")
                else:
                    log.warning("NCBI %s returned %s for %s", endpoint, resp.status_code, merged)
                    return None

            time.sleep(backoff)
            backoff = min(backoff * 2, 60.0)

        log.warning("NCBI request failed for %s: %s", endpoint, last_err)
        return None

    def esearch_nuccore_id(self, accession: str) -> str | None:
        payload = {"db": "nuccore", "term": f"{accession}[Accession]"}
        data = self._request_json("esearch.fcgi", payload)
        if not data:
            return None
        ids = data.get("esearchresult", {}).get("idlist") or []
        return ids[0] if ids else None

    def esearch_nuccore_ids(self, accessions: list[str]) -> list[str]:
        if not accessions:
            return []
        term = " OR ".join(f"{acc}[Accession]" for acc in accessions)
        payload = {
            "db": "nuccore",
            "term": term,
            "retmax": str(len(accessions)),
        }
        data = self._request_json("esearch.fcgi", payload, method="POST")
        if not data:
            return []
        return data.get("esearchresult", {}).get("idlist") or []

    def elink_assembly_ids(self, nuccore_id: str) -> list[str]:
        payload = {"dbfrom": "nuccore", "db": "assembly", "id": nuccore_id}
        data = self._request_json("elink.fcgi", payload)
        if not data:
            return []
        linksets = data.get("linksets") or []
        for linkset in linksets:
            for db in linkset.get("linksetdbs") or []:
                if db.get("dbto") == "assembly":
                    links = db.get("links") or []
                    return [str(link) for link in links]
        return []

    def elink_assembly_map(self, nuccore_ids: list[str]) -> dict[str, list[str]]:
        if not nuccore_ids:
            return {}
        payload = {
            "dbfrom": "nuccore",
            "db": "assembly",
            "id": ",".join(nuccore_ids),
        }
        data = self._request_json("elink.fcgi", payload, method="POST")
        if not data:
            return {}
        mapping: dict[str, list[str]] = {}
        linksets = data.get("linksets") or []
        for linkset in linksets:
            src_ids = [str(val) for val in (linkset.get("ids") or [])]
            if not src_ids:
                continue
            assembly_ids: list[str] = []
            for db in linkset.get("linksetdbs") or []:
                if db.get("dbto") == "assembly":
                    assembly_ids = [str(link) for link in (db.get("links") or [])]
                    break
            for src_id in src_ids:
                mapping[src_id] = assembly_ids
        return mapping

    def esummary_assembly(self, assembly_id: str) -> dict | None:
        payload = {"db": "assembly", "id": assembly_id}
        data = self._request_json("esummary.fcgi", payload)
        if not data:
            return None
        result = data.get("result") or {}
        return result.get(str(assembly_id))

    def esummary_nuccore(self, nuccore_ids: list[str]) -> dict[str, dict]:
        if not nuccore_ids:
            return {}
        payload = {"db": "nuccore", "id": ",".join(nuccore_ids)}
        data = self._request_json("esummary.fcgi", payload, method="POST")
        if not data:
            return {}
        result = data.get("result") or {}
        summaries: dict[str, dict] = {}
        for uid in result.get("uids") or []:
            summary = result.get(uid)
            if summary:
                summaries[str(uid)] = summary
        return summaries

    def esummary_assembly_batch(self, assembly_ids: list[str]) -> dict[str, dict]:
        if not assembly_ids:
            return {}
        payload = {"db": "assembly", "id": ",".join(assembly_ids)}
        data = self._request_json("esummary.fcgi", payload, method="POST")
        if not data:
            return {}
        result = data.get("result") or {}
        summaries: dict[str, dict] = {}
        for uid in result.get("uids") or []:
            summary = result.get(uid)
            if summary:
                summaries[str(uid)] = summary
        return summaries


def _pick_gcf(value: object, *, allow_gca: bool) -> str | None:
    if isinstance(value, str):
        if value.startswith("GCF_"):
            return value
        if allow_gca and value.startswith("GCA_"):
            return value
    elif isinstance(value, Iterable):
        for item in value:
            gcf = _pick_gcf(item, allow_gca=allow_gca)
            if gcf:
                return gcf
    return None


def extract_gcf_from_summary(summary: dict, *, allow_gca: bool) -> str | None:
    for key in ("assemblyaccession", "assembly_accession", "assemblyacc", "assembly_acc"):
        gcf = _pick_gcf(summary.get(key), allow_gca=allow_gca)
        if gcf:
            return gcf

    synonym = summary.get("synonym") or summary.get("synonyms")
    if isinstance(synonym, dict):
        for key in ("refseq", "refseqaccession", "refseq_accession"):
            gcf = _pick_gcf(synonym.get(key), allow_gca=allow_gca)
            if gcf:
                return gcf
        for value in synonym.values():
            gcf = _pick_gcf(value, allow_gca=allow_gca)
            if gcf:
                return gcf
    else:
        gcf = _pick_gcf(synonym, allow_gca=allow_gca)
        if gcf:
            return gcf
    return None


def _summary_accession(summary: dict, *, strip_version: bool) -> str | None:
    for key in ("accessionversion", "accession", "caption"):
        value = summary.get(key)
        if value:
            value = str(value)
            return _strip_version(value) if strip_version else value
    return None


def nuccore_to_gcf(
    client: NcbiClient,
    accession: str,
    *,
    allow_gca: bool,
) -> str | None:
    nuccore_id = client.esearch_nuccore_id(accession)
    if not nuccore_id:
        return None
    assembly_ids = client.elink_assembly_ids(nuccore_id)
    for assembly_id in assembly_ids:
        summary = client.esummary_assembly(assembly_id)
        if not summary:
            continue
        gcf = extract_gcf_from_summary(summary, allow_gca=allow_gca)
        if gcf:
            return gcf
    return None


def resolve_accessions_batch(
    client: NcbiClient,
    accessions: list[str],
    *,
    strip_version: bool,
    allow_gca: bool,
) -> dict[str, str | None]:
    results: dict[str, str | None] = {acc: None for acc in accessions}
    if not accessions:
        return results

    all_ids: list[str] = []
    for chunk in _chunked(accessions, EUTILS_TERM_CHUNK_SIZE):
        ids = client.esearch_nuccore_ids(chunk)
        all_ids.extend(ids)

    all_ids = _unique_preserve(all_ids)
    if not all_ids:
        return results

    summaries: dict[str, dict] = {}
    for chunk in _chunked(all_ids, EUTILS_ID_CHUNK_SIZE):
        summaries.update(client.esummary_nuccore(chunk))

    acc_to_id: dict[str, str] = {}
    for uid, summary in summaries.items():
        acc = _summary_accession(summary, strip_version=strip_version)
        if not acc:
            continue
        target = None
        if acc in results:
            target = acc
        elif not strip_version:
            base = _strip_version(acc)
            if base in results:
                target = base
        if target and target not in acc_to_id:
            acc_to_id[target] = uid

    if not acc_to_id:
        return results

    assembly_map: dict[str, list[str]] = {}
    nuccore_ids = list(acc_to_id.values())
    for chunk in _chunked(nuccore_ids, EUTILS_ID_CHUNK_SIZE):
        assembly_map.update(client.elink_assembly_map(chunk))

    assembly_ids: list[str] = []
    seen_assemblies: set[str] = set()
    for ids_list in assembly_map.values():
        for assembly_id in ids_list:
            if assembly_id in seen_assemblies:
                continue
            seen_assemblies.add(assembly_id)
            assembly_ids.append(assembly_id)

    assembly_summaries: dict[str, dict] = {}
    for chunk in _chunked(assembly_ids, EUTILS_ID_CHUNK_SIZE):
        assembly_summaries.update(client.esummary_assembly_batch(chunk))

    nuccore_to_gcf_map: dict[str, str] = {}
    for nuccore_id, ids_list in assembly_map.items():
        for assembly_id in ids_list:
            summary = assembly_summaries.get(assembly_id)
            if not summary:
                continue
            gcf = extract_gcf_from_summary(summary, allow_gca=allow_gca)
            if gcf:
                nuccore_to_gcf_map[nuccore_id] = gcf
                break

    for acc, uid in acc_to_id.items():
        results[acc] = nuccore_to_gcf_map.get(uid)

    return results


def _needs_header(path: str, append: bool) -> bool:
    if not append:
        return True
    if not os.path.exists(path):
        return True
    return os.path.getsize(path) == 0


def main() -> None:
    args = cli()
    setup_logging("INFO")

    extensions = tuple(ext.strip().lower() for ext in args.extensions.split(",") if ext.strip())
    if not extensions:
        raise SystemExit("No valid extensions specified.")
    if args.batch_size < 1:
        raise SystemExit("batch-size must be >= 1")
    show_progress = args.progress
    if show_progress is None:
        show_progress = sys.stderr.isatty()

    cache = LRUCache(args.cache_size)
    client = NcbiClient(
        api_key=args.api_key,
        email=args.email,
        tool=args.tool,
        rate=args.rate,
        timeout=args.timeout,
        retries=args.retries,
    )

    write_header = _needs_header(args.output, args.append)
    mode = "a" if args.append else "w"

    files_seen = 0
    rows_written = 0
    missing_accessions = 0
    pending: list[str] = []
    pbar = tqdm(desc="GBK files", unit="file", disable=not show_progress)

    with open(args.output, mode, newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        if write_header:
            writer.writerow(["nuccore_accession", "gcf_accession"])

        def flush_pending() -> None:
            nonlocal rows_written
            if not pending:
                return

            unique = _unique_preserve(pending)
            to_lookup: list[str] = []
            for acc in unique:
                if cache.get(acc) is _MISSING:
                    to_lookup.append(acc)

            batch_results: dict[str, str | None] = {}
            if to_lookup:
                batch_results = resolve_accessions_batch(
                    client,
                    to_lookup,
                    strip_version=not args.keep_version,
                    allow_gca=args.allow_gca,
                )
                for acc, gcf in batch_results.items():
                    cache.set(acc, gcf)

            for acc in pending:
                cached = cache.get(acc)
                if cached is _MISSING:
                    gcf = batch_results.get(acc)
                else:
                    gcf = cached
                if gcf or not args.skip_missing:
                    writer.writerow([acc, gcf or ""])
                    rows_written += 1

            pending.clear()

        for path in iter_gbk_files(args.gbk_dir, extensions):
            files_seen += 1
            pbar.update(1)
            accession = extract_nuccore_accession(
                path,
                max_header_lines=args.max_header_lines,
                strip_version=not args.keep_version,
            )
            if not accession:
                missing_accessions += 1
                if files_seen % args.log_every == 0:
                    if show_progress:
                        pbar.set_postfix(rows=rows_written, missing=missing_accessions)
                    else:
                        log.info(
                            "Processed %d files (rows=%d, missing=%d)",
                            files_seen,
                            rows_written,
                            missing_accessions,
                        )
                continue

            pending.append(accession)
            if len(pending) >= args.batch_size:
                flush_pending()

            if files_seen % args.log_every == 0:
                if show_progress:
                    pbar.set_postfix(rows=rows_written, missing=missing_accessions)
                else:
                    log.info(
                        "Processed %d files (rows=%d, missing=%d)",
                        files_seen,
                        rows_written,
                        missing_accessions,
                    )

        flush_pending()

    if show_progress:
        pbar.set_postfix(rows=rows_written, missing=missing_accessions)
    pbar.close()

    log.info(
        "Done. Files=%d, rows=%d, missing_accessions=%d",
        files_seen,
        rows_written,
        missing_accessions,
    )


if __name__ == "__main__":
    main()
