#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import json
import logging
import os
from pathlib import Path
import time
from typing import Dict, Iterable, List, Optional, Tuple
import random
import time
from urllib.parse import quote

import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi  # requires RDKit with InChI support
from tqdm import tqdm

RDLogger.DisableLog('rdApp.*')

NPATLAS_SDF_URL = "https://www.npatlas.org/static/downloads/NPAtlas_download.sdf"
UNICHEM_COMPOUND_API = "https://www.ebi.ac.uk/unichem/api/v1/compounds"
RHEA_SPARQL = "https://sparql.rhea-db.org/sparql"

# How many CHEBI IDs per SPARQL VALUES block (keep modest for reliability)
SPARQL_BATCH = 200
# Gentle pacing for web APIs (seconds)
UNICHEM_SLEEP = 0.02


def setup_logging(workdir: Path, level: str = "INFO") -> Path:
    """Configure console + file logging. Returns log path."""
    log_path = workdir / "parse_chebi.log"
    workdir.mkdir(parents=True, exist_ok=True)

    # Clear existing handlers (if re-running in same interpreter)
    root = logging.getLogger()
    for h in list(root.handlers):
        root.removeHandler(h)

    root.setLevel(getattr(logging, level.upper(), logging.INFO))
    fmt = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setFormatter(fmt)
    fh.setLevel(getattr(logging, level.upper(), logging.INFO))

    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    sh.setLevel(getattr(logging, level.upper(), logging.INFO))

    root.addHandler(fh)
    root.addHandler(sh)

    logging.info("Logging initialized. File: %s", log_path)
    return log_path


def ensure_sdf(workdir: Path, url: str = NPATLAS_SDF_URL) -> Path:
    log = logging.getLogger("ensure_sdf")
    workdir.mkdir(parents=True, exist_ok=True)
    sdf_path = workdir / "NPAtlas_download.sdf"

    if sdf_path.exists() and sdf_path.stat().st_size > 0:
        log.info("SDF already present: %s (%.2f MB)", sdf_path, sdf_path.stat().st_size / (1024 ** 2))
        return sdf_path

    log.info("Downloading NPAtlas SDF from %s ...", url)
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(sdf_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                if chunk:
                    f.write(chunk)
    log.info("Downloaded SDF to %s (%.2f MB)", sdf_path, sdf_path.stat().st_size / (1024 ** 2))
    return sdf_path


def iter_npatlas_mols(sdf_path: Path) -> Iterable[Tuple[int, Chem.Mol]]:
    log = logging.getLogger("iter_npatlas_mols")
    log.info("Opening SDF: %s", sdf_path)
    suppl = Chem.SDMolSupplier(str(sdf_path), sanitize=True, removeHs=False)
    for idx, mol in enumerate(suppl):
        if mol is None:
            continue
        yield idx, mol


def best_npatlas_id(mol: Chem.Mol, idx: int) -> str:
    candidates = [
        "npa_id", "npatlas_id", "npaid", "np_atlas_id", "id", "ID",
        "record_id", "accession", "NAME", "Name"
    ]
    for key in candidates:
        if mol.HasProp(key):
            val = mol.GetProp(key).strip()
            if val:
                return val
    name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
    return name if name else f"npatlas_idx_{idx}"


def mol_inchikey(mol: Chem.Mol) -> Optional[str]:
    try:
        return inchi.MolToInchiKey(mol)
    except Exception:
        logging.getLogger("mol_inchikey").exception("Failed to generate InChIKey for a molecule")
        return None


def unichem_sources_for_inchikey(inchikey: str, session: requests.Session) -> Optional[Dict]:
    """
    Query UniChem 2.0 compound endpoint with type='inchikey'.
    Returns the first compound record (dict) or None.
    """
    log = logging.getLogger("unichem")
    payload = {"type": "inchikey", "compound": inchikey}
    try:
        r = session.post(UNICHEM_COMPOUND_API, json=payload, timeout=30)
        r.raise_for_status()
        data = r.json()
        comps = data.get("compounds") or []
        if not comps:
            log.debug("No UniChem record for %s", inchikey)
            return None
        for c in comps:
            if c.get("standardInchiKey") == inchikey:
                return c
        return comps[0]
    except Exception:
        log.exception("UniChem request failed for %s", inchikey)
        return None


def extract_chebi_ids_from_sources(compound_record: Dict) -> List[str]:
    """
    From UniChem 'sources' list, find ChEBI entries.
    Prefer matching by source name (robust to sourceID changes).
    """
    out: List[str] = []
    sources = compound_record.get("sources") or []
    for s in sources:
        name = (s.get("shortName") or "").lower()
        if name == "chebi":
            cid = str(s.get("compoundId") or "").strip()
            if cid:
                if not cid.upper().startswith("CHEBI:"):
                    cid = f"CHEBI:{cid}"
                out.append(cid)
    return sorted(set(out))


def batched(iterable, n):
    batch = []
    for x in iterable:
        batch.append(x)
        if len(batch) >= n:
            yield batch
            batch = []
    if batch:
        yield batch

BASE = "https://www.ebi.ac.uk/chebi/backend/api/public/compound"
DEFAULT_HEADERS = {
    "User-Agent": "chebi-role-fetcher/1.0 (+https://example.org)",
    "Accept": "application/json",
}

def _normalize_chebi_id(cid: str) -> str:
    cid = cid.strip()
    return cid if cid.upper().startswith("CHEBI:") else f"CHEBI:{cid}"

def _build_url(chebi_id: str) -> str:
    enc = quote(chebi_id, safe="")
    return f"{BASE}/{enc}/?only_ontology_parents=false&only_ontology_children=false"

def _get_with_retries(
    session: requests.Session,
    url: str,
    *,
    max_retries: int = 5,
    timeout: float = 20.0,
    base_sleep: float = 0.8,
    max_sleep: float = 8.0,
) -> requests.Response:
    logger = logging.getLogger(__name__)
    attempt = 0
    while True:
        try:
            resp = session.get(url, headers=DEFAULT_HEADERS, timeout=timeout)
        except requests.RequestException:
            if attempt >= max_retries:
                logger.exception(f"Network error while fetching {url}")
                raise
            sleep = min(max_sleep, base_sleep * (2 ** attempt)) * (1 + random.random() * 0.2)
            logger.warning(f"Retrying {url} after network error, sleeping {sleep:.1f}s")
            time.sleep(sleep)
            attempt += 1
            continue

        if resp.status_code in (429, 503):
            if attempt >= max_retries:
                logger.error(f"Too many retries for {url} (status {resp.status_code})")
                resp.raise_for_status()
            retry_after = resp.headers.get("Retry-After")
            if retry_after:
                try:
                    sleep = float(retry_after)
                except ValueError:
                    sleep = min(max_sleep, base_sleep * (2 ** attempt))
            else:
                sleep = min(max_sleep, base_sleep * (2 ** attempt))
            sleep *= 1 + random.random() * 0.2
            logger.warning(f"Rate limited for {url}, sleeping {sleep:.1f}s")
            time.sleep(sleep)
            attempt += 1
            continue

        if 500 <= resp.status_code < 600:
            if attempt >= max_retries:
                logger.error(f"Server error for {url} (status {resp.status_code}) after retries")
                resp.raise_for_status()
            sleep = min(max_sleep, base_sleep * (2 ** attempt)) * (1 + random.random() * 0.2)
            logger.warning(f"Server error {resp.status_code} for {url}, retrying after {sleep:.1f}s")
            time.sleep(sleep)
            attempt += 1
            continue

        resp.raise_for_status()
        return resp


def get_biological_roles_for_chebis(
    chebi_ids: Iterable[str],
    *,
    per_request_pause: float = 0.2,
    max_retries: int = 5,
    timeout: float = 20.0,
) -> Dict[str, List[Tuple[str, str]]]:
    logger = logging.getLogger("chebi_roles")
    out: Dict[str, List[Tuple[str, str]]] = {}
    ids = [_normalize_chebi_id(cid) for cid in chebi_ids]

    total = len(ids)
    with requests.Session() as sess:
        for idx, cid in enumerate(ids, start=1):
            logger.info(f"[{idx}/{total}] Fetching roles for {cid}")
            url = _build_url(cid)

            try:
                resp = _get_with_retries(
                    sess, url, max_retries=max_retries, timeout=timeout
                )
                data = resp.json()
            except Exception:
                logger.exception(f"Failed to fetch roles for {cid}")
                out[cid] = []
                if per_request_pause:
                    time.sleep(per_request_pause)
                continue

            roles = set()
            for r in data.get("roles_classification", []) or []:
                if r.get("biological_role") is True:
                    rid = r.get("chebi_accession")
                    rlabel = r.get("name")
                    if rid and rlabel:
                        roles.add((rid, rlabel))

            out[cid] = sorted(roles, key=lambda x: (x[1].lower(), x[0]))
            logger.info(f"  -> {len(out[cid])} roles found")

            if per_request_pause:
                time.sleep(per_request_pause)

    return out


def main():
    ap = argparse.ArgumentParser(
        description="Extract ChEBI biological roles for all NPAtlas compounds."
    )
    ap.add_argument(
        "--workdir",
        required=True,
        help="Directory to store NPAtlas SDF, caches, logs, and outputs.",
    )
    ap.add_argument(
        "--sdf-url",
        default=NPATLAS_SDF_URL,
        help="Override NPAtlas SDF URL (defaults to official downloads).",
    )
    ap.add_argument(
        "--output",
        default="npatlas_chebi_roles.csv",
        help="Output CSV filename (inside workdir).",
    )
    ap.add_argument(
        "--cache-unichem",
        default="unichem_cache.jsonl",
        help="JSONL cache file (inside workdir) for UniChem responses.",
    )
    ap.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Console/file log level (default: INFO).",
    )
    args = ap.parse_args()
    workdir = Path(args.workdir).expanduser().resolve()
    log_path = setup_logging(workdir, args.log_level)

    logging.info("Starting NPAtlas → ChEBI roles extraction")
    logging.info("Workdir: %s", workdir)

    try:
        sdf_path = ensure_sdf(workdir, args.sdf_url)
    except Exception:
        logging.exception("Failed to ensure/download NPAtlas SDF")
        return

    # Load / resume UniChem cache
    cache_path = workdir / args.cache_unichem
    unichem_cache: Dict[str, Dict] = {}
    if cache_path.exists():
        try:
            with open(cache_path, "r", encoding="utf-8") as fh:
                for line in fh:
                    try:
                        obj = json.loads(line)
                        unichem_cache[obj["inchikey"]] = obj["record"]
                    except Exception:
                        continue
            logging.info("Loaded UniChem cache: %s (entries: %d)", cache_path, len(unichem_cache))
        except Exception:
            logging.exception("Failed to read UniChem cache: %s", cache_path)

    session = requests.Session()

    # First pass: parse SDF -> (npatlas_id, inchikey), then map to ChEBI
    rows_tmp: List[Tuple[str, str, List[str]]] = []
    seen_chebis: set = set()
    mols_total = 0
    mols_with_ik = 0
    mols_without_ik = 0
    cache_hits = 0
    cache_misses = 0
    chebi_linked = 0

    logging.info("Parsing SDF and querying UniChem ...")
    for idx, mol in tqdm(iter_npatlas_mols(sdf_path), desc="Reading NPAtlas SDF", unit="molecule"):
        mols_total += 1
        ik = mol_inchikey(mol)
        if not ik:
            mols_without_ik += 1
            continue
        mols_with_ik += 1
        npa_id = best_npatlas_id(mol, idx)

        if ik in unichem_cache:
            rec = unichem_cache[ik]
            cache_hits += 1
        else:
            rec = unichem_sources_for_inchikey(ik, session)
            unichem_cache[ik] = rec if rec is not None else {}
            try:
                with open(cache_path, "a", encoding="utf-8") as fh:
                    fh.write(json.dumps({"inchikey": ik, "record": rec}) + "\n")
            except Exception:
                logging.getLogger("cache").exception("Failed to append to cache file")
            cache_misses += 1
            time.sleep(UNICHEM_SLEEP)

        chebis = extract_chebi_ids_from_sources(rec) if rec else []
        if chebis:
            chebi_linked += 1
        rows_tmp.append((npa_id, ik, chebis))
        seen_chebis.update(chebis)

    logging.info("Parsed molecules: total=%d, with_inchikey=%d, without_inchikey=%d",
                 mols_total, mols_with_ik, mols_without_ik)
    logging.info("UniChem cache: hits=%d, misses=%d", cache_hits, cache_misses)
    logging.info("Compounds with at least one ChEBI mapping: %d", chebi_linked)
    logging.info("Unique ChEBI IDs to query for roles: %d", len(seen_chebis))

    # Second pass: SPARQL for roles
    chebi_list = sorted(seen_chebis)
    chebi_to_roles = get_biological_roles_for_chebis(chebi_list)

    # Write output
    out_path = workdir / args.output
    logging.info("Writing output CSV: %s", out_path)
    try:
        with open(out_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["npatlas_id", "inchikey", "chebi_id", "chebi_label", "role_id", "role_label"])
            row_count = 0
            for npa_id, ik, chebis in rows_tmp:
                if not chebis:
                    w.writerow([npa_id, ik, "", "", "", ""])
                    row_count += 1
                    continue
                for c in chebis:
                    roles = chebi_to_roles.get(c, [])
                    if not roles:
                        w.writerow([npa_id, ik, c, "", "", ""])
                        row_count += 1
                    else:
                        for role_id, role_label in roles:
                            w.writerow([npa_id, ik, c, "", role_id, role_label])
                            row_count += 1
        logging.info("Wrote %d rows to %s", row_count, out_path)
    except Exception:
        logging.exception("Failed to write output CSV: %s", out_path)
        return

    logging.info("Done. Output: %s", out_path)
    logging.info("Log file: %s", log_path)


if __name__ == "__main__":
    main()