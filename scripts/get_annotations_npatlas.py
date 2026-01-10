#!/usr/bin/env python3

"""Script for retrieving ChEBI annotations for compounds; consumes a CSV/TSV file containing a SMILES column."""

import argparse
import asyncio
import csv
import json
import logging
import random
import time
from dataclasses import dataclass
from typing import Any

import aiohttp
import requests
from tqdm import tqdm
from rdkit import RDLogger

from retromol.chem.mol import smiles_to_mol, mol_to_inchikey

from bionexus.utils.logging import setup_logging


# Disable RDKit warnings
RDLogger.DisableLog('rdApp.*')


log = logging.getLogger(__name__)


# Configuration
TESTING = False  # limit to 100 entries for testing
CHUNK_SIZE = 1000
NP_CONCURRENCY = 8
NP_RATE = 10.0  # requests per second
NP_TIMEOUT = 10
NP_RETRIES = 3
CHEBI_PAUSE = 0.2 # seconds between ChEBI requests
ENRICH_NP = True
ENRICH_CHEBI = True


# NPClassifier API settings
NP_API_BASES = [
    "https://npclassifier.ucsd.edu/classify",
    "https://npclassifier.gnps2.org/classify",
]

NP_HEADERS = {
    "Accept": "application/json",
    "User-Agent": "BioNexus/0.2 (NPClassifier bulk annotate)",
}

# ChEBI via UniChem
CHEBI_BASE = "https://www.ebi.ac.uk/chebi/backend/api/public/compound"
UNICHEM_COMPOUND_API = "https://www.ebi.ac.uk/unichem/api/v1/compounds"

CHEBI_HEADERS = {
    "User-Agent": "bionexus-chebi-annotator/1.0",
    "Accept": "application/json",
}

_HTTP: requests.Session | None = None
_UNICHEM_CACHE: dict[str, dict | None] = {}
_CHEBI_ROLE_CACHE: dict[str, list[tuple[str, str, str]]] = {}  # CHEBI_ID -> [(role_kind, role_id, role_label)]


def cli() -> argparse.Namespace:
    """
    Command line interface for script.
    
    :return: parsed command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--npatlas", type=str, required=True, help="path to NPAtlas database in JSON format")
    parser.add_argument("--output", type=str, required=True, help="output annotations file path")
    return parser.parse_args()


def _safe_str(x: Any) -> str:
    """
    Safely convert an object to string, returning an empty string if the object is None.

    :param x: object to convert
    :return: string representation of the object or empty string
    """
    return "" if x is None else str(x)


def _taxonomy_species(genus: str | None, species_epithet: str | None) -> str | None:
    """
    Construct species name from genus and species epithet.

    :param genus: genus name
    :param species_epithet: species epithet
    :return: full species name or None
    """
    genus = (genus or "").strip() or None
    epithet = (species_epithet or "").strip() or None

    if genus and epithet:
        return f"{genus} {epithet}"
    
    return genus or epithet


# -------------------------
# NPClassifier fetching
# -------------------------

def parse_npclassifier_rows(blob: dict[str, Any] | None) -> list[tuple[str, str, str]]:
    """
    Parse NPClassifier annotation rows from a given blob.

    :param blob: NPClassifier annotation blob
    :return: list of annotation tuples (scheme, key, value)
    """
    out: list[tuple[str, str, str]] = []
    if not blob or not isinstance(blob, dict):
        return out
    
    isgly = blob.get("isglycoside", blob.get("isgly", blob.get("is_glycoside")))
    if isgly is not None:
        out.append(("NPClassifier", "is_glycoside", str(bool(isgly)).lower()))
    
    for key_json, key_out in [
        ("class_results", "class"),
        ("superclass_results", "superclass"),
        ("pathway_results", "pathway"),
    ]:
        for v in (blob.get(key_json) or []):
            out.append(("NPClassifier", key_out, _safe_str(v)))

    return out


async def _np_fetch_one(
    session: aiohttp.ClientSession,
    smiles: str,
    *,
    timeout_s: int,
    retries: int,
    delay_between: float | None,
    semaphore: asyncio.Semaphore,
) -> dict[str, Any] | None:
    """
    Fetch NPClassifier annotation for a single SMILES string.

    :param session: aiohttp client session
    :param smiles: SMILES string to query
    :param timeout_s: request timeout in seconds
    :param retries: number of retries on failure
    :param delay_between: delay between requests in seconds
    :param semaphore: asyncio semaphore for limiting concurrency
    :return: NPClassifier annotation blob or None
    """
    async with semaphore:
        if delay_between:
            await asyncio.sleep(delay_between)

        backoff = 0.5
        last_err: Exception | None = None

        for attempt in range(retries + 1):
            for base in NP_API_BASES:
                try:
                    async with session.get(
                        base,
                        params={"smiles": smiles},
                        timeout=timeout_s,
                        headers=NP_HEADERS,
                    ) as resp:
                        status = resp.status
                        text = await resp.text()
                        if status == 200:
                            # Header sometimes wrong; try parse anyway
                            try:
                                return json.loads(text)
                            except json.JSONDecodeError:
                                pass

                        if status == 429:
                            ra = resp.headers.get("Retry-After")
                            wait = float(ra) if (ra and ra.isdigit()) else backoff * (1 + random.random())
                            await asyncio.sleep(wait)
                            continue

                        if 500 <= status < 600:
                            last_err = RuntimeError(f"{base} returned {status}")
                            continue

                        if 400 <= status < 500:
                            log.warning("NPClassifier 4xx (%s) for SMILES=%r", status, smiles)
                            return None
                        
                except (aiohttp.ClientError, asyncio.TimeoutError) as e:
                    # Log and retry
                    last_err = e

            if attempt < retries:
                # Exponential backoff with jitter
                await asyncio.sleep(backoff * (1 + random.random()))
                backoff *= 2

        if last_err:
            log.error("NPClassifier giving up for SMILES=%r :: %s", smiles, last_err)
        
        return None


async def np_fetch_many(
    smiles_list: list[str],
    *,
    concurrency: int,
    rate_per_sec: float | None,
    timeout_s: int,
    retries: int,
    pbar: tqdm | None = None,
) -> dict[str, dict[str, Any] | None]:
    """
    Fetch NPClassifier annotations for multiple SMILES strings concurrently.

    :param smiles_list: list of SMILES strings to query
    :param concurrency: maximum number of concurrent requests
    :param rate_per_sec: maximum request rate per second
    :param timeout_s: request timeout in seconds
    :param retries: number of retries on failure
    :param pbar: optional tqdm progress bar
    :return: dictionary mapping SMILES strings to their NPClassifier annotation blobs or None
    """
    semaphore = asyncio.Semaphore(max(1, concurrency))
    delay_between = (1.0 / rate_per_sec) if rate_per_sec and rate_per_sec > 0 else None

    connector = aiohttp.TCPConnector(
        limit=max(4, concurrency),
        limit_per_host=max(2, concurrency // 2),
        enable_cleanup_closed=True,
    )
    timeout = aiohttp.ClientTimeout(total=None)

    out: dict[str, dict[str, Any] | None] = {}
    async with aiohttp.ClientSession(timeout=timeout, connector=connector, headers=NP_HEADERS) as session:
        async def _run(smi: str):
            out[smi] = await _np_fetch_one(
                session, smi,
                timeout_s=timeout_s,
                retries=retries,
                delay_between=delay_between,
                semaphore=semaphore,
            )
            if pbar is not None:
                pbar.update(1)
        tasks = [asyncio.create_task(_run(smi)) for smi in smiles_list]
        await asyncio.gather(*tasks)
    return out


# -------------------------
# ChEBI fetching
# -------------------------

def _http() -> requests.Session:
    """
    Get or create a global HTTP session.
    
    :return: requests session
    """
    global _HTTP
    if _HTTP is None:
        _HTTP = requests.Session()
    return _HTTP


def _normalize_chebi_id(cid: str) -> str:
    """
    Normalize a ChEBI ID to ensure it starts with "CHEBI:".

    :param cid: ChEBI ID
    :return: normalized ChEBI ID
    """
    cid = (cid or "").strip()
    if not cid:
        return cid

    return cid if cid.upper().startswith("CHEBI:") else f"CHEBI:{cid}"


def _chebi_url(chebi_id: str) -> str:
    """
    Construct ChEBI API URL for a given ChEBI ID.

    :param chebi_id: ChEBI ID
    :return: ChEBI API URL
    """
    from urllib.parse import quote
    enc = quote(chebi_id, safe="")
    return f"{CHEBI_BASE}/{enc}/?only_ontology_parents=false&only_ontology_children=false"


def _get_with_retries(url: str, *, max_retries: int = 5, timeout: float = 20.0) -> requests.Response:
    """
    Perform an HTTP GET request with retries and exponential backoff.

    :param url: URL to fetch
    :param max_retries: maximum number of retries
    :param timeout: request timeout in seconds
    :return: HTTP response
    """
    sess = _http()
    base_sleep = 0.8
    max_sleep = 8.0
    attempt = 0

    while True:
        try:
            resp = sess.get(url, headers=CHEBI_HEADERS, timeout=timeout)
        except requests.RequestException:
            if attempt >= max_retries:
                raise
            time.sleep(min(max_sleep, base_sleep * (2**attempt)) * (1 + random.random() * 0.2))
            attempt += 1
            continue

        if resp.status_code in (429, 503) or (500 <= resp.status_code < 600):
            if attempt >= max_retries:
                resp.raise_for_status()
            ra = resp.headers.get("Retry-After")
            sleep = float(ra) if (ra and ra.isdigit()) else min(max_sleep, base_sleep * (2**attempt))
            time.sleep(sleep * (1 + random.random() * 0.2))
            attempt += 1
            continue

        resp.raise_for_status()
        return resp
    

def unichem_record_for_inchikey(inchikey: str) -> dict | None:
    """
    Retrieve UniChem compound record for a given InChIKey.

    :param inchikey: InChIKey of the compound
    :return: UniChem compound record or None
    """
    if inchikey in _UNICHEM_CACHE:
        return _UNICHEM_CACHE[inchikey]
    
    payload = {"type": "inchikey", "compound": inchikey}

    try:
        r = _http().post(UNICHEM_COMPOUND_API, json=payload, timeout=30)
        r.raise_for_status()
        data = r.json()
        comps = data.get("compounds", [])
        rec = None
        for c in comps:
            if c.get("standardInchiKey") == inchikey:
                rec = c
                break
        rec = rec or (comps[0] if comps else None)
        _UNICHEM_CACHE[inchikey] = rec
        return rec
    except Exception:
        log.exception(f"UniChem fetch failed for InChIKey={inchikey}")
        _UNICHEM_CACHE[inchikey] = None
        return None
    

def chebi_ids_from_unichem(rec: dict | None) -> list[str]:
    """
    Extract ChEBI IDs from a UniChem compound record.

    :param rec: UniChem compound record
    :return: list of ChEBI IDs
    """
    if not rec:
        return []
    out: set[str] = set()
    for s in (rec.get("sources") or []):
        if (s.get("shortName") or "").lower() == "chebi":
            cid = _normalize_chebi_id(str(s.get("compoundId") or "").strip())
            if cid:
                out.add(cid)

    return sorted(out)


def fetch_roles_for_chebi_ids(chebi_ids: list[str], *, pause_s: float = 0.2) -> None:
    """
    Fetch and cache roles for a list of ChEBI IDs.

    :param chebi_ids: list of ChEBI IDs
    :param pause_s: pause in seconds between requests
    """
    # Only fetch uncached
    to_fetch = [c for c in chebi_ids if c not in _CHEBI_ROLE_CACHE]

    for cid in tqdm(to_fetch, desc="ChEBI roles", unit="chebi", leave=False):
        try:
            data = _get_with_retries(_chebi_url(cid)).json()

            roles: set[tuple[str, str, str]] = set()
            for r in (data.get("roles_classification") or []):
                rid = r.get("chebi_accession")
                rlabel = r.get("name")
                if not (rid and rlabel):
                    continue

                if r.get("biological_role") is True:
                    roles.add(("biological_role", str(rid), str(rlabel)))
                if r.get("chemical_role") is True:
                    roles.add(("chemical_role", str(rid), str(rlabel)))
                    
            _CHEBI_ROLE_CACHE[cid] = sorted(roles, key=lambda x: (x[0], x[2].lower(), x[1]))
        except Exception:
            log.exception(f"ChEBI fetch failed for ChEBI ID={cid}")
            _CHEBI_ROLE_CACHE[cid] = []
        
        if pause_s:
            time.sleep(pause_s)


def chebi_rows_for_inchikey(inchikey: str) -> list[tuple[str, str, str]]:
    """
    Retrieve ChEBI annotation rows for a given InChIKey.

    :param inchikey: InChIKey of the compound
    :return: list of annotation tuples (scheme, key, value)
    """
    rec = unichem_record_for_inchikey(inchikey)
    chebis = chebi_ids_from_unichem(rec)
    
    rows: list[tuple[str, str, str]] = []    
    for cid in chebis:
        for role_kind, _rid, rlabel in _CHEBI_ROLE_CACHE.get(cid, []):
            rows.append(("ChEBI", role_kind, rlabel))
        
    return rows


# -------------------------
# Main
# -------------------------

@dataclass
class EntryMini:
    """
    Minimal NPAtlas entry representation.

    :var inchikey: InChIKey of the compound
    :var smiles: SMILES string of the compound
    :var npaid: NPAtlas identifier
    :var name: original name of the compound
    :var synonyms: list of synonyms
    :var taxonomy: taxonomy information
    :var npclassifier: NPClassifier annotation blob
    """

    inchikey: str
    smiles: str
    npaid: str
    name: str | None
    synonyms: list[str]
    taxonomy: dict[str, Any]
    npclassifier: dict[str, Any] | None


def main() -> None:
    """
    Main function to retrieve ChEBI annotations.
    """
    args = cli()
    setup_logging(level=logging.INFO)

    delimiter = "\t" if args.output.endswith(".tsv") else ","
    database = "NPAtlas"
    annotation_target = "compound"

    with open(args.npatlas, "r") as infile:
        raw = json.load(infile)

    # First pass: parse, compoute inchikey, collect enrichment needs
    entries: list[EntryMini] = []
    need_np_smiles: set[str] = set()
    need_chebi_inchikeys: set[str] = set()

    for entry in tqdm(raw, desc="Parsing NPAtlas", unit="entry"):
        npaid = entry.get("npaid")
        smiles = entry.get("smiles")

        if not npaid or not smiles:
            log.warning("skipping entry without NPAID or SMILES")
            continue

        try:
            ik = mol_to_inchikey(smiles_to_mol(smiles))
        except Exception:
            log.warning(f"skipping entry with invalid SMILES: {smiles}")
            continue

        npblob = entry.get("npclassifier") or None
        if not npblob or not npblob.get("class_results") or not npblob.get("superclass_results") or not npblob.get("pathway_results"):
            need_np_smiles.add(smiles)

        need_chebi_inchikeys.add(ik)
            
        synonyms = []
        for syn in (entry.get("synonyms") or []):
            nm = syn.get("name") if isinstance(syn, dict) else None
            if nm:
                synonyms.append(nm)

        entries.append(EntryMini(
            inchikey=ik,
            smiles=smiles,
            npaid=str(npaid),
            name=entry.get("original_name"),
            synonyms=synonyms,
            taxonomy=entry.get("origin_organism") or {},
            npclassifier=npblob,
        ))

    # LIMIT TO 1000 FOR TESTING
    if TESTING:
        entries = entries[:100]
        need_np_smiles = need_np_smiles.intersection({e.smiles for e in entries})
        need_chebi_inchikeys = need_chebi_inchikeys.intersection({e.inchikey for e in entries})

    log.info(f"parsed {len(entries)} NPAtlas entries")
    log.info(f"need NPClassifier annotations for {len(need_np_smiles)} unique SMILES")
    log.info(f"need ChEBI annotations for {len(need_chebi_inchikeys)} unique InChIKeys")

    # Enrich NPClassifier annotations
    np_enriched: dict[str, dict[str, Any] | None] = {}
    if ENRICH_NP and need_np_smiles:
        log.info("fetching NPClassifier annotations...")
        
        # Chunk to avoid creating too many tasks at once
        need_list = list(need_np_smiles)

        with tqdm(total=len(need_list), desc="NPClassifier", unit="smiles") as p_np:
            for i in range(0, len(need_list), CHUNK_SIZE):
                chunk = need_list[i : i + CHUNK_SIZE]
                res = asyncio.run(
                    np_fetch_many(
                        chunk,
                        concurrency=NP_CONCURRENCY,
                        rate_per_sec=NP_RATE,
                        timeout_s=NP_TIMEOUT,
                        retries=NP_RETRIES,
                        pbar=p_np,
                    )
                )
                np_enriched.update(res)

    # Enrich ChEBI roles (inchikey -> chebi ids -> roles)
    if ENRICH_CHEBI and need_chebi_inchikeys:
        log.info("fetching ChEBI annotations...")
        
        # inchikey -> chebi ids (cached)
        all_chebis: set[str] = set()

        for ik in tqdm(list(need_chebi_inchikeys), desc="UniChem", unit="inchikey"):
            rec = unichem_record_for_inchikey(ik)
            all_chebis.update(chebi_ids_from_unichem(rec))

        log.info(f"found {len(all_chebis)} unique ChEBI IDs in UniChem for the requested InChIKeys")

        if all_chebis:
            fetch_roles_for_chebi_ids(sorted(all_chebis), pause_s=CHEBI_PAUSE)
        else:
            log.info("no ChEBI IDs found in UniChem for the requested InChIKeys")


    # Second pass: write output
    with open(args.output, "w", newline="") as out:
        w = csv.writer(out, delimiter=delimiter, quoting=csv.QUOTE_MINIMAL)
        w.writerow(["table", "type", "inchikey", "scheme", "key", "value"])

        for e in tqdm(entries, desc="writing", unit="entry"):
            ik = e.inchikey

            # References (NPAtlas)
            if e.name:
                w.writerow(["reference", annotation_target, ik, database, e.npaid, e.name])
            for syn in e.synonyms:
                w.writerow(["reference", annotation_target, ik, database, e.npaid, syn])

            # Taxonomy annotations
            tx = e.taxonomy or {}
            if tx:
                domain = tx.get("type")
                genus = tx.get("genus")
                species = _taxonomy_species(genus, tx.get("species"))

                if domain is not None:
                    w.writerow(["annotation", annotation_target, ik, "taxonomy", "domain", _safe_str(domain)])
                if genus is not None:
                    w.writerow(["annotation", annotation_target, ik, "taxonomy", "genus", _safe_str(genus)])
                if species is not None:
                    w.writerow(["annotation", annotation_target, ik, "taxonomy", "species", _safe_str(species)])

            # NPClassifier annotations (prefer existing, else enriched)
            if ENRICH_NP:
                npblob = e.npclassifier
                if not npblob or not npblob.get("class_results") or not npblob.get("superclass_results") or not npblob.get("pathway_results"):
                    npblob = np_enriched.get(e.smiles) or npblob

                for scheme, key, value in parse_npclassifier_rows(npblob):
                    w.writerow(["annotation", annotation_target, ik, scheme, key, value])

            # ChEBI annotations
            if ENRICH_CHEBI:
                for scheme, key, value in chebi_rows_for_inchikey(ik):
                    w.writerow(["annotation", annotation_target, ik, scheme, key, value])


if __name__ == "__main__":
    main()
