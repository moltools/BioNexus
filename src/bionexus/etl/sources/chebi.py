"""Module for retrieving ChEBI annotations."""

import logging
import random
import time
from typing import Iterable, Optional

import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi
from sqlalchemy import select, delete
from sqlalchemy.dialects.postgresql import insert
from sqlalchemy.exc import SQLAlchemyError
from tqdm import tqdm

from bionexus.db.engine import SessionLocal
from bionexus.db.models import Annotation, Compound


RDLogger.DisableLog("rdApp.*")


logger = logging.getLogger(__name__)


BASE = "https://www.ebi.ac.uk/chebi/backend/api/public/compound"
UNICHEM_COMPOUND_API = "https://www.ebi.ac.uk/unichem/api/v1/compounds"

DEFAULT_HEADERS = {
    "User-Agent": "chebi-db-annotator/1.0 (+https://example.org)",
    "Accept": "application/json",
}

# Global in-memory caches for this process
_UNICHEM_CACHE: dict[str, Optional[dict]] = {}
_CHEBI_ROLES_CACHE: dict[str, list[tuple[str, str]]] = {}
_HTTP_SESSION: Optional[requests.Session] = None


def _get_session() -> requests.Session:
    global _HTTP_SESSION
    if _HTTP_SESSION is None:
        _HTTP_SESSION = requests.Session()
    return _HTTP_SESSION


def _normalize_chebi_id(cid: str) -> str:
    cid = cid.strip()
    return cid if cid.upper().startswith("CHEBI:") else f"CHEBI:{cid}"


def _build_chebi_url(chebi_id: str) -> str:
    from urllib.parse import quote

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
    log = logging.getLogger(__name__)
    attempt = 0
    while True:
        try:
            resp = session.get(url, headers=DEFAULT_HEADERS, timeout=timeout)
        except requests.RequestException:
            if attempt >= max_retries:
                log.exception("Network error while fetching %s", url)
                raise
            sleep = min(max_sleep, base_sleep * (2**attempt)) * (1 + random.random() * 0.2)
            log.warning("Retrying %s after network error, sleeping %.1fs", url, sleep)
            time.sleep(sleep)
            attempt += 1
            continue

        if resp.status_code in (429, 503):
            if attempt >= max_retries:
                log.error("Too many retries for %s (status %s)", url, resp.status_code)
                resp.raise_for_status()
            retry_after = resp.headers.get("Retry-After")
            if retry_after:
                try:
                    sleep = float(retry_after)
                except ValueError:
                    sleep = min(max_sleep, base_sleep * (2**attempt))
            else:
                sleep = min(max_sleep, base_sleep * (2**attempt))
            sleep *= 1 + random.random() * 0.2
            log.warning("Rate limited for %s, sleeping %.1fs", url, sleep)
            time.sleep(sleep)
            attempt += 1
            continue

        if 500 <= resp.status_code < 600:
            if attempt >= max_retries:
                log.error("Server error for %s (status %s)", url, resp.status_code)
                resp.raise_for_status()
            sleep = min(max_sleep, base_sleep * (2**attempt)) * (1 + random.random() * 0.2)
            log.warning("Server error %s for %s, retrying after %.1fs", resp.status_code, url, sleep)
            time.sleep(sleep)
            attempt += 1
            continue

        resp.raise_for_status()
        return resp
    

def mol_inchikey_from_smiles(smiles: str) -> Optional[str]:
    """SMILES -> RDKit Mol -> InChIKey, or None."""
    log = logging.getLogger("chebi_smiles")
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return inchi.MolToInchiKey(mol)
    except Exception:
        log.exception("Failed to get InChIKey for smiles=%r", smiles)
        return None
    

def unichem_sources_for_inchikey(inchikey: str) -> Optional[dict]:
    """
    Query UniChem 2.0 compound endpoint with type='inchikey'.
    Returns the first compound record (dict) or None.
    Uses in-memory cache to avoid repeated calls.
    """
    log = logging.getLogger("unichem")

    if inchikey in _UNICHEM_CACHE:
        return _UNICHEM_CACHE[inchikey]

    payload = {"type": "inchikey", "compound": inchikey}
    sess = _get_session()
    try:
        r = sess.post(UNICHEM_COMPOUND_API, json=payload, timeout=30)
        r.raise_for_status()
        data = r.json()
        comps = data.get("compounds") or []
        if not comps:
            log.debug("No UniChem record for %s", inchikey)
            _UNICHEM_CACHE[inchikey] = None
            return None
        for c in comps:
            if c.get("standardInchiKey") == inchikey:
                _UNICHEM_CACHE[inchikey] = c
                return c
        _UNICHEM_CACHE[inchikey] = comps[0]
        return comps[0]
    except Exception:
        log.exception("UniChem request failed for %s", inchikey)
        _UNICHEM_CACHE[inchikey] = None
        return None


def extract_chebi_ids_from_sources(compound_record: dict) -> list[str]:
    """
    From UniChem 'sources' list, find ChEBI entries.
    Prefer matching by source name (robust to sourceID changes).
    """
    out: list[str] = []
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


def get_biological_roles_for_chebis(
    chebi_ids: Iterable[str],
    *,
    per_request_pause: float = 0.2,
    max_retries: int = 5,
    timeout: float = 20.0,
) -> dict[str, list[tuple[str, str]]]:
    """
    Fetch biological roles for a list of ChEBI IDs, using global role cache.
    Returns mapping CHEBI_ID -> list[(role_id, role_label)].
    """
    log = logging.getLogger("chebi_roles")
    sess = _get_session()
    ids = [_normalize_chebi_id(cid) for cid in chebi_ids]

    # Only query new ones
    to_fetch = [cid for cid in ids if cid not in _CHEBI_ROLES_CACHE]
    total = len(to_fetch)

    for idx, cid in enumerate(to_fetch, start=1):
        log.info("[ChEBI roles %d/%d] Fetching roles for %s", idx, total, cid)
        url = _build_chebi_url(cid)
        try:
            resp = _get_with_retries(
                sess, url, max_retries=max_retries, timeout=timeout
            )
            data = resp.json()
        except Exception:
            log.exception("Failed to fetch roles for %s", cid)
            _CHEBI_ROLES_CACHE[cid] = []
            if per_request_pause:
                time.sleep(per_request_pause)
            continue

        roles = set()
        for r in data.get("roles_classification", []) or []:
            if (r.get("biological_role") is True) or (r.get("chemical_role") is True):
                rid = r.get("chebi_accession")
                rlabel = r.get("name")
                if rid and rlabel:
                    roles.add((rid, rlabel))

        _CHEBI_ROLES_CACHE[cid] = sorted(roles, key=lambda x: (x[1].lower(), x[0]))
        log.info("  -> %d roles found", len(_CHEBI_ROLES_CACHE[cid]))

        if per_request_pause:
            time.sleep(per_request_pause)

    # Build result mapping for all ids (including those previously cached)
    return {cid: _CHEBI_ROLES_CACHE.get(cid, []) for cid in ids}


def fetch_chebi_annotations(
    compound_pairs: list[tuple[int, str]],
) -> dict[int, dict]:
    """
    Given a list of (compound_id, smiles) pairs, fetch ChEBI annotations.

    :param compound_pairs: list of (compound_id, smiles) tuples
    :return: mapping compound_id -> json_blob with ChEBI annotations
    """
    log = logging.getLogger("fetch_chebi_annotations")
    results: dict[int, dict] = {}

    # First pass: SMILES -> inchikey -> unichem -> chebi ids
    comp_to_smiles: dict[int, str] = {}
    comp_to_inchikey: dict[int, Optional[str]] = {}
    comp_to_chebis: dict[int, list[str]] = {}
    seen_chebis: set[str] = set()

    for comp_id, smiles in tqdm(compound_pairs, leave=False):
        smiles = (smiles or "").strip()
        if not smiles:
            continue

        comp_to_smiles[comp_id] = smiles
        ik = mol_inchikey_from_smiles(smiles)
        comp_to_inchikey[comp_id] = ik

        if not ik:
            continue

        rec = unichem_sources_for_inchikey(ik)
        chebis = extract_chebi_ids_from_sources(rec) if rec else []
        comp_to_chebis[comp_id] = chebis
        seen_chebis.update(chebis)

    if not seen_chebis:
        log.info("No ChEBI IDs found for any compounds in this chunk.")
        for comp_id, smiles in comp_to_smiles.items():
            results[comp_id] = {
                "smiles": smiles,
                "inchikey": comp_to_inchikey.get(comp_id),
                "chebi_ids": [],
            }
        return results
    
    log.info(f"Fetching ChEBI annotations for {len(seen_chebis)} unique ChEBI IDs.")

    # Second pass: chebi -> roles (with cache)
    chebi_roles = get_biological_roles_for_chebis(sorted(seen_chebis))

    # Build final json_blob per compound
    for comp_id, smiles in comp_to_smiles.items():
        ik = comp_to_inchikey.get(comp_id)
        chebis = comp_to_chebis.get(comp_id) or []

        chebi_entries = []
        for cid in chebis:
            roles = chebi_roles.get(_normalize_chebi_id(cid)) or []
            chebi_entries.append(
                {
                    "chebi_id": _normalize_chebi_id(cid),
                    "roles": [
                        {"role_id": rid, "role_label": rlabel}
                        for rid, rlabel in roles
                    ],
                }
            )

        results[comp_id] = {
            "smiles": smiles,
            "inchikey": ik,
            "chebi_entries": chebi_entries,
        }

    return results


def _parse_annotations(json_blob: dict) -> list[tuple[str, str, str]]:
    """
    From ChEBI JSON blob, extract annotation tuples (scheme, key, value).

    :param json_blob: ChEBI JSON blob
    :return: list of (scheme, key, value) tuples
    """
    out: set[tuple[str, str, str]] = set()

    for chebi_entry in json_blob.get("chebi_entries", []) or []:
        cid = chebi_entry.get("chebi_id")
        for role in chebi_entry.get("roles", []) or []:
            rid = role.get("role_id")
            rlabel = role.get("role_label")
            if rid and rlabel:
                out.add(("chebi", "role", str(rlabel)))

    return sorted(out)


def annotate_with_chebi(
    recompute: bool,
    *,
    chunk_size: int =  5000,
) -> int:
    """
    Annotate compounds with ChEBI annotations.

    :param recompute: if True, recompute annotations for all compounds
    :param chunk_size: number of compounds to process in each chunk
    :return: total number of annotations inserted
    """
    total_inserted = 0

    with SessionLocal() as s:
        if not recompute:
            # Subquery: does an chebi annotation already exist for this compound?
            chebi_exists = (
                select(1)
                .where((Annotation.compound_id == Compound.id) & (Annotation.scheme == "chebi"))
                .exists()
            )
            # IDs that DO NOT have chebi annotations yet
            comp_ids = s.scalars(select(Compound.id).where(~chebi_exists)).all()
            logger.info(f"Found {len(comp_ids)} compounds without ChEBI annotations.")
        else:
            comp_ids = s.scalars(select(Compound.id)).all()
            logger.info(f"Recomputing ChEBI annotations for {len(comp_ids)} compounds.")

        if not comp_ids:
            logger.info("No compounds to annotate. Exiting.")
            return 0

    # PRocess in manageable chunks
    with tqdm(total=len(comp_ids), desc="Annotatin compounds", unit="cmpd") as pbar:
        for start in range(0, len(comp_ids), chunk_size):
            end = start + chunk_size
            chunk_ids = comp_ids[start:end]

            # Fetch compound structures for this chunk (id and smiles)
            with SessionLocal() as s:
                comps: list[Compound] = s.scalars(select(Compound).where(Compound.id.in_(chunk_ids))).all()

                # Delete existing annotations if recomputing
                if recompute and comps:
                    s.execute(
                        delete(Annotation).where(
                            Annotation.compound_id.in_([c.id for c in comps]),
                            Annotation.scheme == "chebi",
                        )
                    )
                    s.commit()

                pairs: list[tuple[int, str]] = [(c.id, c.smiles) for c in comps if c.smiles]

            if not pairs:
                pbar.update(len(chunk_ids))  # all missing SMILES
                continue

            # Fetch annotations from ChEBI for this chunk
            results = fetch_chebi_annotations(pairs)

            # Build bulk rows and insert
            rows_to_insert = []
            for comp_id, json_blob in results.items():
                for scheme, key, value in _parse_annotations(json_blob):
                    rows_to_insert.append(
                        dict(
                            compound_id=comp_id,
                            scheme=scheme,
                            key=key,
                            value=value,
                            metadata_json=None,
                        )
                    )

            if not rows_to_insert:
                pbar.update(len(chunk_ids))
                continue

            try:
                with SessionLocal() as s:
                    stmt = (
                        insert(Annotation)
                        .values(rows_to_insert)
                        .on_conflict_do_nothing(constraint="uq_annotation_target_scheme_key_value")
                    )
                    res = s.execute(stmt)
                    s.commit()
                    inserted = res.rowcount if (res.rowcount or 0) > 0 else 0
                    if inserted == 0:
                        inserted = len(rows_to_insert)
                    total_inserted += inserted
                    logger.info(f"Inserted {inserted} ChEBI annotations for chunk {start}-{end}.")
            except SQLAlchemyError as e:
                logger.error(f"Failed to insert ChEBI annotations for chunk {start}-{end}: {e}")

    return total_inserted
