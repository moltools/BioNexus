"""NPClassifier annotation source for BioNexus ETL."""

from __future__ import annotations
import logging
import json
import random
from typing import Any, Optional

import asyncio
import aiohttp
from tqdm import tqdm
from sqlalchemy import select, delete
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.dialects.postgresql import insert

from bionexus.db.models import Compound, Annotation
from bionexus.db.engine import SessionLocal


logger = logging.getLogger(__name__)


API_BASES = [
    "https://npclassifier.ucsd.edu/classify",  # primary (usually more stable)
    "https://npclassifier.gnps2.org/classify",  # fallback
]


DEFAULT_HEADERS = {
    "Accept": "application/json",
    "User-Agent": "BioNexus/0.2 (NPClassifier bulk annotate; contact: david.meijer@wur.nl)",
}


def _parse_annotations(annotations: dict[str, Any]) -> list[tuple[str, str, str]]:
    """
    Parse NPClassifier JSON annotations into list of (scheme, key, value) tuples.

    :param annotations: JSON dict from NPClassifier
    :return: list of (scheme, key, value) tuples
    """
    out: list[tuple[str, str, str]] = []

    if not annotations or not isinstance(annotations, dict):
        return out

    # Booleans to lowercase string
    isgly = annotations.get("isglycoside", None)
    if isgly is not None:
        out.append(("npclassifier", "isglycoside", str(bool(isgly)).lower()))

    for key_json, key_db in [
        ("class_results", "class"),
        ("pathway_results", "pathway"),
        ("superclass_results", "superclass"),
    ]:
        values = annotations.get(key_json) or []
        for v in values:
            # Be defensive about v being non-string; coerce to str
            out.append(("npclassifier", key_db, str(v)))

    return out


async def _fetch_one(
    session: aiohttp.ClientSession,
    smiles: str,
    *,
    timeout_s: int = 10,
    retries: int = 3,
    delay_between: Optional[float] = None,
    semaphore: Optional[asyncio.Semaphore] = None,
) -> Optional[dict[str, Any]]:
    """
    Fetch NPClassifier result for a single SMILES string with retries.

    :param session: aiohttp ClientSession
    :param smiles: SMILES string to classify
    :param timeout_s: per-request timeout in seconds
    :param retries: number of retries on failure
    :param delay_between: optional delay before request (for rate limiting)
    :param semaphore: optional asyncio.Semaphore for concurrency limiting
    :return: JSON dict from NPClassifier or None on failure
    """
    guard = semaphore or asyncio.Semaphore(1)
    async with guard:
        if delay_between:
            await asyncio.sleep(delay_between)

        last_err: Optional[Exception] = None
        backoff = 0.5

        for attempt in range(retries + 1):
            for base in API_BASES:
                try:
                    async with session.get(
                        base,
                        params={"smiles": smiles},
                        timeout=timeout_s,
                        headers=DEFAULT_HEADERS,
                    ) as resp:
                        status = resp.status
                        ctype = resp.headers.get("Content-Type", "")
                        text = await resp.text()  # read once

                        # Handle successful JSON the normal way
                        if status == 200:
                            # If server mislabeled JSON as text/html, try manual parse
                            if "application/json" in ctype.lower():
                                return json.loads(text)
                            # Sometimes the body is JSON but header is wrong
                            if text and text.lstrip().startswith(("{", "[")):
                                try:
                                    return json.loads(text)
                                except json.JSONDecodeError:
                                    pass  # fall through to retry

                            # Got HTML (likely rate-limit or WAF). Retry.
                            logger.debug(
                                f"Non-JSON 200 from {base}: {ctype}; first 120 chars: {text[:120]!r}"
                            )

                        # Explicit backoff cases
                        if status == 429:
                            # Respect Retry-After if present
                            ra = resp.headers.get("Retry-After")
                            wait = (
                                float(ra)
                                if ra and ra.isdigit()
                                else backoff * (1 + random.random())
                            )
                            await asyncio.sleep(wait)
                            continue  # try again (same base)

                        if 500 <= status < 600:
                            # Transient server error
                            last_err = RuntimeError(f"{base} returned {status}")
                            continue  # try fallback base or retry

                        if 400 <= status < 500:
                            # Bad request (not retryable except maybe 429 handled above)
                            logger.warning(
                                f"4xx for SMILES (no retry): {status} :: {smiles}"
                            )
                            return None

                except (aiohttp.ClientError, asyncio.TimeoutError) as e:
                    last_err = e
                    # Try next base or fall through to retry loop

            # Retry across bases with jittered backoff
            if attempt < retries:
                sleep_s = backoff * (1 + random.random())
                await asyncio.sleep(sleep_s)
                backoff *= 2
            else:
                if last_err:
                    logger.error(
                        f"Giving up on SMILES after retries: {smiles} :: {last_err}"
                    )
                return None


async def _fetch_many(
    items: list[tuple[int, str]],  # (compound_id, smiles)
    *,
    concurrency: int = 16,
    rate_per_sec: Optional[float] = 10.0,
    timeout_s: int = 10,
    retries: int = 3,
    pbar: Optional[tqdm] = None,
) -> dict[int, Optional[dict[str, Any]]]:
    """
    Concurrently fetch NPClassifier results for multiple SMILES strings.

    :param items: list of (compound_id, smiles) tuples
    :param concurrency: max concurrent requests
    :param rate_per_sec: optional rate limit (requests per second)
    :param timeout_s: per-request timeout in seconds
    :param retries: number of retries per request
    :param pbar: optional tqdm progress bar to update
    :return: dict mapping compound_id to JSON dict or None on failure
    """
    # Bound concurrency
    semaphore = asyncio.Semaphore(max(1, concurrency))
    # Optional simple throttle: delay each request roughly 1/rate_per_sec
    delay_between = (1.0 / rate_per_sec) if rate_per_sec and rate_per_sec > 0 else None

    timeout = aiohttp.ClientTimeout(total=None)  # per-request handled in _fetch_one
    connector = aiohttp.TCPConnector(
        limit=max(4, concurrency),
        limit_per_host=max(2, concurrency // 2),
        enable_cleanup_closed=True,
    )

    results: dict[int, Optional[dict[str, Any]]] = {}

    async with aiohttp.ClientSession(
        timeout=timeout, connector=connector, headers=DEFAULT_HEADERS
    ) as sess:

        async def run_one(comp_id: int, smi: str):
            res = await _fetch_one(
                sess,
                smi,
                timeout_s=timeout_s,
                retries=retries,
                delay_between=delay_between,
                semaphore=semaphore,
            )
            results[comp_id] = res
            if pbar:
                pbar.update(1)

        tasks = [asyncio.create_task(run_one(cid, smi)) for cid, smi in items]
        # Gather and propagate any exceptions
        await asyncio.gather(*tasks)

    return results


def annotate_with_npclassifier(
    recompute: bool,
    *,
    chunk_size: int = 5000,
    concurrency: int = 8,
    rate_per_sec: Optional[float] = 10.0,
    timeout_s: int = 10,
    retries: int = 3,
) -> int:
    """
    Annotate compounds in the database with NPClassifier annotations.

    :param recompute: if True, recompute annotations for all compounds; if False, only for those missing
    :param chunk_size: number of compounds to process per chunk
    :param concurrency: max concurrent requests to NPClassifier API
    :param rate_per_sec: optional rate limit (requests per second)
    :param timeout_s: per-request timeout in seconds
    :param retries: number of retries per request
    :return: total number of new annotations inserted
    """
    total_inserted = 0

    # 1) Collect target compound ids (and SMILES)
    with SessionLocal() as s:
        if not recompute:
            # Subquery: does an npclassifier annotation exist for this compound?
            npclass_exists = (
                select(1)
                .where(
                    (Annotation.compound_id == Compound.id)
                    & (Annotation.scheme == "npclassifier")
                )
                .exists()
            )
            # IDs that DO NOT have npclassifier yet
            comp_ids = s.scalars(select(Compound.id).where(~npclass_exists)).all()
            logger.info(
                f"Found {len(comp_ids)} compounds without NPClassifier annotations"
            )
        else:
            comp_ids = s.scalars(select(Compound.id)).all()
            logger.info(
                f"Recomputing NPClassifier annotations for all {len(comp_ids)} compounds"
            )

    if not comp_ids:
        logger.info("No compounds to annotate.")
        return 0

    # 2) Process in manageable chunks
    with tqdm(total=len(comp_ids), desc="Annotating compounds", unit="cmpd") as pbar:
        for start in range(0, len(comp_ids), chunk_size):
            end = start + chunk_size
            chunk_ids = comp_ids[start:end]

            # Fetch the compounds for this chunk (id + smiles)
            with SessionLocal() as s:
                comps: list[Compound] = s.scalars(
                    select(Compound).where(Compound.id.in_(chunk_ids))
                ).all()

                # Optionally delete existing annotations if recompute
                if recompute and comps:
                    s.execute(
                        delete(Annotation).where(
                            Annotation.compound_id.in_([c.id for c in comps]),
                            Annotation.scheme == "npclassifier",
                        )
                    )
                    s.commit()

                pairs: list[tuple[int, str]] = [
                    (c.id, c.smiles) for c in comps if c.smiles
                ]

            if not pairs:
                pbar.update(len(chunk_ids))  # all missing SMILES; advance bar
                continue

            # 3) Concurrently fetch NPClassifier results for this chunk
            results = asyncio.run(
                _fetch_many(
                    pairs,
                    concurrency=concurrency,
                    rate_per_sec=rate_per_sec,
                    timeout_s=timeout_s,
                    retries=retries,
                    pbar=pbar,
                )
            )

            # 4) Build bulk rows and insert
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
                continue

            try:
                with SessionLocal() as s:
                    stmt = (
                        insert(Annotation)
                        .values(rows_to_insert)
                        .on_conflict_do_nothing(
                            constraint="uq_annotation_target_scheme_key_value"
                        )
                    )
                    res = s.execute(stmt)
                    s.commit()
                    # Rowcount may be -1 depending on backend; defensively derive from attempted size
                    inserted = res.rowcount if (res.rowcount or 0) > 0 else 0
                    if inserted == 0:
                        # Fallback heuristic when rowcount unavailable
                        inserted = len(rows_to_insert)
                    total_inserted += inserted
                    logger.info(
                        f"Chunk {start//chunk_size + 1}: attempted {len(rows_to_insert)} inserts"
                    )
            except SQLAlchemyError as e:
                logger.exception(f"DB error on chunk {start}-{end}: {e}")
                # Continue to next chunk but do not lose progress bar state

    return total_inserted
