"""ETL functions for MIBiG data source: downloading GBK files and loading JSON/GBK data into the database."""

from __future__ import annotations
import json
import hashlib
import logging
import re
import random

import asyncio
import aiohttp
from aiohttp import (
    ClientResponseError,
    ClientConnectorError,
    ClientOSError,
    ClientPayloadError,
)
from pathlib import Path
from sqlalchemy import select
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.dialects.postgresql import insert
from tqdm import tqdm

from bionexus.config import DEFAULT_MAX_BYTES_GBK
from bionexus.etl.chemistry import get_atom_counts, smiles_to_inchikey
from bionexus.db.engine import SessionLocal
from bionexus.db.models import Annotation, Compound, CompoundRecord, GenBankRegion


logger = logging.getLogger(__name__)


MIBIG_URL_TEMPLATE = "https://mibig.secondarymetabolites.org/repository/{accession}.{version}/generated/{accession}.gbk"


async def _fetch_one(
    session: aiohttp.ClientSession,
    sem: asyncio.Semaphore,
    url: str,
    out_path: Path,
    retries: int = 5,
    timeout_s: int = 20,
) -> Path | None:
    """
    Download a single file with retries and exponential backoff.

    :param session: aiohttp ClientSession
    :param sem: asyncio Semaphore for concurrency control
    :param url: URL to download
    :param out_path: Path to save the downloaded file
    :param retries: number of retry attempts
    :param timeout_s: per-request timeout in seconds
    :return: Path to the downloaded file, or None if failed or 404
    """
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path  # already downloaded

    tmp_path = out_path.with_suffix(".part")

    attempt = 0
    while attempt < retries:
        attempt += 1
        try:
            async with sem:
                timeout = aiohttp.ClientTimeout(total=timeout_s)
                async with session.get(url, timeout=timeout) as resp:
                    if resp.status == 404:
                        return None
                    resp.raise_for_status()
                    tmp_path.parent.mkdir(parents=True, exist_ok=True)
                    with open(tmp_path, "wb") as f:
                        async for chunk in resp.content.iter_chunked(1 << 14):
                            f.write(chunk)
                    tmp_path.replace(out_path)
                    return out_path
        except (
            ClientConnectorError,
            ClientOSError,
            ClientPayloadError,
            asyncio.TimeoutError,
        ) as e:
            err = f"{e.__class__.__name__}: {e}"
        except ClientResponseError as e:
            err = f"HTTP {e.status}: {e.message}"
        except Exception as e:
            err = str(e)

        # Backoff
        await asyncio.sleep(min(2 ** (attempt - 1), 16) + random.uniform(0, 0.3))

    logger.error(f"[fail] {url} ({err})")
    if tmp_path.exists():
        tmp_path.unlink(missing_ok=True)
    return None


async def _download_all(
    accessions: list[tuple[str, str]],
    outdir: str | Path,
    concurrency=10,
    retries=5,
    timeout_s=20,
) -> list[Path]:
    """
    Download multiple MIBiG GBK files concurrently.

    :param accessions: list of (accession, version) tuples
    :param outdir: output directory
    :param concurrency: number of concurrent downloads
    :param retries: number of retry attempts per file
    :param timeout_s: per-request timeout in seconds
    :return: list of Paths to successfully downloaded files
    """
    connector = aiohttp.TCPConnector(limit=concurrency, limit_per_host=5)
    sem = asyncio.Semaphore(concurrency)
    headers = {
        "User-Agent": "BionexusDownloader/1.0 (+contact: scripts@local)",
        "Accept": "*/*",
        "Accept-Encoding": "gzip, deflate, br",
    }

    async with aiohttp.ClientSession(connector=connector, headers=headers) as session:
        tasks = []
        for accession, version in accessions:
            url = MIBIG_URL_TEMPLATE.format(accession=accession, version=version)
            out_path = Path(outdir) / f"{accession}.gbk"
            tasks.append(_fetch_one(session, sem, url, out_path, retries, timeout_s))

        results = []
        with tqdm(total=len(tasks), desc="Downloading GBKs") as pbar:
            for coro in asyncio.as_completed(tasks):
                res = await coro
                if res:
                    results.append(res)
                pbar.update(1)

        return results


def download_mibig_gbks(
    accessions: list[tuple[str, str]],
    outdir: str | Path,
    *,
    concurrency: int = 10,
    retries: int = 5,
    timeout: int = 20,
) -> list[Path]:
    """
    Download MIBiG GBK files for given (accession, version) tuples.

    :param accessions: list of (accession, version) tuples
    :param outdir: output directory
    :param concurrency: number of concurrent downloads
    :param retries: number of retry attempts per file
    :param timeout: per-request timeout in seconds
    :return: list of Paths to successfully downloaded files
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    return asyncio.run(_download_all(accessions, outdir, concurrency, retries, timeout))


def _read_gbk_text(path: str) -> str:
    """
    Read GenBank file text with UTF-8 encoding and error replacement.

    :param path: path to GenBank file
    :return: file content as string
    """
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        return fh.read()


def _sha256_bytes(b: bytes) -> str:
    """
    Compute SHA256 hash of bytes.

    :param b: input bytes
    :return: hex digest of SHA256 hash
    """
    return hashlib.sha256(b).hexdigest()


def extract_taxonomy_info(record) -> dict[str, str | list[str] | None]:
    """
    Extract taxonomy information from a Biopython SeqRecord's annotations.

    :param record: Biopython SeqRecord
    :return: dict with keys 'organism', 'taxonomy', 'genus', 'species', 'strain'
    """
    organism = record.annotations.get("organism", "")
    taxonomy = record.annotations.get("taxonomy", [])
    genus = taxonomy[-1] if taxonomy else None

    # Parse organism string more robustly
    parts = organism.split()
    genus_guess = parts[0] if parts else None
    species_guess = None
    strain_guess = None

    if len(parts) > 1:
        # Species is lowercase and not "sp."
        if re.match(r"^[a-z-]+$", parts[1]) and parts[1] != "sp.":
            species_guess = parts[1]
            if len(parts) > 2:
                strain_guess = " ".join(parts[2:])
        else:
            # no clear species, maybe "Streptomyces sp."
            strain_guess = " ".join(parts[1:])

    return {
        "organism": organism or None,
        "taxonomy": taxonomy,
        "genus": genus or genus_guess,
        "species": species_guess,
        "strain": strain_guess,
    }


def extract_region_products(record):
    """
    Extract product annotations from 'region' features in a Biopython SeqRecord.

    :param record: Biopython SeqRecord
    :return: list of unique product strings
    """
    products = []
    for feature in record.features:
        if feature.type == "region":
            prod = feature.qualifiers.get("product", [])
            if prod:
                products.extend(prod)

    # Deduplicate while preserving order
    seen = set()
    unique_products = [p for p in products if not (p in seen or seen.add(p))]

    return unique_products


def load_mibig_files(
    json_paths: list[str] | None, gbk_paths: str | list[str], chunk_size: int = 1000
) -> tuple[int, int, int, int]:
    """
    Load MIBiG compounds from JSON files and GenBank regions from GBK files into the database.

    :param json_paths: list of paths to MIBiG JSON files
    :param gbk_paths: list of paths to MIBiG GenBank files
    :param chunk_size: number of records to process per database commit
    :return: tuple of counts (new_compounds, new_records, new_regions, new_annotations)
    """
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    source = "mibig"

    new_compounds = 0
    new_records = 0
    new_regions = 0
    new_annotations = 0

    batch_compounds: list[Compound] = []
    batch_records: list[CompoundRecord] = []
    batch_regions: list[GenBankRegion] = []

    # Collect annotations as dicts and upsert in bulk per chunk
    batch_ann_dicts: list[
        dict
    ] = []  # each has {"_comp_obj"/"_gbk_obj", "scheme", "key", "value"}
    pending_ann = (
        0  # count staged annotations to trigger chunk commit even if only ann present
    )

    # De-dupe within this load (by inchikey and by (inchikey, source, ext_id) pre-flush)
    seen_inchikey: set[str] = set()
    seen_record_key: set[tuple[str, str, str]] = set()  # (inchikey, source, ext_id)
    seen_region_key: set[str] = set()  # (sha256)
    seen_annotation_key: set[tuple[str, str, str, str]] = (
        set()
    )  # (inchikey/sha256, scheme, key, value)

    with SessionLocal() as s:
        # Compounds from JSON files
        for path in tqdm(json_paths or [], desc="Loading MIBiG JSONs"):
            data = json.load(open(path))
            ext_id = data["accession"]

            for c in data.get("compounds", []):
                name = c.get("name")
                smiles = c.get("structure")
                mol_formula = c.get("formula")
                exact_mass = c.get("mass")

                if not smiles:
                    continue

                if smiles:
                    inchikey: str | None = smiles_to_inchikey(smiles)

                    # Must have an inchikey to unify; skip if missing
                    if not inchikey:
                        continue

                    # Compute atom counts
                    atom_counts = get_atom_counts(smiles)

                    # 1) Find or stage canonical compound by inchikey
                    comp = s.scalars(
                        select(Compound).where(Compound.inchikey == inchikey)
                    ).first()
                    if not comp:
                        # Also avoid creating the same compound twice in one batch before flush
                        if inchikey in seen_inchikey:
                            # If we have already staged it in this batch, fit it from batch list
                            comp = next(
                                (c for c in batch_compounds if c.inchikey == inchikey),
                                None,
                            )
                        if not comp:
                            comp = Compound(
                                inchikey=inchikey,
                                smiles=smiles,
                                mol_formula=mol_formula,
                                exact_mass=exact_mass,
                                **(atom_counts or {}),
                            )
                            batch_compounds.append(comp)
                            seen_inchikey.add(inchikey)
                            new_compounds += 1

                    else:
                        # Optionally update canonical props to fill in gaps only
                        if comp.smiles is None and smiles:
                            comp.smiles = smiles
                        if comp.mol_formula is None and mol_formula:
                            comp.mol_formula = mol_formula
                        if comp.exact_mass is None and exact_mass:
                            comp.exact_mass = exact_mass
                        # Update atom counts if missing
                        if atom_counts:
                            for k, v in atom_counts.items():
                                if getattr(comp, k) is None and v is not None:
                                    setattr(comp, k, v)

                    # 2) Ensure CompoundRecord, batch de-dupe by inchikey+source+ext_id (pre-flush)
                    rk = (inchikey, source, ext_id)
                    if rk in seen_record_key:
                        pass
                    else:
                        seen_record_key.add(rk)
                        rec_exists = s.scalars(
                            select(CompoundRecord).where(
                                CompoundRecord.compound_id
                                == comp.id,  # None for new; query will skip
                                CompoundRecord.source == source,
                                CompoundRecord.ext_id == ext_id,
                            )
                        ).first()

                        if not rec_exists:
                            batch_records.append(
                                CompoundRecord(
                                    compound=comp,  # safe even if comp not flushed
                                    source=source,
                                    ext_id=ext_id,
                                    name=name,
                                )
                            )
                            new_records += 1

                    # Flush if batch full
                    if (
                        len(batch_compounds)
                        + len(batch_records)
                        + len(batch_regions)
                        + pending_ann
                    ) >= chunk_size:
                        try:
                            if batch_compounds:
                                s.add_all(batch_compounds)
                            if batch_records:
                                s.add_all(batch_records)
                            if batch_regions:
                                s.add_all(batch_regions)
                            s.flush()  # assign ids to new objects for FKs in annotations

                            # Bulk upsert annotations
                            # NOTE: currently no annotations from JSONs

                            s.commit()

                        except SQLAlchemyError as e:
                            s.rollback()
                            logger.error(f"Error during MIBiG load chunk commit: {e}")
                            raise

                        batch_compounds.clear()
                        batch_records.clear()
                        batch_regions.clear()
                        batch_ann_dicts.clear()
                        seen_inchikey.clear()
                        seen_record_key.clear()
                        seen_region_key.clear()
                        seen_annotation_key.clear()
                        pending_ann = 0

        # Regions from GenBank files
        for path in tqdm(gbk_paths or [], desc="Loading MIBiG GenBank files"):
            ext_id = path.split("/")[-1].replace(".gbk", "").replace(".gb", "")

            gbk_text = _read_gbk_text(path)
            enc = gbk_text.encode("utf-8", errors="replace")
            size_bytes = len(enc)
            sha256 = _sha256_bytes(enc)

            # Parse taxonomy annotations from GenBank annotations
            gbk_records = list(SeqIO.parse(path, "genbank"))
            if len(gbk_records) != 1:
                logger.warning(
                    f"Skipping MIBiG GenBank {ext_id} due to unexpected record count {len(gbk_records)} != 1"
                )
                continue
            gbk_record: SeqRecord = gbk_records[0]
            annotations: list[tuple[str, str, str]] = []
            tax_info = extract_taxonomy_info(gbk_record)
            if tax_info["genus"]:
                genus = tax_info["genus"]
                annotations.append(("taxonomy", "genus", genus))
            if tax_info["genus"] and tax_info["species"]:
                species = f"{tax_info['genus']} {tax_info['species']}"
                annotations.append(("taxonomy", "species", species))

            # Parse biosynthetic class annotations from GenBank features
            region_products = extract_region_products(gbk_record)
            for rp in region_products:
                annotations.append(("biosynthesis", "product", rp))

            # Skip if too large
            if size_bytes > DEFAULT_MAX_BYTES_GBK:
                logger.warning(
                    f"Skipping MIBiG GenBank {ext_id} due to size {size_bytes} > {DEFAULT_MAX_BYTES_GBK}"
                )
                continue

            # 1) Find or create GenBankRegion by sha256
            region = s.scalars(
                select(GenBankRegion).where(GenBankRegion.sha256 == sha256)
            ).first()
            if not region:
                # Also avoid creating same region twice in one batch before flush
                if sha256 in seen_region_key:
                    # If we have already staged in this batch, fetch it from batch list
                    region = next(
                        (r for r in batch_regions if r.sha256 == sha256), None
                    )
                if not region:
                    region = GenBankRegion(
                        source=source,
                        ext_id=ext_id,
                        gbk_text=gbk_text,
                        size_bytes=size_bytes,
                        sha256=sha256,
                    )
                    batch_regions.append(region)
                    seen_region_key.add(sha256)
                    new_regions += 1
            else:
                # Optionally update props to fill in gaps only
                if region.gbk_text is None and gbk_text:
                    region.gbk_text = gbk_text
                if region.size_bytes is None and size_bytes:
                    region.size_bytes = size_bytes
                if region.ext_id is None and ext_id:
                    region.ext_id = ext_id

            # 2) Annotations from GBK features
            for scheme, key, value in annotations:
                if not value:
                    continue

                # Avoid annotation dupes in this batch
                ann_key = (sha256, scheme, key, value)
                if ann_key in seen_annotation_key:
                    continue
                seen_annotation_key.add(ann_key)

                batch_ann_dicts.append(
                    {
                        "_gbk_obj": region,  # safe even if region not flushed
                        "scheme": scheme,
                        "key": key,
                        "value": value,
                    }
                )
                pending_ann += 1

            # Flush if batch full
            if (
                len(batch_compounds)
                + len(batch_records)
                + len(batch_regions)
                + pending_ann
            ) >= chunk_size:
                try:
                    if batch_compounds:
                        s.add_all(batch_compounds)
                    if batch_records:
                        s.add_all(batch_records)
                    if batch_regions:
                        s.add_all(batch_regions)
                    s.flush()  # assign ids to new objects for FKs in annotations

                    # Bulk upsert annotations
                    if batch_ann_dicts:
                        ann_values = [
                            {
                                "genbank_region_id": d["_gbk_obj"].id,
                                "scheme": d["scheme"],
                                "key": d["key"],
                                "value": d["value"],
                                "metadata_json": None,
                            }
                            for d in batch_ann_dicts
                        ]
                        if ann_values:
                            stmt = (
                                insert(Annotation)
                                .values(ann_values)
                                .on_conflict_do_nothing(
                                    constraint="uq_annotation_target_scheme_key_value"
                                )
                                .returning(Annotation.id)
                            )
                            inserted_ids = s.execute(stmt).scalars().all()
                            new_annotations += len(inserted_ids)

                    s.commit()

                except SQLAlchemyError as e:
                    s.rollback()
                    logger.error(f"Error during MIBiG load chunk commit: {e}")
                    raise

                batch_compounds.clear()
                batch_records.clear()
                batch_regions.clear()
                batch_ann_dicts.clear()
                seen_inchikey.clear()
                seen_record_key.clear()
                seen_region_key.clear()
                seen_annotation_key.clear()
                pending_ann = 0

        # Final flush
        if batch_compounds or batch_records or batch_regions or batch_ann_dicts:
            try:
                if batch_compounds:
                    s.add_all(batch_compounds)
                if batch_records:
                    s.add_all(batch_records)
                if batch_regions:
                    s.add_all(batch_regions)
                s.flush()  # assign ids to new objects for FKs in annotations

                # Bulk upsert annotations
                if batch_ann_dicts:
                    ann_values = [
                        {
                            "genbank_region_id": d["_gbk_obj"].id,
                            "scheme": d["scheme"],
                            "key": d["key"],
                            "value": d["value"],
                            "metadata_json": None,
                        }
                        for d in batch_ann_dicts
                    ]
                    if ann_values:
                        stmt = (
                            insert(Annotation)
                            .values(ann_values)
                            .on_conflict_do_nothing(
                                constraint="uq_annotation_target_scheme_key_value"
                            )
                            .returning(Annotation.id)
                        )
                        inserted_ids = s.execute(stmt).scalars().all()
                        new_annotations += len(inserted_ids)

                s.commit()

            except SQLAlchemyError as e:
                s.rollback()
                logger.error(f"Error during MIBiG final load commit: {e}")
                raise

    return new_compounds, new_records, new_regions, new_annotations
