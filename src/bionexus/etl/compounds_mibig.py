"""ETL for loading MIBiG compound data, references, and annotations into BioNexus."""

import json
import uuid
import logging
from pathlib import Path
from dataclasses import dataclass
from typing import Generator, TypeVar, Iterator

from tqdm import tqdm
from rdkit import RDLogger

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import insert as pg_insert
from sqlalchemy.exc import SQLAlchemyError

from biocracker.utils.download import download_and_prepare

from retromol.model.rules import RuleSet
from retromol.model.result import Result
from retromol.io.streaming import run_retromol_stream
from retromol.chem.mol import smiles_to_mol, mol_to_inchikey
from retromol.fingerprint.fingerprint import FingerprintGenerator

from bionexus.db.engine import SessionLocal
from bionexus.db.models import Compound, Reference, Annotation, compound_reference, compound_annotation
from bionexus.etl.compound import calculate_compound_props


# Disable RDKit warnings
RDLogger.DisableLog("rdApp.*")


MIBIG_URL = r"https://dl.secondarymetabolites.org/mibig/mibig_json_4.0.tar.gz"


log = logging.getLogger(__name__)


T = TypeVar("T")


ruleset = RuleSet.load_default()
generator = FingerprintGenerator(ruleset.matching_rules)


@dataclass(frozen=True)
class MIBiGCompoundJob:
    """
    Single MIBiG compound job parsed from MIBiG JSON data.
    """

    guid: str  # internal job identifier

    # Structure
    inchikey: str
    smiles: str

    # Unique identifier; database reference (multiple compounds might have same MIBiG ID since it references a producing BGC)
    mibig_id: str

    # Name references
    name: str | None

    # Taxonomy annotations
    taxonomy_genus: str | None
    taxonomy_species: str | None


def iter_jobs(jobs: list[MIBiGCompoundJob]) -> Generator[dict[str, str], None, None]:
    """
    Iterate over MIBiG compound jobs and yield data for database loading.
    """
    for j in jobs:
        yield {
            "guid": j.guid,
            "smiles": j.smiles,
            "mibig_id": j.mibig_id,
        }


def iter_chunks(items: list[T], chunk_size: int) -> Iterator[list[T]]:
    for i in range(0, len(items), chunk_size):
        yield items[i : i + chunk_size]


def upsert_compounds_and_get_ids(s, compound_rows: list[dict]) -> dict[str, int]:
    """
    Insert compounds (do nothing on conflict) then return inchikey->id for all rows.
    """
    if not compound_rows:
        return {}

    stmt = (
        pg_insert(Compound)
        .values(compound_rows)
        .on_conflict_do_nothing(index_elements=[Compound.inchikey])
    )
    s.execute(stmt)
    s.flush()

    inchikeys = [r["inchikey"] for r in compound_rows]
    rows = s.execute(
        sa.select(Compound.inchikey, Compound.id).where(Compound.inchikey.in_(inchikeys))
    ).all()
    return {ik: cid for ik, cid in rows}


def upsert_references_and_get_ids(s, ref_rows: list[dict]) -> dict[tuple[str, str, str], int]:
    """
    Upsert Reference rows and return (name, db, dbid)->id for all.
    """
    if not ref_rows:
        return {}

    stmt = (
        pg_insert(Reference)
        .values(ref_rows)
        .on_conflict_do_nothing(
            constraint="ux_reference_name_dbname_dbid"
        )
    )
    s.execute(stmt)
    s.flush()

    keys = {(r["name"], r["database_name"], r["database_identifier"]) for r in ref_rows}
    rows = s.execute(
        sa.select(Reference.name, Reference.database_name, Reference.database_identifier, Reference.id).where(
            sa.tuple_(Reference.name, Reference.database_name, Reference.database_identifier).in_(list(keys))
        )
    ).all()
    return {(n, db, dbid): rid for n, db, dbid, rid in rows}


def upsert_annotations_and_get_ids(s, ann_rows: list[dict]) -> dict[tuple[str, str, str], int]:
    """
    Upsert Annotation rows and return (scheme, key, value)->id for all.
    """
    if not ann_rows:
        return {}

    stmt = (
        pg_insert(Annotation)
        .values(ann_rows)
        .on_conflict_do_nothing(
            constraint="ux_annotation_scheme_key_value"
        )
    )
    s.execute(stmt)
    s.flush()

    keys = {(r["scheme"], r["key"], r["value"]) for r in ann_rows}
    rows = s.execute(
        sa.select(Annotation.scheme, Annotation.key, Annotation.value, Annotation.id).where(
            sa.tuple_(Annotation.scheme, Annotation.key, Annotation.value).in_(list(keys))
        )
    ).all()
    return {(sch, k, v): aid for sch, k, v, aid in rows}


def insert_links_do_nothing(s, table, rows: list[dict], conflict_cols: list):
    """
    Bulk insert into link table with do-nothing on conflict.
    """
    if not rows:
        return
    stmt = (
        pg_insert(table)
        .values(rows)
        .on_conflict_do_nothing(index_elements=conflict_cols)
    )
    s.execute(stmt)


def refs_for_job(job: MIBiGCompoundJob) -> list[dict]:
    """
    Generate reference rows for a given MIBiG compound job.
    """
    db = "MIBiG"
    dbid = job.mibig_id
    out: list[dict] = []

    def norm_name(x) -> str | None:
        if x is None:
            return None
        if isinstance(x, dict):
            x = x.get("name")
        x = str(x).strip()
        if not x:
            return None
        return x[:256]  # Reference.name is varchar(256)

    n = norm_name(job.name)
    if n:
        out.append({"name": n, "database_name": db, "database_identifier": dbid})

    return out


def annotations_for_job(job: MIBiGCompoundJob) -> list[dict]:
    out: list[dict] = []

    def ann(scheme: str, key: str, value: str | None):
        if not value:
            return
        out.append({
            "scheme": scheme[:64],
            "key": key[:64],
            "value": str(value).strip()[:256],  # Annotation.value is varchar(256)
        })

    # Taxonomy strings
    if job.taxonomy_genus:
        ann("taxonomy", "genus", job.taxonomy_genus)
    if job.taxonomy_species:
        ann("taxonomy", "species", job.taxonomy_species)
    return out


def load_compounds_mibig(
    workdir: Path | str,
    workers: int | None = None,
    batch_size: int = 1_000,
    chunk_size: int = 1_000,
) -> None:
    """
    Load MIBiG compounds into the database.

    - Runs RetroMol in batches of `batch_size` jobs (CPU side)
    - Flushes to DB every `chunk_size` processed compounds (I/O side)
    - Inserts:
        * Compound rows
        * Reference rows (guid + MIBiG accession + name) + compound_reference links
        * Annotation rows (taxonomy) + compound_annotation links
    """
    workdir = Path(workdir).expanduser()
    workdir.mkdir(parents=True, exist_ok=True)

    # ----------------------------
    # Download + load MIBiG JSONs
    # ----------------------------
    mibig_path = download_and_prepare(url=MIBIG_URL, cache_dir=workdir)
    mibig_json_iter = mibig_path.rglob("*.json")

    # ----------------------------
    # Parse into jobs + lookup map
    # ----------------------------
    jobs: list[MIBiGCompoundJob] = []
    job_by_guid: dict[str, MIBiGCompoundJob] = {}

    for path in tqdm(mibig_json_iter, desc="Parsing MIBiG JSONs"):

        with open(path, "r", encoding="utf-8") as f:
            item = json.load(f)

        # Parse out MIBiG accession for the producing BGC
        accession = item.get("accession")
        if not accession:
            raise ValueError(f"missing accession in MIBiG JSON: {path}")
        version = item.get("version")
        full_accession = f"{accession}.{version}" if version else accession

        # Parse out taxonomy information for the producing organism that contains this BGC
        taxonomy = item.get("taxonomy", {}) or {}
        genus = None
        species = None
        if taxonomy:
            species = taxonomy.get("name")
            if species:
                genus = species.split(" ")[0]

        # Loop over compounds BGC produces
        compounds = item.get("compounds") or []
        for compound in compounds:

            guid = str(uuid.uuid4())

            name = compound.get("name")
            if not name:
                name = f"Compound produced by {full_accession}"
            smiles = compound.get("structure")
            if not smiles:
                continue
                
            mol = smiles_to_mol(smiles)
            inchikey = mol_to_inchikey(mol)

            job = MIBiGCompoundJob(
                guid=guid,
                inchikey=inchikey,
                smiles=smiles,
                mibig_id=full_accession,
                name=name,
                taxonomy_genus=genus,
                taxonomy_species=species,
            )

            if not job.inchikey or not job.smiles or not job.mibig_id or not job.guid:
                log.warning("skipping malformed MIBiG entry (missing required fields)")
                continue

            jobs.append(job)
            job_by_guid[guid] = job

        continue

    # No filtering, prefer double work and good referencing in database
    # with SessionLocal() as s:
    #     existing = set(s.execute(sa.select(Compound.inchikey)).scalars().all())
    
    # jobs = [j for j in jobs if j.inchikey not in existing]
    # job_by_guid = {j.guid: j for j in jobs}
    # log.info(f"Prepared {len(jobs)} MIBiG compound jobs for processing (after filtering existing).")

    # ----------------------------
    # Progress + counters
    # ----------------------------
    pbar = tqdm(total=len(jobs), desc="Processing MIBiG compounds")
    result_counts = {"successes": 0, "failures": 0, "errors": 0}

    inserted_compounds = 0
    failed_db = 0
    failed_processing = 0

    # ----------------------------
    # Chunk buffers (DB flush unit)
    # ----------------------------
    seen_inchikey: set[str] = set()
    chunk_rows: list[dict] = []
    chunk_jobs: list[MIBiGCompoundJob] = []

    for batch_idx, batch_jobs in enumerate(iter_chunks(jobs, batch_size), start=1):
        batch_iter = iter_jobs(batch_jobs)
        batch_results: list[tuple[MIBiGCompoundJob, Result]] = []

        # ------------------------
        # RetroMol: produce results
        # ------------------------
        for evt in run_retromol_stream(
            ruleset=ruleset,
            row_iter=batch_iter,
            smiles_col="smiles",
            workers=(workers or 1),
        ):
            if evt.error is not None:
                log.error(evt.error)
                result_counts["errors"] += 1

            elif evt.result is not None:
                try:
                    result = Result.from_dict(evt.result)
                    guid = result.submission.props["guid"]
                    job = job_by_guid[guid]
                    batch_results.append((job, result))
                    result_counts["successes"] += 1
                except Exception as e:
                    log.error(f"failed to deserialize/map result: {e}")
                    result_counts["failures"] += 1

            else:
                log.error("received empty result with no error message")
                result_counts["failures"] += 1

            pbar.update(1)

        # ------------------------
        # DB flush for this batch
        # ------------------------
        log.info(f"Flushing RetroMol batch {batch_idx} ({len(batch_results)} results) to DB...")

        with SessionLocal() as s:
            for job, result in batch_results:
                try:
                    inchikey = job.inchikey
                    smiles = job.smiles

                    # chunk-local de-dupe to reduce work (DB enforces global uniqueness anyway)
                    if inchikey in seen_inchikey:
                        continue
                    seen_inchikey.add(inchikey)

                    mol = smiles_to_mol(smiles)
                    props = calculate_compound_props(mol)

                    coverage = result.calculate_coverage()
                    fp_counted = [
                        float(x)
                        for x in generator.fingerprint_from_result(result, num_bits=1024, counted=True)
                    ]
                    fp_binary = [float(int(x > 0)) for x in fp_counted]

                    chunk_rows.append(
                        {
                            "inchikey": inchikey,
                            "smiles": smiles,
                            "mol_weight": props.mol_weight,
                            "c_atom_count": props.c_atom_count,
                            "h_atom_count": props.h_atom_count,
                            "n_atom_count": props.n_atom_count,
                            "o_atom_count": props.o_atom_count,
                            "p_atom_count": props.p_atom_count,
                            "s_atom_count": props.s_atom_count,
                            "f_atom_count": props.f_atom_count,
                            "cl_atom_count": props.cl_atom_count,
                            "br_atom_count": props.br_atom_count,
                            "i_atom_count": props.i_atom_count,
                            "morgan_fp": props.morgan_fp,
                            "retromol_fp_counted": fp_counted,
                            "retromol_fp_binary": fp_binary,
                            "retromol": result.to_dict(),
                            "coverage": coverage,
                        }
                    )
                    chunk_jobs.append(job)

                    # ------------------------------------
                    # Flush chunk: compounds + refs + anns
                    # ------------------------------------
                    if len(chunk_rows) >= chunk_size:
                        try:
                            inchikey_to_cid = upsert_compounds_and_get_ids(s, chunk_rows)
                            inserted_compounds += len(inchikey_to_cid)

                            # collect refs/anns for THIS chunk
                            ref_rows: list[dict] = []
                            ann_rows: list[dict] = []
                            for j in chunk_jobs:
                                ref_rows.extend(refs_for_job(j))
                                ann_rows.extend(annotations_for_job(j))

                            ref_key_to_id = upsert_references_and_get_ids(s, ref_rows)
                            ann_key_to_id = upsert_annotations_and_get_ids(s, ann_rows)

                            # link rows
                            cref_rows: list[dict] = []
                            cann_rows: list[dict] = []

                            for j in chunk_jobs:
                                cid = inchikey_to_cid.get(j.inchikey)
                                if cid is None:
                                    continue

                                for r in refs_for_job(j):
                                    rid = ref_key_to_id[(r["name"], r["database_name"], r["database_identifier"])]
                                    cref_rows.append({"compound_id": cid, "reference_id": rid})

                                for a in annotations_for_job(j):
                                    aid = ann_key_to_id[(a["scheme"], a["key"], a["value"])]
                                    cann_rows.append({"compound_id": cid, "annotation_id": aid})

                            insert_links_do_nothing(
                                s,
                                compound_reference,
                                cref_rows,
                                conflict_cols=[compound_reference.c.compound_id, compound_reference.c.reference_id],
                            )
                            insert_links_do_nothing(
                                s,
                                compound_annotation,
                                cann_rows,
                                conflict_cols=[compound_annotation.c.compound_id, compound_annotation.c.annotation_id],
                            )

                            s.commit()

                        except SQLAlchemyError as e:
                            s.rollback()
                            failed_db += len(chunk_rows)
                            log.error(f"database error during chunk flush: {e}")

                        finally:
                            chunk_rows.clear()
                            chunk_jobs.clear()
                            seen_inchikey.clear()

                except Exception as e:
                    failed_processing += 1
                    log.warning(f"failed to process compound inchikey={job.inchikey}: {e}")
                    continue

            # ------------------------------------
            # Flush remaining chunk rows at end of batch
            # ------------------------------------
            if chunk_rows:
                try:
                    inchikey_to_cid = upsert_compounds_and_get_ids(s, chunk_rows)
                    inserted_compounds += len(inchikey_to_cid)

                    ref_rows: list[dict] = []
                    ann_rows: list[dict] = []
                    for j in chunk_jobs:
                        ref_rows.extend(refs_for_job(j))
                        ann_rows.extend(annotations_for_job(j))

                    ref_key_to_id = upsert_references_and_get_ids(s, ref_rows)
                    ann_key_to_id = upsert_annotations_and_get_ids(s, ann_rows)

                    cref_rows: list[dict] = []
                    cann_rows: list[dict] = []

                    for j in chunk_jobs:
                        cid = inchikey_to_cid.get(j.inchikey)
                        if cid is None:
                            continue

                        for r in refs_for_job(j):
                            rid = ref_key_to_id[(r["name"], r["database_name"], r["database_identifier"])]
                            cref_rows.append({"compound_id": cid, "reference_id": rid})

                        for a in annotations_for_job(j):
                            aid = ann_key_to_id[(a["scheme"], a["key"], a["value"])]
                            cann_rows.append({"compound_id": cid, "annotation_id": aid})

                    insert_links_do_nothing(
                        s,
                        compound_reference,
                        cref_rows,
                        conflict_cols=[compound_reference.c.compound_id, compound_reference.c.reference_id],
                    )
                    insert_links_do_nothing(
                        s,
                        compound_annotation,
                        cann_rows,
                        conflict_cols=[compound_annotation.c.compound_id, compound_annotation.c.annotation_id],
                    )

                    s.commit()

                except SQLAlchemyError as e:
                    s.rollback()
                    failed_db += len(chunk_rows)
                    log.error(f"database error during final chunk flush: {e}")

                finally:
                    chunk_rows.clear()
                    chunk_jobs.clear()
                    seen_inchikey.clear()

    pbar.close()

    log.info(
        "MIBiG compound loading completed: "
        f"{result_counts['successes']} successes, "
        f"{result_counts['failures']} failures, "
        f"{result_counts['errors']} errors. "
        f"Inserted compounds (seen in chunk selects): {inserted_compounds}. "
        f"Failed processing: {failed_processing}. Failed DB flush rows: {failed_db}."
    )
