"""ETL for annotation data."""

import csv
import logging
from pathlib import Path

import sqlalchemy as sa
from sqlalchemy.exc import SQLAlchemyError
from tqdm import tqdm

from bionexus.db.engine import SessionLocal
from bionexus.db.models import (
    Compound,
    CandidateCluster,
    Annotation,
    Reference,
    compound_annotation,
    candidate_cluster_annotation,
    compound_reference,
    candidate_cluster_reference,
)


log = logging.getLogger(__name__)


def resolve_target(
    s: sa.orm.Session,
    row: dict[str, str]
) -> tuple[str, int] | None:
    """
    Resolve the target entity from the row data.

    :param s: database session
    :param row: dictionary representing a row of data
    :return: tuple of (entity type, entity ID) or None if not found
    """
    typ = (row.get("type") or "").strip().lower()

    if typ == "compound":
        inchikey = (row.get("inchikey") or "").strip()
        if not inchikey:
            return None
        cid = s.execute(
            sa.select(Compound.id).where(Compound.inchikey == inchikey)
        ).scalar_one_or_none()
        return ("compound", int(cid)) if cid is not None else None
    
    elif typ == "candidate_cluster":
        try:
            record_name = row["record_name"].strip()
            file_name = row["file_name"].strip()
            start_bp = int(row["start_bp"].strip())
            end_bp = int(row["end_bp"].strip())
        except (KeyError, ValueError):
            return None
        
        ccid = s.execute(
            sa.select(CandidateCluster.id).where(
                CandidateCluster.record_name == record_name,
                CandidateCluster.file_name == file_name,
                CandidateCluster.start_bp == start_bp,
                CandidateCluster.end_bp == end_bp,
            )
        ).scalar_one_or_none()
        return ("candidate_cluster", int(ccid)) if ccid is not None else None
    
    return None


def upsert_annotation_id(
    s: sa.orm.Session,
    scheme: str,
    key: str,
    value: str
) -> int:
    """
    Upsert an annotation and return its ID.
    
    :param s: database session
    :param scheme: annotation scheme
    :param key: annotation key
    :param value: annotation value
    :return: annotation ID
    """
    stmt = (
        sa.dialects.postgresql.insert(Annotation)
        .values(scheme=scheme, key=key, value=value)
        .on_conflict_do_update(
            index_elements=[Annotation.scheme, Annotation.key, Annotation.value],
            set_={"value": value},  # no-op update to avoid changing anything
        )
        .returning(Annotation.id)
    )
    return int(s.execute(stmt).scalar_one())


def upsert_reference_id(
    s: sa.orm.Session,
    name: str,
    database: str,
    database_identifier: str
) -> int:
    """
    Upsert a reference and return its ID.

    :param s: database session
    :param name: reference name
    :param database: reference database name
    :param database_identifier: reference database identifier
    :return: reference ID
    """
    stmt = (
        sa.dialects.postgresql.insert(Reference)
        .values(
            name=name,
            database_name=database,
            database_identifier=database_identifier,
        )
        .on_conflict_do_update(
            index_elements=[
                Reference.name,
                Reference.database_name,
                Reference.database_identifier,
            ],
            set_={"name": name},  # no-op update to avoid changing anything
        )
        .returning(Reference.id)
    )
    return int(s.execute(stmt).scalar_one())


def flush_links(
    s: sa.orm.Session,
    ann_links: list[tuple[str, int, int]],  # (kind, target_id, ann_id)
    ref_links: list[tuple[str, int, int]],  # (kind, target_id, ref_id)
) -> None:
    """
    Flush accumulated annotation and reference links to the database.

    :param s: database session
    :param ann_links: list of annotation links to insert
    :param ref_links: list of reference links to insert
    """
    for kind, tid, aid in ann_links:
        if kind == "compound":
            s.execute(
                sa.dialects.postgresql.insert(compound_annotation)
                .values(compound_id=tid, annotation_id=aid)
                .on_conflict_do_nothing()
            )
        else:
            s.execute(
                sa.dialects.postgresql.insert(candidate_cluster_annotation)
                .values(candidate_cluster_id=tid, annotation_id=aid)
                .on_conflict_do_nothing()
            )

    for kind, tid, rid in ref_links:
        if kind == "compound":
            s.execute(
                sa.dialects.postgresql.insert(compound_reference)
                .values(compound_id=tid, reference_id=rid)
                .on_conflict_do_nothing()
            )
        else:
            s.execute(
                sa.dialects.postgresql.insert(candidate_cluster_reference)
                .values(candidate_cluster_id=tid, reference_id=rid)
                .on_conflict_do_nothing()
            )

    s.commit()


def load_annotations(filepath: Path | str, chunk_size: int = 2_000) -> None:
    """
    Load annotations from a file into the database.

    :param filepath: path to the file containing annotation data
    """
    if isinstance(filepath, str):
        filepath = Path(filepath)

    ann_links = []
    ref_links = []

    processed = 0
    skipped = 0
    failed = 0

    with SessionLocal() as s:
        try:
            with filepath.open("r", newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                
                for row in tqdm(reader):
                    try:
                        tgt = resolve_target(s, row)
                        if tgt is None:
                            skipped += 1
                            log.warning(f"skipping row, target not found: {row}")
                            continue

                        kind, target_id = tgt
                        table = row["table"].strip().lower()
                        scheme = row["scheme"].strip()
                        key = row["key"].strip()
                        value = row["value"].strip()

                        if table == "annotation":
                            ann_id = upsert_annotation_id(s, scheme, key, value)
                            ann_links.append((kind, target_id, ann_id))
                        
                        elif table == "reference":
                            ref_id = upsert_reference_id(
                                s,
                                name=value,
                                database=scheme,
                                database_identifier=key
                            )
                            ref_links.append((kind, target_id, ref_id))
                        else:
                            skipped += 1
                            log.warning(f"skipping row, unknown table '{table}': {row}")
                            continue

                        processed += 1

                        if len(ann_links) + len(ref_links) >= chunk_size:
                            try:
                                flush_links(s, ann_links, ref_links)
                            except Exception as e:
                                s.rollback()
                                failed += len(ann_links) + len(ref_links)
                                log.error(f"failed to flush links: {e}")
                            finally:
                                ann_links.clear()
                                ref_links.clear()

                    except Exception as e:
                        failed += 1
                        s.rollback()
                        ann_links.clear()
                        ref_links.clear()
                        log.error(f"failed to process row {row}: {type(e).__name__}: {e}")

            # Final flush
            if ann_links or ref_links:
                try:
                    flush_links(s, ann_links, ref_links)
                except Exception as e:
                    s.rollback()
                    failed += len(ann_links) + len(ref_links)
                    log.error(f"failed to flush final links: {e}")
        
        except SQLAlchemyError as e:
            s.rollback()
            log.error(f"database error occurred: {e}")
            raise

    log.info(f"processed rows: {processed}")
    log.info(f"skipped rows: {skipped}")
    log.info(f"failed rows: {failed}")
