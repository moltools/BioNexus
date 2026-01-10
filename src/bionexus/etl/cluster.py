"""ETL for candidate cluster data."""

import logging
from pathlib import Path

from tqdm import tqdm

import sqlalchemy as sa
from sqlalchemy.exc import SQLAlchemyError

from retromol.model.rules import RuleSet
from retromol.fingerprint.fingerprint import FingerprintGenerator

from biocracker.utils.json import iter_json
from biocracker.query.modules import LinearReadout

from bionexus.db.engine import SessionLocal
from bionexus.db.models import CandidateCluster


log = logging.getLogger(__name__)


def load_clusters(jsonl: Path | str, chunk_size: int = 10_000) -> None:
    """
    Load candidate clusters from a JSONL file into the database.

    :param jsonl: path to the JSONL file containing candidate cluster data
    """
    if isinstance(jsonl, str):
        jsonl = Path(jsonl)

    ruleset = RuleSet.load_default()
    generator = FingerprintGenerator(ruleset.matching_rules)

    inserted = 0
    duplicates = 0
    failed = 0

    seen_cluster_ids: set[frozenset[str, int, int]] = set()
    batch_rows: list[dict] = []

    def flush_batch(session, rows: list[dict]) -> int:
        """
        Flush a batch of rows to the database.

        :param session: database session
        :param rows: list of row dictionaries to insert
        :return: number of rows inserted
        """
        if not rows:
            return 0

        stmt = (
            sa.dialects.postgresql.insert(CandidateCluster)
            .values(rows)
            .on_conflict_do_nothing(index_elements=[
                CandidateCluster.record_name,
                CandidateCluster.file_name,
                CandidateCluster.start_bp,
                CandidateCluster.end_bp
            ])
            .returning(CandidateCluster.id)
        )
        res = session.execute(stmt)
        session.commit()
        return len(res.fetchall())  # returns one row per inserted record
    
    with SessionLocal() as s:
        for rec in tqdm(iter_json(jsonl, jsonl=True)):
            
            try:
                r = LinearReadout.from_dict(rec)
                
                # Batch level de-dupe of clusters
                cluster_id = frozenset((str(r.id), int(r.start), int(r.end)))
                if cluster_id in seen_cluster_ids:
                    duplicates += 1
                    continue
                seen_cluster_ids.add(cluster_id)

                # retromol_fp_counted_by_orf = [float(x) for x in generator.fingerprint_from_biocracker_readout(r, num_bits=1024, counted=True, by_orf=True)]
                # retromol_fp_binary_by_orf = [float(int(x > 0)) for x in retromol_fp_counted_by_orf]
                retromol_fp_counted_by_region = [float(x) for x in generator.fingerprint_from_biocracker_readout(r, num_bits=512, counted=True, by_orf=False)]
                # retromol_fp_binary_by_region = [float(int(x > 0)) for x in retromol_fp_counted_by_region]

                batch_rows.append({
                    "record_name": str(r.id),
                    "file_name": str(r.file_name),
                    "start_bp": int(r.start),
                    "end_bp": int(r.end),
                    # "retromol_fp_counted_by_orf": retromol_fp_counted_by_orf,
                    # "retromol_fp_binary_by_orf": retromol_fp_binary_by_orf,
                    "retromol_fp_counted_by_region": retromol_fp_counted_by_region,
                    # "retromol_fp_binary_by_region": retromol_fp_binary_by_region,
                    "biocracker": r.to_dict(),
                })

                if len(batch_rows) >= chunk_size:
                    try:
                        n_ins = flush_batch(s, batch_rows)
                        inserted += n_ins
                        duplicates += len(batch_rows) - n_ins
                    except SQLAlchemyError as e:
                        s.rollback()
                        failed += len(batch_rows)
                        log.error(f"database error during batch insert: {e}")
                        exit()
                    finally:
                        batch_rows.clear()
                        seen_cluster_ids.clear()

            except Exception as e:
                failed += 1
                log.error(f"failed to process record: {e}")
                continue

        # Flush any remaining rows
        try:
            n_ins = flush_batch(s, batch_rows)
            inserted += n_ins
            duplicates += len(batch_rows) - n_ins
        except SQLAlchemyError as e:
            s.rollback()
            failed += len(batch_rows)
            log.error(f"database error during final batch insert: {e}")

    log.info(f"total clusters inserted: {inserted}")
    log.info(f"total duplicate clusters skipped: {duplicates}")
    log.info(f"total failed clusters: {failed}")
            