from __future__ import annotations
from typing import List
from sqlalchemy import text, bindparam
from sqlalchemy.dialects.postgresql import BIT
from bionexus.db.engine import SessionLocal
from pgvector.sqlalchemy import Vector

def jaccard_search_exact(bitstr: str, top_k: int = 50):
    """
    Return: id, name (first found in compound_record), jacc
    """
    sql = text("""
    SELECT c.id,
           n.name,
           CASE
             WHEN length(replace((c.fp_morgan_b2048_r2_bit | :qb)::text, '0','')) = 0
             THEN 1.0
             ELSE length(replace((c.fp_morgan_b2048_r2_bit & :qb)::text, '0',''))::float
                  / NULLIF(length(replace((c.fp_morgan_b2048_r2_bit | :qb)::text, '0','')), 0)
           END AS jacc
    FROM compound c
    LEFT JOIN LATERAL (
        SELECT r.name
        FROM compound_record r
        WHERE r.compound_id = c.id AND r.name IS NOT NULL
        ORDER BY r.source, r.ext_id
        LIMIT 1
    ) AS n ON TRUE
    WHERE c.fp_morgan_b2048_r2_bit IS NOT NULL
    ORDER BY jacc DESC
    LIMIT :k
    """).bindparams(
        bindparam("qb", type_=BIT(2048)),
        bindparam("k")
    )

    with SessionLocal() as s:
        rows = s.execute(sql, {"qb": bitstr, "k": top_k}).mappings().all()
        return rows

def jaccard_search_hybrid(bitstr: str, qvec: List[float], top_k: int = 50, cand_k: int = 2000):
    """
    Return: id, name (first found in compound_record), jacc
    """
    sql = text("""
    WITH cand AS (
      SELECT c.id, c.fp_morgan_b2048_r2_bit
      FROM compound c
      WHERE c.fp_morgan_b2048_r2_bit IS NOT NULL
      ORDER BY c.fp_morgan_b2048_r2_vec <#> :qv   -- inner-product distance
      LIMIT :cand_k
    )
    SELECT c.id,
           n.name,
           CASE
             WHEN length(replace((c.fp_morgan_b2048_r2_bit | :qb)::text, '0','')) = 0
             THEN 1.0
             ELSE length(replace((c.fp_morgan_b2048_r2_bit & :qb)::text, '0',''))::float
                  / NULLIF(length(replace((c.fp_morgan_b2048_r2_bit | :qb)::text, '0','')), 0)
           END AS jacc
    FROM cand c
    LEFT JOIN LATERAL (
        SELECT r.name
        FROM compound_record r
        WHERE r.compound_id = c.id AND r.name IS NOT NULL
        ORDER BY r.source, r.ext_id
        LIMIT 1
    ) AS n ON TRUE
    ORDER BY jacc DESC
    LIMIT :k
    """).bindparams(
        bindparam("qb", type_=BIT(2048)),
        bindparam("qv", type_=Vector(2048)),
        bindparam("k"),
        bindparam("cand_k"),
    )

    with SessionLocal() as s:
        return s.execute(sql, {"qb": bitstr, "qv": qvec, "k": top_k, "cand_k": cand_k}).mappings().all()
