from __future__ import annotations
from typing import List
from sqlalchemy import text, bindparam
from sqlalchemy.dialects.postgresql import BIT
from bionexus.db.engine import SessionLocal
from pgvector.sqlalchemy import Vector

def jaccard_search_exact(bitstr: str, top_k: int = 50):
    sql = text("""
    WITH top_compounds AS (
        SELECT 
            c.id,
            c.smiles,
            CASE
                WHEN length(replace((c.fp_morgan_b2048_r2_bit | :qb)::text, '0','')) = 0
                THEN 1.0
                ELSE length(replace((c.fp_morgan_b2048_r2_bit & :qb)::text, '0',''))::float
                    / NULLIF(length(replace((c.fp_morgan_b2048_r2_bit | :qb)::text, '0','')), 0)
            END AS jacc
        FROM compound c
        WHERE c.fp_morgan_b2048_r2_bit IS NOT NULL
        ORDER BY jacc DESC
        LIMIT :k
    )
    SELECT
        t.id,
        r.source,
        r.name,
        t.jacc,
        t.smiles
    FROM top_compounds t
    LEFT JOIN compound_record r
        ON r.compound_id = t.id
    WHERE r.name IS NOT NULL
    ORDER BY t.jacc DESC, r.source, r.name, t.id
    """).bindparams(
        bindparam("qb", type_=BIT(2048)),
        bindparam("k")
    )

    with SessionLocal() as s:
        rows = s.execute(sql, {"qb": bitstr, "k": top_k}).mappings().all()
        return rows

def jaccard_search_hybrid(bitstr: str, qvec: List[float], top_k: int = 50, cand_k: int = 2000):
    sql = text("""
    WITH cand AS (
      SELECT c.id, c.smiles, c.fp_morgan_b2048_r2_bit
      FROM compound c
      WHERE c.fp_morgan_b2048_r2_bit IS NOT NULL
      ORDER BY c.fp_morgan_b2048_r2_vec <#> :qv   -- inner-product distance
      LIMIT :cand_k
    ),
    scored AS (
      SELECT
        c.id,
        c.smiles,
        CASE
          WHEN length(replace((c.fp_morgan_b2048_r2_bit | :qb)::text, '0','')) = 0
          THEN 1.0
          ELSE length(replace((c.fp_morgan_b2048_r2_bit & :qb)::text, '0',''))::float
               / NULLIF(length(replace((c.fp_morgan_b2048_r2_bit | :qb)::text, '0','')), 0)
        END AS jacc
      FROM cand c
    ),
    top_compounds AS (
      SELECT id, smiles, jacc
      FROM scored
      ORDER BY jacc DESC
      LIMIT :k
    )
    SELECT
      t.id,
      t.smiles,
      r.source,
      r.name,
      t.jacc
    FROM top_compounds t
    LEFT JOIN compound_record r
      ON r.compound_id = t.id
    WHERE r.name IS NOT NULL
    ORDER BY t.jacc DESC, r.source, r.ext_id
    """).bindparams(
        bindparam("qb", type_=BIT(2048)),
        bindparam("qv", type_=Vector(2048)),
        bindparam("k"),
        bindparam("cand_k"),
    )

    with SessionLocal() as s:
        return s.execute(sql, {"qb": bitstr, "qv": qvec, "k": top_k, "cand_k": cand_k}).mappings().all()

