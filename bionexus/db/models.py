from __future__ import annotations
from sqlalchemy import ARRAY, BigInteger, Float, SmallInteger, String, Text
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column
from sqlalchemy.dialects.postgresql import BIT
from pgvector.sqlalchemy import Vector

class Base(DeclarativeBase):
    pass

class Compound(Base):
    __tablename__ = "compound"
    id: Mapped[int] = mapped_column(BigInteger, primary_key=True)
    source: Mapped[str] = mapped_column(String(32), default="npatlas")
    ext_id: Mapped[str | None] = mapped_column(String(64))
    name: Mapped[str | None] = mapped_column(String(512))
    synonyms: Mapped[list[str] | None] = mapped_column(ARRAY(String))
    mol_formula: Mapped[str | None] = mapped_column(String(64))
    mol_weight: Mapped[float | None] = mapped_column(Float)
    exact_mass: Mapped[float | None] = mapped_column(Float)
    smiles: Mapped[str | None] = mapped_column(Text)
    inchi: Mapped[str | None] = mapped_column(Text)
    inchikey: Mapped[str | None] = mapped_column(String(27), index=True)
    m_plus_h: Mapped[float | None] = mapped_column(Float)
    m_plus_na: Mapped[float | None] = mapped_column(Float)
    
    # fingerprints
    # exact 2048-bit Morgan fingerprint (for Jaccard/Tanimoto)
    fp_morgan_b2048_r2_bit: Mapped[str | None] = mapped_column(BIT(2048))
    # cached popcount (# of 1-bits) for heuristic filters
    fp_morgan_b2048_r2_pop: Mapped[int | None] = mapped_column(SmallInteger)
    # binary-as-floats vector for fast ANN candidate search
    fp_morgan_b2048_r2_vec: Mapped[list[float] | None] = mapped_column(Vector(2048))
