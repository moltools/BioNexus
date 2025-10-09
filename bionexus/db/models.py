from __future__ import annotations
from sqlalchemy import ARRAY, BigInteger, Float, ForeignKey, SmallInteger, String, Text, UniqueConstraint
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship
from sqlalchemy.dialects.postgresql import BIT
from pgvector.sqlalchemy import Vector

class Base(DeclarativeBase):
    pass

class Compound(Base):
    __tablename__ = "compound"
    id: Mapped[int] = mapped_column(BigInteger, primary_key=True)

    # canonical identifiers/structure (one per molecule)
    inchikey: Mapped[str | None] = mapped_column(String(27), index=True)
    inchi: Mapped[str | None] = mapped_column(Text)
    smiles: Mapped[str | None] = mapped_column(Text)

    # canonical properties (computer/aggregated)
    mol_formula: Mapped[str | None] = mapped_column(String(64))
    mol_weight: Mapped[float | None] = mapped_column(Float)
    exact_mass: Mapped[float | None] = mapped_column(Float)
    m_plus_h: Mapped[float | None] = mapped_column(Float)
    m_plus_na: Mapped[float | None] = mapped_column(Float)
    
    # fingerprints
    fp_morgan_b2048_r2_bit: Mapped[str | None] = mapped_column(BIT(2048))
    fp_morgan_b2048_r2_pop: Mapped[int | None] = mapped_column(SmallInteger)
    fp_morgan_b2048_r2_vec: Mapped[list[float] | None] = mapped_column(Vector(2048))

    # backrefs
    records: Mapped[list["CompoundRecord"]] = relationship(
        back_populates="compound",
        cascade="all, delete-orphan"
    )

class CompoundRecord(Base):
    __tablename__ = "compound_record"
    id: Mapped[int] = mapped_column(BigInteger, primary_key=True)

    compound_id: Mapped[int] = mapped_column(ForeignKey("compound.id", ondelete="CASCADE"), index=True)

    source: Mapped[str] = mapped_column(String(32))
    ext_id: Mapped[str] = mapped_column(String(128))
    name: Mapped[str | None] = mapped_column(String(512))
    synonyms: Mapped[list[str] | None] = mapped_column(ARRAY(String))

    compound: Mapped["Compound"] = relationship(back_populates="records")

    # prevent duplicate mapping of the same external record to a compound
    __table_args__ = (UniqueConstraint("source", "ext_id", name="uq_compound_record_source_ext"),)
