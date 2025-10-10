from __future__ import annotations
from datetime import datetime
from sqlalchemy import (
    ARRAY, BigInteger, DateTime, Float, ForeignKey, SmallInteger, String, Text,
    UniqueConstraint, CheckConstraint, Index, text
)
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
    records: Mapped[list["CompoundRecord"]] = relationship(back_populates="compound", cascade="all, delete-orphan")

    __table_args__ = (UniqueConstraint("inchikey", name="uq_compound_inchikey"),)

class CompoundRecord(Base):
    __tablename__ = "compound_record"
    id: Mapped[int] = mapped_column(BigInteger, primary_key=True)

    compound_id: Mapped[int] = mapped_column(ForeignKey("compound.id", ondelete="CASCADE"), index=True)

    source: Mapped[str] = mapped_column(String(32))
    ext_id: Mapped[str] = mapped_column(String(128))
    name: Mapped[str | None] = mapped_column(String(512))
    synonyms: Mapped[list[str] | None] = mapped_column(ARRAY(String))

    compound: Mapped["Compound"] = relationship(back_populates="records")

    # allow same (source, ext_id) to link to multiple compounds
    __table_args__ = (UniqueConstraint("compound_id", "source", "ext_id", name="uq_compound_record_compound_source_ext"),)

class GenBankRegion(Base):
    __tablename__ = "genbank_region"

    id: Mapped[int] = mapped_column(BigInteger, primary_key=True)

    source: Mapped[str] = mapped_column(String(32))
    ext_id: Mapped[str] = mapped_column(String(128))

    gbk_text: Mapped[str] = mapped_column(Text)             # full GenBank as TEXT
    size_bytes: Mapped[int | None] = mapped_column()        # raw size on ingest
    sha256: Mapped[str | None] = mapped_column(String(64))  # content checksum

    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=text("NOW()"), nullable=False)
    updated_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=text("NOW()"), nullable=False)

    __table_args__ = (
        # uniqueness on (source, ext_id) (matches 0002_rev)
        UniqueConstraint("source", "ext_id", name="uq_genbank_region_source_ext"),
        # partial unique index on sha256 when present (matches 0002_rev)
        Index("ux_genbank_region_sha256_nonnull", "sha256", unique=True, postgresql_where=text("sha256 IS NOT NULL"),),
    )
