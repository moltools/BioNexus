from __future__ import annotations
from datetime import datetime
from typing import Any
from sqlalchemy import (
    ARRAY, Integer, BigInteger, DateTime, Float, ForeignKey, SmallInteger, String, Text,
    CheckConstraint, UniqueConstraint, Index, text
)
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship
from sqlalchemy.dialects.postgresql import BIT, JSONB
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
    # allows to navigate from compound to records like:
    # 
    #   comp = s.get(Compound, 123)
    #   for rec in comp.records:
    #       print(rec.source, rec.ext_id)
    #
    records: Mapped[list["CompoundRecord"]] = relationship(back_populates="compound", cascade="all, delete-orphan")
    # allows to navigate from compound to annotations like:
    # 
    #   comp = s.get(Compound, 123)
    #   for ann in comp.annotations:
    #       print(ann.scheme, ann.key, ann.value)
    #
    annotations: Mapped[list["Annotation"]] = relationship(back_populates="compound", passive_deletes=True)  # no delete-orphan because FK is nullable by design

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

    # backrefs
    # allows to navigate from region to annotations like:
    # 
    #   region = s.get(GenBankRegion, 123)
    #   for ann in region.annotations:
    #       print(ann.scheme, ann.key, ann.value)
    #
    annotations: Mapped[list["Annotation"]] = relationship(back_populates="genbank_region", passive_deletes=True)  # no delete-orphan because FK is nullable by design

    __table_args__ = (
        # uniqueness on (source, ext_id) (matches 0002_rev)
        UniqueConstraint("source", "ext_id", name="uq_genbank_region_source_ext"),
        # partial unique index on sha256 when present (matches 0002_rev)
        Index("ux_genbank_region_sha256_nonnull", "sha256", unique=True, postgresql_where=text("sha256 IS NOT NULL"),),
    )

class Annotation(Base):
    __tablename__ = "annotation"

    id: Mapped[int] = mapped_column(BigInteger, primary_key=True)

    # exactly one of these must be non-null (enforced via CHECK)
    compound_id: Mapped[int | None] = mapped_column(ForeignKey("compound.id", ondelete="CASCADE"), index=True)
    genbank_region_id: Mapped[int | None] = mapped_column(ForeignKey("genbank_region.id", ondelete="CASCADE"), index=True)

    # flexible label triplet
    scheme: Mapped[str] = mapped_column(String(64), nullable=False)
    key: Mapped[str] = mapped_column(String(64), nullable=False)
    value: Mapped[str] = mapped_column(String(256), nullable=False)

    # optional structured payload
    metadata_json: Mapped[dict[str, Any] | None] = mapped_column(JSONB, nullable=True)

    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=text("NOW()"), nullable=False)
    updated_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=text("NOW()"), nullable=False)

    # relationship (no delete-orphan because FK is nullable be design)
    compound: Mapped["Compound | None"] = relationship(back_populates="annotations")
    genbank_region: Mapped["GenBankRegion | None"] = relationship(back_populates="annotations")

    __table_args__ = (
        # enforce exactly one target
        CheckConstraint(
            "("
            "(compound_id IS NOT NULL AND genbank_region_id IS NULL) OR "
            "(compound_id IS NULL AND genbank_region_id IS NOT NULL)"
            ")",
            name="ck_annotation_exactly_one_target",
        ),
        # prevent duplicate facts per target
        UniqueConstraint(
            "compound_id", "genbank_region_id", "scheme", "key", "value",
            name="uq_annotation_target_scheme_key_value",
        ),
        # common lookup indexes
        Index("ix_annotation_scheme_key_value", "scheme", "key", "value"),
        # JSONB GIN index (same as in 0003_rev)
        Index("ix_annotation_metadata_json", metadata_json, postgresql_using="gin"),
    )

class Ruleset(Base):
    __tablename__ = "ruleset"

    id: Mapped[int] = mapped_column(BigInteger, primary_key=True)
    version: Mapped[int] = mapped_column(Integer, nullable=False, unique=True)

    matching_rules_yaml: Mapped[str] = mapped_column(Text, nullable=False)
    matching_rules_sha256: Mapped[str] = mapped_column(String(64), nullable=False, unique=True)
    reaction_rules_yaml: Mapped[str] = mapped_column(Text, nullable=False)
    reaction_rules_sha256: Mapped[str] = mapped_column(String(64), nullable=False, unique=True)
    ruleset_sha256: Mapped[str] = mapped_column(String(64), nullable=False, unique=True)

    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=text("NOW()"), nullable=False)
    updated_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=text("NOW()"), nullable=False)

    # backref to results
    results: Mapped[list["RetroMolCompound"]] = relationship(back_populates="ruleset", cascade="all, delete-orphan")

    __table_args__ = (Index("ix_ruleset_version", "version", unique=True),)

class RetroMolCompound(Base):
    __tablename__ = "retromol_compound"

    id: Mapped[int] = mapped_column(BigInteger, primary_key=True)
    compound_id: Mapped[int] = mapped_column(ForeignKey("compound.id", ondelete="CASCADE"), nullable=False, index=True)
    ruleset_id: Mapped[int] = mapped_column(ForeignKey("ruleset.id", ondelete="RESTRICT"), nullable=False, index=True)

    result_json: Mapped[dict[str, Any]] = mapped_column(JSONB, nullable=False)
    fp_retro: Mapped[str | None] = mapped_column(BIT(2048))

    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=text("NOW()"), nullable=False)
    updated_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=text("NOW()"), nullable=False)

    compound: Mapped["Compound"] = relationship(back_populates="retromol_results")
    ruleset: Mapped["Ruleset"] = relationship(back_populates="results")

    __table_args__ = (
        Index("ix_retromol_compound_compound_id", "compound_id"),
        Index("ix_retromol_compound_ruleset_id", "ruleset_id"),
        # only one result per (compound, ruleset)
        UniqueConstraint(
            "compound_id",
            "ruleset_id",
            name="uq_retromol_compound_compound_ruleset",
        ),
    )

# add the backref to Compound
Compound.retromol_results = relationship(
    "RetroMolCompound",
    back_populates="compound",
    cascade="all, delete-orphan",
)
