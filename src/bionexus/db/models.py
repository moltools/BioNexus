"""Database models for BioNexus."""

from __future__ import annotations

from datetime import datetime
from typing import Any

from pgvector.sqlalchemy import Vector
from sqlalchemy import (
    ARRAY,
    BigInteger,
    CheckConstraint,
    DateTime,
    Float,
    ForeignKey,
    Index,
    Integer,
    SmallInteger,
    String,
    Text,
    UniqueConstraint,
    text,
)
from sqlalchemy.dialects.postgresql import BIT, JSONB
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


class Base(DeclarativeBase):
    """
    Base class for all ORM models.
    """

    pass


class Compound(Base):
    """
    Represents a chemical compound with various properties and identifiers.

    """

    __tablename__ = "compound"

    id: Mapped[int] = mapped_column(BigInteger, primary_key=True)

    # Canonical identifiers/structure (one per molecule)
    inchikey: Mapped[str | None] = mapped_column(String(27), index=True)
    inchi: Mapped[str | None] = mapped_column(Text)
    smiles: Mapped[str | None] = mapped_column(Text)

    # Canonical properties (computer/aggregated)
    mol_formula: Mapped[str | None] = mapped_column(String(64))
    mol_weight: Mapped[float | None] = mapped_column(Float)
    exact_mass: Mapped[float | None] = mapped_column(Float)
    m_plus_h: Mapped[float | None] = mapped_column(Float)
    m_plus_na: Mapped[float | None] = mapped_column(Float)

    # Elemental counts
    c_count: Mapped[int | None] = mapped_column(Integer)
    h_count: Mapped[int | None] = mapped_column(Integer)
    n_count: Mapped[int | None] = mapped_column(Integer)
    o_count: Mapped[int | None] = mapped_column(Integer)
    s_count: Mapped[int | None] = mapped_column(Integer)
    p_count: Mapped[int | None] = mapped_column(Integer)
    f_count: Mapped[int | None] = mapped_column(Integer)
    cl_count: Mapped[int | None] = mapped_column(Integer)
    br_count: Mapped[int | None] = mapped_column(Integer)
    i_count: Mapped[int | None] = mapped_column(Integer)

    # Fingerprints
    fp_morgan_b2048_r2_bit: Mapped[str | None] = mapped_column(BIT(2048))
    fp_morgan_b2048_r2_pop: Mapped[int | None] = mapped_column(SmallInteger)
    fp_morgan_b2048_r2_vec: Mapped[list[float] | None] = mapped_column(Vector(2048))

    # Backrefs:
    # allows to navigate from compound to records like:
    #
    #   comp = s.get(Compound, 123)
    #   for rec in comp.records:
    #       print(rec.source, rec.ext_id)
    #
    records: Mapped[list[CompoundRecord]] = relationship(back_populates="compound", cascade="all, delete-orphan")
    # allows to navigate from compound to annotations like:
    #
    #   comp = s.get(Compound, 123)
    #   for ann in comp.annotations:
    #       print(ann.scheme, ann.key, ann.value)
    #
    annotations: Mapped[list[Annotation]] = relationship(
        back_populates="compound", passive_deletes=True
    )  # no delete-orphan because FK is nullable by design

    __table_args__ = (UniqueConstraint("inchikey", name="uq_compound_inchikey"),)


class CompoundRecord(Base):
    """
    Represents an external record or identifier for a compound.
    """

    __tablename__ = "compound_record"
    id: Mapped[int] = mapped_column(BigInteger, primary_key=True)

    compound_id: Mapped[int] = mapped_column(ForeignKey("compound.id", ondelete="CASCADE"), index=True)

    source: Mapped[str] = mapped_column(String(32))
    ext_id: Mapped[str] = mapped_column(String(128))
    name: Mapped[str | None] = mapped_column(String(512))
    synonyms: Mapped[list[str] | None] = mapped_column(ARRAY(String))

    compound: Mapped[Compound] = relationship(back_populates="records")

    # Allow same (source, ext_id) to link to multiple compounds
    __table_args__ = (
        UniqueConstraint(
            "compound_id",
            "source",
            "ext_id",
            name="uq_compound_record_compound_source_ext",
        ),
    )


class GenBankRegion(Base):
    """
    Represents a GenBank region with associated annotations.
    """

    __tablename__ = "genbank_region"

    id: Mapped[int] = mapped_column(BigInteger, primary_key=True)

    source: Mapped[str] = mapped_column(String(32))
    ext_id: Mapped[str] = mapped_column(String(128))

    gbk_text: Mapped[str] = mapped_column(Text)  # full GenBank as TEXT
    size_bytes: Mapped[int | None] = mapped_column()  # raw size on ingest
    sha256: Mapped[str | None] = mapped_column(String(64))  # content checksum

    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=text("NOW()"), nullable=False)
    updated_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=text("NOW()"), nullable=False)

    # Backrefs:
    # allows to navigate from region to annotations like:
    #
    #   region = s.get(GenBankRegion, 123)
    #   for ann in region.annotations:
    #       print(ann.scheme, ann.key, ann.value)
    #
    annotations: Mapped[list[Annotation]] = relationship(
        back_populates="genbank_region", passive_deletes=True
    )  # no delete-orphan because FK is nullable by design

    __table_args__ = (
        # uniqueness on (source, ext_id) (matches 0002_rev)
        UniqueConstraint("source", "ext_id", name="uq_genbank_region_source_ext"),
        # partial unique index on sha256 when present (matches 0002_rev)
        Index(
            "ux_genbank_region_sha256_nonnull",
            "sha256",
            unique=True,
            postgresql_where=text("sha256 IS NOT NULL"),
        ),
    )


class Annotation(Base):
    """
    Represents an annotation that can be linked to either a Compound or a GenBankRegion.
    """

    __tablename__ = "annotation"

    id: Mapped[int] = mapped_column(BigInteger, primary_key=True)

    # Exactly one of these must be non-null (enforced via CHECK)
    compound_id: Mapped[int | None] = mapped_column(ForeignKey("compound.id", ondelete="CASCADE"), index=True)
    genbank_region_id: Mapped[int | None] = mapped_column(
        ForeignKey("genbank_region.id", ondelete="CASCADE"), index=True
    )

    # Flexible label triplet
    scheme: Mapped[str] = mapped_column(String(64), nullable=False)
    key: Mapped[str] = mapped_column(String(64), nullable=False)
    value: Mapped[str] = mapped_column(String(256), nullable=False)

    # Optional structured payload
    metadata_json: Mapped[dict[str, Any] | None] = mapped_column(JSONB, nullable=True)

    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=text("NOW()"), nullable=False)
    updated_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=text("NOW()"), nullable=False)

    # Relationship (no delete-orphan because FK is nullable be design)
    compound: Mapped[Compound | None] = relationship(back_populates="annotations")
    genbank_region: Mapped[GenBankRegion | None] = relationship(back_populates="annotations")

    __table_args__ = (
        # Enforce exactly one target
        CheckConstraint(
            "("
            "(compound_id IS NOT NULL AND genbank_region_id IS NULL) OR "
            "(compound_id IS NULL AND genbank_region_id IS NOT NULL)"
            ")",
            name="ck_annotation_exactly_one_target",
        ),
        # Prevent duplicate facts per target
        UniqueConstraint(
            "compound_id",
            "genbank_region_id",
            "scheme",
            "key",
            "value",
            name="uq_annotation_target_scheme_key_value",
        ),
        # Common lookup indexes
        Index("ix_annotation_scheme_key_value", "scheme", "key", "value"),
        # JSONB GIN index (same as in 0003_rev)
        Index("ix_annotation_metadata_json", metadata_json, postgresql_using="gin"),
    )


class Ruleset(Base):
    """
    Represents a set of rules used for retrosynthetic analysis.
    """

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

    # Backref to results
    results: Mapped[list[RetroMolCompound]] = relationship(back_populates="ruleset", cascade="all, delete-orphan")

    __table_args__ = (Index("ix_ruleset_version", "version", unique=True),)


class RetroMolCompound(Base):
    """
    Represents the results of retrosynthetic analysis for a compound using a specific ruleset.
    """

    __tablename__ = "retromol_compound"

    id: Mapped[int] = mapped_column(BigInteger, primary_key=True)
    compound_id: Mapped[int] = mapped_column(ForeignKey("compound.id", ondelete="CASCADE"), nullable=False, index=True)
    ruleset_id: Mapped[int] = mapped_column(ForeignKey("ruleset.id", ondelete="RESTRICT"), nullable=False, index=True)

    result_json: Mapped[dict[str, Any]] = mapped_column(JSONB, nullable=False)
    
    # Precomputed properties from result_json
    coverage: Mapped[float | None] = mapped_column(Float)

    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=text("NOW()"), nullable=False)
    updated_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=text("NOW()"), nullable=False)

    compound: Mapped[Compound] = relationship(back_populates="retromol_results")
    ruleset: Mapped[Ruleset] = relationship(back_populates="results")

    retrofingerprints: Mapped[list["RetroFingerprint"]] = relationship(
        back_populates="retromol_compound",
        passive_deletes=True,
    )  # no delete-orphan because FK is nullable by design

    __table_args__ = (
        Index("ix_retromol_compound_compound_id", "compound_id"),
        Index("ix_retromol_compound_ruleset_id", "ruleset_id"),
        # Only one result per (compound, ruleset)
        UniqueConstraint(
            "compound_id",
            "ruleset_id",
            name="uq_retromol_compound_compound_ruleset",
        ),
    )


# Add the backref to Compound
Compound.retromol_results = relationship(
    "RetroMolCompound",
    back_populates="compound",
    cascade="all, delete-orphan",
)


class BioCrackerGenBank(Base):
    """
    Represents the readouts for a GenBank region using BioCracker
    """

    __tablename__ = "biocracker_genbank"

    id: Mapped[int] = mapped_column(BigInteger, primary_key=True)
    genbank_region_id: Mapped[int] = mapped_column(ForeignKey("genbank_region.id", ondelete="CASCADE"), nullable=False, index=True)

    result_json: Mapped[dict[str, Any]] = mapped_column(JSONB, nullable=False)

    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=text("NOW()"), nullable=False)
    updated_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=text("NOW()"), nullable=False)

    retrofingerprints: Mapped[list["RetroFingerprint"]] = relationship(
        back_populates="biocracker_genbank",
        passive_deletes=True,
    )  # no delete-orphan because FK is nullable by design

    __table_args__ = (
        Index("ix_biocracker_genbank_genbank_region_id", "genbank_region_id"),
        # Only one result per (genbank_region)
        UniqueConstraint(
            "genbank_region_id",
            name="uq_biocracker_genbank_genbank_region",
        ),
    )


class RetroFingerprint(Base):
    """
    Represents RetroMol biosynthetic fingerprints for compounds and genbank records.
    """

    __tablename__ = "retrofingerprint"

    id: Mapped[int] =  mapped_column(BigInteger, primary_key=True)

    # Can be linked to either RetroMolCompound or BioCrackerGenBank; multiple RetroFingerprints per target allowed
    retromol_compound_id: Mapped[int | None] = mapped_column(ForeignKey("retromol_compound.id", ondelete="CASCADE"), index=True)
    biocracker_genbank_id: Mapped[int | None] = mapped_column(ForeignKey("biocracker_genbank.id", ondelete="CASCADE"), index=True) 

    # Fingerprints
    fp_retro_b512_bit: Mapped[str] = mapped_column(BIT(512))
    fp_retro_b512_pop: Mapped[int] = mapped_column(SmallInteger)
    fp_retro_b512_vec_binary: Mapped[list[float]] = mapped_column(Vector(512))
    fp_retro_b512_vec_counted: Mapped[list[float]] = mapped_column(Vector(512))

    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=text("NOW()"), nullable=False)
    updated_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=text("NOW()"), nullable=False)

    # Backrefs
    retromol_compound: Mapped["RetroMolCompound | None"] = relationship(back_populates="retrofingerprints")
    biocracker_genbank: Mapped["BioCrackerGenBank | None"] = relationship(back_populates="retrofingerprints")

    __table_args__ = (
        # Enforce exactly one target
        CheckConstraint(
            "("
            "(retromol_compound_id IS NOT NULL AND biocracker_genbank_id IS NULL) OR "
            "(retromol_compound_id IS NULL AND biocracker_genbank_id IS NOT NULL)"
            ")",
            name="ck_retrofingerprint_exactly_one_target",
        ),
    )