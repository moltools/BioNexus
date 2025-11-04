"""
Revision ID: 0003_rev
Revises: 0002_rev
Create Date: 2025-10-11 10:32:00.000000
"""

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB

from alembic import op

revision = "0003_rev"
down_revision = "0002_rev"
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        "annotation",
        sa.Column("id", sa.BigInteger, primary_key=True),
        # polymorphic target: either compound OR genbank_region (exactly one non-null)
        sa.Column(
            "compound_id",
            sa.BigInteger,
            sa.ForeignKey("compound.id", ondelete="CASCADE"),
            nullable=True,
        ),
        sa.Column(
            "genbank_region_id",
            sa.BigInteger,
            sa.ForeignKey("genbank_region.id", ondelete="CASCADE"),
            nullable=True,
        ),
        # flexible classification fields
        # examples:
        #   scheme: "taxonomy", key: "species", value: "Streptomyces coelicolor"
        #   scheme: "bioactivity", key: "antibacterial", value: "true"
        #   scheme: "npclassifier", key: "Superclass", "value": "Polyketides"
        #   scheme: "biosynthetic_class", key: "class", value: "NRPS"
        sa.Column("scheme", sa.String(64), nullable=False),
        sa.Column("key", sa.String(64), nullable=False),
        sa.Column("value", sa.String(256), nullable=False),
        # structured payload (for provenance, scores, ontology IDs, etc.)
        sa.Column("metadata_json", JSONB, nullable=True),
        sa.Column(
            "created_at",
            sa.DateTime(timezone=True),
            server_default=sa.text("NOW()"),
            nullable=False,
        ),
        sa.Column(
            "updated_at",
            sa.DateTime(timezone=True),
            server_default=sa.text("NOW()"),
            nullable=False,
        ),
        sa.CheckConstraint(
            "("
            "(compound_id IS NOT NULL AND genbank_region_id IS NULL) OR "
            "(compound_id IS NULL AND genbank_region_id IS NOT NULL)"
            ")",
            name="ck_annotation_exactly_one_target",
        ),
    )

    # avoid duplicate facts for a given target and label triplet
    op.create_unique_constraint(
        "uq_annotation_target_scheme_key_value",
        "annotation",
        ["compound_id", "genbank_region_id", "scheme", "key", "value"],
    )

    # helpful indexes for common lookups
    op.create_index("ix_annotation_compound_id", "annotation", ["compound_id"])
    op.create_index("ix_annotation_genbank_region_id", "annotation", ["genbank_region_id"])
    op.create_index("ix_annotation_scheme_key_value", "annotation", ["scheme", "key", "value"])

    # JSON metadata indexing (existence/path queries)
    op.execute("""
    CREATE INDEX IF NOT EXISTS ix_annotation_metadata_json_gin
    ON annotation USING GIN (metadata_json);         
    """)

    # timestamp trigger
    op.execute("""
    CREATE OR REPLACE FUNCTION public.set_timestamp_annotation()
    RETURNS TRIGGER AS $$
    BEGIN
      NEW.updated_at = NOW();
      RETURN NEW;
    END;
    $$ LANGUAGE plpgsql;
    """)
    op.execute("""
    CREATE TRIGGER annotation_set_timestamp
    BEFORE UPDATE ON public.annotation
    FOR EACH ROW
    EXECUTE FUNCTION public.set_timestamp_annotation();
    """)


def downgrade():
    # drop triggers and their functions first
    op.execute("DROP TRIGGER IF EXISTS annotation_set_timestamp ON public.annotation;")
    op.execute("DROP FUNCTION IF EXISTS public.set_timestamp_annotation;")

    # drop indexes created via SQL
    op.execute("DROP INDEX IF EXISTS ix_annotation_metadata_json_gin;")

    # drop contraints/indexes created via Alembic helpers
    op.drop_index("ix_annotation_scheme_key_value", table_name="annotation")
    op.drop_index("ix_annotation_genbank_region_id", table_name="annotation")
    op.drop_index("ix_annotation_compound_id", table_name="annotation")
    op.drop_constraint("uq_annotation_target_scheme_key_value", "annotation", type_="unique")

    # drop table
    op.drop_table("annotation")
