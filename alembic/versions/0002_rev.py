"""
Revision ID: 0002_rev
Revises: 0001_init
Create Date: 2025-10-10 21:56:15.797293
"""
import os
from alembic import op
import sqlalchemy as sa

revision = "0002_rev"
down_revision = "0001_init"
branch_labels = None
depends_on = None

def upgrade():
    op.create_table(
        "genbank_region",
        sa.Column("id", sa.BigInteger, primary_key=True),
    
        sa.Column("source", sa.String(32), nullable=False),
        sa.Column("ext_id", sa.String(128), nullable=False),

        sa.Column("gbk_text", sa.Text, nullable=False),         # full GenBank as TEXT
        sa.Column("size_bytes", sa.Integer, nullable=True),     # raw size on ingest
        sa.Column("sha256", sa.String(64), nullable=True),      # content checksum for dedup/integrity

        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.text("NOW()"), nullable=False),
        sa.Column("updated_at", sa.DateTime(timezone=True), server_default=sa.text("NOW()"), nullable=False),
    )

    # uniqueness on (source, ext_id) — this also creates the btree index you need
    op.create_unique_constraint("uq_genbank_region_source_ext", "genbank_region", ["source", "ext_id"],)

    # fast de-dup by content when sha256 is present (partial unique index)
    op.execute("""
    CREATE UNIQUE INDEX IF NOT EXISTS ux_genbank_region_sha256_nonnull
    ON genbank_region(sha256)
    WHERE sha256 IS NOT NULL;
    """)

    # timestamp trigger (table-scoped function name to avoid collisions)
    op.execute("""
    CREATE OR REPLACE FUNCTION public.set_timestamp_genbank_region()
    RETURNS TRIGGER AS $$
    BEGIN
      NEW.updated_at = NOW();  -- use clock_timestamp() if you want wall-clock per statement
      RETURN NEW;
    END;
    $$ LANGUAGE plpgsql;
    """)
    op.execute("""
    CREATE TRIGGER genbank_region_set_timestamp
    BEFORE UPDATE ON public.genbank_region
    FOR EACH ROW
    EXECUTE FUNCTION public.set_timestamp_genbank_region();
    """)

def downgrade():
    # drop trigger & its function first
    op.execute("DROP TRIGGER IF EXISTS genbank_region_set_timestamp ON public.genbank_region;")
    op.execute("DROP FUNCTION IF EXISTS public.set_timestamp_genbank_region;")

    # drop optional partial unique index (if created)
    op.execute("DROP INDEX IF EXISTS ux_genbank_region_sha256_nonnull;")

    op.drop_constraint("uq_genbank_region_source_ext", "genbank_region", type_="unique")
    op.drop_table("genbank_region")
