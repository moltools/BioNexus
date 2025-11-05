import sqlalchemy as sa
from alembic import op
from pgvector.sqlalchemy import Vector
from sqlalchemy.dialects.postgresql import BIT, JSONB

revision = "0007_rev"
down_revision = "0006_rev"
branch_labels = None
depends_on = None


def upgrade():
    op.execute("CREATE EXTENSION IF NOT EXISTS vector;")
    op.create_table(
        "biocracker_genbank",
        sa.Column("id", sa.BigInteger, primary_key=True),
        sa.Column(
            "genbank_region_id",
            sa.BigInteger,
            sa.ForeignKey("genbank_region.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column("result_json", JSONB, nullable=False),
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
        sa.UniqueConstraint(
            "genbank_region_id",
            name="uq_biocracker_genbank_genbank_region"
        ),
    )

    # Explicit index for clarity
    op.create_index(
        "ix_biocracker_genbank_genbank_region_id",
        "biocracker_genbank",
        ["genbank_region_id"]
    )

    op.create_table(
        "retrofingerprint",
        sa.Column("id", sa.BigInteger, primary_key=True),
        sa.Column(
            "retromol_compound_id",
            sa.BigInteger,
            sa.ForeignKey("retromol_compound.id", ondelete="CASCADE"),
            nullable=True,
        ),
        sa.Column(
            "biocracker_genbank_id",
            sa.BigInteger,
            sa.ForeignKey("biocracker_genbank.id", ondelete="CASCADE"),
            nullable=True,
        ),
        sa.Column("fp_retro_b512_bit", BIT(512), nullable=False),
        sa.Column("fp_retro_b512_pop", sa.SmallInteger, nullable=False),
        sa.Column("fp_retro_b512_vec_binary", Vector(512), nullable=False),
        sa.Column("fp_retro_b512_vec_counted", Vector(512), nullable=False),
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
    )

    # Explicit indexes for clarity
    op.create_index(
        "ix_retrofingerprint_retromol_compound_id",
        "retrofingerprint",
        ["retromol_compound_id"],
    )
    op.create_index(
        "ix_retrofingerprint_biocracker_genbank_id",
        "retrofingerprint",
        ["biocracker_genbank_id"],
    )

    # Timestamp trigger (table-scoped function name to avoid collisions)
    op.execute("""
    CREATE OR REPLACE FUNCTION public.set_timestamp_biocracker_genbank()
    RETURNS TRIGGER AS $$
    BEGIN
      NEW.updated_at = NOW();  -- use clock_timestamp() if you want wall-clock per statement
      RETURN NEW;
    END;
    $$ LANGUAGE plpgsql;
    """)
    op.execute("""
    CREATE TRIGGER biocracker_genbank_set_timestamp
    BEFORE UPDATE ON public.biocracker_genbank
    FOR EACH ROW
    EXECUTE FUNCTION public.set_timestamp_biocracker_genbank();
    """)

    op.execute("""
    CREATE OR REPLACE FUNCTION public.set_timestamp_retrofingerprint()
    RETURNS TRIGGER AS $$
    BEGIN
      NEW.updated_at = NOW();  -- use clock_timestamp() if you want wall-clock per statement
      RETURN NEW;
    END;
    $$ LANGUAGE plpgsql;
    """)
    op.execute("""
    CREATE TRIGGER retrofingerprint_set_timestamp
    BEFORE UPDATE ON public.retrofingerprint
    FOR EACH ROW
    EXECUTE FUNCTION public.set_timestamp_retrofingerprint();
    """)


def downgrade():
    # drop trigger and its function first
    op.execute("DROP TRIGGER IF EXISTS retrofingerprint_set_timestamp ON public.retrofingerprint;")
    op.execute("DROP FUNCTION IF EXISTS public.set_timestamp_retrofingerprint;")

    op.execute("DROP TRIGGER IF EXISTS biocracker_genbank_set_timestamp ON public.biocracker_genbank;")
    op.execute("DROP FUNCTION IF EXISTS public.set_timestamp_biocracker_genbank;")

    op.drop_index("ix_retrofingerprint_biocracker_genbank_id", table_name="retrofingerprint")
    op.drop_index("ix_retrofingerprint_retromol_compound_id", table_name="retrofingerprint")
    op.drop_table("retrofingerprint")

    op.drop_index("ix_biocracker_genbank_genbank_region_id", table_name="biocracker_genbank")
    op.drop_table("biocracker_genbank")
