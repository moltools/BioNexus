"""
Revision ID: 0004_rev
Revises: 0003_rev
Create Date: 2025-10-14 12:25:00.000000
"""

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import BIT, JSONB

revision = "0004_rev"
down_revision = "0003_rev"
branch_labels = None
depends_on = None


def upgrade():
    # for digest() to compute SHA-256 in SQL
    op.execute("CREATE EXTENSION IF NOT EXISTS pgcrypto;")

    # --- RULESET (immutable, versioned by DB trigger) -------------------------
    # add table that stores all rule YAMLS with version, and associate version with result
    op.create_table(
        "ruleset",
        sa.Column("id", sa.BigInteger, primary_key=True),
        sa.Column(
            "version", sa.Integer, nullable=False, unique=True
        ),  # auto-assign by trigger
        # store rules as text (YAML) for both monomers and reaction rules, also store sha256 for integrity
        sa.Column("matching_rules_yaml", sa.Text, nullable=False),
        sa.Column("matching_rules_sha256", sa.String(64), nullable=False, unique=True),
        sa.Column("reaction_rules_yaml", sa.Text, nullable=False),
        sa.Column("reaction_rules_sha256", sa.String(64), nullable=False, unique=True),
        # combined content hash (enforces uniqueness of the pair)
        sa.Column("ruleset_sha256", sa.String(64), nullable=False, unique=True),
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

    # timestamps (auto-update updated_at)
    op.execute("""
    CREATE OR REPLACE FUNCTION public.set_timestamp_ruleset()
    RETURNS TRIGGER AS $$
    BEGIN
      NEW.updated_at = NOW();
      RETURN NEW;
    END;
    $$ LANGUAGE plpgsql;
    """)

    # immutable rows: reject UPDATE on ruleset
    op.execute("""
    CREATE OR REPLACE FUNCTION public.prevent_update_ruleset()
    RETURNS TRIGGER AS $$
    BEGIN
      RAISE EXCEPTION 'ruleset rows are immutable (use INSERT to add a new version)';
    END;
    $$ LANGUAGE plpgsql;
    """)

    # BEFORE INSERT: compute combined SHA, assign next version atomically
    # Also sanity-check provided per-section SHA256 strings
    op.execute("""
    CREATE OR REPLACE FUNCTION public.ruleset_before_insert()
    RETURNS TRIGGER AS $$
    DECLARE
    computed_matching_sha text;
    computed_reaction_sha text;
    computed_ruleset_sha  text;
    next_ver integer;
    BEGIN
    computed_matching_sha := encode(digest(coalesce(NEW.matching_rules_yaml, ''), 'sha256'), 'hex');
    computed_reaction_sha := encode(digest(coalesce(NEW.reaction_rules_yaml, ''), 'sha256'), 'hex');

    IF NEW.matching_rules_sha256 IS NULL OR NEW.reaction_rules_sha256 IS NULL THEN
        RAISE EXCEPTION 'matching_rules_sha256 and reaction_rules_sha256 must be provided';
    END IF;

    IF lower(NEW.matching_rules_sha256) <> lower(computed_matching_sha) THEN
        RAISE EXCEPTION 'Provided matching_rules_sha256 (%) does not match computed (%)',
        NEW.matching_rules_sha256, computed_matching_sha;
    END IF;

    IF lower(NEW.reaction_rules_sha256) <> lower(computed_reaction_sha) THEN
        RAISE EXCEPTION 'Provided reaction_rules_sha256 (%) does not match computed (%)',
        NEW.reaction_rules_sha256, computed_reaction_sha;
    END IF;

    computed_ruleset_sha := encode(
        digest(coalesce(NEW.matching_rules_yaml,'') || '||' || coalesce(NEW.reaction_rules_yaml,''), 'sha256'),
        'hex'
    );
    NEW.ruleset_sha256 := computed_ruleset_sha;

    IF EXISTS (
        SELECT 1 FROM public.ruleset r WHERE lower(r.ruleset_sha256) = lower(NEW.ruleset_sha256)
    ) THEN
        RAISE EXCEPTION 'Identical ruleset already exists (ruleset_sha256=%). Insert skipped.', NEW.ruleset_sha256;
    END IF;

    SELECT COALESCE(MAX(version), 0) + 1 INTO next_ver FROM public.ruleset;
    NEW.version := next_ver;

    RETURN NEW;
    END;
    $$ LANGUAGE plpgsql;
    """)

    # triggers for ruleset
    op.execute("""
    CREATE TRIGGER ruleset_set_timestamp
    BEFORE UPDATE ON public.ruleset
    FOR EACH ROW
    EXECUTE FUNCTION public.set_timestamp_ruleset();
    """)
    op.execute("""
    CREATE TRIGGER ruleset_prevent_update
    BEFORE UPDATE ON public.ruleset
    FOR EACH ROW
    EXECUTE FUNCTION public.prevent_update_ruleset();
    """)
    op.execute("""
    CREATE TRIGGER ruleset_before_insert
    BEFORE INSERT ON public.ruleset
    FOR EACH ROW
    EXECUTE FUNCTION public.ruleset_before_insert();
    """)

    # --- RETROMOL results (multiple per compound) -----------------------------
    # one compound can have multiple RetroMol entries associated with it
    # (e.g, different versions, different parameters)
    op.create_table(
        "retromol_compound",
        sa.Column("id", sa.BigInteger, primary_key=True),
        sa.Column(
            "compound_id",
            sa.BigInteger,
            sa.ForeignKey("compound.id", ondelete="CASCADE"),
            nullable=False,
        ),
        # RESTRICT: if ruleset is in use, don't allow deletion
        sa.Column(
            "ruleset_id",
            sa.BigInteger,
            sa.ForeignKey("ruleset.id", ondelete="RESTRICT"),
            nullable=False,
        ),
        # structured payload (for provenance, scores, parameters, etc.)
        sa.Column("result_json", JSONB, nullable=False),
        sa.Column("fp_retro", BIT(512), nullable=True),
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

    # speed lookups by compound or ruleset
    op.create_index(
        "ix_retromol_compound_compound_id",
        "retromol_compound",
        ["compound_id"],
        unique=False,
    )
    op.create_index(
        "ix_retromol_compound_ruleset_id",
        "retromol_compound",
        ["ruleset_id"],
        unique=False,
    )

    # tmestamp tigger (table-scoped function name to avoid collisions)
    op.execute("""
    CREATE OR REPLACE FUNCTION public.set_timestamp_retromol_compound()
    RETURNS TRIGGER AS $$
    BEGIN
      NEW.updated_at = NOW();  -- use clock_timestamp() if you want wall-clock per statement
      RETURN NEW;
    END;
    $$ LANGUAGE plpgsql;
    """)
    op.execute("""
    CREATE TRIGGER retromol_compound_set_timestamp
    BEFORE UPDATE ON public.retromol_compound
    FOR EACH ROW
    EXECUTE FUNCTION public.set_timestamp_retromol_compound();
    """)


def downgrade():
    # drop child triggers/indexes/tables first
    op.execute(
        "DROP TRIGGER IF EXISTS retromol_compound_set_timestamp ON public.retromol_compound;"
    )
    op.execute("DROP FUNCTION IF EXISTS public.set_timestamp_retromol_compound;")
    op.drop_index("ix_retromol_compound_ruleset_id", table_name="retromol_compound")
    op.drop_index("ix_retromol_compound_compound_id", table_name="retromol_compound")
    op.drop_table("retromol_compound")

    # drop ruleset triggers & functions, then table
    op.execute("DROP TRIGGER IF EXISTS ruleset_before_insert ON public.ruleset;")
    op.execute("DROP TRIGGER IF EXISTS ruleset_prevent_update ON public.ruleset;")
    op.execute("DROP TRIGGER IF EXISTS ruleset_set_timestamp ON public.ruleset;")
    op.execute("DROP FUNCTION IF EXISTS public.ruleset_before_insert;")
    op.execute("DROP FUNCTION IF EXISTS public.prevent_update_ruleset;")
    op.execute("DROP FUNCTION IF EXISTS public.set_timestamp_ruleset;")
    op.drop_table("ruleset")
