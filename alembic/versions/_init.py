from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import BIT
from pgvector.sqlalchemy import Vector

revision = "_init"
down_revision = None  # initial migration, no prior revision
branch_labels = None
depends_on = None

def upgrade():
    op.execute("CREATE EXTENSION IF NOT EXISTS vector;")
    op.create_table(
        "compound",
        sa.Column("id", sa.BigInteger, primary_key=True),
        sa.Column("source", sa.String(32), nullable=False, server_default="npatlas"),
        sa.Column("ext_id", sa.String(64)),
        sa.Column("name", sa.String(512)),
        sa.Column("synonyms", sa.ARRAY(sa.String)),
        sa.Column("mol_formula", sa.String(64)),
        sa.Column("mol_weight", sa.Float),
        sa.Column("exact_mass", sa.Float),
        sa.Column("smiles", sa.Text),
        sa.Column("inchi", sa.Text),
        sa.Column("inchikey", sa.String(27)),
        sa.Column("m_plus_h", sa.Float),
        sa.Column("m_plus_na", sa.Float),

        # fingerprints
        sa.Column("fp_morgan_b2048_r2_bit", BIT(2048)),
        sa.Column("fp_morgan_b2048_r2_pop", sa.SmallInteger),
        sa.Column("fp_morgan_b2048_r2_vec", Vector(2048)),
    )
    op.create_index("ix_compound_inchikey", "compound", ["inchikey"])
    op.create_unique_constraint("uq_compound_inchikey", "compound", ["inchikey"])

def downgrade():
    op.drop_constraint("uq_compound_inchikey", "compound", type_="unique")
    op.drop_index("ix_compound_inchikey", table_name="compound")
    op.drop_table("compound")
