# alembic/versions/_init.py
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import BIT
from pgvector.sqlalchemy import Vector

revision = "_init"
down_revision = None
branch_labels = None
depends_on = None

def upgrade():
    op.execute("CREATE EXTENSION IF NOT EXISTS vector;")  # pgvector for fp_morgan_b2048_r2_vec

    op.create_table(
        "compound",
        sa.Column("id", sa.BigInteger, primary_key=True),

        sa.Column("inchikey", sa.String(27), nullable=True),
        sa.Column("inchi", sa.Text, nullable=True),
        sa.Column("smiles", sa.Text, nullable=True),

        sa.Column("mol_formula", sa.String(64), nullable=True),
        sa.Column("mol_weight", sa.Float, nullable=True),
        sa.Column("exact_mass", sa.Float, nullable=True),
        sa.Column("m_plus_h", sa.Float, nullable=True),
        sa.Column("m_plus_na", sa.Float, nullable=True),

        sa.Column("fp_morgan_b2048_r2_bit", BIT(2048), nullable=True),
        sa.Column("fp_morgan_b2048_r2_pop", sa.SmallInteger, nullable=True),
        sa.Column("fp_morgan_b2048_r2_vec", Vector(2048), nullable=True),
    )

    op.create_unique_constraint("uq_compound_inchikey", "compound", ["inchikey"])
    op.create_index("ix_compound_inchikey", "compound", ["inchikey"])

    op.create_table(
        "compound_record",
        sa.Column("id", sa.BigInteger, primary_key=True),
        sa.Column("compound_id", sa.BigInteger, sa.ForeignKey("compound.id", ondelete="CASCADE"), nullable=False),

        sa.Column("source", sa.String(32), nullable=False),
        sa.Column("ext_id", sa.String(128), nullable=False),
        sa.Column("name", sa.String(512), nullable=True),
        sa.Column("synonyms", sa.ARRAY(sa.String), nullable=True),
    )

    # allow same (source, ext_id) to link to multiple compounds
    op.create_unique_constraint("uq_compound_record_source_ext", "compound_record", ["compound_id", "source", "ext_id"])
    op.create_index("ix_compound_record_compound_id", "compound_record", ["compound_id"])
    # speed lookups by accession
    op.create_index("ix_compound_record_source_ext", "compound_record", ["source", "ext_id"], unique=False)

def downgrade():
    op.drop_index("ix_compound_record_source_ext", table_name="compound_record")
    op.drop_index("ix_compound_record_compound_id", table_name="compound_record")
    op.drop_constraint("uq_compound_record_compound_source_ext", "compound_record", type_="unique")
    op.drop_table("compound_record")

    op.drop_constraint("uq_compound_inchikey", "compound", type_="unique")
    op.drop_index("ix_compound_inchikey", table_name="compound")
    op.drop_table("compound")