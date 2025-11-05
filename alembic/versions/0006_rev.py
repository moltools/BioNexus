import sqlalchemy as sa
from alembic import op
from pgvector.sqlalchemy import Vector
from sqlalchemy.dialects.postgresql import BIT

revision = "0006_rev"
down_revision = "0005_rev"
branch_labels = None
depends_on = None


def upgrade():
    op.execute("CREATE EXTENSION IF NOT EXISTS vector")  # ensure pgvector is available (no-op if already installed)
    op.drop_column("retromol_compound", "fp_retro")
    op.add_column("retromol_compound", sa.Column("coverage", sa.Float(), nullable=True))


def downgrade():
    op.add_column("retromol_compound", sa.Column("fp_retro", BIT(512), nullable=True))
    op.drop_column("retromol_compound", "coverage")
