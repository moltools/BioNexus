"""
Revision ID: 0004_rev
Revises: 0003_rev
Create Date: 2025-10-14 12:25:00.000000
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB

revision = "0004_rev"
down_revision = "0003_rev"
branch_labels = None
depends_on = None

def upgrade():
    pass

def downgrade():
    pass
