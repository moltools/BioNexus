"""${message}

Revision ID: ${up_revision}
Revises: ${down_revision | comma,n}
Create Date: ${create_date}
"""

from alembic import op
import sqlalchemy as sa

revision = '${up_revision}'
down_revision = ${'None' if down_revision in (None, 'None') else repr(down_revision)}
branch_labels = ${'None' if branch_labels in (None, 'None') else repr(branch_labels)}
depends_on = ${'None' if depends_on in (None, 'None') else repr(depends_on)}

def upgrade():
    pass

def downgrade():
    pass
