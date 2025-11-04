from alembic import op

revision = "0005_rev"
down_revision = "0004_rev"
branch_labels = None
depends_on = None


def upgrade():
    # 1) remove duplicates (keep the greatest id per (compound_id, ruleset_id); adjust policy if you prefer)
    op.execute("""
        DELETE FROM retromol_compound a
        USING retromol_compound b
        WHERE a.compound_id = b.compound_id
          AND a.ruleset_id  = b.ruleset_id
          AND a.id < b.id;
    """)

    # 2) add the unique constraint to match the ORM
    op.create_unique_constraint(
        "uq_retromol_compound_compound_ruleset",
        "retromol_compound",
        ["compound_id", "ruleset_id"],
    )


def downgrade():
    op.drop_constraint(
        "uq_retromol_compound_compound_ruleset",
        "retromol_compound",
        type_="unique",
    )
