"""Database engine and session setup for BioNexus application."""

import os

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

DB_URL = os.getenv("BIONEXUS_DB_URL", "postgresql+psycopg://bionexus:bionexus@localhost:5432/bionexus")
engine = create_engine(DB_URL, pool_pre_ping=True)
SessionLocal = sessionmaker(bind=engine, autoflush=False, expire_on_commit=False)
