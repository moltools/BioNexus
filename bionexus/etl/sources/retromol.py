from __future__ import annotations
from typing import List, Tuple, Dict, Any, Optional
from tqdm import tqdm
from sqlalchemy import select, delete
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.dialects.postgresql import insert
from bionexus.db.models import Compound, Ruleset, RetroMolCompound
from bionexus.db.engine import SessionLocal

# TODO