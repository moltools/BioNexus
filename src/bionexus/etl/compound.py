"""ETL for compound data."""

import logging
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from tqdm import tqdm
from rdkit import RDLogger
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import rdMolDescriptors

import sqlalchemy as sa
from sqlalchemy.exc import SQLAlchemyError

from retromol.io.json import iter_json
from retromol.model.result import Result
from retromol.model.rules import RuleSet
from retromol.chem.mol import mol_to_inchikey, smiles_to_mol
from retromol.chem.fingerprint import mol_to_morgan_fingerprint
from retromol.chem.tagging import remove_tags
from retromol.fingerprint.fingerprint import FingerprintGenerator

from bionexus.db.engine import SessionLocal
from bionexus.db.models import Compound


RDLogger.DisableLog("rdApp.*")


log = logging.getLogger(__name__)


@dataclass(frozen=True)
class CompoundProps:
    """
    Dataclass to hold computed compound properties.

    :var mol_weight: molecular weight
    :var c_atom_count: number of carbon atoms
    :var h_atom_count: number of hydrogen atoms
    :var n_atom_count: number of nitrogen atoms
    :var o_atom_count: number of oxygen atoms
    :var p_atom_count: number of phosphorus atoms
    :var s_atom_count: number of sulfur atoms
    :var f_atom_count: number of fluorine atoms
    :var cl_atom_count: number of chlorine atoms
    :var br_atom_count: number of bromine atoms
    :var i_atom_count: number of iodine atoms
    :var morgan_fp: Morgan fingerprint as a list of floats
    """
    
    mol_weight: float
    c_atom_count: int
    h_atom_count: int
    n_atom_count: int
    o_atom_count: int
    p_atom_count: int
    s_atom_count: int
    f_atom_count: int
    cl_atom_count: int
    br_atom_count: int
    i_atom_count: int
    morgan_fp: list[float]


def calculate_compound_props(mol: Mol) -> CompoundProps:
    """
    Calculate compound properties from a SMILES string.

    :param smiles: SMILES string of the compound
    :return: CompoundProps dataclass with computed properties
    """
    # Calculate molecular weight
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)

    # Count atom symbols
    atom_counts = Counter()
    h_atom_count = 0
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol().lower()
        atom_counts[symbol] += 1
        h_atom_count += atom.GetTotalNumHs()

    # Get Morgan fingerprint
    morgan_fp = mol_to_morgan_fingerprint(mol, radius=2, num_bits=2048, use_chirality=True)
    morgan_fp_list = [float(x) for x in morgan_fp]

    return CompoundProps(
        mol_weight=mol_weight,
        c_atom_count=atom_counts.get("c", 0),
        h_atom_count=h_atom_count,
        n_atom_count=atom_counts.get("n", 0),
        o_atom_count=atom_counts.get("o", 0),
        p_atom_count=atom_counts.get("p", 0),
        s_atom_count=atom_counts.get("s", 0),
        f_atom_count=atom_counts.get("f", 0),
        cl_atom_count=atom_counts.get("cl", 0),
        br_atom_count=atom_counts.get("br", 0),
        i_atom_count=atom_counts.get("i", 0),
        morgan_fp=morgan_fp_list,
    )


def load_compounds(jsonl: Path | str, chunk_size: int = 1_000) -> None:
    """
    Load compounds from a JSONL file into the database.

    :param jsonl: path to the JSONL file containing compound data
    :param database_name: name of the database for cross-references
    :param name_key: property key for the compound name
    :param idx_key: property key for the database cross-reference
    :param chunk_size: number of records to process in each chunk
    """
    if isinstance(jsonl, str):
        jsonl = Path(jsonl)

    ruleset = RuleSet.load_default()
    generator = FingerprintGenerator(ruleset.matching_rules)

    inserted = 0
    duplicates = 0
    failed = 0

    seen_inchikey: set[str] = set()
    batch_rows: list[dict] = []

    def flush_batch(session, rows: list[dict]) -> int:
        """
        Flush a batch of rows to the database.

        :param session: database session
        :param rows: list of row dictionaries to insert
        :return: number of rows inserted
        """
        if not rows:
            return 0
        
        stmt = (
            sa.dialects.postgresql.insert(Compound)
            .values(rows)
            .on_conflict_do_nothing(index_elements=[Compound.inchikey])
            .returning(Compound.id)
        )
        res = session.execute(stmt)
        session.commit()
        return len(res.fetchall())  # returns one row per inserted record

    with SessionLocal() as s:
        for rec in tqdm(iter_json(jsonl, jsonl=True)):
            smiles = None

            try:
                r = Result.from_dict(rec)
                smiles = r.submission.smiles
                mol = remove_tags(smiles_to_mol(smiles))
                inchikey = mol_to_inchikey(mol)

                # Batch level de-dupe of compounds
                if inchikey in seen_inchikey:
                    duplicates += 1
                    continue
                seen_inchikey.add(inchikey)

                props = calculate_compound_props(mol)
                coverage = r.calculate_coverage()
                retromol_fp_counted = [float(x) for x in generator.fingerprint_from_result(r, num_bits=1024, counted=True)]
                retromol_fp_binary = [float(int(x > 0)) for x in retromol_fp_counted]

                batch_rows.append({
                    "inchikey": inchikey,
                    "smiles": smiles,
                    "mol_weight": props.mol_weight,
                    "c_atom_count": props.c_atom_count,
                    "h_atom_count": props.h_atom_count,
                    "n_atom_count": props.n_atom_count,
                    "o_atom_count": props.o_atom_count,
                    "p_atom_count": props.p_atom_count,
                    "s_atom_count": props.s_atom_count,
                    "f_atom_count": props.f_atom_count,
                    "cl_atom_count": props.cl_atom_count,
                    "br_atom_count": props.br_atom_count,
                    "i_atom_count": props.i_atom_count,
                    "morgan_fp": props.morgan_fp,
                    "retromol_fp_counted": retromol_fp_counted,
                    "retromol_fp_binary": retromol_fp_binary,
                    "retromol": r.to_dict(),
                    "coverage": coverage,
                })

                if len(batch_rows) >= chunk_size:
                    try:
                        n_ins = flush_batch(s, batch_rows)
                        inserted += n_ins
                        duplicates += len(batch_rows) - n_ins
                    except SQLAlchemyError as e:
                        s.rollback()
                        failed += len(batch_rows)
                        log.error(f"database error during batch insert: {e}")
                    finally:
                        batch_rows.clear()
                        seen_inchikey.clear()

            except Exception as e:
                log.warning(f"failed to process compound with SMILES {smiles}: {e}")
                failed += 1
                continue
        
        # Flush any remaining rows
        if batch_rows:
            try:
                n_ins = flush_batch(s, batch_rows)
                inserted += n_ins
                duplicates += len(batch_rows) - n_ins
            except SQLAlchemyError as e:
                s.rollback()
                failed += len(batch_rows)
                log.error(f"database error during final batch insert: {e}")

    log.info(f"total compounds inserted: {inserted}")
    log.info(f"total duplicate compounds skipped: {duplicates}")
    log.info(f"total failed compounds: {failed}")
