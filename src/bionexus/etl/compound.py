"""ETL for compound data."""

import logging
from collections import Counter
from dataclasses import dataclass

from rdkit import RDLogger
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import rdMolDescriptors

from retromol.chem.fingerprint import mol_to_morgan_fingerprint


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
