#!/usr/bin/env python3

from dotenv import load_dotenv

try:
    load_dotenv(".env")
except Exception:
    exit("could not load .env file")

import argparse
import logging
import warnings
import time
from dataclasses import dataclass

import numpy as np
import sqlalchemy as sa
from Bio import BiopythonDeprecationWarning
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from tqdm import tqdm

from retromol.model.rules import RuleSet
from retromol.model.submission import Submission
from retromol.model.reaction_graph import MolNode
from retromol.pipelines.parsing import run_retromol_with_timeout
from retromol.chem.mol import smiles_to_mol
from retromol.chem.fingerprint import mol_to_morgan_fingerprint, calculate_tanimoto_similarity
from retromol.fingerprint.fingerprint import FingerprintGenerator

from biocracker.query.modules import LinearReadout, PKSModule, NRPSModule, PKSExtenderUnit, linear_readout as biocracker_linear_readout
from biocracker.io.readers import load_regions
from biocracker.io.options import AntiSmashOptions
from biocracker.inference.registry import register_domain_model
from biocracker.inference.model_paras import ParasModel
from biocracker.pipelines.annotate_region import annotate_region

from bionexus.utils.logging import setup_logging
from bionexus.db.engine import SessionLocal
from bionexus.db.models import CandidateCluster

from versalign.aligner import setup_aligner
from versalign.scoring import create_substituion_matrix_dynamically
from versalign.docking import dock_against_target


log = logging.getLogger(__name__)


warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)


def cli() -> argparse.Namespace:
    """
    Command line interface for querying GBK files.

    :return: parsed command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles", type=str, required=False)
    parser.add_argument("--gbk", type=str, required=False)
    # PARAS model path needed when gbk is chosen
    parser.add_argument("--paras-model-path", type=str, required=False)
    return parser.parse_args()


def tanimoto_bits(a: np.ndarray, b: np.ndarray) -> float:
    """
    Calculate the Tanimoto similarity between two bit vectors.
    
    :param a: first bit vector
    :param b: second bit vector
    :return: Tanimoto similarity
    """
    inter = int(np.logical_and(a, b).sum())
    union = int(np.logical_or(a, b).sum())
    return inter / union if union else 0.0



@dataclass(frozen=True)
class SequenceItem:
    """
    Represents an item in a biosynthetic sequence (e.g., NRPS/PKS module or monomer).
    
    :var name: name of the item
    :var morgan_fp: Morgan fingerprint of the item's structure (if applicable)
    """
    
    name: str
    morgan_fp: ExplicitBitVect | None = None

    def __hash__(self) -> int:
        return hash((self.name, self.morgan_fp.ToBitString() if self.morgan_fp else None))

    @classmethod
    def from_nrps_module(cls, mod: NRPSModule) -> "SequenceItem":
        """
        Create a SequenceItem from an NRPS module.
        """
        if mod.substrate.smiles is not None:
            name = mod.substrate.name
            smiles = mod.substrate.smiles
            if smiles == "O=NN(O)CCC[C@H](N)(C(=O)O":  # graminine fix (fixed in >=2.0.1 versions of BioCracker)
                smiles = "O=NN(O)CCC[C@H](N)(C(=O)O)"
            mol = smiles_to_mol(smiles)
            morgan_fp = mol_to_morgan_fingerprint(mol, radius=2, num_bits=2048, use_chirality=False)
            return cls(name, morgan_fp)
        else:
            return cls("UNKNOWN")

    @classmethod
    def from_pks_module(cls, mod: PKSModule) -> "SequenceItem":
        """
        Create a SequenceItem from a PKS module.
        """
        match mod.substrate.extender_unit:
            case PKSExtenderUnit.PKS_A: name = "PKS_A"
            case PKSExtenderUnit.PKS_B: name = "PKS_B"
            case PKSExtenderUnit.PKS_C: name = "PKS_C"
            case PKSExtenderUnit.PKS_D: name = "PKS_D"
            case _: name = "PKS_A"
        return cls(name)

    @classmethod
    def from_molnode(cls, node: MolNode) -> "SequenceItem":
        """
        Create a SequenceItem from a MolNode.
        """
        if node.is_identified:
            rule = node.identity.matched_rule
            name = rule.name
            mol = smiles_to_mol(rule.smiles)
            morgan_fp = mol_to_morgan_fingerprint(mol, radius=2, num_bits=2048, use_chirality=False)
            return cls(name, morgan_fp)
        else:
            return cls("UNKNOWN")


def item_compare(a: SequenceItem | str, b: SequenceItem | str) -> float:
    """
    Compare two SequenceItems or gap representations.
    """
    if a == "-" or b == "-":
        return 0.0  # gap penalty
    
    elif isinstance(a, SequenceItem) and isinstance(b, SequenceItem):
        pks_mod_names = {"PKS_A", "PKS_B", "PKS_C", "PKS_D", "B1", "B12", "B5"}
        if a.name in pks_mod_names and b.name in pks_mod_names:
            # all PKS extender units treated as equal
            return 1.0

        elif a.name == "UNKNOWN" or b.name == "UNKNOWN":
            # could be correct, but we have no info
            return 0.5

        elif a.morgan_fp is not None and b.morgan_fp is not None:
            return calculate_tanimoto_similarity(a.morgan_fp, b.morgan_fp)
    
    return -2.0

def label_fn (r: SequenceItem | str) -> str:
    """
    """
    return str(hash(r)) if isinstance(r, SequenceItem) else r


def main() -> None:
    """
    Main function for querying GBK files.
    """
    args = cli()
    setup_logging(logging.INFO)

    ruleset = RuleSet.load_default()
    generator = FingerprintGenerator(ruleset.matching_rules)

    if args.smiles:
        smiles = args.smiles
        submission = Submission(smiles)
        result = run_retromol_with_timeout(submission, ruleset)
        log.info(f"coverage: {result.calculate_coverage() * 100:.2f}%")
        retromol_fp_counted = generator.fingerprint_from_result(result, num_bits=512, counted=True)
        retromol_fp_counted = retromol_fp_counted.astype(float).tolist()
        retromol_fp_binary = [float(int(x > 0)) for x in retromol_fp_counted]
        
        linear_readouts = result.linear_readout.paths
        linear_readout = max(linear_readouts, key=lambda x: len(x))
        log.info(f"best linear readout has {len(linear_readout)} module(s)")
        
        # Parse linear_readout into sequence of SequenceItems
        seq1: list[SequenceItem] = [SequenceItem.from_molnode(n) for n in linear_readout]
    
    elif args.gbk:
        gbk = args.gbk
        pm = ParasModel(threshold=0.1, keep_top=3, cache_dir=".", model_path=args.paras_model_path)
        register_domain_model(pm)
        options = AntiSmashOptions(readout_level="cand_cluster")
        regions = load_regions(gbk, options)
        for region in regions:
            annotate_region(region)
        log.info(f"loaded {len(regions)} region(s) from {gbk}")
        readout = biocracker_linear_readout(regions[0])
        retromol_fp_counted = generator.fingerprint_from_biocracker_readout(readout, num_bits=512, counted=True, by_orf=False)
        retromol_fp_counted = retromol_fp_counted.astype(float).tolist()
        sum(retromol_fp_counted)
        retromol_fp_binary = [float(int(x > 0)) for x in retromol_fp_counted]

        seq1: list[SequenceItem] = []
        for mod in readout.modules:
            if isinstance(mod, NRPSModule): seq1.append(SequenceItem.from_nrps_module(mod))
            elif isinstance(mod, PKSModule): seq1.append(SequenceItem.from_pks_module(mod))
            else: raise ValueError(f"unknown module type: {type(mod)}")

    # Settings
    keep_top = 10
    best_clusters = []

    t0 = time.time()

    # First extract 1000 nearest neighbors on retromol_fp
    with SessionLocal() as s:

        # Works for pgvector 0.8.0+
        s.execute(sa.text("SET LOCAL hnsw.iterative_scan = strict_order"))
        # increase how far it is allowed to scan
        s.execute(sa.text("SET LOCAL hnsw.max_scan_tuples = 1000000"))
        # optional: allow more memory for scanning
        s.execute(sa.text("SET LOCAL hnsw.scan_mem_multiplier = 2"))
        # increase ef_search for better accuracy
        s.execute(sa.text("SET LOCAL hnsw.ef_search = 5000"))

        dist = CandidateCluster.retromol_fp_counted_by_region.cosine_distance(retromol_fp_counted).label("dist")
        stmt = (
            sa.select(CandidateCluster, dist)
            .where(
                CandidateCluster.retromol_fp_counted_by_region.is_not(None),
                # CandidateCluster.file_name.ilike("BGC%"),
            )
            .order_by(dist.asc())
            .limit(100)
        )
        rows = s.execute(stmt).all()
        log.info(f"found {len(rows)} candidate clusters from ANN search")
        log.info(f"min cosine distance: {rows[0][1]:.4f}, max cosine distance: {rows[-1][1]:.4f}")

    # Rerank through docking alignment
    scored = []
    for cluster, cosine_dist in tqdm(rows):
        rec = LinearReadout.from_dict(cluster.biocracker)

        # Assembly seq2 from rec
        seq2: list[SequenceItem] = []
        for mod in rec.modules:
            if isinstance(mod, NRPSModule): seq2.append(SequenceItem.from_nrps_module(mod))
            elif isinstance(mod, PKSModule): seq2.append(SequenceItem.from_pks_module(mod))
            else: raise ValueError(f"unknown module type: {type(mod)}")
        
        if len(seq2):
            # Dynamically create scoring matrix
            items = ["-"]
            items.extend(seq1)
            items.extend(seq2)
            unique_items = list(set(items))
            sm, _ = create_substituion_matrix_dynamically(unique_items, compare=item_compare, label_fn=label_fn)
            aligner = setup_aligner(
                sm,
                "global",
                target_internal_open_gap_score=-5.0,
                target_left_open_gap_score=-5.0,
                target_right_open_gap_score=-5.0,
                query_internal_open_gap_score=-5.0,
                query_left_open_gap_score=-5.0,
                query_right_open_gap_score=-5.0,
                label_fn=label_fn,
            )
            aln = dock_against_target(
                aligner=aligner,
                target=seq1,
                candidates=[seq2],
                gap_repr="-",
                allow_block_reverse=True,
                strategy="nonoverlap",
            )
            s = aln.total_score  # the higher the score the bigger and stronge the match between the two; favors long matches

            # Penalize unmatched parts; if we are using shorter blocks
            # TODO

            # Penalize big size differences
            cum_size_seq1 = len(seq1)
            cum_size_seq2 = len(seq2)
            size_diff = abs(cum_size_seq1 - cum_size_seq2)
            size_penalty = size_diff * 0.5  # adjust penalty factor as needed
            s -= size_penalty

            if len(best_clusters) < keep_top or s > best_clusters[-1][0]:
                best_clusters.append((s, cluster, cosine_dist))
                best_clusters.sort(key=lambda x: x[0], reverse=True)
                if len(best_clusters) > keep_top:
                    best_clusters.pop()
        else:
            scored.append((cluster, -float("inf"))) 

    for i, (score, cluster, cosine_dist) in enumerate(best_clusters, 1):
        print(i, cluster.file_name, cluster.start_bp, cluster.end_bp, f"score: {score:.4f}", f"cosine: {1.0 - cosine_dist:.4f}", sep="\t")

    te = time.time()
    log.info(f"total query time: {te - t0:.2f} seconds")


if __name__ == "__main__":
    main()
