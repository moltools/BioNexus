import argparse
import os
import time
import statistics
from collections import defaultdict
import logging
import json
import random
import numpy as np

from sqlalchemy import select, func, and_

from bionexus.db.engine import SessionLocal
from bionexus.db.models import (
    GenBankRegion,
    BioCrackerGenBank,
    RetroFingerprint,
    RetroMolCompound,
    CompoundRecord,
    Annotation,
)
from tqdm import tqdm

import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

# set font to Helvetica
plt.rcParams["font.family"] = "Helvetica"

FONT_SIZE_AXIS_TITLE = 18
FONT_SIZE_AXIS_LABELS = 16
FONT_SIZE_TITLE = 20
FONT_SIZE_LEGEND = 16
PLOT_PALETTE = [
    "#e69f00",
    "#56b4e9",
]


def _bitvalue_to_int(v: object) -> int:
    """
    Convert a PostgreSQL BIT(512) value to an integer bitset.
    """
    if v is None:
        return 0

    if isinstance(v, memoryview):
        v = v.tobytes()
    if isinstance(v, (bytes, bytearray)):
        try:
            s = v.decode("ascii")
            if set(s) <= {"0", "1"}:
                return int(s, 2)
            return int.from_bytes(v, byteorder="big")
        except UnicodeDecodeError:
            return int.from_bytes(v, byteorder="big")
    if isinstance(v, str):
        return int(v, 2)

    return int(v)


def _tanimoto_bits(a: int, b: int) -> float:
    """
    Tanimoto coefficient between two bitset integers.
    """
    if a == 0 and b == 0:
        return 0.0
    inter = (a & b).bit_count()
    union = (a | b).bit_count()
    return inter / union if union else 0.0


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, help="output dir")
    return parser.parse_args()


def compute_similar_truth_compounds(
    compound_bit_fps: dict[int, int],
    true_cids: set[int],
    tc_threshold: float,
) -> set[int]:
    """
    Expand true_cids with chemically similar compounds using Tanimoto >= tc_threshold.
    Similarity is computed in compound-fingerprint bit space (compound vs compound).
    """
    expanded = set(true_cids)
    for cid in true_cids:
        fp1 = compound_bit_fps.get(cid)
        if fp1 is None:
            continue
        for other_cid, fp2 in compound_bit_fps.items():
            if other_cid in expanded:
                continue
            tc = _tanimoto_bits(fp1, fp2)
            if tc >= tc_threshold:
                expanded.add(other_cid)
    return expanded


def evaluate_mibig_bgc_to_compound_retrieval(
    top_ks: list[int],
    tc_thresholds: list[float],
    limit_bgc: int | None = None,
    metric: str = "tanimoto",  # "tanimoto" or "cosine"
    ann_scheme: str | None = None,
    ann_key: str | None = None,
    ann_values: list[str] | None = None,
    min_coverage: float | None = None,
) -> dict:
    """
    Evaluate retrieval for MIBiG GenBankRegions against all compounds.

    metric == "tanimoto":
        - BGC and compound fingerprints are bitsets (BIT(512))
        - similarity = Tanimoto(bits)

    metric == "cosine":
        - BGC and compound fingerprints use fp_retro_b512_vec_binary
        - similarity = cosine, computed in the DB via pgvector (cosine_distance)
        - we query only the top-N nearest compounds per BGC (N = max(top_ks))

    If ann_scheme, ann_key, ann_values are given, only evaluate BGCs whose
    GenBankRegion has Annotation rows such that:

        scheme == ann_scheme
        key == ann_key
        and the *set* of values == set(ann_values)

    If min_coverage is not None, only compounds with RetroMolCompound.coverage
    >= min_coverage are considered (others are ignored).
    """
    metric = metric.lower()
    assert metric in {"tanimoto", "cosine"}, "metric must be 'tanimoto' or 'cosine'"

    top_ks = sorted(set(top_ks))
    tc_thresholds = sorted(set(tc_thresholds))

    with SessionLocal() as s:
        # --- 1. Build ext_id -> set(compound_id) for MIBiG compounds ---
        mibig_comp_rows = s.execute(
            select(CompoundRecord.ext_id, CompoundRecord.compound_id)
            .where(CompoundRecord.source == "mibig")
        ).all()

        ext_to_compound_ids: dict[str, set[int]] = defaultdict(set)
        for ext_id, compound_id in mibig_comp_rows:
            ext_to_compound_ids[ext_id].add(compound_id)

        # --- 2. Load all compound RetroFingerprints (bits; vectors stay in DB) ---
        comp_rf_rows = s.execute(
            select(
                RetroFingerprint.fp_retro_b512_bit,
                RetroMolCompound.compound_id,
                RetroMolCompound.coverage,
            )
            .join(
                RetroMolCompound,
                RetroFingerprint.retromol_compound_id == RetroMolCompound.id,
            )
        ).all()

        compound_bit_fps: dict[int, int] = {}
        for fp_bits, compound_id, coverage in comp_rf_rows:
            if min_coverage is not None:
                if coverage is None or coverage < min_coverage:
                    continue

            fp_int = _bitvalue_to_int(fp_bits)
            if fp_int != 0:
                compound_bit_fps[compound_id] = fp_int

        if not compound_bit_fps:
            logger.warning("No compound bit RetroFingerprints found; aborting evaluation.")
            return {}

        # --- 3. Load MIBiG BGC RetroFingerprints (bits + vec presence flag) ---
        bgc_rows = s.execute(
            select(
                GenBankRegion.id,
                GenBankRegion.ext_id,
                RetroFingerprint.fp_retro_b512_bit,
                RetroFingerprint.fp_retro_b512_vec_binary,
            )
            .join(
                BioCrackerGenBank,
                BioCrackerGenBank.genbank_region_id == GenBankRegion.id,
            )
            .join(
                RetroFingerprint,
                RetroFingerprint.biocracker_genbank_id == BioCrackerGenBank.id,
            )
            .where(GenBankRegion.source == "mibig")
            .order_by(GenBankRegion.id)
        ).all()

        bgc_bit_fps: dict[int, list[int]] = defaultdict(list)  # region_id -> list[int]
        bgc_has_vec: dict[int, bool] = defaultdict(bool)        # region_id -> has_vec
        bgc_ext_id: dict[int, str] = {}

        for region_id, ext_id, fp_bits, fp_vec in bgc_rows:
            fp_int = _bitvalue_to_int(fp_bits)
            if fp_int != 0:
                bgc_bit_fps[region_id].append(fp_int)
            if fp_vec is not None:
                bgc_has_vec[region_id] = True
            bgc_ext_id[region_id] = ext_id

        region_ids = sorted(bgc_ext_id.keys())

        # --- 3b. Optional: restrict BGCs by annotation (scheme/key/values) ---
        if ann_scheme is not None and ann_key is not None and ann_values:
            target_value_set = set(ann_values)

            ann_rows = s.execute(
                select(Annotation.genbank_region_id, Annotation.value)
                .where(
                    Annotation.genbank_region_id.in_(region_ids),
                    Annotation.scheme == ann_scheme,
                    Annotation.key == ann_key,
                )
            ).all()

            region_to_values: dict[int, set[str]] = defaultdict(set)
            for region_id, val in ann_rows:
                region_to_values[region_id].add(val)

            region_ids = [
                rid
                for rid in region_ids
                if region_to_values.get(rid, set()) == target_value_set
            ]

        if limit_bgc is not None:
            region_ids = region_ids[:limit_bgc]

        # Precompute which ext_ids have at least one compound with a bit fp
        ext_with_valid_compounds: dict[str, set[int]] = {}
        for ext, cids in ext_to_compound_ids.items():
            valid = {cid for cid in cids if cid in compound_bit_fps}
            if valid:
                ext_with_valid_compounds[ext] = valid

        # Precompute expanded truth sets for each Tc threshold: T -> ext -> set[cid]
        expanded_truth: dict[float, dict[str, set[int]]] = {T: {} for T in tc_thresholds}
        for T in tc_thresholds:
            for ext, raw_cids in ext_with_valid_compounds.items():
                expanded_truth[T][ext] = compute_similar_truth_compounds(
                    compound_bit_fps, raw_cids, tc_threshold=T
                )

        # Aggregates and per-BGC results
        aggregates = {
            T: {
                "evaluated_bgcs": 0,
                "hit_at_k": {k: 0 for k in top_ks},
            }
            for T in tc_thresholds
        }
        per_bgc: dict[int, dict] = {}

        # Timing per BGC (ms)
        per_bgc_times_ms: list[float] = []

        logger.info(
            f"Evaluating retrieval for {len(region_ids)} MIBiG BGCs "
            f"(top_ks={top_ks}, tc_thresholds={tc_thresholds}, metric={metric}, "
            f"filter={ann_scheme}:{ann_key}:{ann_values}, min_coverage={min_coverage})"
        )

        max_k_search = max(top_ks)

        for region_id in tqdm(region_ids, desc=f"MIBiG BGC retrieval ({metric})", unit="bgc"):
            ext = bgc_ext_id.get(region_id)
            if not ext:
                continue

            if metric == "tanimoto":
                bgc_fps_bits = bgc_bit_fps.get(region_id, [])
                if not bgc_fps_bits:
                    continue
            else:  # cosine
                if not bgc_has_vec.get(region_id, False):
                    continue

            # --- 4. Compute best similarity per compound across all BGC fingerprints ---
            scores: dict[int, float] = {}

            start_t = time.perf_counter()

            if metric == "tanimoto":
                for cid, c_fp in compound_bit_fps.items():
                    best_sim = 0.0
                    for bg_fp in bgc_fps_bits:
                        sim = _tanimoto_bits(bg_fp, c_fp)
                        if sim > best_sim:
                            best_sim = sim
                    scores[cid] = best_sim
            else:
                # cosine: use pgvector in the DB; do not pull all vectors into Python
                bgc_vec_rows = s.execute(
                    select(
                        RetroFingerprint.fp_retro_b512_vec_binary,
                        RetroFingerprint.id,
                    )
                    .join(
                        BioCrackerGenBank,
                        RetroFingerprint.biocracker_genbank_id == BioCrackerGenBank.id,
                    )
                    .where(BioCrackerGenBank.genbank_region_id == region_id)
                ).all()

                for bg_vec, _rf_id in bgc_vec_rows:
                    if bg_vec is None:
                        continue

                    distance_expr = RetroFingerprint.fp_retro_b512_vec_binary.cosine_distance(bg_vec)

                    q = (
                        select(
                            RetroMolCompound.compound_id,
                            (1 - distance_expr).label("sim"),
                        )
                        .join(
                            RetroMolCompound,
                            RetroFingerprint.retromol_compound_id == RetroMolCompound.id,
                        )
                        .where(RetroFingerprint.fp_retro_b512_vec_binary.isnot(None))
                        .order_by(distance_expr)
                        .limit(max_k_search)
                    )

                    rows = s.execute(q).all()
                    for cid, sim in rows:
                        if cid not in scores or sim > scores[cid]:
                            scores[cid] = sim

            end_t = time.perf_counter()
            elapsed_ms = (end_t - start_t) * 1000.0

            if not scores:
                continue

            ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
            if not ranked or ranked[0][1] <= 0.0:
                continue

            per_bgc_times_ms.append(elapsed_ms)

            top1_cid, top1_score = ranked[0]

            per_bgc_entry = per_bgc.setdefault(
                region_id,
                {
                    "ext_id": ext,
                    "top1": {"cid": top1_cid, "sim": top1_score},
                },
            )
            per_bgc_entry.setdefault("best_true_sim", {})
            per_bgc_entry.setdefault("hit_at_k", {})

            # For each Tc threshold, compute per-BGC stats
            for T in tc_thresholds:
                true_cids = expanded_truth[T].get(ext, set())
                if not true_cids:
                    continue

                aggregates[T]["evaluated_bgcs"] += 1

                best_true_sim = max(scores.get(cid, 0.0) for cid in true_cids)
                per_bgc_entry["best_true_sim"][T] = best_true_sim

                hits_for_T = per_bgc_entry["hit_at_k"].setdefault(T, {})
                for k in top_ks:
                    topk_cids = [cid for cid, _ in ranked[:k]]
                    hit = any(cid in true_cids for cid in topk_cids)
                    hits_for_T[k] = hit
                    if hit:
                        aggregates[T]["hit_at_k"][k] += 1

        # Compute rates
        for T in tc_thresholds:
            eval_n = aggregates[T]["evaluated_bgcs"]
            if eval_n == 0:
                continue
            aggregates[T]["hit_at_k_rate"] = {
                k: aggregates[T]["hit_at_k"][k] / eval_n
                for k in top_ks
            }

        # Timing summary
        if per_bgc_times_ms:
            mean_ms = statistics.mean(per_bgc_times_ms)
            std_ms = statistics.pstdev(per_bgc_times_ms) if len(per_bgc_times_ms) > 1 else 0.0
        else:
            mean_ms = 0.0
            std_ms = 0.0

        timing = {
            "mean_ms": mean_ms,
            "std_ms": std_ms,
            "n_bgcs_timed": len(per_bgc_times_ms),
        }

        logger.info(
            "Timing (%s, min_coverage=%s): mean=%.2f ms, std=%.2f ms over %d BGCs",
            metric, str(min_coverage), mean_ms, std_ms, len(per_bgc_times_ms),
        )

        return {
            "metric": metric,
            "top_ks": top_ks,
            "tc_thresholds": tc_thresholds,
            "aggregates": aggregates,
            "per_bgc": per_bgc,
            "filter": {
                "scheme": ann_scheme,
                "key": ann_key,
                "values": ann_values,
            },
            "min_coverage": min_coverage,
            "timing": timing,
        }


def evaluate_mibig_bgc_for_biosynthesis_products(
    product_values: list[str],
    top_ks: list[int],
    tc_thresholds: list[float],
    limit_bgc: int | None = None,
    metric: str = "tanimoto",
    min_coverage: float | None = None,
) -> dict:
    """
    Convenience wrapper to restrict to BGCs with biosynthesis/product annotations
    having exactly the given set of values.
    """
    return evaluate_mibig_bgc_to_compound_retrieval(
        top_ks=top_ks,
        tc_thresholds=tc_thresholds,
        limit_bgc=limit_bgc,
        metric=metric,
        ann_scheme="biosynthesis",
        ann_key="product",
        ann_values=product_values,
        min_coverage=min_coverage,
    )


def load_results(path: str) -> dict:
    with open(path, "r") as f:
        return json.load(f)


def extract_lines_with_counts(results: dict):
    """
    Convert benchmark output into:
        ks = [...]
        lines = { T (float) : [hit@k for k in ks] }
        counts = { T (float) : evaluated_bgcs }

    Handles JSON stringification of dict keys.
    """
    ks = results["top_ks"]
    aggr = results["aggregates"]

    lines: dict[float, list[float]] = {}
    counts: dict[float, int] = {}

    for T_key, data in aggr.items():
        hit_at_k_rate = data.get("hit_at_k_rate")
        if not hit_at_k_rate:
            continue

        T = float(T_key)
        hit_rates = [hit_at_k_rate.get(str(k), 0.0) for k in ks]
        lines[T] = hit_rates
        counts[T] = data.get("evaluated_bgcs", 0)

    return ks, lines, counts


def extract_timing(path: str) -> tuple[float, float, int]:
    """
    Return (mean_ms, std_ms, n_bgcs_timed) from a benchmark JSON.
    """
    res = load_results(path)
    timing = res.get("timing", {}) or {}
    mean_ms = timing.get("mean_ms", 0.0)
    std_ms = timing.get("std_ms", 0.0)
    n = timing.get("n_bgcs_timed", 0)
    return mean_ms, std_ms, n


def plot_three_panels_compare_coverage(
    all_any_path: str,
    all_cov_path: str,
    t1_any_path: str,
    t1_cov_path: str,
    nrps_any_path: str,
    nrps_cov_path: str,
    out_path: str,
):
    """
    3 panels: All / T1PKS / NRPS.
    In each panel: 4 lines:
        - Tc ≥ 0.4, coverage any
        - Tc ≥ 0.4, coverage ≥ 0.9
        - Tc ≥ 0.675, coverage any
        - Tc ≥ 0.675, coverage ≥ 0.9
    """
    res_all_any = load_results(all_any_path)
    res_all_cov = load_results(all_cov_path)
    res_t1_any = load_results(t1_any_path)
    res_t1_cov = load_results(t1_cov_path)
    res_nrps_any = load_results(nrps_any_path)
    res_nrps_cov = load_results(nrps_cov_path)

    ks_all, lines_all_any, counts_all_any = extract_lines_with_counts(res_all_any)
    ks_all2, lines_all_cov, counts_all_cov = extract_lines_with_counts(res_all_cov)

    ks_t1, lines_t1_any, counts_t1_any = extract_lines_with_counts(res_t1_any)
    ks_t12, lines_t1_cov, counts_t1_cov = extract_lines_with_counts(res_t1_cov)

    ks_nrps, lines_nrps_any, counts_nrps_any = extract_lines_with_counts(res_nrps_any)
    ks_nrps2, lines_nrps_cov, counts_nrps_cov = extract_lines_with_counts(res_nrps_cov)

    assert ks_all == ks_all2
    assert ks_t1 == ks_t12
    assert ks_nrps == ks_nrps2

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)

    panels = [
        ("All BGCs", ks_all, lines_all_any, lines_all_cov, counts_all_any, counts_all_cov),
        ("T1PKS BGCs", ks_t1, lines_t1_any, lines_t1_cov, counts_t1_any, counts_t1_cov),
        ("NRPS BGCs", ks_nrps, lines_nrps_any, lines_nrps_cov, counts_nrps_any, counts_nrps_cov),
    ]

    colors = {0.4: PLOT_PALETTE[0], 0.675: PLOT_PALETTE[1]}
    linestyles = {"any": "-", "cov": "--"}
    yticks = [i / 10 for i in range(0, 11)]

    for ax, (title, ks, lines_any, lines_cov, counts_any, counts_cov) in zip(axes, panels):
        thresholds = sorted(set(lines_any.keys()) | set(lines_cov.keys()))
        for T in thresholds:
            color = colors.get(T, None)

            if T in lines_any:
                hr_any = lines_any[T]
                n_any = counts_any.get(T, 0)
                ax.plot(
                    ks,
                    hr_any,
                    marker="o",
                    linestyle=linestyles["any"],
                    color=color,
                    linewidth=2,
                    label=f"Tc ≥ {T}, cov any (N={n_any})",
                )

            if T in lines_cov:
                hr_cov = lines_cov[T]
                n_cov = counts_cov.get(T, 0)
                ax.plot(
                    ks,
                    hr_cov,
                    marker="s",
                    linestyle=linestyles["cov"],
                    color=color,
                    linewidth=2,
                    label=f"Tc ≥ {T}, cov ≥ 0.9 (N={n_cov})",
                )

        ax.set_title(title, fontsize=FONT_SIZE_TITLE)
        ax.set_xlabel("Top-K", fontsize=FONT_SIZE_AXIS_TITLE)
        ax.set_xticks(ks)
        ax.set_ylim(0.0, 1.0)
        ax.set_yticks(yticks, [f"{ytick:.1f}" for ytick in yticks], fontsize=FONT_SIZE_AXIS_LABELS)
        ax.set_xticklabels([str(k) for k in ks], fontsize=FONT_SIZE_AXIS_LABELS)
        ax.grid(True, axis="y", alpha=0.3)
        ax.legend(fontsize=FONT_SIZE_LEGEND)

    axes[0].set_ylabel("Hit rate", fontsize=FONT_SIZE_AXIS_TITLE)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    print(f"Saved figure → {out_path}")


def plot_three_panels_timing(
    metric: str,
    all_any_path: str,
    all_cov_path: str,
    t1_any_path: str,
    t1_cov_path: str,
    nrps_any_path: str,
    nrps_cov_path: str,
    out_path: str,
):
    """
    3 panels: All / T1PKS / NRPS.
    In each panel: 2 bars:
        - coverage any
        - coverage ≥ 0.9

    Bars: mean ms, error bars = std ms.
    """
    mean_all_any, std_all_any, n_all_any = extract_timing(all_any_path)
    mean_all_cov, std_all_cov, n_all_cov = extract_timing(all_cov_path)
    mean_t1_any, std_t1_any, n_t1_any = extract_timing(t1_any_path)
    mean_t1_cov, std_t1_cov, n_t1_cov = extract_timing(t1_cov_path)
    mean_nrps_any, std_nrps_any, n_nrps_any = extract_timing(nrps_any_path)
    mean_nrps_cov, std_nrps_cov, n_nrps_cov = extract_timing(nrps_cov_path)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)

    panels = [
        ("All BGCs", (mean_all_any, std_all_any, n_all_any), (mean_all_cov, std_all_cov, n_all_cov)),
        ("T1PKS BGCs", (mean_t1_any, std_t1_any, n_t1_any), (mean_t1_cov, std_t1_cov, n_t1_cov)),
        ("NRPS BGCs", (mean_nrps_any, std_nrps_any, n_nrps_any), (mean_nrps_cov, std_nrps_cov, n_nrps_cov)),
    ]

    x_labels = ["cov any", "cov ≥ 0.9"]
    x_pos = [0, 1]

    for ax, (title, (m_any, s_any, n_any), (m_cov, s_cov, n_cov)) in zip(axes, panels):
        means = [m_any, m_cov]
        stds = [s_any, s_cov]
        ns = [n_any, n_cov]

        bars = ax.bar(
            x_pos,
            means,
            yerr=stds,
            capsize=5,
            align="center",
        )

        # annotate N above each bar
        for xpos, bar, n in zip(x_pos, bars, ns):
            height = bar.get_height()
            ax.text(
                xpos,
                height,
                f"N={n}",
                ha="center",
                va="bottom",
                fontsize=9,
            )

        ax.set_title(title, fontsize=15)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(x_labels, rotation=0)
        ax.set_ylabel("Time per BGC (ms)" if ax is axes[0] else "")
        ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(f"Retrieval time per BGC ({metric})", fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(out_path, dpi=300)
    print(f"Saved timing figure → {out_path}")


def sample_true_vs_random_similarities(
    *,
    metric: str,
    ann_scheme: str | None = None,
    ann_key: str | None = None,
    ann_values: list[str] | None = None,
    min_coverage: float | None = None,
    limit_bgc: int | None = None,
    max_pairs: int = 20000,
    random_per_true: int = 1,
    seed: int = 0,
) -> dict:
    """
    Sample similarity distributions for:
      - true pairs: (BGC, compound in its MIBiG mapping)
      - random pairs: (same BGC, random compound not in its truth set)

    Currently implemented for metric == "tanimoto" only.
    Returns dict with lists: true_sims, random_sims, plus some counts.
    """
    metric = metric.lower()
    assert metric in {"tanimoto", "cosine"}
    if metric != "tanimoto":
        raise NotImplementedError("This sampler currently supports metric='tanimoto' only.")

    rng = random.Random(seed)

    with SessionLocal() as s:
        # --- 1) ext_id -> set(compound_id) for MIBiG compounds ---
        mibig_comp_rows = s.execute(
            select(CompoundRecord.ext_id, CompoundRecord.compound_id)
            .where(CompoundRecord.source == "mibig")
        ).all()

        ext_to_compound_ids: dict[str, set[int]] = defaultdict(set)
        for ext_id, compound_id in mibig_comp_rows:
            ext_to_compound_ids[ext_id].add(compound_id)

        # --- 2) load compound bit fps (optionally filter by coverage) ---
        comp_rf_rows = s.execute(
            select(
                RetroFingerprint.fp_retro_b512_bit,
                RetroMolCompound.compound_id,
                RetroMolCompound.coverage,
            )
            .join(
                RetroMolCompound,
                RetroFingerprint.retromol_compound_id == RetroMolCompound.id,
            )
        ).all()

        compound_bit_fps: dict[int, int] = {}
        for fp_bits, compound_id, coverage in comp_rf_rows:
            if min_coverage is not None:
                if coverage is None or coverage < min_coverage:
                    continue
            fp_int = _bitvalue_to_int(fp_bits)
            if fp_int != 0:
                compound_bit_fps[compound_id] = fp_int

        if not compound_bit_fps:
            logger.warning("No compound bit RetroFingerprints found; aborting sampling.")
            return {}

        all_compound_ids = list(compound_bit_fps.keys())

        # --- 3) load BGC bit fps + ext_id ---
        bgc_rows = s.execute(
            select(
                GenBankRegion.id,
                GenBankRegion.ext_id,
                RetroFingerprint.fp_retro_b512_bit,
            )
            .join(
                BioCrackerGenBank,
                BioCrackerGenBank.genbank_region_id == GenBankRegion.id,
            )
            .join(
                RetroFingerprint,
                RetroFingerprint.biocracker_genbank_id == BioCrackerGenBank.id,
            )
            .where(GenBankRegion.source == "mibig")
            .order_by(GenBankRegion.id)
        ).all()

        bgc_bit_fps: dict[int, list[int]] = defaultdict(list)
        bgc_ext_id: dict[int, str] = {}

        for region_id, ext_id, fp_bits in bgc_rows:
            fp_int = _bitvalue_to_int(fp_bits)
            if fp_int != 0:
                bgc_bit_fps[region_id].append(fp_int)
            bgc_ext_id[region_id] = ext_id

        region_ids = sorted(bgc_ext_id.keys())

        # --- 3b) optional annotation filtering (exact set match, same as evaluation) ---
        if ann_scheme is not None and ann_key is not None and ann_values:
            target_value_set = set(ann_values)
            ann_rows = s.execute(
                select(Annotation.genbank_region_id, Annotation.value)
                .where(
                    Annotation.genbank_region_id.in_(region_ids),
                    Annotation.scheme == ann_scheme,
                    Annotation.key == ann_key,
                )
            ).all()

            region_to_values: dict[int, set[str]] = defaultdict(set)
            for region_id, val in ann_rows:
                region_to_values[region_id].add(val)

            region_ids = [
                rid for rid in region_ids
                if region_to_values.get(rid, set()) == target_value_set
            ]

        if limit_bgc is not None:
            region_ids = region_ids[:limit_bgc]

        # Precompute ext -> valid true cids (must exist in compound_bit_fps)
        ext_with_valid_compounds: dict[str, set[int]] = {}
        for ext, cids in ext_to_compound_ids.items():
            valid = {cid for cid in cids if cid in compound_bit_fps}
            if valid:
                ext_with_valid_compounds[ext] = valid

        true_sims: list[float] = []
        random_sims: list[float] = []

        def best_sim_to_bgc(bgc_fps: list[int], comp_fp: int) -> float:
            best = 0.0
            for bg_fp in bgc_fps:
                sim = _tanimoto_bits(bg_fp, comp_fp)
                if sim > best:
                    best = sim
            return best

        # --- 4) sample pairs, capped by max_pairs ---
        for rid in region_ids:
            bgc_fps = bgc_bit_fps.get(rid, [])
            if not bgc_fps:
                continue

            ext = bgc_ext_id.get(rid)
            if not ext:
                continue

            true_cids = ext_with_valid_compounds.get(ext, set())
            if not true_cids:
                continue

            # True pairs
            for cid in true_cids:
                if len(true_sims) >= max_pairs:
                    break
                sim = best_sim_to_bgc(bgc_fps, compound_bit_fps[cid])
                true_sims.append(sim)

                # For each true pair, sample random_per_true negatives from global pool
                for _ in range(random_per_true):
                    if len(random_sims) >= max_pairs:
                        break

                    # try a few times to avoid picking a true cid
                    for _attempt in range(10):
                        rcid = rng.choice(all_compound_ids)
                        if rcid not in true_cids:
                            break
                    else:
                        # fallback: allow true (rare if truth set is large)
                        rcid = rng.choice(all_compound_ids)

                    rsim = best_sim_to_bgc(bgc_fps, compound_bit_fps[rcid])
                    random_sims.append(rsim)

            if len(true_sims) >= max_pairs and len(random_sims) >= max_pairs:
                break

        return {
            "metric": metric,
            "filter": {"scheme": ann_scheme, "key": ann_key, "values": ann_values},
            "min_coverage": min_coverage,
            "n_bgcs_considered": len(region_ids),
            "n_true_pairs": len(true_sims),
            "n_random_pairs": len(random_sims),
            "true_sims": true_sims,
            "random_sims": random_sims,
            "seed": seed,
            "max_pairs": max_pairs,
            "random_per_true": random_per_true,
        }


def _summarize_dist(x: list[float]) -> dict:
    if not x:
        return {"n": 0}
    arr = np.array(x, dtype=float)
    return {
        "n": int(arr.size),
        "mean": float(arr.mean()),
        "p50": float(np.percentile(arr, 50)),
        "p90": float(np.percentile(arr, 90)),
        "p95": float(np.percentile(arr, 95)),
        "p99": float(np.percentile(arr, 99)),
        "max": float(arr.max()),
    }


def plot_true_vs_random_similarity_panels(
    *,
    metric: str,
    out_path: str,
    all_data: dict,
    t1_data: dict,
    nrps_data: dict,
    bins: int = 50,
):
    """
    3 panels: All / T1PKS / NRPS.
    Each panel overlays histograms of true vs random similarity.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)

    panels = [
        ("All BGCs", all_data),
        ("T1PKS BGCs", t1_data),
        ("NRPS BGCs", nrps_data),
    ]

    for ax, (title, dat) in zip(axes, panels):
        true_sims = dat.get("true_sims", [])
        rand_sims = dat.get("random_sims", [])
        if not true_sims or not rand_sims:
            ax.set_title(f"{title}\n(no data)")
            continue

        # consistent binning across [0,1]
        edges = np.linspace(0.0, 1.0, bins + 1)

        ax.hist(rand_sims, bins=edges, density=True, alpha=0.5, label=f"random (n={len(rand_sims)})", color=PLOT_PALETTE[0])
        ax.hist(true_sims, bins=edges, density=True, alpha=0.5, label=f"true (n={len(true_sims)})", color=PLOT_PALETTE[1])

        t_med = np.median(true_sims)
        r_med = np.median(rand_sims)
        ax.axvline(r_med, linewidth=2, linestyle="--", label=f"median random={r_med:.3f}", color="black")
        ax.axvline(t_med, linewidth=2, linestyle="-", label=f"median true={t_med:.3f}", color="black")

        ax.set_title(title, fontsize=FONT_SIZE_TITLE)
        ax.set_xlabel("Similarity", fontsize=FONT_SIZE_AXIS_TITLE)
        ax.set_xlim(0.0, 1.0)
        ax.set_yticklabels([f"{ytick:.2f}" for ytick in ax.get_yticks()], fontsize=FONT_SIZE_AXIS_LABELS)
        ax.set_xticklabels([f"{xtick:.1f}" for xtick in ax.get_xticks()], fontsize=FONT_SIZE_AXIS_LABELS)
        ax.grid(True, axis="y", alpha=0.3)
        ax.legend(fontsize=FONT_SIZE_LEGEND)

    axes[0].set_ylabel("Density", fontsize=FONT_SIZE_AXIS_TITLE)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(out_path, dpi=300)
    print(f"Saved similarity distribution figure → {out_path}")


def plot_compound_coverage_distribution(
    *,
    out_path: str,
    bin_width: float = 0.01,  # bars from 0.00..1.00 in steps of 0.01
) -> dict:
    """
    For each compound_id, take the *best* (max) RetroMolCompound.coverage across rows.
    Plot a bar histogram over coverage in [0, 1] and report counts >= 0.9 and >= 0.5.

    Returns:
        {
          "n_compounds_with_coverage": int,
          "count_ge_0_9": int,
          "count_ge_0_5": int,
          "fraction_ge_0_9": float,
          "fraction_ge_0_5": float,
        }
    """
    assert 0 < bin_width <= 0.1, "bin_width should be reasonable (e.g., 0.01 or 0.02)"

    with SessionLocal() as s:
        # Best coverage per compound_id (max over rows)
        rows = s.execute(
            select(
                RetroMolCompound.compound_id,
                func.max(RetroMolCompound.coverage).label("best_cov"),
            )
            .where(RetroMolCompound.coverage.isnot(None))
            .group_by(RetroMolCompound.compound_id)
        ).all()

    best_covs = [float(cov) for _cid, cov in rows if cov is not None]

    if not best_covs:
        logger.warning("No compound coverage values found; nothing to plot.")
        return {
            "n_compounds_with_coverage": 0,
            "count_ge_0_9": 0,
            "count_ge_0_5": 0,
            "fraction_ge_0_9": 0.0,
            "fraction_ge_0_5": 0.0,
        }

    # Clamp into [0, 1] just in case
    best_covs = [min(1.0, max(0.0, x)) for x in best_covs]

    n = len(best_covs)
    count_ge_0_9 = sum(1 for x in best_covs if x >= 0.9)
    count_ge_0_5 = sum(1 for x in best_covs if x >= 0.5)

    # Histogram bins from 0..1
    bins = int(round(1.0 / bin_width))
    edges = np.linspace(0.0, 1.0, bins + 1)
    hist, _ = np.histogram(best_covs, bins=edges)

    # Plot as bars
    centers = (edges[:-1] + edges[1:]) / 2.0
    width = (edges[1] - edges[0])

    plt.figure(figsize=(5, 4))
    plt.bar(centers, hist, width=width, align="center", color=PLOT_PALETTE[1], edgecolor="black", linewidth=1)
    plt.xlim(0.0, 1.0)
    plt.xlabel("coverage", fontsize=FONT_SIZE_AXIS_TITLE)
    plt.ylabel("count", fontsize=FONT_SIZE_AXIS_TITLE)
    plt.grid(True, axis="y", alpha=0.3)

    # annotate the two threshold counts
    plt.text(
        0.02, 0.98,
        f"N compounds={n}\n"
        f"coverage ≥ 0.9: {count_ge_0_9} ({count_ge_0_9/n:.1%})\n"
        f"coverage ≥ 0.5: {count_ge_0_5} ({count_ge_0_5/n:.1%})",
        transform=plt.gca().transAxes,
        va="top",
        fontsize=12,
    )

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    print(f"Saved coverage distribution → {out_path}")

    # Also print counts plainly (as you asked)
    print(f"Compounds with best coverage ≥ 0.9: {count_ge_0_9} / {n} ({count_ge_0_9/n:.2%})")
    print(f"Compounds with best coverage ≥ 0.5: {count_ge_0_5} / {n} ({count_ge_0_5/n:.2%})")

    return {
        "n_compounds_with_coverage": n,
        "count_ge_0_9": count_ge_0_9,
        "count_ge_0_5": count_ge_0_5,
        "fraction_ge_0_9": count_ge_0_9 / n,
        "fraction_ge_0_5": count_ge_0_5 / n,
    }


def benchmark(out_dir: str):
    os.makedirs(out_dir, exist_ok=True)

    plot_compound_coverage_distribution(
        out_path=os.path.join(out_dir, "compound_coverage_distribution_best_per_compound.svg"),
        bin_width=0.05,
    )

    metrics = ["tanimoto", "cosine"]
    top_ks = [1, 5, 10, 50, 100]
    tc_thresholds = [0.4, 0.675]

    coverage_settings = [
        (None, "any"),
        (0.9, "cov_0.9"),
    ]

    subsets = [
        ("all", None),
        ("T1PKS", ["T1PKS"]),
        ("NRPS", ["NRPS"]),
    ]

    for metric in metrics:
        # Run benchmarks for each subset + coverage setting
        for subset_name, product_filter in subsets:
            for min_cov, cov_label in coverage_settings:
                filename = f"benchmark_{metric}_{subset_name}_cov_{cov_label}.json"
                out_path = os.path.join(out_dir, filename)

                if os.path.exists(out_path):
                    logger.info("Skipping existing benchmark output: %s", filename)
                    continue

                if product_filter is None:
                    results = evaluate_mibig_bgc_to_compound_retrieval(
                        top_ks=top_ks,
                        tc_thresholds=tc_thresholds,
                        limit_bgc=None,
                        metric=metric,
                        min_coverage=min_cov,
                    )
                else:
                    results = evaluate_mibig_bgc_for_biosynthesis_products(
                        product_values=product_filter,
                        top_ks=top_ks,
                        tc_thresholds=tc_thresholds,
                        limit_bgc=None,
                        metric=metric,
                        min_coverage=min_cov,
                    )

                with open(out_path, "w") as f:
                    json.dump(results, f, indent=2)
                logger.info(
                    "Wrote %s benchmark results for subset=%s, cov=%s to %s",
                    metric, subset_name, cov_label, out_path
                )

        # Plot hit-rate comparison (cov-any vs cov-0.9) for this metric
        plot_three_panels_compare_coverage(
            all_any_path=os.path.join(out_dir, f"benchmark_{metric}_all_cov_any.json"),
            all_cov_path=os.path.join(out_dir, f"benchmark_{metric}_all_cov_cov_0.9.json"),
            t1_any_path=os.path.join(out_dir, f"benchmark_{metric}_T1PKS_cov_any.json"),
            t1_cov_path=os.path.join(out_dir, f"benchmark_{metric}_T1PKS_cov_cov_0.9.json"),
            nrps_any_path=os.path.join(out_dir, f"benchmark_{metric}_NRPS_cov_any.json"),
            nrps_cov_path=os.path.join(out_dir, f"benchmark_{metric}_NRPS_cov_cov_0.9.json"),
            out_path=os.path.join(out_dir, f"benchmark_hit_rates_cov_compare_{metric}.svg"),
        )

        # Plot timing comparison for this metric
        plot_three_panels_timing(
            metric=metric,
            all_any_path=os.path.join(out_dir, f"benchmark_{metric}_all_cov_any.json"),
            all_cov_path=os.path.join(out_dir, f"benchmark_{metric}_all_cov_cov_0.9.json"),
            t1_any_path=os.path.join(out_dir, f"benchmark_{metric}_T1PKS_cov_any.json"),
            t1_cov_path=os.path.join(out_dir, f"benchmark_{metric}_T1PKS_cov_cov_0.9.json"),
            nrps_any_path=os.path.join(out_dir, f"benchmark_{metric}_NRPS_cov_any.json"),
            nrps_cov_path=os.path.join(out_dir, f"benchmark_{metric}_NRPS_cov_cov_0.9.json"),
            out_path=os.path.join(out_dir, f"benchmark_timing_cov_compare_{metric}.svg"),
        )

        if metric != "tanimoto":
            # similarity distribution sampling only implemented for tanimoto
            continue

        for min_cov, cov_label in coverage_settings:
            # All
            all_dat = sample_true_vs_random_similarities(
                metric=metric,
                min_coverage=min_cov,
                ann_scheme=None,
                ann_key=None,
                ann_values=None,
                max_pairs=20000,
                random_per_true=1,
                seed=0,
            )
            # T1PKS
            t1_dat = sample_true_vs_random_similarities(
                metric=metric,
                min_coverage=min_cov,
                ann_scheme="biosynthesis",
                ann_key="product",
                ann_values=["T1PKS"],
                max_pairs=20000,
                random_per_true=1,
                seed=0,
            )
            # NRPS
            nrps_dat = sample_true_vs_random_similarities(
                metric=metric,
                min_coverage=min_cov,
                ann_scheme="biosynthesis",
                ann_key="product",
                ann_values=["NRPS"],
                max_pairs=20000,
                random_per_true=1,
                seed=0,
            )

            # optional: print quick summaries
            logger.info("Similarity dist summary (%s, cov=%s) ALL true=%s random=%s",
                        metric, cov_label, _summarize_dist(all_dat.get("true_sims", [])), _summarize_dist(all_dat.get("random_sims", [])))
            logger.info("Similarity dist summary (%s, cov=%s) T1PKS true=%s random=%s",
                        metric, cov_label, _summarize_dist(t1_dat.get("true_sims", [])), _summarize_dist(t1_dat.get("random_sims", [])))
            logger.info("Similarity dist summary (%s, cov=%s) NRPS true=%s random=%s",
                        metric, cov_label, _summarize_dist(nrps_dat.get("true_sims", [])), _summarize_dist(nrps_dat.get("random_sims", [])))

            plot_true_vs_random_similarity_panels(
                metric=metric,
                out_path=os.path.join(out_dir, f"benchmark_true_vs_random_similarity_{metric}_cov_{cov_label}.svg"),
                all_data=all_dat,
                t1_data=t1_dat,
                nrps_data=nrps_dat,
                bins=60,
            )