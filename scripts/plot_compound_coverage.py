#!/usr/bin/env python3

from dotenv import load_dotenv

try:
    load_dotenv(".env")
except Exception:
    exit("could not load .env file")

import argparse
import logging

import sqlalchemy as sa
import matplotlib.pyplot as plt
from matplotlib_venn import venn3  # pip install matplotlib-venn     

from bionexus.utils.logging import setup_logging
from bionexus.db.engine import SessionLocal
from bionexus.db.models import Compound


log = logging.getLogger(__name__)


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", "-o", type=str, required=True, help="Output file path for the plot")
    return parser.parse_args()


def main() -> None:
    setup_logging(level="DEBUG")
    args = cli()

    with SessionLocal() as s:
        coverages = s.scalars(sa.select(Compound.coverage).where(Compound.coverage.isnot(None))).all()
        log.info(f"retrieved {len(coverages)} compound coverages from the database")

        total = s.scalar(sa.select(sa.func.count(Compound.id)))

        ge05 = s.scalar(
            sa.select(sa.func.count(Compound.id))
            .where(Compound.coverage.isnot(None))
            .where(Compound.coverage >= 0.5)
        )

        ge09 = s.scalar(
            sa.select(sa.func.count(Compound.id))
            .where(Compound.coverage.isnot(None))
            .where(Compound.coverage >= 0.9)
        )
        log.info(f"out of {total} compounds:")
        log.info(f"  - {ge05} ({ge05/total:.2%}) have coverage >= 0.5")
        log.info(f"  - {ge09} ({ge09/total:.2%}) have coverage >= 0.9")

    
    fig, (ax_hist, ax_venn) = plt.subplots(
        2, 1,
        figsize=(3.5, 3.9),
        gridspec_kw={"height_ratios": [2, 3]},
        constrained_layout=True
    )

    # --- Histogram ---
    ax_hist.hist(
        coverages,
        bins=20,
        color="#56b4e9",
        edgecolor="black",
        linewidth=1.0,
    )
    ax_hist.set_xlabel("coverage", fontsize=9)
    ax_hist.set_ylabel("count", fontsize=9)
    ax_hist.tick_params(axis="both", labelsize=8)
    ax_hist.set_axisbelow(True)
    ax_hist.grid(axis="y", alpha=0.6)

    # --- Venn ---
    subsets = (
        total - ge05,
        0,
        ge05 - ge09,
        0,
        0,
        0,
        ge09,
    )

    v = venn3(
        subsets=subsets,
        set_labels=("All", "≥0.5", "≥0.9"),
        ax=ax_venn,
    )

    colors = {
        "100": "#56b4e9", 
        "110": "#039e73", 
        "111": "#f0e442",
    }

    for region, color in colors.items():
        patch = v.get_patch_by_id(region)
        if patch:
            patch.set_facecolor(color)
            patch.set_edgecolor("black")
            patch.set_alpha(0.9)

    # Make venn text readable but compact
    for txt in v.set_labels:
        if txt:
            txt.set_fontsize(8)
    for txt in v.subset_labels:
        if txt:
            txt.set_fontsize(8)

    # ax_venn.set_title(
    #     f"Total={total} | ≥0.5={ge05} | ≥0.9={ge09}",
    #     fontsize=9,
    # )

    plt.savefig(args.out, dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
