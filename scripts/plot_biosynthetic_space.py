"""Plot biosynthetic space embeddings."""

import argparse
import pandas as pd

import matplotlib.pyplot as plt


COLOR_BLUE = "#5db4e5"
COLOR_GREEN = "#029e73"


def cli() -> argparse.Namespace:
    """Command line interface."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, required=True, help="Path to TSV data embedding space.")
    parser.add_argument("-o", type=str, required=True, help="Output path for plot (.svg for best quality).")
    return parser.parse_args()


def main() -> None:
    args = cli()

    data = pd.read_csv(args.i, sep="\t")
    kind = data["kind"].values
    x = data["x"].values
    y = data["y"].values

    plt.figure(figsize=(6, 6))
    for k in set(kind):
        mask = kind == k
        color = COLOR_BLUE if k == "gene_cluster" else COLOR_GREEN
        # give black line
        plt.scatter(x[mask], y[mask], label=k, s=40, color=color, edgecolor="black", linewidth=1)

    plt.xlabel("dimension 1")
    plt.ylabel("dimension 2")
    # plt.legend()

    # no axis ticks and labels
    plt.xticks([])
    plt.yticks([])

    plt.savefig(args.o, bbox_inches="tight")



if __name__ == "__main__":
    main()
