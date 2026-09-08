"""Rule `qc_summary`: one row per section from the per-sample stats JSONs,
plus a kept/dropped bar plot."""


import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from scst_common import log, redirect_log, write_tsv  # noqa: E402

LEAD_COLS = [
    "sample_id",
    "platform",
    "condition",
    "n_cells_raw",
    "n_cells_kept",
    "n_cells_dropped",
    "pct_kept",
]


def summary_table(stats_files):
    rows = []
    for p in stats_files:
        with open(p) as fh:
            d = json.load(fh)
        d.pop("thresholds", None)
        rows.append(d)
    df = pd.DataFrame(rows)
    drop_cols = sorted(c for c in df.columns if c.startswith("drop_"))
    rest = [c for c in df.columns if c not in LEAD_COLS and c not in drop_cols]
    df = df[[c for c in LEAD_COLS if c in df.columns] + drop_cols + rest]
    for c in drop_cols:
        df[c] = df[c].fillna(0).astype(int)
    return df


def plot_summary(df, path):
    n = len(df)
    fig, axes = plt.subplots(1, 2, figsize=(max(6, 1.2 * n + 4), 4.5))
    x = np.arange(n)
    axes[0].bar(x, df["n_cells_kept"], color="steelblue", label="kept")
    axes[0].bar(x, df["n_cells_dropped"], bottom=df["n_cells_kept"], color="tomato", label="dropped")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(df["sample_id"], rotation=45, ha="right", fontsize=8)
    axes[0].set_ylabel("cells / spots")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].set_title("QC verdict per section", fontsize=10)
    axes[1].bar(x, df["median_counts_kept"], color="slategrey")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(df["sample_id"], rotation=45, ha="right", fontsize=8)
    axes[1].set_ylabel("median counts (kept)")
    axes[1].set_yscale("log")
    axes[1].set_title("Depth per section", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def run(stats_files, out_tsv, out_png):
    df = summary_table(stats_files)
    write_tsv(df, out_tsv, index=False)
    plot_summary(df, out_png)
    log(f"summarised {len(df)} sections → {out_tsv}")


if "snakemake" in globals():  # pragma: no cover
    redirect_log(snakemake.log[0])  # noqa: F821
    run(list(snakemake.input.stats), snakemake.output.tsv, snakemake.output.png)  # noqa: F821
