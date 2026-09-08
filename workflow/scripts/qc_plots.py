"""Rule `qc_plots`: QC figures from the per-cell table alone.

qc_metrics.png  histograms of the platform's QC metrics (kept vs dropped),
                with the configured thresholds drawn as lines
spatial_qc.png  the section coloured by log10 counts, detected genes, and
                kept / drop reason
"""


import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from scst_common import (  # noqa: E402
    COORD_UNIT,
    UNIT_NAME,
    XENIUM,
    categorical_palette,
    log,
    redirect_log,
    spatial_scatter,
)

# (column, axis label, log10?, [threshold keys drawn as lines])
METRIC_PANELS = {
    XENIUM: [
        ("total_counts", "transcripts / cell", True, ["min_counts", "max_counts"]),
        ("n_genes_by_counts", "genes / cell", False, ["min_genes"]),
        ("cell_area", "cell area (µm²)", False, ["min_area", "max_area"]),
        ("control_frac", "control-feature fraction", False, ["max_control_frac"]),
    ],
    "visium": [
        ("total_counts", "UMIs / spot", True, ["min_counts", "max_counts"]),
        ("n_genes_by_counts", "genes / spot", False, ["min_genes"]),
        ("pct_counts_mt", "% mitochondrial", False, ["max_pct_mt"]),
        ("pct_counts_hb", "% hemoglobin", False, ["max_pct_hb"]),
    ],
}


def plot_metrics(df, platform, thresholds, path, title):
    panels = METRIC_PANELS[platform]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    kept = df["qc_keep"].astype(bool).values
    for ax, (col, label, use_log, keys) in zip(axes.ravel(), panels):
        if col not in df.columns:
            ax.axis("off")
            ax.set_title(f"{label}: not available", fontsize=9)
            continue
        vals = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(vals)
        x = np.log10(vals + 1) if use_log else vals
        bins = np.histogram_bin_edges(x[finite], bins=50) if finite.any() else 10
        ax.hist(x[finite & kept], bins=bins, color="steelblue", alpha=0.8, label="kept")
        if (~kept).any():
            ax.hist(x[finite & ~kept], bins=bins, color="tomato", alpha=0.7, label="dropped")
        for key in keys:
            limit = float(thresholds.get(key, 0) or 0)
            if limit > 0:
                ax.axvline(np.log10(limit + 1) if use_log else limit, color="black", ls="--", lw=1)
        ax.set_xlabel(f"log10({label} + 1)" if use_log else label, fontsize=9)
        ax.set_ylabel(f"{UNIT_NAME[platform]}s", fontsize=9)
        ax.legend(fontsize=8, frameon=False)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def plot_spatial(df, platform, path, title):
    xy = df[["x", "y"]].to_numpy(dtype=float)
    unit = COORD_UNIT[platform]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    spatial_scatter(
        axes[0], xy, np.log10(df["total_counts"].to_numpy(dtype=float) + 1),
        title="log10(counts + 1)", cmap="viridis", unit=unit,
    )
    spatial_scatter(
        axes[1], xy, df["n_genes_by_counts"].to_numpy(dtype=float),
        title="detected genes", cmap="magma", unit=unit,
    )
    status = np.where(df["qc_keep"].astype(bool), "kept", df["qc_drop_reason"].astype(str))
    cats = ["kept"] + sorted({s for s in status if s != "kept"})
    palette = {"kept": "lightgrey"}
    palette.update({c: col for c, col in categorical_palette(cats[1:]).items()})
    spatial_scatter(
        axes[2], xy, pd.Categorical(status, categories=cats),
        title="QC verdict", categorical=True, palette=palette, unit=unit,
    )
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def run(cells_tsv, out_metrics, out_spatial, platform, thresholds, sample):
    df = pd.read_csv(cells_tsv, sep="\t", index_col=0)
    n_keep = int(df["qc_keep"].astype(bool).sum())
    title = f"{sample} ({platform}): {n_keep} / {len(df)} {UNIT_NAME[platform]}s kept"
    plot_metrics(df, platform, thresholds, out_metrics, title)
    plot_spatial(df, platform, out_spatial, title)
    log(f"wrote {out_metrics} and {out_spatial}")


if "snakemake" in globals():  # pragma: no cover
    redirect_log(snakemake.log[0])  # noqa: F821
    run(
        snakemake.input.cells,  # noqa: F821
        snakemake.output.metrics,  # noqa: F821
        snakemake.output.spatial,  # noqa: F821
        snakemake.params.platform,  # noqa: F821
        dict(snakemake.params.thresholds),  # noqa: F821
        snakemake.wildcards.sample,  # noqa: F821
    )
