"""Shared helpers for the workflow's Python ``script:`` rules.

Snakemake puts this directory on ``sys.path`` for local scripts, so every rule
script can ``import scst_common``. The module is deliberately free of
scanpy/squidpy: the loader and QC modules (and their unit tests) need only
anndata / h5py / numpy / pandas / scipy, so the static CI job can exercise
them with a pip install.
"""


import json
import os
import sys

import numpy as np
import pandas as pd

XENIUM = "xenium"
VISIUM = "visium"
PLATFORMS = (XENIUM, VISIUM)
COORD_UNIT = {XENIUM: "µm", VISIUM: "px"}
UNIT_NAME = {XENIUM: "cell", VISIUM: "spot"}
GENE_FEATURE_TYPE = "Gene Expression"

# Cell/spot drop rules per platform, applied IN THIS ORDER: the first rule a
# cell fails is recorded as its drop reason. Each entry is
# (reason, obs column, config key, comparison), with a config value of 0
# meaning "rule disabled".
DROP_RULES = {
    XENIUM: [
        ("low_counts", "total_counts", "min_counts", "<"),
        ("high_counts", "total_counts", "max_counts", ">"),
        ("low_genes", "n_genes_by_counts", "min_genes", "<"),
        ("small_area", "cell_area", "min_area", "<"),
        ("large_area", "cell_area", "max_area", ">"),
        ("high_control_frac", "control_frac", "max_control_frac", ">"),
    ],
    VISIUM: [
        ("low_counts", "total_counts", "min_counts", "<"),
        ("high_counts", "total_counts", "max_counts", ">"),
        ("low_genes", "n_genes_by_counts", "min_genes", "<"),
        ("high_mito", "pct_counts_mt", "max_pct_mt", ">"),
        ("high_hb", "pct_counts_hb", "max_pct_hb", ">"),
    ],
}


# ── logging ─────────────────────────────────────────────────────────────
def redirect_log(path):
    """Send stdout + stderr of a script rule into its Snakemake log file."""
    fh = open(path, "w", buffering=1)
    sys.stdout = fh
    sys.stderr = fh
    return fh


def log(msg):
    print(msg, flush=True)


# ── 10x / vendor file readers ───────────────────────────────────────────
def _decode(arr):
    """h5py returns fixed-length strings as bytes arrays and variable-length
    ones as object arrays of bytes; normalise both to str."""
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr])


def read_10x_h5(path):
    """Read a 10x Genomics HDF5 feature-barcode matrix (Cell Ranger v3+ /
    Xenium / Space Ranger layout) into an AnnData (cells x features, CSR).

    Kept independent of scanpy so the loader is testable without it.
    """
    import anndata as ad
    import h5py
    import scipy.sparse as sp

    with h5py.File(path, "r") as f:
        if "matrix" not in f:
            raise ValueError(f"{path}: not a v3 10x HDF5 matrix (no 'matrix' group)")
        g = f["matrix"]
        n_features, n_barcodes = (int(x) for x in g["shape"][()])
        data = g["data"][()]
        indices = g["indices"][()]
        indptr = g["indptr"][()]
        barcodes = _decode(g["barcodes"][()])
        feat = g["features"]
        var = pd.DataFrame(
            {
                "gene_ids": _decode(feat["id"][()]),
                "feature_types": _decode(feat["feature_type"][()]),
            },
            index=pd.Index(_decode(feat["name"][()]), name="gene"),
        )
        if "genome" in feat:
            var["genome"] = _decode(feat["genome"][()])
    # 10x stores a CSC matrix of shape (features, barcodes): one column per cell.
    X = sp.csc_matrix((data, indices, indptr), shape=(n_features, n_barcodes)).T.tocsr()
    X = X.astype(np.float32)
    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=pd.Index(barcodes, name="barcode")), var=var)
    adata.var_names_make_unique()
    return adata


def read_cells_table(dirpath):
    """Xenium per-cell table (cells.parquet, else cells.csv.gz / cells.csv),
    indexed by cell_id."""
    for name in ("cells.parquet", "cells.csv.gz", "cells.csv"):
        p = os.path.join(dirpath, name)
        if os.path.exists(p):
            df = pd.read_parquet(p) if name.endswith(".parquet") else pd.read_csv(p)
            if "cell_id" not in df.columns:
                raise ValueError(f"{p}: no cell_id column")
            df["cell_id"] = df["cell_id"].astype(str)
            return df.set_index("cell_id")
    raise FileNotFoundError(f"{dirpath}: no cells.parquet / cells.csv.gz / cells.csv")


POSITION_COLS = ["in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]


def read_positions(spatial_dir):
    """Space Ranger tissue positions, indexed by barcode. Handles the >= 2.0
    `tissue_positions.csv` (with header) and the older
    `tissue_positions_list.csv` (no header); both share the column order."""
    p = os.path.join(spatial_dir, "tissue_positions.csv")
    if os.path.exists(p):
        df = pd.read_csv(p, header=0)
    else:
        p = os.path.join(spatial_dir, "tissue_positions_list.csv")
        if not os.path.exists(p):
            raise FileNotFoundError(f"{spatial_dir}: no tissue_positions.csv / tissue_positions_list.csv")
        df = pd.read_csv(p, header=None)
    if df.shape[1] < 6:
        raise ValueError(f"{p}: expected 6 columns (barcode + {POSITION_COLS}), found {df.shape[1]}")
    df = df.iloc[:, :6]
    df.columns = ["barcode"] + POSITION_COLS
    df["barcode"] = df["barcode"].astype(str)
    return df.set_index("barcode")


def read_scalefactors(spatial_dir):
    with open(os.path.join(spatial_dir, "scalefactors_json.json")) as fh:
        return json.load(fh)


def read_images(spatial_dir):
    """The hires/lowres tissue images when present (arrays), for plotting."""
    images = {}
    for key in ("hires", "lowres"):
        p = os.path.join(spatial_dir, f"tissue_{key}_image.png")
        if os.path.exists(p):
            import matplotlib.image as mpimg

            images[key] = mpimg.imread(p)
    return images


# ── QC threshold logic ──────────────────────────────────────────────────
def apply_thresholds(obs, platform, thresholds):
    """Apply the platform drop rules to a per-cell metrics table.

    Returns a DataFrame (same index) with ``keep`` (bool) and ``drop_reason``
    (the first failing rule's name, '' when kept). A threshold of 0 disables
    its rule; a rule whose metric column is missing is skipped.
    """
    if platform not in DROP_RULES:
        raise ValueError(f"unknown platform {platform!r}; expected one of {PLATFORMS}")
    reason = pd.Series("", index=obs.index, dtype=object)
    for name, col, key, op in DROP_RULES[platform]:
        limit = float(thresholds.get(key, 0) or 0)
        if limit <= 0 or col not in obs.columns:
            continue
        vals = pd.to_numeric(obs[col], errors="coerce")
        fail = (vals < limit) if op == "<" else (vals > limit)
        fail = fail.fillna(False).astype(bool) & (reason == "")
        reason[fail] = name
    return pd.DataFrame({"keep": reason == "", "drop_reason": reason}, index=obs.index)


def drop_reasons(platform):
    return [r[0] for r in DROP_RULES[platform]]


# ── plotting helpers (matplotlib imported lazily) ───────────────────────
def point_size(n):
    """Marker area that keeps a section readable from ~100 to ~1M points."""
    return float(np.clip(20000.0 / max(int(n), 1), 0.3, 40.0))


def categorical_palette(categories):
    import matplotlib

    n = len(categories)
    if n <= 10:
        cmap = matplotlib.colormaps["tab10"]
        colors = [cmap(i) for i in range(n)]
    elif n <= 20:
        cmap = matplotlib.colormaps["tab20"]
        colors = [cmap(i) for i in range(n)]
    else:
        cmap = matplotlib.colormaps["nipy_spectral"]
        colors = [cmap(i / max(n - 1, 1)) for i in range(n)]
    return dict(zip(list(categories), colors))


def spatial_scatter(
    ax,
    xy,
    values=None,
    *,
    title="",
    cmap="viridis",
    categorical=False,
    palette=None,
    size=None,
    invert_y=True,
    colorbar=True,
    legend=True,
    unit="",
):
    """Draw a section on ``ax``: points at ``xy`` coloured by ``values``.

    Image coordinates (both Xenium µm and Visium pixels) grow downwards, so
    the y axis is inverted by default to match the tissue image.
    """
    import matplotlib.pyplot as plt

    xy = np.asarray(xy, dtype=float)
    s = point_size(len(xy)) if size is None else size
    if values is None:
        ax.scatter(xy[:, 0], xy[:, 1], s=s, c="steelblue", linewidths=0, rasterized=True)
    elif categorical:
        cats = pd.Categorical(values)
        pal = palette or categorical_palette(cats.categories)
        for cat in cats.categories:
            m = np.asarray(cats == cat)
            ax.scatter(xy[m, 0], xy[m, 1], s=s, color=pal[cat], label=str(cat), linewidths=0, rasterized=True)
        if legend:
            handles, labels = ax.get_legend_handles_labels()
            if handles:
                leg = ax.legend(
                    handles, labels, loc="center left", bbox_to_anchor=(1.01, 0.5),
                    fontsize=7, frameon=False, markerscale=max(1.0, 8.0 / max(s, 0.3) ** 0.5),
                    ncol=1 if len(labels) <= 20 else 2,
                )
                for h in leg.legend_handles:
                    h.set_sizes([30])
    else:
        vals = np.asarray(values, dtype=float)
        sc = ax.scatter(xy[:, 0], xy[:, 1], c=vals, s=s, cmap=cmap, linewidths=0, rasterized=True)
        if colorbar:
            plt.colorbar(sc, ax=ax, shrink=0.7, pad=0.02)
    ax.set_aspect("equal")
    if invert_y:
        ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    if unit:
        ax.set_xlabel(unit, fontsize=7)
    ax.set_title(title, fontsize=9)


def text_figure(path, message, figsize=(6, 3)):
    """A placeholder PNG carrying a message (e.g. 'fewer than 2 clusters')."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=figsize)
    ax.axis("off")
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=10, wrap=True)
    fig.savefig(path, dpi=100, bbox_inches="tight")
    plt.close(fig)


# ── misc ────────────────────────────────────────────────────────────────
def sanitize_uns(obj):
    """Make a nested dict h5ad-writable: None → '', tuples → lists, numpy
    scalars → Python scalars."""
    if obj is None:
        return ""
    if isinstance(obj, dict):
        return {str(k): sanitize_uns(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [sanitize_uns(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


def write_tsv(df, path, index=True):
    df.to_csv(path, sep="\t", index=index, float_format="%.6g")
