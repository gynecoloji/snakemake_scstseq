"""Rule `qc_filter`: per-cell/spot metrics, platform thresholds, gene filter.

Metrics: total_counts, n_genes_by_counts, pct_counts_mt / _ribo / _hb (from
the gene-name regexes in config.genome). Xenium cells additionally carry
cell_area, nucleus_area and control_frac from the loader.

Thresholds come from config.qc.<platform>; 0 disables a rule; the FIRST rule
a cell fails is recorded as its drop reason (see scst_common.DROP_RULES).
Genes detected in fewer than `min_cells` kept cells are removed afterwards.
"""


import json
import re

import anndata as ad
import numpy as np
import pandas as pd

from scst_common import (
    UNIT_NAME,
    XENIUM,
    apply_thresholds,
    drop_reasons,
    log,
    redirect_log,
    write_tsv,
)

QC_PATTERNS = (("mt", "mito_pattern"), ("ribo", "ribo_pattern"), ("hb", "hb_pattern"))


def compute_metrics(adata, patterns):
    """Add total_counts, n_genes_by_counts and pct_counts_<subset> to obs
    (and the boolean subset flags to var). Works on any sparse/dense X."""
    X = adata.X
    total = np.asarray(X.sum(axis=1)).ravel()
    n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
    adata.obs["total_counts"] = total
    adata.obs["n_genes_by_counts"] = n_genes.astype(int)
    names = pd.Index(adata.var_names.astype(str))
    for key, cfg_key in QC_PATTERNS:
        pat = patterns.get(cfg_key, "") or ""
        mask = names.str.contains(pat, regex=True) if pat else pd.Series(False, index=names)
        mask = np.asarray(mask, dtype=bool)
        adata.var[key] = mask
        sub = np.asarray(X[:, mask].sum(axis=1)).ravel() if mask.any() else np.zeros(adata.n_obs)
        adata.obs[f"pct_counts_{key}"] = np.where(total > 0, 100.0 * sub / np.maximum(total, 1e-12), 0.0)
    return adata


def filter_genes(adata, min_cells):
    """Keep genes detected in >= min_cells cells (0 = keep all). Returns the
    number of genes removed."""
    if not min_cells or min_cells <= 0:
        return 0
    detected = np.asarray((adata.X > 0).sum(axis=0)).ravel()
    keep = detected >= int(min_cells)
    removed = int((~keep).sum())
    if keep.sum() == 0:
        raise ValueError(f"gene filter min_cells={min_cells} removed every gene")
    if removed:
        adata._inplace_subset_var(keep)
    return removed


CELL_TABLE_COLS = [
    "sample_id",
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt",
    "pct_counts_ribo",
    "pct_counts_hb",
    "cell_area",
    "nucleus_area",
    "transcript_counts",
    "control_counts",
    "control_frac",
    "in_tissue",
    "array_row",
    "array_col",
]


def cell_table(adata, verdict):
    df = pd.DataFrame(index=adata.obs_names)
    xy = np.asarray(adata.obsm["spatial"], dtype=float)
    df["x"] = xy[:, 0]
    df["y"] = xy[:, 1]
    for col in CELL_TABLE_COLS:
        if col in adata.obs.columns:
            df[col] = adata.obs[col].values
    df["qc_keep"] = verdict["keep"].values
    df["qc_drop_reason"] = verdict["drop_reason"].values
    df.index.name = "cell"
    return df


def summarize(adata_raw, verdict, adata_kept, platform, thresholds, n_genes_removed):
    obs = adata_raw.obs
    kept = verdict["keep"].values
    stats = {
        "sample_id": str(adata_raw.uns["scst"].get("sample_id", "")),
        "platform": platform,
        "condition": str(adata_raw.uns["scst"].get("condition", "")),
        "n_cells_raw": int(adata_raw.n_obs),
        "n_cells_kept": int(kept.sum()),
        "n_cells_dropped": int((~kept).sum()),
        "pct_kept": float(100.0 * kept.sum() / max(adata_raw.n_obs, 1)),
    }
    counts = verdict["drop_reason"].value_counts()
    for reason in drop_reasons(platform):
        stats[f"drop_{reason}"] = int(counts.get(reason, 0))
    stats.update(
        {
            "n_genes_raw": int(adata_raw.n_vars),
            "n_genes_kept": int(adata_kept.n_vars),
            "n_genes_removed_min_cells": int(n_genes_removed),
            "median_counts_raw": float(np.median(obs["total_counts"])),
            "median_counts_kept": float(np.median(obs["total_counts"][kept])),
            "median_genes_raw": float(np.median(obs["n_genes_by_counts"])),
            "median_genes_kept": float(np.median(obs["n_genes_by_counts"][kept])),
            "median_pct_mt_kept": float(np.median(obs["pct_counts_mt"][kept])),
        }
    )
    if platform == XENIUM:
        if "cell_area" in obs:
            stats["median_cell_area_kept"] = float(np.median(obs["cell_area"][kept]))
        if "control_frac" in obs:
            stats["median_control_frac_kept"] = float(np.median(obs["control_frac"][kept]))
    stats["thresholds"] = dict(thresholds)
    return stats


def run(in_h5ad, out_h5ad, out_cells, out_stats, platform, thresholds, patterns):
    adata = ad.read_h5ad(in_h5ad)
    compute_metrics(adata, patterns)
    verdict = apply_thresholds(adata.obs, platform, thresholds)
    adata.obs["qc_keep"] = verdict["keep"].values
    adata.obs["qc_drop_reason"] = verdict["drop_reason"].values
    write_tsv(cell_table(adata, verdict), out_cells)

    unit = UNIT_NAME[platform]
    n_keep = int(verdict["keep"].sum())
    log(f"{platform}: {adata.n_obs} {unit}s, {n_keep} kept, {adata.n_obs - n_keep} dropped")
    log("drop reasons: " + json.dumps(verdict["drop_reason"].value_counts().to_dict()))
    if n_keep == 0:
        raise SystemExit(
            f"qc_filter: every {unit} of {adata.uns['scst'].get('sample_id')} was dropped by the "
            f"{platform} thresholds {json.dumps(thresholds)} — relax them in config.qc.{platform}"
        )
    kept = adata[verdict["keep"].values].copy()
    n_removed = filter_genes(kept, thresholds.get("min_cells", 0))
    log(f"gene filter (min_cells={thresholds.get('min_cells', 0)}): removed {n_removed}, kept {kept.n_vars}")

    stats = summarize(adata, verdict, kept, platform, thresholds, n_removed)
    with open(out_stats, "w") as fh:
        json.dump(stats, fh, indent=2)
    kept.uns["scst"] = {**kept.uns.get("scst", {}), "qc": {"thresholds": dict(thresholds), "n_cells_raw": int(adata.n_obs)}}
    kept.write_h5ad(out_h5ad)
    log(f"wrote {kept.n_obs} x {kept.n_vars} to {out_h5ad}")
    return stats


if "snakemake" in globals():  # pragma: no cover
    redirect_log(snakemake.log[0])  # noqa: F821
    run(
        snakemake.input.h5ad,  # noqa: F821
        snakemake.output.h5ad,  # noqa: F821
        snakemake.output.cells,  # noqa: F821
        snakemake.output.stats,  # noqa: F821
        snakemake.params.platform,  # noqa: F821
        dict(snakemake.params.thresholds),  # noqa: F821
        dict(snakemake.params.patterns),  # noqa: F821
    )
