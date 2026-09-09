"""Shared helpers for the imputation rules (MAGIC / ALRA in impute_sample.py,
scVI in impute_scvi.py). Scanpy-free so both conda envs can import it."""

import json
import time

import numpy as np
import pandas as pd

from scst_common import COORD_UNIT, log, sanitize_uns, spatial_scatter, text_figure, write_tsv

ZERO_TOL = 1e-6  # imputed values below this count as zero in the statistics


def select_genes(adata, mode):
    """Boolean mask over var: the HVGs when `mode` is hvg and any exist,
    otherwise every gene. Returns (mask, mode actually used)."""
    if mode == "hvg" and "highly_variable" in adata.var.columns and bool(adata.var["highly_variable"].any()):
        return adata.var["highly_variable"].to_numpy(dtype=bool), "hvg"
    if mode == "hvg":
        log("no HVG column (targeted panel): imputing every gene")
    return np.ones(adata.n_vars, dtype=bool), "all"


def dense(X):
    return X.toarray() if hasattr(X, "toarray") else np.asarray(X)


def gene_stats(observed, imputed, genes):
    """Per gene: zero fraction before/after, means, Pearson r observed vs
    imputed across cells (nan for constant columns)."""
    obs = dense(observed).astype(np.float64)
    imp = np.asarray(imputed, dtype=np.float64)
    obs_c = obs - obs.mean(axis=0)
    imp_c = imp - imp.mean(axis=0)
    denom = np.sqrt((obs_c**2).sum(axis=0) * (imp_c**2).sum(axis=0))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(denom > 0, (obs_c * imp_c).sum(axis=0) / denom, np.nan)
    df = pd.DataFrame(
        {
            "pct_zero_observed": 100.0 * (obs <= ZERO_TOL).mean(axis=0),
            "pct_zero_imputed": 100.0 * (imp <= ZERO_TOL).mean(axis=0),
            "mean_observed": obs.mean(axis=0),
            "mean_imputed": imp.mean(axis=0),
            "pearson_r": r,
        },
        index=pd.Index(genes, name="gene"),
    )
    return df


def plot_observed_vs_imputed(adata, imputed, genes, path, title, platform, method):
    """Top genes on the section: observed (top row) vs imputed (bottom row)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    genes = [g for g in genes if g in adata.var_names]
    if not genes:
        text_figure(path, f"{title}: no genes to draw")
        return
    xy = np.asarray(adata.obsm["spatial"], dtype=float)
    n = len(genes)
    fig, axes = plt.subplots(2, n, figsize=(4.6 * n, 8.6), squeeze=False)
    col = {g: i for i, g in enumerate(adata.var_names)}
    for j, g in enumerate(genes):
        obs = dense(adata[:, g].X).ravel()
        imp = np.asarray(imputed[:, col[g]]).ravel()
        vmax = float(max(np.nanmax(obs), np.nanmax(imp), 1e-9))
        for ax, vals, label in ((axes[0, j], obs, "observed"), (axes[1, j], imp, f"{method} imputed")):
            spatial_scatter(ax, xy, np.clip(vals, 0, vmax), title=f"{g} — {label}", cmap="magma", unit=COORD_UNIT[platform])
    fig.suptitle(f"{title} — observed vs {method}-imputed (log-normalized)", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def top_genes(moran_tsv, candidates, n):
    """The top-n Moran's I genes among `candidates` (the imputed genes)."""
    if n <= 0:
        return []
    df = pd.read_csv(moran_tsv, sep="\t")
    cand = set(candidates)
    return [g for g in df.sort_values("rank")["gene"] if g in cand][: int(n)]


def finish(adata, mask, imputed, method, cfg, genes_mode, extras, started, outputs, moran_tsv, plot_top, platform, sample):
    """Common tail of every imputation rule: restrict the object to the
    imputed genes, attach layers["imputed"], write stats / run JSON / figure."""
    out = adata[:, mask].copy()
    imputed = np.asarray(imputed, dtype=np.float32)
    if imputed.shape != (out.n_obs, out.n_vars):
        raise ValueError(f"imputed matrix has shape {imputed.shape}, expected {(out.n_obs, out.n_vars)}")
    out.layers["imputed"] = imputed
    stats = gene_stats(out.X, imputed, out.var_names)
    write_tsv(stats, outputs["stats"])
    runtime = time.time() - started
    run = {
        "sample_id": sample,
        "platform": platform,
        "method": method,
        "genes_mode": genes_mode,
        "n_cells": int(out.n_obs),
        "n_genes_imputed": int(out.n_vars),
        "pct_zero_observed_mean": float(stats["pct_zero_observed"].mean()),
        "pct_zero_imputed_mean": float(stats["pct_zero_imputed"].mean()),
        "median_gene_pearson_r": float(np.nanmedian(stats["pearson_r"])) if stats["pearson_r"].notna().any() else float("nan"),
        "runtime_s": round(runtime, 2),
        "params": dict(cfg),
        **extras,
    }
    with open(outputs["run"], "w") as fh:
        json.dump(run, fh, indent=2)
    genes = top_genes(moran_tsv, out.var_names, plot_top)
    plot_observed_vs_imputed(out, imputed, genes, outputs["png"], f"{sample} ({platform})", platform, method)
    out.uns["scst"] = {**out.uns.get("scst", {}), "imputation": sanitize_uns({k: v for k, v in run.items() if k != "params"} | {"params": dict(cfg)})}
    out.write_h5ad(outputs["h5ad"])
    log(
        f"{method}: {out.n_obs} cells x {out.n_vars} genes in {runtime:.1f}s; zero fraction "
        f"{run['pct_zero_observed_mean']:.1f}% -> {run['pct_zero_imputed_mean']:.1f}%, median r {run['median_gene_pearson_r']:.3f}"
    )
    return run
