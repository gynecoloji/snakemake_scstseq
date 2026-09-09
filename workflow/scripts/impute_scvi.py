"""Rule `impute_scvi`: scVI denoising of one section (Lopez et al. 2018).

Trains an scvi-tools SCVI model on the raw counts (layers["counts"]) of the
selected genes and exports `get_normalized_expression` at a fixed library
size, log1p-transformed, so the result sits on the same log-normalized scale
as the other methods. The latent space is kept in obsm["X_scVI"].
Runs in the py-scvi (CPU) or py-scvi-gpu env.
"""

import time

import anndata as ad
import numpy as np

from scst_common import log, redirect_log
from scst_impute import dense, finish, select_genes


def library_size(adata, cfg):
    """Counts per cell the imputed matrix is scaled to before log1p: the config
    value, else the platform's normalize_total target, else the median depth."""
    if float(cfg["library_size"]) > 0:
        return float(cfg["library_size"])
    target = float(adata.uns.get("scst", {}).get("processing", {}).get("target_sum", 0) or 0)
    if target > 0:
        return target
    totals = np.asarray(adata.layers["counts"].sum(axis=1)).ravel()
    return float(np.median(totals))


def run(in_h5ad, moran_tsv, outputs, platform, cfg, genes_mode, plot_top, seed, threads, sample):
    import scvi
    import torch

    started = time.time()
    scvi.settings.seed = int(seed)
    torch.set_num_threads(max(1, int(threads)))
    adata = ad.read_h5ad(in_h5ad)
    mask, genes_used = select_genes(adata, genes_mode)
    lib = library_size(adata, cfg)
    sub = adata[:, mask].copy()
    sub.layers["counts"] = sub.layers["counts"].astype(np.float32)
    log(f"scvi: {sub.n_obs} cells, {sub.n_vars} genes ({genes_used}), library size {lib:.1f}")

    scvi.model.SCVI.setup_anndata(sub, layer="counts")
    model = scvi.model.SCVI(sub, n_latent=int(cfg["n_latent"]), n_layers=int(cfg["n_layers"]), n_hidden=int(cfg["n_hidden"]))
    max_epochs = int(cfg["max_epochs"]) or None
    accelerator = "gpu" if cfg["gpu"] else "cpu"
    model.train(max_epochs=max_epochs, accelerator=accelerator, devices=1, early_stopping=False, enable_progress_bar=False)
    epochs = int(len(model.history["elbo_train"])) if "elbo_train" in model.history else (max_epochs or 0)
    imputed = model.get_normalized_expression(library_size=lib, return_numpy=True)
    imputed = np.log1p(np.asarray(imputed, dtype=np.float32))
    latent = np.asarray(model.get_latent_representation(), dtype=np.float32)
    adata.obsm["X_scVI"] = latent
    extras = {"epochs_trained": epochs, "library_size_used": lib, "accelerator": accelerator, "n_latent": int(cfg["n_latent"])}
    return finish(adata, mask, imputed, "scvi", cfg, genes_used, extras, started, outputs, moran_tsv, plot_top, platform, sample)


if "snakemake" in globals():  # pragma: no cover
    redirect_log(snakemake.log[0])  # noqa: F821
    run(
        snakemake.input.h5ad,  # noqa: F821
        snakemake.input.moran,  # noqa: F821
        dict(snakemake.output),  # noqa: F821
        snakemake.params.platform,  # noqa: F821
        dict(snakemake.params.cfg),  # noqa: F821
        snakemake.params.genes,  # noqa: F821
        int(snakemake.params.plot_top),  # noqa: F821
        int(snakemake.params.seed),  # noqa: F821
        int(snakemake.threads),  # noqa: F821
        snakemake.wildcards.sample,  # noqa: F821
    )
