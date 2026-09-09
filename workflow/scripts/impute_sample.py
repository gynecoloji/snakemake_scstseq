"""Rule `impute` (methods magic | alra): impute the dropouts of the measured
genes of one section on the log-normalized X.

magic  magic-impute: kNN graph on (PCA of) the log-normalized values, alpha-
       decay kernel, diffusion for `t` steps (van Dijk et al. 2018).
alra   the numpy re-implementation in alra.py (Linderman et al. 2022).
"""

import time

import anndata as ad
import numpy as np

from scst_common import log, redirect_log
from scst_impute import dense, finish, select_genes


def run_magic(X, cfg, n_cells, n_genes, threads, seed):
    import magic

    n_pca = None if n_genes <= 100 else int(min(cfg["n_pca"], n_genes - 1, n_cells - 1))
    knn = int(min(cfg["knn"], n_cells - 1))
    op = magic.MAGIC(
        knn=knn, t=int(cfg["t"]), n_pca=n_pca, decay=float(cfg["decay"]), n_jobs=int(threads),
        random_state=seed, verbose=0,
    )
    imputed = np.asarray(op.fit_transform(X), dtype=np.float32)
    return imputed, {"knn_used": knn, "n_pca_used": n_pca if n_pca is not None else 0}


def run_alra(X, cfg, seed):
    from alra import alra

    imputed, k = alra(X, k=int(cfg["k"]), quantile=float(cfg["quantile"]), k_max=int(cfg["k_max"]), seed=seed)
    return imputed, {"rank_used": int(k)}


def run(in_h5ad, moran_tsv, outputs, method, platform, cfg, genes_mode, plot_top, seed, threads, sample):
    started = time.time()
    adata = ad.read_h5ad(in_h5ad)
    mask, genes_used = select_genes(adata, genes_mode)
    X = adata[:, mask].X
    log(f"{method}: {adata.n_obs} cells, {int(mask.sum())} genes ({genes_used})")
    if method == "magic":
        imputed, extras = run_magic(X, cfg, adata.n_obs, int(mask.sum()), threads, seed)
    elif method == "alra":
        imputed, extras = run_alra(X, cfg, seed)
    else:
        raise ValueError(f"impute_sample.py handles magic | alra, not {method!r}")
    imputed = np.clip(imputed, 0, None)
    return finish(adata, mask, imputed, method, cfg, genes_used, extras, started, outputs, moran_tsv, plot_top, platform, sample)


if "snakemake" in globals():  # pragma: no cover
    redirect_log(snakemake.log[0])  # noqa: F821
    run(
        snakemake.input.h5ad,  # noqa: F821
        snakemake.input.moran,  # noqa: F821
        dict(snakemake.output),  # noqa: F821
        snakemake.wildcards.method,  # noqa: F821
        snakemake.params.platform,  # noqa: F821
        dict(snakemake.params.cfg),  # noqa: F821
        snakemake.params.genes,  # noqa: F821
        int(snakemake.params.plot_top),  # noqa: F821
        int(snakemake.params.seed),  # noqa: F821
        int(snakemake.threads),  # noqa: F821
        snakemake.wildcards.sample,  # noqa: F821
    )
