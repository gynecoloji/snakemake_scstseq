"""Rule `process_sample`: normalization → (HVG) → PCA → neighbours → Leiden →
UMAP for one section, with the platform block of config.processing.

X ends up log-normalized (used by the marker and SVG steps); the raw counts
are kept in layers["counts"]. Scaling, when enabled, happens on a temporary
copy so X is not overwritten. n_pcs and n_neighbors are capped to the data
size so tiny sections (and the executable test case) still run.
"""


import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402

from scst_common import (  # noqa: E402
    COORD_UNIT,
    categorical_palette,
    log,
    redirect_log,
    sanitize_uns,
    spatial_scatter,
    write_tsv,
)


def select_hvg(adata, n_top):
    """seurat_v3 on the counts layer; falls back to the log-data 'seurat'
    flavour when the loess fit cannot be made (very small panels)."""
    try:
        sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_top, layer="counts")
    except Exception as exc:  # pragma: no cover - depends on data size
        log(f"seurat_v3 HVG failed ({exc}); falling back to flavor='seurat' on log data")
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=n_top)
    return adata.var["highly_variable"].to_numpy(dtype=bool)


def process(adata, cfg, seed, threads=1):
    sc.settings.n_jobs = int(threads)
    adata.layers["counts"] = adata.X.copy()
    target = float(cfg["target_sum"]) if cfg["target_sum"] else None
    sc.pp.normalize_total(adata, target_sum=target)
    sc.pp.log1p(adata)

    n_top = int(cfg["n_top_genes"])
    use_hvg = bool(cfg["hvg"]) and 0 < n_top < adata.n_vars
    if use_hvg:
        mask = select_hvg(adata, n_top)
        log(f"HVG: {int(mask.sum())} of {adata.n_vars} genes")
    else:
        mask = np.ones(adata.n_vars, dtype=bool)
        if cfg["hvg"]:
            log(f"HVG requested but n_top_genes={n_top} >= n_genes={adata.n_vars}: using every gene")
        else:
            log(f"using every gene ({adata.n_vars}) for PCA")

    sub = adata[:, mask].copy()
    if cfg["scale"]:
        sc.pp.scale(sub, max_value=10)
    n_pcs = int(max(2, min(int(cfg["n_pcs"]), sub.n_vars - 1, sub.n_obs - 1)))
    sc.tl.pca(sub, n_comps=n_pcs, random_state=seed)
    adata.obsm["X_pca"] = sub.obsm["X_pca"]
    adata.uns["pca"] = {
        "variance_ratio": np.asarray(sub.uns["pca"]["variance_ratio"]),
        "variance": np.asarray(sub.uns["pca"]["variance"]),
    }

    n_neighbors = int(max(2, min(int(cfg["n_neighbors"]), adata.n_obs - 1)))
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep="X_pca", random_state=seed)
    sc.tl.leiden(
        adata,
        resolution=float(cfg["resolution"]),
        random_state=seed,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        key_added="leiden",
    )
    sc.tl.umap(adata, random_state=seed)
    n_clusters = adata.obs["leiden"].nunique()
    log(f"PCA {n_pcs} PCs, neighbours {n_neighbors}, Leiden res {cfg['resolution']} → {n_clusters} clusters")
    adata.uns["scst"] = {
        **adata.uns.get("scst", {}),
        "processing": sanitize_uns(
            {**cfg, "n_pcs_used": n_pcs, "n_neighbors_used": n_neighbors, "n_genes_for_pca": int(mask.sum()), "n_clusters": int(n_clusters)}
        ),
    }
    return adata


def cluster_table(adata):
    xy = np.asarray(adata.obsm["spatial"], dtype=float)
    umap = np.asarray(adata.obsm["X_umap"], dtype=float)
    df = pd.DataFrame(
        {
            "leiden": adata.obs["leiden"].astype(str).values,
            "umap_1": umap[:, 0],
            "umap_2": umap[:, 1],
            "x": xy[:, 0],
            "y": xy[:, 1],
        },
        index=adata.obs_names,
    )
    df.index.name = "cell"
    return df


def plot_clusters(adata, out_umap, out_spatial, title, platform):
    cats = adata.obs["leiden"].astype("category")
    palette = categorical_palette(cats.cat.categories)
    umap = np.asarray(adata.obsm["X_umap"], dtype=float)
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    spatial_scatter(ax, umap, cats, title=f"{title} — UMAP", categorical=True, palette=palette, invert_y=False)
    ax.set_xlabel("UMAP 1", fontsize=8)
    ax.set_ylabel("UMAP 2", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_umap, dpi=110, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    spatial_scatter(
        ax, adata.obsm["spatial"], cats, title=f"{title} — Leiden clusters on the section",
        categorical=True, palette=palette, unit=COORD_UNIT[platform],
    )
    fig.tight_layout()
    fig.savefig(out_spatial, dpi=110, bbox_inches="tight")
    plt.close(fig)


def run(in_h5ad, out_h5ad, out_clusters, out_umap, out_spatial, platform, cfg, seed, threads, sample):
    adata = sc.read_h5ad(in_h5ad)
    process(adata, cfg, seed, threads)
    write_tsv(cluster_table(adata), out_clusters)
    plot_clusters(adata, out_umap, out_spatial, f"{sample} ({platform})", platform)
    adata.write_h5ad(out_h5ad)
    log(f"wrote {adata.n_obs} x {adata.n_vars} to {out_h5ad}")


if "snakemake" in globals():  # pragma: no cover
    redirect_log(snakemake.log[0])  # noqa: F821
    run(
        snakemake.input.h5ad,  # noqa: F821
        snakemake.output.h5ad,  # noqa: F821
        snakemake.output.clusters,  # noqa: F821
        snakemake.output.umap,  # noqa: F821
        snakemake.output.spatial,  # noqa: F821
        snakemake.params.platform,  # noqa: F821
        dict(snakemake.params.cfg),  # noqa: F821
        int(snakemake.params.seed),  # noqa: F821
        int(snakemake.threads),  # noqa: F821
        snakemake.wildcards.sample,  # noqa: F821
    )
