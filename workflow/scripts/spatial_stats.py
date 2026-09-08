"""Rule `spatial_stats`: spatial neighbours graph + neighborhood enrichment +
Moran's I spatially variable genes for one section (squidpy).

Graph:  Xenium — generic coordinates (µm): kNN (`n_neighs`), Delaunay, or a
        fixed `radius`; Visium — the hexagonal spot grid (`n_neighs` per ring,
        `n_rings`).
Stats:  nhood_enrichment on the Leiden clusters (permutation z-scores);
        spatial_autocorr (Moran's I) over the HVGs when they exist, else all
        genes, with permutation p-values (n_perms = 0 → analytic only) and
        BH-FDR; the top SVGs drawn on the section.
"""


import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import seaborn as sns  # noqa: E402
import squidpy as sq  # noqa: E402

from scst_common import (  # noqa: E402
    COORD_UNIT,
    XENIUM,
    log,
    redirect_log,
    sanitize_uns,
    spatial_scatter,
    text_figure,
    write_tsv,
)

CLUSTER_KEY = "leiden"
NHOOD_KEY = f"{CLUSTER_KEY}_nhood_enrichment"


def graph_kwargs(platform, g):
    if platform == XENIUM:
        kw = {"coord_type": "generic", "delaunay": bool(g["delaunay"]), "n_neighs": int(g["n_neighs"])}
        if float(g["radius"]) > 0:
            kw["radius"] = float(g["radius"])
        return kw
    return {"coord_type": "grid", "n_neighs": int(g["n_neighs"]), "n_rings": int(g["n_rings"])}


def build_graph(adata, platform, g):
    """Build obsp['spatial_connectivities'] with the platform block of
    config.spatial. Uses squidpy's per-method builders (1.8+) when present;
    the classic `spatial_neighbors` (deprecated in 1.8, removed in 1.9) is the
    fallback for older releases."""
    kw = graph_kwargs(platform, g)
    if hasattr(sq.gr, "spatial_neighbors_grid"):
        if kw["coord_type"] == "grid":
            sq.gr.spatial_neighbors_grid(adata, spatial_key="spatial", n_neighs=kw["n_neighs"], n_rings=kw["n_rings"])
        elif "radius" in kw:
            sq.gr.spatial_neighbors_radius(adata, spatial_key="spatial", radius=kw["radius"])
        elif kw["delaunay"]:
            sq.gr.spatial_neighbors_delaunay(adata, spatial_key="spatial")
        else:
            sq.gr.spatial_neighbors_knn(adata, spatial_key="spatial", n_neighs=kw["n_neighs"])
    else:  # pragma: no cover - squidpy < 1.8
        sq.gr.spatial_neighbors(adata, spatial_key="spatial", **kw)
    conn = adata.obsp["spatial_connectivities"]
    degree = np.asarray(conn.sum(axis=1)).ravel()
    log(f"spatial graph {kw}: mean degree {degree.mean():.2f}, {int((degree == 0).sum())} isolated")
    return kw


def nhood_enrichment(adata, n_perms, seed, n_jobs, out_tsv, out_png, title):
    cats = list(adata.obs[CLUSTER_KEY].cat.categories)
    if len(cats) < 2:
        log(f"{len(cats)} cluster(s): neighborhood enrichment needs >= 2; writing empty outputs")
        write_tsv(pd.DataFrame(index=pd.Index([], name="cluster")), out_tsv)
        text_figure(out_png, f"{title}: only {len(cats)} cluster — no neighborhood enrichment")
        return None
    sq.gr.nhood_enrichment(
        adata, cluster_key=CLUSTER_KEY, n_perms=int(n_perms), seed=seed, n_jobs=n_jobs, show_progress_bar=False,
    )
    z = pd.DataFrame(adata.uns[NHOOD_KEY]["zscore"], index=cats, columns=cats)
    z.index.name = "cluster"
    write_tsv(z, out_tsv)
    fig, ax = plt.subplots(figsize=(max(4, 0.45 * len(cats) + 2.5), max(3.5, 0.45 * len(cats) + 2)))
    vmax = float(np.nanmax(np.abs(z.values))) if np.isfinite(z.values).any() else 1.0
    sns.heatmap(
        z, ax=ax, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax, annot=len(cats) <= 15, fmt=".1f",
        annot_kws={"size": 7}, cbar_kws={"label": "z-score", "shrink": 0.8}, square=True,
    )
    ax.set_title(f"{title} — neighborhood enrichment", fontsize=10)
    ax.set_xlabel("cluster")
    ax.set_ylabel("cluster")
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)
    return z


def moran_genes(adata):
    if "highly_variable" in adata.var.columns and bool(adata.var["highly_variable"].any()):
        genes = adata.var_names[adata.var["highly_variable"].to_numpy(dtype=bool)]
        log(f"Moran's I over {len(genes)} HVGs")
    else:
        genes = adata.var_names
        log(f"Moran's I over all {len(genes)} genes")
    return list(genes)


def morans_i(adata, cfg, seed, n_jobs, out_tsv):
    genes = moran_genes(adata)
    n_perms = int(cfg["n_perms"]) or None
    sq.gr.spatial_autocorr(
        adata, mode="moran", genes=genes, n_perms=n_perms, n_jobs=n_jobs, seed=seed, show_progress_bar=False,
    )
    df = adata.uns["moranI"].copy()
    df.index.name = "gene"
    df = df.sort_values("I", ascending=False)
    df.insert(0, "rank", np.arange(1, len(df) + 1))
    write_tsv(df, out_tsv)
    log(f"Moran's I: top gene {df.index[0]} (I = {df['I'].iloc[0]:.3f})")
    return df


def plot_top_svg(adata, df, n, platform, out_png, title):
    n = int(min(n, len(df)))
    if n <= 0:
        text_figure(out_png, f"{title}: no SVG panels requested")
        return
    genes = list(df.index[:n])
    ncols = 3 if n >= 3 else n
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 4.6 * nrows), squeeze=False)
    xy = np.asarray(adata.obsm["spatial"], dtype=float)
    for ax, gene in zip(axes.ravel(), genes):
        vals = adata[:, gene].X
        vals = vals.toarray().ravel() if hasattr(vals, "toarray") else np.asarray(vals).ravel()
        spatial_scatter(
            ax, xy, vals, title=f"{gene}  (Moran's I = {df.loc[gene, 'I']:.2f})", cmap="magma", unit=COORD_UNIT[platform],
        )
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.suptitle(f"{title} — top spatially variable genes (log-normalized)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)


def run(in_h5ad, outputs, platform, graph, nhood, moran, seed, threads, sample):
    adata = sc.read_h5ad(in_h5ad)
    adata.obs[CLUSTER_KEY] = adata.obs[CLUSTER_KEY].astype(str).astype("category")
    title = f"{sample} ({platform})"
    kw = build_graph(adata, platform, graph)
    nhood_enrichment(adata, nhood["n_perms"], seed, threads, outputs["nhood"], outputs["nhood_png"], title)
    df = morans_i(adata, moran, seed, threads, outputs["moran"])
    plot_top_svg(adata, df, moran["plot_top"], platform, outputs["svg_png"], title)
    adata.uns["scst"] = {
        **adata.uns.get("scst", {}),
        "spatial": sanitize_uns(
            {"graph": kw, "n_clusters": len(adata.obs[CLUSTER_KEY].cat.categories), "moran_genes": len(df), "moran_n_perms": moran["n_perms"], "nhood_n_perms": nhood["n_perms"]}
        ),
    }
    if "spatial_neighbors" in adata.uns:
        adata.uns["spatial_neighbors"] = sanitize_uns(dict(adata.uns["spatial_neighbors"]))
    adata.write_h5ad(outputs["h5ad"])
    log(f"wrote {outputs['h5ad']}")


if "snakemake" in globals():  # pragma: no cover
    redirect_log(snakemake.log[0])  # noqa: F821
    run(
        snakemake.input.h5ad,  # noqa: F821
        dict(snakemake.output),  # noqa: F821
        snakemake.params.platform,  # noqa: F821
        dict(snakemake.params.graph),  # noqa: F821
        dict(snakemake.params.nhood),  # noqa: F821
        dict(snakemake.params.moran),  # noqa: F821
        int(snakemake.params.seed),  # noqa: F821
        int(snakemake.threads),  # noqa: F821
        snakemake.wildcards.sample,  # noqa: F821
    )
