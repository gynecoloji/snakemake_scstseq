"""Rule `cluster_markers`: per-cluster marker genes (scanpy rank_genes_groups
on the log-normalized X) → top-N table + a dotplot of the leading markers."""


import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402

from scst_common import log, redirect_log, text_figure, write_tsv  # noqa: E402

COLUMNS = ["cluster", "rank", "gene", "score", "log2fc", "pval", "padj", "pct_in", "pct_out"]
RENAME = {
    "group": "cluster",
    "names": "gene",
    "scores": "score",
    "logfoldchanges": "log2fc",
    "pvals": "pval",
    "pvals_adj": "padj",
    "pct_nz_group": "pct_in",
    "pct_nz_reference": "pct_out",
}


def marker_table(adata, method, n_top, key="leiden"):
    sc.tl.rank_genes_groups(adata, groupby=key, method=method, use_raw=False, pts=True)
    df = sc.get.rank_genes_groups_df(adata, group=None).rename(columns=RENAME)
    df["rank"] = df.groupby("cluster", observed=True).cumcount() + 1
    df = df[df["rank"] <= int(n_top)]
    for col in COLUMNS:
        if col not in df.columns:
            df[col] = float("nan")
    return df[COLUMNS].reset_index(drop=True)


def leading_genes(df, per_cluster=3):
    genes = []
    for _, sub in df.groupby("cluster", sort=False, observed=True):
        for g in sub.sort_values("rank")["gene"].head(per_cluster):
            if g not in genes:
                genes.append(g)
    return genes


def dotplot(adata, genes, path, key="leiden"):
    try:
        dp = sc.pl.dotplot(adata, var_names=genes, groupby=key, use_raw=False, return_fig=True, show=False)
        dp.savefig(path, bbox_inches="tight", dpi=110)
        plt.close("all")
    except Exception as exc:  # pragma: no cover - plotting fallback
        log(f"dotplot failed ({exc}); writing a placeholder")
        text_figure(path, f"dotplot unavailable: {exc}")


def run(in_h5ad, out_tsv, out_png, cfg, sample):
    adata = sc.read_h5ad(in_h5ad)
    adata.obs["leiden"] = adata.obs["leiden"].astype(str).astype("category")
    n_clusters = len(adata.obs["leiden"].cat.categories)
    if n_clusters < 2:
        log(f"{sample}: {n_clusters} cluster(s) — marker detection needs >= 2; writing empty outputs")
        write_tsv(pd.DataFrame(columns=COLUMNS), out_tsv, index=False)
        text_figure(out_png, f"{sample}: only {n_clusters} cluster — no markers")
        return
    df = marker_table(adata, cfg["method"], cfg["n_top"])
    write_tsv(df, out_tsv, index=False)
    genes = leading_genes(df)
    log(f"{sample}: {n_clusters} clusters, {len(df)} marker rows, dotplot of {len(genes)} genes")
    dotplot(adata, genes, out_png)


if "snakemake" in globals():  # pragma: no cover
    redirect_log(snakemake.log[0])  # noqa: F821
    run(
        snakemake.input.h5ad,  # noqa: F821
        snakemake.output.tsv,  # noqa: F821
        snakemake.output.dotplot,  # noqa: F821
        dict(snakemake.params.cfg),  # noqa: F821
        snakemake.wildcards.sample,  # noqa: F821
    )
