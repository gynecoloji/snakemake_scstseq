#!/usr/bin/env python3
"""Assert the NUMERIC outputs of the executable test case.

Run after `snakemake -d .test` (the default target). Exits non-zero on any
failed assertion.

The expected values are known BY CONSTRUCTION from .test/make_testdata.py —
how many low-quality cells and high-mito spots were planted, which genes are
domain markers, which cells belong to which domain — not recorded from a
previous run. File existence would not catch a QC rule that silently drops
nothing, a cluster that mixes the two domains, or a Moran's I that ranks
housekeeping genes first; this does.

Only the standard library, numpy, h5py and pandas are needed, so it runs in
the driver environment (no scanpy).
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from make_testdata import N_GENES, SAMPLES, expected_kept, n_cells  # noqa: E402

TEST = Path(__file__).resolve().parent
RESULTS = TEST / "results"


def imputation_methods():
    """The methods listed in .test/config/config.yaml (stdlib parse: the
    `methods:` line inside the imputation block)."""
    for line in (TEST / "config" / "config.yaml").read_text().splitlines():
        if line.strip().startswith("methods:"):
            inner = line.split("[", 1)[1].split("]", 1)[0]
            return [m.strip().strip('"').strip("'") for m in inner.split(",") if m.strip()]
    return []

FAILURES = []


def check(cond, msg):
    status = "ok  " if cond else "FAIL"
    print(f"[{status}] {msg}")
    if not cond:
        FAILURES.append(msg)


def truth(sample):
    df = pd.read_csv(TEST / "truth" / f"{sample}_truth.tsv", sep="\t", dtype=str)
    df["cell"] = sample + "_" + df["barcode"]
    return df.set_index("cell")


def check_summary():
    df = pd.read_csv(RESULTS / "qc" / "qc_summary_all.tsv", sep="\t")
    check(len(df) == len(SAMPLES), f"summary has exactly {len(SAMPLES)} rows (got {len(df)})")
    check(set(df["sample_id"]) == set(SAMPLES), "summary lists every sample once")
    for s, spec in SAMPLES.items():
        row = df[df["sample_id"] == s]
        if row.empty:
            check(False, f"{s}: missing from summary")
            continue
        row = row.iloc[0]
        check(row["platform"] == spec["platform"], f"{s}: platform {row['platform']}")
        check(int(row["n_cells_raw"]) == n_cells(s), f"{s}: n_cells_raw == {n_cells(s)} (got {row['n_cells_raw']})")
        check(int(row["n_cells_kept"]) == expected_kept(s), f"{s}: n_cells_kept == {expected_kept(s)} (got {row['n_cells_kept']})")
        check(int(row["drop_low_counts"]) == spec["n_low"], f"{s}: drop_low_counts == {spec['n_low']} (got {row['drop_low_counts']})")
        if spec["platform"] == "visium":
            check(int(row["drop_high_mito"]) == spec["n_mito"], f"{s}: drop_high_mito == {spec['n_mito']} (got {row['drop_high_mito']})")
        check(int(row["n_genes_raw"]) == N_GENES, f"{s}: {N_GENES} genes after loading (controls removed; got {row['n_genes_raw']})")


def check_qc_cells(sample):
    df = pd.read_csv(RESULTS / "qc" / f"{sample}_qc_cells.tsv", sep="\t", index_col=0)
    t = truth(sample)
    df = df.join(t, how="left")
    check(df["quality"].notna().all(), f"{sample}: every QC row maps to a planted cell")
    low = df[df["quality"] == "low"]
    check((low["qc_drop_reason"] == "low_counts").all(), f"{sample}: all {len(low)} planted low-quality cells dropped as low_counts")
    normal = df[df["quality"] == "normal"]
    check(normal["qc_keep"].all(), f"{sample}: all {len(normal)} normal cells kept")
    if SAMPLES[sample]["platform"] == "visium":
        mito = df[df["quality"] == "mito"]
        check((mito["qc_drop_reason"] == "high_mito").all(), f"{sample}: all {len(mito)} planted high-mito spots dropped as high_mito")
        check((mito["pct_counts_mt"] > 40).all(), f"{sample}: high-mito spots have > 40 % mito")
    else:
        kept = df[df["qc_keep"].astype(bool)]
        check("control_frac" in df.columns and (kept["control_frac"] < 0.05).all(), f"{sample}: control fraction < 5 % in every kept cell (max {kept['control_frac'].max():.3f})")
        check("cell_area" in df.columns, f"{sample}: cell_area carried from cells.csv.gz")


def check_clusters(sample):
    df = pd.read_csv(RESULTS / "processed" / f"{sample}_clusters.tsv", sep="\t", index_col=0, dtype={"leiden": str})
    t = truth(sample)
    df = df.join(t[["domain"]], how="left")
    check(len(df) == expected_kept(sample), f"{sample}: clusters table has {expected_kept(sample)} kept cells (got {len(df)})")
    n_clusters = df["leiden"].nunique()
    check(n_clusters >= 2, f"{sample}: >= 2 Leiden clusters (got {n_clusters})")
    majority = df.groupby("leiden")["domain"].agg(lambda d: d.value_counts().idxmax())
    df["majority"] = df["leiden"].map(majority)
    agreement = float((df["majority"] == df["domain"]).mean())
    check(agreement >= 0.9, f"{sample}: {agreement:.1%} of cells sit in a cluster dominated by their own domain (>= 90 %)")
    check(set(majority) == {"A", "B"}, f"{sample}: both domains are recovered by some cluster (got {sorted(set(majority))})")
    return majority


def check_markers(sample, majority):
    df = pd.read_csv(RESULTS / "processed" / f"{sample}_markers.tsv", sep="\t", dtype={"cluster": str})
    check(len(df) > 0, f"{sample}: marker table is non-empty")
    for cluster, dom in majority.items():
        top = df[(df["cluster"] == cluster) & (df["rank"] <= 5)]["gene"].tolist()
        hits = sum(g.startswith(f"DOM{dom}_") for g in top)
        check(hits >= 3, f"{sample}: cluster {cluster} (domain {dom}) has >= 3 DOM{dom} genes in its top 5 markers ({top})")


def check_spatial(sample):
    moran = pd.read_csv(RESULTS / "spatial" / sample / "morans_i.tsv", sep="\t")
    top10 = moran.sort_values("rank")["gene"].head(10).tolist()
    check(all(g.startswith(("DOMA_", "DOMB_")) for g in top10), f"{sample}: top-10 Moran's I genes are all domain markers ({top10})")
    check(not any(g.startswith("HK_") for g in top10), f"{sample}: no housekeeping gene in the top-10 SVGs")
    check(float(moran["I"].max()) > 0.5, f"{sample}: best Moran's I > 0.5 (got {moran['I'].max():.3f})")
    hk = moran[moran["gene"].str.startswith("HK_")]["I"]
    check(len(hk) == 0 or float(hk.max()) < 0.3, f"{sample}: housekeeping Moran's I < 0.3 (max {hk.max() if len(hk) else float('nan'):.3f})")

    z = pd.read_csv(RESULTS / "spatial" / sample / "nhood_enrichment.tsv", sep="\t", index_col=0)
    diag = np.diag(z.values.astype(float)) if len(z) else np.array([])
    check(len(diag) >= 2 and (diag > 0).all(), f"{sample}: every cluster is self-enriched in the neighborhood z-scores (diag {np.round(diag, 1).tolist()})")

    with h5py.File(RESULTS / "spatial" / f"{sample}.h5ad", "r") as f:
        check("leiden" in f["obs"], f"{sample}: final h5ad carries obs/leiden")
        check("spatial" in f["obsm"], f"{sample}: final h5ad carries obsm/spatial")
        check("spatial_connectivities" in f["obsp"], f"{sample}: final h5ad carries the spatial graph")


def _h5ad_names(group):
    """obs/var names from an h5ad group: a plain string dataset, or (pandas >= 3
    string dtype) a nullable-string-array group with `values` + `mask`."""
    node = group[group.attrs["_index"]]
    arr = node["values"][()] if isinstance(node, h5py.Group) else node[()]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr])


def check_imputation(sample, method):
    stats = pd.read_csv(RESULTS / "imputed" / f"{sample}.{method}_gene_stats.tsv", sep="\t", index_col=0)
    markers = stats[stats.index.str.startswith(("DOMA_", "DOMB_"))]
    check(len(markers) >= 10, f"{sample}/{method}: >= 10 domain markers among the imputed genes ({len(markers)})")
    check(
        markers["pct_zero_imputed"].mean() < markers["pct_zero_observed"].mean(),
        f"{sample}/{method}: zero fraction of domain markers drops ({markers['pct_zero_observed'].mean():.1f}% -> {markers['pct_zero_imputed'].mean():.1f}%)",
    )
    r = float(markers["pearson_r"].median())
    check(r > 0.3, f"{sample}/{method}: median observed-vs-imputed r over domain markers > 0.3 (got {r:.3f})")

    t = truth(sample)
    with h5py.File(RESULTS / "imputed" / f"{sample}.{method}.h5ad", "r") as f:
        check("imputed" in f["layers"], f"{sample}/{method}: h5ad carries layers/imputed")
        genes = _h5ad_names(f["var"])
        cells = _h5ad_names(f["obs"])
        imp = f["layers"]["imputed"][()]
    dom = t.reindex(cells)["domain"].to_numpy()
    is_a = dom == "A"
    a_genes = np.array([g.startswith("DOMA_") for g in genes])
    b_genes = np.array([g.startswith("DOMB_") for g in genes])
    if a_genes.any() and b_genes.any():
        sep_a = imp[is_a][:, a_genes].mean() - imp[~is_a][:, a_genes].mean()
        sep_b = imp[~is_a][:, b_genes].mean() - imp[is_a][:, b_genes].mean()
        check(sep_a > 0 and sep_b > 0, f"{sample}/{method}: imputed domain markers still separate the domains (A {sep_a:.2f}, B {sep_b:.2f})")
    check(np.isfinite(imp).all() and (imp >= 0).all(), f"{sample}/{method}: imputed values are finite and non-negative")


def check_imputation_summary(methods):
    p = RESULTS / "imputed" / "imputation_summary.tsv"
    check(p.exists(), "imputation summary exists")
    if p.exists():
        df = pd.read_csv(p, sep="\t")
        check(len(df) == len(SAMPLES) * len(methods), f"imputation summary has {len(SAMPLES) * len(methods)} rows (got {len(df)})")
        check(set(df["method"]) == set(methods), f"imputation summary covers {sorted(methods)}")


def check_report():
    p = RESULTS / "report" / "scstseq_report.html"
    check(p.exists(), "report exists")
    if p.exists():
        html = p.read_text(encoding="utf-8")
        for s in SAMPLES:
            check(f"id='{s}'" in html, f"report has a section for {s}")
        check("data:image/png;base64," in html, "report embeds its images")


def main():
    check_summary()
    for s in SAMPLES:
        check_qc_cells(s)
        majority = check_clusters(s)
        check_markers(s, majority)
        check_spatial(s)
    check_report()
    methods = imputation_methods()
    if methods and (RESULTS / "imputed").exists():
        for s in SAMPLES:
            for m in methods:
                check_imputation(s, m)
        check_imputation_summary(methods)
    elif methods:
        print("[skip] imputation outputs not present (run `snakemake -d .test impute_all` to test them)")
    if FAILURES:
        print(f"\n{len(FAILURES)} assertion(s) FAILED:")
        for m in FAILURES:
            print(f"  - {m}")
        sys.exit(1)
    print("\nall assertions passed")


if __name__ == "__main__":
    main()
