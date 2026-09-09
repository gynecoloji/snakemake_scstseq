"""Imputation helpers shared by the MAGIC / ALRA / scVI rules (scanpy-free)."""
import json
import pathlib
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "workflow" / "scripts"))

import impute_summary  # noqa: E402
import scst_impute  # noqa: E402


def small_adata(with_hvg=True):
    X = sp.csr_matrix(np.array([[1.0, 0.0, 2.0], [0.0, 0.0, 1.0], [3.0, 1.0, 0.0], [2.0, 0.0, 1.0]], dtype=np.float32))
    a = ad.AnnData(X=X, obs=pd.DataFrame(index=list("abcd")), var=pd.DataFrame(index=["G1", "G2", "G3"]))
    if with_hvg:
        a.var["highly_variable"] = [True, False, True]
    a.obsm["spatial"] = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
    a.uns["scst"] = {"platform": "visium", "sample_id": "s"}
    return a


def test_select_genes_hvg_and_fallback():
    mask, mode = scst_impute.select_genes(small_adata(), "hvg")
    assert mask.tolist() == [True, False, True] and mode == "hvg"
    mask, mode = scst_impute.select_genes(small_adata(with_hvg=False), "hvg")
    assert mask.all() and mode == "all"
    mask, mode = scst_impute.select_genes(small_adata(), "all")
    assert mask.all() and mode == "all"


def test_gene_stats():
    obs = np.array([[1.0, 0.0], [0.0, 0.0], [3.0, 0.0]])
    imp = np.array([[1.0, 0.5], [0.5, 0.5], [3.0, 0.5]])
    df = scst_impute.gene_stats(sp.csr_matrix(obs), imp, ["A", "B"])
    assert df.loc["A", "pct_zero_observed"] == pytest.approx(100 / 3)
    assert df.loc["A", "pct_zero_imputed"] == 0
    assert df.loc["A", "pearson_r"] == pytest.approx(np.corrcoef(obs[:, 0], imp[:, 0])[0, 1])
    assert np.isnan(df.loc["B", "pearson_r"])  # constant columns
    assert df.loc["B", "pct_zero_observed"] == 100 and df.loc["B", "pct_zero_imputed"] == 0


def test_top_genes_respects_candidates(tmp_path):
    p = tmp_path / "moran.tsv"
    pd.DataFrame({"gene": ["Z", "G3", "G1", "G2"], "rank": [1, 2, 3, 4]}).to_csv(p, sep="\t", index=False)
    assert scst_impute.top_genes(p, ["G1", "G3"], 5) == ["G3", "G1"]
    assert scst_impute.top_genes(p, ["G1", "G3"], 1) == ["G3"]
    assert scst_impute.top_genes(p, ["G1"], 0) == []


def test_finish_writes_every_output(tmp_path):
    pytest.importorskip("matplotlib")
    a = small_adata()
    mask = np.array([True, False, True])
    imputed = np.array([[1.0, 2.0], [0.4, 1.0], [3.0, 0.3], [2.0, 1.0]], dtype=np.float32)
    moran = tmp_path / "moran.tsv"
    pd.DataFrame({"gene": ["G3", "G1", "G2"], "rank": [1, 2, 3]}).to_csv(moran, sep="\t", index=False)
    outputs = {k: str(tmp_path / f"{k}.{ext}") for k, ext in (("h5ad", "h5ad"), ("stats", "tsv"), ("run", "json"), ("png", "png"))}
    run = scst_impute.finish(a, mask, imputed, "alra", {"k": 3}, "hvg", {"rank_used": 3}, 0.0, outputs, str(moran), 2, "visium", "s")
    out = ad.read_h5ad(outputs["h5ad"])
    assert list(out.var_names) == ["G1", "G3"]
    assert np.allclose(out.layers["imputed"], imputed)
    assert out.uns["scst"]["imputation"]["method"] == "alra"
    stats = pd.read_csv(outputs["stats"], sep="\t", index_col=0)
    assert list(stats.index) == ["G1", "G3"]
    meta = json.loads(pathlib.Path(outputs["run"]).read_text())
    assert meta["n_genes_imputed"] == 2 and meta["rank_used"] == 3 and meta["params"] == {"k": 3}
    assert run["pct_zero_imputed_mean"] == 0
    assert pathlib.Path(outputs["png"]).stat().st_size > 0


def test_finish_rejects_shape_mismatch(tmp_path):
    a = small_adata()
    with pytest.raises(ValueError):
        scst_impute.finish(a, np.array([True, True, True]), np.zeros((4, 2)), "magic", {}, "all", {}, 0.0, {}, "", 0, "visium", "s")


def test_impute_summary(tmp_path):
    files = []
    for i, m in enumerate(["magic", "alra"]):
        p = tmp_path / f"{m}.json"
        p.write_text(json.dumps({"sample_id": "s", "platform": "xenium", "method": m, "n_cells": 10, "runtime_s": i, "params": {"x": 1}, "rank_used": 3 if m == "alra" else None}))
        files.append(str(p))
    df = impute_summary.summary_table(files)
    assert list(df["method"]) == ["magic", "alra"]
    assert df.columns[0] == "sample_id" and df.columns[-1] == "params"
    assert json.loads(df["params"].iloc[0]) == {"x": 1}
