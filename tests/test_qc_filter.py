"""QC threshold logic and metric computation (no scanpy needed)."""
import pathlib
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import qc_filter  # noqa: E402
import scst_common as sc_  # noqa: E402


def xenium_obs():
    return pd.DataFrame(
        {
            "total_counts": [100, 5, 100, 100, 100, 100, 5000],
            "n_genes_by_counts": [40, 3, 2, 40, 40, 40, 40],
            "cell_area": [80, 80, 80, 5, 80, 900, 80],
            "control_frac": [0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0],
        },
        index=[f"c{i}" for i in range(7)],
    )


def test_apply_thresholds_xenium_reasons_and_precedence():
    thr = {"min_counts": 10, "max_counts": 1000, "min_genes": 5, "min_area": 10, "max_area": 500, "max_control_frac": 0.1}
    v = sc_.apply_thresholds(xenium_obs(), "xenium", thr)
    assert v["keep"].tolist() == [True, False, False, False, False, False, False]
    assert v["drop_reason"].tolist() == ["", "low_counts", "low_genes", "small_area", "high_control_frac", "large_area", "high_counts"]


def test_zero_disables_a_rule():
    thr = {"min_counts": 10, "max_counts": 0, "min_genes": 0, "min_area": 0, "max_area": 0, "max_control_frac": 0}
    v = sc_.apply_thresholds(xenium_obs(), "xenium", thr)
    # only the low-counts cell is dropped; high counts / few genes / areas / controls pass
    assert v["drop_reason"].tolist() == ["", "low_counts", "", "", "", "", ""]


def test_first_failing_rule_wins():
    obs = pd.DataFrame({"total_counts": [1], "n_genes_by_counts": [1], "cell_area": [1], "control_frac": [1.0]}, index=["c"])
    thr = {"min_counts": 10, "min_genes": 5, "min_area": 10, "max_control_frac": 0.1}
    v = sc_.apply_thresholds(obs, "xenium", thr)
    assert v["drop_reason"].iloc[0] == "low_counts"


def test_missing_metric_column_is_skipped():
    obs = pd.DataFrame({"total_counts": [100], "n_genes_by_counts": [50]}, index=["s"])
    v = sc_.apply_thresholds(obs, "visium", {"min_counts": 10, "min_genes": 5, "max_pct_mt": 20})
    assert v["keep"].iloc[0]


def test_visium_mito_and_hb():
    obs = pd.DataFrame(
        {"total_counts": [1000, 1000, 1000], "n_genes_by_counts": [500, 500, 500], "pct_counts_mt": [5, 40, 5], "pct_counts_hb": [0, 0, 30]},
        index=["a", "b", "c"],
    )
    v = sc_.apply_thresholds(obs, "visium", {"min_counts": 500, "min_genes": 250, "max_pct_mt": 30, "max_pct_hb": 10})
    assert v["drop_reason"].tolist() == ["", "high_mito", "high_hb"]


def test_unknown_platform_raises():
    with pytest.raises(ValueError):
        sc_.apply_thresholds(pd.DataFrame(), "cosmx", {})


def test_nan_metric_never_fails_a_rule():
    obs = pd.DataFrame({"total_counts": [np.nan], "n_genes_by_counts": [50]}, index=["s"])
    v = sc_.apply_thresholds(obs, "visium", {"min_counts": 10})
    assert v["keep"].iloc[0]


def small_adata():
    X = sp.csr_matrix(np.array([[5, 0, 5, 0], [0, 0, 0, 0], [2, 2, 2, 2]], dtype=np.float32))
    a = ad.AnnData(X=X, obs=pd.DataFrame(index=["a", "b", "c"]), var=pd.DataFrame(index=["MT-CO1", "RPL3", "HBB", "GENE"]))
    a.uns["scst"] = {"sample_id": "t", "condition": ""}
    return a


def test_compute_metrics():
    a = qc_filter.compute_metrics(small_adata(), {"mito_pattern": "^MT-", "ribo_pattern": "^RP[SL]", "hb_pattern": "^HB[^P]"})
    assert a.obs["total_counts"].tolist() == [10, 0, 8]
    assert a.obs["n_genes_by_counts"].tolist() == [2, 0, 4]
    assert a.obs["pct_counts_mt"].tolist() == pytest.approx([50, 0, 25])
    assert a.obs["pct_counts_ribo"].tolist() == pytest.approx([0, 0, 25])
    assert a.obs["pct_counts_hb"].tolist() == pytest.approx([50, 0, 25])
    assert a.var["mt"].tolist() == [True, False, False, False]


def test_filter_genes():
    a = small_adata()
    removed = qc_filter.filter_genes(a, 2)
    assert removed == 2  # MT-CO1 and HBB are detected in cells a + c; RPL3 and GENE only in c
    assert list(a.var_names) == ["MT-CO1", "HBB"]
    assert qc_filter.filter_genes(small_adata(), 0) == 0


def test_filter_genes_removing_everything_raises():
    with pytest.raises(ValueError):
        qc_filter.filter_genes(small_adata(), 10)


def test_run_end_to_end(tmp_path):
    a = small_adata()
    a.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]])
    a.obs["sample_id"] = "t"
    a.uns["scst"]["platform"] = "visium"
    a.write_h5ad(tmp_path / "raw.h5ad")
    thr = {"min_counts": 5, "max_counts": 0, "min_genes": 1, "max_pct_mt": 30, "max_pct_hb": 0, "min_cells": 1}
    stats = qc_filter.run(
        str(tmp_path / "raw.h5ad"), str(tmp_path / "f.h5ad"), str(tmp_path / "cells.tsv"), str(tmp_path / "stats.json"),
        "visium", thr, {"mito_pattern": "^MT-", "ribo_pattern": "^RP[SL]", "hb_pattern": "^HB[^P]"},
    )
    # a: 50 % mito -> high_mito; b: 0 counts -> low_counts; c kept
    assert stats["n_cells_kept"] == 1 and stats["drop_high_mito"] == 1 and stats["drop_low_counts"] == 1
    cells = pd.read_csv(tmp_path / "cells.tsv", sep="\t", index_col=0)
    assert cells["qc_drop_reason"].fillna("").tolist() == ["high_mito", "low_counts", ""]
    assert list(cells.columns[:2]) == ["x", "y"]
    kept = ad.read_h5ad(tmp_path / "f.h5ad")
    assert kept.n_obs == 1 and kept.n_vars == 4  # c expresses every gene; min_cells=1


def test_run_dropping_everything_exits(tmp_path):
    a = small_adata()
    a.obsm["spatial"] = np.zeros((3, 2))
    a.uns["scst"]["platform"] = "visium"
    a.write_h5ad(tmp_path / "raw.h5ad")
    with pytest.raises(SystemExit):
        qc_filter.run(
            str(tmp_path / "raw.h5ad"), str(tmp_path / "f.h5ad"), str(tmp_path / "c.tsv"), str(tmp_path / "s.json"),
            "visium", {"min_counts": 1000, "min_genes": 0, "max_pct_mt": 0, "min_cells": 0},
            {"mito_pattern": "^MT-", "ribo_pattern": "^RP[SL]", "hb_pattern": "^HB[^P]"},
        )
