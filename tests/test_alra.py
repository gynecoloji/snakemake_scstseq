"""The ALRA re-implementation on a planted low-rank matrix with dropouts."""
import pathlib
import sys

import numpy as np
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "workflow" / "scripts"))

import alra  # noqa: E402


def planted(n_cells=300, n_genes=80, rank=3, dropout=0.4, seed=0):
    """Two cell programs (rank-3 signal incl. a baseline) + Poisson-ish noise,
    then random dropouts. Returns (observed, truth, dropout mask)."""
    rng = np.random.default_rng(seed)
    programs = rng.uniform(0.5, 3.0, size=(rank, n_genes))
    loadings = np.abs(rng.normal(size=(n_cells, rank)))
    truth = loadings @ programs
    noisy = np.clip(truth + rng.normal(scale=0.15, size=truth.shape), 0, None)
    drop = rng.random(truth.shape) < dropout
    observed = np.where(drop, 0.0, noisy)
    return observed.astype(np.float32), truth, drop


def test_randomized_svd_matches_numpy_on_small_matrix():
    rng = np.random.default_rng(1)
    A = rng.normal(size=(40, 25))
    U, s, Vt = alra.randomized_svd(A, k=5, n_iter=6, seed=1)
    s_ref = np.linalg.svd(A, compute_uv=False)[:5]
    assert np.allclose(s, s_ref, rtol=1e-3)
    assert np.allclose(U.T @ U, np.eye(5), atol=1e-6)


def test_randomized_svd_accepts_sparse():
    obs, _, _ = planted()
    U, s, Vt = alra.randomized_svd(sp.csr_matrix(obs), k=3, seed=0)
    assert U.shape == (300, 3) and s.shape == (3,) and Vt.shape == (3, 80)


def test_choose_k_finds_the_planted_rank():
    obs, _, _ = planted(rank=3)
    k, s = alra.choose_k(obs, k_max=40, seed=0)
    assert 2 <= k <= 6


def test_choose_k_on_tiny_matrix_does_not_crash():
    k, _ = alra.choose_k(np.ones((5, 4)), k_max=100)
    assert k >= 1


def test_alra_recovers_dropouts_and_keeps_observed_nonzeros():
    obs, truth, drop = planted()
    imputed, k = alra.alra(obs, k=0, k_max=40, seed=0)
    assert imputed.shape == obs.shape and imputed.dtype == np.float32
    assert (imputed >= 0).all()
    # observed non-zeros are never zeroed out
    assert (imputed[obs > 0] > 0).all()
    # dropouts with a real underlying value are mostly filled in
    real_dropouts = drop & (truth > 1.0)
    assert (imputed[real_dropouts] > 0).mean() > 0.8
    # and the filled values track the truth
    r = np.corrcoef(imputed[real_dropouts], truth[real_dropouts])[0, 1]
    assert r > 0.6
    # the zero fraction went down
    assert (imputed == 0).mean() < (obs == 0).mean()


def test_alra_explicit_rank_and_sparse_input_agree():
    obs, _, _ = planted()
    a, k_a = alra.alra(obs, k=3, seed=0)
    b, k_b = alra.alra(sp.csr_matrix(obs), k=3, seed=0)
    assert k_a == k_b == 3
    assert np.allclose(a, b, atol=1e-4)


def test_alra_rank_is_capped_to_the_matrix():
    obs, _, _ = planted(n_cells=20, n_genes=10)
    _, k = alra.alra(obs, k=50, seed=0)
    assert k <= 9


@pytest.mark.parametrize("quantile", [0.001, 0.01])
def test_alra_quantile_runs(quantile):
    obs, _, _ = planted()
    imputed, _ = alra.alra(obs, k=3, quantile=quantile, seed=0)
    assert np.isfinite(imputed).all()
