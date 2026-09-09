"""ALRA — adaptively-thresholded low-rank approximation for imputation.

A numpy/scipy re-implementation of Linderman et al. (2022, Nature
Communications 13:192; github.com/KlugerLab/ALRA) so the workflow does not
need R:

1. rank-k randomized SVD of the log-normalized cells x genes matrix,
   with k chosen from the singular-value noise floor when not given;
2. per-gene thresholding: reconstructed values at or below the magnitude
   of the gene's `quantile` (most negative) value are set to 0 — the
   symmetry of the noise around zero tells biological zeros from dropouts;
3. per-gene rescaling so the non-zero mean / sd match the observed
   non-zero mean / sd, negatives clipped at 0;
4. observed non-zeros that the approximation zeroed are restored.

Only numpy and scipy.sparse are used, so it is unit-testable in the static CI
job and runs inside the py-scst env without extra packages.
"""

import numpy as np
import scipy.sparse as sp


def randomized_svd(A, k, n_oversamples=10, n_iter=4, seed=0):
    """Halko-Martinsson-Tropp randomized SVD of a dense or sparse matrix.

    Returns U (m x k), s (k,), Vt (k x n) for the top-k singular triplets.
    """
    rng = np.random.default_rng(seed)
    m, n = A.shape
    l = int(min(k + n_oversamples, min(m, n)))
    omega = rng.standard_normal((n, l)).astype(np.float32 if A.dtype == np.float32 else np.float64)
    Y = A @ omega
    Y = np.asarray(Y)
    for _ in range(n_iter):
        Q, _ = np.linalg.qr(Y)
        Z = np.asarray(A.T @ Q)
        Q, _ = np.linalg.qr(Z)
        Y = np.asarray(A @ Q)
    Q, _ = np.linalg.qr(Y)
    B = np.asarray((A.T @ Q)).T  # (l x n)
    Ub, s, Vt = np.linalg.svd(B, full_matrices=False)
    U = Q @ Ub
    return U[:, :k], s[:k], Vt[:k]


def choose_k(A, k_max=100, noise_start=None, n_sds=6.0, seed=0):
    """Rank selection as in ALRA's `choose_k`: singular values are computed up
    to k_max; successive differences beyond `noise_start` define the noise
    level; k is the last index whose drop stands more than `n_sds` standard
    deviations above that noise.

    Returns (k, singular_values). Falls back to a small rank when the spectrum
    has no clear elbow (e.g. tiny synthetic panels).
    """
    m, n = A.shape
    K = int(min(k_max, min(m, n) - 1))
    if K < 3:
        return max(K, 1), np.array([])
    if noise_start is None:
        noise_start = max(int(round(0.8 * K)), 2)
    noise_start = int(min(noise_start, K - 2))
    _, s, _ = randomized_svd(A, K, seed=seed)
    diffs = s[:-1] - s[1:]
    noise = diffs[noise_start - 1 :]
    mu, sigma = float(np.mean(noise)), float(np.std(noise, ddof=1)) if len(noise) > 1 else 0.0
    if sigma <= 0:
        return max(2, min(10, K)), s
    num_sds = (diffs - mu) / sigma
    above = np.flatnonzero(num_sds > n_sds)
    k = int(above[-1] + 1) if len(above) else max(2, min(10, K))
    return max(k, 2), s


def _nonzero_stats(M):
    """Per-column mean and sd of the non-zero entries (nan when none)."""
    nz = M != 0
    cnt = nz.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(cnt > 0, M.sum(axis=0) / np.maximum(cnt, 1), np.nan)
        sq = (M * M).sum(axis=0)
        var = np.where(cnt > 1, (sq - cnt * mean**2) / np.maximum(cnt - 1, 1), np.nan)
    sd = np.sqrt(np.clip(var, 0, None))
    return mean, sd


def alra(A, k=0, quantile=0.001, k_max=100, seed=0):
    """Impute a log-normalized cells x genes matrix (dense ndarray or scipy
    sparse). Returns (imputed dense float32 array, rank used)."""
    A_sp = sp.csr_matrix(A) if not sp.issparse(A) else A.tocsr()
    A_sp = A_sp.astype(np.float32)
    if k <= 0:
        k, _ = choose_k(A_sp, k_max=k_max, seed=seed)
    k = int(min(k, min(A_sp.shape) - 1))
    U, s, Vt = randomized_svd(A_sp, k, seed=seed)
    A_k = (U * s) @ Vt  # dense reconstruction
    A_k = A_k.astype(np.float32)

    # 2. per-gene threshold at the magnitude of the most negative values
    thresh = np.abs(np.quantile(A_k, quantile, axis=0))
    A_cor = np.where(A_k <= thresh[None, :], 0.0, A_k).astype(np.float32)

    # 3. rescale non-zeros per gene to the observed non-zero mean / sd
    A_dense = A_sp.toarray()
    mu_1, sd_1 = _nonzero_stats(A_cor)
    mu_2, sd_2 = _nonzero_stats(A_dense)
    toscale = np.isfinite(sd_1) & np.isfinite(sd_2) & (sd_1 > 0)
    scale = np.ones(A_cor.shape[1], dtype=np.float32)
    shift = np.zeros(A_cor.shape[1], dtype=np.float32)
    scale[toscale] = sd_2[toscale] / sd_1[toscale]
    shift[toscale] = -mu_1[toscale] * sd_2[toscale] / sd_1[toscale] + mu_2[toscale]
    nz = A_cor != 0
    A_sc = np.where(nz, A_cor * scale[None, :] + shift[None, :], 0.0).astype(np.float32)
    A_sc[A_sc < 0] = 0.0

    # 4. restore observed non-zeros the approximation zeroed out
    restore = (A_dense > 0) & (A_sc == 0)
    A_sc[restore] = A_dense[restore]
    return A_sc, k
