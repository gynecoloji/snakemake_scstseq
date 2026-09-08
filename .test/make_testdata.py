#!/usr/bin/env python3
"""Generate the miniature Xenium + Visium sections for the executable test case.

Everything here is SYNTHETIC and DETERMINISTIC (fixed seed), so the expected
outputs are known by construction rather than by having run the pipeline once
and recorded whatever came out. That is what lets `.test/assert_outputs.py`
assert numbers instead of just file existence.

Design
------
Genes       60 genes: 20 domain-A markers (DOMA_*), 20 domain-B markers
            (DOMB_*), 15 housekeeping genes (HK_*), 5 mitochondrial genes
            (MT-*). The Xenium matrix additionally carries 7 control features
            (negative-control probes/codewords, an unassigned codeword) that
            the loader must count and remove.

Layout      xenium_1: 30 x 30 cells on a jittered grid in a 1000 x 1000 µm
            field; the left half is domain A, the right half domain B.
            visium_1 / visium_2: hexagonal spot grids; the top half of the
            array is domain A, the bottom half domain B. Each Visium section
            also lists a ring of out-of-tissue positions (absent from the
            matrix) plus 4 barcodes that ARE in the matrix but flagged
            in_tissue = 0, which the loader must drop.

Expression  Poisson counts. Domain-A cells: A markers λ=6, B markers λ=0.2,
            housekeeping λ=3, mito λ=0.6 (Visium) / 0.05 (Xenium). Domain-B
            cells mirror this.

Truth       Planted and later asserted:
              1. LOW-QUALITY cells/spots (all rates x 0.05, ~9 counts) —
                 exactly `n_low` per section, all dropped as `low_counts`.
              2. HIGH-MITO spots (Visium only; mito λ=40, > 50 % mito) —
                 exactly `n_mito` per section, all dropped as `high_mito`.
              3. Two spatial DOMAINS with disjoint marker programs — Leiden
                 must recover them, Moran's I must rank the domain markers
                 first, and neighborhood enrichment must be self-enriched.

Everything is written where the workflow expects it:
    .test/data/<sample>/...   and the per-cell truth in .test/truth/<sample>_truth.tsv

Usage:  python .test/make_testdata.py [--outdir .test]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import random
from pathlib import Path

import h5py
import numpy as np

SEED = 20260908

N_A, N_B, N_HK = 20, 20, 15
MT_GENES = ["MT-CO1", "MT-CO2", "MT-ND1", "MT-ND2", "MT-CYB"]
GENES = (
    [f"DOMA_{i:02d}" for i in range(1, N_A + 1)]
    + [f"DOMB_{i:02d}" for i in range(1, N_B + 1)]
    + [f"HK_{i:02d}" for i in range(1, N_HK + 1)]
    + MT_GENES
)
N_GENES = len(GENES)  # 60
A_IDX = np.arange(0, N_A)
B_IDX = np.arange(N_A, N_A + N_B)
HK_IDX = np.arange(N_A + N_B, N_A + N_B + N_HK)
MT_IDX = np.arange(N_A + N_B + N_HK, N_GENES)

XENIUM_CONTROLS = [
    ("NegControlProbe_00001", "Negative Control Probe"),
    ("NegControlProbe_00002", "Negative Control Probe"),
    ("NegControlProbe_00003", "Negative Control Probe"),
    ("NegControlProbe_00004", "Negative Control Probe"),
    ("NegControlCodeword_0500", "Negative Control Codeword"),
    ("NegControlCodeword_0501", "Negative Control Codeword"),
    ("UnassignedCodeword_0002", "Unassigned Codeword"),
]

LAMBDA = {"marker_hi": 6.0, "marker_lo": 0.2, "hk": 3.0, "mt_visium": 0.6, "mt_xenium": 0.05, "mt_high": 40.0, "control": 0.02}
LOW_QUALITY_SCALE = 0.05

# Sample specs. `n_low` low-quality cells/spots and (Visium) `n_mito` high-
# mito spots are planted; everything else is a normal domain cell.
SAMPLES = {
    "xenium_1": {"platform": "xenium", "n_side": 30, "n_low": 30, "n_mito": 0, "condition": "tumor"},
    "visium_1": {"platform": "visium", "rows": 24, "cols": 24, "n_low": 20, "n_mito": 15, "condition": "tumor"},
    "visium_2": {"platform": "visium", "rows": 20, "cols": 24, "n_low": 12, "n_mito": 10, "condition": "normal"},
}
N_OUT_OF_TISSUE_IN_MATRIX = 4  # Visium barcodes present in the matrix but flagged in_tissue = 0


def n_cells(sample):
    """Number of analysis cells/spots the loader hands to QC."""
    s = SAMPLES[sample]
    return s["n_side"] ** 2 if s["platform"] == "xenium" else s["rows"] * s["cols"]


def expected_kept(sample):
    s = SAMPLES[sample]
    return n_cells(sample) - s["n_low"] - s["n_mito"]


# ── expression model ────────────────────────────────────────────────────
def rates(domain, platform, quality):
    lam = np.empty(N_GENES)
    lam[A_IDX] = LAMBDA["marker_hi"] if domain == "A" else LAMBDA["marker_lo"]
    lam[B_IDX] = LAMBDA["marker_hi"] if domain == "B" else LAMBDA["marker_lo"]
    lam[HK_IDX] = LAMBDA["hk"]
    lam[MT_IDX] = LAMBDA["mt_visium"] if platform == "visium" else LAMBDA["mt_xenium"]
    if quality == "low":
        lam *= LOW_QUALITY_SCALE
    elif quality == "mito":
        lam[MT_IDX] = LAMBDA["mt_high"]
    return lam


def simulate(rng, domains, platform, qualities):
    X = np.zeros((len(domains), N_GENES), dtype=np.int32)
    for i, (d, q) in enumerate(zip(domains, qualities)):
        X[i] = rng.poisson(rates(d, platform, q))
    return X


# ── 10x HDF5 writer (no scipy: the CSC arrays are built by hand) ────────
def write_10x_h5(path, X, barcodes, names, ids, feature_types, genome="GRCh38", library_id=None):
    """Write a Cell Ranger v3-style feature-barcode HDF5: a CSC matrix of
    shape (features, barcodes) under /matrix, one column per cell."""
    X = np.asarray(X)
    n_cells_, n_feat = X.shape
    data, indices, indptr = [], [], [0]
    for c in range(n_cells_):
        nz = np.flatnonzero(X[c])
        data.append(X[c, nz])
        indices.append(nz)
        indptr.append(indptr[-1] + len(nz))
    data = np.concatenate(data).astype(np.int32) if data else np.zeros(0, np.int32)
    indices = np.concatenate(indices).astype(np.int64) if indices else np.zeros(0, np.int64)
    with h5py.File(path, "w") as f:
        f.attrs["filetype"] = "matrix"
        f.attrs["version"] = 2
        if library_id is not None:
            f.attrs["library_ids"] = np.array([library_id.encode()])
        g = f.create_group("matrix")
        g.create_dataset("barcodes", data=np.array(barcodes, dtype="S"))
        g.create_dataset("data", data=data)
        g.create_dataset("indices", data=indices)
        g.create_dataset("indptr", data=np.array(indptr, dtype=np.int64))
        g.create_dataset("shape", data=np.array([n_feat, n_cells_], dtype=np.int32))
        feat = g.create_group("features")
        feat.create_dataset("id", data=np.array(ids, dtype="S"))
        feat.create_dataset("name", data=np.array(names, dtype="S"))
        feat.create_dataset("feature_type", data=np.array(feature_types, dtype="S"))
        feat.create_dataset("genome", data=np.array([genome] * n_feat, dtype="S"))
        feat.create_dataset("_all_tag_keys", data=np.array(["genome"], dtype="S"))


def write_truth(path, barcodes, domains, qualities):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["barcode", "domain", "quality"])
        for b, d, q in zip(barcodes, domains, qualities):
            w.writerow([b, d, q])


# ── Xenium section ──────────────────────────────────────────────────────
def make_xenium(outdir, truth_dir, sample, spec, rng):
    n_side = spec["n_side"]
    n = n_side * n_side
    step = 1000.0 / n_side
    gx, gy = np.meshgrid(np.arange(n_side), np.arange(n_side), indexing="ij")
    x = (gx.ravel() + 0.5) * step + rng.uniform(-0.25 * step, 0.25 * step, n)
    y = (gy.ravel() + 0.5) * step + rng.uniform(-0.25 * step, 0.25 * step, n)
    domains = np.where(x < 500.0, "A", "B")
    qualities = np.array(["normal"] * n, dtype=object)
    qualities[rng.choice(n, spec["n_low"], replace=False)] = "low"

    X = simulate(rng, domains, "xenium", qualities)
    controls = rng.poisson(LAMBDA["control"], size=(n, len(XENIUM_CONTROLS))).astype(np.int32)
    barcodes = [f"cell_{i + 1:05d}" for i in range(n)]
    cell_area = np.clip(rng.normal(80.0, 15.0, n), 20.0, None)

    outdir.mkdir(parents=True, exist_ok=True)
    full = np.hstack([X, controls])
    names = GENES + [c[0] for c in XENIUM_CONTROLS]
    ids = [f"ENSG{i:011d}" for i in range(N_GENES)] + [c[0] for c in XENIUM_CONTROLS]
    types = ["Gene Expression"] * N_GENES + [c[1] for c in XENIUM_CONTROLS]
    write_10x_h5(outdir / "cell_feature_matrix.h5", full, barcodes, names, ids, types)
    with gzip.open(outdir / "cells.csv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["cell_id", "x_centroid", "y_centroid", "transcript_counts", "cell_area", "nucleus_area", "nucleus_count"])
        for i in range(n):
            w.writerow([barcodes[i], f"{x[i]:.3f}", f"{y[i]:.3f}", int(full[i].sum()), f"{cell_area[i]:.3f}", f"{0.3 * cell_area[i]:.3f}", 1])
    write_truth(truth_dir / f"{sample}_truth.tsv", barcodes, domains, qualities)
    return n


# ── Visium section ──────────────────────────────────────────────────────
def random_barcode(rng):
    return "".join(rng.choice(list("ACGT"), 16)) + "-1"


def make_visium(outdir, truth_dir, sample, spec, rng):
    rows, cols = spec["rows"], spec["cols"]
    # In-tissue spots on a hexagonal grid: array_col parity follows the row.
    in_r, in_c = np.meshgrid(np.arange(rows), np.arange(cols), indexing="ij")
    in_r, in_c = in_r.ravel(), in_c.ravel()
    n = rows * cols
    domains = np.where(in_r < rows / 2, "A", "B")
    qualities = np.array(["normal"] * n, dtype=object)
    picks = rng.choice(n, spec["n_low"] + spec["n_mito"], replace=False)
    qualities[picks[: spec["n_low"]]] = "low"
    qualities[picks[spec["n_low"]:]] = "mito"
    X = simulate(rng, domains, "visium", qualities)

    # 4 extra barcodes in the matrix flagged out of tissue (loader drops them),
    # placed on the row below the tissue.
    n_extra = N_OUT_OF_TISSUE_IN_MATRIX
    ex_r = np.full(n_extra, rows)
    ex_c = np.arange(n_extra)
    X_extra = simulate(rng, ["A"] * n_extra, "visium", ["normal"] * n_extra)

    seen = set()
    barcodes = []
    while len(barcodes) < n + n_extra:
        b = random_barcode(rng)
        if b not in seen:
            seen.add(b)
            barcodes.append(b)

    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "spatial").mkdir(exist_ok=True)
    ids = [f"ENSG{i:011d}" for i in range(N_GENES)]
    write_10x_h5(
        outdir / "filtered_feature_bc_matrix.h5", np.vstack([X, X_extra]), barcodes, GENES, ids,
        ["Gene Expression"] * N_GENES, library_id=sample,
    )

    def pixel(r, c):
        c2 = 2 * c + (r % 2)  # Visium array_col parity
        return 500 + r * 87, 500 + c2 * 50  # (pxl_row, pxl_col): ~100 px hex spacing

    with open(outdir / "spatial" / "tissue_positions.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"])
        for i in range(n):
            pr, pc = pixel(int(in_r[i]), int(in_c[i]))
            w.writerow([barcodes[i], 1, int(in_r[i]), 2 * int(in_c[i]) + int(in_r[i]) % 2, pr, pc])
        for j in range(n_extra):
            pr, pc = pixel(int(ex_r[j]), int(ex_c[j]))
            w.writerow([barcodes[n + j], 0, int(ex_r[j]), 2 * int(ex_c[j]) + int(ex_r[j]) % 2, pr, pc])
        # a ring of out-of-tissue positions that are NOT in the matrix
        for c in range(cols):
            pr, pc = pixel(rows + 1, c)
            w.writerow([f"OUT{c:013d}-1", 0, rows + 1, 2 * c + (rows + 1) % 2, pr, pc])
    with open(outdir / "spatial" / "scalefactors_json.json", "w") as fh:
        json.dump(
            {"spot_diameter_fullres": 45.0, "tissue_hires_scalef": 0.1, "fiducial_diameter_fullres": 72.0, "tissue_lowres_scalef": 0.03},
            fh,
        )
    write_truth(truth_dir / f"{sample}_truth.tsv", barcodes[:n], domains, qualities)
    return n


def main(outdir):
    outdir = Path(outdir)
    rng = np.random.default_rng(SEED)
    random.seed(SEED)
    truth_dir = outdir / "truth"
    truth_dir.mkdir(parents=True, exist_ok=True)
    for sample, spec in SAMPLES.items():
        target = outdir / "data" / sample
        if spec["platform"] == "xenium":
            n = make_xenium(target, truth_dir, sample, spec, rng)
        else:
            n = make_visium(target, truth_dir, sample, spec, rng)
        print(f"{sample}: {spec['platform']}, {n} cells/spots, {spec['n_low']} low-quality, {spec['n_mito']} high-mito → {target}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=str(Path(__file__).resolve().parent))
    main(ap.parse_args().outdir)
