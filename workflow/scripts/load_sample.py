"""Rule `load_sample`: read one vendor output directory into an AnnData with
spatial coordinates.

xenium  cell_feature_matrix.h5 + cells.parquet (or cells.csv.gz): centroids
        (µm) → obsm["spatial"]; cell/nucleus area and transcript counts →
        obs; control probes/codewords are counted per cell (control_counts,
        control_frac) and then REMOVED from the gene matrix.
visium  filtered_feature_bc_matrix.h5 + spatial/: full-resolution pixel
        positions → obsm["spatial"] (x = column, y = row); array row/col and
        in_tissue → obs; scale factors + images → uns["spatial"]; spots
        outside the tissue are dropped.

Barcodes are prefixed `{sample_id}_` so later merges never collide; the
original barcode is kept in obs["barcode"].
"""


import os

import numpy as np

from scst_common import (
    COORD_UNIT,
    GENE_FEATURE_TYPE,
    PLATFORMS,
    VISIUM,
    XENIUM,
    log,
    read_10x_h5,
    read_cells_table,
    read_images,
    read_positions,
    read_scalefactors,
    redirect_log,
)

XENIUM_OBS_COLS = (
    "x_centroid",
    "y_centroid",
    "transcript_counts",
    "cell_area",
    "nucleus_area",
    "nucleus_count",
    "segmentation_method",
)


def _row_sums(X):
    return np.asarray(X.sum(axis=1)).ravel()


def read_xenium(path):
    adata = read_10x_h5(os.path.join(path, "cell_feature_matrix.h5"))
    cells = read_cells_table(path)
    missing = adata.obs_names.difference(cells.index)
    if len(missing):
        raise ValueError(
            f"{path}: {len(missing)} cell(s) in cell_feature_matrix.h5 are absent from the "
            f"cells table (e.g. {list(missing[:3])})"
        )
    cells = cells.loc[adata.obs_names]
    for col in XENIUM_OBS_COLS:
        if col in cells.columns:
            adata.obs[col] = cells[col].values
    for col in ("x_centroid", "y_centroid"):
        if col not in cells.columns:
            raise ValueError(f"{path}: cells table has no {col} column")
    adata.obsm["spatial"] = cells[["x_centroid", "y_centroid"]].to_numpy(dtype=float)

    is_gene = (adata.var["feature_types"] == GENE_FEATURE_TYPE).to_numpy()
    control_types = sorted(set(adata.var["feature_types"][~is_gene]))
    control = _row_sums(adata.X[:, ~is_gene]) if (~is_gene).any() else np.zeros(adata.n_obs)
    gene = _row_sums(adata.X[:, is_gene])
    total = control + gene
    adata.obs["control_counts"] = control
    adata.obs["control_frac"] = np.where(total > 0, control / np.maximum(total, 1), 0.0)
    log(
        f"xenium: {adata.n_obs} cells, {int(is_gene.sum())} genes, "
        f"{int((~is_gene).sum())} control features {control_types} removed"
    )
    adata = adata[:, is_gene].copy()
    adata.uns["scst"] = {"control_feature_types": control_types}
    return adata


def read_visium(path, library_id):
    adata = read_10x_h5(os.path.join(path, "filtered_feature_bc_matrix.h5"))
    is_gene = (adata.var["feature_types"] == GENE_FEATURE_TYPE).to_numpy()
    if not is_gene.all():
        log(f"visium: dropping {int((~is_gene).sum())} non-gene features")
        adata = adata[:, is_gene].copy()
    spatial_dir = os.path.join(path, "spatial")
    pos = read_positions(spatial_dir)
    missing = adata.obs_names.difference(pos.index)
    if len(missing):
        raise ValueError(
            f"{path}: {len(missing)} barcode(s) in the matrix are absent from the tissue "
            f"positions (e.g. {list(missing[:3])})"
        )
    pos = pos.loc[adata.obs_names]
    adata.obs["in_tissue"] = pos["in_tissue"].astype(int).values
    adata.obs["array_row"] = pos["array_row"].astype(int).values
    adata.obs["array_col"] = pos["array_col"].astype(int).values
    # (x, y) = (column, row) in full-resolution image pixels.
    adata.obsm["spatial"] = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)
    adata.uns["spatial"] = {
        library_id: {
            "images": read_images(spatial_dir),
            "scalefactors": read_scalefactors(spatial_dir),
            "metadata": {"source": os.path.abspath(path)},
        }
    }
    n_out = int((adata.obs["in_tissue"] == 0).sum())
    if n_out:
        log(f"visium: dropping {n_out} spots outside the tissue")
        adata = adata[adata.obs["in_tissue"] == 1].copy()
    log(f"visium: {adata.n_obs} in-tissue spots, {adata.n_vars} genes")
    adata.uns["scst"] = {"n_out_of_tissue": n_out}
    return adata


def load_sample(platform, path, sample_id, condition=""):
    if platform == XENIUM:
        adata = read_xenium(path)
    elif platform == VISIUM:
        adata = read_visium(path, sample_id)
    else:
        raise ValueError(f"unknown platform {platform!r}; expected one of {PLATFORMS}")
    if adata.n_obs == 0:
        raise ValueError(f"{sample_id}: no cells/spots after loading {path}")
    adata.obs["barcode"] = adata.obs_names.astype(str)
    adata.obs_names = [f"{sample_id}_{b}" for b in adata.obs_names]
    adata.obs["sample_id"] = sample_id
    adata.obs["platform"] = platform
    if condition:
        adata.obs["condition"] = condition
    adata.uns["scst"] = {
        **adata.uns.get("scst", {}),
        "platform": platform,
        "sample_id": sample_id,
        "condition": condition or "",
        "coord_unit": COORD_UNIT[platform],
        "source_dir": os.path.abspath(path),
    }
    return adata


if "snakemake" in globals():  # pragma: no cover
    redirect_log(snakemake.log[0])  # noqa: F821
    _adata = load_sample(
        snakemake.params.platform,  # noqa: F821
        snakemake.params.path,  # noqa: F821
        snakemake.wildcards.sample,  # noqa: F821
        snakemake.params.condition,  # noqa: F821
    )
    _adata.write_h5ad(snakemake.output.h5ad)  # noqa: F821
    log(f"wrote {_adata.n_obs} x {_adata.n_vars} to {snakemake.output.h5ad}")  # noqa: F821
