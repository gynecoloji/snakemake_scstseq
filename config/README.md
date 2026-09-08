# Configuration

This workflow is configured through two files in this directory:

- `config.yaml` — all workflow parameters (nested by section)
- `samples.csv` — the sample sheet

plus the per-sample vendor output directories you place under `data/`.

## Sample sheet (`config/samples.csv`)

CSV with one row per tissue section:

| column      | required | description                                                                 |
|-------------|----------|-----------------------------------------------------------------------------|
| `sample_id` | yes      | Section name; used in every output path and as the barcode prefix. Must be unique. |
| `platform`  | yes      | `xenium` (single cells, imaging-based) or `visium` (55 µm spots, sequencing-based). |
| `path`      | yes      | The vendor output directory (see below). Relative to `paths.data_dir` unless absolute. |
| `condition` | no       | Free-text group label, carried into the QC summary and report (reserved for future comparative stages). |

Example:

```csv
sample_id,platform,path,condition
xenium_tumor_1,xenium,xenium_tumor_1/outs,tumor
visium_tumor_1,visium,visium_tumor_1/outs,tumor
visium_normal_1,visium,visium_normal_1/outs,normal
```

### What each `path` must contain

| platform | files read from `path` | produced by |
|---|---|---|
| `xenium` | `cell_feature_matrix.h5` and `cells.parquet` (or `cells.csv.gz`) | Xenium Onboard Analysis `outs/` |
| `visium` | `filtered_feature_bc_matrix.h5`, `spatial/tissue_positions.csv` (or the older `tissue_positions_list.csv`), `spatial/scalefactors_json.json`; `spatial/tissue_hires_image.png` / `tissue_lowres_image.png` are used when present | Space Ranger `outs/` |

Everything else in those directories (transcripts, boundaries, morphology
images, zarr stores) is ignored — nothing is copied, so you can point `path`
straight at the vendor directory. Missing required files fail at DAG-build
time with the missing path named.

Xenium **control features** (negative-control probes/codewords, unassigned
and genomic controls) are counted per cell into `control_counts` /
`control_frac` and then removed from the gene matrix. Visium spots outside
the tissue (`in_tissue == 0`) are dropped at load time.

## Parameters (`config/config.yaml`)

Every parameter — type, default, description — is defined in the config schema,
[`workflow/schemas/config.schema.yaml`](../workflow/schemas/config.schema.yaml).
The workflow validates `config.yaml` against it on every run (and fills in
defaults for anything omitted), and the Snakemake Workflow Catalog renders it as
a parameter table. `config.yaml` ships with working defaults and an inline
comment on each parameter.

Highlights:

- `genome.*_pattern` — gene-name regexes for the mito / ribo / hemoglobin
  QC subsets (defaults are human; switch to `^mt-` / `^Rp[sl]` / `^Hb[^p]`
  for mouse).
- `qc.xenium.*` / `qc.visium.*` — cell/spot thresholds per platform; every
  threshold uses `0` = disabled. Cells are dropped by the first failing rule,
  and that reason is recorded.
- `processing.xenium` / `processing.visium` — normalization target, HVG
  selection, scaling, PCs, neighbours, Leiden resolution.
- `spatial.xenium` / `spatial.visium` — how the spatial neighbours graph is
  built (kNN / Delaunay / radius for cells; hexagonal grid rings for spots);
  `spatial.nhood_enrichment.n_perms` and `spatial.moran.*` control the
  permutation tests.
- `report.enabled` — build the self-contained HTML report as part of `all`.
- `threads.*` — threads per rule.
