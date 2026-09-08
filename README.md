# Spatial Transcriptomics Analysis Pipeline (Xenium + Visium)

[![CI](https://github.com/gynecoloji/snakemake_scstseq/actions/workflows/ci.yml/badge.svg)](https://github.com/gynecoloji/snakemake_scstseq/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/gynecoloji/snakemake_scstseq?label=release)](https://github.com/gynecoloji/snakemake_scstseq/releases/latest)
[![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A58.0-brightgreen.svg)](https://snakemake.github.io)
[![Docker Pulls](https://img.shields.io/docker/pulls/gynecoloji/scstseq-pipeline)](https://hub.docker.com/r/gynecoloji/scstseq-pipeline)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
<!-- After the first GitHub release with the Zenodo integration enabled, add:
[![DOI](https://zenodo.org/badge/<repo-id>.svg)](https://doi.org/<concept-doi>) -->

A Snakemake workflow for **spatial transcriptomics** that starts from the
vendor output directories of both **single-cell, imaging-based** sections
(10x Genomics **Xenium**) and **spot-level, sequencing-based** sections
(10x Genomics **Visium**), and produces per-section QC with platform-specific
thresholds, normalization / clustering / marker genes, a spatial neighbours
graph with neighborhood enrichment and Moran's I spatially variable genes,
and a self-contained cross-sample HTML report.

The repository mirrors the structure and automation of its siblings
[snakemake_ATACseq_spikein](https://github.com/gynecoloji/snakemake_ATACseq_spikein)
and [snakemake_scRNAseq](https://github.com/gynecoloji/snakemake_scRNAseq):
standard [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/)
layout, schema-validated configuration, per-rule conda environments,
Docker/Apptainer containers, CI with an executable test case that asserts
planted values, and release-please versioning.

## Overview

One `snakemake --use-conda` run builds every stage in dependency order
(unified DAG); each stage is also a target of its own.

| stage | target | per section |
|---|---|---|
| **QC** | `qc_all` | platform reader → per-cell/spot metrics (counts, genes, % mito/ribo/hemoglobin; Xenium: cell/nucleus area, control-feature fraction) → platform thresholds with the **first failing rule recorded as the drop reason** → QC plots → cross-sample summary |
| **Processing** | `process_all` | `normalize_total` + `log1p` → HVGs (Visium; Xenium panels use every gene) → optional scaling → PCA → neighbours → **Leiden** → UMAP → **marker genes** per cluster (dotplot + table) |
| **Spatial** | `spatial_all` | spatial neighbours graph (Xenium: kNN / Delaunay / radius over µm centroids; Visium: hexagonal grid rings) → **neighborhood enrichment** between clusters → **Moran's I** spatially variable genes (permutation p-values, BH-FDR) → top SVGs drawn on the section |
| **Report** | `all` (default) | one self-contained HTML page: summary table + every figure and top table per section, images embedded, no JavaScript |

Everything is Python: [scanpy](https://scanpy.readthedocs.io) for I/O, QC and
clustering, [squidpy](https://squidpy.readthedocs.io) for the spatial
statistics, in one conda environment (`workflow/envs/py-scst.yaml`).

**Roadmap** (not yet implemented): multi-sample integration, cell-type
annotation (marker scoring / reference mapping), Visium deconvolution, Visium
HD (its positions file is parquet), CosMx / MERSCOPE readers.

## Workflow diagram

![Rule graph](images/rulegraph.svg)

Regenerate it from the executable test case with:

```bash
snakemake -s workflow/Snakefile -d .test --forceall --rulegraph | dot -Tsvg > images/rulegraph.svg
```

## Requirements

- [Snakemake](https://snakemake.github.io) ≥ 8.0 with conda/mamba
  (`snakemake-minimal` + `pandas` is enough as the driver), **or** just
  Apptainer/Docker using the provided image (see below).
- The per-rule tool environment (`workflow/envs/py-scst.yaml`: Python 3.12,
  scanpy 1.11, squidpy 1.8, anndata, h5py, pyarrow, leidenalg, igraph,
  matplotlib, seaborn) is created automatically on the first `--use-conda`
  run.

## Input files

Point each `samples.csv` row at the **vendor output directory**; nothing is
copied or converted beforehand.

| `platform` | files read from `path` | produced by |
|---|---|---|
| `xenium` | `cell_feature_matrix.h5`, `cells.parquet` (or `cells.csv.gz`) | Xenium Onboard Analysis `outs/` |
| `visium` | `filtered_feature_bc_matrix.h5`, `spatial/tissue_positions.csv` (or `tissue_positions_list.csv`), `spatial/scalefactors_json.json`; tissue images used when present | Space Ranger `outs/` |

Xenium control features (negative-control probes/codewords, unassigned and
genomic controls) are counted per cell into `control_frac` and removed from
the gene matrix; Visium spots outside the tissue are dropped at load time.
See [`config/README.md`](config/README.md) for the column reference.

```csv
sample_id,platform,path,condition
xenium_tumor_1,xenium,xenium_tumor_1/outs,tumor
visium_tumor_1,visium,visium_tumor_1/outs,tumor
visium_normal_1,visium,visium_normal_1/outs,normal
```

## Installation

```bash
git clone https://github.com/gynecoloji/snakemake_scstseq.git
cd snakemake_scstseq

# Driver environment (the tool env is created on the first --use-conda run).
# The containers and CI pin snakemake-minimal 9.3.2; any Snakemake >= 8 works.
mamba create -n scstseq -c conda-forge -c bioconda snakemake-minimal=9.3.2 pandas
conda activate scstseq
```

> **No local install?** Use the container image instead — see
> [Container execution](#container-execution-docker--apptainer).

## Configuration

All parameters live in [`config/config.yaml`](config/config.yaml), documented
inline and validated on every run against
[`workflow/schemas/config.schema.yaml`](workflow/schemas/config.schema.yaml)
(the single source of truth: types, defaults, descriptions; the Snakemake
Workflow Catalog renders it as a parameter table). The sections are
per-platform where the platforms differ:

| section | Xenium default | Visium default |
|---|---|---|
| `qc.<platform>` | `min_counts: 10`, `min_genes: 5`, `min_area` / `max_area` / `max_control_frac` off, `min_cells: 5` | `min_counts: 500`, `min_genes: 250`, `max_pct_mt: 30`, `max_pct_hb` off, `min_cells: 10` |
| `processing.<platform>` | `target_sum: 100`, every gene, no scaling, 30 PCs, 15 neighbours, Leiden 0.5 | median target, 2000 `seurat_v3` HVGs, scaling, 30 PCs, 15 neighbours, Leiden 0.5 |
| `spatial.<platform>` | `generic` coordinates, kNN `n_neighs: 10` (or `delaunay` / `radius` in µm) | `grid`, `n_neighs: 6`, `n_rings: 1` |

Every threshold uses `0` = disabled. Shared keys: `genome.*_pattern`
(mito / ribo / hemoglobin gene-name regexes; switch for mouse),
`processing.markers`, `processing.seed`, `spatial.nhood_enrichment.n_perms`,
`spatial.moran.{n_perms,n_top,plot_top}`, `report.enabled`, `threads.*`.

## Running the pipeline

The workflow is the standard-layout `workflow/Snakefile`; Snakemake finds it
automatically from the repo root, so `-s workflow/Snakefile` is optional.

### Dry run

```bash
snakemake -s workflow/Snakefile -n            # validates config + builds the DAG
```

### Local execution

```bash
snakemake -s workflow/Snakefile --use-conda --cores 8                # everything (default target)
snakemake -s workflow/Snakefile --use-conda --cores 8 qc_all         # one stage
snakemake -s workflow/Snakefile --use-conda --cores 8 process_all
snakemake -s workflow/Snakefile --use-conda --cores 8 spatial_all
```

> **Running from a batch script?** Invoke Snakemake by its **absolute path**
> rather than putting the driver environment on `PATH`, so the per-rule
> environment's own Python is the one the scripts run with.

### Cluster execution

For a SLURM cluster, the standard executor plugin works unchanged (not
benchmarked here):

```bash
pip install snakemake-executor-plugin-slurm
snakemake -s workflow/Snakefile --use-conda --executor slurm --jobs 20 \
  --default-resources slurm_partition=<partition> mem_mb=32000 runtime=240
```

Large Xenium sections (hundreds of thousands of cells) are the memory-hungry
case: `spatial_stats` holds the kNN graph and the permutation tests, so give
that rule the most memory and threads (`threads.spatial_stats`).

### Container execution (Docker / Apptainer)

One prebuilt image covers the whole workflow — you install nothing except
Docker or Apptainer:

| image | contents | used for |
|---|---|---|
| **`gynecoloji/scstseq-pipeline`** | Snakemake + the pre-built `py-scst` conda env | every stage |

**Download** — no build required. On HPC (Apptainer / Singularity) a single
command downloads it and writes a ready-to-run `.sif` in your current
directory:

```bash
apptainer pull scstseq-pipeline.sif docker://gynecoloji/scstseq-pipeline:latest   # HPC
docker pull gynecoloji/scstseq-pipeline:latest                                     # Docker
```

The image's entrypoint is
`snakemake --use-conda --conda-frontend mamba --conda-prefix /opt/wf-conda`,
so anything after the image name goes straight to `snakemake`. Data are not
baked in; run from your project directory (which holds `config/` and `data/`).

```bash
# Apptainer (auto-mounts $HOME, /tmp and the current directory; runs as you)
module load apptainer                                   # if your cluster uses modules
apptainer run scstseq-pipeline.sif -s workflow/Snakefile --cores 16
apptainer run scstseq-pipeline.sif -s workflow/Snakefile --cores 16 qc_all
apptainer run scstseq-pipeline.sif -s workflow/Snakefile --cores 1 -n

# Docker
docker run --rm -v "$(pwd)":/workflow -e HOME=/tmp --user "$(id -u):$(id -g)" \
    gynecoloji/scstseq-pipeline:latest -s workflow/Snakefile --cores 16
./run_pipeline.sh --cores 16                            # convenience wrapper (builds locally if needed)
docker compose run --rm scstseq -s workflow/Snakefile --cores 16 qc_all
```

If your vendor directories live outside the project (e.g. on scratch), bind
them in — `--bind /scratch/xenium:/scratch/xenium` — and use absolute paths
in `samples.csv`. Build the image yourself with `docker compose build` or
`apptainer build --fakeroot scstseq-pipeline.sif apptainer.def`; details in
[`DOCKER.md`](DOCKER.md).

## Deploying with snakedeploy

This repository follows the Snakemake Workflow Catalog structure, so it can be
deployed into another project without cloning it by hand:

```bash
pip install snakedeploy
snakedeploy deploy-workflow https://github.com/gynecoloji/snakemake_scstseq . --tag main
```

Then edit the `config/` copy it writes, place your data, and run with
`snakemake --use-conda`.

## Outputs (`results/`)

```
results/
├── raw/{sample}.h5ad                       # loaded section: counts, obsm["spatial"], platform metadata
├── qc/
│   ├── {sample}.filtered.h5ad              # cells/spots + genes that passed QC
│   ├── {sample}_qc_cells.tsv               # per cell: x, y, every metric, qc_keep, qc_drop_reason
│   ├── {sample}_qc_stats.json              # counts before/after, per-reason drops, medians, thresholds
│   ├── plots/{sample}/qc_metrics.png       # metric histograms with threshold lines
│   ├── plots/{sample}/spatial_qc.png       # counts / genes / verdict on the section
│   └── qc_summary_all.{tsv,png}            # one row per section
├── processed/
│   ├── {sample}.h5ad                       # log-normalized X, layers["counts"], X_pca, X_umap, obs["leiden"]
│   ├── {sample}_clusters.tsv               # cell, leiden, UMAP, x, y
│   ├── {sample}_markers.tsv                # top-N markers per cluster (score, log2FC, padj, % expressing)
│   └── plots/{sample}/{umap_clusters,spatial_clusters,markers_dotplot}.png
├── spatial/
│   ├── {sample}.h5ad                       # + obsp["spatial_connectivities"], uns["moranI"], nhood enrichment
│   └── {sample}/{nhood_enrichment.tsv,nhood_enrichment.png,morans_i.tsv,top_svg_spatial.png}
└── report/scstseq_report.html              # self-contained cross-sample report
```

Logs go to `logs/<rule>/<sample>.log`.

### Directory structure

```
snakemake_scstseq/                     # Snakemake Workflow Catalog layout
├── config/                            # config.yaml, samples.csv, README.md
├── workflow/
│   ├── Snakefile                      # entry point (unified DAG; targets: qc_all, process_all, spatial_all)
│   ├── rules/                         # common.smk, qc.smk, process.smk, spatial.smk, report.smk
│   ├── scripts/                       # one Python script per rule + scst_common.py helpers
│   ├── envs/py-scst.yaml              # the tool environment
│   └── schemas/config.schema.yaml     # parameter definitions (single source of truth)
├── .test/                             # executable test case (synthetic Xenium + Visium sections)
├── tests/                             # unit tests (pytest)
├── .snakemake-workflow-catalog.yml    # catalog metadata (enables snakedeploy)
├── Dockerfile, docker-compose.yml, apptainer.def, run_pipeline.sh, create_envs.smk, DOCKER.md
├── data/                              # vendor output directories (you provide)
├── results/, logs/                    # outputs
└── images/rulegraph.svg
```

## Tests & CI

- **Unit tests** (`pytest tests/ -q`): schema validity and shipped configs,
  sample-sheet checks, the QC threshold logic (drop-reason precedence,
  `0` = off), and both platform readers on in-memory vendor directories.
  They need only anndata / numpy / pandas / scipy / h5py / pyyaml / jsonschema.
- **Executable test case** ([`.test/`](.test)): a deterministic generator
  plants the truth — two spatial domains with disjoint marker programs, a
  known number of low-quality cells and high-mito spots — for one Xenium and
  two Visium sections, and `assert_outputs.py` asserts the **values** that
  come out (exact kept counts and drop reasons, domain-pure Leiden clusters,
  domain markers first in Moran's I, self-enriched neighborhoods).

  ```bash
  python .test/make_testdata.py
  snakemake -s workflow/Snakefile -d .test --sdm conda --cores 2
  python .test/assert_outputs.py
  ```
- **CI** (`.github/workflows/ci.yml`): a `static` job (snakefmt, pytest, dry
  runs of every target over stubbed inputs, lint) and an `integration` job
  that runs the executable test case end to end with a cached conda env and
  gates on the assertions.

## Versioning & releases

Commits follow [Conventional Commits](https://www.conventionalcommits.org);
[release-please](https://github.com/googleapis/release-please) turns them
into a release PR that bumps `version.txt`, `CITATION.cff` and
`CHANGELOG.md`, and tags a GitHub Release on merge. See
[`CONTRIBUTING.md`](CONTRIBUTING.md).

## Troubleshooting

- **Every cell/spot dropped** — `qc_filter` stops with the thresholds it
  used; relax `qc.<platform>` (`0` disables a rule). Xenium counts are
  transcripts per cell (tens to hundreds), Visium counts are UMIs per spot
  (thousands); do not reuse one platform's thresholds for the other.
- **Missing input file at DAG build** — the loader lists the exact vendor
  files it needs as rule inputs; the error names the missing path. Check
  `platform` and `path` in `samples.csv`.
- **One Leiden cluster** — markers and neighborhood enrichment need ≥ 2; raise
  `processing.<platform>.resolution`. Tiny pilot sections are handled (PCs
  and neighbours are capped to the data size).
- **Slow `spatial_stats` on large Xenium sections** — lower
  `spatial.moran.n_perms` (`0` keeps the analytic p-values only) and
  `spatial.nhood_enrichment.n_perms`, or restrict Moran's I to HVGs by
  enabling `processing.xenium.hvg`.

## Citation

If you use this workflow in your research, please cite it via the
**"Cite this repository"** button on GitHub (generated from
[`CITATION.cff`](CITATION.cff)).

**Please also cite the tools it builds on:**
- **Snakemake**: Mölder, F. et al. (2021). Sustainable data analysis with Snakemake. F1000Research, 10, 33.
- **scanpy**: Wolf, F.A., Angerer, P. and Theis, F.J. (2018). SCANPY: large-scale single-cell gene expression data analysis. Genome Biology, 19, 15.
- **squidpy**: Palla, G. et al. (2022). Squidpy: a scalable framework for spatial omics analysis. Nature Methods, 19, 171–178.
- **Leiden**: Traag, V.A., Waltman, L. and van Eck, N.J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. Scientific Reports, 9, 5233.
- **UMAP**: McInnes, L., Healy, J. and Melville, J. (2018). UMAP: Uniform Manifold Approximation and Projection for dimension reduction. arXiv:1802.03426.
- **anndata**: Virshup, I. et al. (2024). anndata: Access and store annotated data matrices. Journal of Open Source Software, 9(101), 4371.

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file.

## Contact

**Author**: gynecoloji
**Project repository**: [https://github.com/gynecoloji/snakemake_scstseq](https://github.com/gynecoloji/snakemake_scstseq)

For questions, issues, or feature requests, please check the existing
[Issues](https://github.com/gynecoloji/snakemake_scstseq/issues), then open a
new one with the exact `snakemake` command, the platform, and the relevant
`logs/` file.
