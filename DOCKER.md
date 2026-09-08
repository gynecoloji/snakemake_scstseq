# Running the spatial-transcriptomics workflow in a container

The pipeline runs its rules in **one conda environment** (`workflow/envs/py-scst.yaml`:
scanpy + squidpy). The image therefore ships **Snakemake + that pre-built env**
and runs Snakemake with `--use-conda`.

Vendor output directories are **not** baked into the image — you mount your
project directory at run time. **Apptainer is the primary path** (this workflow
is meant to run on HPC from a `.sif`); Docker is the fallback for workstations.

## 1. What you need on the host (in `data/`)

The container reads these from the mounted project directory:

| Path | What |
|---|---|
| `data/<sample>/…` (Xenium) | `cell_feature_matrix.h5` + `cells.parquet` / `cells.csv.gz` from the Xenium `outs/` |
| `data/<sample>/…` (Visium) | `filtered_feature_bc_matrix.h5` + `spatial/` from the Space Ranger `outs/` |
| `config/config.yaml`, `config/samples.csv` | config + sample sheet (tracked in the repo) |

List your sections (with `platform` and `path`) in `config/samples.csv`.

## 2. Get the image

### Pull the published image (no build)

```bash
apptainer pull scstseq-pipeline.sif docker://gynecoloji/scstseq-pipeline:latest   # HPC
docker pull gynecoloji/scstseq-pipeline:latest                                     # Docker
```

### Or build it yourself

**Apptainer (HPC):**

```bash
module load apptainer                       # if your cluster uses modules
cd /path/to/snakemake_scstseq               # must run from the repo root (%files context)
apptainer build --fakeroot scstseq-pipeline.sif apptainer.def
```

Because this is a long, memory-hungry build, run it as a batch job rather than
on a login node, and point Apptainer's scratch at a large filesystem:

```bash
export APPTAINER_TMPDIR=/big/scratch/apptainer_tmp
export APPTAINER_CACHEDIR=/big/scratch/apptainer_cache
```

**Docker:**

```bash
docker compose build
# or:  docker build -t scstseq-pipeline:latest .
```

This pre-builds the conda env into the image (a few GB; squidpy pulls in
spatialdata/dask). For a reproducible image, pin the base tag in the
`Dockerfile` / `apptainer.def` (`condaforge/miniforge3:<version>`).

## 3. Run

The image `ENTRYPOINT` is `snakemake --use-conda --conda-frontend mamba
--conda-prefix /opt/wf-conda`; everything after the image name goes to
`snakemake`. A single run builds every stage (unified DAG).

**Apptainer** (auto-mounts `$HOME`, `/tmp` and the current directory; runs as you):

```bash
apptainer run scstseq-pipeline.sif -s workflow/Snakefile --cores 1 -n         # dry run first
apptainer run scstseq-pipeline.sif -s workflow/Snakefile --cores 16           # everything
apptainer run scstseq-pipeline.sif -s workflow/Snakefile --cores 16 qc_all    # one stage
```

**Docker** — the helper script or compose:

```bash
./run_pipeline.sh -n                          # dry run
./run_pipeline.sh --cores 16                  # everything
./run_pipeline.sh --cores 16 spatial_all      # one stage
docker compose run --rm scstseq -s workflow/Snakefile --cores 16
```

Notes:
- **Data outside the project dir:** bind it in —
  `apptainer run --bind /scratch/xenium:/scratch/xenium …` /
  `docker run -v /scratch/xenium:/scratch/xenium …` — and use absolute paths
  in `config/samples.csv`.
- **Pre-built env:** the conda environment is baked at `/opt/wf-conda`
  (read-only in the SIF) and reused via `--conda-prefix`. If Apptainer
  reports a read-only error writing there (an edited `py-scst.yaml` forces a
  rebuild), add `--writable-tmpfs`, or pass a writable
  `--conda-prefix /path/on/scratch` after the image name.
- **Outputs as your user:** the Docker wrapper passes `--user "$(id -u):$(id -g)"`
  and `HOME=/tmp` so `results/` is not root-owned.

## 4. Publish (maintainers)

```bash
docker build -t gynecoloji/scstseq-pipeline:latest .
docker push gynecoloji/scstseq-pipeline:latest
```

Apptainer users pull that same image with `apptainer pull docker://…`; there
is no separate SIF upload.
