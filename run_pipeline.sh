#!/usr/bin/env bash
# Build (once) and run the spatial-transcriptomics Snakemake workflow in Docker.
#
# The workflow is the standard-layout workflow/Snakefile: a single run builds
# every stage (QC → processing → spatial → report) in dependency order. Pass
# extra snakemake args (cores, a target, -n, ...) through to the container.
#
# Usage:
#   ./run_pipeline.sh                          # everything, 4 cores
#   ./run_pipeline.sh --cores 16               # everything, 16 cores
#   ./run_pipeline.sh --cores 16 qc_all        # QC stage only
#   ./run_pipeline.sh --cores 16 spatial_all   # spatial statistics (after processing)
#   ./run_pipeline.sh -n                       # dry run: check the DAG first
#
# The current directory (code + data/ vendor directories) is mounted at
# /workflow; results/ are written back here. The pre-built conda env lives in
# the image at /opt/wf-conda and is reused via --conda-prefix. The image
# ENTRYPOINT runs `snakemake --use-conda ...`.
set -euo pipefail

IMAGE="scstseq-pipeline:latest"

if [ "$#" -eq 0 ]; then set -- --cores 4; fi   # default snakemake args

if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
    echo ">> Building $IMAGE (first time only; pre-builds the conda env, ~10-20 min)..."
    docker build -t "$IMAGE" .
fi

echo ">> snakemake -s workflow/Snakefile $*"
docker run --rm \
    -v "$(pwd)":/workflow \
    -e HOME=/tmp \
    --user "$(id -u):$(id -g)" \
    "$IMAGE" -s workflow/Snakefile "$@"
