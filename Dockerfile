# Spatial transcriptomics (Xenium + Visium) — Snakemake workflow container.
#
# The workflow runs one conda env PER RULE (workflow/envs/*.yaml). This image
# ships Snakemake + the pre-built per-rule conda env (py-scst: scanpy + squidpy)
# and runs --use-conda.
#
# Vendor output directories are NOT baked in — mount your project directory at
# runtime (see docker-compose.yml / run_pipeline.sh / DOCKER.md).

# For a fully reproducible build, pin the base, e.g. FROM condaforge/miniforge3:24.11.3-2
FROM condaforge/miniforge3:latest

LABEL org.opencontainers.image.title="scstseq-snakemake"
LABEL org.opencontainers.image.description="Spatial transcriptomics (Xenium + Visium) QC, clustering and spatial statistics Snakemake workflow, --use-conda"
LABEL org.opencontainers.image.source="https://github.com/gynecoloji/snakemake_scstseq"

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    WF_CONDA_PREFIX=/opt/wf-conda

# Disable apt's privilege-dropping sandbox so the image also builds under a rootless
# engine (e.g. podman without /etc/subuid), mirroring apptainer.def. No-op under root
# docker (apt would run as root anyway). git helps Snakemake's conda handling; procps
# for subprocess management.
RUN printf 'APT::Sandbox::User "root";\n' > /etc/apt/apt.conf.d/01-no-sandbox && \
    apt-get update && apt-get install -y --no-install-recommends \
        git procps ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# FLEXIBLE channel priority (safe for the conda-forge + bioconda env spec);
# best-effort accept the Anaconda defaults ToS so a non-interactive solve never stalls.
RUN conda config --system --set channel_priority flexible && \
    ( conda tos accept --override-channels \
        --channel https://repo.anaconda.com/pkgs/main \
        --channel https://repo.anaconda.com/pkgs/r 2>/dev/null || true )

# Snakemake driver in its OWN env (miniforge base pins Python 3.13; snakemake-minimal
# 9.3.2 needs <3.13). pandas is imported by common.smk at parse time. Prepend the
# driver env to PATH.
RUN mamba create -y -n driver -c conda-forge -c bioconda \
        python=3.12 snakemake-minimal=9.3.2 pandas && \
    mamba clean -afy
ENV PATH=/opt/conda/envs/driver/bin:$PATH

WORKDIR /workflow

# Env spec first, then pre-build the per-rule conda env INTO the image. Doing this
# BEFORE copying the workflow code means later edits to the rules/scripts don't
# invalidate the (slow) conda-env layer.
COPY workflow/envs/ ./workflow/envs/
COPY create_envs.smk ./
RUN snakemake -s create_envs.smk --use-conda --conda-create-envs-only \
        --conda-frontend mamba --conda-prefix "${WF_CONDA_PREFIX}" --cores 1 && \
    mamba clean -afy && \
    rm -rf build .snakemake

# Workflow code + config (data are mounted at runtime, not baked)
COPY workflow/ ./workflow/
COPY config/ ./config/
COPY tests/ ./tests/
COPY version.txt ./

# ENTRYPOINT fixes the conda settings; pass the target/cores at `docker run`.
ENTRYPOINT ["snakemake", "--use-conda", "--conda-frontend", "mamba", "--conda-prefix", "/opt/wf-conda"]
CMD ["-s", "workflow/Snakefile", "--cores", "4"]
