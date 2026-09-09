# Build-time helper: pre-create the per-rule conda envs (workflow/envs/*.yaml)
# so they are baked into the Docker/Apptainer image and reused at runtime via
# --conda-prefix.
#
# Used ONLY by the Dockerfile / apptainer.def / apptainer-gpu.def:
#   snakemake -s create_envs.smk --use-conda --conda-create-envs-only \
#       --conda-frontend mamba --conda-prefix /opt/wf-conda --cores 1
#
# It has no external inputs, so the DAG resolves without any data. Snakemake
# keys each conda env by the CONTENT of its workflow/envs/*.yaml file, so the
# envs built here are reused by workflow/Snakefile at runtime (same env files,
# same --conda-prefix). If content differs, Snakemake just rebuilds that one
# env at runtime — no hard failure.
#
# Default env set = the CPU image (py-scst + the CPU scVI env used by the
# opt-in imputation stage). apptainer-gpu.def overrides via WF_CONDA_ENVS to
# swap py-scvi for its CUDA build (py-scvi-gpu).
import os

ENVS = os.environ.get("WF_CONDA_ENVS", "py-scst,py-scvi").split(",")


rule all:
    input:
        expand("build/conda_env_{env}.ready", env=ENVS),


rule create_env:
    output:
        "build/conda_env_{env}.ready",
    conda:
        "workflow/envs/{env}.yaml"
    shell:
        "touch {output}"
