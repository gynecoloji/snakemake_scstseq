# Imputation stage (OPT-IN: `snakemake impute_all`): denoise / fill in the
# dropouts of the measured genes, one output set per method listed in
# config.imputation.methods. Every method writes on the log-normalized scale:
#
#   magic  diffusion over the expression kNN graph (magic-impute, py-scst env)
#   alra   adaptively-thresholded low-rank approximation (numpy re-implementation)
#   scvi   scvi-tools model per section, get_normalized_expression -> log1p
#          (own conda env: py-scvi, or py-scvi-gpu with imputation.scvi.gpu)
#
# Inputs are the processed object (clusters, HVGs, counts layer) and the
# Moran's I table (to pick which genes to draw observed vs imputed). Outputs
# under results/imputed/: {sample}.{method}.h5ad restricted to the imputed
# genes with layers["imputed"], a per-gene statistics table, a run JSON and
# the comparison figure; impute_summary collects one row per section x method.


rule impute_all:
    input:
        expand(
            f"{IMPUTE_DIR}/{{sample}}.{{method}}.h5ad",
            sample=SAMPLES,
            method=IMPUTE_METHODS,
        ),
        expand(
            f"{IMPUTE_DIR}/{{sample}}.{{method}}_gene_stats.tsv",
            sample=SAMPLES,
            method=IMPUTE_METHODS,
        ),
        expand(
            f"{IMPUTE_DIR}/plots/{{sample}}.{{method}}_observed_vs_imputed.png",
            sample=SAMPLES,
            method=IMPUTE_METHODS,
        ),
        [f"{IMPUTE_DIR}/imputation_summary.tsv"] if IMPUTE_METHODS else [],


_IMPUTE_OUTPUTS = dict(
    h5ad=f"{IMPUTE_DIR}/{{sample}}.{{method}}.h5ad",
    stats=f"{IMPUTE_DIR}/{{sample}}.{{method}}_gene_stats.tsv",
    run=f"{IMPUTE_DIR}/{{sample}}.{{method}}_run.json",
    png=f"{IMPUTE_DIR}/plots/{{sample}}.{{method}}_observed_vs_imputed.png",
)


# Rule 8a: MAGIC and ALRA share the py-scst env and one dispatching script.
rule impute:
    input:
        h5ad=f"{RESULTS}/processed/{{sample}}.h5ad",
        moran=f"{RESULTS}/spatial/{{sample}}/morans_i.tsv",
    output:
        **_IMPUTE_OUTPUTS,
    log:
        f"{LOGS}/impute/{{sample}}.{{method}}.log",
    wildcard_constraints:
        method="magic|alra",
    conda:
        "../envs/py-scst.yaml"
    threads: config["threads"]["impute"]
    params:
        platform=lambda w: sample_platform(w.sample),
        cfg=lambda w: config["imputation"][w.method],
        genes=config["imputation"]["genes"],
        plot_top=config["imputation"]["plot_top"],
        seed=config["processing"]["seed"],
    script:
        "../scripts/impute_sample.py"


# The scVI env comes in a CPU and a CUDA flavour; imputation.scvi.gpu picks
# the latter (baked into scstseq-pipeline-gpu.sif; run with `apptainer --nv`).
_SCVI_ENV = (
    "../envs/py-scvi-gpu.yaml"
    if config["imputation"]["scvi"]["gpu"]
    else "../envs/py-scvi.yaml"
)


# Rule 8b: scVI, in its own environment.
rule impute_scvi:
    input:
        h5ad=f"{RESULTS}/processed/{{sample}}.h5ad",
        moran=f"{RESULTS}/spatial/{{sample}}/morans_i.tsv",
    output:
        **_IMPUTE_OUTPUTS,
    log:
        f"{LOGS}/impute/{{sample}}.{{method}}.log",
    wildcard_constraints:
        method="scvi",
    conda:
        _SCVI_ENV
    threads: config["threads"]["impute"]
    params:
        platform=lambda w: sample_platform(w.sample),
        cfg=config["imputation"]["scvi"],
        genes=config["imputation"]["genes"],
        plot_top=config["imputation"]["plot_top"],
        seed=config["processing"]["seed"],
    script:
        "../scripts/impute_scvi.py"


# Rule 9: one row per section x method (genes imputed, zero fraction before /
# after, median observed-vs-imputed gene correlation, runtime, method extras).
rule impute_summary:
    input:
        runs=expand(
            f"{IMPUTE_DIR}/{{sample}}.{{method}}_run.json",
            sample=SAMPLES,
            method=IMPUTE_METHODS,
        ),
    output:
        tsv=f"{IMPUTE_DIR}/imputation_summary.tsv",
    log:
        f"{LOGS}/impute/summary.log",
    conda:
        "../envs/py-scst.yaml"
    threads: config["threads"]["impute_summary"]
    script:
        "../scripts/impute_summary.py"
