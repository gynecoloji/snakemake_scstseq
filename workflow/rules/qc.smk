# QC stage: per section, read the vendor output directory into an AnnData
# with spatial coordinates, compute per-cell/spot metrics, apply the
# platform-specific thresholds (recording the first failing rule as the drop
# reason), draw the QC plots, then aggregate a cross-sample summary.


rule qc_all:
    input:
        expand(f"{RESULTS}/qc/{{sample}}.filtered.h5ad", sample=SAMPLES),
        expand(f"{RESULTS}/qc/plots/{{sample}}/qc_metrics.png", sample=SAMPLES),
        expand(f"{RESULTS}/qc/plots/{{sample}}/spatial_qc.png", sample=SAMPLES),
        f"{RESULTS}/qc/qc_summary_all.tsv",
        f"{RESULTS}/qc/qc_summary_all.png",


# Rule 1: platform reader. Xenium: cell_feature_matrix.h5 + cells.parquet/csv
# (control probes/codewords are counted per cell, then removed from the gene
# matrix). Visium: filtered_feature_bc_matrix.h5 + spatial/ (positions, scale
# factors, images when present; out-of-tissue spots dropped). Barcodes are
# prefixed {sample}_ so later merges never collide.
rule load_sample:
    input:
        platform_input_files,
    output:
        h5ad=f"{RESULTS}/raw/{{sample}}.h5ad",
    log:
        f"{LOGS}/load_sample/{{sample}}.log",
    conda:
        "../envs/py-scst.yaml"
    threads: config["threads"]["load_sample"]
    params:
        platform=lambda w: sample_platform(w.sample),
        path=lambda w: sample_dir(w.sample),
        condition=lambda w: sample_condition(w.sample),
    script:
        "../scripts/load_sample.py"


# Rule 2: per-cell/spot metrics (total counts, detected genes, % mito / ribo /
# hemoglobin; Xenium adds cell/nucleus area and the control fraction), the
# platform thresholds from config (0 = off; first failing rule = drop reason),
# and the gene filter. Fails with a clear message if nothing survives.
rule qc_filter:
    input:
        h5ad=f"{RESULTS}/raw/{{sample}}.h5ad",
    output:
        h5ad=f"{RESULTS}/qc/{{sample}}.filtered.h5ad",
        cells=f"{RESULTS}/qc/{{sample}}_qc_cells.tsv",
        stats=f"{RESULTS}/qc/{{sample}}_qc_stats.json",
    log:
        f"{LOGS}/qc_filter/{{sample}}.log",
    conda:
        "../envs/py-scst.yaml"
    threads: config["threads"]["qc_filter"]
    params:
        platform=lambda w: sample_platform(w.sample),
        thresholds=platform_cfg("qc"),
        patterns=config["genome"],
    script:
        "../scripts/qc_filter.py"


# Rule 3: QC plots from the per-cell table alone (no matrix reload): metric
# histograms with the threshold lines, and the section coloured by counts,
# genes and kept/dropped.
rule qc_plots:
    input:
        cells=f"{RESULTS}/qc/{{sample}}_qc_cells.tsv",
    output:
        metrics=f"{RESULTS}/qc/plots/{{sample}}/qc_metrics.png",
        spatial=f"{RESULTS}/qc/plots/{{sample}}/spatial_qc.png",
    log:
        f"{LOGS}/qc_plots/{{sample}}.log",
    conda:
        "../envs/py-scst.yaml"
    threads: config["threads"]["qc_plots"]
    params:
        platform=lambda w: sample_platform(w.sample),
        thresholds=platform_cfg("qc"),
    script:
        "../scripts/qc_plots.py"


# Rule 4: cross-sample summary — one row per section (platform, condition,
# cells before/after, per-reason drop counts, medians) + a kept/dropped bar
# plot.
rule qc_summary:
    input:
        stats=expand(f"{RESULTS}/qc/{{sample}}_qc_stats.json", sample=SAMPLES),
    output:
        tsv=f"{RESULTS}/qc/qc_summary_all.tsv",
        png=f"{RESULTS}/qc/qc_summary_all.png",
    log:
        f"{LOGS}/qc_summary/all.log",
    conda:
        "../envs/py-scst.yaml"
    threads: config["threads"]["qc_summary"]
    script:
        "../scripts/qc_summary.py"
