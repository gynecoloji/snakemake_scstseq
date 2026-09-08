# Processing stage: per section, normalize (platform-specific target), select
# HVGs when configured, PCA → neighbours → Leiden → UMAP, then per-cluster
# marker genes. Counts are kept in layers["counts"]; X holds the log-
# normalized values used for markers and plots.


rule process_all:
    input:
        expand(f"{RESULTS}/processed/{{sample}}.h5ad", sample=SAMPLES),
        expand(f"{RESULTS}/processed/{{sample}}_clusters.tsv", sample=SAMPLES),
        expand(f"{RESULTS}/processed/{{sample}}_markers.tsv", sample=SAMPLES),
        expand(
            f"{RESULTS}/processed/plots/{{sample}}/umap_clusters.png", sample=SAMPLES
        ),
        expand(
            f"{RESULTS}/processed/plots/{{sample}}/spatial_clusters.png",
            sample=SAMPLES,
        ),
        expand(
            f"{RESULTS}/processed/plots/{{sample}}/markers_dotplot.png", sample=SAMPLES
        ),


# Rule 5: normalize_total(target_sum) + log1p; HVGs (seurat_v3 on the counts
# layer) when `hvg`; optional scaling on a temporary copy so X stays log-
# normalized; PCA / neighbours / Leiden (igraph flavour) / UMAP with n_pcs and
# n_neighbors capped to the data size so tiny sections still run.
rule process_sample:
    input:
        h5ad=f"{RESULTS}/qc/{{sample}}.filtered.h5ad",
    output:
        h5ad=f"{RESULTS}/processed/{{sample}}.h5ad",
        clusters=f"{RESULTS}/processed/{{sample}}_clusters.tsv",
        umap=f"{RESULTS}/processed/plots/{{sample}}/umap_clusters.png",
        spatial=f"{RESULTS}/processed/plots/{{sample}}/spatial_clusters.png",
    log:
        f"{LOGS}/process_sample/{{sample}}.log",
    conda:
        "../envs/py-scst.yaml"
    threads: config["threads"]["process_sample"]
    params:
        platform=lambda w: sample_platform(w.sample),
        cfg=platform_cfg("processing"),
        seed=config["processing"]["seed"],
    script:
        "../scripts/process_sample.py"


# Rule 6: rank_genes_groups over the Leiden clusters on the log-normalized
# values; top-N table per cluster and a dotplot of the leading markers.
rule cluster_markers:
    input:
        h5ad=f"{RESULTS}/processed/{{sample}}.h5ad",
    output:
        tsv=f"{RESULTS}/processed/{{sample}}_markers.tsv",
        dotplot=f"{RESULTS}/processed/plots/{{sample}}/markers_dotplot.png",
    log:
        f"{LOGS}/cluster_markers/{{sample}}.log",
    conda:
        "../envs/py-scst.yaml"
    threads: config["threads"]["cluster_markers"]
    params:
        cfg=config["processing"]["markers"],
    script:
        "../scripts/cluster_markers.py"
