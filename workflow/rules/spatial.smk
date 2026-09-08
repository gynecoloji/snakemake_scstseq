# Spatial stage: per section, build the spatial neighbours graph (kNN /
# Delaunay / radius over Xenium cell centroids in µm; hexagonal grid rings
# over Visium spots), then neighborhood enrichment between the Leiden
# clusters and Moran's I spatially variable genes (squidpy).


rule spatial_all:
    input:
        expand(f"{RESULTS}/spatial/{{sample}}.h5ad", sample=SAMPLES),
        expand(f"{RESULTS}/spatial/{{sample}}/nhood_enrichment.tsv", sample=SAMPLES),
        expand(f"{RESULTS}/spatial/{{sample}}/nhood_enrichment.png", sample=SAMPLES),
        expand(f"{RESULTS}/spatial/{{sample}}/morans_i.tsv", sample=SAMPLES),
        expand(f"{RESULTS}/spatial/{{sample}}/top_svg_spatial.png", sample=SAMPLES),


# Rule 7: sq.gr.spatial_neighbors with the platform block of config.spatial;
# sq.gr.nhood_enrichment on `leiden` (z-score matrix → TSV + heatmap);
# sq.gr.spatial_autocorr (Moran's I) over the HVGs when they exist, else all
# genes, with permutation p-values (n_perms = 0 → analytic only) and BH-FDR;
# the top SVGs drawn on the section. The final AnnData carries the graph.
rule spatial_stats:
    input:
        h5ad=f"{RESULTS}/processed/{{sample}}.h5ad",
    output:
        h5ad=f"{RESULTS}/spatial/{{sample}}.h5ad",
        nhood=f"{RESULTS}/spatial/{{sample}}/nhood_enrichment.tsv",
        nhood_png=f"{RESULTS}/spatial/{{sample}}/nhood_enrichment.png",
        moran=f"{RESULTS}/spatial/{{sample}}/morans_i.tsv",
        svg_png=f"{RESULTS}/spatial/{{sample}}/top_svg_spatial.png",
    log:
        f"{LOGS}/spatial_stats/{{sample}}.log",
    conda:
        "../envs/py-scst.yaml"
    threads: config["threads"]["spatial_stats"]
    params:
        platform=lambda w: sample_platform(w.sample),
        graph=platform_cfg("spatial"),
        nhood=config["spatial"]["nhood_enrichment"],
        moran=config["spatial"]["moran"],
        seed=config["processing"]["seed"],
    script:
        "../scripts/spatial_stats.py"
