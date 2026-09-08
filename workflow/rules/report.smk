# Report: one self-contained HTML page (PNGs base64-embedded, tables inline,
# no external assets, no JavaScript) with the cross-sample QC summary and,
# per section, the QC plots, cluster UMAP/spatial plots, marker dotplot,
# neighborhood-enrichment heatmap and the top spatially variable genes.


rule build_report:
    input:
        summary_tsv=f"{RESULTS}/qc/qc_summary_all.tsv",
        summary_png=f"{RESULTS}/qc/qc_summary_all.png",
        qc_metrics=expand(
            f"{RESULTS}/qc/plots/{{sample}}/qc_metrics.png", sample=SAMPLES
        ),
        qc_spatial=expand(
            f"{RESULTS}/qc/plots/{{sample}}/spatial_qc.png", sample=SAMPLES
        ),
        umap=expand(
            f"{RESULTS}/processed/plots/{{sample}}/umap_clusters.png", sample=SAMPLES
        ),
        clusters=expand(
            f"{RESULTS}/processed/plots/{{sample}}/spatial_clusters.png",
            sample=SAMPLES,
        ),
        dotplot=expand(
            f"{RESULTS}/processed/plots/{{sample}}/markers_dotplot.png", sample=SAMPLES
        ),
        markers=expand(f"{RESULTS}/processed/{{sample}}_markers.tsv", sample=SAMPLES),
        nhood_png=expand(
            f"{RESULTS}/spatial/{{sample}}/nhood_enrichment.png", sample=SAMPLES
        ),
        moran=expand(f"{RESULTS}/spatial/{{sample}}/morans_i.tsv", sample=SAMPLES),
        svg_png=expand(
            f"{RESULTS}/spatial/{{sample}}/top_svg_spatial.png", sample=SAMPLES
        ),
    output:
        html=f"{RESULTS}/report/scstseq_report.html",
    log:
        f"{LOGS}/build_report/all.log",
    conda:
        "../envs/py-scst.yaml"
    threads: config["threads"]["build_report"]
    params:
        samples=SAMPLES,
        platforms={s: sample_platform(s) for s in SAMPLES},
        conditions={s: sample_condition(s) for s in SAMPLES},
        n_top_svg=config["spatial"]["moran"]["n_top"],
        version_file=os.path.join(workflow.basedir, "..", "version.txt"),
    script:
        "../scripts/build_report.py"
