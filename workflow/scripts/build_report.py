"""Rule `build_report`: one self-contained HTML page over every section.

PNGs are base64-embedded and tables inlined, so the file has no external
assets and no JavaScript — it can be mailed or archived as-is.
"""


import base64
import datetime as dt
import html
import os

import pandas as pd

from scst_common import log, redirect_log

CSS = """
body { font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; margin: 0; padding: 0 0 4rem 0;
       color: #1f2933; background: #fafbfc; }
main { max-width: 1240px; margin: 0 auto; padding: 1.5rem; }
h1 { font-size: 1.7rem; margin: 0.4rem 0 0.2rem 0; }
h2 { font-size: 1.3rem; border-bottom: 2px solid #d9e2ec; padding-bottom: 0.25rem; margin-top: 2.4rem; }
h3 { font-size: 1.05rem; margin: 1.4rem 0 0.4rem 0; color: #334e68; }
.meta { color: #52606d; font-size: 0.9rem; }
.grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: 1rem; }
figure { margin: 0; background: #fff; border: 1px solid #e4e7eb; border-radius: 6px; padding: 0.5rem; }
figure img { max-width: 100%; height: auto; display: block; }
figcaption { font-size: 0.8rem; color: #52606d; padding-top: 0.3rem; }
table { border-collapse: collapse; font-size: 0.8rem; background: #fff; margin: 0.4rem 0 1rem 0; }
th, td { border: 1px solid #e4e7eb; padding: 0.25rem 0.5rem; text-align: right; }
th { background: #f0f4f8; }
td:first-child, th:first-child { text-align: left; }
.wide { overflow-x: auto; }
.toc a { margin-right: 1rem; }
"""


def img(path, caption=""):
    with open(path, "rb") as fh:
        b64 = base64.b64encode(fh.read()).decode("ascii")
    cap = f"<figcaption>{html.escape(caption)}</figcaption>" if caption else ""
    return f'<figure><img src="data:image/png;base64,{b64}" alt="{html.escape(caption)}">{cap}</figure>'


def table(df, max_rows=None):
    if max_rows is not None:
        df = df.head(max_rows)
    return '<div class="wide">' + df.to_html(index=False, border=0, float_format=lambda x: f"{x:.4g}") + "</div>"


def markers_block(path, per_cluster=5):
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return "<p class='meta'>No markers (fewer than 2 clusters).</p>"
    top = df[df["rank"] <= per_cluster][["cluster", "rank", "gene", "score", "log2fc", "padj", "pct_in", "pct_out"]]
    return table(top)


def moran_block(path, n_top):
    df = pd.read_csv(path, sep="\t")
    cols = [c for c in ["rank", "gene", "I", "pval_norm", "pval_norm_fdr_bh", "pval_sim", "pval_sim_fdr_bh"] if c in df.columns]
    return table(df[cols], max_rows=n_top)


def read_version(path):
    try:
        with open(path) as fh:
            return fh.read().strip()
    except OSError:
        return "dev"


def build(inputs, samples, platforms, conditions, n_top_svg, version, out_html):
    summary = pd.read_csv(inputs["summary_tsv"], sep="\t")
    n_xen = sum(1 for s in samples if platforms[s] == "xenium")
    n_vis = len(samples) - n_xen
    parts = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        "<title>snakemake_scstseq report</title>",
        f"<style>{CSS}</style></head><body><main>",
        "<h1>Spatial transcriptomics report</h1>",
        f"<p class='meta'>snakemake_scstseq v{html.escape(version)} · generated {dt.datetime.now():%Y-%m-%d %H:%M} · "
        f"{len(samples)} sections ({n_xen} Xenium, {n_vis} Visium)</p>",
        "<p class='toc'>" + " ".join(f"<a href='#{html.escape(s)}'>{html.escape(s)}</a>" for s in samples) + "</p>",
        "<h2 id='summary'>Cross-sample QC summary</h2>",
        table(summary),
        img(inputs["summary_png"], "Kept / dropped cells or spots per section and median depth of the kept ones."),
    ]
    for i, s in enumerate(samples):
        platform = platforms[s]
        cond = conditions.get(s, "")
        unit = "cells" if platform == "xenium" else "spots"
        row = summary[summary["sample_id"] == s]
        kept = f"{int(row['n_cells_kept'].iloc[0])} / {int(row['n_cells_raw'].iloc[0])} {unit} kept" if len(row) else ""
        parts.append(f"<h2 id='{html.escape(s)}'>{html.escape(s)} <span class='meta'>· {platform}"
                     f"{' · ' + html.escape(cond) if cond else ''} · {kept}</span></h2>")
        parts.append("<h3>Quality control</h3><div class='grid'>")
        parts.append(img(inputs["qc_metrics"][i], "QC metric distributions (kept vs dropped) with the configured thresholds."))
        parts.append(img(inputs["qc_spatial"][i], "Counts, detected genes and the QC verdict on the section."))
        parts.append("</div><h3>Clusters</h3><div class='grid'>")
        parts.append(img(inputs["umap"][i], "Leiden clusters in UMAP space."))
        parts.append(img(inputs["clusters"][i], "Leiden clusters on the section."))
        parts.append("</div>")
        parts.append(img(inputs["dotplot"][i], "Leading marker genes per cluster (mean expression and fraction expressing)."))
        parts.append("<h3>Top markers per cluster</h3>")
        parts.append(markers_block(inputs["markers"][i]))
        parts.append("<h3>Spatial statistics</h3><div class='grid'>")
        parts.append(img(inputs["nhood_png"][i], "Neighborhood enrichment between clusters (permutation z-scores)."))
        parts.append(img(inputs["svg_png"][i], "Top spatially variable genes by Moran's I."))
        parts.append("</div><h3>Spatially variable genes (Moran's I)</h3>")
        parts.append(moran_block(inputs["moran"][i], n_top_svg))
    parts.append("</main></body></html>")
    os.makedirs(os.path.dirname(out_html) or ".", exist_ok=True)
    with open(out_html, "w", encoding="utf-8") as fh:
        fh.write("\n".join(parts))
    log(f"wrote {out_html} ({os.path.getsize(out_html) / 1e6:.1f} MB)")


if "snakemake" in globals():  # pragma: no cover
    redirect_log(snakemake.log[0])  # noqa: F821
    _inputs = {k: (list(v) if isinstance(v, (list, tuple)) else v) for k, v in dict(snakemake.input).items()}  # noqa: F821
    build(
        _inputs,
        list(snakemake.params.samples),  # noqa: F821
        dict(snakemake.params.platforms),  # noqa: F821
        dict(snakemake.params.conditions),  # noqa: F821
        int(snakemake.params.n_top_svg),  # noqa: F821
        read_version(snakemake.params.version_file),  # noqa: F821
        snakemake.output.html,  # noqa: F821
    )
