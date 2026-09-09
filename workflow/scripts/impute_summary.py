"""Rule `impute_summary`: one row per section x method from the run JSONs."""

import json

import pandas as pd

from scst_common import log, redirect_log, write_tsv

LEAD = [
    "sample_id",
    "platform",
    "method",
    "genes_mode",
    "n_cells",
    "n_genes_imputed",
    "pct_zero_observed_mean",
    "pct_zero_imputed_mean",
    "median_gene_pearson_r",
    "runtime_s",
]


def summary_table(run_files):
    rows = []
    for p in run_files:
        with open(p) as fh:
            d = json.load(fh)
        d["params"] = json.dumps(d.get("params", {}), sort_keys=True)
        rows.append(d)
    df = pd.DataFrame(rows)
    rest = [c for c in df.columns if c not in LEAD and c != "params"]
    return df[[c for c in LEAD if c in df.columns] + rest + ["params"]]


def run(run_files, out_tsv):
    df = summary_table(run_files)
    write_tsv(df, out_tsv, index=False)
    log(f"summarised {len(df)} imputation runs → {out_tsv}")


if "snakemake" in globals():  # pragma: no cover
    redirect_log(snakemake.log[0])  # noqa: F821
    run(list(snakemake.input.runs), snakemake.output.tsv)  # noqa: F821
