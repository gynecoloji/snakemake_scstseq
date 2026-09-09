# Shared setup for every rule file: config validation, the sample sheet,
# config-derived shortcuts and the per-platform input-file helpers. Included
# first by workflow/Snakefile, so every name defined here is visible to the
# rules in the other included files.

import os
import re

import pandas as pd
from snakemake.utils import validate

# ── Config validation (also fills defaults) ─────────────────────────────
# Path is relative to this file (workflow/rules/).
validate(config, "../schemas/config.schema.yaml")

# ── Samples (from the sample sheet) ─────────────────────────────────────
PLATFORMS = ("xenium", "visium")

samples_df = pd.read_csv(config["samples_table"], dtype=str).fillna("")
samples_df.columns = [c.strip() for c in samples_df.columns]
for col in ("sample_id", "platform", "path"):
    if col not in samples_df.columns:
        raise ValueError(f"{config['samples_table']}: missing required column '{col}'")
for col in samples_df.columns:
    samples_df[col] = samples_df[col].astype(str).str.strip()
if "condition" not in samples_df.columns:
    samples_df["condition"] = ""

SAMPLES = samples_df["sample_id"].tolist()
if not SAMPLES:
    raise ValueError(f"{config['samples_table']}: no samples listed")
if len(set(SAMPLES)) != len(SAMPLES):
    dups = sorted({s for s in SAMPLES if SAMPLES.count(s) > 1})
    raise ValueError(f"{config['samples_table']}: duplicate sample_id(s): {dups}")
_bad = samples_df[~samples_df["platform"].isin(PLATFORMS)]
if len(_bad):
    raise ValueError(
        f"{config['samples_table']}: unknown platform for "
        + ", ".join(f"{r.sample_id} ({r.platform!r})" for r in _bad.itertuples())
        + f"; expected one of {PLATFORMS}"
    )
_BY_SAMPLE = samples_df.set_index("sample_id")

# ── Config shortcuts ────────────────────────────────────────────────────
DATA = config["paths"]["data_dir"]
RESULTS = config["paths"]["results_dir"]
LOGS = config["paths"]["logs_dir"]

# Helper scripts live next to the Snakefile, NOT under the working directory,
# so `snakemake -d <dir>` (the executable test case, the catalog) and
# module/snakedeploy imports resolve them.
SCRIPTS = os.path.join(workflow.basedir, "scripts")

REPORT_TARGETS = (
    [f"{RESULTS}/report/scstseq_report.html"] if config["report"]["enabled"] else []
)

# ── Imputation (opt-in stage; see rules/impute.smk) ──────────────────────
IMPUTE_DIR = f"{RESULTS}/imputed"
IMPUTE_METHODS = list(config["imputation"]["methods"])


# ── Per-sample helpers ──────────────────────────────────────────────────
def sample_platform(sample):
    return _BY_SAMPLE.loc[sample, "platform"]


def sample_condition(sample):
    return _BY_SAMPLE.loc[sample, "condition"]


def sample_dir(sample):
    """Vendor output directory; sample-sheet paths are relative to the data dir unless absolute."""
    path = _BY_SAMPLE.loc[sample, "path"]
    return path if os.path.isabs(path) else os.path.join(DATA, path)


def _first_existing(base, candidates):
    """The first candidate that exists, else the first one (so the missing-input
    error names the canonical file)."""
    for rel in candidates:
        if os.path.exists(os.path.join(base, rel)):
            return os.path.join(base, rel)
    return os.path.join(base, candidates[0])


def platform_input_files(wildcards):
    """The concrete vendor files a sample's loader reads. Listing them as the
    rule's input makes a missing file fail at DAG-build time, with the path
    named, instead of deep inside a job."""
    base = sample_dir(wildcards.sample)
    platform = sample_platform(wildcards.sample)
    if platform == "xenium":
        return [
            os.path.join(base, "cell_feature_matrix.h5"),
            _first_existing(base, ["cells.parquet", "cells.csv.gz", "cells.csv"]),
        ]
    return [
        os.path.join(base, "filtered_feature_bc_matrix.h5"),
        os.path.join(base, "spatial", "scalefactors_json.json"),
        _first_existing(
            base,
            [
                os.path.join("spatial", "tissue_positions.csv"),
                os.path.join("spatial", "tissue_positions_list.csv"),
            ],
        ),
    ]


def platform_cfg(section):
    """`params:` helper: the platform-specific block of a config section."""
    return lambda wildcards: config[section][sample_platform(wildcards.sample)]


# Constrain wildcards to known values so no stray files are matched.
wildcard_constraints:
    sample="|".join(re.escape(s) for s in SAMPLES) if SAMPLES else "a^",
