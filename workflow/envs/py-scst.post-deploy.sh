#!/usr/bin/env bash
# Snakemake post-deployment script for workflow/envs/py-scst.yaml: runs once,
# right after the conda env is created, with $CONDA_PREFIX pointing at it.
#
# Installs MAGIC (magic-impute) and its pure-Python chain from PyPI WITHOUT
# dependency resolution. Reason: the bioconda magic-impute build depends on a
# conda `scprep` that is years behind, and PyPI scprep 1.2.3 pins
# pandas < 2.1 (no Python 3.12 wheel), so a normal `pip install` tries to
# compile pandas and fails — while MAGIC itself runs fine on the current
# numpy / pandas / scikit-learn (smoke-tested on numpy 2.5, pandas 3.0).
# The runtime deps these wheels need (numpy, scipy, scikit-learn, pandas,
# matplotlib, decorator, packaging, future, deprecated) come from the conda
# spec.
#
# NOTE: Snakemake hashes only the .yaml to decide whether an env is current;
# editing this script alone does not trigger a rebuild (bump a comment in the
# yaml to force one).
set -euo pipefail
"${CONDA_PREFIX}/bin/pip" install --no-deps --no-cache-dir \
    magic-impute==3.0.0 scprep==1.2.3 graphtools==2.1.0 tasklogger==1.2.0 PyGSP==0.6.1
"${CONDA_PREFIX}/bin/python" -c "import magic, scprep, graphtools; print('magic-impute', magic.__version__, 'installed')"
