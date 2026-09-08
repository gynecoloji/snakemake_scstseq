# Contributing

Thanks for your interest in improving this workflow! This guide covers how to
report problems, set up a development environment, run the tests, and propose
changes.

By participating, you agree to abide by the [Code of Conduct](CODE_OF_CONDUCT.md).

## Reporting issues

Open a [GitHub issue](https://github.com/gynecoloji/snakemake_scstseq/issues)
with:

- what you ran (the exact `snakemake` command / target),
- what you expected vs. what happened,
- the relevant log(s) from `logs/`, your platform (Xenium / Visium), and your
  OS + Snakemake version.

## Development setup

```bash
git clone https://github.com/gynecoloji/snakemake_scstseq.git
cd snakemake_scstseq

# Driver environment (the per-rule tool env is created on the first --use-conda run)
mamba create -n scstseq -c conda-forge -c bioconda snakemake-minimal=9.3.2 pandas pytest pyyaml jsonschema
conda activate scstseq
```

Vendor output directories (Xenium `outs/`, Space Ranger `outs/`) are not
tracked; place them under `data/` and list them in `config/samples.csv` as
described in [`config/README.md`](config/README.md) and the top-level
`README.md`.

## Running the workflow

```bash
snakemake -s workflow/Snakefile -n                        # dry run (validates config + DAG)
snakemake -s workflow/Snakefile --use-conda --cores 8     # everything (QC → processing → spatial → report)
snakemake -s workflow/Snakefile --use-conda --cores 8 qc_all   # one stage
```

Configuration is validated against [`workflow/schemas/config.schema.yaml`](workflow/schemas/config.schema.yaml)
on every run — that schema is the single source of truth for parameters.

## Tests

Unit tests need only `anndata`, `numpy`, `pandas`, `scipy`, `h5py`, `pyyaml`
and `jsonschema`, and run in seconds:

```bash
python -m pytest tests/ -q
```

The executable test case in [`.test/`](.test) runs the whole workflow on a
synthetic dataset and asserts the numbers that come out:

```bash
python .test/make_testdata.py
snakemake -s workflow/Snakefile -d .test --sdm conda --cores 4
python .test/assert_outputs.py
```

CI (`.github/workflows/ci.yml`) runs both, plus dry runs over stubbed inputs.
Please make sure they pass before opening a PR.

## Commit messages & releases

This repo uses **[Conventional Commits](https://www.conventionalcommits.org)**
and [release-please](https://github.com/googleapis/release-please) to automate
versioning, the changelog, and releases. Prefix your commits:

| Prefix | Effect |
|---|---|
| `feat: …` | new feature → minor version bump, listed under *Added* |
| `fix: …` | bug fix → patch bump, under *Fixed* |
| `feat!: …` or a `BREAKING CHANGE:` footer | major bump |
| `docs:` / `refactor:` / `test:` / `chore:` | no release on their own |

On merge to `main`, release-please opens/updates a "release PR"; merging that PR
tags the version, publishes a GitHub Release, and updates `CHANGELOG.md` and
`CITATION.cff` (and, once the Zenodo integration is enabled, mints a DOI). You
do **not** edit the changelog or version numbers by hand.

## Pull requests

1. Branch from `main`, make your change, and add/adjust tests where relevant.
2. Run `pytest tests/ -q` and `snakemake -s workflow/Snakefile -n`.
3. Open a PR with a clear, Conventional-Commit-style title and description.
