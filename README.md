[![CI](https://github.com/bihealth/snappy-pipeline/workflows/CI/badge.svg?branch=main)](https://github.com/bihealth/snappy-pipeline/actions/workflows/main.yml)
[![Coverage Status](https://coveralls.io/repos/github/bihealth/snappy-pipeline/badge.svg?branch=master)](https://coveralls.io/github/bihealth/snappy-pipeline?branch=master)
[![Documentation Status](https://readthedocs.org/projects/snappy-pipeline/badge/?version=latest)](https://snappy-pipeline.readthedocs.io/en/latest/?badge=latest)

# SNAPPY - SNAPPY Nucleic Acid Processing Pipeline

A Snakemake-based pipeline for processing NGS data — read mapping, variant calling, annotation, filtration, and more.

## Quick Install

**Prerequisites**: Install [pixi](https://pixi.sh) (see https://pixi.sh/latest/#installation).

```bash
git clone git@github.com:bihealth/snappy-pipeline.git
cd snappy-pipeline
pixi install
```

After installation the `snappy` command is available via `pixi run snappy ...` or by activating the pixi environment with `eval "$(pixi shell-hook)"`.

## Quick Start

```bash
snappy init --directory my_project
cd my_project
# edit config.yaml, add samplesheet.tsv
snappy run
```

See [user quickstart](docs/quickstart.rst) for details.

See [developer setup](docs/installation.rst) for testing, linting, and building documentation.

## Development Notes

Here, you can find the required layout for post-PR commit messages:

- https://github.com/amannn/action-semantic-pull-request#configuration
