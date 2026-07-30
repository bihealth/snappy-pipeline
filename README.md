[![CI](https://github.com/bihealth/snappy-pipeline/workflows/CI/badge.svg?branch=main)](https://github.com/bihealth/snappy-pipeline/actions/workflows/main.yml)
[![Coverage Status](https://coveralls.io/repos/github/bihealth/snappy-pipeline/badge.svg?branch=master)](https://coveralls.io/github/bihealth/snappy-pipeline?branch=master)
[![Documentation Status](https://readthedocs.org/projects/snappy-pipeline/badge/?version=latest)](https://snappy-pipeline.readthedocs.io/en/latest/?badge=latest)

# SNAPPY - SNAPPY Nucleic Acid Processing Pipeline

## Installation

Installation should be complete in 10 to 15 minutes.

**Prerequisites**: Install [pixi](https://pixi.sh) (see https://pixi.sh/latest/#installation).

**In a nutshell**:

```bash
git clone git@github.com:bihealth/snappy-pipeline.git
cd snappy-pipeline
pixi install
```

After installation the `snappy` command is available via `pixi run snappy ...` or by activating the pixi environment with `eval "$(pixi shell-hook)"`.

See [user quickstart](docs/quickstart.rst) if you just want to use the pipeline.

See [developer setup](docs/installation.rst) for testing, linting, and building documentation.

## Using GATK3

Some wrappers rely on GATK 3.
GATK v3 is not free software and cannot be redistributed.
If you are a member of CUBI, you can use the central GATK download.
Alternatively, you can download the tarball [from the Broad archive](https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2).

To use GATKv3 with the pipeline, run:

```bash
$ cd variant_calling
$ snappy run -- --conda-create-envs-only
```

Find which conda environments use GATK v3 and register the tarball as described in the Snakemake documentation.

## Development Notes

Here, you can find the required layout for post-PR commit messages:

- https://github.com/amannn/action-semantic-pull-request#configuration
