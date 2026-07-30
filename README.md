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
Earlier, we had an internal CUBI conda server but this limits use of the wrapper for the general public.
Now, the using pipeline steps must be activated as follows.

If you are a member of CUBI, you can use the central GATK download.
Alternatively, you can download the tarball [from the Broad archive](https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2).

```bash
$ ls -lh /fast/groups/cubi/work/projects/biotools/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
-rw-rw---- 1 holtgrem_c hpc-ag-cubi 14M Dec 19  2019 /fast/groups/cubi/work/projects/biotools/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
```

First, go to the pipeline directory where you want to run:

```bash
$ cd variant_calling
```

Explicitely create any missing conda environment

```bash
$ snappy-snake --conda-create-envs-only
[...]
12-27 17:18 snakemake.logging WARNING  Downloading and installing remote packages.
[...]
```

Find out which conda environments use GATK v3

```bash
$ grep 'gatk.*3' .snakemake/conda/*.yaml
.snakemake/conda/d76b719b718c942f8e49e55059e956a6.yaml:  - gatk =3
```

Activate each conda environment and register

```bash
$ for yaml in $(grep -l 'gatk.*3' .snakemake/conda/*.yaml); do
        environ=${yaml%.yaml};
        conda activate $environ
        gatk3-register /fast/groups/cubi/work/projects/biotools/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
        conda deactivate
    done
Moving GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 to /home/holtgrem_c/miniconda3/envs/gatk3/opt/gatk-3.8
```

You are now ready to run GATK v3 from this environment.

## Development Notes

Here, you can find the required layout for post-PR commit messages:

- https://github.com/amannn/action-semantic-pull-request#configuration
