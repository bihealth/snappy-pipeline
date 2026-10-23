(installation)=

# Installation

## Prerequisites

Install [pixi](https://pixi.sh) (see <https://pixi.sh/latest/#installation>).
pixi manages all dependencies of the pipeline, the Python packages and the conda packages (BWA, STAR, samtools, ...).

For cluster execution, a Snakemake profile of the system or the user selects the executor, for example SLURM; see {ref}`usage`.

PureCN and pVACtools run from container images, which need [apptainer](https://apptainer.org) on the machines that run their jobs.
On Ubuntu 24.04 and newer, apptainer also needs an AppArmor profile that allows unprivileged user namespaces.

## User installation

To run the pipeline, clone the repository and let pixi create the environment:

```shell
git clone git@github.com:bihealth/snappy-pipeline.git
cd snappy-pipeline
pixi install
```

The `snappy` command is then available as `pixi run snappy <subcommand> ...`.
Alternatively, activate the environment with `eval "$(pixi shell-hook)"`.

`pixi install --frozen` installs the exact versions in `pixi.lock`.
`pixi update` updates the lock file from `pyproject.toml`.

The tools of each step are installed in conda environments that Snakemake creates on the first run.

## Developer installation

Use the `dev` environment, which adds the test, lint and documentation tools:

```shell
pixi install -e dev
pixi run -e dev <command>   # one-off commands
pixi shell -e dev           # activate the environment
```

### Tests

```shell
pixi run -e dev test
```

### Style checks

```shell
pixi run -e dev lint       # ruff check, ruff format --check and snakefmt --check
pixi run -e dev fmt        # format with ruff
pixi run -e dev snakefmt   # format the Snakefiles
pixi run -e dev srcfmt     # both formatters
```

### Documentation

```shell
pixi run -e docs docs
```

The HTML output is in `docs/_build/html/`.
Format Markdown pages with `uvx rumdl fmt --flavor myst <file>`.

### Developer documentation

Read the "Pipeline Developers Docs" as well, starting with {ref}`dev_intro`.

## Configuring GATK3

Some wrappers use GATK 3.
GATK 3 is not free software and cannot be redistributed.
Download the tarball [from the Broad archive](https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2).

Create the conda environments first.
Then register the tarball in each environment that has GATK 3:

```shell
cd /path/to/project
snappy run -- --deploy-software-envs
for meta in .snakemake/software/deployments/conda/*/conda-meta/gatk-3.*.json; do
    environ=$(dirname $(dirname $meta))
    conda activate $environ
    gatk3-register /path/to/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
    conda deactivate
done
```
