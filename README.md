[![CI](https://github.com/bihealth/snappy-pipeline/workflows/CI/badge.svg?branch=main)](https://github.com/bihealth/snappy-pipeline/actions/workflows/main.yml)
[![Coverage Status](https://coveralls.io/repos/github/bihealth/snappy-pipeline/badge.svg?branch=master)](https://coveralls.io/github/bihealth/snappy-pipeline?branch=master)
[![Documentation Status](https://readthedocs.org/projects/snappy-pipeline/badge/?version=latest)](https://snappy-pipeline.readthedocs.io/en/latest/?badge=latest)

# SNAPPY - SNAPPY Nucleic Acid Processing Pipeline

## Installation

Installation should be complete in 10 to 15 minutes.

**In a nutshell**:

```bash
# Download & preparation
git clone git@github.com:bihealth/snappy-pipeline.git
cd snappy-pipeline

# If you want to select a given branch, uncomment the following:
# git checkout <branch_name>

# WARNING- make sure that you are in your conda base environment

# Create conda environment "snappy_env" with the minimal requirements:
mamba env create --file environment.yml
conda activate snappy_env

# Add testing & development requirements:
pip install -r requirements/test.txt
pip install -r requirements/dev.txt

# Optionally add "pytest-pdb" missing from anaconda
pip install pytest-pdb

# Install snappy in snappy_env environment
pip install -e .
```

**Note:** To create the environment under another name, replace the commands for the environment creation & activation of the correct environment by:

```
mamba env create --file environment.yml --name <other_environment_name>
conda activate <other_environment_name>
```


See [user installation](docs/quickstart.rst) if you just want to use the pipeline.

See [developer installation)[docs/installation.rst) for getting started with working on the pipeline code and also building the documentation.

