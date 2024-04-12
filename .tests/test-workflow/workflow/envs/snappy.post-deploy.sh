#!env bash
set -o pipefail
conda init
conda activate "$CONDA_PREFIX"
cd snappy-pipeline || exit
pip install -e .
