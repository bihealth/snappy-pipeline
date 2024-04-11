#!env bash
set -o pipefail
conda init
conda activate "$CONDA_PREFIX"
git clone https://github.com/bihealth/snappy-pipeline.git
cd snappy-pipeline || exit
git checkout e260347ebed0b54cc8c84e48c932581ac59cd667
pip install -e .
