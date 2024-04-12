#!env bash
set -o pipefail
conda init
conda activate "$CONDA_PREFIX"
if [[ -z "${GITHUB_WORKSPACE}" ]]; then
  cd "${GITHUB_WORKSPACE}" || exit
else
  cd snappy-pipeline || exit
fi
pip install -e .
