#!env bash
set -o pipefail
conda init
conda activate "$CONDA_PREFIX"

if [[ -n "${GITHUB_WORKSPACE}" ]]; then
  cd "${GITHUB_WORKSPACE}" || exit
else
  if [[ -d /github/workspace ]]; then
    cd /github/workspace/ || exit
  else
    cd snappy-pipeline || exit
  fi
fi
pip install -e .
