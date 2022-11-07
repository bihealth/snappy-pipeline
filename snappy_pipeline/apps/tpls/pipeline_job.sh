#!/bin/bash

# SNAPPY best practice pipeline_job.sh
#
# Version: 3
# Date: 2017-02-02

# The medium project/queue is a sensible default.
#SBATCH --partition %(partition)s
# Set a required running time for the master job.
#SBATCH --time 3-00
# Reserve some resources
#SBATCH --mem=6G
# Keep current environment variables
#SBATCH --export=all
# Send a mail upon job completion and error
%(line_m)s
%(line_M)s
# Logs should be written into "slurm_log" sub directory.
#SBATCH --output slurm_log/%%x-%%J.log
# Use more descriptive name in Slurm.
#SBATCH --job-name %(step_name)s

# Enable the official bash strict mode (fail early, fail often)
set -euo pipefail

# Fix the umask.
umask ug=rwx,o=

# Configuration variables ---------------------------------------------------

# Maximal number of jobs to execute at the same time
MAX_JOBS=500
# Maximal number of jobs per second
MAX_JOBS_PER_SECOND=10
# Number of times to restart jobs
RESTART_TIMES=0

# Check preconditions -------------------------------------------------------

# Ensure slurm_log is a directory
test -d slurm_log || { >&2 echo "${PWD}/slurm_log does not exist"; exit 1; }

# Enforce existence of TMPDIR -----------------------------------------------

export TMPDIR=${HOME}/scratch/tmp
mkdir -p ${TMPDIR}

# Create one log directory per Snakemake run --------------------------------

test -z "${SLURM_JOB_ID-}" && SLURM_JOB_ID=$(date +%%Y-%%m-%%d_%%H-%%M)
LOGDIR=slurm_log/${SLURM_JOB_ID}
mkdir -p ${LOGDIR}
export SBATCH_DEFAULTS=" --output=${LOGDIR}/%%x-%%j.log"

# Activate appropriate Miniconda3 installation ------------------------------

# 1. If CONDA_PATH is set, use this.
# 2. Look into parent directories for miniconda3 (owned by current user)
# 3. Look whether there is a conda in $PATH and use it.
# 4. Look for ~/miniconda3 and use it
# 5. If all fails, bail out.

conda-in-parent()
{
    current=$PWD
    while [[ -n "$current" ]] && [[ "$current" != "/" ]]; do
        if [[ -e "$current/miniconda3.$USER" ]] && \
                [[ $(stat -c %%u $current/miniconda3.$USER) == $UID ]]; then
            echo "$current/miniconda3.$USER"
            return 0
        fi
        if [[ -e "$current/miniconda3" ]] && \
                [[ $(stat -c %%u $current/miniconda3) == $UID ]]; then
            echo "$current/miniconda3"
            return 0
        fi
        current=$(dirname $current)
    done

    return 1
}

if [[ -n "${CONDA_PATH-}" ]] || CONDA_PATH=$(conda-in-parent); then
    :
elif which conda >/dev/null; then
    CONDA_PATH=$(dirname $(dirname $(which conda)))
elif [[ -e $HOME/miniconda3 ]]; then
    CONDA_PATH=$HOME/miniconda3
elif [[ -e $HOME/work/miniconda3 ]]; then
    CONDA_PATH=$HOME/work/miniconda3
else
    >&2 echo "Could not determine a suitable CONDA_PATH."
    exit 1
fi

>&2 echo "Using conda installation in $CONDA_PATH"
>&2 echo "+ conda activate %(conda)s"
set +euo pipefail
conda deactivate &>/dev/null || true  # disable any existing
source $CONDA_PATH/etc/profile.d/conda.sh
conda activate %(conda)s # enable found
set -euo pipefail

# Activate bash cmd printing, debug info ------------------------------------

set -x
>&2 hostname
>&2 date

# Kick off Snakemake --------------------------------------------------------

# Using the medium project/queue is a sensible default.
snappy-snake --printshellcmds \
    --snappy-pipeline-use-profile "cubi-v1" \
    --snappy-pipeline-jobs $MAX_JOBS \
    --restart-times ${RESTART_TIMES} \
    --default-partition="medium" \
    --rerun-incomplete \
    -- \
    $*

# Print date after finishing, for good measure ------------------------------

>&2 date
>&2 echo "All done. Have a nice day."
