#!/bin/bash

# SNAPPY best practice pipeline_job.sh
#
# Submit a snakemake job through the snappy pipeline wrapper "snappy-snake"
# This script is written for snappy pipeline >=0.6, snakemake >= 9.13 & SLURM
# A job can be submitted either using salloc or sbatch (HPC rules don't allow srun on login nodes).
# * salloc:
#     +: recommended in the snakemake slurm plugin documentation
#        (https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html#how-to-run-snakemake-in-an-interactive-slurm-job)
#     -: replaces current process, so better to put it in the background
#        SBATCH directives are ignored (no jobname, partition, time or mem).
#        They must be explicitely added by the user on the command line:
#        salloc --job-name snappy_step --partition medium --time 3-00 --mem 6G pipeline_job.sh &
# * sbatch:
#     +: SBATCH directives are followed
#     -: not recommended in the official snakemake slurm plugin documentation
#
#
# Version: 4
# Date: 2025-10-31

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
# Use more descriptive name in Slurm.
#SBATCH --job-name %(step_name)s

# Enable the official bash strict mode (fail early, fail often)
set -euo pipefail

# Fix the umask.
umask ug=rwx,o=

# Create jobname from pipeline step name ------------------------------------
JOBNAME="%(step_name)s"

# Create one log directory per Snakemake run --------------------------------

test -z "${SLURM_JOB_ID-}" && SLURM_JOB_ID=$(date +%%Y-%%m-%%d_%%H-%%M)
LOGDIR=slurm_log/${SLURM_JOB_ID}
mkdir -p ${LOGDIR}
LOGFILE=slurm_log/snakemake.${JOBNAME}-${SLURM_JOB_ID}.log

# Find SLURM account --------------------------------------------------------

ME=$(whoami)
ACCOUNT=$(sacctmgr show associations -P user=$ME | tail -n +2 | head -n 1 | tr '|' '\t' | cut -f 2)

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
        if [[ -e "$current/miniforge3.$USER" ]] && \
                [[ $(stat -c %%u $current/miniforge3.$USER) == $UID ]]; then
            echo "$current/miniforge3.$USER"
            return 0
        fi
        if [[ -e "$current/miniforge3" ]] && \
                [[ $(stat -c %%u $current/miniforge3) == $UID ]]; then
            echo "$current/miniforge3"
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
elif [[ -e $HOME/miniforge3 ]]; then
    CONDA_PATH=$HOME/miniforge3
elif [[ -e $HOME/work/miniforge3 ]]; then
    CONDA_PATH=$HOME/work/miniforge3
else
    echo "Could not determine a suitable CONDA_PATH."
    exit 1
fi

echo "Using conda installation in $CONDA_PATH"
echo "+ conda activate snappy_dev"
set +euo pipefail
conda deactivate &>/dev/null || true  # disable any existing
source $CONDA_PATH/etc/profile.d/conda.sh
conda activate snappy_dev # enable found
set -euo pipefail

# Activate bash cmd printing, debug info ------------------------------------

hostname > $LOGFILE
date    >> $LOGFILE

# Kick off Snakemake --------------------------------------------------------

# The slurm account value cannot be defined on the command line, so the job is dispached without account
snappy-snake \
    --profile-snappy-pipeline \
    --printshellcmds --jobname "snakemake.$JOBNAME.{jobid}" --slurm-logdir $LOGDIR --slurm-keep-successful-logs --slurm-no-account \
    $* \
1>> $LOGFILE 2>&1

# Print date after finishing, for good measure ------------------------------

date >> $LOGFILE
echo "All done. Have a nice day." >> $LOGFILE
