# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for computing PureCN coverage
"""

import os

from snakemake import shell

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]["purecn"]

shell.executable("/bin/bash")

# Prepare files and directories that must be accessible by the container
bound_files = {
    "bam": os.path.realpath(snakemake.input.bam),
}

keys = list(bound_files.keys())
bindings = []
for i in range(len(keys)):
    k = keys[i]
    if bound_files[k]:
        # Binding directory to /bindings/d<i>
        bindings.append(" -B {}:/bindings/d{}:ro".format(os.path.dirname(bound_files[k]), i))
        bound_files[k] = "/bindings/d{}/{}".format(i, os.path.basename(bound_files[k]))
bindings = " ".join(bindings)

shell(
    r"""
set -x

# Also pipe everything to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec &> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

# Setup auto-cleaned tmpdir
export tmpdir=$(mktemp -d)
trap "rm -rf $tmpdir" EXIT

# Create coverage
cmd="/usr/local/bin/Rscript /opt/PureCN/Coverage.R --force \
    --seed {config[seed]} \
    --out-dir $(dirname {snakemake.output.coverage}) \
    --bam {bound_files[bam]} \
    --intervals {snakemake.input.intervals}
"
mkdir -p $(dirname {snakemake.output.coverage})
apptainer exec --home $PWD {bindings} {snakemake.input.container} $cmd

pushd $(dirname {snakemake.output.coverage})
md5sum $(basename {snakemake.output.coverage}) > $(basename {snakemake.output.coverage}).md5
popd
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
