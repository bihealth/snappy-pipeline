# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for computing PureCN coverage
"""

from snakemake import shell

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]["purecn"]

shell.executable("/bin/bash")

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

export R_LIBS_USER=$(dirname {snakemake.input.packages})

# Find PureCN scripts in extdata
pkg_folders=$(R --quiet --vanilla -e 'cat(.libPaths(), sep="\n")' | grep -v '^>')
PURECN=$(for folder in $pkg_folders ; do ls -1 $folder | grep -E '^PureCN$' | sed -e "s#^#$folder/#" ; done)
PURECN="$PURECN/extdata"

# Create coverage
Rscript $PURECN/Coverage.R --force \
    --seed {config[seed]} \
    --out-dir $(dirname {snakemake.output.coverage}) \
    --bam {snakemake.input.bam} \
    --intervals {snakemake.input.intervals}
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
