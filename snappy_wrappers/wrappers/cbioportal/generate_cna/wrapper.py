# -*- coding: utf-8 -*-
"""Wrapper for merging multiple tables in R on shared columns"""

import os

from snakemake import shell

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

r_script = os.path.abspath(os.path.join(os.path.dirname(__file__), "script.R"))
helper_functions = os.path.join(os.path.dirname(r_script), "..", "helper_functions.R")

args = getattr(snakemake.params, "args", {})

shell(
    r"""
set -x

# Write files for reproducibility -----------------------------------------------------------------

conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
pushd $(dirname {snakemake.log.conda_list}) ; md5sum $(basename {snakemake.log.conda_list}) > $(basename {snakemake.log.conda_list_md5}) ; popd
pushd $(dirname {snakemake.log.conda_info}) ; md5sum $(basename {snakemake.log.conda_info}) > $(basename {snakemake.log.conda_info_md5}) ; popd

# Also pipe stderr to log file --------------------------------------------------------------------

if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Create auto-cleaned temporary directory ---------------------------------------------------------

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Run the R script --------------------------------------------------------------------------------

R --vanilla --slave << __EOF
source("{helper_functions}")
source("{r_script}")
write.table(
    cns_to_cna("{snakemake.input.DNAcopy}", "{snakemake.input.features}", "{args[pipeline_id]}"),
    file="{snakemake.output}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
)
__EOF

pushd $(dirname {snakemake.output}) ; md5sum $(basename {snakemake.output}) > $(basename {snakemake.output}).md5 ; popd
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
sleep 1s  # try to wait for log file flush
pushd $(dirname {snakemake.log.log}) ; md5sum $(basename {snakemake.log.log}) > $(basename {snakemake.log.log_md5}) ; popd
"""
)
