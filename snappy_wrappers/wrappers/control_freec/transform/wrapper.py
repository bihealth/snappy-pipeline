# -*- coding: utf-8 -*-
"""Wrapper for merging multiple tables in R on shared columns"""

from snakemake import shell
import os

shell.executable("/bin/bash")

config = snakemake.config["step_config"]["somatic_wgs_cnv_calling"]["control_freec"]["convert"]
rscript = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "snappy-convert-control_freec.R"
)

shell(
    r"""
set -x

# Write out information about conda installation --------------------------------------------------

conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}

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

# Setup auto-cleaned TMPDIR -----------------------------------------------------------------------

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

R --vanilla -e "source(\"{rscript}\") ; library(magrittr) ; \
    control_freec_write_files( \
    sample_name = \"{snakemake.wildcards.cancer_library}\", \
    ratios_fn = \"{snakemake.input.ratio}\", \
    log2_fn = \"{snakemake.output.log2}\", \
    call_fn = \"{snakemake.output.call}\", \
    segments_fn = \"{snakemake.output.segments}\", \
    cns_fn = \"{snakemake.output.cns}\", \
    cnr_fn = \"{snakemake.output.cnr}\", \
    org_obj={config[org_obj]}, \
    tx_obj={config[tx_obj]}, \
    bs_obj={config[bs_obj]})"

for f in {snakemake.output.log2} {snakemake.output.call} {snakemake.output.segments} \
    {snakemake.output.cns} {snakemake.output.cnr}; do
    md5sum $f >$f.md5
done

"""
)
