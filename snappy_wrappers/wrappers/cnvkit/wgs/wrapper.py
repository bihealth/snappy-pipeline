# -*- coding: utf-8 -*-
"""Wrapper vor cnvkit.py batch --method wgs
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bihealth.de"

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Also pipe stderr to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

# Run cvnkit.py

cnvkit.py batch \
    {snakemake.input.tumor_bam} \
    -n {snakemake.input.normal_bam} \
    --method wgs \
    --annotate {snakemake.config[step_config][somatic_wgs_cnv_calling][cnvkit][path_annotate_refflat]} \
    -f {snakemake.config[static_data_config][reference][path]} \
    --output-dir $(dirname {snakemake.output.segment})

pushd $(dirname {snakemake.output.segment})
file=$(basename {snakemake.output.segment})
mapper=$(echo $file | cut -f1 -d".")
for i in $mapper.*; do ln -sr -T $i $(echo $i | sed "s/$mapper/$mapper\.cnvkit/") ; done
popd

unset DISPLAY

cnvkit.py scatter \
    {snakemake.output.bins} \
    -o {snakemake.output.scatter} \
    -s {snakemake.output.segment}

pushd $(dirname {snakemake.output.segment})
md5sum $(basename {snakemake.output.segment}) > $(basename {snakemake.output.segment}).md5
md5sum $(basename {snakemake.output.bins}) > $(basename {snakemake.output.bins}).md5
md5sum $(basename {snakemake.output.scatter}) > $(basename {snakemake.output.scatter}).md5
popd
"""
)
