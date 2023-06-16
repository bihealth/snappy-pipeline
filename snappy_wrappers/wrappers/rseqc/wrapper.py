# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for infer_experiment: Snakemake wrapper.py
"""

import os

from snakemake import shell

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

shell.executable("/bin/bash")

current_step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][current_step]["strandedness"]
out_link_dir = (
    os.path.dirname(snakemake.output.output) if "output" in snakemake.output.keys() else ""
)
log_link_dir = os.path.dirname(snakemake.output.log) if "output" in snakemake.output.keys() else ""

shell(
    r"""
set -euo pipefail
set -x

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}


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

# Run rseqc to infer strandedness
infer_experiment.py \
    -r "{config[path_exon_bed]}" \
    -i "{snakemake.input.bam}" \
    > "{snakemake.output.tsv}"


# Set strandedness based on input parameter or inferred value
strand={config[strand]}
if [ ${{strand}} -eq -1 ]
then
    pattern="^This is (Pair|Single)End Data$"
    endedness=$(grep -E "$pattern" {snakemake.output.tsv} | sed -E "s/$pattern/\1/")

    pattern="^Fraction of reads failed to determine: *([+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([EeDd][+-]?[0-9]+)?)$"
    failed=$(grep -E "$pattern" {snakemake.output.tsv} | sed -E "s/$pattern/\1/")

    if [ "$endedness" == "Single" ]
    then
        pattern="^Fraction of reads explained by \"\\+\\+,\\-\\-\": *([+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([EeDd][+-]?[0-9]+)?)$"
        forward=$(grep -E "$pattern" {snakemake.output.tsv} | sed -E "s/$pattern/\1/")
        pattern="^Fraction of reads explained by \"\\+\\-,\\-\\+\": *([+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([EeDd][+-]?[0-9]+)?)$"
        reverse=$(grep -E "$pattern" {snakemake.output.tsv} | sed -E "s/$pattern/\1/")
    else
        pattern="^Fraction of reads explained by \"1\\+\\+,1\\-\\-,2\\+\\-,2\\-\\+\": *([+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([EeDd][+-]?[0-9]+)?)$"
        forward=$(grep -E "$pattern" {snakemake.output.tsv} | sed -E "s/$pattern/\1/")
        pattern="^Fraction of reads explained by \"1\\+\\-,1\\-\\+,2\\+\\+,2\\-\\-\": *([+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([EeDd][+-]?[0-9]+)?)$"
        reverse=$(grep -E "$pattern" {snakemake.output.tsv} | sed -E "s/$pattern/\1/")
    fi

    reverse=$(echo $reverse | tr '[Dde]' 'E')

    strand=0
    if [ $(echo "$forward > {config[threshold]}" | bc -l) -gt 0 ]
    then
        strand=1
    fi
    if [ $(echo "$reverse > {config[threshold]}" | bc -l) -gt 0 ]
    then
        strand=2
    fi
fi

# Write strandedness
echo ${{strand}} > {snakemake.output.decision}

pushd $(dirname {snakemake.output.decision})
md5sum $(basename {snakemake.output.decision}) > $(basename {snakemake.output.decision}).md5
md5sum $(basename {snakemake.output.tsv}) > $(basename {snakemake.output.tsv}).md5
popd

if [[ -n "{out_link_dir}" ]];
then
    ln -sr {snakemake.output.decision} {out_link_dir}/.
    ln -sr {snakemake.output.decision}.md5 {out_link_dir}/.
fi
if [[ -n "{log_link_dir}" ]];
then
    ln -sr {snakemake.log.log} {log_link_dir}/.
    ln -sr {snakemake.log.log}.md5 {log_link_dir}/.
    ln -sr {snakemake.log.conda_list} {log_link_dir}/.
    ln -sr {snakemake.log.conda_list}.md5 {log_link_dir}/.
    ln -sr {snakemake.log.conda_info} {log_link_dir}/.
    ln -sr {snakemake.log.conda_info}.md5 {log_link_dir}/.
fi
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
