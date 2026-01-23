# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for STAR: Snakemake wrapper.py"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

args = getattr(snakemake.params, "args", {})

# Input fastqs are passed through snakemake.params.
# snakemake.input is a .done file touched after linking files in.
reads_left = args["input"]["reads_left"]
reads_right = args["input"].get("reads_right", "")

# salmon flag for first reads changes for single-end data.
if reads_right:
    read_flag = "-1"
else:
    read_flag = "-r"

this_file = __file__

shell(
    r"""
set -x

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

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/tmp.d $TMPDIR/pre.d

# Define left and right reads as Bash arrays
declare -a reads_left=({reads_left})
# declared but never used
declare -a reads_right=({reads_right})

left_files=$(IFS=" " ; echo "${{reads_left[*]}}")
left_files_prefixed=" ${{left_files}}"

right_files_prefixed=""
if [[ "{reads_right}" != "" ]]; then
    right_files=$(IFS=" " ; echo "${{reads_right[*]}}")
    right_files_prefixed=" -2 ${{right_files}}"
fi

libraryType="A"
if [[ {args[strand]} -ge 0 ]]
then
    libraryType="I"
    if [[ {args[strand]} -gt 0 ]]
    then
        libraryType="${{libraryType}}S"
        if [[ {args[strand]} -eq 1 ]]
        then
            libraryType="${{libraryType}}F"
        else
            libraryType="${{libraryType}}R"
        fi
    else
        libraryType="${{libraryType}}U"
    fi
fi

salmon quant \
    -i $(dirname {input.indices}) \
    -l $libraryType \
    {read_flag} ${{left_files_prefixed}} ${{right_files_prefixed}} \
    -g {input.features} \
    -o $TMPDIR \
    -p {args[num_threads]} \
    --auxDir aux \
    {args[salmon_params]}

# Copy over the output files
cp $TMPDIR/quant.sf {snakemake.output.transcript_sf}
md5sum {snakemake.output.transcript_sf} > {snakemake.output.transcript_sf_md5}
if [[ "${{t2g_cmd}}" != "" ]]
then
    cp $TMPDIR/quant.genes.sf {snakemake.output.gene_sf}
    md5sum {snakemake.output.gene_sf} > {snakemake.output.gene_sf_md5}
fi

# Copy log files
log=$(dirname {snakemake.log.log})
cp $TMPDIR/cmd_info.json $log/cmd_info.json
cp $TMPDIR/lib_format_counts.json $log/lib_format_counts.json
cp $TMPDIR/logs/salmon_quant.log $log/salmon_quant.log
md5sum $log/cmd_info.json > $log/cmd_info.json.md5
md5sum $log/lib_format_counts.json > $log/lib_format_counts.json.md5
md5sum $log/salmon_quant.log > $log/salmon_quant.log.md5

# Copy extra files
aux=$(dirname {snakemake.output.transcript_sf})/aux
mkdir -p $aux
cp -R $TMPDIR/aux/* $aux/.

# Logging: Save a copy this wrapper (with the pickle details in the header)
cp {this_file} $(dirname $log)/wrapper_star.py
# Logging: Save a permanent copy of the environment file used
# Commented out: this_file is the temp copy, not the original
# cp $(dirname {this_file})/environment.yaml $(dirname $log)/environment_wrapper_star.yaml
"""
)
