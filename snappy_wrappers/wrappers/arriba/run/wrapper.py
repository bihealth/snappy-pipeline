# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for arriba: Snakemake wrapper.py"""

from typing import TYPE_CHECKING

from snakemake.shell import shell

if TYPE_CHECKING:
    from snakemake.script import snakemake

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

shell.executable("/bin/bash")

# Input fastqs are passed through snakemake.params.
# snakemake.input is a .done file touched after linking files in.
args = getattr(snakemake.params, "args", {})
reads_left = args["input"]["reads_left"]
reads_right = (
    args["input"]["reads_right"]
    if args["input"]["reads_right"]
    else ""
)

trim_adapters = args["trim_adapters"]
num_threads_trimming = args["num_threads_trimming"]
trim_cmd = "trimadap-mt -p {num_threads_trimming}" if trim_adapters else "zcat"

num_threads = args["num_threads"]
arriba_index = args["path_index"]
star_parameters = args["star_parameters"]

reference_path = args["reference_path"]
features_path = args["features_path"]

blacklist = args["blacklist"]
blacklist_param = f"-b {blacklist}" if blacklist else ""

known_fusions = args["known_fusions"]
known_fusions_param = f"-k {known_fusions}" if known_fusions else ""

tags = args["tags"]
tags_param = f"-t {tags}" if tags else ""

structural_variants = args["structural_variants"]
structural_variants_param = f"-d {structural_variants}" if structural_variants else ""

protein_domains = args["protein_domains"]
protein_domains_param = f"-p {protein_domains}" if protein_domains else ""

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
mkdir -p $TMPDIR

# Define left and right reads as Bash arrays
declare -a reads_left=({reads_left})
# declared but never used
declare -a reads_right=({reads_right})

left_files=$(IFS="," ; echo "${{reads_left[*]}}")

right_files=""
if [[ "{reads_right}" != "" ]]; then
    right_files=$(IFS="," ; echo "${{reads_right[*]}}")
fi

STAR \
    --runThreadN {num_threads} \
    --genomeDir {arriba_index} --genomeLoad NoSharedMemory \
    --readFilesIn ${{left_files}} ${{right_files}} --readFilesCommand {trim_cmd} \
    --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
    --outFileNamePrefix $TMPDIR/ \
    {star_parameters} \
| arriba \
    -x /dev/stdin \
    -o $TMPDIR/fusions.tsv -O $TMPDIR/fusions.discarded.tsv \
    -a {reference_path} \
    -g {features_path} \
    {blacklist_param} \
    {known_fusions_param} \
    {tags_param} \
    {structural_variants_param} \
    {protein_domains_param}

cp $TMPDIR/fusions.tsv {snakemake.output.fusions}
pushd $(dirname {snakemake.output.fusions})
md5sum $(basename {snakemake.output.fusions}) > $(basename {snakemake.output.fusions}).md5
popd

gzip -c $TMPDIR/fusions.discarded.tsv > {snakemake.output.discarded}
pushd $(dirname {snakemake.output.discarded})
md5sum $(basename {snakemake.output.discarded}) > $(basename {snakemake.output.discarded}).md5
popd

log_dir=$(dirname {snakemake.log.log})
star_logs=$(echo "Log.out Log.std.out Log.final.out SJ.out.tab" | tr ' ' '\n')
for star_log in $star_logs
do
    cp $TMPDIR/$star_log $log_dir/$star_log
done

pushd $log_dir
for star_log in $star_logs
do
    md5sum $star_log > $star_log.md5
done
popd

touch {snakemake.output.done}
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
