# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for arriba: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

shell.executable("/bin/bash")

# Input fastqs are passed through snakemake.params.
# snakemake.input is a .done file touched after linking files in.
reads_left = snakemake.params.args["input"]["reads_left"]
reads_right = (
    snakemake.params.args["input"]["reads_right"]
    if snakemake.params.args["input"]["reads_right"]
    else ""
)

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

trim_cmd=""
if [[ "{snakemake.config[step_config][somatic_gene_fusion_calling][arriba][trim_adapters]}" == "True" ]]; then
    trim_cmd="\"trimadap-mt -p {snakemake.config[step_config][somatic_gene_fusion_calling][arriba][num_threads_trimming]}\""
else
    trim_cmd="zcat"
fi

STAR \
    --runThreadN {snakemake.config[step_config][somatic_gene_fusion_calling][arriba][num_threads]} \
    --genomeDir {snakemake.config[step_config][somatic_gene_fusion_calling][arriba][path_index]} --genomeLoad NoSharedMemory \
    --readFilesIn ${{left_files}} ${{right_files}} --readFilesCommand ${{trim_cmd}} \
    --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
    --outFileNamePrefix $TMPDIR/ \
    {snakemake.config[step_config][somatic_gene_fusion_calling][arriba][star_parameters]} \
| arriba \
    -x /dev/stdin \
    -o $TMPDIR/fusions.tsv -O $TMPDIR/fusions.discarded.tsv \
    -a {snakemake.config[static_data_config][reference][path]} \
    -g {snakemake.config[step_config][somatic_gene_fusion_calling][arriba][features]} \
    $(if [[ -n "{snakemake.config[step_config][somatic_gene_fusion_calling][arriba][blacklist]}" ]]; then \
        echo -b {snakemake.config[step_config][somatic_gene_fusion_calling][arriba][blacklist]}
    fi) \
    $(if [[ -n "{snakemake.config[step_config][somatic_gene_fusion_calling][arriba][known_fusions]}" ]]; then \
        echo -k {snakemake.config[step_config][somatic_gene_fusion_calling][arriba][known_fusions]}
    fi) \
    $(if [[ -n "{snakemake.config[step_config][somatic_gene_fusion_calling][arriba][tags]}" ]]; then \
        echo -t {snakemake.config[step_config][somatic_gene_fusion_calling][arriba][tags]}
    fi) \
    $(if [[ -n "{snakemake.config[step_config][somatic_gene_fusion_calling][arriba][structural_variants]}" ]]; then \
        echo -d {snakemake.config[step_config][somatic_gene_fusion_calling][arriba][structural_variants]}
    fi) \
    $(if [[ -n "{snakemake.config[step_config][somatic_gene_fusion_calling][arriba][protein_domains]}" ]]; then \
        echo -p {snakemake.config[step_config][somatic_gene_fusion_calling][arriba][protein_domains]}
    fi)

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
