# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for Epitoper: Snakemake wrapper.py"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

info_fields = " ".join(snakemake.params.args["vcf_info_fields"])

shell(
    r"""
set -x

# Also pipe stderr to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec 2> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/tmp.d

# Call epitoper
_JAVA_OPTIONS="-Xmx5G" \
epitoper \
    --blastp-blastdb {snakemake.config[step_config][epitope_prediction][epitoper][path_blastdb]} \
    --jv-db {snakemake.config[step_config][epitope_prediction][epitoper][path_jannovar_ser]} \
    --binding-prediction-method {snakemake.params.args[binding_prediction_method]} \
    --blastp-num-threads {snakemake.params.args[blastp_num_threads]} \
    --hla-type {snakemake.params.args[hla_type]} \
    --input-vcf {snakemake.input.vcf} \
    $(if [[ -n "{info_fields}" ]]; then
        for pair in {info_fields}; do \
            echo --vcf-info-field $pair; \
        done; \
    fi) \
    $(if [[ -n "{snakemake.params.args[tumor_dna_sample]}" ]]; then \
        echo --tumor-dna-sample {snakemake.params.args[tumor_dna_sample]}; \
    fi) \
    $(if [[ -n "{snakemake.params.args[tumor_rna_sample]}" ]]; then \
        echo --tumor-rna-sample {snakemake.params.args[tumor_rna_sample]}; \
    fi) \
    $(if [[ "{snakemake.params.args[verbose]}" == "True" ]]; then \
        echo -vv; \
    fi) \
| gzip -c \
> {snakemake.output.txt_gz}

pushd $(dirname {snakemake.output.txt_gz}) && \
    md5sum $(basename {snakemake.output.txt_gz}) >$(basename {snakemake.output.txt_gz}).md5
"""
)
