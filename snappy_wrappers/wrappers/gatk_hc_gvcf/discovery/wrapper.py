# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for GATK HC GVCF Calling: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

import sys

print("PARAMS", snakemake.params, file=sys.stderr)

arg_annotations = " ".join(
    [
        "--annotation {}".format(anno)
        for anno in snakemake.config["step_config"]["variant_calling"]["gatk_hc"]["annotations"]
    ]
)
arg_intervals = " ".join(
    ["--intervals {}".format(interval) for interval in snakemake.params["args"]["intervals"]]
)

shell.executable("/bin/bash")

shell(
    r"""
set -x

# TODO: add through shell.prefix
export TMPDIR=$HOME/scratch/tmp

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/tmp.d

# TODO: Fix link problems of tabix. What the link?
export JAVA_HOME=$(dirname $(which gatk_nonfree))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

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

# Prepare arguments
if [[ "{snakemake.config[step_config][variant_calling][gatk_hc_gvcf][allow_seq_dict_incompatibility]}" == "True" ]]; then
    arg_seq_dict="-U ALLOW_SEQ_DICT_INCOMPATIBILITY"
else
    arg_seq_dict=
fi

# Call GATK through the bioconda wrapper
gatk_nonfree -Xmx6g -Djava.io.tmpdir=$TMPDIR \
    --analysis_type HaplotypeCaller \
    $arg_seq_dict \
    -nct 1 \
    --emitRefConfidence GVCF \
    --variant_index_type LINEAR \
    --variant_index_parameter 128000 \
    --out {snakemake.output.vcf} \
    --reference_sequence {snakemake.config[static_data_config][reference][path]} \
    --sample_ploidy 2 \
    --pair_hmm_implementation VECTOR_LOGLESS_CACHING \
    -minPruning 3 \
    --dbsnp {snakemake.config[static_data_config][dbsnp][path]} \
    {arg_intervals} \
    --input_file $(echo {snakemake.input} | tr ' ' '\n' | grep '\.bam$')

# compute tabix index of the resulting VCF file
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf}) && \
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5 && \
    md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi}).md5
"""
)
