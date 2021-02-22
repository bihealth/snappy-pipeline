# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for GATK GenotypeGVCF: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

arg_intervals = " ".join(
    ["--intervals {}".format(interval) for interval in snakemake.params["args"]["intervals"]]
)

shell.executable("/bin/bash")

shell(
    r"""
set -x

# TODO: add through shell.prefix
export TMPDIR=$HOME/scratch/tmp

# TODO: Fix link problems of tabix. What the link?
export JAVA_HOME=$(dirname $(which gatk_nonfree))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/tmp.d

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
#
# We break at multiples of 10'000 as we use this for window overlap by default
if [[ $(echo {snakemake.input} | grep '\.combine_gvcf\.') ]]; then
    mem=20g
else
    mem=8g
fi
gatk_nonfree -Xmx$mem -Djava.io.tmpdir=$TMPDIR \
    --analysis_type GenotypeGVCFs \
    $arg_seq_dict \
    --out {snakemake.output.vcf} \
    --reference_sequence {snakemake.config[static_data_config][reference][path]} \
    {arg_intervals} \
    $(for fname in $(echo {snakemake.input} | tr ' ' '\n' | grep '\.vcf.gz$'); do \
        echo -V $fname; \
    done)

# compute tabix index of the resulting VCF file
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf}) && \
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5 && \
    md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi}).md5
"""
)
