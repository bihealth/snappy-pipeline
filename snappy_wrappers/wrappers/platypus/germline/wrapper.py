# -*- coding: utf-8 -*-
"""Wrapper for running Platypus in germline
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

shell(
    r"""
set -x

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}

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

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Platypus Variant Calling ----------------------------------------------------

# TODO: interpret documentation for cores and ignore chroms

platypus callVariants \
    --logFileName=$(dirname {snakemake.log})/platypus.log \
    --bamFiles=$(echo "{snakemake.input}" \
                 | tr ' ' '\n' \
                 | grep -v '\.bai'$ \
                 | tr '\n' ',' \
                 | sed 's/,$//g') \
    --nCPU=16 \
    --refFile={snakemake.config[static_data_config][reference][path]} \
    --output={snakemake.output.vcf}
    # --output=$TMPDIR/temp_platypus.vcf

# # Fix Resulting VCF File ------------------------------------------------------

# snappy-fix_vcf \
#     --fix-platypus-gq \
#     --faidx {snakemake.config[static_data_config][reference][path]}.fai \
#     --input-vcf $TMPDIR/temp_platypus.vcf \
# | bgzip -c \
# > {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) > $(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.tbi}) > $(basename {snakemake.output.tbi}).md5
"""
)
