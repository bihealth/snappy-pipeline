# -*- coding: utf-8 -*-
"""Wrapper for merging the ERDS call files before calling SV2 on these files.
"""

import os

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bihealth.de"

base_dir = os.path.dirname(os.path.realpath(__file__))

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec &> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

export LC_ALL=C
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

time python3 {base_dir}/merge_erds.py \
    --fai {snakemake.config[static_data_config][reference][path]}.fai \
    --output {snakemake.output.vcf} \
    --input {snakemake.input}

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf}).tbi >$(basename {snakemake.output.vcf}).tbi.md5
"""
)
