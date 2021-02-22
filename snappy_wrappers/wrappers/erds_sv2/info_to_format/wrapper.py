# -*- coding: utf-8 -*-
"""Move SV2-specific INFO fields to FORMAT.
"""

import os.path

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bihealth.de"

# Get path to this file's (wrapper.py) directory.
base_dir = os.path.dirname(os.path.realpath(__file__))

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
# exec &> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

export LC_ALL=C

# Move FILTER/{{FILTER,DENOVO_FILTER,REF_GTL,AF}} to FORMAT/FT
python3 {base_dir}/info_to_format.py \
    --svmethod "ERDSv1.1+SV2v1.4.3.4" \
    --input {snakemake.input} \
    --output {snakemake.output.vcf}

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf}).tbi >$(basename {snakemake.output.vcf}).tbi.md5
"""
)
