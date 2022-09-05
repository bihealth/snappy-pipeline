# -*- coding: utf-8 -*-
# isort:skip_file
from snappy_wrappers.tools import vcf_merge_exome_cnvs
import os

from snakemake.shell import shell


vcf_merge_exome_cnvs.main(
    [snakemake.output.vcf] + ["--sv-method", "GATK-4.1.7.0+gcnvkernel-0.7"] + list(snakemake.input)
)

shell(
    r"""
set -x
set -euo pipefail

cd $(dirname {snakemake.output.vcf})
tabix -f $(basename {snakemake.output.vcf})

md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf}.tbi) >$(basename {snakemake.output.vcf}).tbi.md5
"""
)
